using BenchmarkTools

# --- Component Structs (our data) ---
struct P1; val::Int; end
struct P2; val::Int; end
struct P3; val::Int; end
struct P4; val::Int; end
struct P5; val::Int; end # This is the one we'll search for

# Sentinel for "not found"
struct NoP5 end
const NO_P5 = NoP5()

# --- Setup: Create the data containers ---

# 1. The Tuple container
# This mimics the `optionals` tuple in our proposed design.
const TUPLE_DATA = (P1(1), P2(2), P3(3), P4(4), P5(5))

# 2. The Dictionary container
# This mimics the `pdict` from the "Property Bag" design.
const DICT_DATA = Dict(
    :p1 => P1(1),
    :p2 => P2(2),
    :p3 => P3(3),
    :p4 => P4(4),
    :p5 => P5(5)  # The key we'll search for
)

# --- Method 1: Dictionary Lookup ---
function get_from_dict(d, key, default)
    return get(d, key, default)
end

# --- Method 2: Iterative Tuple Search ---
@inline function get_iterative(tup::Tuple, ::Type{T}, default) where {T}
    for component in tup
        if component isa T
            return component
        end
    end
    return default
end

# --- Method 3: Recursive Tuple Search ---
@inline function get_recursive(tup::Tuple{}, ::Type{T}, default) where {T}
    return default # Base case
end

@inline function get_recursive(tup::Tuple, ::Type{T}, default) where {T}
    if first(tup) isa T
        return first(tup)
    else
        # The @inline hint here is important for the compiler
        return @inline get_recursive(Base.tail(tup), T, default)
    end
end

# --- Benchmarking ---

# To ensure a fair comparison, we wrap the calls in a function
# so the compiler can specialize the code for the specific inputs.
run_dict() = get_from_dict(DICT_DATA, :p5, NO_P5)
run_iterative() = get_iterative(TUPLE_DATA, P5, NO_P5)
run_recursive() = get_recursive(TUPLE_DATA, P5, NO_P5)


println("--- Benchmarking Property Access ---")
println("\n1. Dictionary Method:")
@btime run_dict()

println("\n2. Iterative Tuple Search Method:")
@btime run_iterative()

println("\n3. Recursive Tuple Search Method:")
@btime run_recursive()


# --- Let's also check the generated code (for advanced users) ---
# The @code_native macro shows the actual machine code.
# For the recursive version, you'll see it's just a few simple instructions
# with no loops or branches.

println("\n--- Sample of Generated Machine Code for Recursive Method ---")
@code_native debuginfo=:none syntax=:intel run_recursive()
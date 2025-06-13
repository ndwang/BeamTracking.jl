using BenchmarkTools

# --- Define some simple structs to act as components ---
struct ComponentA; x::Int; end
struct ComponentB; y::Float64; end
struct ComponentC; z::Bool; end
struct ComponentD; w::Char; end
struct ComponentE; v::String; end

# --- The Tuple we will search through ---
const TUPLE_TO_SEARCH = (
    ComponentA(1), 
    ComponentB(2.0), 
    ComponentC(true), 
    ComponentD('d'), 
    ComponentE("hello")
)

# A default/sentinel value to return if nothing is found
const NOT_FOUND = nothing

# ==========================================================
# Method 1: Iterative `for` loop
# ==========================================================
@inline function get_component_iterative(tup::Tuple, ::Type{T}) where {T}
    for component in tup
        if component isa T
            return component
        end
    end
    return NOT_FOUND
end

# ==========================================================
# Method 2: Compile-time Recursive function
# ==========================================================
@inline function get_component_recursive(tup::Tuple, ::Type{T}) where {T}
    # This is the user-facing function that calls the internal recursive helper.
    # It provides the initial default value.
    _get_component_recursive_impl(tup, T, NOT_FOUND)
end

# The internal recursive implementation
@inline function _get_component_recursive_impl(tup::Tuple{}, ::Type{T}, default) where {T}
    return default # Base case for an empty tuple
end

@inline function _get_component_recursive_impl(tup::Tuple, ::Type{T}, default) where {T}
    # Compile-time check: is the first element's type what we want?
    if first(tup) isa T
        return first(tup)
    else
        # If not, recurse on the rest of the tuple. The compiler unrolls this.
        return _get_component_recursive_impl(Base.tail(tup), T, default)
    end
end


# ==========================================================
# Benchmarking
# ==========================================================

println("--- BENCHMARKING ---")
println("Tuple to search: ", typeof(TUPLE_TO_SEARCH))
println("-"^20)

# --- Best Case Scenario: Find the first element (ComponentA) ---
println("Best Case: Searching for ComponentA (first element)")
print("Iterative:  ")
@btime get_component_iterative($TUPLE_TO_SEARCH, ComponentA)
print("Recursive:  ")
@btime get_component_recursive($TUPLE_TO_SEARCH, ComponentA)
println()

# --- Worst Case Scenario: Find the last element (ComponentE) ---
println("Worst Case: Searching for ComponentE (last element)")
print("Iterative:  ")
@btime get_component_iterative($TUPLE_TO_SEARCH, ComponentE)
print("Recursive:  ")
@btime get_component_recursive($TUPLE_TO_SEARCH, ComponentE)
println()

# --- Not Found Scenario: Search for a type not in the tuple ---
println("Not Found Case: Searching for Float32 (not present)")
print("Iterative:  ")
@btime get_component_iterative($TUPLE_TO_SEARCH, Float32)
print("Recursive:  ")
@btime get_component_recursive($TUPLE_TO_SEARCH, Float32)
println()
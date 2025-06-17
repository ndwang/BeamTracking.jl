using BeamTracking
using BeamTracking: get_N_particle, runkernels!, MAX_TEMPS, KernelCall, KernelChain, BunchView
using BenchmarkTools
using SciMLBase, OrdinaryDiffEq
using StaticArrays

# Add CUDA support for GPU benchmarking
using CUDA

"""
    evaluate_kernel_performance(bunch, kernel, args...; n_runs=10)

Evaluate the performance of any tracking kernel and return a dictionary of metrics.

# Arguments
- `bunch`: Initial particle bunch
- `kernel`: The kernel function to evaluate
- `args...`: Arguments to pass to the kernel
- `n_runs`: Number of runs for performance evaluation (default: 10)

# Returns
A dictionary containing the following metrics:
- `min_time`: Minimum tracking time (in nanoseconds)
- `min_memory`: Minimum memory allocation (in bytes)
- `min_allocs`: Minimum number of allocations
- `success`: Boolean whether the tracking was successful

"""
function evaluate_kernel_performance(bunch, kernel, args...; n_runs=10, kwargs...)
    try
        # Create kernel chain
        kc = (KernelCall(kernel, args),)
        
        # Benchmark the tracking with specified sample size and time budget
        result = @benchmark begin
            runkernels!(nothing, $bunch, $kc; $(kwargs)...)
        end samples=n_runs seconds=10
        
        metrics = Dict(
            "min_time" => time(minimum(result)),
            "min_memory" => memory(minimum(result)),
            "min_allocs" => allocs(minimum(result)),
            "success" => true
        )

        return metrics
    catch e
        @warn "Tracking failed: $e"
        return Dict(
            "min_time" => NaN,
            "min_memory" => NaN,
            "min_allocs" => NaN,
            "success" => false
        )
    end 
end

"""
    create_gpu_bunch(n_particles; species=ELECTRON, Brho_ref=NaN)

Create a particle bunch with GPU arrays for benchmarking.

# Arguments
- `n_particles`: Number of particles in the bunch
- `species`: Particle species (default: ELECTRON)
- `Brho_ref`: Reference magnetic rigidity (default: NaN)

# Returns
A Bunch with GPU arrays, or nothing if CUDA is not available
"""
function create_gpu_bunch(n_particles)
    if !CUDA.functional()
        @warn "CUDA is not available, cannot create GPU bunch"
        return nothing
    end
        return Bunch(CUDA.zeros(n_particles, 6))
end

"""
    evaluate_gpu_kernel_performance(bunch, kernel, args...; n_runs=10, groupsize=256)

Evaluate the performance of any tracking kernel on GPU and return a dictionary of metrics.

# Arguments
- `bunch`: Initial particle bunch (should be on GPU)
- `kernel`: The kernel function to evaluate
- `args...`: Arguments to pass to the kernel
- `n_runs`: Number of runs for performance evaluation (default: 10)
- `groupsize`: GPU workgroup size (default: 256)

# Returns
A dictionary containing the following metrics:
- `min_time`: Minimum tracking time (in nanoseconds)
- `min_memory`: Minimum memory allocation (in bytes)
- `min_allocs`: Minimum number of allocations
- `success`: Boolean whether the tracking was successful
- `gpu_info`: Dictionary with GPU device information

"""
function evaluate_gpu_kernel_performance(bunch, kernel, args...; n_runs=10, groupsize=256, kwargs...)
    if !CUDA.functional()
        @warn "CUDA is not available"
        return Dict(
            "min_time" => NaN,
            "min_memory" => NaN,
            "min_allocs" => NaN,
            "success" => false,
            "gpu_info" => Dict("available" => false)
        )
    end
    
    try
        # Create kernel chain
        kc = (KernelCall(kernel, args),)
        
        # Ensure we're using KernelAbstractions for GPU
        kwargs_dict = Dict(kwargs...)
        kwargs_dict[:use_KA] = true
        kwargs_dict[:groupsize] = groupsize
        
        # Benchmark the tracking with specified sample size and time budget
        result = @benchmark begin
            runkernels!(nothing, $bunch, $kc; $(kwargs_dict)...)
        end samples=n_runs seconds=10
        
        # Get GPU device information
        device = CUDA.device()
        gpu_info = Dict(
            "available" => true,
            "name" => CUDA.name(device),
            "compute_capability" => string(CUDA.capability(device)),
            "memory_total" => CUDA.total_memory(),
            "memory_free" => CUDA.free_memory()
        )
        
        metrics = Dict(
            "min_time" => time(minimum(result)),
            "min_memory" => memory(minimum(result)),
            "min_allocs" => allocs(minimum(result)),
            "success" => true,
            "gpu_info" => gpu_info
        )

        return metrics
    catch e
        @warn "GPU tracking failed: $e"
        return Dict(
            "min_time" => NaN,
            "min_memory" => NaN,
            "min_allocs" => NaN,
            "success" => false,
            "gpu_info" => Dict("available" => false, "error" => string(e))
        )
    end 
end

function benchmark_kernel_cpu_vs_gpu(kernel, args...; n_particles=1000, n_runs=10, groupsize=256, kwargs...)
    # CPU benchmark
    cpu_bunch = Bunch(n_particles)
    cpu_metrics = evaluate_kernel_performance(BunchView(cpu_bunch), kernel, args...; n_runs=n_runs, kwargs...)
    
    # GPU benchmark
    gpu_bunch = create_gpu_bunch(n_particles)
    gpu_metrics = evaluate_gpu_kernel_performance(BunchView(gpu_bunch), kernel, args...; n_runs=n_runs, groupsize=groupsize, kwargs...)
    
    # Calculate speedup
    speedup = cpu_metrics["success"] && gpu_metrics["success"] ? 
              cpu_metrics["min_time"] / gpu_metrics["min_time"] : NaN
    
    return Dict("cpu" => cpu_metrics, "gpu" => gpu_metrics, "speedup" => speedup)
end

function evaluate_field_track_performance(; n_runs=10, n_particles=1000, solver=Tsit5(), solver_params=(save_everystep=false,save_start=false,save_end=true,dense=false,calck=false))
    bunch = Bunch(n_particles)
    L = 1.0
    field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
    params = nothing
    return evaluate_kernel_performance(BunchView(bunch), FieldTracking.field_track!, L, field_func, params, solver, solver_params; n_runs=n_runs)
end

function evaluate_linear_track_performance(;n_runs=10, n_particles=1000)
    # suggest good default values for bunch, L, r56
    bunch = Bunch(n_particles)
    L = 1.0
    r56 = 1.0
    return evaluate_kernel_performance(BunchView(bunch), LinearTracking.linear_drift!, L, r56; n_runs=n_runs)
end

function evaluate_rk4_track_performance(;n_runs=10, n_particles=1000)
    bunch = Bunch(n_particles)
    t_span = (0.0, 1.0)
    field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
    params = nothing
    return evaluate_kernel_performance(BunchView(bunch), RungeKuttaTracking.rk4_track!, t_span, field_func, params, 10; n_runs=n_runs)
end

# GPU-specific benchmark functions
"""
    evaluate_field_track_gpu_performance(; n_runs=10, n_particles=1000, groupsize=256, solver=Tsit5(), solver_params=(save_everystep=false,save_start=false,save_end=true,dense=false,calck=false))
"""
function evaluate_field_track_gpu_performance(; n_runs=10, n_particles=1000, groupsize=256, solver=Tsit5(), solver_params=(save_everystep=false,save_start=false,save_end=true,dense=false,calck=false))
    gpu_bunch = create_gpu_bunch(n_particles)
    if isnothing(gpu_bunch)
        return Dict("success" => false, "error" => "Failed to create GPU bunch")
    end
    
    L = 1.0
    field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
    params = nothing
    return evaluate_gpu_kernel_performance(BunchView(gpu_bunch), FieldTracking.field_track!, L, field_func, params, solver, solver_params; n_runs=n_runs, groupsize=groupsize)
end

"""
    evaluate_linear_track_gpu_performance(; n_runs=10, n_particles=1000, groupsize=256)
"""
function evaluate_linear_track_gpu_performance(; n_runs=10, n_particles=1000, groupsize=256)
    gpu_bunch = create_gpu_bunch(n_particles)
    if isnothing(gpu_bunch)
        return Dict("success" => false, "error" => "Failed to create GPU bunch")
    end
    
    L = 1.0
    r56 = 1.0
    return evaluate_gpu_kernel_performance(BunchView(gpu_bunch), LinearTracking.linear_drift!, L, r56; n_runs=n_runs, groupsize=groupsize)
end

"""
    evaluate_rk4_track_gpu_performance(; n_runs=10, n_particles=1000, groupsize=256)
"""
function evaluate_rk4_track_gpu_performance(; n_runs=10, n_particles=1000, groupsize=256)
    gpu_bunch = create_gpu_bunch(n_particles)
    if isnothing(gpu_bunch)
        return Dict("success" => false, "error" => "Failed to create GPU bunch")
    end
    
    t_span = (0.0, 1.0)
    field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
    params = nothing
    return evaluate_gpu_kernel_performance(BunchView(gpu_bunch), RungeKuttaTracking.rk4_track!, t_span, field_func, params, 10; n_runs=n_runs, groupsize=groupsize)
end

"""
    run_comprehensive_gpu_benchmark(; n_particles_range=[100, 1000, 10000], n_runs=10, groupsize=256)

Run a comprehensive benchmark comparing CPU vs GPU performance across different particle counts.

# Arguments
- `n_particles_range`: Array of particle counts to test (default: [100, 1000, 10000])
- `n_runs`: Number of benchmark runs per test (default: 10)
- `groupsize`: GPU workgroup size (default: 256)

# Returns
A dictionary containing benchmark results for all kernels and particle counts
"""
function run_comprehensive_gpu_benchmark(; n_particles_range=[100, 1000, 10000], n_runs=10, groupsize=256)
    results = Dict()
    
    # Test linear drift kernel
    results["linear_drift"] = Dict()
    for n_particles in n_particles_range
        results["linear_drift"][n_particles] = benchmark_kernel_cpu_vs_gpu(
            LinearTracking.linear_drift!, 1.0, 1.0; 
            n_particles=n_particles, n_runs=n_runs, groupsize=groupsize
        )
    end
    
    # Test RK4 kernel
    results["rk4_track"] = Dict()
    t_span = (0.0, 1.0)
    field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
    params = nothing
    for n_particles in n_particles_range
        results["rk4_track"][n_particles] = benchmark_kernel_cpu_vs_gpu(
            RungeKuttaTracking.rk4_track!, t_span, field_func, params, 10; 
            n_particles=n_particles, n_runs=n_runs, groupsize=groupsize
        )
    end
    
    #=
    # Test field tracking kernel (if available)
    try
        results["field_track"] = Dict()
        solver = Tsit5()
        solver_params = (save_everystep=false, save_start=false, save_end=true, dense=false, calck=false)
        field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
        params = nothing
        for n_particles in n_particles_range
            results["field_track"][n_particles] = benchmark_kernel_cpu_vs_gpu(
                FieldTracking.field_track!, 1.0, field_func, params, solver, solver_params; 
                n_particles=n_particles, n_runs=n_runs, groupsize=groupsize
            )
        end
    catch e
        @warn "Field tracking kernel not available: $e"
    end
    =#
    
    return results
end

"""
    print_benchmark_results(results)

Print benchmark results in a formatted way.

# Arguments
- `results`: Results from run_comprehensive_gpu_benchmark
"""
function print_benchmark_results(results)
    println("="^80)
    println("GPU vs CPU Kernel Performance Benchmark Results")
    println("="^80)
    
    for (kernel_name, particle_results) in results
        println("\n$(uppercase(kernel_name)) KERNEL:")
        println("-"^50)
        
        for (n_particles, result) in particle_results
            println("\nParticles: $n_particles")
            
            if result["cpu"]["success"]
                cpu_time = result["cpu"]["min_time"] / 1e6  # Convert to ms
                println("  CPU:   $(round(cpu_time, digits=3)) ms")
            else
                println("  CPU:   Failed")
            end
            
            if result["gpu"]["success"]
                gpu_time = result["gpu"]["min_time"] / 1e6  # Convert to ms
                println("  GPU:   $(round(gpu_time, digits=3)) ms")
                
                if !isnan(result["speedup"])
                    println("  Speedup: $(round(result["speedup"], digits=2))x")
                end
                
                if haskey(result["gpu"]["gpu_info"], "name")
                    println("  GPU Device: $(result["gpu"]["gpu_info"]["name"])")
                end
            else
                println("  GPU:   Failed")
            end
        end
    end
    
    println("\n" * "="^80)
end

#=
function example()
    # Basic GPU vs CPU comparison
    result = benchmark_kernel_cpu_vs_gpu(
        LinearTracking.linear_drift!, 1.0, 1.0; 
        n_particles=1000, n_runs=10
    )

    # Comprehensive benchmark across multiple particle counts
    results = run_comprehensive_gpu_benchmark(
        n_particles_range=[100, 1000, 10000],
        n_runs=10,
        groupsize=256
    )
end
=#
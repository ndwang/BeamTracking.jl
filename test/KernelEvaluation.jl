using BeamTracking: get_N_particle, runkernel!, MAX_TEMPS, soaview
using BenchmarkTools

"""
    evaluate_kernel_performance(bunch, kernel, args...; n_runs=10, kwargs...)

Evaluate the performance of any tracking kernel and return a dictionary of metrics.

# Arguments
- `bunch`: Initial particle bunch
- `kernel`: The kernel function to evaluate
- `args...`: Arguments to pass to the kernel
- `n_runs`: Number of runs for performance evaluation (default: 10)
- `kwargs...`: Additional keyword arguments to pass to runkernel!

# Returns
A dictionary containing the following metrics:
- `min_time`: Minimum tracking time per particle
- `min_memory`: Minimum memory allocation per particle
- `min_allocs`: Minimum number of allocations per particle
- `success`: Boolean whether the tracking was successful

"""
function evaluate_kernel_performance(bunch, kernel, args...; n_runs=10, kwargs...)
    n_particles = get_N_particle(bunch)

    # Get the tracking method from the kernel's module
    tracking_method = parentmodule(kernel).TRACKING_METHOD()
    n_temps = MAX_TEMPS(tracking_method)
    work = zeros(eltype(bunch.v), n_particles, n_temps)
    v = soaview(bunch)

    try
        # Benchmark the tracking with specified sample size and time budget
        result = @benchmark begin
            runkernel!($kernel, nothing, $v, $work, $(args...); $(kwargs...))
        end samples=n_runs seconds=10
        
        metrics = Dict(
            "min_time" => time(minimum(result)) / n_particles,
            "min_memory" => memory(minimum(result)) / n_particles,
            "min_allocs" => allocs(minimum(result)) / n_particles,
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
    evaluate_field_track_performance(bunch, L, field_func, params, solver; n_runs=10)

Evaluate the performance of field-based particle tracking and return detailed metrics.

# Arguments
- `bunch`: Initial particle bunch to be tracked
- `L`: Drift length for the tracking simulation
- `field_func`: Function that returns the field at a given position (x, y, z)
- `params`: Additional parameters for the field function
- `solver`: ODE solver to use for the integration
- `n_runs`: Number of runs for performance evaluation (default: 10)

# Returns
A dictionary containing the following metrics:
- `min_time`: Minimum tracking time per particle
- `min_memory`: Minimum memory allocation per particle
- `min_allocs`: Minimum number of allocations per particle
- `success`: Boolean whether the tracking was successful

"""
function evaluate_field_track_performance(bunch, L, field_func, params, solver; n_runs=10)
    return evaluate_kernel_performance(bunch, field_track!, L, field_func, params, solver; n_runs=n_runs)
end
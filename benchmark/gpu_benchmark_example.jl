#!/usr/bin/env julia

"""
GPU Benchmark Example for BeamTracking

This script demonstrates how to use the GPU benchmark functions to compare
CPU vs GPU performance for various tracking kernels.

Usage:
    julia gpu_benchmark_example.jl

Requirements:
    - CUDA.jl package installed
    - NVIDIA GPU with CUDA support
"""

using BeamTracking
using BeamTracking: LinearTracking, RungeKuttaTracking, FieldTracking
using BenchmarkTools
using StaticArrays
using CUDA

# Include the benchmark functions
include("KernelEvaluation.jl")

function main()
    println("BeamTracking GPU Benchmark Example")
    println("="^50)
    
    # Check if CUDA is available
    if !CUDA.functional()
        println("CUDA is not available on this system")
        println("This example requires an NVIDIA GPU with CUDA support")
        return
    end
    
    # Display GPU information
    device = CUDA.device()
    println("CUDA is available")
    println("GPU: $(CUDA.name(device))")
    println("Compute Capability: $(CUDA.capability(device))")
    println("Memory: $(round(CUDA.totalmem(device) / 1024^3, digits=1)) GB")
    println()
    
    # Test individual kernels
    println("Testing individual kernels...")
    println("-"^30)
    
    # Test linear drift kernel
    println("1. Linear Drift Kernel")
    linear_result = benchmark_kernel_cpu_vs_gpu(
        LinearTracking.linear_drift!, 1.0, 1.0; 
        n_particles=1000, n_runs=5
    )
    
    if linear_result["cpu"]["success"] && linear_result["gpu"]["success"]
        cpu_time = linear_result["cpu"]["min_time"] / 1e6
        gpu_time = linear_result["gpu"]["min_time"] / 1e6
        speedup = linear_result["speedup"]
        println("CPU: $(round(cpu_time, digits=3)) ms")
        println("GPU: $(round(gpu_time, digits=3)) ms")
        println("Speedup: $(round(speedup, digits=2))x")
    else
        println("Benchmark failed")
    end
    println()
    
    # Test RK4 kernel
    println("2. RK4 Tracking Kernel")
    t_span = (0.0, 1.0)
    field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
    params = nothing
    rk4_result = benchmark_kernel_cpu_vs_gpu(
        RungeKuttaTracking.rk4_track!, t_span, field_func, params, 10; 
        n_particles=1000, n_runs=5
    )
    
    if rk4_result["cpu"]["success"] && rk4_result["gpu"]["success"]
        cpu_time = rk4_result["cpu"]["min_time"] / 1e6
        gpu_time = rk4_result["gpu"]["min_time"] / 1e6
        speedup = rk4_result["speedup"]
        println("CPU: $(round(cpu_time, digits=3)) ms")
        println("GPU: $(round(gpu_time, digits=3)) ms")
        println("Speedup: $(round(speedup, digits=2))x")
    else
        println("Benchmark failed")
    end
    println()
    
    #=
    # Test field tracking kernel (if available)
    println("3. Field Tracking Kernel")
    try
        solver = Tsit5()
        solver_params = (save_everystep=false, save_start=false, save_end=true, dense=false, calck=false)
        field_func = (u, t, params) -> SVector(u[2], 0.0, u[4], 0.0, u[6], 0.0)
        params = nothing
        field_result = benchmark_kernel_cpu_vs_gpu(
            FieldTracking.field_track!, 1.0, field_func, params, solver, solver_params; 
            n_particles=1000, n_runs=5
        )
        
        if field_result["cpu"]["success"] && field_result["gpu"]["success"]
            cpu_time = field_result["cpu"]["min_time"] / 1e6
            gpu_time = field_result["gpu"]["min_time"] / 1e6
            speedup = field_result["speedup"]
            println("CPU: $(cpu_time:.3f) ms")
            println("GPU: $(gpu_time:.3f) ms")
            println("Speedup: $(speedup:.2f)x")
        else
            println("Benchmark failed")
        end
    catch e
        println("Field tracking kernel not available: $e")
    end
    println()
    =#
    # Run comprehensive benchmark
    println("Running comprehensive benchmark...")
    println("-"^30)
    
    # Use smaller particle counts for faster demonstration
    comprehensive_results = run_comprehensive_gpu_benchmark(
        n_particles_range=[100, 1000, 10000, 100000], 
        n_runs=3,  # Fewer runs for faster execution
        groupsize=256
    )
    
    # Print results
    print_benchmark_results(comprehensive_results)
    
    println("\nGPU benchmark completed successfully!")
end

# Run the example if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end 
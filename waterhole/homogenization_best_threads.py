#!/usr/bin/env python3
"""
Benchmark script to find optimal number of threads for FITS image processing.
This script tests different thread counts for loading and processing multiple FITS images.
"""

import os
import sys
import time
import psutil
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
from astropy.io import fits
import glob
import gc

def msg(message):
    """Print timestamped message with flush"""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"{timestamp} | {message}", flush=True)

def get_image(fitsfile):
    """Load FITS image data, handling degenerate axes properly"""
    try:
        input_hdu = fits.open(fitsfile)[0]
        # Always use the last 2 dimensions regardless of number of axes
        data = input_hdu.data
        while len(data.shape) > 2:
            data = data[0]  # Remove first axis until we have 2D
        return np.array(data, dtype=np.float32)
    except Exception as e:
        msg(f"Error loading {fitsfile}: {e}")
        return None

def get_system_info():
    """Get system resource information"""
    cpu_physical = psutil.cpu_count(logical=False)
    cpu_logical = psutil.cpu_count(logical=True)
    memory = psutil.virtual_memory()
    
    msg("System Information:")
    msg(f"  Physical CPU cores: {cpu_physical}")
    msg(f"  Logical CPU cores: {cpu_logical}")
    msg(f"  Total RAM: {memory.total / (1024**3):.1f} GB")
    msg(f"  Available RAM: {memory.available / (1024**3):.1f} GB")
    
    return {
        'cpu_physical': cpu_physical,
        'cpu_logical': cpu_logical,
        'memory_total_gb': memory.total / (1024**3),
        'memory_available_gb': memory.available / (1024**3)
    }

def process_chunk_threaded(args):
    """Process a single chunk using multiple images"""
    chunk_id, images, y_start, y_end, x_start, x_end = args
    
    try:
        chunk_height = y_end - y_start
        chunk_width = x_end - x_start
        chunk_data = np.zeros((len(images), chunk_height, chunk_width), dtype=np.float32)
        
        # Load chunk from all images
        for img_idx, img_path in enumerate(images):
            img = get_image(img_path)
            if img is not None:
                chunk_data[img_idx] = img[y_start:y_end, x_start:x_end]
            del img
        
        # Compute median
        result_chunk = np.median(chunk_data, axis=0)
        del chunk_data
        
        return chunk_id, result_chunk, y_start, y_end, x_start, x_end
        
    except Exception as e:
        msg(f"Error processing chunk {chunk_id}: {e}")
        return chunk_id, None, y_start, y_end, x_start, x_end

def benchmark_thread_count(images, n_threads, chunk_size=512, n_test_chunks=20):
    """
    Benchmark processing with specific number of threads
    
    Parameters
    ----------
    images : list
        List of image file paths
    n_threads : int
        Number of threads to test
    chunk_size : int
        Size of chunks to process
    n_test_chunks : int
        Number of test chunks to process (for speed)
    
    Returns
    -------
    dict
        Timing and performance results
    """
    # Get image dimensions from first image
    first_img = get_image(images[0])
    if first_img is None:
        return None
        
    ny, nx = first_img.shape
    del first_img
    
    # Create test chunks (just process a subset for benchmarking)
    chunk_tasks = []
    chunk_id = 0
    
    for i in range(0, min(ny, n_test_chunks * chunk_size), chunk_size):
        for j in range(0, min(nx, n_test_chunks * chunk_size), chunk_size):
            y_start, y_end = i, min(i + chunk_size, ny)
            x_start, x_end = j, min(j + chunk_size, nx)
            
            chunk_tasks.append((chunk_id, images, y_start, y_end, x_start, x_end))
            chunk_id += 1
            
            if len(chunk_tasks) >= n_test_chunks:
                break
        if len(chunk_tasks) >= n_test_chunks:
            break
    
    # Benchmark processing with specified thread count
    start_time = time.time()
    memory_start = psutil.virtual_memory()
    
    completed_chunks = 0
    
    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = [executor.submit(process_chunk_threaded, task) for task in chunk_tasks]
        
        for future in as_completed(futures):
            result = future.result()
            if result[1] is not None:  # Check if chunk was processed successfully
                completed_chunks += 1
    
    end_time = time.time()
    memory_end = psutil.virtual_memory()
    
    elapsed_time = end_time - start_time
    chunks_per_second = completed_chunks / elapsed_time if elapsed_time > 0 else 0
    memory_used_mb = (memory_end.used - memory_start.used) / (1024**2)
    
    # Force garbage collection
    gc.collect()
    
    return {
        'n_threads': n_threads,
        'elapsed_time': elapsed_time,
        'completed_chunks': completed_chunks,
        'chunks_per_second': chunks_per_second,
        'memory_used_mb': memory_used_mb,
        'efficiency': chunks_per_second / n_threads  # Chunks per second per thread
    }

def find_optimal_threads_for_image_count(images, image_counts, max_threads=None, chunk_size=512, n_test_chunks=20):
    """
    Find optimal number of threads for different image counts
    
    Parameters
    ----------
    images : list
        Full list of image file paths
    image_counts : list
        Different numbers of images to test with
    max_threads : int, optional
        Maximum number of threads to test
    chunk_size : int
        Size of chunks to process
    n_test_chunks : int
        Number of test chunks for benchmarking
    
    Returns
    -------
    dict
        Results for all tested combinations
    """
    system_info = get_system_info()
    
    if max_threads is None:
        max_threads = min(32, system_info['cpu_logical'])
    
    # Test different thread counts
    thread_counts = [1, 2, 4, 6, 8, 12, 16, 20, 24, 28, 32]
    thread_counts = [t for t in thread_counts if t <= max_threads]
    
    msg(f"Benchmarking thread counts: {thread_counts}")
    msg(f"Testing with image counts: {image_counts}")
    msg(f"Chunk size: {chunk_size}x{chunk_size}, Test chunks: {n_test_chunks}")
    msg("=" * 80)
    
    all_results = {}
    
    for num_images in image_counts:
        if num_images > len(images):
            msg(f"Warning: Requested {num_images} images but only {len(images)} available")
            num_images = len(images)
        
        test_images = images[:num_images]
        msg(f"\nTesting with {num_images} images")
        msg("-" * 40)
        
        results = {}
        
        for n_threads in thread_counts:
            msg(f"Testing {n_threads} threads with {num_images} images...")
            
            result = benchmark_thread_count(test_images, n_threads, chunk_size, n_test_chunks)
            
            if result is not None:
                results[n_threads] = result
                
                msg(f"  Time: {result['elapsed_time']:.2f}s, "
                    f"Rate: {result['chunks_per_second']:.2f} chunks/s, "
                    f"Efficiency: {result['efficiency']:.3f} chunks/s/thread")
            else:
                msg(f"  Failed to benchmark {n_threads} threads")
            
            # Wait between tests to let system settle
            time.sleep(1)
        
        all_results[num_images] = results
    
    return all_results

def analyze_results_by_image_count(all_results):
    """Analyze benchmark results across different image counts"""
    msg("\n" + "=" * 80)
    msg("COMPREHENSIVE BENCHMARK RESULTS")
    msg("=" * 80)
    
    recommendations = {}
    
    for num_images in sorted(all_results.keys()):
        results = all_results[num_images]
        if not results:
            continue
            
        msg(f"\nResults for {num_images} images:")
        msg("-" * 50)
        
        # Print detailed results table
        msg(f"{'Threads':<8} {'Time(s)':<8} {'Rate(c/s)':<10} {'Efficiency':<12} {'Memory(MB)':<12} {'Speedup':<8}")
        msg("-" * 70)
        
        baseline_time = None
        best_rate = 0
        best_threads = 1
        best_efficiency = 0
        best_efficiency_threads = 1
        
        for n_threads in sorted(results.keys()):
            result = results[n_threads]
            
            if baseline_time is None:
                baseline_time = result['elapsed_time']
                speedup = 1.0
            else:
                speedup = baseline_time / result['elapsed_time']
            
            # Track best performance
            if result['chunks_per_second'] > best_rate:
                best_rate = result['chunks_per_second']
                best_threads = n_threads
            
            if result['efficiency'] > best_efficiency:
                best_efficiency = result['efficiency']
                best_efficiency_threads = n_threads
            
            msg(f"{n_threads:<8} {result['elapsed_time']:<8.2f} {result['chunks_per_second']:<10.2f} "
                f"{result['efficiency']:<12.3f} {result['memory_used_mb']:<12.1f} {speedup:<8.2f}")
        
        # Store recommendation for this image count
        recommendations[num_images] = {
            'best_performance_threads': best_threads,
            'best_efficiency_threads': best_efficiency_threads,
            'best_rate': best_rate,
            'best_efficiency': best_efficiency
        }
        
        msg(f"\nBest performance: {best_threads} threads ({best_rate:.2f} chunks/s)")
        msg(f"Best efficiency: {best_efficiency_threads} threads ({best_efficiency:.3f} chunks/s/thread)")
    
    # Overall analysis
    msg("\n" + "=" * 80)
    msg("OVERALL ANALYSIS")
    msg("=" * 80)
    
    # Check if optimal thread count changes with image count
    perf_threads = [rec['best_performance_threads'] for rec in recommendations.values()]
    eff_threads = [rec['best_efficiency_threads'] for rec in recommendations.values()]
    
    if len(set(perf_threads)) == 1:
        msg(f"Optimal performance thread count is consistent: {perf_threads[0]} threads")
    else:
        msg("Optimal performance thread count varies with image count:")
        for num_images, rec in recommendations.items():
            msg(f"  {num_images} images: {rec['best_performance_threads']} threads")
    
    if len(set(eff_threads)) == 1:
        msg(f"Optimal efficiency thread count is consistent: {eff_threads[0]} threads")
    else:
        msg("Optimal efficiency thread count varies with image count:")
        for num_images, rec in recommendations.items():
            msg(f"  {num_images} images: {rec['best_efficiency_threads']} threads")
    
    # Determine if image count affects optimal threading
    if len(set(perf_threads)) > 1 or len(set(eff_threads)) > 1:
        msg("\nCONCLUSION: Optimal thread count DOES depend on number of images")
        msg("Recommendation: Use adaptive threading based on image count")
    else:
        msg(f"\nCONCLUSION: Optimal thread count is INDEPENDENT of image count")
        msg(f"Recommendation: Use {perf_threads[0]} threads for best performance")
        msg(f"Alternative: Use {eff_threads[0]} threads for best efficiency")
    
    return recommendations

def main():
    """Main function"""
    if len(sys.argv) not in [2, 3, 4, 5, 6]:
        msg("Usage: python3 homogenization_best_threads.py <image_pattern> [max_images] [max_threads] [chunk_size] [test_counts]")
        msg("")
        msg("Arguments:")
        msg("  image_pattern: Glob pattern for FITS images (e.g., '*.fits')")
        msg("  max_images: Maximum number of images to use (default: 512)")
        msg("  max_threads: Maximum number of threads to test (default: auto-detect)")
        msg("  chunk_size: Size of chunks for testing (default: 512)")
        msg("  test_counts: Comma-separated list of image counts to test (default: '50,100,200,512')")
        msg("")
        msg("Examples:")
        msg("  python3 homogenization_best_threads.py '*.fits'")
        msg("  python3 homogenization_best_threads.py '*.fits' 512 16 1024 '50,100,200,400'")
        sys.exit(1)
    
    # Parse arguments
    image_pattern = sys.argv[1]
    max_images = int(sys.argv[2]) if len(sys.argv) > 2 else 512
    max_threads = int(sys.argv[3]) if len(sys.argv) > 3 else None
    chunk_size = int(sys.argv[4]) if len(sys.argv) > 4 else 512
    test_counts_str = sys.argv[5] if len(sys.argv) > 5 else "50,100,200,512"
    
    # Parse test counts
    try:
        test_counts = [int(x.strip()) for x in test_counts_str.split(',')]
    except ValueError:
        msg(f"ERROR: Invalid test_counts format: {test_counts_str}")
        msg("Use comma-separated integers like: '50,100,200,512'")
        sys.exit(1)
    
    # Find images
    msg(f"Searching for images matching: {image_pattern}")
    images = glob.glob(f'{image_pattern}-[!MFS]*-I-image.fits')
    
    if not images:
        msg(f"ERROR: No images found matching pattern: {image_pattern}")
        sys.exit(1)
    
    # Limit number of images
    if len(images) > max_images:
        msg(f"Found {len(images)} images, using first {max_images}")
        images = images[:max_images]
    else:
        msg(f"Found {len(images)} images")
    
    # Filter test counts to available images
    test_counts = [count for count in test_counts if count <= len(images)]
    if not test_counts:
        msg(f"ERROR: No valid test counts for {len(images)} available images")
        sys.exit(1)
    
    # Verify first image can be loaded
    test_img = get_image(images[0])
    if test_img is None:
        msg(f"ERROR: Cannot load test image: {images[0]}")
        sys.exit(1)
    
    msg(f"Image dimensions: {test_img.shape}")
    del test_img
    
    # Run comprehensive benchmark
    msg(f"\nStarting comprehensive benchmark:")
    msg(f"  Total images available: {len(images)}")
    msg(f"  Image counts to test: {test_counts}")
    msg(f"  Max threads: {max_threads or 'auto-detect'}")
    msg(f"  Chunk size: {chunk_size}x{chunk_size}")
    
    all_results = find_optimal_threads_for_image_count(
        images, test_counts, max_threads, chunk_size
    )
    
    # Analyze and display results
    if all_results:
        recommendations = analyze_results_by_image_count(all_results)
        
        # Save results to file
        import json
        output_file = "thread_benchmark_comprehensive.json"
        
        # Convert results for JSON serialization
        json_results = {}
        for img_count, results in all_results.items():
            json_results[str(img_count)] = {}
            for thread_count, result in results.items():
                json_results[str(img_count)][str(thread_count)] = {
                    k: float(v) if isinstance(v, (np.floating, np.integer)) else v 
                    for k, v in result.items()
                }
        
        # Add recommendations
        json_results['recommendations'] = {
            str(k): v for k, v in recommendations.items()
        }
        
        with open(output_file, 'w') as f:
            json.dump(json_results, f, indent=2)
        
        msg(f"\nDetailed results saved to: {output_file}")
    else:
        msg("No successful benchmark results obtained")

if __name__ == "__main__":
    main()
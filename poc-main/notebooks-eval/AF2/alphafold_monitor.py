#!/usr/bin/env python3
"""
AlphaFold Resource Monitor

This script monitors and logs the resource usage (CPU, RAM, GPU, VRAM) of an AlphaFold
Docker container during inference. It records metrics at regular intervals and
saves them to a CSV file for later analysis.

Requirements:
- docker Python package
- nvidia-ml-py3 (for GPU monitoring)
- matplotlib and pandas (for optional visualization)

Usage:
python alphafold_monitor.py --container_name alphafold --interval 5 --output_file alphafold_metrics.csv
"""

import argparse
import csv
import datetime
import json
import time
import os
import sys
import signal
import subprocess
from threading import Event

# Try to import optional dependencies
try:
    import docker
except ImportError:
    print("Warning: docker package not found. Install with: pip install docker")
    docker = None

try:
    import pynvml
except ImportError:
    print("Warning: nvidia-ml-py3 package not found. GPU monitoring will be disabled.")
    print("Install with: pip install nvidia-ml-py3")
    pynvml = None

# Optional visualization dependencies
try:
    import pandas as pd
    import matplotlib.pyplot as plt
    visualization_available = True
except ImportError:
    visualization_available = False
    print("Warning: pandas or matplotlib not found. Visualization features disabled.")

# Global variable to handle graceful shutdown
stop_monitoring = Event()

class ResourceMonitor:
    """Monitor resource usage of Docker containers and GPUs"""
    
    def __init__(self, container_name, interval=5, output_file="resource_metrics.csv"):
        """
        Initialize the resource monitor
        
        Args:
            container_name (str): Name of the Docker container to monitor
            interval (int): Monitoring interval in seconds
            output_file (str): Path to the output CSV file
        """
        self.container_name = container_name
        self.interval = interval
        self.output_file = output_file
        self.docker_client = None
        self.container = None
        self.gpu_available = False
        self.metrics = []
        
        # Check if Docker is available
        if docker:
            try:
                self.docker_client = docker.from_env()
                print(f"Connected to Docker daemon: {self.docker_client.version()['Version']}")
            except Exception as e:
                print(f"Error connecting to Docker: {e}")
        
        # Check if NVIDIA GPU monitoring is available
        if pynvml:
            try:
                pynvml.nvmlInit()
                self.gpu_count = pynvml.nvmlDeviceGetCount()
                self.gpu_available = self.gpu_count > 0
                if self.gpu_available:
                    print(f"Found {self.gpu_count} NVIDIA GPUs")
                else:
                    print("No NVIDIA GPUs found")
            except Exception as e:
                print(f"Error initializing NVIDIA Management Library: {e}")
        
        # Setup CSV file header
        self.setup_output_file()
    
    def setup_output_file(self):
        """Initialize the output CSV file with headers"""
        file_exists = os.path.isfile(self.output_file)
        
        with open(self.output_file, 'a', newline='') as f:
            writer = csv.writer(f)
            if not file_exists:
                headers = [
                    'timestamp', 'container_name', 'cpu_percent', 'memory_usage_mb', 
                    'memory_limit_mb', 'memory_percent'
                ]
                
                # Add GPU headers if available
                if self.gpu_available:
                    for i in range(self.gpu_count):
                        headers.extend([
                            f'gpu{i}_utilization', f'gpu{i}_memory_used_mb', 
                            f'gpu{i}_memory_total_mb', f'gpu{i}_memory_percent', 
                            f'gpu{i}_temperature'
                        ])
                
                writer.writerow(headers)
    
    def find_container(self):
        """Find the container by name"""
        if not self.docker_client:
            print("Docker client not available")
            return False
        
        try:
            containers = self.docker_client.containers.list()
            for container in containers:
                if self.container_name in container.name:
                    self.container = container
                    print(f"Found container: {container.name} (ID: {container.short_id})")
                    return True
            
            print(f"Container '{self.container_name}' not found. Available containers:")
            for container in containers:
                print(f"  - {container.name}")
            return False
        except Exception as e:
            print(f"Error finding container: {e}")
            return False
    
    def get_container_stats(self):
        """Get current container metrics"""
        metrics = {}
        
        if not self.container:
            if not self.find_container():
                return None
        
        try:
            # Get container stats
            stats = self.container.stats(stream=False)
            
            # Extract CPU metrics
            cpu_delta = stats['cpu_stats']['cpu_usage']['total_usage'] - \
                        stats['precpu_stats']['cpu_usage']['total_usage']
            system_delta = stats['cpu_stats']['system_cpu_usage'] - \
                          stats['precpu_stats']['system_cpu_usage']
            online_cpus = stats['cpu_stats'].get('online_cpus', len(stats['cpu_stats']['cpu_usage']['percpu_usage']))
            
            if system_delta > 0 and online_cpus > 0:
                cpu_percent = (cpu_delta / system_delta) * online_cpus * 100.0
            else:
                cpu_percent = 0.0
            
            # Extract memory metrics
            memory_usage = stats['memory_stats']['usage']
            memory_limit = stats['memory_stats']['limit']
            memory_percent = (memory_usage / memory_limit) * 100.0
            
            # Convert to MB
            memory_usage_mb = memory_usage / (1024 * 1024)
            memory_limit_mb = memory_limit / (1024 * 1024)
            
            metrics = {
                'timestamp': datetime.datetime.now().isoformat(),
                'container_name': self.container.name,
                'cpu_percent': round(cpu_percent, 2),
                'memory_usage_mb': round(memory_usage_mb, 2),
                'memory_limit_mb': round(memory_limit_mb, 2),
                'memory_percent': round(memory_percent, 2)
            }
            
            return metrics
        
        except Exception as e:
            print(f"Error getting container stats: {e}")
            # Try to reconnect
            self.container = None
            return None
    
    def get_gpu_stats(self):
        """Get current GPU metrics"""
        gpu_metrics = {}
        
        if not self.gpu_available or not pynvml:
            return gpu_metrics
        
        try:
            for i in range(self.gpu_count):
                handle = pynvml.nvmlDeviceGetHandleByIndex(i)
                
                # Get GPU utilization
                utilization = pynvml.nvmlDeviceGetUtilizationRates(handle)
                gpu_util = utilization.gpu
                
                # Get memory information
                memory_info = pynvml.nvmlDeviceGetMemoryInfo(handle)
                memory_used = memory_info.used / (1024 * 1024)  # Convert to MB
                memory_total = memory_info.total / (1024 * 1024)  # Convert to MB
                memory_percent = (memory_used / memory_total) * 100.0
                
                # Get temperature
                temperature = pynvml.nvmlDeviceGetTemperature(handle, pynvml.NVML_TEMPERATURE_GPU)
                
                gpu_metrics.update({
                    f'gpu{i}_utilization': gpu_util,
                    f'gpu{i}_memory_used_mb': round(memory_used, 2),
                    f'gpu{i}_memory_total_mb': round(memory_total, 2),
                    f'gpu{i}_memory_percent': round(memory_percent, 2),
                    f'gpu{i}_temperature': temperature
                })
        except Exception as e:
            print(f"Error getting GPU stats: {e}")
        
        return gpu_metrics
    
    def collect_metrics(self):
        """Collect all metrics and return combined dictionary"""
        container_metrics = self.get_container_stats()
        if not container_metrics:
            return None
        
        gpu_metrics = self.get_gpu_stats()
        container_metrics.update(gpu_metrics)
        
        return container_metrics
    
    def log_metrics(self, metrics):
        """Log metrics to CSV file"""
        if not metrics:
            return
        
        with open(self.output_file, 'a', newline='') as f:
            writer = csv.writer(f)
            
            # Create a row with values in the same order as headers
            row = []
            with open(self.output_file, 'r') as read_file:
                headers = next(csv.reader(read_file))
            
            for header in headers:
                row.append(metrics.get(header, ''))
            
            writer.writerow(row)
        
        # Also store metrics in memory for live reporting
        self.metrics.append(metrics)
    
    def print_current_metrics(self, metrics):
        """Print current metrics to console"""
        if not metrics:
            return
        
        print(f"\n--- {metrics['timestamp']} ---")
        print(f"Container: {metrics['container_name']}")
        print(f"CPU: {metrics['cpu_percent']}%")
        print(f"Memory: {metrics['memory_usage_mb']:.2f}MB / {metrics['memory_limit_mb']:.2f}MB ({metrics['memory_percent']:.2f}%)")
        
        # Print GPU metrics if available
        if self.gpu_available:
            for i in range(self.gpu_count):
                print(f"GPU {i}:")
                print(f"  Utilization: {metrics.get(f'gpu{i}_utilization', 'N/A')}%")
                print(f"  Memory: {metrics.get(f'gpu{i}_memory_used_mb', 'N/A'):.2f}MB / " + 
                      f"{metrics.get(f'gpu{i}_memory_total_mb', 'N/A'):.2f}MB " + 
                      f"({metrics.get(f'gpu{i}_memory_percent', 'N/A'):.2f}%)")
                print(f"  Temperature: {metrics.get(f'gpu{i}_temperature', 'N/A')}°C")
    
    def start_monitoring(self):
        """Start monitoring resources at the specified interval"""
        print(f"Starting monitoring of '{self.container_name}' every {self.interval} seconds...")
        print(f"Press Ctrl+C to stop monitoring")
        print(f"Logging to {self.output_file}")
        
        try:
            while not stop_monitoring.is_set():
                metrics = self.collect_metrics()
                if metrics:
                    self.log_metrics(metrics)
                    self.print_current_metrics(metrics)
                time.sleep(self.interval)
        
        except KeyboardInterrupt:
            print("\nMonitoring stopped by user")
        
        finally:
            self.cleanup()
    
    def cleanup(self):
        """Clean up resources"""
        print("\nCleaning up...")
        if pynvml and self.gpu_available:
            try:
                pynvml.nvmlShutdown()
            except:
                pass
    
    def generate_report(self):
        """Generate a report and visualizations from the collected data"""
        if not visualization_available:
            print("Visualization packages not available. Install pandas and matplotlib.")
            return
        
        if not os.path.exists(self.output_file):
            print(f"Output file {self.output_file} not found")
            return
        
        try:
            # Load data
            df = pd.read_csv(self.output_file)
            
            if len(df) == 0:
                print("No data collected")
                return
            
            # Convert timestamp to datetime for plotting
            df['timestamp'] = pd.to_datetime(df['timestamp'])
            
            # Create output directory for plots
            report_dir = "resource_report"
            os.makedirs(report_dir, exist_ok=True)
            
            # Plot CPU usage
            plt.figure(figsize=(12, 6))
            plt.plot(df['timestamp'], df['cpu_percent'])
            plt.title('CPU Usage')
            plt.xlabel('Time')
            plt.ylabel('CPU (%)')
            plt.grid(True)
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(os.path.join(report_dir, 'cpu_usage.png'))
            
            # Plot Memory usage
            plt.figure(figsize=(12, 6))
            plt.plot(df['timestamp'], df['memory_usage_mb'])
            plt.title('Memory Usage')
            plt.xlabel('Time')
            plt.ylabel('Memory (MB)')
            plt.grid(True)
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(os.path.join(report_dir, 'memory_usage.png'))
            
            # Plot GPU metrics if available
            gpu_columns = [col for col in df.columns if col.startswith('gpu') and 'utilization' in col]
            if gpu_columns:
                # Plot GPU utilization
                plt.figure(figsize=(12, 6))
                for col in gpu_columns:
                    gpu_id = col.split('_')[0][3:]
                    plt.plot(df['timestamp'], df[col], label=f'GPU {gpu_id}')
                plt.title('GPU Utilization')
                plt.xlabel('Time')
                plt.ylabel('Utilization (%)')
                plt.legend()
                plt.grid(True)
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig(os.path.join(report_dir, 'gpu_utilization.png'))
                
                # Plot GPU memory
                plt.figure(figsize=(12, 6))
                for i in range(self.gpu_count):
                    col = f'gpu{i}_memory_used_mb'
                    if col in df.columns:
                        plt.plot(df['timestamp'], df[col], label=f'GPU {i}')
                plt.title('GPU Memory Usage')
                plt.xlabel('Time')
                plt.ylabel('Memory (MB)')
                plt.legend()
                plt.grid(True)
                plt.xticks(rotation=45)
                plt.tight_layout()
                plt.savefig(os.path.join(report_dir, 'gpu_memory.png'))
            
            # Generate summary statistics
            summary = {
                'avg_cpu': df['cpu_percent'].mean(),
                'max_cpu': df['cpu_percent'].max(),
                'avg_memory_mb': df['memory_usage_mb'].mean(),
                'max_memory_mb': df['memory_usage_mb'].max(),
                'memory_limit_mb': df['memory_limit_mb'].iloc[0]
            }
            
            for i in range(self.gpu_count):
                util_col = f'gpu{i}_utilization'
                mem_col = f'gpu{i}_memory_used_mb'
                if util_col in df.columns:
                    summary[f'avg_gpu{i}_util'] = df[util_col].mean()
                    summary[f'max_gpu{i}_util'] = df[util_col].max()
                if mem_col in df.columns:
                    summary[f'avg_gpu{i}_memory_mb'] = df[mem_col].mean()
                    summary[f'max_gpu{i}_memory_mb'] = df[mem_col].max()
            
            # Save summary to file
            with open(os.path.join(report_dir, 'summary.txt'), 'w') as f:
                f.write("Resource Usage Summary\n")
                f.write("=====================\n\n")
                f.write(f"Monitoring period: {df['timestamp'].min()} to {df['timestamp'].max()}\n")
                f.write(f"Total samples: {len(df)}\n\n")
                
                f.write("CPU Usage:\n")
                f.write(f"  Average: {summary['avg_cpu']:.2f}%\n")
                f.write(f"  Maximum: {summary['max_cpu']:.2f}%\n\n")
                
                f.write("Memory Usage:\n")
                f.write(f"  Average: {summary['avg_memory_mb']:.2f} MB\n")
                f.write(f"  Maximum: {summary['max_memory_mb']:.2f} MB\n")
                f.write(f"  Limit: {summary['memory_limit_mb']:.2f} MB\n\n")
                
                if self.gpu_available:
                    f.write("GPU Usage:\n")
                    for i in range(self.gpu_count):
                        f.write(f"  GPU {i}:\n")
                        if f'avg_gpu{i}_util' in summary:
                            f.write(f"    Average Utilization: {summary[f'avg_gpu{i}_util']:.2f}%\n")
                            f.write(f"    Maximum Utilization: {summary[f'max_gpu{i}_util']:.2f}%\n")
                        if f'avg_gpu{i}_memory_mb' in summary:
                            f.write(f"    Average Memory: {summary[f'avg_gpu{i}_memory_mb']:.2f} MB\n")
                            f.write(f"    Maximum Memory: {summary[f'max_gpu{i}_memory_mb']:.2f} MB\n")
            
            print(f"Report generated in {report_dir} directory")
            
        except Exception as e:
            print(f"Error generating report: {e}")

def signal_handler(sig, frame):
    """Handle SIGINT (Ctrl+C) to gracefully stop monitoring"""
    print("\nReceived signal to stop monitoring")
    stop_monitoring.set()

def main():
    """Main function"""
    # Parse arguments
    parser = argparse.ArgumentParser(description='Monitor resources of AlphaFold docker container')
    parser.add_argument('--container_name', default='alphafold', 
                        help='Name of the AlphaFold docker container')
    parser.add_argument('--interval', type=int, default=5, 
                        help='Monitoring interval in seconds')
    parser.add_argument('--output_file', default='alphafold_metrics.csv', 
                        help='Output CSV file path')
    parser.add_argument('--report', action='store_true', 
                        help='Generate report from existing metrics file')
    
    args = parser.parse_args()
    
    # Setup signal handler for graceful shutdown
    signal.signal(signal.SIGINT, signal_handler)
    
    # Create monitor instance
    monitor = ResourceMonitor(
        container_name=args.container_name,
        interval=args.interval,
        output_file=args.output_file
    )
    
    # Generate report if requested
    if args.report:
        monitor.generate_report()
        return
    
    # Start monitoring
    monitor.start_monitoring()

if __name__ == "__main__":
    main()
import yaml
import argparse
import logging
import time
import os
import sys
import pandas as pd
import psutil
import threading
import time
import numpy as np

from tabulate import tabulate
from datetime import datetime
from subprocess import Popen, PIPE, STDOUT
from pathlib import Path
from pydantic import BaseModel
from typing import List

class Compressor(BaseModel):
    name: str
    ext: str
    executable: Path
    use_mpi: bool
    n_threads: List[int]
    compress_params: List[str]
    decompress_params: List[str]

class Dataset(BaseModel):
    name: str
    dimension: List[int]
    ext: str
    fileNames: List[str]
    folder: Path
    ebs: List[float]
    depths: List[int]
    data_type: str

class CompressionStats():
    compressor_name: str = ""
    n_threads: int = 1
    layer_depth: int = 1
    dataset_name: str = ""
    data_file_name: str = ""
    error_bound: float = 0
    num_of_elements: int = 0
    compress_wall_time: float = 0
    compress_cpu_time: float = 0
    compressed_size: int = 0
    original_size: int = 0
    compression_ratio: float = 0
    decompress_wall_time: float = 0
    decompress_cpu_time: float = 0

formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')


def setup_logger(name, log_file, level=logging.INFO):
    """To setup as many loggers as you want"""

    handler = logging.FileHandler(log_file)     
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

# current_time = datetime.now()
# log_prefix = f"logs/{current_time.month}-{current_time.day}-{current_time.year}_{current_time.hour}-{current_time.minute}-{current_time.second}"
# os.makedirs(log_prefix)

# logger = setup_logger('app_logger', f'{log_prefix}/app.log')
# verboseLogger = setup_logger('verbose_logger', f'{log_prefix}/verbose.log', logging.DEBUG)

class ResourceUsage(threading.Thread):
    def __init__(self, pid, filename, earlyStop=False):
        threading.Thread.__init__(self)
        self.p = psutil.Process(pid)
        self.memory_percents_list = []
        self.cpu_percents_list = []
        self.time_counts = []
        self.total_memory_list = []
        self.used_memory_list = []
        self.filename = filename
        self.earlyStop = earlyStop
    
    def run(self):
        time_count = 0
        df = pd.DataFrame({"time_count" : [0],
                           "memory_percents" : [0],
                           "total_memory": [0],
                            "used_memory": [0],
                           "cpu_percents": [0]})
        df.to_csv(self.filename, index=False)
        counter = 0
        while self.p.is_running():
            self.cpu_percents_list.append(psutil.cpu_percent())
            virtual_memory = psutil.virtual_memory()
            self.memory_percents_list.append(virtual_memory.percent)
            self.total_memory_list.append(virtual_memory.total)
            self.used_memory_list.append(virtual_memory.total * virtual_memory.percent / (1024 * 1024))
            self.time_counts.append(time_count)

            if counter >= 3:
                df = pd.DataFrame({"time_count" : self.time_counts,
                                        "memory_percents" : self.memory_percents_list,
                                        "total_memory": self.total_memory_list,
                                        "used_memory": self.used_memory_list,
                                        "cpu_percents": self.cpu_percents_list})
                df.to_csv(self.filename, mode='a', index=False, header=False)
                self.cpu_percents_list = []
                self.memory_percents_list = []
                self.time_counts = []
                self.total_memory_list = []
                self.used_memory_list = []
                counter = 0
            if self.earlyStop and time_count >= 500: # exit early for only 500s
                break
            time_count += 5
            counter += 1
            time.sleep(5)
            
        if len(self.cpu_percents_list) > 0:
            df = pd.DataFrame({"time_count" : self.time_counts,
                               "memory_percents" : self.memory_percents_list,
                               "total_memory": self.total_memory_list,
                               "used_memory": self.used_memory_list,
                               "cpu_percents": self.cpu_percents_list})
            df.to_csv(self.filename, mode='a', index=False, header=False)

        print("<ResourceUsage> Memory and CPU usage collection finished!")
        logger.info(f"<ResourceUsage> Memory and CPU usage collection finished! Data saved to {self.filename}")

def run_compress_and_collect(command_line_param_str: str, stats: CompressionStats, dataset:Dataset, data_file, compressor: Compressor, compressed_file, eb, dimension, depth=1):
    logger.info(f"start running compression for dataset <{dataset.name}> with compressor <{compressor.name}>")
    logger.info(command_line_param_str)
    start = time.perf_counter()
    process = Popen(' '.join(["time", command_line_param_str]), shell=True, stdout = PIPE, stderr=PIPE, text=True, executable='/bin/bash')
    resource_thread = ResourceUsage(pid=process.pid,
                                    filename=str(Path(log_prefix) / f"{dataset.name}-{eb}-depth-{depth}-{compressor.name}-compress-resources.csv"))
    resource_thread.start()
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            verboseLogger.debug(f"[compress]<{dataset.name}><{compressor.name}><{eb}> {output.strip()}")
    process.wait()
    resource_thread.join()
    stdout, stderr = process.communicate()
    logger.debug(f"[compress]<{dataset.name}><{compressor.name}> {stderr.strip()}")
    end = time.perf_counter()
    elapsed_time = end - start
    cpu_times = stderr.strip().split('\n')[-3:]  # real_time, user_time, sys_time
    real_time = cpu_times[1].split()[1]
    logger.info(f"finished running compression for dataset <{dataset.name}> with compressor <{compressor.name}> of error bound <{eb}>")
    # collect stats
    try:
        stats.compress_wall_time = elapsed_time
        stats.compress_cpu_time = real_time
        stats.error_bound = eb
        stats.data_file_name = Path(data_file).name
        stats.num_of_elements = np.prod(dimension)
        stats.compressed_size = os.path.getsize(compressed_file)
        stats.original_size = os.path.getsize(data_file)
        stats.compression_ratio = stats.original_size / stats.compressed_size
        stats.layer_depth = depth
    except Exception as error:
        logger.error("Cannot finish collecting the stats for compression!")
        logger.error(error)

def run_decompress_collect(command_line_param_str: str, stats: CompressionStats, dataset: Dataset, compressor: Compressor, eb, depth=1):
    logger.info(f"start running decompression for dataset <{dataset.name}> with compressor <{compressor.name}>")
    logger.info(command_line_param_str)
    start = time.perf_counter()
    process = Popen(' '.join(["time", command_line_param_str]), shell=True, stdout = PIPE, stderr=PIPE, text=True, executable='/bin/bash')
    resource_thread = ResourceUsage(pid=process.pid,
                                    filename=str(Path(log_prefix) / f"{dataset.name}-{eb}-depth-{depth}-{compressor.name}-decompress-resources.csv"))
    resource_thread.start()
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            verboseLogger.debug(f"[decompress]<{dataset.name}><{compressor.name}> {output.strip()}")
    process.wait()
    resource_thread.join()
    _, stderr = process.communicate()
    logger.debug(f"[decompress]<{dataset.name}><{compressor.name}> {stderr.strip()}")
    end = time.perf_counter()
    elapsed_time = end - start
    print("stderr splited: ", stderr.strip().split('\n'))
    cpu_times = stderr.strip().split('\n')[-3:]  # real_time, user_time, sys_time
    real_time = cpu_times[1].split()[1]
    stats.decompress_wall_time = elapsed_time
    stats.decompress_cpu_time = real_time

def sz3_split_compression(stats: CompressionStats, dataset: Dataset, data_file, compressor: Compressor, compressed_file: str, eb, dimension, depth, n_threads:int=16):
    params = compressor.compress_params.copy()

    for i, param in enumerate(params):
        if param == '$fileName':
            params[i] = data_file
        elif param == '$compressedFileName':
            params[i] = compressed_file
        elif param == '$eb':
            params[i] = str(eb)
    params = [str(compressor.executable)] + params + [f'-d'] + [str(dim) for dim in dimension] + ["--depth", str(depth)]
    params = params + [f"--data_type {dataset.data_type}"]
    if compressor.use_mpi:
        io_threads = max(2, n_threads / 8)
        params = ["mpiexec -n", str(n_threads)] + params + ["--threads", str(io_threads), "--mpi"]
    else:
        params = params + ["--threads", str(n_threads), "--mpi"]
    command_line_param_str = " ".join(params)
    print(command_line_param_str)
    run_compress_and_collect(command_line_param_str, stats, dataset, data_file, compressor, compressed_file, eb, dimension, depth)
   

def sz3_split_decompression(stats: CompressionStats, dataset: Dataset, compressor: Compressor, compressed_file: str, decompressed_file, eb, dimension, depth, n_threads:int=16):
    params = compressor.decompress_params.copy()
    for i, param in enumerate(params):
        if param == '$compressedFileName':
            params[i] = compressed_file
        elif param == '$decompressedFileName':
            params[i] = decompressed_file
        elif param == '$eb':
            params[i] = str(eb)
    # num_threads = 16 # min(64 / depth, 16)
    params = [str(compressor.executable)] + params + [f'-d'] + [str(dim) for dim in dimension] + ["--depth", str(depth), "--threads", str(n_threads)]
    params = params + [f"--data_type {dataset.data_type}"]
    if compressor.use_mpi:
        io_threads = max(2, n_threads / 8)
        params = ["mpiexec -n", str(n_threads)] + params + ["--threads", str(io_threads), "--mpi"]
    else:
        params = params + ["--threads", str(n_threads), "--mpi"]
    command_line_param_str = " ".join(params)
    print(command_line_param_str)
    run_decompress_collect(command_line_param_str, stats, dataset, compressor, eb, depth)



def sz3_compression(stats: CompressionStats, dataset: Dataset, data_file, compressor: Compressor, compressed_file: str, eb, dimension):
    params = compressor.compress_params.copy()

    for i, param in enumerate(params):
        if param == '$fileName':
            params[i] = data_file
        elif param == '$compressedFileName':
            params[i] = compressed_file
        elif param == '$eb':
            params[i] = str(eb)
    n_dim = len(dimension)
    params = [str(compressor.executable)] + params + [f'-{n_dim}'] + [str(dim) for dim in dimension]
    if dataset.data_type == "float64":
        params = params + ["-d"]
    else:
        params = params + ["-f"]
    command_line_param_str = " ".join(params)
    run_compress_and_collect(command_line_param_str, stats, dataset, data_file, compressor, compressed_file, eb, dimension)

def sz3_decompression(stats: CompressionStats, dataset: Dataset, compressor: Compressor, compressed_file, decompressed_file, eb, dimension):
    params = compressor.decompress_params.copy()
    for i, param in enumerate(params):
        if param == '$compressedFileName':
            params[i] = compressed_file
        elif param == '$decompressedFileName':
            params[i] = decompressed_file
        elif param == '$eb':
            params[i] = str(eb)
    n_dim = len(dimension)
    params = [str(compressor.executable)] + params + [f'-{n_dim}'] + [str(dim) for dim in dimension]
    if dataset.isFloat64:
        params = params + ["-d"]
    else:
        params = params + ["-f"]
    command_line_param_str = " ".join(params)
    run_decompress_collect(command_line_param_str, stats, dataset, compressor, eb)


def benchmark(config, do_compression: bool, do_decompression: bool):
    datasets = [Dataset(**dataset) for dataset in config["datasets"]]
    compressors = [Compressor(**compressor) for compressor in config["compressors"]]
    stats_file_path = f"{log_prefix}/benchmark_stats.csv"
    global_stats = {
        "compressor_name": [],
        "n_threads": [],
        "layer_depth": [],
        "dataset_name": [],
        "data_file_name": [],
        "error_bound": [],
        "num_of_elements": [],
        "compress_wall_time": [],
        "compress_cpu_time": [],
        "compressed_size": [],
        "original_size": [],
        "compression_ratio": [],
        "decompress_wall_time": [],
        "decompress_cpu_time": []
    }
    logger.info(f"start running benchmark on #{len(datasets)} datasets with #{len(compressors)} compressors!")
    for compressor in compressors:
        logger.info(f"Will test on the following number of threads: {compressor.n_threads}")
        for n_thread in compressor.n_threads:
            for dataset in datasets:
                logger.info(f"ebs: {dataset.ebs}")
                stats = CompressionStats()
                stats.n_threads = n_thread
                stats.compressor_name = compressor.name
                stats.dataset_name = dataset.name
                dataset_files = [dataset.folder / filename for filename in dataset.fileNames]
                for data_file_path in dataset_files:
                    for eb in dataset.ebs:
                        for depth in dataset.depths:
                            data_file = str(data_file_path)
                            filename = data_file_path.name
                            compressed_file = str(Path(config["global"]["large_file_output_folder"]) / (filename + '-' + str(eb) + '-' + compressor.name + '-' + str(n_thread) + '-threads-' + compressor.ext))
                            verboseLogger.info(f"compressed file: {compressed_file}")
                            if do_compression:     
                                try:
                                    sz3_split_compression(stats, dataset, data_file, compressor, compressed_file, eb, dataset.dimension, depth, n_thread)
                                except Exception as error:
                                    logger.error("Compression Failed:", error)
                                    sys.exit(1)
                            decompressed_file = compressed_file + ".dp"
                            if do_decompression:
                                try:
                                    sz3_split_decompression(stats, dataset, compressor, compressed_file, decompressed_file, eb, dataset.dimension, depth, n_thread)
                                except Exception as error:
                                    logger.error("Decompression Failed:", error)
                                    sys.exit(2)
                            for key in global_stats.keys():
                                global_stats[key].append(getattr(stats, key))
                            print(global_stats)
                            tmp_df = pd.DataFrame(global_stats)
                            print(tabulate(tmp_df, headers='keys', tablefmt='psql'))
                            logger.info(f"Tabulated partial results\n{tabulate(tmp_df, headers='keys', tablefmt='psql')}")
                            tmp_df.to_csv(stats_file_path, index=False)
                logger.info(f"finished benchmarking dataset <{dataset.name}> with compressor <{compressor.name}>")

def main():
    parser = argparse.ArgumentParser(
        description="Benchmarking lossy compression algorithms such as SZ3 and its variants",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "-c",
        "--config",
        type=str,
        default="anvil_benchmark_config.yml",
        help="The config file for benchmarking"
    )

    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        default="benchmark",
        help="There are three modes: benchmark, compress, decompress. `benchmark` mode will do compression and decopmression and collect all data, while the other two only collects partial data."
    )

    args = parser.parse_args()
    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    global log_prefix, verboseLogger, logger
    current_time = datetime.now()
    metrics_folder = config["global"]["metrics_output_folder"]
    log_prefix = Path(metrics_folder) / "logs" / f"{current_time.month}-{current_time.day}-{current_time.year}_{current_time.hour}-{current_time.minute}-{current_time.second}"
    log_prefix.mkdir(parents=True)
    large_file_folder = Path(config["global"]["large_file_output_folder"])
    large_file_folder.mkdir(parents=True, exist_ok=True)

    logger = setup_logger('app_logger', log_prefix / "app.log")
    verboseLogger = setup_logger('verbose_logger', log_prefix / "verbose.log", logging.DEBUG)
    logger.info(f"saving all metrics in {metrics_folder}")

    do_compression = True
    do_decompression = True
    if(args.mode == 'compress'):
        do_decompression = False
        logger.info("This program does compression ONLY and collects partial data")
    elif(args.mode == 'decompress'):
        do_compression = False
        logger.info("This program does decompression ONLY and collects partial data")
    else:
        logger.info("This program does full compression/decompression and will collect all data")

    logger.info("finished loading YAML config, ready for benchmarking!")
    benchmark(config, do_compression, do_decompression)
    logger.info("Congratulations! You have finished all benchmarking!")

if __name__ == '__main__':
    main()
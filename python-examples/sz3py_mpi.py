import sz3py
import numpy as np
import argparse
import os

from pathlib import Path
from PIL import Image
from matplotlib import pyplot as plt

def get_partial_preview_data(dimension: list[int], data_file: str, layer_number: int, is_float64: bool=False):
    import numpy as np
    from matplotlib.figure import Figure
    from matplotlib import pyplot as plt
    from io import BytesIO
    if len(dimension) != 3: 
        # error, this function only deals with 3D tensor
        return (-1, -1, -1)
    
    layer_size = dimension[0] * dimension[1]
    data_type = np.float32 if not is_float64 else np.float64
    
    with open(data_file, 'rb') as f:
        f.seek(layer_number * layer_size * np.dtype(data_type).itemsize)
        layer_data = np.fromfile(f, dtype=data_type, count=layer_size)
        layer_data = layer_data.reshape(dimension[0], dimension[1])

    fig = Figure(dpi=100)
    fig.subplots_adjust(bottom=0, top=1, left=0, right=1)
    ax = fig.add_subplot(111)
    data_min, data_max = np.min(layer_data), np.max(layer_data)
    ax.imshow(layer_data, cmap=plt.get_cmap('rainbow'), norm=plt.Normalize(vmin=data_min, vmax=data_max), aspect='auto')
    buf = BytesIO()
    fig.savefig(buf, format='png')

    return (buf, data_min, data_max)

def parse_arguments():
    """ Parse command-line arguments """
    parser = argparse.ArgumentParser(description="Run SZ3 Compression and Decompression with optional MPI support.")
    parser.add_argument("--use-mpi", action="store_true", help="Enable MPI for parallel execution.")
    parser.add_argument("--threads", type=int, default=2, help="Number of threads of processes depending on whether MPI is used.")
    return parser.parse_args()

def get_mpi_rank():
    return int(os.environ.get("OMPI_COMM_WORLD_RANK", 0))

def main():

    args = parse_arguments()
    use_mpi = args.use_mpi

    # comm = MPI.COMM_WORLD if use_mpi else None
    rank = 0 if not use_mpi else get_mpi_rank()
    size = args.threads
    # rank = comm.Get_rank() if use_mpi else 0
    # size = comm.Get_size() if use_mpi else args.threads

    if rank == 0:
        print(f"Running with MPI: {use_mpi}, Processes/Threads: {size}")
    
    DATA_PATH = Path("/media/briannas/Research/compression-data/nyx/")

    COMPRESSED_DATA_PATH = Path("/media/briannas/Research/compression-data/compressed/")

    temperature_raw_file_path = DATA_PATH / "temperature.f32"
    temperature_compressed_file_path = COMPRESSED_DATA_PATH / "temperature.f32.szsplit"
    temperature_decompressed_file_path = COMPRESSED_DATA_PATH  / "temperature.f32.dp"

    compress_args = [
        "compress",
        "-i", str(temperature_raw_file_path),
        "-o", str(temperature_compressed_file_path),
        "-d", "512", "512", "512",
        "-e", "0.01",
        "--mode", "layer",
        "--depth", "32",
        "--threads", str(args.threads)
    ]

    if use_mpi:
        compress_args.append("--mpi")

    # Compress the data
    # sz3py.compress(compress_args)

    decompress_args = [
        "decompress",
        "-i", str(temperature_compressed_file_path),
        "-o", str(temperature_decompressed_file_path),
        "-d", "512", "512", "512",
        "-e", "0.01",
        "--mode", "layer",
        "--depth", "32",
        "--threads", str(args.threads)
    ]

    if use_mpi:
        decompress_args.append("--mpi")

    # Decompress the data
    sz3py.decompress(decompress_args)

    sz3py.safeMPIFinalize()

    if rank == 0:
        # Visualize the raw data and the decompressed data
        layer_number = 400
        dimension = [512, 512, 512]
        is_float64 = False
        raw_image_data, raw_data_min, raw_data_max = get_partial_preview_data(dimension, temperature_raw_file_path, layer_number, is_float64)
        dp_image_data, dp_data_min, dp_data_max = get_partial_preview_data(dimension, temperature_decompressed_file_path, layer_number, is_float64)
        raw_image = Image.open(raw_image_data)
        raw_image_array = np.array(raw_image)
        dp_image = Image.open(dp_image_data)
        dp_image_array = np.array(dp_image)

        figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 5))
        axes[0].imshow(raw_image_array)
        axes[0].set_title(f"raw image at layer {layer_number}")

        axes[1].imshow(dp_image_array)
        axes[1].set_title(f"dp image at layer {layer_number}")
        print("Saving the figure of comparsion...")
        plt.tight_layout()
        plt.savefig("raw-dp-image.png", dpi=300)

if __name__ == '__main__':
    main()
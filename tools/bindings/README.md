# The Python binding of sz3_split

The package is named `sz3py`. To install the python package, go to the project's root directory. Note that it is recommended that you install the package into a Python virtual environment like a custom conda environment.

```bash
mkdir build && cd build
cmake ..
cmake --build . --target sz3py
cmake --install .
cd ..  # to the project root directory
pip install tools/bindings
```

To use the Python package, it is currently a direct mirror of the C++ version of the sz3_split command. We can use it this way. Note that the command options are exactly the same as the C++ version but we pass it in a Python way.

```c++
"Usage: sz3_split (de)compress [options]\n"
"options:  --threads/-t     INT   number of threads, default is 1, for MPI mode, this specifies the number of I/O processes\n"
"          --input/-i       STR   the RAW file/compressed file\n"
"          --output/-o      STR   the compressed file/decompressed file location\n"
"          --help/-h              print this help information\n"
"          --dimension/-d   STR   the data dimension of the file, e.g., 256 256 512\n"
"          --errorbound/-e  FLOAT the error bound to use in compression\n"
"          --float64              the default is float32 for each datapoint, this param changes it to float64\n"
"          --mode           STR   select 'layer' for layer-by-layer compression, 'direct' for direct compression\n"
"          --depth          INT   select the layer depth in layer-by-layer compression\n"
"          --mpi                  use MPI to run multiple processes";
```

```python
import sz3py
from pathlib import Path

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
    # "--mpi",
    # running a mpi program requires the `mpirun` prorgam instead of inside Python
    "--threads", "4"
]

sz3py.compress(compress_args)

decompress_args = [
    "decompress",
    "-i", str(temperature_compressed_file_path),
    "-o", str(temperature_decompressed_file_path),
    "-d", "512", "512", "512",
    "-e", "0.01",
    "--mode", "layer",
    "--depth", "32",
    # "--mpi",
    # running a mpi program requires the `mpirun` prorgam instead of inside Python
    "--threads", "4"
]

sz3py.decompress(decompress_args)
```

There are two example Python files under `${PROJECT_ROOT}/python-examples` folder. The [Jupyter Notebook](../../python-examples/sz3pytest.ipynb) is self-explanatory. You can run the multi-threaded sz3_split compression in there.

To run the MPI version of the sz3_split, we have to use the `mpirun` command instead of running a normal python file or jupyter notebook. You can check the provided [sz3py_mpi.py](../../python-examples/sz3py_mpi.py) file to see how it works.

```bash
mpirun -n 4 python sz3py_mpi.py --use-mpi --threads 2
```

You can also run the multi-threading version with the same python file but do not provide `--use-mpi` option. Note that the `--threads` option means the number of I/O threads when using `--use-mpi` option, while it means the total number of worker threads when not using `--use-mpi` option.

```bash
python sz3py_mpi.py --use-mpi --threads 4
```

We also provide a convenient tool to collect data for predicting compression performance. The following function runs the prediction method and the standard sz3 compression and decompression, so the collected data contains both data for prediction and actual compression performance.

```python
import sz3py
import pandas as pd

collected_data = sz3py.collect_prediction_data(
    inputFilePath=str(temperature_raw_file_path),
    dims=[512, 512, 512],
    error_bound=0.1
)

print(collected_data)
columns = collected_data[0]
data = collected_data[1:]
df = pd.DataFrame(data, columns=columns)
print(df)
```

SZ3.0: A Modular Compression Framework for Scientific Datasets

This version supports multiple error bounds in different data intervals or regions.
The idea can be found in the paper
[Optimizing Error-Bounded Lossy Compression for Scientific Data With Diverse Constraints](https://ieeexplore.ieee.org/document/9844293).

The parameters for the program `sz_region` is shown below. Both range and region based compression is implemented in the `sz_region` executable.

```text
-i            [input] input file name
-c            [output] the compressed file name
-q            [output] the decompressed file name
-r/--range    [param] specify error bounds for each range
-w/--region   [param] specify error bounds for each region
-d            [param] define the dimension of the data
-b            [param] use the background removal algorithm during compression
-m            [param] select the mode of compression, can be "test", "compress", "decompress", "bg_pre"
```

We provide some examples on how to use the program.

To run the range based compression, we use the `-r/--range` parameter. The format is
`range_start range_end error_bound;`, and you can add multiple ranges. The start of the next range must be the end of the previous range.
In the `IntegerQuantizer.hpp`, we specify a LARGE_ENOUGH_NUMBER to be 1,000,000, please make sure the data is within
-LARGE_ENOUGH_NUMBER to LARGE_ENOUGH_NUMBER, otherwise the program will render a incorrect result.

For instance,
we have the Katrina Hurricane simulation data, the dimension of which is (417642, 162). It has a background data -99999 which is very far from the normal data range.
The normal data range for this dataset is about from -2 to 10.

The compression program will still work if you don't handle the background data with a lower compression ratio. You can handle the background data by setting `-b -99999`,
and the program will use the background removal algorithm during compression. You can also set mode to be "bg_pre", and the
program will clean the data beforehand by setting all background data to 0 and then use the same algorithms to do compression.


If we hope to set error bound 0.2 for range [-1, 0.98), error bound 0.001 for range [0.98, 1.02),
and error bound 0.01 for range [1.02, 10), we can run the program this way.

```bash
sz_region -i katrina.bin -c katrina.bin.sz -q katrina.bin.out \
-r "-1 0.98 0.2; 0.98 1.02 0.001; 1.02 10 0.01;" -d "417642 162;" -m test
```

Similarly, for the QMCPACK dataset of dimension (33120, 69, 69), we can run the program this way.

```bash
sz_region -i QMCPACK.f32 -c QMCPACK.sz -q QMCPACK.dp \
-r "-17 -8 0.01; -8 -5 0.1; -5 17 1;" -d "33120 69 69;" -m test
```

To run the region based compression, we use the `-w/--region` parameter. The format is as follows:
`default_error_bound>x_lower y_lower z_lower:x_length y_length z_length error_bound;`

```bash
sz_region -i CLDHGH_1_1800_3600.dat -c CLDHGH.sz -q CLDHGH.dp \
--region "0.001>5 5 5: 8 10 12 0.1; 20 20 20: 90 90 90 0.01;" -d "1800 3600;" -l log.txt -m test
```

To clarify, the "20 20 20: 90 90 90 0.01;" means inside the bounding box 20<=x<=110, 20<=y<=110, 20<=z<=110,
the error bound is 0.01.


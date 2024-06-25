# Parallel-Marching-Squares-Algorithm

A parallel implementation of the Marching Squares algorithm using pthreads in C.

## Short Description

If the input image dimensions exceed 2048x2048, the image is scaled down using bicubic interpolation.
The image is divided into a grid, and each cell in the grid is assigned a binary value based on the average color intensity compared to a threshold (sigma).
The algorithm identifies the type of contour for each subgrid and updates the image with the corresponding contour pixels.
The workload is divided among multiple threads, each processing a portion of the image. Barriers are used to synchronize the threads at critical points.

Information about the marching squares algorithm can be found [here](https://www.baeldung.com/cs/marching-squares).

To run the program, use the following command: `gcc p_ms_algorithm.c helpers.c -o p_ms_algorithm -lm -lpthread -Wall -Wextra`

Execute using `./p_ms_algorithm <in_file> <out_file> <P>` where `in_file` is the input file name (ppm image) found in the `countours` directory, `out_file` is your custom output file name, and `P` is the number of threads to use.
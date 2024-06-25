
#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }

typedef struct {
    int thread_id;
    int thread_count;
    ppm_image *image;
    unsigned char **grid;
    ppm_image **contour_map;
    int step_x;
    int step_y;
    ppm_image *new_image;
    pthread_barrier_t *barrier;
} thread_data;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

// Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
// Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
// pixel values compare to the `sigma` reference value. The points are taken at equal distances
// in the original image, based on the `step_x` and `step_y` arguments.
void sample_grid(thread_data *data, unsigned char sigma, int start_index, int end_index, unsigned char **grid) {
    int p = data->image->x / data->step_x;
    int q = data->image->y/ data->step_y;
    start_index= start_index/data->step_y;
    end_index=end_index/data->step_y;


    for (int i = 0; i < p; i++) {
        for (int j = start_index; j < end_index; j++) {
            ppm_pixel curr_pixel = data->image->data[i * data->step_x * data->image->y + j * data->step_y];
            unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;
            if (curr_color > sigma) {
                grid[i][j] = 0;
            } else {
                grid[i][j] = 1;
            }
        }
    }
    grid[p][q] = 0;
    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    if(data->thread_id==data->thread_count-1){
    for (int i =  0; i < p; i++) {
        ppm_pixel curr_pixel = data->image->data[i * data->step_x * data->image->y + data->image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[i][q] = 0;
        } else {
            grid[i][q] = 1;
        }
    }
    for (int j = 0; j < q; j++) {
        ppm_pixel curr_pixel = data->image->data[(data->image->x - 1) * data->image->y + j * data->step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > sigma) {
            grid[p][j] = 0;
        } else {
            grid[p][j] = 1;
        }
    }
}
}

// Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
// type of contour which corresponds to each subgrid. It determines the binary value of each
// sample fragment of the original image and replaces the pixels in the original image with
// the pixels of the corresponding contour image accordingly.
void march( thread_data *data, int start, int end) {
    int p = data->image->x / data->step_x;
    int q = data->image->y/ data->step_y;
    start= start/data->step_y;
    end=end/data->step_y;
    

    for (int i = 0; i < p; i++) {
        for (int j = start; j < end; j++) {
            unsigned char k = 8 * data->grid[i][j] + 4 * data->grid[i][j + 1] + 2 *data-> grid[i + 1][j + 1] + 1 * data->grid[i + 1][j];
            update_image(data->image, data->contour_map[k], i * data->step_x, j * data->step_y);
        }
    }
}


// Calls `free` method on the utilized resources.
void free_resources(thread_data *data) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        free(data->contour_map[i]->data);
        free(data->contour_map[i]);
    }
    free(data->contour_map);

    for (int i = 0; i <= data->image->y/ data->step_y; i++) {
        free(data->grid[i]);
    }
    free(data->grid);

    free(data->image->data);
    free(data->image);

}

void rescale_image(thread_data *data, int start, int end) {
    uint8_t sample[3];

    // use bicubic interpolation for scaling
    for (int i = 0; i < data->new_image->x; i++) {
        for (int j = start; j < end; j++) {
            float u = (float)i / (float)(data->new_image->x - 1);
            float v = (float)j / (float)(data->new_image->y - 1);
            sample_bicubic(data->image, u, v, sample);

            data->new_image->data[i * data->new_image->y + j].red = sample[0];
            data->new_image->data[i * data->new_image->y + j].green = sample[1];
            data->new_image->data[i * data->new_image->y + j].blue = sample[2];
        }
    }

}

void *thread_func(void *arg) {
    thread_data *data = (thread_data *)arg;

    int n = data->image->y;
    int start_index = data->thread_id * (double) data->image->y / data->thread_count;
    int end_index =fmin((data->thread_id + 1) * (double) data->image->y/ data->thread_count, data->image->y);
    if (data->thread_id == data->thread_count - 1) {
        end_index = n;
    }

    if (data->image->x >RESCALE_X || data->new_image->y > RESCALE_Y) {
        start_index = data->thread_id * (double) data->new_image->y / data->thread_count;
        end_index =fmin((data->thread_id + 1) * (double) data->new_image->y/ data->thread_count, data->new_image->y);
        rescale_image(data, start_index, end_index);
        pthread_barrier_wait(data->barrier);
        data->image=data->new_image;
    }
        sample_grid(data, SIGMA, start_index, end_index, data->grid);
        pthread_barrier_wait(data->barrier);
        march(data, start_index, end_index);
        pthread_barrier_wait(data->barrier);
    
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./p_ms_algorithm <in_file> <out_file> <P>\n");
        return 1;
    }

    ppm_image *image = read_ppm(argv[1]);
    ppm_image *new_image = (ppm_image *)malloc(sizeof(ppm_image));

    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    int THREAD_COUNT = atoi(argv[3]);

    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        new_image = image;
    }
    else{
    new_image->x = RESCALE_X;
    new_image->y = RESCALE_Y;
    
    new_image->data = (ppm_pixel*)malloc(new_image->x * new_image->y * sizeof(ppm_pixel));
    if (!new_image) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
        }
    }

    ppm_image **contour_map = init_contour_map();

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, THREAD_COUNT);

    thread_data* data;
    data=(thread_data*)malloc(THREAD_COUNT*sizeof(thread_data));
    data->contour_map = contour_map;
    data->image = image;
    data->new_image = new_image;
    data->step_x = STEP;
    data->step_y = STEP;
    data->barrier = &barrier;

    int p,q;
    p = new_image->x / data->step_x;
    q = new_image->y / data->step_y;
    
    pthread_t threads[THREAD_COUNT];
    unsigned char **grid = (unsigned char **)malloc((p + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= p; i++) {
        grid[i] = (unsigned char *)malloc((q + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }
    int r;
    for (int i = 0; i < THREAD_COUNT; ++i) {
        data[i].grid=grid;
        data[i].new_image=new_image;
        data[i].thread_id = i;
        data[i].thread_count = THREAD_COUNT;
        data[i].image = image;
        data[i].contour_map = contour_map;
        data[i].step_x = STEP;
        data[i].step_y = STEP;
        data[i].barrier = &barrier;
        r= pthread_create(&threads[i], NULL, thread_func, &data[i]);
     if(r){
        printf("Error joining thread %d\n", i);
        exit(-1);
    }
    }
    for (int i = 0; i < THREAD_COUNT; ++i) {
        if (pthread_join(threads[i], NULL)) {
            fprintf(stderr, "Error joining thread\n");
            exit(-1);
        }
    }

    pthread_barrier_destroy(&barrier);

    write_ppm(data->image, argv[2]);
    free_resources(data);
    free(data); 
    return 0;
}

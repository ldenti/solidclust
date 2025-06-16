#ifndef OPTION_H
#define OPTION_H

#define TMP_FILE_SIZE 20

typedef struct {
    char *input_fastq;
    char *output_mapping;
    char *tmp_filename;
    unsigned char k;
    unsigned char w;
    unsigned char canonical;
    unsigned char weighted;
    double quality_threshold;
    double similarity_threshold;
    double post_cluster;
    unsigned long long seed;
    unsigned char tmp_filename_allocated;
} option_t;

void init_options(option_t *const opts);
void destroy_options(option_t *const opts);

#endif /* OPTION_H */

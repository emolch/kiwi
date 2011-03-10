


void gfdb_io_init_(bool *ok) {
}


void gfdb_io_deinit_() {
}


void gfdb_io_read_index_(char *filename, float *dt, float *dx, float *dz, 
                         float *firstx, float *firstz, int *nchunks, int *nx, 
                         int *nxc, int *nz, int *ng, bool *ok ) {
}


void gfdb_io_create_index_(char *filename, float *dt, float *dx, float *dz, 
                          float *firstx, float *firstz, int *nchunks, int *nx, 
                          int *nxc, int *nz, int *ng, bool *ok ) {
}


void gfdb_io_create_chunk_(char *filename, int *ng, int *nz, int *nxc, 
                          int *size_hint, bool *ok) {
}


void gfdb_io_save_trace(char *filename, FILE *file, int *dataset_index, 
            int *ixc, int *iz, int *ig, 
            float *packed, int *packed_size, int *pofs, int *ofs, int *nstrips, 
            int *size_hint1, int *size_hint2, 
            void *reference, bool *ok ) {
}


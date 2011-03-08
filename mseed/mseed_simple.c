
#include <math.h>
#include <stdio.h>
#include <libmseed.h>

static void record_handler (char *record, int reclen, void *outfile) {
    if ( fwrite(record, reclen, 1, outfile) != 1 ) {
      fprintf(stderr, "Error writing mseed record to output file\n");
    }
}

int writemseed(char *filename, float *seismogram, int length, double toffset, double deltat, 
    char *network, char *station, char *location, char *channel  ) {
    MSTrace     *trace;
    char        srcname[50];
    int         psamples, precords;
    FILE        *outfile;
    double      sec_em5;
    hptime_t    hptime;
    
    outfile = fopen(filename, "w" );
    if (outfile == NULL) {
        fprintf (stderr, "Error opening file %s\n", filename);
        return -1;
    }
    
    trace = mst_init (NULL);

    strncpy( trace->network, network, 2 );
    strncpy( trace->station, station, 4 );
    strncpy( trace->location, location, 2 );
    strncpy( trace->channel, channel, 3 );
    trace->network[3]  = '\0';
    trace->station[5]  = '\0';
    trace->location[3] = '\0';
    trace->channel[4]  = '\0';


    /* double is only precise up to something slightly below 1e-5 s in this 
       century, so last digit in  hptime with HPTMODULUS==1000000 is not reliable.
       I try to set it to zero here... */
    
    sec_em5 = rint(toffset*1e5)*1e-5;
    hptime = (hptime_t)(sec_em5 * HPTMODULUS);
        
    trace->starttime = hptime;
    trace->samprate = 1./deltat;
    trace->sampletype = 'f';
    
    /* The datasamples pointer and numsamples counter will be adjusted by
       the packing routine, the datasamples array must be dynamic memory
       allocated by the malloc() family of routines. */
    
    trace->datasamples = seismogram;
    trace->numsamples = length;
    trace->samplecnt = length;
    trace->sampletype = 'f';
    mst_srcname (trace, srcname, 0);
    
    precords = mst_pack (trace, &record_handler, outfile, 4096, DE_FLOAT32,
                                     1, &psamples, 1, 0, NULL);

    mst_free( &trace );                                 
    fclose( outfile );
    return 0;
}


int readmseed(char *filename, float **seismogram, int *length,  double *toffset, double *deltat) {

    MSTraceGroup *mstg = NULL;
    MSTrace      *trace;
    int retcode;
    
    char *source_char;
    int *source_int;
    float *source_float;
    double *source_double;
        
    float *data;
    
    double dtime;
    hptime_t hptime;
    
    int i;
        
    retcode = ms_readtraces (&mstg, filename, 0, -1.0, -1.0, 0, 1, 1, 0);
    if ( retcode < 0 )
        fprintf (stderr, "Cannot read %s: %s\n", filename, ms_errorstr(retcode));

    if ( ! mstg ) {
        fprintf (stderr, "Error reading file\n");
        return -1;
    }
    
    trace = mstg->traces;
    
    if (trace->datasamples == NULL) {
        fprintf (stderr, "Error reading file\n");
        return -1;
    }

    data = (float*) calloc(sizeof(float), trace->numsamples);
    
    if (trace->sampletype == 'i') {
        source_int = (int*)trace->datasamples;
        for (i=0; i<trace->numsamples; i++) {
            data[i] = (float)source_int[i];
        }
    }
    if (trace->sampletype == 'a') {
        source_char = (char*)trace->datasamples;
        for (i=0; i<trace->numsamples; i++) {
            data[i] = (float)source_char[i];
        }
    }
    if (trace->sampletype == 'f') {
        source_float = (float*)trace->datasamples;
        for (i=0; i<trace->numsamples; i++) {
            data[i] = source_float[i];
        }
    }
    if (trace->sampletype == 'd') {
        source_double = (double*)trace->datasamples;
        for (i=0; i<trace->numsamples; i++) {
            data[i] = (float)source_double[i];
        }
    }
    
    hptime = trace->starttime ;
    
    /* loosing last digit of precision here... */
    dtime = (double)hptime / ((double)HPTMODULUS);

    *length = trace->numsamples;
    *toffset = dtime;
    *deltat = 1./trace->samprate;
    *seismogram = data;
    
    mst_freegroup (&mstg);

    return 0;
}

/*  Fortran wrappers

    Tested with ifort and g95, which follow the convention, that length of fortran strings are
    appended to the end of the argument list.

    readmseed1_ reads seismogram into *databuffer and returns length and other infos
    readmseed2_ copies the seismogram to a a given fortran array
    
    This way a fortran array of the required length can be allocated between the calls. 

    Example usage:
    
    program mseedio
        
        implicit none
        
        integer                         :: length
        character(len=1024)             :: filename
        real, dimension(:), allocatable :: seismogram
        real(kind=8)                    :: toffset, deltat
        integer                         :: nerr
        integer                         :: i
        
      ! read the file
        filename = 'test.mseed';
        call readmseed1(trim(filename)//char(0), length, toffset, deltat, nerr )
        if (nerr .ne. 0) stop 'failed to read file test.mseed'
    
      ! make seimsmogram large enough and fill in data
        allocate( seismogram(length) )
        call readmseed2( seismogram )
        
      ! print seismogram to stdout
        do i=1, length
            print *, toffset+(i-1)*deltat, seismogram(i)
        end do
        
      ! output seismogram to a different file
        filename = 'abcd.mseed'
        call writemseed(trim(filename)//char(0), seismogram, length, toffset, deltat, nerr )
        if (nerr .ne. 0) stop 'failed to write file abcd.mseed'
    
    end program

*/
    
float *databuffer;
int databufferlength;



/* readmseed1_()
 
   *** not thread safe ***

   Fortran wrapper, to read a single mseed trace from file.
   It is not directly returned to the caller, but it's length, so that
   the fortran has a chance to allocate memory for the trace.
   The actual data can then be aquired with a call to readmseed2_().
   
   IN:
     filename:  Name of file, must be null-terminated externally.
   
   OUT;
     length:    Number of samples in the trace.
     toffset:   Time of first sample. Seconds since epoch.
                Double precision gives accuracy of up to 1e-5 s.
     deltat:    Sampling interval.
     nerr:      Zero if everything went well.
*/
     
void readmseed1_(char *filename, int *length, double *toffset, double *deltat, int *nerr ) {

    databuffer = NULL;
    databufferlength = 0;
    *nerr = readmseed(filename, &databuffer, length, toffset, deltat);
    databufferlength = *length;
}

/* readmseed2_()

   Get data read by previous call to readmseed1_().
   
   OUT;
     seismogram: Array long enough to hold the trace data. (See readmseed1_())
*/

void readmseed2_( float *seismogram ) {
    int i;
    if (databuffer != NULL) {
        for (i=0;i<databufferlength; i++) {
            seismogram[i] = databuffer[i];
        }
        free( databuffer );
    }
}


/* writemseed_()

   Dump a seismogram trace to file.
   
   IN:
     filename:    Name of file, must be null-terminated externally.
     seismogram:  Seismogram data.
     length:      Number of samples in seismogram.
     toffset:     Time of first sample. Seconds since epoch.
                  Double precision gives accuracy of up to 1e-5 s.
     deltat:      Sampling interval.
   
   OUT;
     nerr:        Zero if all went nice.
       
*/
   
void writemseed_( char *filename, float *seismogram, int *length, double *toffset,
                  double *deltat, 
                char *network, char *station, char *location, char *channel, int *nerr ) {

    float *data;
    int i;
    
    data = calloc(*length,sizeof(float));
    for (i=0;i<(*length); i++) {
        data[i] = seismogram[i];
    }
    *nerr = writemseed( filename, data, *length,  *toffset, *deltat, network, station, location, channel ); /* data is deallocated inside of this */
}


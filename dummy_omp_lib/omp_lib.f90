module omp_lib

    implicit none
    
    private
    public :: omp_get_num_threads
    public :: omp_get_thread_num
    public :: omp_get_max_threads

  contains
    
    pure function omp_get_num_threads()
        integer :: omp_get_num_threads 
        omp_get_num_threads = 1        
    end function

    pure function omp_get_max_threads()
        integer :: omp_get_max_threads 
        omp_get_max_threads = 1        
    end function

    pure function omp_get_thread_num()
        integer :: omp_get_thread_num
        omp_get_thread_num = 0        
    end function

end module

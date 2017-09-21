      integer fakesize, fakerank
      parameter (fakesize = 1)
      parameter (fakerank = 1)
c
c     Other variables used by MPI.  These are best treated as
c      integers:
c
      integer MPI_COMM_WORLD,MPI_DOUBLE_PRECISION,MPI_INTEGER,
     $     MPI_LOGICAL
c
c     These are actually MPI operations in MPI_REDUCE and 
c      MPI_ALLREDUCE, but we might as well treate them as integers
c      here:
c
      integer MPI_MIN,MPI_MAX,MPI_SUM,MPI_PROD

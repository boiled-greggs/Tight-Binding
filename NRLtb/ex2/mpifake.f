c
c     This is a collection of fake MPI subroutines which allow us
c      to run the code on a single processor without too much
c      conversion.  Note that this way is preferable when we
c      don't have a consistent set of compiler directives from
c      machine to machine.
c
c=======================================================================
c
      subroutine MPI_INIT (ierr)
c
c     Should be OK:
c
      ierr = 0
      return
      end
c
c=======================================================================
c
      subroutine MPI_COMM_SIZE (MPI_COMM_WORLD, MPISIZE, ierr)
c
c     It's a single processor:
c
      MPI_COMM_WORLD = 0
      MPISIZE = 1
      ierr = 0
      return
      end
c
c=======================================================================
c
      subroutine MPI_COMM_RANK (MPI_COMM_WORLD, MPIRANK, ierr)
c
c     MPI_COMM_WORLD should already be defined.  This is a single
c      processor machine:
c
      MPIRANK = 0
      ierr = 0
      return
      end
c
c=======================================================================
c
      subroutine MPI_BCAST(x,istart,itype,iplace,MPI_COMM_WORLD,ierr)
c
c     We don't want to broadcast anything here, especially since
c      we don't know if x is really a double precision variable:
c
      ierr = 0
      return
      end
c
c=======================================================================
c
      subroutine MPI_FINALIZE(ierr)
      ierr = 0
      return
      end
c
c=======================================================================
c
      subroutine MPI_GATHERV(array,icount,itype,aarray,jcount,ioff,
     $     jtype,jplace,MPI_COMM_WORLD,ierr)
c
c     Don't even think about doing anything:
c
      ierr = 0
      return
      end
c
c=======================================================================
c
      subroutine MPI_ALLREDUCE(xlocal,xfinal,idim,mptype,mpop,
     $     MPI_COMM_WORLD,ierr)
c
c     In an ideal world, this would copy the idim elements of xlocal
c      into xfinal, but we don't know what xlocal and xfinal
c      really are, so we can't do much about that:
c
      ierr = 0
      return
      end
c=======================================================================
c
      subroutine MPI_REDUCE(xlocal,xfinal,idim,mptype,mpop,inode,
     $     MPI_COMM_WORLD,ierr)
c
      ierr = 0
      return
      end

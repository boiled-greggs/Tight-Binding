      program bandplt
c
c     Given the band structure produced by the LAPW program in the
c      Binary QLMT file format, produces the output in a file which can
c      be read by plotting programs such as gnuplot to produce a band
c      structure plot
c
c     Note that our definition of reciprocal lattice vectors is
c      without the 2 pi factor.
c
c
      implicit double precision (a-h,o-z)
      double precision kold(3),k(3)
c
c     If you want to plot more than 1000 bands at a given k-point
c      you are going to have to redimension this.
c
      dimension bande(2000),avec(3,3),bvec(3,3),vold(3),vnew(3)
      dimension kp_directions(20),ene_kp(20)
      character*25  infile,outfile,outfileg
      character*20  title
	real*8 knum
c=======================================================================
c
      parameter (zero = 0d0)
      parameter (shift = 0.2854d0)
	parameter (twopi = 6.28318530717958647693d0)
c
c=======================================================================
c
c     Input and output file information:
c
        open(9,file='dosdat.in')

c         read(9,*) ymin,dummy,ymax
         read(9,1050) ymin,dummy,ymax
1050	format(/////,3f10.5)
        open(10,file='bandplot.in')
c100  write(*,'('' Name of input file?'')')
      read(10,1051,end=1000) title
1051	format(A20)
      read(10,*) structure
c        structure: (1=bcc, 2=fcc, 3=sc,...)
      read(10,'(a)',end=1000) infile
      open(unit=11,file=infile,status='unknown')
c     write(*,'('' Enter direct space primitive lattice vectors:  '')')
      read(10,*) 
      read(10,*) ((avec(i,j),i=1,3),j=1,3)
c
c     Find reciprocal space vectors:
c
c     First just get them orthogonal to the avec:
c
      do j=1,3
         j1=j+1
         if(j1.gt.3) j1=j1-3
         j2=j1+1
         if(j2.gt.3) j2=j2-3
         do i=1,3
            i1=i+1
            if(i1.gt.3) i1=i1-3
            i2=i1+1
            if(i2.gt.3) i2=i2-3
            bvec(i,j)=avec(i1,j1)*avec(i2,j2)-avec(i2,j1)*avec(i1,j2)
         end do
      end do
c
c     Compute the unit cell volume
c
      vol=zero
      do i=1,3
         vol=vol+avec(i,1)*bvec(i,1)
      end do
c
c     Now get the vectors properly scaled.  Note that
c
c     Sum[ avec(i,j)*bvec(i,k), i = 1,3] = delta(j,k)
c
      do i=1,3
         do j=1,3
            bvec(i,j)=bvec(i,j)/vol
         end do
      end do
c
      write(*,1225) ((avec(i,j),i=1,3),j=1,3),vol,
     $     ((bvec(i,j),i=1,3),j=1,3)
 1225 format(/' Primitive vectors of the direct lattice:'/
     $     3(/1X,3F12.6)//' Volume of the unit cell=',F18.8/
     $     /' Primitive vectors of the reciprocal lattice:'/
     $     3(/1X,3F12.6))
c
c     write(*,1231)
 1231 format(' Enter number of bands to be output, and maximum'/
     $     ' eigenvalue to be printed, which is also the number that'/
     $     ' will be printed if there are not enough bands:')
c
      read(10,*) 
      read(10,*) nbout,eigdum
      read(10,*) 
      read(10,*)efermi
      write(*,*)efermi
c
c -----------------------------------------------------------------
c        read the number of points in each k-point direction
c        write(*,*)' Enter the number of k-directions in the Brillouin zone '
        read(10,*)
        read(10,*)n_dirkp
        if(n_dirkp.gt.20)then
                stop 'Please increase the dimension of directions_kp(20) 
     &          in bandplot.f'
        else
                ksum=0.
                do i=1,n_dirkp
c                write(*,*)'Enter the number of k-point in the ',i,' direction'
                read(10,*)kp_tot
                 ksum=ksum+kp_tot
                 kp_directions(i)=ksum
                enddo
        endif
        write(*,*)'Number of kp-directions = ',n_dirkp
        write(*,*)'Finall k-point in each direction  = ',
     &             (kp_directions(i),i=1,n_dirkp)
c -----------------------------------------------------------------
c     write(*,'('' Name of output file?  '')')
c
      read(10,'(a)') outfile
      open(unit=21,file=outfile,status='unknown')
      read(10,'(a)') outfileg
      open(unit=31,file=outfileg,status='unknown')
c
c        read the ymin and ymax of the band structure plot
c       read(10,*)ymin,ymax
c
c     Read in the number of windows to be plotted
c      and the number of atoms in the QLMT file.  Note that "0" is
c      an allowed value for natoms
c
      read(11,*) nspins,nwind,nwsc,natoms
      write(6,145) nspins,nwind,nwsc,natoms
 145  format(/1x,i5,' Spins',i5,' Windows',i5,' Semi-core windows'/
     $     ' The system contains ',i5,' atoms')
c
      if(nspins.ne.1) stop 'Only Spin=1 allowed for now'
c
c     Plot each window separately:
c
      do 400 nw=1,nwind
c
c        How many points in this window? (Plus extraneous stuff)
c
         read(11,*) npts
c        read(11,*) npts,nat1,zfill,tkb
c
         write(*,165) npts,nw
 165     format(//1x,i5,' points in window number ',i5/)
c
        if(npts.lt.1) go to 400
c
c       We'll plot the information on a linear graph.  The y coordinates
c         will be the band energies.  The x coordinate will be the sum of
c         the distances between the neighboring points in k space.  This
c         means that the first point has distance zero.  It also means
c         that we'll have to handle the first point separately.
c         The information is read in by a subroutine, and the
c         band energies are stored in BANDE.  KNUM is the number of points
c         found in this band at this K.
c
        dis=zero
        call rdband(kold,knum,bande,natoms,nbout,eigdum)
c
c       Convert to Cartesian coordinates:
c
        do l=1,3
           vold(l)=zero
           do j=1,3
              vold(l)=vold(l)+kold(j)*bvec(l,j)
           end do
        end do
c
c$$$        WRITE(6,175) KNUM,KOLD,VOLD,DIS,(BANDE(I),I=1,KNUM)
c$$$175     FORMAT(/1X,I5,' Points at K = (',3F10.4,' ) = ('
c$$$     1     3F10.4,')'/1X,7F10.5/
c$$$     2         (10X,6F10.5))
c
c       Write NBOUT eigenvalues, including dummy values, if any:
c
c       do ii=1,nbout-1
c          do jj=1,nbout-ii
c             if(bande(jj+1).lt.bande(jj)) then
c                temp=bande(jj)
c                bande(jj)=bande(jj+1)
c                bande(jj+1)=temp
c             end if
c          end do
c       end do

c	Modified by lwn 01/09 to adjust TB x-axis to LAPW x-axis for plot
c         write(21,185) dis,(bande(i)-efermi,i=1,nbout)
c shift efermi by a known constant 
         write(21,185) dis,(bande(i),i=1,nbout)
c       write(21,185) dis*1.7071/0.224,(bande(i)-efermi,i=1,nbout)
185     format(120f13.6)
         write(*,*) 'Gamma  point at = ',0.0
c

        do i=2,npts
           call rdband(k,knum,bande,natoms,nbout,eigdum)
c
c          Convert to Cartesian:
c
           do j=1,3
              vnew(j)=0d0
              do l=1,3
                 vnew(j)=vnew(j)+k(l)*bvec(j,l)
              end do
           end do
c
c          Find distance between this point and the last and
c           add it to DIS:
c
           dsum=zero
           do j=1,3
              dsum=dsum+(vnew(j)-vold(j))**2
           end do
           dis=dis+sqrt(dsum)
c       do ii=1,nbout-1
c          do jj=1,nbout-ii
c             if(bande(jj+1).lt.bande(jj)) then
c                temp=bande(jj)
c                bande(jj)=bande(jj+1)
c                bande(jj+1)=temp
c             end if
c          end do
c       end do

c	Modified by lwn 01/09 to adjust TB x-axis to LAPW x-axis for plot
c           write(21,185) dis,(bande(l)-efermi,l=1,nbout)
           write(21,185) dis,(bande(l),l=1,nbout)
c         write(21,185) dis*1.7071/0.224,(bande(l)-efermi,l=1,nbout)
             do ikkk=1,n_dirkp
                if(kp_directions(ikkk).eq.i)then
                ene_kp(ikkk)=dis
                write(*,*)i,' High symmetry points at = ',ene_kp(ikkk)
                endif
             enddo
c
c         Update the "last" vector:
c
          do j=1,3
             vold(j)=vnew(j)
             kold(j)=k(j)
          end do
       end do
 400  continue
 1000 continue
c  ------------------------ Make the gnuplot file ---------------------
c
        write(31,*)'set title "Band structure of  ',title,'"'
c
        write(31,*)' set ylabel "Energy (Ry)" '

c
c
  176 format(A10,I6,A6,F7.3,A13,F7.3,A20)
c       if(structure.eq.2) then
c          ene_kp(5)=ene_kp(4)+(ene_kp(5)-ene_kp(4))*3./4.
c       end if

c       if(structure.ne.2) then
           do ii=1,n_dirkp
              write(31,176)'set arrow', ii,' from ',ene_kp(ii),
     $        ', graph 0 to ',ene_kp(ii),', graph 1 nohead '
           enddo
c       else
c          do ii=1,n_dirkp-1
c             write(31,176)'set arrow', ii,' from ',ene_kp(ii) ,
c    $        ', graph 0 to ',ene_kp(ii),', graph 1 nohead '
c          enddo
c          ii=n_dirkp-1
c          write(31,176)'set arrow', ii+1,' from ',(ene_kp(ii+1)-
c    $     ene_kp(ii))*6.*sqrt(2.)/(sqrt(40.)+6.*sqrt(2.))+ene_kp(ii),
c    $     ', graph 0 to ',(ene_kp(ii+1)-ene_kp(ii))*6.*sqrt(2.)/(
c    $     sqrt(40.)+6.*sqrt(2.))+ene_kp(ii),', graph 1 nohead '
c          write(31,176)'set arrow', ii+2,' from ',ene_kp(ii+1) ,
c    $     ', graph 0 to ',ene_kp(ii+1),', graph 1 nohead '

c          write(31,176)'set arrow', ii+1,' from ',(ene_kp(ii+1)-
c    $     ene_kp(ii))*73./101.+ene_kp(ii),
c    $     ', graph 0 to ',(ene_kp(ii+1)-ene_kp(ii))*73./
c    $     101.+ene_kp(ii),', graph 1 nohead '
c          write(31,176)'set arrow', ii+2,' from ',ene_kp(ii+1) ,
c    $     ', graph 0 to ',ene_kp(ii+1),', graph 1 nohead '
c       end if

       write(31,188)'set encoding iso_8859_1'
c
c                     bcc --------------------------------------------
        if(structure.eq.1)then
c 188 format(A10,A20,A4,1X,F7.3,A5,1X,F7.3,A15,1X,F7.3,A5,1X,F7.3,
c    $     A5,1X,F7.3,a5,1X,f7.3,a5,1X,f7.3,a3)
  188 format(A10,A20,14(A15,1X,F7.3),A3)
       write(31,188)'set xtics','("{/Symbol G}" 0.,',
     $'"{/Symbol D}"',ene_kp(1)/2., 
     $',"H"',ene_kp(1), 
     $',"G"',ene_kp(1)+(ene_kp(2)-ene_kp(1))/2.,
     $',"N"', ene_kp(2), 
     $',"{/Symbol S}"',ene_kp(2)+(ene_kp(3)-ene_kp(2))/2.,
     $',"{/Symbol G}"',ene_kp(3),
     $',"{/Symbol L}"',ene_kp(3)+(ene_kp(4)-ene_kp(3))/2.,
     $',"P"', ene_kp(4),
     $',"D"',ene_kp(4)+(ene_kp(5)-ene_kp(4))/2., 
     $',"N"',ene_kp(5),
     $',"D"',ene_kp(4)+(ene_kp(6)-ene_kp(5))/2.,
     $',"P"', ene_kp(6),
     $',"F"',ene_kp(6)+(ene_kp(7)-ene_kp(6))/2.,
     $',"H"',ene_kp(7),
     $')'
        endif
c                     fcc --------------------------------------------
        if(structure.eq.2)then
c 189 format(A10,A20,A4,1X,F7.3,A5,1X,F7.3,A5,1X,F7.3,
c    $A20,1X,F7.3,A5,1X,F7.3,a3)
  189 format(A10,A20,14(A15,1X,F7.3),A3)
       write(31,189)'set xtics','("{/Symbol G}" 0.,',
     $'"{/Symbol D}"',(ene_kp(1)/2.),
     $',"X"',ene_kp(1),
     $',"Z"',(ene_kp(1)+(ene_kp(2)-ene_kp(1))/2.), 
     $',"W"', ene_kp(2),
     $',"Q"',(ene_kp(2)+(ene_kp(3)-ene_kp(2))/2.), 
     $',"L"', ene_kp(3), 
     $',"{/Symbol L}"',(ene_kp(3)+(ene_kp(4)-ene_kp(3))/2.),
     $',"{/Symbol G}"',ene_kp(4),
     $',"{/Symbol S}"',(ene_kp(4)+(ene_kp(5)-ene_kp(4))/2.),
c    $',"K"',ene_kp(4)+(ene_kp(5)-ene_kp(4))*6.*sqrt(2.)/
c    $ (8.+6.*sqrt(2.)),
c    $',"X"', ene_kp(5), 
c    $',"K"', ene_kp(4)+(ene_kp(5)-ene_kp(4))*3./4., 
     $',"K"', ene_kp(5), 
     $')'
        endif
c
c                     sc --------------------------------------------
        if(structure.eq.3)then
c 187 format(A10,A20,A4,1X,F7.3,A5,1X,F7.3,A15,1X,F7.3,A5,1X,F7.3,
c    $     A5,1X,F7.3,a5,1X,f7.3,a5,1X,f7.3,a3)
  187 format(A10,A20,14(A15,1X,F7.3),A3)
       write(31,187)'set xtics','("{/Symbol G}" 0.,',
     $'"{/Symbol D}"',ene_kp(1)/2.,
     $',"X"',ene_kp(1), 
     $',"Z"',ene_kp(1)+(ene_kp(2)-ene_kp(1))/2.,
     $',"M"', ene_kp(2), 
     $',"{/Symbol S}"',ene_kp(2)+(ene_kp(3)-ene_kp(2))/2.,
     $',"{/Symbol G}"',ene_kp(3),
     $',"{/Symbol L}"',ene_kp(3)+(ene_kp(4)-ene_kp(3))/2.,
     $',"R"', ene_kp(4), 
     $',"S"',ene_kp(4)+(ene_kp(5)-ene_kp(4))/2.,
     $',"X"', ene_kp(5),
     $',"S"',ene_kp(5)+(ene_kp(6)-ene_kp(5))/2., 
     $',"R"', ene_kp(6),
     $',"T"',ene_kp(6)+(ene_kp(7)-ene_kp(6))/2., 
     $',"M"', ene_kp(7), 
     $')'
        endif

c
c        write(31,*)'efermi =',efermi-efermi
        write(31,*)'efermi =',efermi
  177 format(A10,I6,A6,F5.3,A3,F10.5,A4,F7.3,A3,F10.5,A9)
        write(31,177)'set arrow ', n_dirkp+1 ,' from ',0.0,','
     &  ,efermi, ' to ', ene_kp(5),',',efermi,'  nohead'
c     &  ,efermi, ' to ', ene_kp(n_dirkp),',',efermi ,'  nohead'
c
        write(31,*)'set yrange [',ymin,':',ymax,']'
        write(31,*)'set xrange [ 0.0    :',ene_kp(5),']'
c
        write(31,*)'set term postscript eps enhanced solid'
        write(31,*)'set output "band.eps" '
c       write(31,*)'replot'

        write(31,"(a)",advance="NO") ' plot "banddata" u 1:2 t "" w l, '
        do ib = 3, nbout
          ib1 = mod(ib,10)
          ib2 = int((ib-ib1)/10)
          if(ib2.eq.0) then
            write(31,"(a)",advance="NO")            
     $    ' "banddata" u 1:'//char(48+ib1)//' t "" w l, '
          else
          write(31,"(a)",advance="NO")
     $    ' "banddata" u 1:'//char(48+ib2)//char(48+ib1)//' t "" w l, '
          endif
        end do
        ib1 = mod(nbout+1,10)
        ib2 = int((nbout+1-ib1)/10)
        write(31,"(a)",advance="NO")
     $    ' "banddata" u 1:'//char(48+ib2)//char(48+ib1)//' t "" w l '


        write(31,*)'set output '
c       write(31,*)'set term x11 '
        

c
c  ------------------------  gnuplot file ---------------------
      end
c_______________________________________________________________________
c
      subroutine rdband(k,knum,bande,natoms,nbout,eigdum)
      implicit real*8 (a-h,o-z)
      dimension bande(*)
      real*8 k(3),dummy(4),knum
c
c     Reads the information about one K point in a window from
c       the binary QLMT file in file descriptor 11
c
c     K is the point in momentum space in reciprocal lattice coordinates.
c     KNUM is the number of energy bands found in this window at this
c       point.  Everything else is a dummy variable
c
      read(11,*) k,wt0,knum,wght
ctemp
c     write(*,*) k,knum
cend temp
c
c     Now read the energies for this band.  There is a lot of information
c       which is not relevant to this plot, and this information is read
c       into dummy:
c
      do i=1,int(knum)
         read(11,*) bande(i)
ctemp
cend temp
c
c        There is one line of unneeded information per "atom".
c        Note that natoms=0 skips this entirely.
c
         do j = 1,natoms
            read(11,*) dummy
         end do
      end do
       do ii=1,int(knum)-1
          do jj=1,int(knum)-ii
             if(bande(jj+1).lt.bande(jj)) then
                temp=bande(jj)
                bande(jj)=bande(jj+1)
                bande(jj+1)=temp
             end if
          end do
       end do

c
c     If there are not nbout eigenstates, fill in a dummy value:
c
      if(knum.lt.nbout) then
         do i=int(knum)+1,nbout
            bande(i)=eigdum
         end do
      end if
      
c
c     If an eigenvalue is > eigdum, set it equal to eigdum:
c
      do i = 1,nbout
         bande(i) = min (bande(i),eigdum)
      end do
      return
      end

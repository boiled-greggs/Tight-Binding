      subroutine setlat(alat,vol,ierr)
c
c     Given the lattice type and a set of lattice parameters, this
c      routine computes the primitive vectors (alat) and
c      volume (vol) of the desired primitive cell.  The lattice
c      parameters may include:
c
c     (1) A choice of lattice types (see below) specified by one
c         of several means,
c     or
c
c     (2) The complete specification of the lattice, with scaling
c         parameters
c
c     which may be modified by
c
c     (3) The traditional elastic strain parameters (e_i, i=[1,6]),
c         which may be specified directly or in terms a convenient
c         volume conserving strain parameter.
c
c=======================================================================
c
c     REVISION HISTORY
c
c-----------------------------------------------------------------------
c
c     Switched from left-handed to right-handed coordinate systems
c      for hexagonal lattices.              --  mjm  26 March 1998
c
c-----------------------------------------------------------------------
c
c     Pulled from findmin, rewritten to provide lattice and volume
c      based on data from the static SKIN file.  In addition, a
c      menu of choices for the e_i has been added.
c
c     Note that to keep things in the same order as static, we have
c      reversed the meaning of the alat indices from the way they
c      appeared in findmin.  Now, alat(i,j) is the jth Cartesian
c      component of the ith primitive lattice vector.
c                                           --  mjm  21 April 1998
c
c-----------------------------------------------------------------------
c
c     Inadvertantly defined the fraction two3rd as "1/3".  Fixed
c      this.                                --  mjm  27 May   1998
c
c-----------------------------------------------------------------------
c
c     Added support for the "magic strain" lattices for centered cubic
c      systems.  These are in options 15 and 16. -- mjm  26 June 2003
c
c=======================================================================
c
      implicit none
c
c=======================================================================
c
c     Arguments
c
c     ierr = 0 if the lattice was successfully constructed, otherwise
c      ierr has a value corresponding to an error number given
c      at the bottom of this subroutine, and the primitive vectors
c      return with strange values
c
      integer ierr
c
      real*8 vol,alat(3,3)
c
c=======================================================================
c
c
c     Elastic strain parameters.  If strain is non-zero, then after
c      we construct the lattice we strain it by the parameters e(i).
c      Note that this may interfer with c/a and b/a options.  However,
c      this ability is so useful that I want to put it in, and the
c      user will have to beware of the problems.
c     We have modified the calculation to provide a "menu" of
c      predefined strain patterns which depend on only one variable,
c      either a strain x or on x2 = x*x.  These applied strains
c      may be used, for example, in elastic constant calculations.
c
      integer strain
      real*8 x,x2
c
c     And here are the real e(i):
c
      real*8 e(6)
c
c=======================================================================
c
c     Reals
c
c     a, b, c, boa, coa, and theta are used for lattice setup
c
      real*8 a,b,c,boa,coa,theta
c
c     afcc, abcc, and rmagic are used with vol and theta in
c      the setup of the "magic strain" lattices
c     w, x, y, and z are auxillary quantities needed there.  x is
c      already declared.
c
      real*8 afcc,abcc,rmagic,w,y,z
c
c     scale(i) is the scaling parameter for the ith direction, used
c      only when we specify latnum = 0.  If scale(i) = 0 we will
c      reset it to the default of 1.
c
      real*8 scale(3)
c
c     est is the strain matrix.  alats temporarily holds the
c      strained lattice, and sum holds temporary variables
c
      real*8 est(3,3),alats(3,3),sum
c
c=======================================================================
c
c     Integers
c
c     latnum is the lattice input index
c
      integer latnum
c
c     Do loop variables and counters
c
      integer i,j,k
c
c=======================================================================
c
c     Numerical factors:
c
      real*8 zero,one,two,three,four
c
      parameter (zero  = 0d0)
      parameter (one   = 1d0)
      parameter (two   = 2d0)
      parameter (three = 3d0)
      parameter (four  = 4d0)
c
      real*8 third,thirdm,half,halfm,quarter,two3rd
c
      parameter (third  = one/three)
      parameter (thirdm = -third)
      parameter (half   = 5d-1)
      parameter (halfm  = -half)
      parameter (quarter = one/four)
      parameter (two3rd = two/three)
c
c     sqrt(3)/2 and its inverse
c
      real*8 fac1,fac1m,fac2
      parameter (fac1 = 0.866025403784438647d0)
      parameter (fac1m = -fac1)
      parameter (fac2 = one/fac1)
c
c     Circular parameters
c
      real*8 pi,tpi3
      parameter (pi   = 3.14159265358979323846d0)
      parameter (tpi3 = two*pi/three)
c
c=======================================================================
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     There are several ways to input the lattice, based on the flag
c      latnum
c
c     latnum                Method
c
c         0                 Read vectors in Cartesian Coordinates
c         1                 standard fcc lattice, input a
c        -1                 standard fcc lattice, input vol
c         2                 standard bcc lattice, input a
c        -2                 standard bcc lattice, input vol
c         3                 standard hexagonal lattice, input a,c
c        -3                 standard hexagonal lattice, input vol,c/a
c         4                 standard sc  lattice, input a
c        -4                 standard sc  lattice, input vol
c         5                 simple tetragonal lattice, input a,c
c        -5                 simple tetragonal lattice, input vol,c/a
c         6                 body-centered tetragonal lattice, a,c
c        -6                 body-centered tetragonal lattice, vol,c/a
c         7                 orthorhombic lattice, a,b,c
c        -7                 orthorhombic lattice, vol,b/a,c/a
c         8                 base-centered orthorhombic lattice,a,b,c
c        -8                 base-centered orthorhombic lattice, vol,b/a,c/a
c         9                 simple tetragonal lattice, rotated 45 degrees
c        -9                 simple tetragonal lattice, rotated, V,c/a
c        10                 face-centered tetragonal,
c       -10                 face-centered tetragonal, V,c/a
c        11                 face-centered tetragonal 'fcc-like', input a,c
c       -11                 face-centered tetragonal 'fcc-like, input V,c/a
c        12                 rhombohedral lattice, read a,b (see below)
c       -12                 rhombohedral lattice, read V, angle in radians
c        13                 body-centered tetragonal, 'bcc-like', a,c
c       -13                 body-centered tetragonal, 'bcc-like', vol,c/a
c        14                 base-centered orthorhombic lattice a_plane,c,theta
c       -14                 base-centered orthorhombic lattice V,c/a_plane,theta
c        15                 "magic strains" around fcc lattice, a_bcc,r,theta
c       -15                 "magic strains" around fcc lattice, V,r,theta
c        16                 "magic strains" around bcc lattice, a_fcc,r,theta
c       -16                 "magic strains" around bcc lattice, V,r,theta
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     At the start we have no error
c
      ierr = 0
c
c     Zero all of the primitive vectors.  We'll fill in the
c      ones we want to change:
c
      do i = 1,3
         do j = 1,3
            alat(j,i) = zero
         end do
      end do
c
      read(10,*,err=11000) latnum
      write( *,'(/''Lattice type '',i5)') latnum
      write(15,'(/''Lattice type '',i5)') latnum
c
c     Set up the primitive vectors:
c
      if(latnum.eq.0) then
c
c        Arbitrary lattice.  The lattice vectors are specified
c         directly.  Note the input1.f compatible indices
c
         write(15,'(''Arbitrary lattice'')')
         read(10,*,err=12000) ((alat(i,j),j=1,3),i=1,3)
         write(15,'(''a('',i1,'') = '',3f12.5)')
     $        (i,(alat(i,j),j=1,3),i=1,3)
c
c        For this lattice, and this lattice only, we allow
c         scale parameters.  Read these in now, and invoke
c         after we apply the elastic strains.  These will default
c         to 1 if scale(i) = 0
c
         read(10,*,err=12000) scale
         do i = 1,3
            if(scale(i).eq.zero) scale(i) = one
         end do
         write(15,'(''Scale factors = '',3f10.5)') scale
      else if(abs(latnum).eq.1) then
c
c        fcc lattice
c
         write(15,'(''fcc lattice'')')
         if(latnum.eq.1) then
            read(10,*,err=12000) a
            write(15,'('' a = '',f12.5)') a
         else
            read(10,*,err=12000) vol
            write(15,'('' V = '',f12.5)') vol
            a = (four*abs(vol))**third
         end if
c
c        Note that we don't have to reverse the indicies here,
c         as the primitive lattice "matrix" is symmetric
c
         do i = 1,3
            do j = 1,3
               if(j.ne.i) alat(j,i) = half*a
            end do
         end do
      else if(abs(latnum).eq.2) then
c
c        bcc lattice
c
         write(15,'(''bcc lattice'')')
         if(latnum.eq.2) then
            read(10,*,err=12000) a
            write(15,'('' a = '',f12.5)') a
         else
            read(10,*,err=12000) vol
            write(15,'('' V = '',f12.5)') vol
            a = (two*abs(vol))**third
         end if
         do i = 1,3
            do j = 1,3
               if(j.eq.i) then
                  alat(j,i) = -half*a
               else
                  alat(j,i) = half*a
               end if
            end do
         end do
      else if(abs(latnum).eq.3) then
c
c        hexagonal lattice
c
         write(15,'(''hexagonal lattice'')')
         if(latnum.eq.3) then
            read(10,*,err=12000) a,c
            write(15,'('' a = '',f12.5,'' c = '',f12.5)') a,c
            coa = c/a
         else
            read(10,*,err=12000) vol,coa
            write(15,'('' V = '',f12.5,'' c/a = '',f10.6)') vol,coa
            a = (fac2*abs(vol)/abs(coa))**third
            c = abs(coa)*a
         end if
c
c        Formerly this was a left-hand lattice, but let's change it
c         to a right-hand lattice for consistency with our space
c         Crystal Lattice page
c
c        Here we do have to reverse the indices from the findmin
c         order
c
         alat(1,1) = half*a
         alat(2,1) = alat(1,1)
         alat(1,2) = fac1m*a
         alat(2,2) = -alat(1,2)
         alat(3,3) = c
      else if(abs(latnum).eq.4) then
c
c        sc lattice
c
         write(15,'(''simple cubic lattice'')')
         if(latnum.eq.4) then
            read(10,*,err=12000) a
            write(15,'('' a = '',f12.5)') a
         else
            read(10,*,err=12000) vol
            write(15,'('' V = '',f12.5)') vol
            a = (abs(vol))**third
         end if
         do i = 1,3
            alat(i,i) = a
         end do
      else if(abs(latnum).eq.5) then
c
c        Simple tetragonal lattice
c
         write(15,'(''tetragonal lattice'')')
         if(latnum.eq.5) then
            read(10,*,err=12000) a,c
            write(15,'('' a = '',f12.5,'' c = '',f12.5)') a,c
         else
            read(10,*,err=12000) vol,coa
            write(15,'('' V = '',f12.5,'' c/a = '',f10.6)') vol,coa
            a = (abs(vol)/abs(coa))**third
            c = coa*a
         end if
         alat(1,1) = a
         alat(2,2) = a
         alat(3,3) = c
      else if(abs(latnum).eq.6) then
c
c        Body-centered tetragonal lattice
c
         write(15,'(''body-centered tetragonal lattice'')')
         if(latnum.eq.6) then
            read(10,*,err=12000) a,c
            write(15,'('' a = '',f12.5,'' c = '',f12.5)') a,c
         else
            read(10,*,err=12000) vol,coa
            write(15,'('' V = '',f12.5,'' c/a = '',f10.6)') vol,coa
            a = (two*abs(vol)/abs(coa))**third
            c = coa*a
         end if
c
c        Here we again must reverse the indices from the findmin
c         order
c
         alat(1,1) = a
         alat(2,2) = a
         alat(3,1) = half*a
         alat(3,2) = half*a
         alat(3,3) = half*c
      else if(abs(latnum).eq.7) then
c
c        Orthorhombic lattice
c
         write(15,'(''orthorhombic lattice'')')
         if(latnum.eq.7) then
            read(10,*,err=12000) a,b,c
            write(15,'('' a = '',f12.5,'' b = '',f12.5,
     $           '' c = '',f12.5)') a,b,c
         else
            read(10,*,err=12000) vol,boa,coa
            write(15,'('' V = '',f12.5,'' b/a = '',f10.6,
     $           '' c/a = '',f10.6)') vol,boa,coa
            a = (abs(vol)/(abs(boa*coa)))**third
            b = boa*a
            c = coa*a
         end if
         alat(1,1) = a
         alat(2,2) = b
         alat(3,3) = c
      else if(abs(latnum).eq.8) then
c
c        Base-Centered Orthorhombic lattice.  Primitive vectors are
c         Again we reverse the indices from findmin's setlat order:
c
c        a1 = ( a/2 , -b/2 , 0 )
c        a2 = ( a/2 ,  b/2 , 0 )
c        a3 = (  0  ,   0  , c )
c
         write(15,'(''base-centered orthorhombic lattice'')')
         if(latnum.eq.8) then
            read(10,*,err=12000) a,b,c
            write(15,'('' a = '',f12.5,'' b = '',f12.5,
     $           '' c = '',f12.5)') a,b,c
         else
            read(10,*,err=12000) vol,boa,coa
            write(15,'('' V = '',f12.5,'' b/a = '',f10.6,
     $           '' c/a = '',f10.6)') vol,boa,coa
            a = (abs(two*vol)/(abs(boa*coa)))**third
            b = boa*a
            c = coa*a
         end if
         alat(1,1) = half*a
         alat(2,1) = alat(1,1)
c
c        This coordinate system will be right-handed:
c
         alat(2,2) = half*b
         alat(1,2) = -alat(2,2)
         alat(3,3) = c
      else if(abs(latnum).eq.9) then
c
c        Simple tetragonal lattice again, except that the x-y
c         plane primitive vectors are along (110) and (1-10).  This
c         is for compatibility with L1_0 using the t2g-eg split
c         parametrization
c
c        E.g., the primitive vectors are
c
c        a1 = ( a/2 , -a/2 , 0 )
c        a2 = ( a/2 ,  a/2 , 0 )
c        a3 = (  0  ,   0  , c )
c
c        This is analogous to the base-centered orthorhombic lattice
c         above, with b = a.
c
         write(15,'(''rotated tetragonal lattice'')')
         if(latnum.eq.9) then
            read(10,*,err=12000) a,c
            write(15,'('' a = '',f12.5,'' c = '',f12.5)') a,c
         else
            read(10,*,err=12000) vol,coa
            write(15,'('' V = '',f12.5,'' c/a = '',f10.6)') vol,coa
            a = abs(two*vol/coa)**third
            c = coa*a
         end if
         alat(1,1) =  half*a
         alat(1,2) = -alat(1,1)
         alat(2,1) =  alat(1,1)
         alat(2,2) =  alat(1,1)
         alat(3,3) =  c
      else if(abs(latnum).eq.10) then
c
c        face-centered tetragonal lattice, with the x-y
c         plane primitive vectors are along (110) and (1-10).  This
c         is for compatibility with L1_0 using the t2g-eg split
c         parametrization
c
         write(15,'(''face-centered tetragonal lattice'')')
         if(latnum.eq.10) then
            read(10,*,err=12000) a,c
            write(15,'('' a = '',f12.5,'' c = '',f12.5)') a,c
         else
            read(10,*,err=12000) vol,coa
            write(15,'('' V = '',f12.5,'' c/a = '',f10.6)') vol,coa
            a = abs(four*vol/coa)**third
            c = coa*a
         end if
         alat(1,1) =  half*a
         alat(1,2) = -alat(1,1)
         alat(2,1) =  alat(1,1)
         alat(2,2) =  alat(1,1)
         alat(3,1) =  alat(1,1)
         alat(3,3) =  half*c
      else if(abs(latnum).eq.11) then
c
c        Face-centered tetragonal lattice -- an alternate form, with
c        "fcc-like" primitive vectors.
c
         write(15,'(''fcc-like face-centered tetragonal lattice'')')
         if(latnum.eq.11) then
            read(10,*,err=12000) a,c
            write(15,'('' a = '',f12.5,'' c = '',f12.5)') a,c
         else
            read(10,*,err=12000) vol,coa
            write(15,'('' V = '',f12.5,'' c/a = '',f10.6)') vol,coa
            a = (four*abs(vol)/abs(coa))**third
            c = coa*a
         end if
         alat(2,1) = half*a
         alat(3,1) = half*a
         alat(1,2) = half*a
         alat(3,2) = half*a
         alat(1,3) = half*c
         alat(2,3) = half*c
      else if(abs(latnum).eq.12) then
c
c        Rhombohedral lattice.  Use the form
c
c        a1 = (a,b,b)
c        a2 = (b,a,b)
c        a3 = (b,b,a)
c
c        Or specify the volume V and the angle theta in radians
c
         write(15,'(''rhombohedral lattice'')')
         if(latnum.eq.12) then
            read(10,*,err=12000) a,b
            write(15,'('' a = '',f12.5,'' b = '',f12.5)') a,b
         else
            read(10,*,err=12000) vol,theta
            write(15,'('' V = '',f12.5,'' theta = '',f10.6)') vol,theta
c
c           Use coa to store the cosine of the angle
c
            coa = cos(theta)
c
c           If coa < -.5, we've got a bad angle
c
            if(coa.lt.halfm) stop 'Angle to large for lattice 12'
c
c           boa holds the auxillary lattice parameter
c
            boa = third*(sqrt((one+two*coa)/(one-coa))-one)
c
c           We'll let c do some work, too:
c
            c = (abs(vol)/(one+three*boa))**third
c
c           Finally
c
            a = c*(one+boa)
            b = c*boa
         end if
         do i = 1,3
            do j = 1,3
               if(j.eq.i) then
                  alat(j,i) = a
               else
                  alat(j,i) = b
               end if
            end do
         end do
      else if(abs(latnum).eq.13) then
c
c        This is Yet Another Version of the Body-Centered Tetragonal
c         lattice (YAVBCT).  In this version the primitive vectors
c         are in a form reminiscent of the standard bcc primitive
c         vectors:
c
c        a1 = (-a/2, a/2, c/2)
c        a2 = ( a/2,-a/2, c/2)
c        a3 = ( a/2, a/2,-c/2)
c
         write(15,'(''bcc-like bct lattice'')')
         if(latnum.eq.13) then
            read(10,*,err=12000) a,c
            write(15,'('' a = '',f12.5,'' c = '',f12.5)') a,c
         else
            read(10,*,err=12000) vol,coa
            write(15,'('' V = '',f12.5,'' c/a = '',f10.6)') vol,coa
            a = (two*abs(vol)/abs(coa))**third
            c = coa*a
         end if
         alat(1,1) = -half*a
         alat(1,2) = half*a
         alat(1,3) = half*c
         alat(2,1) = half*a
         alat(2,2) = -half*a
         alat(2,3) = half*c
         alat(3,1) = half*a
         alat(3,2) = half*a
         alat(3,3) = -half*c
      else if (abs(latnum).eq.14) then
c
c        Base-centered orthorhombic lattice.  "a" is the length
c         of one of the two primitive vectors in the x-y plane,
c         and "theta" is the angle between them (in radians).
c         c is the length of the primitive vector along the (001)
c         direction.
c
c        Note that V = a^2 c sin(theta)
c
         if(latnum.eq.14) then
            read(10,*,err=12000) a,c,theta
            write(15,'(''a = '',f12.5,'' c = '',f12.5,
     $           '' theta = '',f12.5)') a,c,theta
            vol = a*a*c*sin(theta)
         else
            read(10,*,err=12000) vol,coa,theta
            write(15,'(''V = '',f12.5,'' c/a = '',f12.5,
     $           '' theta = '',f12.5)') vol,coa,theta
            a = (vol/(coa*sin(theta)))**third
            c = a*coa
         end if
c
c        For the rest of the calculation we'll need theta/2 rather
c         than theta:
c
         theta = half*theta
c
         alat(1,1) =  a*cos(theta)
         alat(1,2) = -a*sin(theta)
         alat(2,1) =  alat(1,1)
         alat(2,2) = -alat(1,2)
         alat(3,3) =  c
      else if (abs(latnum).eq.15) then
c
c        This is the "magic strain" lattice _based_ on the fcc lattice.
c         That is, the central point around which everything is based,
c         (rmagic = 0, theta = whatever) is an fcc lattice.  However,
c         this is most useful in describing magic strains for materials
c         with bcc ground states.
c
c        The input is either the volume or the equivalent bcc lattice
c         constant (V = abcc^3/2), rmagic and theta.  The primitive
c         vectors are:
c
c         (a1)               ( 0 y z )
c         (a2) = (V/2)^(1/3) ( x 0 z )
c         (a3)               ( x y 0 )
c
c        where
c
c        x = exp( - rmagic cos(theta + 2 pi/3)
c        y = exp( - rmagic cos(theta - 2 pi/3)
c        z = exp( - rmagic cos(theta))
c
c        So the lattice has a period of 2 pi/3.
c
c        Note that we get a bcc lattice when theta = 0 (or +/- 2 pi/3)
c         and rmagic = 1/3 ln 2 = 0.2310490601866 ..
c
c        Enter theta in radians
c
         if(latnum.eq.15) then
            read(10,*,err=12000) abcc,rmagic,theta
            vol = half*(abcc**3)
         else
            read(10,*,err=12000) vol,rmagic,theta
            abcc = (two*vol)**third
         end if
         write(15,'('' V = '',f12.5,'' abcc = '',f12.5,
     $        '' rmagic = '',f12.7,'' theta = '',f12.7)')
     $        vol, abcc, rmagic, theta
c
c        w is the prefactor in front of every lattice vector:
c
         w = (half*vol)**third
c
c        Define x, y, and z including the w:
c
         x = w*exp(-rmagic*cos(theta + tpi3))
         y = w*exp(-rmagic*cos(theta - tpi3))
         z = w*exp(-rmagic*cos(theta))
c
c        Define the non-zero components of the primitive lattice
c         vectors:
c
         alat(1,2) = y
         alat(1,3) = z
         alat(2,1) = x
         alat(2,3) = z
         alat(3,1) = x
         alat(3,2) = y
      else if (abs(latnum).eq.16) then
c
c        This is the "magic strain" lattice _based_ on the bcc lattice.
c         That is, the central point around which everything is based,
c         (rmagic = 0, theta = whatever) is a bcc lattice.  However,
c         this is most useful in describing magic strains for materials
c         with fcc ground states.
c
c        The input is either the volume or the equivalent bcc lattice
c         constant (V = afcc^3/2), rmagic and theta.  The primitive
c         vectors are:
c
c         (a1)               ( -x  y  z )
c         (a2) = (V/4)^(1/3) (  x -y  z )
c         (a3)               (  x  y -z )
c
c        where
c
c        x = exp( rmagic cos(theta + 2 pi/3)
c        y = exp( rmagic cos(theta - 2 pi/3)
c        z = exp( rmagic cos(theta))
c
c        So the lattice has a period of 2 pi/3.
c
c        Note that we get an fcc lattice when theta = 0 (or +/- 2 pi/3)
c         and rmagic = 1/3 ln 2 = 0.2310490601866 ...
c
c        In addition, we get a bct10 (10-fold coordinated body-centered
c         tetragonal lattice) when rmagic = 1/3 ln 3/2 = 0.135155036036 ...
c         and theta = pi/3.
c
c        Enter theta in radians
c
         if(latnum.eq.16) then
            read(10,*,err=12000) afcc,rmagic,theta
            vol = quarter*(afcc**3)
         else
            read(10,*,err=12000) vol,rmagic,theta
            afcc = (four*vol)**third
         end if
         write(15,'('' V = '',f12.5,'' afcc = '',f12.5,
     $        '' rmagic = '',f12.7,'' theta = '',f12.7)')
     $        vol, afcc, rmagic, theta
c
c        w is the prefactor in front of every lattice vector:
c
         w = (quarter*vol)**third
c
c        Define x, y, and z including the w:
c
         x = w*exp(rmagic*cos(theta + tpi3))
         y = w*exp(rmagic*cos(theta - tpi3))
         z = w*exp(rmagic*cos(theta))
c
c        Define the components of the primitive lattice vectors:
c
         alat(1,1) = -x
         alat(1,2) =  y
         alat(1,3) =  z
         alat(2,1) =  x
         alat(2,2) = -y
         alat(2,3) =  z
         alat(3,1) =  x
         alat(3,2) =  y
         alat(3,3) = -z
      else
         write(0,*) 'Unknown lattice type'
         go to 10000
      end if
c
c=======================================================================
c
c
c     Strain the lattice if desired
c
c     We do this by entering a "strain option".  Note that some options
c      have negative as well as positive values.  Options with the same
c      absolute value have the same symmetry and small-strain behavior
c      of the energy, but differ at large strains.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     0 -- No strain
c
c     1 -- Read the six e(i) independently
c
c     2 -- Read in parameter x2 = x^2.  Construct volume conserving
c            strain:
c            e(1) = sqrt((1+x)/(1-x)) - 1 and
c            e(2) = sqrt((1-x)/(1+x)) - 1
c          Then
c            E(x) = E(0) + 1/2 V (C11 + C22 - 2 C12) x^2 + O(x^4)
c    -2 -- "Traditional" strain, using e(3) to control the volume.  Then
c            e(1) =  x
c            e(2) = -x
c            e(3) = x^2/(1-x^2)
c            Same Elastic constant representation.
c     21 -- Similar to above, but for use when the x and y axis behave
c           differently so that E(x) <> E(-x):
c            e(1) = sqrt((1+x)/(1-x)) - 1 and
c            e(2) = sqrt((1-x)/(1+x)) - 1
c          Then
c            E(x) = E(0) + 1/2 V (C11 + C22 - 2 C12) x^2 + O(x^3)
c    -21 -- The analog of -2:
c            e(1) =  x
c            e(2) = -x
c            e(3) = x^2/(1-x^2)
c            Same Elastic constant representation.
c
c     3 -- Read in x2 = e(6)^2.  Use e(6), e(1), and e(2) to construct
c            a volume conserving strain.  The energy then goes as
c            E(x) = E(0) + 1/2 V C66 x^2 + O(x^4).
c    -3 -- "Traditional" strain, using e(3) to control the volume
c            rather than e(1) and e(2).  Same Elastic constant
c            representation.
c
c     4 -- Let 1+x be the ratio of change of the length of the
c            z-coordinate to the change in the x-y plane coordinates.
c            If this sounds confusing, think of an [001] strain on
c            a cubic lattice.  Then the strained ratio c/a = 1+x
c            In this case the energy goes as
c            E(x) = E(0) + 1/18 V (C11 + C22 + 4 C33
c                              + 2 C12 - 4 C13 - 4 C23) x^2 + O(x^3)
c
c     5 -- Considering the standard picture of the e(i) representing
c            strains of the original Cartesian unit vectors, let x
c            represent the strain of these vectors along the [111]
c            direction.  The cosine between these vectors is
c            (2x + 3 x^2)/(1 + 2x + 3 x^2).  If we start with a
c            simple cubic lattice, then x = 1/3 is an fcc lattice
c            and x = -1/6 is a bcc lattice.  Note that at x = -1/3
c            the three strained unit vectors become coplanar, and
c            the problem blows up.
c            In the general case the energy goes as
c            E(x) = E(0) + 2 V (C44 + C55 + C66) x^2 + O(x^3)
c
c     6 -- One possible strain to determine C_{44} in a hexagonal
c            lattice is the volume conserving strain
c            e1 = e2 = e5 = e6 = 0
c            e6 = x ; e3 = 1 - x^2/4
c          This turns the hexagonal lattice into a (true) monoclinic
c            lattice.  Both the simple hexagonal and hcp lattice
c            are now reduced to C2/m (#12) symmetry.  In the hcp
c            lattice the atoms are allowed to relax by an arbitrary
c            vector of the form (0,y,z).  It is fairly easy to show
c            that the +x and -x lattices are equivalent, so the energy goes as
c            E(x) = E(0) + 1/2 V C44 x^2 + O(x^4)
c
c-----------------------------------------------------------------------
c
c     Note the following relationships between elastic constants
c
c     In a cubic system,
c
c        C11 = C22 = C33
c        C12 = C13 = C23
c        C44 = C55 = C66
c
c        so there are three independent lattice constants
c
c     In a tetragonal system
c
c        C11 = C22
c        C13 = C23
c        C44 = C55
c
c        so there are six independent lattice constants
c
c     In a hexagonal system
c
c        C11 = C22
c        C13 = C23
c        C44 = C55
c        C66 = 1/2 (C11-C12)
c
c        so there are five independent lattice constants
c
c-----------------------------------------------------------------------
c
c     For more details on the relationship between strains and elastic
c      constants, see "First-Principles Calculation of Elastic
c      Properties," by M.J. Mehl, B.M. Klein, and D.A. Papaconstantopoulos,
c      in "Intermetallic Compounds:  Vol. 1, Principles," edited by
c      J.H. Westbrook and R.L. Fleischer, (John Wiley & Sons, London,1994),
c      Chapter 9, pp. 195-210.  In particular, see equation (24) and
c      Table 1.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      read(10,*,err=13000) strain
      write( *,'(/''Elastic strain field '',i5)') strain
      write(15,'(/''Elastic strain field '',i5)') strain
c
      if(strain.ne.0) then
c
c        First zero the strain components:
c
         do i = 1,6
            e(i) = zero
         end do
c
c        Now read in strains as necessary:
c
         if(strain.eq.1) then
            read(10,*,err=13000) (e(i),i=1,6)
            write(15,'(''e_i = '',6f10.6)') (e(i),i=1,6)
         else if(abs(strain).eq.2) then
            read(10,*,err=13000) x2
            x = sqrt(x2)
            write(15,'(''x^2 = '',f12.8)') x2
            if(strain.eq.2) then
               e(1) = sqrt((one+x)/(one-x)) - one
               e(2) = sqrt((one-x)/(one+x)) - one
            else
               e(1) = x
               e(2) = -x
               e(3) = x2/(one-x2)
            end if
         else if(abs(strain).eq.21) then
            read(10,*,err=13000) x
            x2 = x*x
            write(15,'(''x = '',f12.8)') x
            if(strain.eq.21) then
               e(1) = sqrt((one+x)/(one-x)) - one
               e(2) = sqrt((one-x)/(one+x)) - one
            else
               e(1) = x
               e(2) = -x
               e(3) = x2/(one-x2)
            end if
         else if(abs(strain).eq.3) then
            read(10,*,err=13000) x2
            write(15,'(''x^2 = '',f12.8)') x2
            e(6) = sqrt(x2)
            if(strain.eq.3) then
               e(1) = sqrt(one+quarter*x2) - one
               e(2) = e(1)
            else
               e(3) = x2/(four-x2)
            end if
         else if(strain.eq.4) then
            read(10,*,err=13000) x
            write(15,'(''x = '',f12.8)') x
            e(1) = (one+x)**thirdm - one
            e(2) = e(1)
            e(3) = (one+x)**two3rd - one
         else if(strain.eq.5) then
            read(10,*,err=13000) x
            write(15,'(''x = '',f12.8)') x
            e(1) = (one+x)*((one+three*x)**thirdm) - one
            e(2) = e(1)
            e(3) = e(1)
            e(4) = two*x*((one+three*x)**thirdm)
            e(5) = e(4)
            e(6) = e(5)
         else if(strain.eq.6) then
c
c           Strained hexagonal lattice.  Normally the energy will
c            be independent of the sign of x, but well allow
c            negative x strains generated if x2 < 0
c
            read(10,*,err=13000) x2
            if(x2.lt.zero) then
               x2 = abs(x2)
               x = -sqrt(x2)
            else
               x = sqrt(x2)
            end if
            write(15,'(''x^2 = '',f12.8)') x2
            e(4) = x
            e(3) = quarter*x*x
         else
            write(0,*) 'Unknown strain type'
            go to 13000
         end if
c
c        Construct the strain matrix from the e(i):
c
         do i = 1,3
            est(i,i) = one + e(i)
         end do
         est(2,3) = half*e(4)
         est(3,2) = est(2,3)
         est(1,3) = half*e(5)
         est(3,1) = est(1,3)
         est(1,2) = half*e(6)
         est(2,1) = est(1,2)
c
c        We have reversed the strain indicies on alat and alats
c         (compare to the findmin version of setlat):
c
         do i = 1,3
            do j = 1,3
               sum = zero
               do k = 1,3
                  sum = sum + alat(j,k)*est(k,i)
               end do
               alats(j,i) = sum
            end do
         end do
         do j = 1,3
            do i = 1,3
               alat(i,j) = alats(i,j)
            end do
         end do
      end if
c
c=======================================================================
c
c     Now that we have the new lattice vectors, apply the scale
c      field if we have latnum = 0
c
      if(latnum.eq.0) then
         do j = 1,3
            do i = 1,3
               alat(i,j) = alat(i,j)*scale(j)
            end do
         end do
      end if
c
c=======================================================================
c
c     Calculate the volume:
c
      vol = abs(
     $     alat(1,1)*(alat(2,2)*alat(3,3) - alat(3,2)*alat(2,3)) +
     $     alat(2,1)*(alat(3,2)*alat(1,3) - alat(1,2)*alat(3,3)) +
     $     alat(3,1)*(alat(1,2)*alat(2,3) - alat(2,2)*alat(1,3)))
      if(vol.eq.zero) then
         write(0,*) 'Zero Volume'
         go to 10000
      end if
c
c     Standard return:
c
      return
c
c     If we get here, some kind of error has occured (probably
c      an EOF), and the calling program should be warned:
c
10000 write(*,*) 'setlat:  bad lattice'
      ierr = 1
      return
11000 write(*,*) 'setlat:  error in latnum input'
      ierr = 2
      return
12000 write(*,*) 'setlat:  error in lattice parameter input'
      ierr = 3
      return
13000 write(*,*) 'setlat:  error in strain input'
      ierr = 4
      return
      end

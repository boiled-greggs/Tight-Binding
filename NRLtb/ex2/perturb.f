      subroutine perturb(mxd,n,dhmat,dsmat,z,eig,eigcut,ntot,
     $     deig,work)
c
c     Calculate perturbed eigenvalues deig(i) on a non-orthogonal
c      Hamiltonian/Overlap matrix system specified
c      by eigenvalues eig(i) and states z(i,j).
c
c     Note that perturb and dgen form a package (now all in one routine)
c
c     Perturbations are calculated only if eig(k) <= eigcut
c
c=======================================================================
c
c     mxd is size of the matricies
c     n is the actual dimension of the secular equation
c
c     ntot is the number of perturbed eigenenergies returned, i.e.,
c      the number of eigenvalues < eigcut
c
c     dh is the perturbation on the Hamiltonian
c     ds is the perturbation on the Overlap Matrix
c
c     work is a temporary array used to hold intermediate values
c
c     The entire system should be Hermitian
c
c     Use first order perturbation theory.  (Outlined below)
c
c=======================================================================
c
      implicit none
c
      integer mxd,n,ntot
      real*8 eigcut
c$$$      real*8 dhr(mxd*(mxd+1)/2),dhi(mxd*(mxd+1)/2),
c$$$     $     dsr(mxd*(mxd+1)/2),dsi(mxd*(mxd+1)/2)
      complex*16 dhmat(mxd*(mxd+1)/2),dsmat(mxd*(mxd+1)/2)
      complex*16 z(mxd,mxd),work(mxd)
      real*8 eig(mxd),deig(mxd)
c
c=======================================================================
c
c     Assume the maximum degeneracy is mdg.
c     num is the actual degeneracy at any level
c
      integer mdg,num
      parameter (mdg = 9)
c
c     zd should hold the eigenvectors of hd, and should
c      be dimensioned to (mdg,mdg).  However, we don't need the
c      eigenvectors of hd, so we can save a (trivial) bit of
c      storage
c
      complex*16 zd(1,1),zwork(2*mdg-1)
      complex*16 hd(mdg*(mdg+1)/2)
c$$$      complex*16 dcmplx,dconjg
      real*8 seig(mdg),rwork(3*mdg-2)
      integer info
c
      integer i,j,is,js,k,m,start
      complex*16 csum
      real*8 sqsum
c
c     mjk is used to index the Hamiltonian matrix
c
      integer mjk
c
c=======================================================================
c
      real*8 zero,half,four
      parameter (zero = 0d0)
      parameter (half = 5d-1)
      parameter (four = 4d0)
c
c     Tolerance.  Determines whether two states are degenerate or
c      not
c
      real*8 tol
      parameter (tol = 1d-8)
c
c=======================================================================
c
c     First determine the degenerate and non-degenerate states
c
c     Degenerate and non-degenerate
c
      start = 0
 100  start = start + 1
c
c     If the current eigenvalue is too large, then we are finished:
c
      if(eig(start).gt.eigcut) then
         ntot = start - 1
         return
      end if
c
c     How many degenerate states are there at this energy (if any)
      do j = start+1,n
         if(abs(eig(j)-eig(start)).gt.tol) then
            num = j - start
            go to 200
         end if
      end do
      num = n - start + 1
c
 200  if(num.gt.mdg) stop 'dgen:  Degeneracy is too large'
c
c     The first order perturbation energies are the
c      eigenvalues of the matrix
c
c     <a| dh - eig ds |b>
c
c     where {|a>} are the (<a|S_0|b> = delta(a,b)) eigenvectors
c      which span the degenerate subspace of the unperturbed Hamiltonian
c
c     Form the matrix.  Use the packed form
c
      m = 0
      do i = 1,num
c
c        The states are |start>, |start+1>, ... |start+num-1>
c
c        The unperturbed eigenvalue is always eig(start)
c
         is = start + i - 1
c
         do j = 1,n
            csum = dcmplx(zero,zero)
c           Note that (dhr, dhi) = hmat, etc. are now packed, so
c            we have have this little overhead to fix
c            things up.  In particular, we have to respect the
c            Hermitian nature of the matrices.
c
c$$$            mjk = j*(j-1)/2
            do k = 1,j-1
c
c              Case j > k: mjk = j*(j-1)/2+k
c
c$$$               mjk = mjk + 1
               mjk = j*(j-1)/2+k
c
c              Complex conjugation added:
c
c$$$               csum = csum + (dcmplx(dhr(mjk),-dhi(mjk))
c$$$     $              - eig(start)*dcmplx(dsr(mjk),-dsi(mjk)))*z(k,is)
               csum = csum +
     $              dconjg(dhmat(mjk)-eig(start)*dsmat(mjk))*z(k,is)
            end do
            do k = j,n
c
c              Case k >= j:  mjk = k*(k-1)/2+j
c
c$$$               mjk = mjk + k-1
               mjk = k*(k-1)/2+j
c
c              No conjugation here:
c
c$$$               csum = csum + (dcmplx(dhr(mjk),dhi(mjk))
c$$$     $              - eig(start)*dcmplx(dsr(mjk),dsi(mjk)))*z(k,is)
               csum = csum + (dhmat(mjk)-eig(start)*dsmat(mjk))*z(k,is)
            end do
            work(j) = csum
         end do
         do j = 1,i
            m = m + 1
            js = start + j - 1
            csum = dcmplx(zero,zero)
            do k = 1,n
               csum = csum + dconjg(z(k,js))*work(k)
            end do
            hd(m) = csum
ctemp
c$$$            write(15,*) m,hd(m)
cend temp
         end do
      end do
c
c     Find the eigenvalues.  Note that it's a standard Hermetian problem
c      (no overlap matrix), and the diagonal elements are supposed
c      to be real
c
      if(num.eq.1) then
         deig(start) = dble(hd(1))
      else if(num.eq.2) then
c
c        This is a simple 2x2 matrix, so we might as well do it ourselves
c
         sqsum = sqrt(dble(hd(1)-hd(3))**2 + four*
     $        dble(dconjg(hd(2))*hd(2)))
         deig(start)   = half*(dble(hd(1)+hd(3)) - sqsum)
         deig(start+1) = half*(dble(hd(1)+hd(3)) + sqsum)
ctemp
c$$$         write(15,*) '2X2 perturbation matrix'
c$$$         do i = 1,3
c$$$            write(15,*) i,hd(i)
c$$$         end do
cend temp
      else
         call zhpev('N','U',num,hd,seig,zd,mdg,zwork,rwork,info)
         if(info.ne.0) then
            write(0,*) 'zhpev INFO = ',info,' in perturb for num = ',num
            stop 'PERTURB:  INFO'
         end if
         do i = 1,num
            deig(start+i-1) = seig(i)
         end do
      end if
c
c     Search through the next level of states
c
ctemp
c$$$      write(15,15) (m,deig(m),m=start,start+num-1)
c$$$ 15   format(i5,f15.8)
cend temp
      start = start + num - 1
      if (start.lt.n) go to 100
c
c     That should do it.  All states are counted if we get here:
c
      ntot = n
      return
      end

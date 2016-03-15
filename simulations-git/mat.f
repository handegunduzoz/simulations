c-----------------------------------------------------------------------
c     MATERIAL MODEL FOR CMP ---------------------> neohooke, PF
c     (c) C. Miehe, University of Stuttgart
c
c-----------------------------------------------------------------------
      subroutine mate01(d,g, s,aa, nnn,l,isw)
      implicit double precision (a-h,o-z)
c
c.... FEAP common blocks
      include 'iofile.h'
      include 'hdata.h'
      include 'comblk.h'
      include 'tdata.h'
      dimension d(*),g(24),s(24),aa(24,24)
c
c.... go to correct processor
      go to(1,2,3,3,2, 3,2,3,2,2, 2,2,2,3,2, 2,2,2,2,2), isw
c
c.... [nd0] general material-element-model information
1     nd0   = d(1)
      ndele = d(2)
      nhele = d(3)
c.... ndmat - number of material parameters
      ndmat = 13
c.... nhmat - number of history parameters
      nhmat = 1
      d(4)  = ndmat
      d(5)  = nhmat
c.... [ndmat] individual material-model information
      call dinput(d(nd0+ndele+1),7)
      call dinput(d(nd0+ndele+8),6)
c.... output of material data information
      write(*  ,2000) ndmat,nhmat
      write(iow,2000) ndmat,nhmat
      write(iow,2001) (d(nd0+ndele+i),i=1,7)
      write(  *,2001) (d(nd0+ndele+i),i=1,7)
      write(iow,2002) (d(nd0+ndele+i),i=8,13)
      write(  *,2002) (d(nd0+ndele+i),i=8,13)
2     return
c
c.... EXECUTION PHASE
3     nd0   = d(1)
      ndele = d(2)
      nhele = d(3)
      ndmat = d(4)
      nhmat = d(5)
      lmat = nd0 + ndele + 1
      lnh1 = nh1 + (l-1)*nhmat + nhele
      lnh2 = nh2 + (l-1)*nhmat + nhele
c.. set history zero at begin of process (for parameter identification)
      if (ttim.le.dt) then
         do  i = 1,nhmat
            hr(lnh1-1+i) = 0.0d0
         enddo
      endif
c
      call modl01(d(lmat),hr(lnh1),hr(lnh2),nhmat,g,nnn,l,isw, s,aa)
c
      return
c
 2000 format(
     * 10x,''/
     * 10x,'MATERIAL MODEL FOR CMP ----- (c) C. Miehe, Stuttgart'/
     * 10x,'[  ] [           ] 2fgradient large strain neo hooke',/
     * 10x,'[  ] [           ] SC form vol+iso:neo-noniso.f.....',/
     * 10x,'[G4] [ ndmat     ] material parameters .............',i12/
     * 10x,'[G5] [ nhmat     ] material history variables ......',i12)
 2001 format(
     * 10x,'[ 1] [ mu        ] shear modulus....................',e12.5/
     * 10x,'[ 2] [ nu        ] querkontraktion..................',e12.5/
     * 10x,'[ 3] [ scrit     ] critical stress level... ........',e12.5/
     * 10x,'[ 4] [ zeta      ] stiffness of d-evolution.........',e12.5/
     * 10x,'[ 5] [ l         ] length-scale parameter ..........',e12.5/
     * 10x,'[ 6] [ eta       ] viscosity .......................',e12.5/
     * 10x,'[ 7] [   k       ] artificial rest stiffness........',e12.5)
 2002 format(
     * 10x,'[ 8] [    c_v] specific heat .......................',e12.5/
     * 10x,'[ 9] [   xkon] conductivity  .......................',e12.5/
     * 10x,'[10] [   expa] thermal expansion coefficient .......',e12.5/
     * 10x,'[11] [temp_oo] ambient temperature surrounding mate.',e12.5/
     * 10x,'[12] [    h  ] surface convection coefficinent......',e12.5/
     * 10x,'[13] [penalty] penalty to prescribe temp in damage..',e12.5)
      end
c
c***********************************************************************
c
      subroutine modl01(d,hold,hnew,nhmat,g,nnn,lll,isw, s,aa)
      implicit double precision (a-h,o-z)
c
      dimension hold(*), hnew(*)
      dimension d(*),g(24),s(24),aa(24,24)
      dimension fi(3,3),f(3,3),fn(3,3),ci(3,3),xn(2),xm(2)
      dimension xs(3),xc(3,3),sign(3),sig(3), bm(3,3),ss(4),s0(4)
      dimension e(6),eeig(4),c(3,3),fe(3,3)
      dimension ii(9), jj(9),xii(3,3),ha(3,3,3,3),stau(9),sp(3,3)
      dimension ij(3,3)
c
      include 'tdata.h'
      include 'eldata.h'
      include 'hdata.h'
      include 'comblk.h'
      include 'pointer.h'
      include 'cdata.h'
      include 'iofile.h'
      include 'pdata6.h'
      include 'part0.h'
      save
      data ii/1,1,1,2,2,2,3,3,3/, jj/1,2,3,1,2,3,1,2,3/
      data xii/1.d0, 0.d0, 0.d0,
     &         0.d0, 1.d0, 0.d0,
     &         0.d0, 0.d0, 1.d0/
      data ij/1,4,6,7,2,5,9,8,3/

c.... material parameters
      p_mu = d(1) ! shear modulus
      p_nu = d(2) ! querkontraktion
      p_cr = d(3) ! critical fracture stress
      p_ze = d(4) ! unused parameter
      p_ls = d(5) ! length scale
      p_ep = d(6) ! viscosity
      p_k =  d(7) ! rest stiffness
      p_be = 2.d0*p_nu/(1.d0-2.d0*p_nu)
      p_c  = d(8)  ! specific heatr
      p_kon= d(9)  ! konductivity
      p_aexp= d(10) ! thermal expansion coefficient alpha
      t_amb0 = d(11) ! ambient temperature
      t_amb = d(11) ! ambient temperature
c      if(ttim.le.3.d0) then  ! special trick for circle example
c      t_amb = t_amb0*ttim*1.d2
c      else
c      t_amb = t_amb0*3.d0*1.d2
c      endif
      p_h   = d(12)  ! surface convection coefficient
      p_pen = d(13)  ! unused
c
c.... damage function
      dam    =         (1.d0-g(13))**2.d0 + p_k
      p_kon_d =     p_kon   !*(1.d0-g(13))**2.d0
c.... compute deformation gradient
      call pzero(f,9)
      call pzero(fe,9)
      call pzero(c,9)
      f(1,1)    =  g(1) + 1.d0
      f(1,2)    =  g(2)
      f(1,3)    =  g(3)
      f(2,1)    =  g(4)
      f(2,2)    =  g(5) + 1.d0
      f(2,3)    =  g(6)
      f(3,1)    =  g(7)
      f(3,2)    =  g(8)
      f(3,3)    =  g(9) + 1.d0
c.... volumetric expansion
      vol = dexp(-p_aexp*(g(18)-t_amb0))
       do i = 1,3
       do j = 1,3
       fe(i,j) = vol*f(i,j)
       enddo
       enddo

c.... compute determinante of deformation gradient
      call ComputeDeterminant(fe,3, detf)
       djb = detf**(-p_be)
c.... compute inverse of deformation gradient
      call pmove(fe,fi,9)
      call invert(fi,3,3)
      do i=1,3
      do j= 1,3
      do k=1,3
      c(i,j)=c(i,j)+ f(k,i)*f(k,j)
      enddo
       enddo
      enddo
      call pmove(c,ci,9)
      call invert(ci,3,3)
ccc
ccc       call pzero(ci,9)
ccc        ci(1,1) = 1.d0
ccc        ci(2,2) = 1.d0
ccc        ci(3,3) = 1.d0
ccc        
      if(npart.eq.2.or.npart.eq.3) goto 666
c     ---------------------
c.... C) GENERALIZED STRESS ROUTINE
c.... mechanic stresses and tangent
       do 10 i = 1,3
       do 10 j = 1,3
 10    sp(i,j) =  p_mu * ( fe(i,j) - djb * fi(j,i) )
       do 15 i=1,9
 15       s(i) = vol*dam*sp(ii(i),jj(i))  
c
c      kirchhoff spannung tau
       call pzero(stau,9)
       do 20 i = 1,3
       do 20 j = 1,3
       do 20 k = 1,3
 20    stau(ij(i,j)) = stau(ij(i,j))+ sp(i,k)*fe(j,k)*vol
       call spec_cmp(6,stau,eeig,bm)
c
       do 30 i=1,3
       do 30 j=1,3
       do 30 k=1,3 
       do 30 l=1,3
 30       ha(i,j,k,l) = p_mu * xii(i,k)*xii(j,l) 
     A            + p_mu*p_be*djb*fi(j,i)*fi(l,k)
     A            + p_mu*djb*fi(j,k)*fi(l,i)
      do 40 i=1,9
      do 40 j=1,9
 40   aa(i,j) = vol*vol*dam*ha(ii(i),jj(i),ii(j),jj(j))
c.... phase field driving force
      hist = 0.d0
      do i=1,3
      if(eeig(i).gt.0.d0) then
      hist = hist  + (eeig(i)/p_cr)**2.d0
      endif
      enddo
      hist = p_ze*(hist - 1.d0)
      if(hist.gt.1.d3) hist = 1.d3
       
c.... update history field
      if(hist.ge.hold(1)) then
        hnew(1) = hist
      else
        hnew(1) = hold(1)
      endif
c
      s(20)   = eeig(1)
      s(21)   = eeig(2)
      s(22)   = eeig(3)
      s(23)   = hist
      s(24)   = ttim
c
 666  continue
c.... phase field stress and tangent
      s(10) = g(10)*p_ls**2.d0
      s(11) = g(11)*p_ls**2.d0
      s(12) = g(12)*p_ls**2.d0
c
      s(13) = g(13) - p_ls*(1.d0-g(13))*hold(1)
     A     + (p_ep/dt)*(g(13)-g(14))*p_ls
c
      aa(10,10) = p_ls**2.d0
      aa(11,11) = p_ls**2.d0
      aa(12,12) = p_ls**2.d0
      aa(13,13) = 1.d0 + p_ls*hold(1) + (p_ep/dt)*p_ls
c
c===================================thermal======================
c .... surface density function
      gamma = g(13)**2.d0/2.d0/p_ls + p_ls/2.d0*(
     1     g(10)**2.d0 + g(11)**2.d0+   g(12)**2.d0)  
      s(15) = p_kon_d*(ci(1,1) * g(15) + ci(1,2)*g(16) + ci(1,3)*g(17))
      s(16) = p_kon_d*(ci(2,1) * g(15) + ci(2,2)*g(16) + ci(2,3)*g(17))
      s(17) = p_kon_d*(ci(3,1) * g(15) + ci(3,2)*g(16) + ci(3,3)*g(17))
      s(18) = p_c/dt*(g(18)-g(19)) + 2.d0*p_h*gamma*(g(18) - t_amb)
      do i=1,3
      do j=1,3
      aa(14+i,14+j) = p_kon_d*ci(i,j)
      enddo
      enddo
      aa(18,18) = p_c/dt  + 2.d0*p_h*gamma
cc
      s(24) = gamma     
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ComputeDeterminant(a,dim, det)

      implicit none

      integer dim
      real*8  a(dim,dim),det

      if (dim.eq.2) then
        det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      elseif (dim.eq.3) then
        det = a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2))
     1      - a(1,2) * (a(2,1)*a(3,3) - a(2,3)*a(3,1))
     2      + a(1,3) * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
      endif

      end
c

c-----------------------------------------------------------------------
c     ELEMENT FOR FEAP -----------------------
c
c-----------------------------------------------------------------------
      subroutine elmt01(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
c-----------------------------------------------------------------------
c
      implicit double precision (a-h,o-z)
c
      include 'pointer.h'
      include 'hdata.h'
      include 'comblk.h'
      include 'cdata.h'
      include 'tdata.h'
      include 'eldata.h'
      include 'iofile.h'
      include 'strnum.h'
c
      dimension sh(24),aah(24,24),s3(24),aa3(24,24), hh(5,18)
      dimension g(24)
      dimension d(*),ul(ndf,1),xl(ndm,1),ix(1),tl(*),s(nst,nst),p(*)
      dimension shp(4,8),sg(5,8)
      save
c
c ... go to correct processor
      go to(1,2,3,3,2, 3,2,3,2,2, 2,2,2,3,2, 2,2,2,2,2), isw
      return
c
c.... input [nd0] general informations
1     nd0   = 5
      ndele = 0
      nhele = 0
      d(1)  = nd0
      d(2)  = ndele
      d(3)  = nhele
c.... output of element data information
      write(*  ,2000) iel,nd0,ndele,nhele
      write(iow,2000) iel,nd0,ndele,nhele
c.... input [ndmat] individual material-model information
      call mate01(d,g, s3,aa3, n,l,isw)
      ndmat = d(4)
      nhmat = d(5)
c.... set number of element history variables
      ngaus = 8
      mct = nhele + ngaus*nhmat
c.... set number of projected stresses (iste=nr,istv=max nr)
        iste = 24
        istv = iste
2     return
c.... get number of parameters
3     nd0   = d(1)
      ndele = d(2)
      nhele = d(3)
      ndmat = d(4)
      nhmat = d(5)

c.... get gauss integration information (coords and weighting factors)
      igaus = 2
      if(nel.eq.4) then  
      call tint3d(igaus, lint,sg)
      elseif (nel.eq.8) then
      call int3d(igaus, lint,sg)
      endif
c
      do 300 l = 1,lint
c....   get shape functions and derivatives of shape functions
      if (nel.eq.4) then      
      call tetshp(sg(1,l),xl,ndm,nel,detj,shp)
      elseif (nel.eq.8) then
      call shp3d(sg(1,l),detj,shp,xl,ndm,nel)
      endif
        dvol = detj*sg(5,l)
        
c.... generalized strains
      call pzero(g,24)
c 
      do 310 i = 1,nel
c.... displacement gradient (grad u)
      g(1) = g(1) + ul(1,i)*shp(1,i)
      g(2) = g(2) + ul(1,i)*shp(2,i)
      g(3) = g(3) + ul(1,i)*shp(3,i)
      g(4) = g(4) + ul(2,i)*shp(1,i)
      g(5) = g(5) + ul(2,i)*shp(2,i)
      g(6) = g(6) + ul(2,i)*shp(3,i)
      g(7) = g(7) + ul(3,i)*shp(1,i)
      g(8) = g(8) + ul(3,i)*shp(2,i)
      g(9) = g(9) + ul(3,i)*shp(3,i)
c.... damage gradient (grad d)
      g(10) = g(10) + ul(4,i)*shp(1,i)
      g(11) = g(11) + ul(4,i)*shp(2,i)
      g(12) = g(12) + ul(4,i)*shp(3,i)
c.... damage field d
      g(13) = g(13) + ul(4,i)*shp(4,i)
c.... damage field d at time t_n
      g(14) = g(14) + (ul(4,i)-ul(4,i+nel))*shp(4,i)
c.... Temperature gradient (grad theta)
      g(15) = g(15) + ul(5,i)*shp(1,i)
      g(16) = g(16) + ul(5,i)*shp(2,i)
      g(17) = g(17) + ul(5,i)*shp(3,i)
c.... temperature
      g(18) = g(18) + ul(5,i)*shp(4,i)
c.... temperature at tn
      g(19) =  g(19) + (ul(5,i)-ul(5,i+nel))*shp(4,i)
  310 continue
c
c       --------------------------------------------------
c....   compute stresses and moduli in material subroutine
c       --------------------------------------------------
      call mate01(d,g,sh,aah,n,l,isw)

c.... move local stresses for plotting
      if(isw.eq.8) then
      sh(2) = sh(5) !P22
      sh(3) = sh(9) !P33
      sh(4) = sh(20) !s_max
      sh(5) = sh(23) !hist
        call slcnxd(sh,shp,dvol,p,s)
           goto 999
      endif

c....   multiply stress and tangent moduli by volume element
      do 330 i = 1,13
      sh(i) = sh(i)*dvol
      do 330 j = 1,13
  330 aah(i,j) = aah(i,j)*dvol
c
c....   coupled residual
        do 350 i = 1,nel

      i1 = (i-1)*ndf
c....   residual (displacements)
      p(i1+1) = p(i1+1)-shp(1,i)*sh(1)-shp(2,i)*sh(2)-shp(3,i)*sh(3)
      p(i1+2) = p(i1+2)-shp(1,i)*sh(4)-shp(2,i)*sh(5)-shp(3,i)*sh(6)
      p(i1+3) = p(i1+3)-shp(1,i)*sh(7)-shp(2,i)*sh(8)-shp(3,i)*sh(9)
c....   residual  (damage)
      p(i1+4)=p(i1+4) -shp(1,i)*sh(10)-shp(2,i)*sh(11)-shp(3,i)*sh(12)
     A                                -shp(4,i)*sh(13)
c....   residual  (temp)
      p(i1+5)=p(i1+5) -shp(1,i)*sh(15)-shp(2,i)*sh(16)-shp(3,i)*sh(17)
     A                                -shp(4,i)*sh(18)
      call pzero(hh,5*18)
        do 340 ii = 1,9
      hh(1,ii)=shp(1,i)*aah(ii,1)+shp(2,i)*aah(ii,2)+shp(3,i)*aah(ii,3)
      hh(2,ii)=shp(1,i)*aah(ii,4)+shp(2,i)*aah(ii,5)+shp(3,i)*aah(ii,6)
      hh(3,ii)=shp(1,i)*aah(ii,7)+shp(2,i)*aah(ii,8)+shp(3,i)*aah(ii,9)
 340  continue
        do 341 ii = 10,13
      hh(4,ii)=shp(1,i)*aah(ii,10)+shp(2,i)*aah(ii,11)
     A        +shp(3,i)*aah(ii,12)+shp(4,i)*aah(ii,13)
 341  continue
       do 342 ii = 15,18
      hh(5,ii)=shp(1,i)*aah(ii,15)+shp(2,i)*aah(ii,16)
     A        +shp(3,i)*aah(ii,17)+shp(4,i)*aah(ii,18)
 342  continue
c
	do 350 j = i,nel
        j1 = (j-1)*ndf
	do 350 ii = 1,5
      s(i1+ii,j1+1)= s(i1+ii,j1+1)+hh(ii,1)*shp(1,j)+hh(ii,2)*shp(2,j)
     &                            +hh(ii,3)*shp(3,j)
      s(i1+ii,j1+2)= s(i1+ii,j1+2)+hh(ii,4)*shp(1,j)+hh(ii,5)*shp(2,j)
     &                            +hh(ii,6)*shp(3,j)
      s(i1+ii,j1+3)= s(i1+ii,j1+3)+hh(ii,7)*shp(1,j)+hh(ii,8)*shp(2,j)
     &                            +hh(ii,9)*shp(3,j)
c
      s(i1+ii,j1+4)= s(i1+ii,j1+4)+hh(ii,10)*shp(1,j)+hh(ii,11)*shp(2,j)
     &                            +hh(ii,12)*shp(3,j)+hh(ii,13)*shp(4,j)
      s(i1+ii,j1+5)= s(i1+ii,j1+5)+hh(ii,15)*shp(1,j)+hh(ii,16)*shp(2,j)
     &                            +hh(ii,17)*shp(3,j)+hh(ii,18)*shp(4,j)
  350   continue
c
  300   continue

c.... form lower part of tangent by symmetry
      do 360 i = 1,nst-1
      do 360 j = i+1,nst
  360 s(j,i) = s(i,j)

  999 continue
      return

c.... formats for input-output
 2000 format(
     * 10x,'ELEMENT MODEL FOR CMP ------ (c) C. Miehe, Stuttgart'/
     * 10x,'[  ] [     ] 3-/6-node-triangle or 4-node-quad',/
     * 10x,'[  ] [iel  ] 2F-TRIQUAD. Small Strain Vers 24/03/09 ',i12/
     * 10x,'[G1] [nd0  ] general parameters ....................',i12/
     * 10x,'[G2] [ndele] element parameters ....................',i12/
     * 10x,'[G3] [nhele] element history variables .............',i12)
      end

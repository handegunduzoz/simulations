C-----------------------------------------------------------------------
      subroutine umacr0(lct,ctl,prt)
      implicit  none
c-----------------------------------------------------------------------
c
c.... Write tecplot output files, check that folder "plots" exists!
c
c.... Usage: tecp,,ctl(1)
c     tecp,,-1 - initialization and write grid file
c     tecp     - write solution file
c
c-----------------------------------------------------------------------
      include 'umac1.h'  ! uct
      include 'pointer.h'
      include 'comblk.h'
      include 'sdata.h'
      include 'cdata.h'
      include 'iofile.h'
      include 'pdata3.h'
      include 'prstrs.h'
      include 'strnum.h'
      include 'tdata.h'
      character lct*15
      character ugplotname*18
      logical   setvar, palloc
      logical   pcomp,prt
      integer   icount,itec
      real*8    ctl(3)
      common    /solcount/ icount
      save
c
c.... set command word
      if (pcomp(uct,'mac0',4)) then      ! Usual    form
        uct = 'tecp'                     ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation
c
c....   project Gauss-Values to Nodes
        setvar = palloc( 58,'NDNP',numnp*npstr,2)
        setvar = palloc( 57,'NDER',numnp*8    ,2)
        setvar = palloc( 60,'NDNS',max(nen*npstr,nst*nst),2)
        setvar = palloc(207,'NSCR',numel      ,2)
        nph    = np(58)
        ner    = np(57)
        call pjstrs(.false.)
c
c....   INITIALIZATION AND GRID FILE (tecp,,-1)
        if(ctl(1).lt.-0.5d0) then
c....     initialize counter for solutiontime command in tec360
          icount = 0
c....     write grid file
          call writegrid(numnp,numel,mr(np(33)))
          return
        endif
c
c....   SOLUTION FILE (tecp)
        icount = icount + 1
c....   determine current filename
        ugplotname = 'plots/sol00000.dat'
        if (icount.lt.10) then
          write(ugplotname(14:14),'(i1)') icount
        elseif (icount.lt.100) then
          write(ugplotname(13:14),'(i2)') icount
        elseif (icount.lt.1000) then
          write(ugplotname(12:14),'(i3)') icount
        elseif (icount.lt.10000) then
          write(ugplotname(11:14),'(i4)') icount
        elseif (icount.lt.100000) then
          write(ugplotname(10:14),'(i5)') icount
        elseif (icount.gt.99999) then
          write(*,*) ' ->ERROR: Output limited to 999999 files!!!'
        endif
c....   open solution file
        itec = 17
        open(itec,file=ugplotname,status='unknown')
c....   write tecplot header
        call tpheader(numnp,numel,itec)
c....   write nodal information
        call tpnode(hr(np(43)),hr(np(40)),hr(np(58)+numnp),
     A   	    ndm,numnp,ndf,itec)
c....   close file
         write(11,*) 'Sol.Nr.:',icount,'time',ttim
        write(*,*) 'writing tecplot solution file: ',ugplotname(7:14)
        close (itec)
      endif
c
      end
c
c
c ----------------------------------------------------------------------
      subroutine tpheader(numnp,numel,itec)
c ----------------------------------------------------------------------
c.... Ausgabe des Headers fuer Tecplot
      implicit double precision (a-h,o-z)
      include 'eldata.h'
      include 'tdata.h'
      common /solcount/ icount
      save
c
c.... for 3d brick-elements (8-noded)
      if (nel.eq.8) then
        write(itec,1001)
        write(itec,1002)
        write(itec,1003)
        write(itec,1004) numnp,numel,icount
c.... for 2d triangular elements (3- or 6-noded)
      elseif (nel.eq.3.or.nel.eq.6) then
        write(itec,1001)
        write(itec,1002)
        write(itec,1005)
        write(itec,1006) numnp,numel,icount
c.... for 2d quadrilateral elements (4-noded)
      elseif (nel.eq.4) then
        write(itec,1001)
        write(itec,1002)
        write(itec,1008)
        write(itec,1007) numnp,numel,icount
      endif
c
 1001 format('TITLE = "tecplot solution file by FEAP"')
 1002 format('FILETYPE = SOLUTION')
c.... 3d brick-elements
 1003 format('VARIABLES = "xs","ys","zs",',
     A                   '"ux","uy","uz","d", "temp",',
     A   '"s1","s2","s3","smax","hist"')
 1004 format('ZONE N=',I5,', E=',I5,', F=FEPOINT, ET=BRICK '/
     A       'SOLUTIONTIME=',I5)
c.... 2d triangular elements
 1005 format('VARIABLES = "xs","ys",',
     A                   '"ux","uy","-pot","p1","p2","s1","s2","s3",',
     A                   '"s4","s5","etotx","etoty","gradP","ex","ey"')
 1006 format('ZONE N=',I5,', E=',I5,', F=FEPOINT, ET=Triangle '/
     A       'SOLUTIONTIME=',I5)
c.... 2d quadrilateral elements
 1008 format('VARIABLES = "xs","ys",',
     A'"ux","uy","d","temp","s1","s2","s3",')
 1007 format('ZONE N=',I5,', E=',I5,', F=FEPOINT, ET=Quadrilateral '/
     A       'SOLUTIONTIME=',I5)
c
      return
      end
c
c ----------------------------------------------------------------------
      subroutine tpnode(x,u,s,ndm,numnp,ndf,itec)
c ----------------------------------------------------------------------
c.... Ausgabe der Koordinaten,Verschiebungen, Spannungen
      implicit double precision (a-h,o-z)
      include 'pdata3.h'
      include 'eldata.h'
      dimension x(ndm,numnp),xs(ndm,numnp),u(ndf,numnp),s(numnp,npstr)
c
c.... for 3d brick-elements (8-noded)
      if (nel.eq.8) then
        do i = 1,numnp
          do j = 1,ndm
          xs(j,i) = x(j,i) + u(j,i)
          enddo
        write(itec,'(13e12.4)')(xs(j,i),j=1,ndm),
     A                         ( u(j,i),j=1,ndf),
     A    s(i,1),s(i,2),s(i,3),s(i,4),s(i,5)
         enddo
      endif
c.... for 2d triangular elements (3- or 6-noded)
      if (nel.eq.3.or.nel.eq.6) then
        do i = 1,numnp
          do j = 1,ndm
          xs(j,i) = x(j,i) + u(j,i)
          enddo
        write(itec,'(19e12.4)')(xs(j,i),j=1,ndm),
     A                         ( u(j,i),j=1,ndf)
        enddo
      endif
c
c.... for 2d quadrilateral elements (4-noded)
      if (nel.eq.4) then
        do i = 1,numnp
          do j = 1,ndm
          xs(j,i) = x(j,i) + u(j,i)
          enddo
        write(itec,'(9e12.4)')(xs(j,i),j=1,ndm),
     A                         ( u(j,i),j=1,ndf),
     A s(i,1),s(i,2),s(i,3)
        enddo
      endif
c
      return
      end
c
c
c ----------------------------------------------------------------------
      subroutine writegrid(numnp,numel,ix)
      implicit none
c-----------------------------------------------------------------------
c.... Output of grid file for tecplot
      include 'sdata.h'   ! ndm,nen1
      include 'eldata.h'  ! nel
      include 'comblk.h'
      include 'pointer.h'
      integer i,j,numnp,numel,ix(nen1,numel)
      save
c
c.... INITIALIZATION AND TITLE
c.... open grid file
      open(2,file='plots/grid.dat',status='unknown')
c.... writing title etc. for grid file
      write(2,*) 'TITLE="tecplot grid file created by FEAP"'
      write(2,*) 'FILETYPE = GRID'
      if (ndm.eq.2) write(2,*) 'VARIABLES="xm","ym"'
      if (ndm.eq.3) write(2,*) 'VARIABLES="xm","ym","zm"'
c
c.... WRITE ELEMENT HEADER
c....   for 2d linear triangular elements (3-noded)
        if (nel.eq.3) write(2,1005) numnp,numel
c....   for 2d quadratic triangular elements (6-noded)
        if (nel.eq.6) write(2,1005) numnp,4*numel
c....   for 2d quadrilateral elements (4-noded)
        if (nel.eq.4) write(2,1006) numnp,numel
c....   for 3d brick-elements (8-noded)
        if (nel.eq.8) write(2,1003) numnp,numel
c
c.... WRITE NODAL COORDINATES
      if (ndm.eq.2) then
        do i = 1,numnp
        write(2,'(2e12.4)') hr(np(43)+(i-1)*ndm),
     &                      hr(np(43)+(i-1)*ndm+1)
        enddo
      elseif (ndm.eq.3) then
        do i = 1,numnp
        write(2,'(3e12.4)') hr(np(43)+(i-1)*ndm),
     &                      hr(np(43)+(i-1)*ndm+1),
     &                      hr(np(43)+(i-1)*ndm+2)
        enddo
      endif
c
c.... WRITE ELEMENT TOPOLOGY
c.... for 3d brick elements (8-noded)
      if (nel.eq.8) then
        do i = 1,numel
        write(2,'(8i8)') (ix(j,i),j=1,8)
        enddo
c.... for 2d triangular elements (3- or 6-noded)
      elseif (nel.eq.3.or.nel.eq.6) then
        do i = 1,numel
        write(2,'(3i8)') (ix(j,i),j=1,3)
        enddo
c.... for 2d quadrilateral elements (4-noded)
      elseif (nel.eq.4) then
        do i = 1,numel
        write(2,'(4i8)') (ix(j,i),j=1,4)
        enddo
      endif
c.... close grid file
      close(2)
c
c.... FORMAT STATEMENTS
c.... 3d brick elements
 1003 format('ZONE N=',I9,', E=',I9,', F=FEPOINT,',' ET=BRICK')
c.... 2d triangle-elements
 1005 format('ZONE N=',I9,', E=',I9,', F=FEPOINT,',' ET=TRIANGLE')
c.... 2d quadrilateral elements
 1006 format('ZONE N=',I9,', E=',I9,', F=FEPOINT,',' ET=QUADRILATERAL')
c
      end

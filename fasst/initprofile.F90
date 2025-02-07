      subroutine initprofile(ttest,nt0,mtest,nm0,str_flag,nstr_flag)

      use fasst_global

! This places the nodes at given locations, then determines the dt needed for 
! convergence. It also initializes all of the soil property profiles.

! calls the following subroutines:
!     sort, locate (appended to this subroutine)
!     thcond

! uses the functions: dense,head,soilhumid,vap_press

      implicit none

      character(len=1),intent(in):: ttest
      integer(kind=4),intent(inout):: nt0,nm0,str_flag(maxl)
      character(len=1),intent(inout):: mtest
      integer(kind=4),intent(out):: nstr_flag(maxn)

! local variables
      integer(kind=4):: i,j,k,l,jj,nnodesi,nmtest,nttest,fr_test,tcount
      integer(kind=4):: nnodesf,icount,nnodesi2,rocount,lcount,ntemp
      integer(kind=4):: deltaz,d1i,d2i,io,sid,ics(30)
      real(kind=8):: test_thick,delzi,delzz,botz,w,mindeltat,deltatm,f1
      real(kind=8):: w1,c1,logalpha,logn,rh,dense,delz,p,a,b,t1,vpress
      real(kind=8):: sumthick,totthick,aslope,head,htot,newsnow,modz
      real(kind=8):: vp,hatm,pres,rhoa,stemp,deltatt,mindeltatm,rhotot
      real(kind=8):: dens,pors,ssemis,ssalb,shc,smin,smax,salpha,svgn
      real(kind=8):: sspheat,sorgan,spsand,spsilt,spclay,temp,soilhumid
      real(kind=8):: rstep,rdeltaz,rnnodes,deep,smmax,smb,smt,mixr,rhow
      real(kind=8):: nz1(maxn),nz3(maxn),nz4(maxn),smi(maxn),vap_press
      real(kind=8):: rhoda(maxn),rhov(maxn),nzt(maxn)
      real(kind=8):: nz2(maxl),sph(maxl)
      character(len=1):: rotest
      character(len=2):: ssname
      character(len=200):: header

      real(kind=8),parameter:: mindepth = 2d0                            !minimum allowed soil depth
      real(kind=8),parameter:: wltpoint = -1.5d4                         !wilting point (cm)
      real(kind=8),parameter:: fldcap   = -3.4d2                         !field capacity (cm)
      real(kind=8),parameter:: minwat   = -3.1d4                         !minimum allowed water (cm)
      real(kind=8),parameter:: airdry   = -1.0d6                         !minimum head due to air drying (cm)
! Fredlund et al (2002): soil_moist(i) = 0d0 -> phead(i) ~ -10^9 Pa -> -10^7 cm absolute minimum head

!              1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16   17   18 
!              GW,GP,GM,GC,SW,SP,SM,SC,ML,CL,OL,CH,MH,OH,PT,SMSC,CLML,EV,
      data ics/2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 2,   1,   2,&
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
!                 CO,AS,         RO,   AI,      SN
!                 20 21          25    27       30
 
! zero out variables
      jj = 0
      nnodesi = 0
      nmtest = 0
      nttest = 0
      fr_test = 0
      tcount = 0
      nnodesf = 0
      icount = 0
      nnodesi2 = 0
      rocount = 0
      lcount = 0
      ntemp = 0
      d1i = 0
      d2i = 0
      deltaz = 0
      test_thick = 0d0
      delzi = 0d0
      delzz = 0d0
      botz = 0d0
      deltatm = 0d0
      f1 = 0d0
      totthick = 0d0
      w1 = 0d0
      c1 = 0d0
      logalpha = 0d0
      logn = 0d0
      rh = 0d0
      delz = 0d0
      p = 0d0
      a = 0d0
      b = 0d0
      t1 = 0d0
      vpress = 0d0
      sumthick = 0d0
      vp = 0d0
      hatm = 0d0
      pres = 0d0
      rhoa = 0d0
      stemp = 0d0
      deltatt = 0d0
      rhotot = 0d0
      rstep = 0d0
      rdeltaz = 0d0
      rnnodes = 0d0
      deep = 0d0
      smmax = 0d0
      smb = 0d0
      smt = 0d0
      mixr = 0d0
      rhow = 0d0

      do i=1,maxn
        nstr_flag(i) = 0
        nz1(i) = 0d0
        nz3(i) = 0d0
        nz4(i) = 0d0
        smi(i) = 0d0
        rhoda(i) = 0d0
        rhov(i) = 0d0
        nzt(i) = 0d0
      end do

      do i=1,maxl
        nz2(i) = 0d0
        sph(i) = 0d0
      end do

      rotest = ' '

! initialize variables
      nmtest = nm0 + 1
      nttest = nt0 + 1
      toptemp = 0d0

! initialize the layer properties
      do i=1,nlayers
        io = 0
        do k=1,20
          read(30,'(a)') header
        end do

        jj = 0
        sid = 0
        ssname = '  '
        dens = 0d0
        pors = 0d0
        ssalb = 0d0
        ssemis = 0d0
        shc = 0d0
        smin = 0d0
        smax = 0d0
        salpha = 0d0
        svgn = 0d0
        sspheat = 0d0
        sorgan = 0d0
        spsand = 0d0
        spsilt = 0d0
        spclay = 0d0

        do while(io /= -1.and.jj == 0)
          read(30,*,iostat=io) sid,ssname,dens,pors,ssemis,ssalb,shc,   &
                               smin,smax,salpha,svgn,sspheat,sorgan,    &
                               spsand,spsilt,spclay

          if(soiltype(i) == sid) then
            jj = 1
            soiltype(i) = sid
            stype(i) = ssname
            soilp(i,1) = dens
            soilp(i,2) = pors
            soilp(i,3) = ssalb
            soilp(i,4) = ssemis
            soilp(i,7) = shc
            soilp(i,8) = smin
            soilp(i,9) = smax
            soilp(i,11) = svgn
            soilp(i,13) = sspheat
            soilp(i,14) = sorgan
            soilp(i,18) = spsand
            soilp(i,19) = spsilt
            soilp(i,20) = spclay
 
            if(soilp(i,2) <= eps) soilp(i,2) = 1d-3

            if(sid <= 2.or.(sid == 5.or.sid == 6)) then
              soilp(i,23) = 2.5d0
            else if(sid == 3.or.sid == 4) then
              soilp(i,23) = 3d1
            else
              soilp(i,23) = 1d2 - soilp(i,18)
            end if

            if(dabs(rho_fac(i)-1d0) > eps) then
              soilp(i,1) = soilp(i,1)*((1d0 - soilp(i,2)/rho_fac(i))/   &
                                                    (1d0 - soilp(i,2)))
              soilp(i,2) = soilp(i,2)/rho_fac(i)
              soilp(i,7) = soilp(i,7)/rho_fac(i)
              soilp(i,9) = soilp(i,2)
              soilp(i,11) = soilp(i,11)/rho_fac(i)
              salpha = salpha*rho_fac(i)*5d-1
            end if

            if(dabs(soilp(i,5)-spflag) <= eps)                          &
                                          soilp(i,5) = 1d-2*soilp(i,18)
            if(dabs(salpha-spflag) > eps) soilp(i,10) = 1d0/salpha
            if(dabs(soilp(i,11)-spflag) > eps)                          &
                                    soilp(i,12) = 1d0 - 1d0/soilp(i,11)
            if(dabs(soilp(i,18)-spflag) > eps)                          &
                                         soilp(i,18) = 1d-2*soilp(i,18)
            if(dabs(soilp(i,19)-spflag) > eps)                          &
                                         soilp(i,19) = 1d-2*soilp(i,19)
            if(dabs(soilp(i,20)-spflag) > eps)                          &
                                         soilp(i,20) = 1d-2*soilp(i,20)
            if(dabs(soilp(i,21)-spflag) > eps)                          &
                                         soilp(i,21) = 5d-1*soilp(i,14)
            if(dabs(soilp(i,23)-spflag) > eps)                          &
                                         soilp(i,23) = 1d-2*soilp(i,23)
          end if
        end do
        rewind(30)

        if(soiltype(i) >= 19) str_flag(i) = 4

        if(soiltype(i) /= 26) then
          if(soilp(i,9) > soilp(i,2)) soilp(i,9) = soilp(i,2)
        else
          if(soilp(i,9) > soilp(i,1)) soilp(i,9) = soilp(i,1)
        end if
        
        do j=1,maxp
          f1 = anint(isoilp(i,j)*1d10)*1d-10
          if(dabs(f1-spflag) > eps) soilp(i,j) = isoilp(i,j)
        end do

        if(soilp(i,2) <= eps) soilp(i,2) = 1d-3

        if(dabs(isoilp(i,1)-spflag) <= eps.and.                         &
          dabs(isoilp(i,2)-spflag) > eps) soilp(i,1) = soilp(i,2)*2.7d0

        soilp(i,25) = dmax1(0d0,(1d0 - soilp(i,23)) - soilp(i,18))       !percent gravel

! minimum pressure head (cm)
! wilting point, field capacity and minimum allowed water values
        if((soiltype(i) /= 26.and.soiltype(i) /= 27)                     &
                                      .and.(soilp(i,10) > eps.and.       &
                                  dabs(soilp(i,10)-spflag) > eps)) then
          c1 = 1d0 + (soilp(i,10)*dabs(airdry))**soilp(i,11)
          w1 = 1d0/(c1**soilp(i,12))
          temp = soilp(i,8) + w1*(soilp(i,9) - soilp(i,8))               !minimum air-dry water content
          soilp(i,8) = dmax1(1d-3, dmin1(temp,soilp(i,8)))
          if(dabs(soilp(i,8)-soilp(i,9)) <= eps)                        &
                                           soilp(i,8) = 1d-1*soilp(i,9)

          c1 = 1d0 + (soilp(i,10)*dabs(minwat))**soilp(i,11)
          w1 = 1d0/(c1**soilp(i,12))
          soilp(i,15) = soilp(i,8) + w1*(soilp(i,9) - soilp(i,8))        !minimum allowed water content
          soilp(i,15) = dmax1(soilp(i,15),soilp(i,8)*1.001d0)

          w = (soilp(i,15) - soilp(i,8))/(soilp(i,9) - soilp(i,8))
          sph(i) = -(1d0/soilp(i,10))*                                  &
                      (((w**(-1d0/soilp(i,12)))-1d0)**(1d0/soilp(i,11))) !minimum pressure head (cm)
          sph(i) = 1d-2*sph(i)                                           !m
          sph(i) = anint(sph(i)*1d10)*1d-10

          c1 = 1d0 + (soilp(i,10)*dabs(wltpoint))**soilp(i,11)
          w1 = 1d0/(c1**soilp(i,12))
          soilp(i,16) = soilp(i,8) + w1*(soilp(i,9) - soilp(i,8))        !wilting point
          soilp(i,16) = dmax1(soilp(i,16),soilp(i,15)*1.01d0)

          c1 = -6d-1*(2d0 + dlog10(soilp(i,7)*8.64d4))                   !Twarakavi et al., WWR,V45,W10410,doi:10.1029/2009WR007944
          w1 = soilp(i,11)**c1
          soilp(i,17) = soilp(i,8) + w1*(soilp(i,9) - soilp(i,8))        !field capacity

          c1 = 1d0 + (soilp(i,10)*dabs(fldcap))**soilp(i,11)             !Fredlund
          w1 = 1d0/(c1**soilp(i,12))
          temp = soilp(i,8) + w1*(soilp(i,9) - soilp(i,8))               !field capacity
          soilp(i,17) = dmax1(7.5d-1*soilp(i,9),soilp(i,17),temp)
          soilp(i,24) = 0.999d0*soilp(i,9)
        else
          sph(i) = -1d-5                                                 !based on Fredlund
          soilp(i,15) = soilp(i,8)                                       !minimum allowed water content
          soilp(i,16) = soilp(i,15)                                      !wilting point
          soilp(i,17) = soilp(i,9)                                       !field capacity
          soilp(i,24) = soilp(i,9)
        end if

        if(dabs(soilp(i,8)-soilp(i,9)) <= eps) then
          write(*,'('' minimum and maximum water contents are equal'')')
          stop
        end if

! thermal conductivity of dry soil within +/- 20% (W/m*K)
        if(soiltype(i) /= 26.and.soiltype(i) /= 27) then
          if(dabs(soilp(i,6)-spflag) <= eps) then
            if(soiltype(i) <= 4) then
              soilp(i,6) = 0.039d0*(soilp(i,2)**(-2.2d0))
              soilp(i,6) = -0.56d0*soilp(i,2) + 0.51d0
            else if(soiltype(i) > 4.and.soiltype(i) <= 18) then
              f1 = 1d0/(1d0 - soilp(i,2))
              soilp(i,6) = (0.135d0*soilp(i,1)*1d3 + 64.7d0)/           &
                                        (soilp(i,1)*1d3*(f1 - 0.947d0))
!              soilp(i,6) = -5.6d-1*nsoilp(i,2) + 5.1d-1                  !Lu et al (2007)          
              if(soiltype(i) /= 15) then
                soilp(i,6) = -0.56d0*soilp(i,2) + 0.51d0
              else
                soilp(i,6) = 2d0*soilp(i,6)
              end if
            end if
          end if
        end if

        do j=1,maxp
          soilp(i,j) = anint(soilp(i,j)*1d10)*1d-10
        end do
      end do

! determine if there is snow/ice on the ground and set the albedo and emissivity
      sgralbedo = soilp(1,3)
      sgremis = soilp(1,4)
      albedo = 1d0    !initialize
      if((dabs(hi) <= eps.and.dabs(hsaccum) <= eps).and.                &
         (aint(met(1,ip_pt)) /= 3.or.aint(met(1,ip_pt2)) /= 3)) then     !no snow on ground
        albedo = dmin1(albedo,sgralbedo)
        emis = sgremis
        hsaccum = 0d0
        hi = 0d0
        newsnow = 0d0  
      else if(hsaccum > eps) then                                        !snow
        emis = semis
        if(aint(met(1,ip_pt)) /= 3.or.aint(met(1,ip_pt2)) /= 3) then     !old
          albedo = dmin1(albedo,soalbedo)
          newsnow = 0d0
        else                                                             !new
          albedo = dmin1(albedo,snalbedo)
          newsnow = (met(1,ip_prec) + met(1,ip_prec2))*1d-3
        end if
      else if(hi > eps.and.dabs(hsaccum) <= eps) then                    !ice
        albedo = dmin1(albedo,ialbedo)
        emis = iemis
        hsaccum = 0d0
        newsnow = 0d0
      end if
      htot = hsaccum + hi + newsnow

! determine node placement based on # of layers, yes/no to temp and/or moist meas.
      do i=1,nlayers
        totthick = totthick + lthick(i)
        nz2(i) = totthick
      end do
      sumthick = totthick
      if(totthick < mindepth) totthick = mindepth

      if((ttest == 'y'.or.ttest == 'Y').and.                            &
               (mtest == 'y'.or.mtest == 'Y')) then
        nnodesi = nt0
        if(dabs(zti(1)) > eps) then
          do j=2,nnodesi+1
            nz1(j) = zti(j-1)
          end do
          nnodesi = nt0 + 1
        else
          do j=1,nnodesi
            nz1(j) = zti(j)
          end do
          nnodesi = nt0
        end if

        nnodesi2 = nm0
        if(dabs(zm(1)) > eps) then
          do j=2,nnodesi2+1
              nz4(j) = zm(j-1)
          end do
          nnodesi2 = nm0 + 1
        else
          do j=1,nnodesi2
              nz4(j) = zm(j)
          end do
          nnodesi2 = nm0
        end if

        nnodesf = nnodesi + nnodesi2 + nlayers
        do i=1,nnodesf
          if(i <= nnodesi) then
            nz3(i) = nz1(i)
          else if(i > nnodesi.and.i <= nnodesi2+nnodesi) then
            nz3(i) = nz4(i-nnodesi)
          else
            nz3(i) = nz2(i-(nnodesi+nnodesi2))
          end if
        end do

      else if((ttest == 'y'.or.ttest == 'Y').and.                       &
               (mtest == 'n'.or.mtest == 'N')) then
        nnodesi = nt0
        if(dabs(zti(1)) > eps) then
          do j=2,nnodesi+1
            nz1(j) = zti(j-1)
          end do
          nnodesi = nt0 + 1
        else
          do j=1,nnodesi
            nz1(j) = zti(j)
          end do
          nnodesi = nt0
        end if

        nnodesf = nnodesi + nlayers
        do i=1,nnodesf
          if(i <= nnodesi) then
            nz3(i) = nz1(i)
          else
            nz3(i) = nz2(i-nnodesi)
          end if
        end do

      else if((ttest == 'n'.or.ttest == 'N').and.                       &
               (mtest == 'y'.or.mtest == 'Y')) then
        nnodesi = nm0
        if(dabs(zm(1)) > eps) then
          do j=2,nnodesi+1
              nz1(j) = zm(j-1)
          end do
          nnodesi = nm0 + 1
        else
          do j=1,nnodesi
              nz1(j) = zm(j)
          end do
          nnodesi = nm0
        end if

        nnodesf = nnodesi + nlayers
        do i=1,nnodesf
          if(i <= nnodesi) then
            nz3(i) = nz1(i)
          else
            nz3(i) = nz2(i-nnodesi)
          end if
        end do

      else if((ttest == 'n'.or.ttest == 'N').and.                       &
               (mtest == 'n'.or.mtest == 'N')) then
         nz1(2) = totthick
         nnodesi = 2

        nnodesf = nnodesi + nlayers
        do i=1,nnodesf
          if(i <= nnodesi) then
            nz3(i) = nz1(i)
          else
            nz3(i) = nz2(i-nnodesi)
          end if
        end do

      end if

      call sort(nnodesf,nz3,icount,nz)
      nnodesi = icount

! add more nodes if spacing is too large
      nz1(1) = nz(1)
      icount = 1
      if(water_flag == 0) then                                           !not water or air
        do while(nz1(icount) < nz(nnodesi)) !nz(nnodesi))  !1d0)
          icount = icount + 1
          if(nz1(icount-1) <= 1d0) then
            nz1(icount) = nz1(icount-1) + 1.5d-2*(2d0*                  &
                                                   float(icount) - 3d0)
          else
            nz1(icount) = nz1(icount-1) + 3d-2*(2d0*                    &
                                                   float(icount) - 3d0)
          end if
        end do
        if(nz1(icount) > nz(nnodesi)) icount = icount - 1  !1d0

        do i=2,icount
          do k=2,nnodesi
            if(dabs(nz1(i)-nz(k)) <= 0.5d0*(1.5d-2*(2d0*float(i)-3d0))) then
              nz1(i) = nz(k)
            end if
          end do
        end do

        do i=2,nnodesi
          icount = icount + 1
          nz1(icount) = nz(i)
        end do

        call sort(icount,nz1,nnodesi,nz)

      else if(water_flag == 1.or.water_flag == 2) then                  !water
        do while(nz1(icount) < nz(nnodesi))
          icount = icount + 1
          if(nz1(icount-1) <= 1d0) then
            nz1(icount) = nz1(icount-1) + 1.5d-2*(2d0*                  &
                                                   float(icount) - 3d0)
          else
            nz1(icount) = nz1(icount-1) + 3d-2*(2d0*                    &
                                                   float(icount) - 3d0)
          end if
        end do
        if(nz1(icount) > nz(nnodesi)) icount = icount - 1  !1d0

        icount = icount + 1
        if(vegh_type > 0) then
          nz1(icount) = 3d0
        else if(vegl_type > 0) then
          nz1(icount) = 1.5d0
        else
          nz1(icount) = 10d0
        end if
        icount = icount + 1
        nz1(icount) = nz1(icount-1) + 1d0

        do i=2,icount
          do k=2,nnodesi
            if(dabs(nz1(i)-nz(k)) <= 1.5d-2*(2d0*float(i)-3d0)) then
              nz1(i) = nz(k)
            end if
          end do
        end do

        do i=2,nnodesi
          icount = icount + 1
          nz1(icount) = nz(i)
        end do

        call sort(icount,nz1,nnodesi,nz)
      end if

! add more nodes if total thickness < 1m
      if(nz(nnodesi)-totthick < 0d0) then
        delzz = dmin1(3d-1,2d0*(nz(nnodesi) - nz(nnodesi-1)))
        modz = (totthick - nz(nnodesi))                                 &
                         - anint((totthick - nz(nnodesi))/delzz)*delzz
!        if(dabs(dmod(totthick - nz(nnodesi),delzz)) <= eps) then
        if(dabs(modz) < 1d0) then
          rnnodes = (totthick - nz(nnodesi))/delzz
          nnodes = nnodesi + int(rnnodes) + 1
        else
          rnnodes = (totthick - nz(nnodesi))/delzz
          nnodes = nnodesi + int(rnnodes) + 1
        end if

        do j=1,nnodesi
          nz1(j) = nz(j)
        end do

        do j=nnodesi+1,nnodes
          nz1(j) = nz1(j-1) + delzz
        end do

        do while(nz1(nnodes) > totthick) 
          nnodes = nnodes - 1
        end do

        do while(nz1(nnodes) < totthick) 
          nnodes = nnodes + 1
          nz1(nnodes) = nz1(nnodes-1) + delzz
        end do

        call sort(nnodes,nz1,nnodesf,nz)
        nnodes = nnodesf
      else
        nnodes = nnodesi
      end if

      ntot = nnodes + 2
      ntemp = nnodes
      if(veg_flagl == 0) then
        if(hsaccum+hi > 0d0) then
          ntemp = nnodes + 1
          node_type(ntemp) = 'HM'
        end if
      else
        if(dabs(sigfl) <= eps) then
          if(hsaccum+hi > 0d0) then
            ntemp = nnodes + 1
            node_type(ntemp) = 'HM'
          end if
        else
          ntemp = nnodes + 1
          node_type(ntemp) = 'HM'
          if(hsaccum+hi > eps.and.(ihfol > hsaccum+hi                   &
                                   .and.dabs(ihfol-spflag) > eps)) then
            ntemp = nnodes + 2
            node_type(ntemp) = 'VG'
          end if
        end if 
      end if

      totthick = nz(nnodes)
      test_thick = totthick

! reverse nodes so that 1 is at the bottom and nnodes at the surface
      do i=1,nnodes
        j = nnodes - i + 1
        nzt(i) = nz(j)
      end do
      if(nz(nnodes) < eps) nz(nnodes) = 0d0

      refn = 0
      do i=1,nnodes
        nz(i) = anint(nzt(i)*1d5)*1d-5
        if(nz(i) > 1d-1.and.nz(i) < 3d-1) refn = i
      end do

      if(sumthick < totthick)                                           &
        lthick(nlayers) = lthick(nlayers) + (totthick - sumthick)

! initialize soil type at a node
      delzi = totthick
      rotest = 'n'
      do i=nlayers,1,-1
        botz = anint((delzi - lthick(i))*1d5)*1d-5
        do j=1,nnodes
          if(j == 1.and.i == nlayers) then
            node_type(j) = stype(nlayers)
            ntype(j) = soiltype(nlayers)
            do k=1,maxp
              nsoilp(j,k) = soilp(nlayers,k)
            end do
            pheadmin(j) = sph(nlayers)
            nstr_flag(j) = str_flag(nlayers)
            if(ntype(j) >= 19) rotest = 'y'
            icourse(j) = ics(ntype(j))
!          else if(nz(j) > botz.and.nz(j) <= delzi) then
          else if(nz(j)-botz > eps.and.nz(j)-delzi <= eps) then
            node_type(j) = stype(i)
            ntype(j) = soiltype(i)
            do k=1,maxp
              nsoilp(j,k) = soilp(i,k)
            end do
            pheadmin(j) = sph(i)
            nstr_flag(j) = str_flag(i)
            if(ntype(j) >= 19) rotest = 'y'
            icourse(j) = ics(ntype(j))
          else if(j == nnodes.or.nz(j) < botz) then
            node_type(j) = stype(1)
            ntype(j) = soiltype(1)
            do k=1,maxp
              nsoilp(j,k) = soilp(1,k)
            end do
            pheadmin(j) = sph(1)
            nstr_flag(j) = str_flag(1)
            if(ntype(j) >= 19) rotest = 'y'
            icourse(j) = ics(ntype(1))
          end if
        end do

        do l=1,nm0
          if(zm(l) >= botz.and.zm(l) <= delzi) then
            if(str_flag(i) <= 2) sm(l) = sm(l)*soilp(i,9)
            if(sm(l) < 1.1d0*soilp(i,8)) sm(l) = 1.1d0*soilp(i,8)
            if(sm(l) > 0.99d0*soilp(i,9)) sm(l) = 0.99d0*soilp(i,9)
          end if
        end do

        delzi = botz
      end do

! initialize temperature and head above the "soil"; determine location of air nodes
      sid = nnodes
      do i=1,nnodes
        if(0.9d0-nz(i) <= eps) sid = i
      end do

      pres = met(1,ip_ap)                                                !mbar
      stt(nnodes+1) = met(1,ip_tmp) + Tref
      vp = vap_press(nnodes+1,1d-2*met(1,ip_rh),pres)                    !Pa (vapor pressure)
      rhoa = (pres*1d2 - vp)/(Rd*(met(1,ip_tmp) + Tref))
      hatm = -vp/(rhoa*grav)                                             !head due to vapor pressure (m)
      hatm = anint(hatm*1d10)*1d-10

      if(ntemp == nnodes) then
        stt(nnodes+1) = met(1,ip_tmp) + Tref
        if(aint(dabs(met(1,ip_tsoil)-mflag)*1d5)*1d-5 > eps              &
                     .and.mstflag == 1) stt(nnodes+1) = met(1,ip_tsoil)
        stt(nnodes+2) = stt(nnodes+1)
        ftemp = stt(nnodes+2)

        phead(nnodes+1) = hatm
        phead(nnodes+2) = hatm
      else if(ntemp == nnodes+1.and.(dabs(ihfol) <= eps.and.            &
                                           hsaccum+hi > eps)) then       !no foliage, yes snow/ice
        stt(nnodes+1) = met(1,ip_tmp) + Tref
        stt(nnodes+2) = dmin1(Tref,met(1,ip_tmp)+Tref)
        ftemp = stt(nnodes+2)

        phead(nnodes+1) = hatm
        phead(nnodes+2) = 0d0                                            !assumes snow is saturated
      else if(ntemp == nnodes+1.and.(ihfol > eps.and.                   &
                                    dabs(hsaccum+hi) <= eps)) then       !yes foliage, no snow/ice
        stt(nnodes+1) = met(1,ip_tmp) + Tref
        stt(nnodes+2) = met(1,ip_tmp) + Tref
        ftemp = stt(nnodes+1)

        phead(nnodes+1) = hatm
        phead(nnodes+2) = hatm
      else if(ntemp == nnodes+1.and.(ihfol > eps.and.                   &
                                           hsaccum+hi > eps)) then       !foliage buried by snow/ice
       stt(nnodes+1) = dmin1(Tref,met(1,ip_tmp)+Tref)
       stt(nnodes+2) = dmin1(Tref,met(1,ip_tmp)+Tref)
       ftemp = stt(nnodes+1)

        phead(nnodes+1) = 0d0                                            !assumes snow is saturated
        phead(nnodes+2) = 0d0                                            !assumes snow is saturated
      else if(ntemp == nnodes+2) then                                    !yes foliage, yes snow/ice
        stt(nnodes+1) = dmin1(Tref,met(1,ip_tmp)+Tref)
        stt(nnodes+2) = met(1,ip_tmp) + Tref
        ftemp = stt(nnodes+2)

        phead(nnodes+1) = hatm
        phead(nnodes+2) = 0d0                                            !assumes snow is saturated
      end if    

! initialize soil temperature
      if(ttest == 'y'.or.ttest == 'Y') then
        if(zti(nt0) < test_thick) then
          nt0 = nttest
          zti(nt0) = test_thick
          tm(nt0) = tm(nt0-1)
        end if

        if(dabs(zti(1)) > eps) then
          stt(nnodes) = ((met(1,ip_tmp) + Tref)*iheight                 &
                         + nz(nnodes-1)*tm(1))/(nz(nnodes-1) + iheight)
          if(nz(nnodes-1) < 1d-2) stt(nnodes) = tm(1)
!         stt(nnodes) = met(1,ip_tmp) + Tref
        else
          stt(nnodes) = tm(1)
        end if
        if(stt(nnodes) <= Tref) fr_test = 1
        do i = 1,nnodes-1
          call locate(zti,nt0,nz(i),jj)
          stt(i) = ((nz(i)-zti(jj))/(zti(jj)-zti(jj+1)))*               &
                                             (tm(jj)-tm(jj+1)) + tm(jj)
          if(stt(i) <= Tref.and.ntype(i) /= 27) fr_test = 1
          if(node_type(i) == 'AI') stt(i) = met(1,ip_tmp) + Tref
        end do
      else
        if(aint(dabs(met(1,ip_tsoil)-mflag)*1d5)*1d-5 > eps)  then
          stt(nnodes) = met(1,ip_tsoil)
        else if(dabs(htot) <= eps.and.node_type(nnodes) /= 'SN') then
          stt(nnodes) = met(1,ip_tmp) + Tref                             !no snow
        else if(dabs(htot) > eps.or.node_type(nnodes) == 'SN') then
          stt(nnodes) = dmin1(Tref,met(1,ip_tmp)+Tref)                   !snow/ice on grd
        end if
        if(stt(nnodes) <= Tref.and.ntype(nnodes) /= 27) fr_test = 1
        stt(sid) = stt(nnodes)*(1d0 + 2.5d-2*nz(sid))
          if(htot > eps) stt(sid) = stt(nnodes)*(1d0 + 1.5d-2*nz(sid))

        do i=1,nnodes-1
          stt(i) = stt(nnodes)*(1d0 + 2.5d-2*nz(i))
          if(htot > eps) stt(i) = stt(nnodes)*(1d0 + 1.5d-2*nz(i))
          if(nz(i) >= nz(sid)) stt(i) = stt(sid)

          if(node_type(i) == 'SN')                                      &
            stt(i) = dmin1(Tref,met(1,ip_tmp)+Tref)
          if(node_type(i) == 'AI') stt(i) = met(1,ip_tmp) + Tref
          if(stt(i) <= Tref.and.ntype(i) /= 27) fr_test = 1
        end do
      end if

      if(dabs(hsaccum+hi) <= eps.and.dabs(zti(1)) <= eps) then
        toptemp = stt(nnodes)
      else
        toptemp = dmin1(Tref,met(1,ip_tmp)+Tref)
      end if

! initialize soil moisture
!     determine if a buried rock layer is present
      i = nnodes
      do while (i >= 1)
        if(ntype(i) >= 19) then
          tcount = i                                                     !node of top of rock or air layer
          lcount = max(lcount,tcount)
        end if
        i = i - 1
      end do
      tcount = lcount

      lcount = 0
      i = 1
      do while (i < nnodes)
        if(ntype(i) < 19.and.ntype(i+1) >= 19) lcount = i                !node of bottom of rock layer
        i = i + 1
      end do
      if(lcount /= 0) then
        lcount = lcount + 1
      else if(lcount == 0) then
        if(tcount /= 0) then
          lcount = 1
        else
          lcount = tcount
        end if
      end if

!     see if more bured layers; extent of rock layer
      rocount = 0                                                        !number of impermeable layers
      if(rotest == 'y') then
        do i=1,nnodes
          if(ntype(i) >= 19) rocount = rocount + 1
        end do
      end if
      if(rocount == nnodes.or.(water_flag == 1.or.water_flag == 2))     &
        mtest = '3'

      if(tcount == 0) then
        if(rocount == nnodes) then
          tcount = nnodes
        else
          tcount = 1
        end if
      end if      

! determine position of groundwater table if not given
! based on Beven, K., R. Lamb, P. Quinn, R. Romanowicz, J. Freer (1995)
! in Singh, V.P. (Ed.), Computer Models of Watershed Hydrology. WRR Publications, PP. 627-668.
      if(gwl < 0d0) then
        i = max(1,lcount-1)
        if(slope > eps.and.slope < 90d0) then
          aslope = dmax1(0d0,dmin1(1.57d0,slope*pi/1.8d2))
          gwl = 1d1*(1d-1/(nsoilp(i,9) - nsoilp(i,8)))                  &
                   *dlog(dabs((-pheadmin(i)/totthick + dcos(aslope))/   &
                                              (totthick*dtan(aslope))))
          if(gwl < eps) gwl = 1d1
        else
          gwl = 1d1*(1d-1/(nsoilp(i,9) - nsoilp(i,8)))                  &
                                     *dlog(dabs(-pheadmin(i)/totthick)) 
          if(gwl < eps) gwl = 1d1
        end if
        if(vegl_type == 8.or.vegl_type == 11) gwl = 5d1
      end if

      if(mtest == 'y'.or.mtest == 'Y') then                              !some measured moistures
        if(tcount /= nnodes) then
          if(zm(nm0) < nz(tcount)) then
            nm0 = nmtest
            zm(nm0) = test_thick
            sm(nm0) = (nsoilp(tcount,8) + nsoilp(tcount,9))*5d-1
            if(zm(nm0) >= gwl) sm(nm0) = nsoilp(tcount,9)
          end if
        end if

        stemp = 0d0
        do i=1,nnodes
          smmax = nsoilp(i,9) - nsoilp(i,8)

          call locate(zm,nm0,nz(i),jj)

          if(dabs(zm(jj)-zm(jj+1)) > eps.and.zm(jj+1) > zm(jj)) then
            soil_moist(i) = ((nz(i) - zm(jj))/(zm(jj) - zm(jj+1)))*     &
                                (sm(jj) - sm(jj+1)) + sm(jj)

            if(soil_moist(i) < nsoilp(i,8))                             &
                                      soil_moist(i) = 1.1d0*nsoilp(i,8)
            if(nz(i) >= gwl) soil_moist(i) = nsoilp(i,9)
            if(soil_moist(i) > nsoilp(i,9))                             & 
                                      soil_moist(i) = 0.9d0*nsoilp(i,9)
            stemp = soil_moist(i)
            if(ntype(i) > 18.and.(soil_moist(i) < nsoilp(i,8).or.       &
                                          soil_moist(i) > nsoilp(i,9))) &
              soil_moist(i) = nsoilp(i,8) + 0.3d0*smmax
          else if(dabs(zm(jj)-nz(i)) <= eps) then
            soil_moist(i) = sm(jj)
            if(soil_moist(i) < nsoilp(i,8))                             & 
                                      soil_moist(i) = 1.1d0*nsoilp(i,8)
            if(nz(i) >= gwl) soil_moist(i) = nsoilp(i,9)
            if(soil_moist(i) > nsoilp(i,9))                             & 
                                      soil_moist(i) = 0.9d0*nsoilp(i,9)
          else
            if(ntype(i) <= 18.and.dabs(stemp) > eps) then
              soil_moist(i) = stemp + sm(nm0)*smmax
            else
              soil_moist(i) = nsoilp(i,8) + 0.3d0*smmax
            end if
          end if

          if(ntype(i) == 26) soil_moist(i) = 1d0
          if(nz(i) >= gwl.or.((water_flag == 1.or.water_flag == 2)      &
                       .and.ntype(i) /= 26))soil_moist(i) = nsoilp(i,9)

          if(ntype(i) == 27) soil_moist(i) = 1d-2*met(1,ip_rh)

          if(soil_moist(i) < nsoilp(i,8).or.                            &
                                      soil_moist(i) > nsoilp(i,9)) then
            if(single_multi_flag == 0) write(*,'('' Initial soil '',    &
              &''moistures out of range for user soil type. Fix input'',&
              &'' file.'')')
            stop
          end if
        end do
      else                                                               !no measured moistures
        do i=1,nnodes
          smmax = nsoilp(i,9) - nsoilp(i,8)

          if(node_type(i) /= 'SN') then
            soil_moist(i) = nsoilp(i,8) + 0.3d0*smmax
            if(ntype(i) == 27) soil_moist(i) = 1d-2*met(1,ip_rh)
          else
            soil_moist(i) = nsoilp(i,8) + 0.98d0*smmax
          end if
          if(ntype(i) == 26) soil_moist(i) = 1d0
          if(nz(i) >= gwl.or.((water_flag == 1.or.water_flag == 2)      &
                       .and.ntype(i) /= 26))soil_moist(i) = nsoilp(i,9)
        end do
      end if

! determine the head, adjust ksat for depth
! initialize total moisture and ice content (fractions), node state
      ice(nnodes+1) = 0d0
      ice(nnodes+2) = 0d0
      do i=1,nnodes
        ice(i) = 0d0
        if(i == nnodes.and.(ntype(i) == 26.and.hm > eps)) then
          ice(i) = 1d0
          soil_moist(i) = 0d0
        end if

        soil_moist(i) = anint(soil_moist(i)*1d10)*1d-10

        if(ntype(i) < 19) then
          deep = 1d0
          if(elev-nz(i) < 1d0) then
            if(dabs(2d-2*(nz(i) - elev)) < 5d1)                           &
              deep = dmax1(5d-1,dexp(2d-2*(nz(i) - elev)))  !1d-1
            nsoilp(i,7) = deep*nsoilp(i,7)
          end if
!          nsoilp(i,7) = nsoilp(i,7)                                     &
!                        *dexp((-1d-1/(nsoilp(i,9) - nsoilp(i,8)))*nz(i))
        end if

        soil_moist(i) = anint(soil_moist(i)*1d10)*1d-10
        stt(i) = anint(stt(i)*1d10)*1d-10
        ice(i) = anint(ice(i)*1d10)*1d-10

        if(ntype(i) /= 26.and.ntype(i) /= 27) then
          phead(i) = head(i,soil_moist(i))
        else
          phead(i) = 0d0
        end if
      end do

! determine the water vapor content; adjust nz(i) for elevation
      gwl = elev - gwl
      wvc(nnodes+1) = 0d0
      wvc(nnodes+2) = 0d0
      do i=1,nnodes
        if(ntype(i) /= 26.and.ntype(i) /= 27) then
          rh = soilhumid(i,phead(i),soil_moist(i),stt(i))
          d1i = 1
          p = met(1,ip_ap) + 1d-2*(nz(i) + dabs(phead(i)))              &
                                          *dense(stt(i),0d0,d1i)*grav    !mbar
        else if(ntype(i) == 26) then
          rh = 1d0
          p = met(1,ip_ap)
        else if(ntype(i) == 27) then
          rh = met(istart,ip_rh)*1d-2
          p = met(1,ip_ap)
        end if

        rhow = dense(stt(i),0d0,1)
        vpress = vap_press(i,rh,p)                                       !Pa
        mixr = 0.622d0*vpress/(p*1d2 - vpress)                           !unitless (kg/kg)
        rhov(i) = 0.622d0*vpress/(Rv*stt(i))                             !kg/m^3
        rhoda(i) = dmax1(0.95d0,dmin1(2.8d0,(p*1d2 - vpress)/           &
                                                          (Rd*stt(i))))  !kg/m^3

!        wvc(i) = (rhoda(i)/rhov(i))*mr*(nsoilp(i,2) - (soil_moist(i)    &
!                                                    + ice(i)))
!        t1 = dmin1(1d0,dmax1(mixr*rhoda(i)/(rhov(i) + mixr*rhoda(i))    &
!                                                                 ,0d0))
        t1 = dmin1(1d0,dmax1(mixr*rhoda(i)/(rhow + mixr*rhoda(i)),0d0))
        wvc(i) = t1*(nsoilp(i,2) - (soil_moist(i) + ice(i)))
        if(ntype(i) == 27) wvc(i) = t1

        wvc(i) = dmax1(0d0,dmin1(wvc(i),nsoilp(i,2) - (soil_moist(i)    &
                                                            + ice(i))))
        if(ntype(i) == 26) wvc(i) = 0d0
        wvc(i) = anint(wvc(i)*1d15)*1d-15

        nz(i) = elev - nz(i)                                             !place z-axis so that m.s.l.= 0 and z + upward from msl
        nz(i) = anint(nz(i)*1d5)*1d-5
        nzi(i) = nz(i)
      end do

! diffusion*dt/dz^2 < 0.5 for stability in finite diff calculations
      newsd = 0d0
      d1i = 1
      call th_param(d1i,maxn,0d0,0d0,rhov,rhoda,rhotot)

      mindeltat = 99999d0                                                !seconds
      mindeltatm = 99999d0                                               !seconds
      do i=1,nnodes
        if(i == 1) then
          delzs(i) = (nz(i+1) - nz(i))  !5d-1*                            !m
        else if(i == nnodes) then
          delzs(i) = (nz(i) - nz(i-1))  !5d-1*                            !m
        else
          delzs(i) = 5d-1*(nz(i) - nz(i-1)) +  5d-1*(nz(i+1) - nz(i))    !m
        end if
        delzs(i) = anint(delzs(i)*1d5)*1d-5
        delzsi(i) = delzs(i)

        if(dabs(grthcond(i)) > eps)                                     &
         deltatt = 0.45d0*delzs(i)*delzs(i)*grspheat(i)/grthcond(i)
        mindeltat = dmin1(mindeltat,deltatt)                             !s

        if(dabs(nsoilp(i,7)) > eps.and.ntype(i) < 20) then
          deltatm = 0.45d0*delzs(i)*delzs(i)*nsoilp(i,9)                &
                          /(nsoilp(i,7)*1d-2*dabs(head(i,nsoilp(i,17)))) !s
          mindeltatm = dmin1(mindeltatm,deltatm)
        end if
      end do

      deltat = mindeltat
      if(timstep >= 1d0) then
        deltat = dmin1(dmax1(deltat,3d2),timstep*3.6d3)  !1.2d2
      else 
        deltat = dmin1(dmin1(deltat,1d1),timstep*3.6d3)  !1.5d1
      end if

      rstep = timstep*3.6d3/deltat
      step = int(rstep)
      step = max(1,step)

      deltat = timstep*3.6d3/float(step)
      deltati = deltat
      stepi = step

      end subroutine initprofile

! ******************************************************************************
      subroutine sort(n,arr,ncount,b)

      use fasst_global
! this subroutine sorts the nodes in ascending order

      implicit none

      integer(kind=4),intent(in):: n
      integer(kind=4),intent(out):: ncount
      real(kind=8),intent(inout):: arr(maxn)
      real(kind=8),intent(out):: b(maxn)

! local variables
      integer(kind=4):: i,j
      real(kind=8):: a

      a = 0d0
      do i=1,maxn
        b(i) = 0d0
      end do

      do j=2,n
        a = arr(j)
        i = j - 1
        do while(i >= 1.and.arr(i) > a)
          arr(i+1) = arr(i)
          i = i - 1
        end do
        if(i == 1.and.arr(i) > a) i = 0
        arr(i+1) = a
      end do

      ncount = 1
      b(ncount) = arr(1)
      do j=2,n
        a = arr(j)*1000d0
        if(dabs(a-arr(j-1)*1d3) > 1d-4) then !eps) then
          ncount = ncount + 1
          b(ncount) = arr(j)
        end if
      end do

      end subroutine sort

! ******************************************************************************
      subroutine locate(zmi,n0,z,jj)
! this subroutine finds the nearest neighbor to a point

      use fasst_global

      implicit none

      integer(kind=4),intent(in):: n0
      integer(kind=4),intent(out):: jj
      real(kind=8),intent(in):: z,zmi(maxn)

! local variables
      integer(kind=4):: jl,jm,ju

      jj = 0
      jl = 1
      jm = 0
      ju = 0

      ju = n0 + 1
      do while (ju-jl > 1)
        jm = aint((ju+jl)*0.5)

        if((zmi(n0) > zmi(1)).eqv.(z > zmi(jm))) then
           jl = jm
        else 
           ju = jm
        end if
      end do
     
      jj = jl

      end subroutine locate

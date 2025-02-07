      subroutine sub_divide(ii,delvar,delvar1,dmet)

      use fasst_global

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: ii
      real(kind=8),intent(in):: delvar(maxcol),delvar1(12)
      real(kind=8),intent(out):: dmet(13)

! local variables
      integer(kind=4):: i
      real(kind=8):: ltest,solzen,saz,sdir,sdif,delaz,difcoeff
      real(kind=8):: sdircorr,sdifcorr1,sdifcorr2,pre_slope,f2

! zero-out variables
      ltest = 0d0
      solzen = 0d0
      saz = 0d0
      sdir = 0d0
      sdif = 0d0
      delaz = 0d0
      difcoeff = 0d0
      sdircorr = 0d0
      sdifcorr1 = 0d0
      sdifcorr2 = 0d0

      do i=1,13
        dmet(i) = 0d0
      end do

! sectorize the meterological parameters
      ltest = float(ii)  !-1
      if(ltest < 0) ltest = 0d0
      if(iw >= 2) then
        solzen = met(iw-1,ip_zen) + ltest*delvar(ip_zen)
        saz = met(iw-1,ip_az) + ltest*delvar(ip_az)

        dmet(1) = dmet1(iw-1,1) + ltest*delvar1(1)                       !total incoming solar
        dmet(2) = dmet1(iw-1,2) + ltest*delvar1(2)                       !incoming ir
        dmet(3) = dmet1(iw-1,3) + ltest*delvar1(3)                       !wind speed
        dmet(4) = dmet1(iw-1,4) + ltest*delvar1(4)                       !air temperature (K)
        dmet(5) = dmet1(iw-1,5) + ltest*delvar1(5)                       !relative humidity

        dmet(6) = delvar1(6)
        dmet(7) = delvar1(7)

        if(aint(dabs(delvar1(8)-mflag)*1d5)*1d-5 > eps) then             !reflected solar
          dmet(8) = dmet1(iw-1,8) + ltest*delvar1(8)
        else
          dmet(8) = mflag
        end if

        dmet(9) = dmet1(iw-1,9) + ltest*delvar1(9)                       !direct solar
        dmet(10) = dmet1(iw-1,10) + ltest*delvar1(10)                    !diffuse solar
        dmet(11) = dmet1(iw-1,11) + ltest*delvar1(11)                    !air pressure

        if(aint(dabs(delvar1(12)-mflag)*1d5)*1d-5 > eps) then            !upwelling ir
          dmet(12) = met(iw-1,ip_irup) + ltest*delvar1(12)
        else
          dmet(12) = mflag
        end if

        if(mstflag == 1.and.                                            &
                aint(dabs(delvar(ip_tsoil)-mflag)*1d5)*1d-5 > eps) then  !measured soil surface temp
          dmet(13) = met(iw-1,ip_tsoil) + ltest*delvar(ip_tsoil)
        else
          dmet(13) = mflag
        end if

      else
        solzen = met(iw,ip_zen)
        saz = met(iw,ip_az)

        dmet(1) = dmet1(iw,1)                                              !total incoming solar
        dmet(2) = dmet1(iw,2)                                              !incoming IR
        dmet(3) = dmet1(iw,3)                                              !wind speed
        dmet(4) = dmet1(iw,4)                                              !air temperature
        dmet(5) = dmet1(iw,5)                                              !relative humidity

!        dmet(6) = dmet1(iw,6)
!        dmet(7) = dmet1(iw,7)
        dmet(6) = delvar1(6)
        dmet(7) = delvar1(7)

        dmet(8) = dmet1(iw,8)                                              !reflected solar
        dmet(9) = dmet1(iw,9)                                              !direct incoming solar
        dmet(10) = dmet1(iw,10)                                            !diffuse incoming solar
        dmet(11) = dmet1(iw,11)                                            !air pressure
        dmet(12) = dmet1(iw,12)                                            !reflected & emitted IR
        dmet(13) = dmet1(iw,13)                                            !measured soil surface temp
      end if
      if(dmet(1) <= 1d-2) dmet(1) = 0d0
      dmet(2) = dmax1(0d0,dmet(2))

      sdir = dmax1(0d0,dmet(9))
      sdif = dmax1(0d0,dmet(10))


! CALCULATE THE SHORT and LONG WAVE RADIATION TERM
! correct for a surface slope based on Jordan,R., CRREL Rep. 91-16, p.25.
! solar zenith angle[solzen(degrees)]
! solar azimuthal angle; + = clockwise [saz(degrees from north)]
! aspect: 180 = true South; 360 = true North

      if(slope > eps) then
        f2 = pi/180d0
        solzen = solzen*f2

        if(solzen >= 90d0*f2) then
          dmet(1) = 0d0
        else if(solzen < 90d0*f2.and.dmet(1) > 0d0) then
          if(sloper > eps.and.sloper < 90d0*f2) then
            delaz = dabs(aspect - saz)*f2
! short wave
            if(solzen <= 85d0*f2) then
              sdircorr =(dsin(solzen)*dsin(sloper)*dcos(delaz))         &
                        /dcos(solzen)
            else
              sdircorr =(dsin(solzen)*dsin(sloper)*dcos(delaz))         &
                        /dcos(85d0)
            end if

            sdir = dmax1(0d0,sdir*(dcos(sloper) + sdircorr))
            dmet(9) = sdir

!          difcoeff = 1d0 + 0.5d0*sin(solzen) + 2d0*sin(2d0*solzen)
!          sdifcorr1 = (pi - delaz)*(difcoeff + dcos(sloper))
!          sdifcorr2 = delaz*(1d0 + difcoeff*dcos(sloper))
!          sdif = dmax1(0d0,sdif*((sdifcorr1 + sdifcorr2)/(pi*(1d0 + difcoeff))))

            pre_slope = dmet(1)
            dmet(1) = sdir + sdif                                       &
                         + dmet(1)*albedo*(1d0 - dcos(sloper))*5d-1
            if(dmet(1) < 1d-2) dmet(1) = 0d0

            dmet(10) = dmet(1) - dmet(9)

            if(aint(dabs(dmet(8)-mflag)*1d5)*1d-5 > eps.and.           &
                 pre_slope > eps) dmet(8) = dmet(8)*dmet(1)/pre_slope

! long wave (Sellers et al. (1995), JGR 100(D12), p.25607-25629)
!         dmet(2) = dmax1(0d0,dmet(2)*((pi - sloper)/pi)/dcos(sloper))
          end if
        end if
      end if

      if(dabs(dmet(1)) <= eps) then                                      !no up-solar if none coming in
        dmet(8) = 0d0
      else if(dmet(8) >= dmet(1).and.                                   &
                        aint(dabs(dmet(8)-mflag)*1d5)*1d-5 > eps) then   !can't have more going out than coming in
        dmet(8) = albedo*dmet(1)
      end if

      do i=1,13
        dmet(i) = anint(dmet(i)*1d20)*1d-20
      end do

      end subroutine sub_divide
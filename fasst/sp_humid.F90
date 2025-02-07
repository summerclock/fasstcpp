      subroutine sp_humid(ic,apres,temp1,rh2,ph,mixr,dmrdt,vpress,      &
                          wetbulb,rhov,rhoda,vpsat,thvc,dthvdt,dthvdh)

      use fasst_global

! formulation for vpress comes from the WES program teten.c
! no subroutines called

      implicit none

      integer(kind=4),intent(in):: ic
      real(kind=8),intent(in):: apres,temp1,rh2,ph
      real(kind=8),intent(out):: mixr,dmrdt,vpress,wetbulb,vpsat
      real(kind=8),intent(out):: rhov,rhoda,thvc,dthvdt,dthvdh

! local variables
      integer(kind=4):: di0
      real(kind=8):: t1,t2,dedt,a,b,desdt,c1,c2,t3,t4,dewpt,delta,deltad
      real(kind=8):: deltap,ap,dqdt,sphumid,dRHdt,dRHdh,dedh,drvdh,dradh
      real(kind=8):: dmrdh,t5,drvdt,drvsdt,dradt,rhovs,d0,rhow,dense,tc

! zero-out parameters
      di0 = 0
      t1 = 0d0
      t2 = 0d0
      t3 = 0d0
      t4 = 0d0
      t5 = 0d0
      a = 0d0
      b = 0d0
      c1  = 0d0
      c2 = 0d0
      ap = 0d0
      dRHdt = 0d0
      dRHdh = 0d0
      vpsat = 0d0
      vpress = 0d0
      desdt = 0d0
      dedt = 0d0
      dedh = 0d0
      rhovs = 0d0
      rhov = 0d0
      drvsdt = 0d0
      drvdt = 0d0
      drvdh = 0d0
      rhoda = 0d0
      dradt = 0d0
      dradh = 0d0
      sphumid = 0d0
      mixr = 0d0
      dqdt = 0d0
      dmrdt = 0d0
      dmrdh = 0d0
      delta = 0d0
      deltad = 0d0
      deltap = 0d0
      dewpt = 0d0
      wetbulb = 0d0
      thvc = 0d0
      dthvdt = 0d0
      dthvdh = 0d0
      d0 = 0d0
      rhow = 0d0
      tc = 0d0

      dRHdt = -ph*grav/(Rv*temp1*temp1)                                  !1/K (dRH/dtemp1)
      dRHdh = grav/(Rv*temp1)                                            !1/m (dRH/dh)

      if(ic == 0) then                                                   !over water
        a = 17.269d0
        b = 35.86d0
      else                                                               !over ice/snow
        a = 21.8745d0
        b = 7.66d0
      end if

      ap = apres*1d2                                                     !Pa
      t1 = temp1 - b                                                     !K
      tc = temp1 - Tref                                                  !C

! vapor pressures and their derivatives
      if(dabs(a*tc/t1) > 5d1.or.dabs(t1) <= eps) then
        vpress = ap                                                      !boiling point
        if(rh2 > eps) then
          vpsat = vpress/rh2
        else
          vpsat = vpress
        end if
      else
        vpsat = dmin1(ap,610.78d0*dexp(a*tc/t1))                         !Pa (saturation vapor pressure)
        vpress = dmin1(vpsat,rh2*vpsat)                                  !Pa (vapor pressure)

        desdt = vpsat*a*(1d0/t1 - tc/(t1*t1))                            !Pa/K (dvpsat/dttemp1)
        dedt = vpress*(dRHdt + a*(1d0/t1 - tc/(t1*t1)))                  !Pa/K (dvpress/dtemp1)
        dedh = dRHdh*vpress                                              !Pa/m (dvpress/dh)
      end if

! dry air and vapor densities and their derivatives
      if(dabs(temp1) > eps) then
        c1 = 1d0/(Rv*temp1)                                              !kg/J
        c2 = 1d0/(Rd*temp1)                                              !kg/J
        t2 = ap - vpress                                                 !Pa

        rhovs = c1*vpsat                                                 !kg/m^3 (saturated water vapor density)
        rhov = c1*vpress                                                 !kg/m^3 (water vapor density)

        drvsdt = c1*(desdt - vpsat/temp1)                                !kg/m^3*K (drhovs/dtemp1)
        drvdt = c1*(dedt - vpress/temp1)                                 !kg/m^3*K (drhov/dtemp1)
        drvdh = c1*dedh                                                  !kg/m^4 (drhov/dh)

        rhoda = dmax1(0.95d0,dmin1(2.8d0,t2*c2))                         !kg/m^3 (dry air density)
        dradt = -c2*(dedt + t2/temp1)                                    !kg/m^3*K (drhoa/dtemp1)
        dradh = -c2*dedh                                                 !kg/m^4 (drhoa/dh)
      else
        rhoda = 0.95d0
      end if

! specific humidity and mixing ratio and their derivatives
! NOTE: 0.622 = Rd/Rv
      if(dabs(t2) > eps) then
        mixr = 0.622d0*vpress/t2                                         !unitless [kg/kg] (mixing ratio)
        sphumid = mixr/(1d0 + mixr)                                      !unitless [kg/kg] (specific humidity)

        dqdt = (0.622d0/ap)*dedt                                         !1/K (dsphumid/dtemp1)
        dmrdt = dedt*0.622d0*(1d0/t2 + vpress/(t2*t2))                   !1/K (dmr/dtemp1)
        dmrdh = dedh*0.622d0*(1d0/t2 + vpress/(t2*t2))                   !1/m (dmr/dh)
      end if

! dew point
      if(vpress <= 1d-15) vpress = 0d0 !1d-15                            !prevent program crashes
      if(dabs(vpress) > eps.and.rh2 > eps) then
        t3 = dlog(vpress/(rh2*610.78d0))                                 !unitless
      else
        t3 = 0d0
      end if

      if(dabs(t3-a) > eps) then
        dewpt = dmax1(temp1,(t3*b - Tref*a)/(t3 - a))                    !K (dewpt temperature)
      else
        dewpt = temp1
      end if

! wet bulb temperature
      if(dabs(t1) > eps) delta = 4.099d6*vpsat/(t1*t1)                   !Pa
      t4 = dewpt - b                                                     !K
      if(dabs(t4) > eps) then
        deltad = 4.099d6*vpress/(t4*t4)                                  !Pa
      else
        deltad = 0d0
      end if
      deltap = (delta + deltad)*5d-1                                     !Pa

      wetbulb = temp1 - ((vpsat - vpress)/(deltap + 6.6d0*apres))        !K (wetbulb temperature)

! dthetav/dT and dthetav/dh
!      t5 = rhov + mixr*rhoda
      di0 = 1
      rhow = dense(temp1,d0,di0)
      t5 = rhow + mixr*rhoda
      if(dabs(t5) > eps) then
        thvc = mixr*rhoda/t5
!        dthvdt = (dmrdt*rhoda + dradt*mixr -                            &
!                            thvc*(drvdt + dmrdt*rhoda + dradt*mixr))/t5
!        dthvdh = (dmrdh*rhoda + dradh*mixr -                            &
!                            thvc*(drvdh + dmrdh*rhoda + dradh*mixr))/t5
        dthvdt = (dmrdt*rhoda + dradt*mixr -                            &
                                    thvc*(dmrdt*rhoda + dradt*mixr))/t5
        dthvdh = (dmrdh*rhoda + dradh*mixr -                            &
                                    thvc*(dmrdh*rhoda + dradh*mixr))/t5
      end if

      wetbulb = anint(wetbulb*1d20)*1d-20
      vpsat = anint(vpsat*1d20)*1d-20
      sphumid = anint(sphumid*1d20)*1d-20
      dqdt = anint(dqdt*1d20)*1d-20
      vpress = anint(vpress*1d20)*1d-20
      drvdt = anint(drvdt*1d20)*1d-20
      rhovs = anint(rhovs*1d20)*1d-20
      rhov = anint(rhov*1d20)*1d-20
      rhoda = anint(rhoda*1d20)*1d-20
      drvsdt = anint(drvsdt*1d20)*1d-20
      dradt = anint(dradt*1d20)*1d-20
      mixr = anint(mixr*1d20)*1d-20
      dmrdt = anint(dmrdt*1d20)*1d-20
      thvc = anint(thvc*1d20)*1d-20
      dthvdt = anint(dthvdt*1d20)*1d-20
      dthvdh = anint(dthvdh*1d20)*1d-20

      end subroutine sp_humid
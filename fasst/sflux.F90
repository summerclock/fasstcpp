      subroutine sflux(vflag,s,sn,swu,swd,lwu,lwd,sfac,temp1,cc1,taf,   &
                       ph,sh,lh,sphumidaf,sph,kth,temp2,h,r1,r2,iceo,   &
                       wvco,smo,dqdT1,dqdT2,d1,radswd,radswu,radlwd,    &
                       radlwu,phw,shw,lhw,chw,chw1,net,dnet1,dnet2)

      use fasst_global

      implicit none

! no subroutines called
! uses the function: dense,spheats

      integer(kind=4),intent(in):: vflag,s,sn
      real(kind=8),intent(in):: swu,swd,lwu,lwd,sfac,temp1,cc1
      real(kind=8),intent(in):: ph,sh,lh,sphumidaf,sph,kth,taf,d1
      real(kind=8),intent(in):: temp2,h,r1,r2,iceo,wvco,smo,dqdT1,dqdT2
      real(kind=8),intent(out):: radswd,radswu,radlwd,radlwu,phw,shw
      real(kind=8),intent(out):: lhw,chw,chw1,net,dnet1,dnet2
 
! local variables
      integer(kind=4):: vf,d1i
      real(kind=8):: alb,sgl,sgl1,sfc,tt1,ems,fchw,fheat,f1,f2,f3
      real(kind=8):: cheati,cheatv,cheatw,W,gamma,dgdT,t1,dense,spheats
      real(kind=8):: dradlwu1,dphw
      real(kind=8):: dshw,dlhw,dchw,dchw1,dshwfg,dlhwfg,dradlwu2


! if vflag = 1: T1 = ftemp, T2 = toptemp
! if vflag = 0: T1 = toptmep, T2 = ftemp
      vf = 0
      d1i = 0
      alb = 0d0                    !surface albedo
      sgl = 0d0                    !low veg density if vflag = 1; 1-(low veg density) if vflag = 0
      sgl1 = 0d0
      sfc = 0d0                    !amount of short wave radiation penetrating surface
      tt1 = 0d0                    !T1 (K)
      ems = 0d0                    !surface emissivity
      fchw = 0d0                   !conductive heat from veg to ground surface (+/-), W/m^2
      fheat = 0d0                  !heat flux due to water flow through the snow (+/-), W/m^2
      cheati = 0d0                 !heat flux due to ice/water phase change (+/-), W/m^2
      cheatv = 0d0                 !heat flux due to vapor/water phase change (+/-), W/m^2
      cheatw = 0d0                 !heat flux due to water flow (+/-), W/m^2
      W = 0d0
      gamma = 0d0                  !surface tension (N/m)
      dgdT = 0d0
      t1 = 0d0
      f1 = 0d0
      f2 = 0d0
      dradlwu1 = 0d0               !d(radlwu)/dT1
      dradlwu2 = 0d0               !d(radlwu)/dT2
      dphw = 0d0                   !d(phw)/dT1
      dshw = 0d0                   !d(shw)/dT1
      dlhw = 0d0                   !d(lhw)/dT1
      dchw = 0d0                   !d(chw)/dT1
      dchw1 = 0d0                  !d(chw1)/dT1

      radswd = 0d0                 !energy from incoming/downwelling solar radiation (+), W/m^2
      radswu = 0d0                 !energy from upwelling/reflected solar radiation (-), W/m^2
      radlwd = 0d0                 !energy from incoming longwave/IR radiation (+), W/m^2
      radlwu = 0d0                 !energy from upwelling/emitted longwave/IR radiation (-), W/m^2
      phw = 0d0                    !precipitation heat flux (+/-), W/m^2
      shw = 0d0                    !sensible heat flux (+/-), W/m^2
      lhw = 0d0                    !latent heat flux (+/-), W/m^2
      chw = 0d0                    !concuctivie heat flux (+/-), W/m^2
      chw1 = 0d0                   !net heat flux due to water/vapor flow, phase changes (+/-), W/m^2
      net = 0d0                    !sum of heat fluxes (+/-), W/m^2
      dnet1 = 0d0
      dnet2 = 0d0
      dshwfg = 0d0                 !d(shw)/dT2
      dlhwfg = 0d0                 !d(lwh)/dT2

! sign convention: > 0 -> towards the surface
!                  < 0 -> away from the surface
!                  positive z is upwards relative to sea level

      if(vflag == 1) then                                                !low veg
        vf = 1
        alb = albf
        sgl = sigfl
        sgl1 = sigfl
        sfc = 1d0
        tt1 = ftemp
        ems = epf
        f1 = 6d-1
        f2 = 1d-1
      else                                                               !bare soil/snow
        vf = -1
        alb = albedo
        sgl = 1d0 - sigfl
        sgl1 = dmax1(0d0,1d0 - sigfl*5d-1)
        sfc = sfac
        tt1 = temp1
        ems = emis
        f1 = 1d-1
        f2 = 6d-1

! foliage conduction
        if(hfol_tot > eps) then
          fchw = sgl*kveg*(ftemp - temp1)/hfol_tot
          fchw = anint(fchw*1d20)*1d-20
        end if

! water + vapor flow conduction
        fheat = (flowu(sn) - flowl(sn))*dabs(temp1 - Tref)
        fheat = anint(fheat*1d20)*1d-20

! phase change; warming/cooling of pore-space materials
        if(hsaccum <= eps) then
          f3 = 0d0
          f3 = h/(deltat*s)
! ice
          if(ice(sn)+iceo > eps) then
            d1i = 2
            cheati = (lhfus + spheats(tt1,d1i)*dabs(temp1 - Tref))     & 
                             *dense(temp1,0d0,d1i)*(ice(sn) - iceo)*f3  !W/m^2 (+/-)
          end if
! vapor
          if(wvc(sn)+wvco > eps) then
            d1i = 1
            cheatv = (2500775.6d0 + spheats(tt1,d1i)*dabs(temp1 - Tref))&
                               *dense(temp1,0d0,d1i)*((wvc(sn) - wvco)  &
                                                         *f3) + fv1(sn)  !W/m^2 (+/-)
          end if
! water
          if(soil_moist(sn)+smo > eps) then
            d1i = 1
            cheatw = dense(temp1,0d0,d1i)*(spheats(temp1,d1i)           &
                                                   *dabs(temp1 - Tref)  &
                            + grav*phead(sn))*(soil_moist(sn) - smo)*f3  !W/m^2 (+/-)
          end if
        end if
        chw1 = -cheati + cheatv + cheatw + fheat + fchw
        chw1 = anint(chw1*1d20)*1d-20
      end if

! solar/short-wave
!      radswd = sgl*swd*sfc
      radswd = sgl1*swd*sfc
      radswd = anint(radswd*1d20)*1d-20                                  !W/m^2 (+)
      if(aint(dabs(swu-mflag)*1d5)*1d-5 > eps) then
!        radswu = -sgl*swu*sfc
        radswu = -sgl1*swu*sfc
      else
!        radswu = -sgl*swd*alb*sfc
        radswu = -sgl1*swd*alb*sfc
      end if
      radswu = anint(radswu*1d20)*1d-20                                  !W/m^2 (-)

! ir/long-wave
      radlwd = sgl*lwd
!      radlwd = sgl1*lwd
      radlwd = anint(radlwd*1d20)*1d-20                                  !W/m^2 (+)
      if(aint(dabs(lwu-mflag)*1d5)*1d-5 > eps) then
        radlwu = -sgl*lwu
!        radlwu = -sgl1*lwu
        dradlwu1 = 0d0
        dradlwu2 = 0d0
      else
        if(temp1 > eps.and.ftemp > eps)                                 &
          radlwu = -sgl*(1d0 - ems)*lwd                                 &
                                      - sgl*ems*sigma*(tt1*tt1*tt1*tt1) &
                                    + vf*cc1*(temp1**4d0 - ftemp**4d0)
!          radlwu = -sgl1*(1d0 - ems)*lwd                                 &
!                                      - sgl1*ems*sigma*(tt1*tt1*tt1*tt1) &
!                                    + vf*cc1*(temp1**4d0 - ftemp**4d0)
        if(vflag == 0) then                                              !ground
          dradlwu1 = (-sgl*ems*sigma + vf*cc1)*(4d0*tt1*tt1*tt1)
!          dradlwu1 = (-sgl1*ems*sigma + vf*cc1)*(4d0*tt1*tt1*tt1)
          dradlwu2 = -vf*cc1*4d0*ftemp*ftemp*ftemp
        else                                                             !vegetation
          dradlwu1 = (-sgl*ems*sigma - vf*cc1)*(4d0*tt1*tt1*tt1)
!          dradlwu1 = (-sgl1*ems*sigma - vf*cc1)*(4d0*tt1*tt1*tt1)
          dradlwu2 = vf*cc1*4d0*temp2*temp2*temp2
        end if
      end if
      radlwu = anint(radlwu*1d20)*1d-20                                  !W/m^2 (-)
      dradlwu1 = anint(dradlwu1*1d20)*1d-20
      dradlwu2 = anint(dradlwu2*1d20)*1d-20

! precipitation
      phw = ph*(ptemp - tt1)
      phw = anint(phw*1d20)*1d-20                                        !W/m^2 (+/-)
      dphw = anint(ph*1d20)*1d-20

! sensible
      shw = sh*(taf - tt1)
      shw = anint(shw*1d20)*1d-20                                        !W/m^2 (+/-)
      dshw = sh*(f1*sigfl - 1d0)                                         !w.r.t. T1
      dshw = anint(dshw*1d20)*1d-20
      dshwfg = sh*f2*sigfl                                               !w.r.t. T2
      dshwfg = anint(dshwfg*1d20)*1d-20

!latent
      lhw = lh*(sphumidaf - r1*sph)
      lhw = anint(lhw*1d20)*1d-20                                        !W/m^2 (+ = condensation/- = evaporation)
      dlhw = lh*(f1*sigfl - d1)*r1*dqdT1/d1                              !w.r.t. T1
      dlhw = anint(dlhw*1d20)*1d-20
      dlhwfg = lh*f2*sigfl*r2*dqdT2/d1                                   !w.r.t. T2
      dlhwfg = anint(dlhwfg*1d20)*1d-20

! ground conduction
      if(h /= 0d0) then
        if(hsaccum > eps) then
          chw = sgl*kth*(temp2 - tt1)/h
          chw = anint(chw*1d20)*1d-20                                    !W/m^2 (+/-)
          dchw = -anint(sgl*(kth/h)*1d20)*1d-20
        else
          chw = kth*(temp2 - tt1)/h
          chw = anint(chw*1d20)*1d-20                                    !W/m^2 (+/-)
          dchw = -anint((kth/h)*1d20)*1d-20
        end if
      end if

      net = radswd + radswu + radlwd + radlwu + shw + lhw + chw + chw1  &
            + phw
      net = aint(net*1d18)*1d-18

      dnet1 = dradlwu1 + dshw + dlhw - dphw                              !w.r.t. layer of interest; used in soil_tmp.f90
      dnet1 = anint(dnet1*1d18)*1d-18
      dnet2 = dradlwu2 + dshwfg + dlhwfg                                 !w.r.t. non layer of interest; used in soil_tmp.f90
      dnet2 = anint(dnet2*1d18)*1d-18

      end subroutine sflux

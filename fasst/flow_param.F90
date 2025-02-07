      subroutine flow_param(isn,ntemp,qtop,qtopv,simeltsm,zt,delz,qbot, &
                            klhtop,kvhtop,kvttop)

      use fasst_global

! This subroutine calculates:
!     hydraulic conductivity at a node based on van Genuchten (1980) - klh
!      and Niu & Yang (2006) JHM 7(5), pp.937-952 [ice affected]
!     hdr. cond. based on temp gradients based on Hansson et al. (2004) - klt
!     vapor conductivity with pressure gradient based on UNSAT-H, Fayer (2000) - kvh
!     vapor conductivity with temp gradient based on UNSAT-H, Fayer (2000) - kvt
!     d(soil moisture)/d(pressure head) = dsmdh

!     W = J/s; J/m = N

! no subroutines called

! uses the functions: dense,spheats,soilhumid,head

      implicit none

      integer(kind=4),intent(in):: isn,ntemp
      real(kind=8),intent(in):: qtop,qtopv,simeltsm
      real(kind=8),intent(in):: zt(ntot),delz(ntot)
      real(kind=8),intent(out):: qbot,klhtop,kvhtop,kvttop

! internal variables
      integer(kind=4):: i,d1i,d2i,impflag(ntot)
      real(kind=8):: dense,spheats,rhow,w,hc1,c1,c2,h0,fvl,fvu,fll,flu
      real(kind=8):: kl,D,eta,Hr,a,b,t1,vpsat,desdt,rhovs,drvsdt,cv,kv
      real(kind=8):: kdenom,vinlw,vinuw,vinlv,vinuv,dh,dT,lhe,sph1,sph2
      real(kind=8):: knum,cw,gamma,dgdt,soilhumid,G,ptemp1,rd1
      real(kind=8):: khll,khlv,khul,khuv,klbot,kvbot,phtemp,slpr,tsign
      real(kind=8):: fvl1,fvu1,klt(ntot),kvt(ntot),klh(ntot),kvh(ntot)
      real(kind=8):: head,Frz


! initialize variables
      d1i = 0
      d1i = 1
      d2i = 0
      w = 0d0
      hc1 = 0d0
      c1 = 0d0
      c2 = 0d0
      h0 = 0d0
      fvl = 0d0
      fvu = 0d0
      fll = 0d0
      flu = 0d0
      cw = 0d0
      cv = 0d0
      kv = 0d0
      kl = 0d0
      D = 0d0
      eta = 0d0
      Hr = 0d0
      a = 0d0
      b = 0d0
      t1 = 0d0
      vpsat = 0d0
      desdt = 0d0
      rhovs = 0d0
      drvsdt = 0d0
      kdenom = 0d0
      vinlw = 0d0
      vinuw = 0d0
      vinlv = 0d0
      vinuv = 0d0
      dh = 0d0
      dT = 0d0
      lhe = 0d0
      qbot = 0d0
      fvl1 = 0d0
      fvu1 = 0d0
      khll = 0d0
      khlv = 0d0
      khul = 0d0
      khuv = 0d0
      phtemp = 0d0
      ptemp1 = 0d0
      slpr = 0d0
      sph1 = 0d0
      sph2 = 0d0
      rd1 = 0d0
      tsign = 0d0
      Frz = 0d0

      do i=1,ntot
        impflag(i) = 0
        dsmdh(i) = 0d0
        klh(i) = 0d0
        klt(i) = 0d0
        kvh(i) = 0d0
        kvt(i) = 0d0
        vin(i) = 0d0                                                     !m/s
        fv1(i) = 0d0                                                     !W/m^2
        flowu(i) = 0d0                                                   !m/s
        flowl(i) = 0d0                                                   !m/s
        khu(i) = 0d0                                                     !m/s
        khl(i) = 0d0                                                     !m/s
      end do

      slpr = 1d0 !dmax1(0d0,dmin1(1d0,dcos(slope*pi/1.8d2)))

! solve for dsmdh, kvh, kvt, klh, klt
      d2i = 2
      do i=isn,nnodes
        if(dabs(stt(i)) <= eps) stt(i) = 1d0                             !safety catch to prevent /0

        if(((ntype(i) == 26.or.ntype(i) == 27).or.ntype(i) == 0).or.  &
          (nsoilp(i,11) <= eps.or.nsoilp(i,12) <= eps)) impflag(i) = 1

        if(impflag(i) == 0) then
          rhow = dense(stt(i),0d0,d1i)                                   !kg/m^3 (water density)

          w = (soil_moist(i) - nsoilp(i,8))/(nsoilp(i,9) - nsoilp(i,8))  !unitless

          if(w > eps) then
            if(dabs(phead(i)) <= eps) then
              ptemp1 = dabs(head(i,0.99d0*nsoilp(i,9)))*1d2 !1d0 !0
            else
              ptemp1 = dabs(phead(i))*1d2                                !cm
            end if

            c1 = (1d0 + (nsoilp(i,10)*ptemp1)**nsoilp(i,11))            &
                                                **(-nsoilp(i,12) - 1d0)  !unitless
            c2 =  nsoilp(i,12)*nsoilp(i,11)*(nsoilp(i,10)**nsoilp(i,11)) !unitless

            dsmdh(i) = (nsoilp(i,9) - nsoilp(i,8))*c2*c1*               &
                                     (ptemp1**(nsoilp(i,11) - 1d0))*1d2  !1/m
          end if

          if(w <= 0.99d0.and.w > 0d0) then
            hc1 = 1d0 - (1d0 - (w**(1d0/nsoilp(i,12))))**nsoilp(i,12)    !unitless
            klh(i) = nsoilp(i,7)*(dsqrt(w))*hc1*hc1*1d-2                 !m/s (hydraulic conductivity due to waterflow)
            klh(i) = dmin1(klh(i),nsoilp(i,7)*1d-2)
            if(ice(i) > eps) then
              Frz = dexp(-3d-1*(1d0 - ice(i)/nsoilp(i,9))) - dexp(-3d-1)
              klh(i) = klh(i)*(1d0 - Frz)
!              klh(i) = klh(i)*(1d0 + 8d0*ice(i))**2d0                    !m/s
            else if(ice(i) <= eps.and.stt(i) > Tref) then
              klh(i) = klh(i)*(5.3888d-1 + 2.096d-2*(stt(i) - Tref))
            end if
          else if(w > 0.99d0) then
            klh(i) = nsoilp(i,7)*1d-2
            if(ice(i) > eps) then
              Frz = dexp(-3d-1*(1d0 - ice(i)/nsoilp(i,9))) - dexp(-3d-1)
              klh(i) = klh(i)*(1d0 - Frz)
!              klh(i) = klh(i)*(1d0 + 8d0*ice(i))**2d0                    !m/s
            else if(ice(i) <= eps.and.stt(i) > Tref) then
              klh(i) = klh(i)*(5.3888d-1 + 2.096d-2*(stt(i) - Tref))
            end if
          end if

          klh(i) = dmax1(0d0,klh(i))

!          if(ntype(i) <= 4) then
!            G = 5d0  !2d0
!          else if(ntype(i) >= 5.and.ntype(i) <= 8) then
!            G = 7d0
!          else if(ntype(i) == 15) then
!            G = 10d0
!          else
!            G = 9d0
!          end if
          if(ntype(i) /= 15) then
            G = 2.5d1*(nsoilp(i,20) + 1d0)
            if(G < 5d0) G = 5d0
            if(G > 1d1) G = 1d1
          else
            G = 1d1
          end if

          t1 = stt(i) - Tref
          if(t1 <= eps) then
            klt(i) = 0d0
          else
            gamma = 1d-3*(7.56d1 - 1.425d-1*t1 - 2.38d-4*t1*t1)          !surface tension (N/m)
            dgdt = 1d-3*(-1.425d-1 - 2d0*2.38d-4*t1)                     !d(gamma)/dT (N/m*K)
            if(dabs(gamma) > eps)                                       &
              klt(i) = dmax1(0d0,klh(i)*(G*phead(i)*dgdt/gamma))         !m^2/s*K (hydraulic conductivity due to tempflow)
          end if

          if(wvc(i) >= 1d-15) then
            Hr = soilhumid(i,phead(i),soil_moist(i),stt(i))

            if(dabs(ice(i)) <= eps) then                                 !over water
              a = 17.269d0
              b = 35.86d0
            else                                                         !over ice/snow
              a = 21.8745d0
              b = 7.66d0
            end if

            t1 = stt(i) - b                                              !K
            if(dabs(a*(stt(i) - Tref)/t1) > 5d1) then
              vpsat = 0d0
            else
              vpsat = 610.78d0*dexp(a*(stt(i) - Tref)/t1)                !Pa (saturation vapor pressure)
            end if

            rhovs = vpsat/(Rv*stt(i))                                    !kg/m^3 (saturated water vapor density)
            desdt = vpsat*a*(1d0 - (stt(i) - Tref)/t1)/t1                !Pa/K (saturated)
            drvsdt = (desdt - vpsat/stt(i))/(Rv*stt(i))                  !kg/m^3*K (drhovs/dtemp1)

            t1 = dmax1(0d0,nsoilp(i,2) - (soil_moist(i) + ice(i)))
            D = 0d0
            if(t1 > 0d0) then
              D = (t1**(5d0/3d0))*(2.12d-5*(stt(i)/Tref)**2d0)           !m^2/s (diff. coeff. for water vapor in air)
              D = dmax1(0d0,D)
            end if

            kvh(i) = D*rhovs*grav*Hr/(rhow*Rv*stt(i))                    !m/s
            kvh(i) = dmax1(0d0,kvh(i))

            if(nsoilp(i,20) > eps) then
              t1 = ((1d0 + 2.6d0/dsqrt(nsoilp(i,20)*1d2))               &
                                                   *soil_moist(i))**4d0
              if(t1 > 5d1) t1 = 5d1
              eta = 9.5d0 + 6d0*soil_moist(i) - 8.5d0*dexp(-t1)          !unitless
            else
              eta = 0d0
            end if

            kvt(i) = dmax1(0d0,D*eta*Hr*drvsdt/rhow)                     !m^2/s*K
          end if   ! wvc(i) >= 1d-15 
        end if  !else (i <= nnodes)

        dsmdh(i) = anint(dsmdh(i)*1d20)*1d-20                            !1/m
        klh(i) = anint(klh(i)*slpr*1d20)*1d-20                           !m/s
        kvh(i) = anint(kvh(i)*slpr*1d20)*1d-20                           !m/s
        klt(i) = anint(klt(i)*slpr*1d20)*1d-20                           !m^2/s*K
        kvt(i) = anint(kvt(i)*slpr*1d20)*1d-20                           !m^2/s*K
      end do

      if(dabs(zt(1)-gwl) > eps) then
        dh = (phead(1) - 0d0)/(zt(1) - gwl)                              !head at gwl = 0
        dT = 0d0 !dexp(-2.08d0*(zt(1) - gwl))/(zt(1) - gwl)
      else
        dh = 0d0
        dT = 0d0
      end if
      klhtop = klh(nnodes)
      kvhtop = kvh(nnodes)
      kvttop = kvt(nnodes)

      klbot = -(klh(1)*(dh + 1d0) + klt(1)*dT)                           !m/s
      kvbot = -(kvh(1)*dh + kvt(1)*dT)                                   !m/s
      qbot = (klbot + kvbot)

      klbot = anint(klbot*1d20)*1d-20
      kvbot = anint(kvbot*1d20)*1d-20

! calculate the averages of klh and kvh based on Cherry and Freeze assuming perpendicular flow
      if(isn == 1) then
      d2i = 0
      do i=1,ntemp  !ntot
        if((node_type(i) /= 'AI'.and.node_type(i) /= '  ')              &
                                             .and.impflag(i) == 0) then  !not air and not impervious
          lhe = 2500775.6d0 - 2369.729d0*(stt(i) -Tref)                  !latent heat of evaporation, J/kg = (m/s)^2

! top node, snow, veg
          if(i >= nnodes) then
            if(node_type(i) == 'HM') then                               !snow on ground
              fvl1 = 0d0
              fvu1 = 0d0
              fvl = 0d0
              fvu = 0d0
              vinlw = 0d0
              vinuw = 0d0
              vinlv = 0d0
              vinuv = 0d0

              fll = -sphm*simeltsm/deltat                                !W/m^2*K
              flu = -sphm*atopf/deltat/float(step)                       !W/m^2*K
            else if(node_type(i) == 'VG') then
              fvl1 = 0d0
              fvu1 = 0d0
              fvl = 0d0
              fvu = 0d0
              vinlw = 0d0
              vinuw = 0d0
              vinlv = 0d0
              vinuv = 0d0
              fll = 0d0
              flu = 0d0
            else
              rd1 = dense(stt(i),0d0,d1i)
              sph1 = spheats(stt(i),d1i)*rd1
              sph2 = spheats(stt(i),d2i)*dense(stt(i),0d0,d2i)

              knum = delz(i-1) + delz(i)                                 !m
              kdenom = klh(i-1)*delz(i) + klh(i)*delz(i-1)               !m^2/s
              if(dabs(kdenom) /= 0d0) then
                khll = klh(i)*klh(i-1)*knum/kdenom                       !m/s
              else
                khll = 0d0
              end if

              kdenom = kvh(i-1)*delz(i) + kvh(i)*delz(i-1)               !m^s/s
              if(dabs(kdenom) /= 0d0) then
                khlv = kvh(i)*kvh(i-1)*knum/kdenom                       !m/s
              else
                khlv = 0d0
              end if
              khl(i) = (khll + khlv)/(zt(i) - zt(i-1))                   !1/s

              kv = 5d-1*(kvt(i-1) + kvt(i))                              !m^2/s*K
              kl = 5d-1*(klt(i-1) + klt(i))                              !m^2/s*K
              cw = 5d-1*(sph1                                           & 
                       + spheats(stt(i-1),d1i)*dense(stt(i-1),0d0,d1i))  !J/m^3*K
              cv = 5d-1*(sph2                                           &
                       + spheats(stt(i-1),d2i)*dense(stt(i-1),0d0,d2i))  !J/m^3*K
     
              dh = (phead(i) - phead(i-1))/(zt(i) - zt(i-1))             !m/m
              dT = (stt(i) - stt(i-1))/(zt(i) - zt(i-1))                 !K/m

              vinlw = -(khll*(dh + 1d0) + kl*dT)                         !m/s
              fll = cw*vinlw                                             !W/m^2*K
              vinlv = -(khlv*dh + kv*dT)                                 !m/s
              fvl = cv*vinlv                                             !W/m^2*K
              fvl1 = lhe*rd1*vinlv                                       !W/m^2

              khul = 0d0                                                 !m/s
              khuv = 0d0                                                 !m/s
              khu(i) = 0d0                                               !1/s

              vinuw = -qtop - qtopv                                             !m/s
              vinuv = 0d0 !-qtopv                                             !m/s

              if(hm > eps) vinuv = 0d0
              flu = sph1*vinuw                                           !W/m^2*K
              fvu = sph2*vinuv                                           !W/m^2*K
              fvu1 = lhe*rd1*vinuv                                       !W/m^2
            end if

! bottom node
          else if(i == 1) then
            rd1 = dense(stt(i),0d0,d1i)
            sph1 = spheats(stt(i),d1i)*rd1
            sph2 = spheats(stt(i),d2i)*dense(stt(i),0d0,d2i)

            khll = 0d0                                                   !m/s
            khlv = 0d0                                                   !m/s
            khl(i) = 0d0                                                 !1/s

            vinlw = klbot                                                !m/s
            fll = sph1*vinlw                                             !W/m^2*K
            vinlv = kvbot                                                !m/s
            fvl =  sph2*vinlv                                            !W/m^2*K
            fvl1 = lhe*rd1*vinlv                                         !W/m^2

            knum = delz(i) + delz(i+1)                                   !m
            kdenom = klh(i)*delz(i+1) + klh(i+1)*delz(i)                 !m^2/s
            if(dabs(kdenom) /= 0d0) then
              khul = klh(i)*klh(i+1)*knum/kdenom                         !m/s
            else
              khul = 0d0
            end if
            kdenom = kvh(i)*delz(i+1) + kvh(i+1)*delz(i)                 !m^2/s
            if(dabs(kdenom) /= 0d0) then
              khuv = kvh(i)*kvh(i+1)*knum/kdenom                         !m/s
            else
              khuv = 0d0
            end if
            khu(i) = (khul + khuv)/(zt(i+1) - zt(i))                     !1/s

            kv = 5d-1*(kvt(i) + kvt(i+1))                                !m^2/s*K
            kl = 5d-1*(klt(i) + klt(i+1))                                !m^2/s*K
            cw = 5d-1*(sph1                                             &
                         + spheats(stt(i+1),1)*dense(stt(i+1),0d0,d1i))  !J/m^3*K 
            cv = 5d-1*(sph2                                             &
                         + spheats(stt(i+1),0)*dense(stt(i+1),0d0,d2i))  !J/m^3*K

            dh = (phead(i+1) - phead(i))/(zt(i+1) - zt(i))               !m/m
            dT = (stt(i+1) - stt(i))/(zt(i+1) - zt(i))                   !K/m

            vinuw = -(khul*(dh + 1d0) + kl*dT)                           !m/s
            flu = cw*vinuw                                               !W/m^2*K
            vinuv = -(khuv*dh + kv*dT)                                   !m/s
            fvu = cv*vinuv                                               !W/m^2*K
            fvu1 = lhe*rd1*vinuv                                         !W/m^2

! interior nodes
          else if(i > 1.and.i < nnodes) then
            rd1 = dense(stt(i),0d0,d1i)
            sph1 = spheats(stt(i),d1i)*rd1
            sph2 = spheats(stt(i),d2i)*dense(stt(i),0d0,d2i)

            knum = delz(i-1) + delz(i)                                   !m
            kdenom = klh(i-1)*delz(i) + klh(i)*delz(i-1)                 !m^2/s
            if(dabs(kdenom) /= 0d0) then
              khll = klh(i)*klh(i-1)*knum/kdenom                         !m/s
            else
              khll = 0d0
            end if
            kdenom = kvh(i-1)*delz(i) + kvh(i)*delz(i-1)                 !m^s/s
            if(dabs(kdenom) /= 0d0) then
              khlv = kvh(i)*kvh(i-1)*knum/kdenom                         !m/s
            else
              khlv = 0d0
            end if
            khl(i) = (khll + khlv)/(zt(i) - zt(i-1))                     !1/s

            kv = 5d-1*(kvt(i-1) + kvt(i))                                !m^2/s*K
            kl = 5d-1*(klt(i-1) + klt(i))                                !m^2/s*K
            cw = 5d-1*(sph1                                             & 
                       + spheats(stt(i-1),d1i)*dense(stt(i-1),0d0,d1i))  !J/m^3*K
            cv = 5d-1*(sph2                                             &
                       + spheats(stt(i-1),d2i)*dense(stt(i-1),0d0,d2i))  !J/m^3*K

            dh = (phead(i) - phead(i-1))/(zt(i) - zt(i-1))               !m/m
            dT = (stt(i) - stt(i-1))/(zt(i) - zt(i-1))                   !K/m

            vinlw = -(khll*(dh + 1d0) + kl*dT)                           !m/s
            fll = cw*vinlw                                               !W/m^2*K
            vinlv = -(khlv*dh + kv*dT)                                   !m/s
            fvl = cv*vinlv                                               !W/m^2*K
            fvl1 = lhe*rd1*vinlv                                         !W/m^2

            knum = delz(i) + delz(i+1)                                   !m
            kdenom = klh(i)*delz(i+1) + klh(i+1)*delz(i)                 !m^2/s
            if(dabs(kdenom) /= 0d0) then
              khul = klh(i)*klh(i+1)*knum/kdenom                         !m/s
            else
              khul = 0d0
            end if
            kdenom = kvh(i)*delz(i+1) + kvh(i+1)*delz(i)                 !m^2/s
            if(dabs(kdenom) /= 0d0) then
              khuv = kvh(i)*kvh(i+1)*knum/kdenom                         !m/s
            else
              khuv = 0d0
            end if
            khu(i) = (khul + khuv)/(zt(i+1) - zt(i))                     !1/s

            kv = 5d-1*(kvt(i) + kvt(i+1))                                !m^2/s*K
            kl = 5d-1*(klt(i) + klt(i+1))                                !m^2/s*K
            cw = 5d-1*(sph1                                             &
                         + spheats(stt(i+1),1)*dense(stt(i+1),0d0,d1i))  !J/m^3*K 
            cv = 5d-1*(sph2                                             &
                         + spheats(stt(i+1),0)*dense(stt(i+1),0d0,d2i))  !J/m^3*K
            dh = (phead(i+1) - phead(i))/(zt(i+1) - zt(i))               !m/m
            dT = (stt(i+1) - stt(i))/(zt(i+1) - zt(i))                   !K/m

            vinuw = -(khul*(dh + 1d0) + kl*dT)                           !m/s
            flu = cw*vinuw                                               !W/m^2*K
            vinuv = -(khuv*dh + kv*dT)                                   !m/s
            fvu = cv*vinuv                                               !W/m^2*K
            fvu1 = lhe*rd1*vinuv                                         !W/m^2
          end if

          vin(i) = (vinuw + vinuv) - (vinlw + vinlv)                     !m/s

          fv1(i) = (fvu1 - fvl1)                                         !W/m^2
          flowu(i) = (flu + fvu)
          flowl(i) = (fll + fvl)                                         !W/m^2*K

          khu(i) = anint(khu(i)*1d20)*1d-20
          khl(i) = anint(khl(i)*1d20)*1d-20
          vin(i) = anint(vin(i)*1d20)*1d-20
          fv1(i) = anint(fv1(i)*1d20)*1d-20
          flowu(i) = anint(flowu(i)*1d20)*1d-20
          flowl(i) = anint(flowl(i)*1d20)*1d-20
        end if
      end do
      end if

      end subroutine flow_param

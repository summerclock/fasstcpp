      subroutine th_param(iselect,n,pdens,pdensnew,rhov,rhoda,rhotot)

      use fasst_global

! calcultes the thermal properties of both the soil and snow/ice layers

! no subroutines called

! uses the function: dense,spheats,thconds

      implicit none

      integer(kind=4),intent(in):: iselect,n
      real(kind=8),intent(in):: pdens,pdensnew
      real(kind=8),intent(in):: rhov(n),rhoda(n)
      real(kind=8),intent(out):: rhotot

! local variables
      integer(kind=4):: i,d0i,d1i,d2i,d3i
      real(kind=8):: sr,ke,ksat,ks,kw,ki,gs,gw,gi,ga,gv,sphmi,sphms
      real(kind=8):: thconds,dense,spheats,kmi,kms,kmsn,tice,hms,tdens
      real(kind=8):: sphsnow,ttemp,sphmsn,densgc,ka,ns,kv,na,cs,ksq,kso
      real(kind=8):: ksn,pq,np

! zero-out variables
      d0i = 0
      d1i = 0
      d2i = 0
      d3i = 0
      sr = 0d0
      ke = 0d0
      ksat = 0d0
      ks = 0d0
      kw = 0d0
      ki = 0d0
      ka = 0d0
      kv = 0d0
      gs = 0d0
      gw = 0d0
      gi = 0d0
      ga = 0d0
      gv = 0d0
      ns = 0d0
      na = 0d0
      sphmi = 0d0
      sphms = 0d0
      sphmsn = 0d0
      kmi = 0d0
      kms = 0d0
      kmsn = 0d0
      tice = 0d0
      hms = 0d0
      tdens = 0d0
      densgc = 0d0
      ttemp = 0d0
      sphsnow = 0d0
      rhotot = 0d0
      cs = 0d0
      ksq = 0d0
      kso = 0d0
      ksn = 0d0
      pq = 0d0
      np = 0d0

! ke = Kersten number
! sr = saturation ratio (vol. water/vol. voids)
! nsoilp(i,2) = porosity
! nsoilp(i,5) = quartz content of the soil by fraction
! nsoilp(i,13) = specific heat of dry soil (J/kg*K)
! nsoilp(i,14) = organic fraction of the soil
! grspheat(i) = heat capacity (specific heat)/unit volume of soil (J/m^3*K)

      d1i = 1
      d2i = 2
      d3i = 3

      select case (iselect)
        case(1)  !soils, etc....

! equations for soil thcond are based on Farouki, p.112-116, CRREL Monograph 81-1
          do i=1,nnodes
            if(ice(i)+soil_moist(i) > nsoilp(i,2)) then
              np = ice(i) + soil_moist(i)
            else
              np = nsoilp(i,2)
            end if
            ns = (1d0 - np)

            if(ntype(i) /= 27.and.ntype(i) /= 26) then                   !not water or air
              na = dmax1(0d0,np - (soil_moist(i) + ice(i) + wvc(i)))
            else if(ntype(i) == 27) then
              na = dmax1(0d0,1d0 - wvc(i))
            else
              na = 0d0
            end if

            if(ntype(i) <= 18) then
              kw = thconds(stt(i),d1i)                                   !W/m*K, thermal conductivity of water
              ka = thconds(stt(i),d3i)                                   !W/m*K, thermal conductivity of air
              kv = thconds(stt(i),d0i)                                   !W/m*K, thermal conductivity of water vapor
              sr = dmax1(eps,dmin1(1d0,soil_moist(i)/nsoilp(i,9)))       !unitless

! thermal conductivity of soil solids (W/m*K)
! NOTE: thermal cond of quartz [nsoilp(i,5)]: 7.7 - 8.7 (Johansen); 2-3 (engineeringtoolbox.com); 1.3 (Wikipedia);
!                                             1.4 (www.quartz.com)
!       thermal cond of non-quartz: 2.0 - 3.0 (Johansen); 2-7 (engineeringtoolbox.com)
!       thermal cond of organics [nsoilp(i,14)], ics: ~10%(non-quartz); 0.25 (many sources)
              ksq = 1d0                                                  !quartz
!              if(nsoilp(i,5) > 0d0) ksq = 8.4d0**nsoilp(i,5)  !3d0
              pq = 3.39d-1 + 4.17d-1*(nsoilp(i,18) + nsoilp(i,25))       !Tarnawski et al. (2009)
              if(pq > 0d0) ksq = 8.4d0*pq 
              kso = 1d0                                                  !organics
              if(nsoilp(i,14) > 0d0) kso = 0.25d0**nsoilp(i,14)
              ksn = 1d0                                                  !all other solids
!              if(dabs(1d0 - (nsoilp(i,5) + nsoilp(i,14))) > eps)        &
!                ksn = 2.9d0**(1d0 - (nsoilp(i,5) + nsoilp(i,14)))  !5d0
              if(dabs(1d0 - (pq + nsoilp(i,14))) > eps)                 &
                ksn = 2.9d0**(1d0 - (pq + nsoilp(i,14)))  !5d0
              ks = ksq*kso*ksn                                           !total solids

!              ks = (7.7d0**nsoilp(i,5))*(0.25d0**nsoilp(i,14))          &
!                           *(2d0**(1d0 - (nsoilp(i,5) + nsoilp(i,14))))
!!              if(icourse(i) == 1.and.nsoilp(i,5) < 2d-1)                   &
!              if(nsoilp(i,5) < 2d-1)                                    &
!                ks = (7.7d0**nsoilp(i,5))*(0.25d0**nsoilp(i,14))        &
!                          *(3d0**(1d0 - (nsoilp(i,5) + nsoilp(i,14))))

! saturated thermal conductivity (ksat); wetness factor (ke)
              ke = sr
              if(ntype(i) == 15) then
                if(sr > 5d-2) ke = 0.7d0*dlog10(sr) + 1d0
              else
!                if(sr > 0d0) then
!                  if(nsoilp(i,5) <= 4d-1) then
!                    ke = dexp(0.27d0*(1d0 - sr**(0.27d0 - 1.33d0)))
!                  else
!                    ke = dexp(0.96d0*(1d0 - sr**(0.96d0 - 1.33d0)))
!                  end if
!                end if

                if(sr > 0d0) then
                  if(ntype(i) <= 8.or.ntype(i) == 18) then                 !Lu et al (2007) with Tarnawski (2009) coeff.
                    ke = dexp(7.28d-1*(1d0 - sr**(7.28d-1 - 1.165d0)))
                  else
                    ke = dexp(3.7d-1*(1d0 - sr**(3.7d-1 - 1.29d0)))
                  end if
                end if
              end if
              ke = dmin1(dmax1(0d0,ke),1d0)

              if(ice(i) > eps) then                                      !frozen
!                ke = sr
!                if(ke <= eps) ke = 1d0

                ki = thconds(stt(i),d2i)
                if(soil_moist(i) > eps) then
                  ksat = (ks**ns)*(ki**ice(i))*(kw**(np - ice(i)))       !W/m*K
                else if(soil_moist(i) <= eps) then
                  ksat = (ks**ns)*(ki**ice(i))
                end if
              else                                                       !unfrozen
                ksat = (kw**np)*(ks**ns)                                 !W/m*K
              end if

! effective thermal conductivity (W/m*K)
              grthcond(i) = nsoilp(i,6) + (ksat - nsoilp(i,6))*ke        !W/m*K
            end if

! concrete
            if(ntype(i) == 20) then
              if(dabs(nsoilp(i,6)-spflag) <= eps) then
                grthcond(i) = 0.875d0                                    ! range is 0.05 - 1.7 W/m*K
              else
                grthcond(i) = nsoilp(i,6)
              end if
            end if

! asphalt
            if(ntype(i) == 21) then
              if(dabs(nsoilp(i,6)-spflag) <= eps) then
                grthcond(i) = 0.7d0                                      !range is 0.15 - 1.4 W/m*K
              else
                grthcond(i) = nsoilp(i,6)
              end if
            end if

! bed rock - assumed to be granite
            if(ntype(i) == 25) then
              if(dabs(nsoilp(i,6)-spflag) <= eps) then
                grthcond(i) = 2.0d0                                      !range is 2.0 - 3.5 W/m*K
              else
                grthcond(i) = nsoilp(i,6)
              end if
            end if

! water
            if(ntype(i) == 26) then
              if(dabs(nsoilp(i,6)-spflag) <= eps) then
                grthcond(i) = thconds(stt(i),d1i)                        !W/m*K
              else
                grthcond(i) = nsoilp(i,6)
              end if
            end if
! air
            if(ntype(i) == 27) then
              if(dabs(nsoilp(i,6)-spflag) <= eps) then
                grthcond(i) = thconds(stt(i),d3i)                        !W/m*K
              else
                grthcond(i) = nsoilp(i,6)
              end if
            end if

! glaciers +/- permanent snow
            if(ntype(i) == 30) then
              if(dabs(nsoilp(i,6)-spflag) <= eps) then
                grthcond(i) = thconds(stt(i),d2i)                        !W/m*K
              else
                grthcond(i) = nsoilp(i,6)
              end if
            end if

! specific heat (J/m^3*K)
            if(node_type(i) /= 'WA'.and.dabs(nsoilp(i,13)-spflag)       &
                                                            > eps) then
              gs = ns*(nsoilp(i,1)*1d3)*nsoilp(i,13)                     !solids
            else
              gs = 0d0
            end if

            if(ntype(i) /= 27.and.ntype(i) /= 26) then
              gw = soil_moist(i)*dense(stt(i),0d0,d1i)                  &
                                                    *spheats(stt(i),d1i) !water
              gi = ice(i)*dense(stt(i),0d0,d2i)*spheats(stt(i),d2i)      !ice
              ga = na*rhoda(i)*spheats(stt(i),d3i)                       !air
              gv = wvc(i)*rhov(i)*spheats(stt(i),d0i)                    !water vapor
            else if(ntype(i) == 27) then
              gw = 0d0
              gi = 0d0
              ga = na*rhoda(i)*spheats(stt(i),d3i)
              gv = wvc(i)*rhov(i)*spheats(stt(i),d0i)
            else if(ntype(i) == 26) then
              gw = soil_moist(i)*dense(stt(i),0d0,d1i)                  &
                                                    *spheats(stt(i),d1i)
              gi = 0d0
              ga = 0d0
              gv = 0d0
            end if

            grspheat(i) = gs + gw + gi + ga + gv                         !J/m^3*K

            grthcond(i) = anint(grthcond(i)*1d20)*1d-20
            grspheat(i) = anint(grspheat(i)*1d20)*1d-20

! thermal diffusivity (m^2/s) = grthcond(i)/grspheat(i)
! thermal inertia (J/m^2Ks^0.5) = dsqrt(grspheat(i)*grthcond(i))
          end do

! ******************************************************************************

        case(2)   !snow, ice

! calculate snow/ice thermal conductivity
! determine snow +/- ice thermal conductivy
!     kmi = 2.29d0/((-13.3d0 + 7.8d0*ttemp)*idens)    !m^2/s, th. diff. ice
!     kms = 2.3d-2 + (7.75d-5*sdensw + 1.105d-6*sdensw*sdensw)*(2.29d0 - 2.3d-2)  !W/m*K, th. cond. snow (SNTHERM)
!     kms = kms/((-13.3d0 + 7.8d0*ttemp)*sdensw)      !m^2/s, th. diff. snow

! Sturm et al. (1997), J. Glaciology 43(143),pp.26-41

          km = 0d0
          kmi = 0d0
          kms = 0d0
          kmsn = 0d0
          sphm = 0d0
          sphmi = 0d0
          sphms = 0d0
          sphmsn = 0d0
          ttemp = toptemp
          sphsnow = 2.09d3                                               !J/kg*K

          if(ttemp > Tref) ttemp = Tref

          if(hsaccum > eps) then
            if(iw /= 1) then
              tdens = sdens(iw-1)                                        !kg/m^3
              if(dabs(sdens(iw-1)) <= eps) tdens = pdensnew
            else
              tdens = sdensw
            end if
            densgc = tdens*1d-3                                          !g/cm^3

            if(densgc < 1.56d-1) then                                    !W/m*K, th. cond. (Sturm et al)
              kms = 2.3d-2 + 2.34d-1*densgc
              kms = dmin1(dmax1(2.3d-2,kms),1d0)
            else
              kms = 1.38d-1 - 1.01d0*densgc + 3.233d0*densgc*densgc
              kms = dmin1(dmax1(1.38d-1,kms),1d0)
            end if
            sphms = sphsnow*tdens                                        !J/m^3*K, specific heat of snow
          end if

          if(newsd > eps) then
            densgc = pdens*1d-3
            if(densgc < 1.56d-1) then                                    !W/m*K, th. cond. (Sturm et al)
              kmsn = 2.3d-2 + 2.34d-1*densgc
              kmsn = dmin1(dmax1(2.3d-2,kmsn),1d0)
            else
              kmsn = 1.38d-1 - 1.01d0*densgc + 3.233d0*densgc*densgc
              kmsn = dmin1(dmax1(1.38d-1,kmsn),1d0)
            end if
            sphmsn = sphsnow*pdens                                       !J/m^3*K, specific heat of snow
          end if

          if(hi > eps) then
            if(dabs(hsaccum+newsd) <= eps) then
              km = thconds(ttemp,d2i)                                    !W/m*K, thermal conductivity ice
              sphm = spheats(ttemp,d2i)*dense(ttemp,0d0,d2i)             !J/m^3*K, specific heat of ice
            else
              if(hsaccum > eps.and.dabs(newsd) <= eps) then
                tice = (kms/hsaccum + 2.29d0/hi)/                       &
                       ((2.29d0/hi)*stt(nnodes) + (kms/hsaccum)*ttemp)
                kmi = thconds(tice,d2i)                                  !W/m*K, thermal conductivity ice
                sphmi = spheats(tice,d2i)*dense(tice,0d0,d2i)            !J/m^3*K, specific heat of ice

                rhotot = (tdens*hsaccum + dense(tice,0d0,d2i)*hi)       &
                                                         /(hsaccum + hi)
                km = (hsaccum + hi)/(hi/kmi + hsaccum/kms)
                sphm = (sphmi*dense(tice,0d0,d2i) + sphms*tdens)/rhotot  !J/m^3*K, snow/ice layer
              else if(dabs(hsaccum) <= eps.and.newsd > eps) then
                tice = (kmsn/newsd + 2.29d0/hi)/                        &
                        ((2.29d0/hi)*stt(nnodes) + (kmsn/newsd)*ttemp)
                kmi = thconds(tice,d2i)                                  !W/m*K, thermal conductivity ice
                sphmi = spheats(tice,d2i)*dense(tice,0d0,d2i)            !J/m^3*K, specific heat of ice

                rhotot = (pdens*newsd + dense(tice,0d0,d2i)*hi)         &
                                                           /(newsd + hi)
                km = (hi + newsd)/(hi/kmi + newsd/kmsn)
                sphm = (sphmi*dense(tice,0d0,d2i) + sphmsn*pdens)/rhotot !J/m^3*K, snow/ice layer
              else   !both old and new snow
                hms = hsaccum + newsd
                kms = hms/(hsaccum/kms + newsd/kmsn)

                tice = (kms/hms + 2.29d0/hi)/                           &
                        ((2.29d0/hi)*stt(nnodes) + (kms/hms)*ttemp)
                kmi = thconds(tice,d2i)                                  !W/m*K, thermal conductivity ice
                sphmi = spheats(tice,d2i)*dense(tice,0d0,d2i)            !J/m^3*K, specific heat of ice

                rhotot = (tdens*hsaccum + dense(tice,0d0,d2i)*hi        &
                              + pdens*newsd)/(hsaccum + hi + newsd)
                km = (hi + hms)/(hi/kmi + hsaccum/kms + newsd/kmsn)
                sphm = (sphmi*dense(tice,0d0,d2i) + sphms*tdens         &
                              + sphmsn*pdens)/rhotot                     !J/m^3*K, snow/ice layer
              end if
            end if
          else   !no ice layer
            if(hsaccum > eps.and.dabs(newsd) <= eps) then
              rhotot = tdens
              km = kms
              sphm = sphms
            else if(dabs(hsaccum) <= eps.and.newsd > eps) then
              rhotot = pdens
              km = kmsn
              sphm = sphmsn
            else
              rhotot = (pdens*newsd + tdens*hsaccum)/(newsd + hsaccum)
              km = (hsaccum + newsd)/(hsaccum/kms + newsd/kmsn)
              sphm = (sphms*tdens + sphmsn*pdens)/rhotot                 !J/m^3*K, snow/newsnow layer
            end if
          end if

          km = anint(km*1d20)*1d-20
          sphm = anint(sphm*1d20)*1d-20

! ******************************************************************************

      end select

      end subroutine th_param
      subroutine soil_tmp(ntemp,iflag,kmo,sphmo,isurfg,isurff,disurfg,  &
                          disurff,disurffg,disurfgf,sigflo,rhotot,airo, &
                          dmet,grthcondo,grspheato,zt,delzt,told,iceo,  &
                          wvco,flowuo,flowlo,fv1o,sinko,sourceo,ph_old, &
                          sm_old,thvc,dthvdt,rhov,rhoda,rhs_errort,     &
                          errort,sert,ii,iter)
      
      use fasst_global

! this subroutine calcultes the soil temperture profile
! this version uses Crank-Nicholson

! average deep earth heat flux (qb) of 75 mW/m^2.
!   ref.: Baxter, D.O. (1997) "A Comparison of deep soil temperature: Tenessee 
!           versus other locations",Trans. of the ASAE, 40(3), p.727-738

! no subroutines called

! uses the function: dense,spheats,head,soilhumid,vap_press

      implicit none

      integer(kind=4),intent(in):: ntemp,ii,iter
      integer(kind=4),intent(out):: iflag
      real(kind=8),intent(in):: kmo,sphmo,isurfg,isurff,disurfg,disurff
      real(kind=8),intent(in):: sigflo,disurffg,disurfgf,rhotot,airo
      real(kind=8),intent(in):: dmet(13),grthcondo(ntot),grspheato(ntot)
      real(kind=8),intent(in):: wvco(ntot),zt(ntot),fv1o(ntot)
      real(kind=8),intent(in):: delzt(ntot),iceo(ntot),flowuo(ntot)
      real(kind=8),intent(in):: flowlo(ntot),told(ntot),sinko(ntot)
      real(kind=8),intent(in):: sourceo(ntot),ph_old(ntot),sm_old(ntot)
      real(kind=8),intent(in):: thvc(ntot),dthvdt(ntot)
      real(kind=8),intent(inout):: rhov(ntot),rhoda(ntot)
      real(kind=8),intent(out):: rhs_errort,errort,sert


! local variables
      integer(kind=4):: i,j,sn,d1i,d2i,d3i,d0i,d4i,k,fflag,cflag,sn2
      integer(kind=4):: sn3,is1,if1,is2,if2
      real(kind=8):: dense,spheats,head,soilhumid,f2,f3,f3o
      real(kind=8):: kthu,kthl,kthuo,kthlo,bet,f1,rhsv,sph,out,in,lhe
      real(kind=8):: rhow,rhowo,rhsi,rhsl,rhsu,qb,q1rhs,q1lhsa,q1lhsb
      real(kind=8):: q1lhsc,q2rhs,q2lhsa,q2lhsb,q2lhsc,rhsDv,rhsDi,tf2
      real(kind=8):: rhsDl,W,Wo,didT,rhsw,rhoi,ftempt
      real(kind=8):: maxwat,maxice,rhv,rhovs,pres,rh,dtop,maxdel,tsign
      real(kind=8):: deep,mxd1,d1,d2,t1,t2,t3,dfwi,dfiw,y,fs,minice
      real(kind=8):: A(ntot),B(ntot),C(ntot),D(ntot),lhs1
      real(kind=8):: lhs2,gam(ntot),delstt(ntot),tcheck(ntot)
      real(kind=8):: icheck(ntot),wcheck(ntot),scheck(ntot)
      real(kind=8):: pcheck(ntot),lhs3,lhs4(ntot),lhs5(ntot)
      real(kind=8):: smm,icm,stc,dzl,dzu,dhdT,dhdTo,rhovo,vap_press
      real(kind=8):: rho,rhovi,tl,tu,tlo,tuo,po,vpress,mixr,dtopi,tmp !,qm

      real(kind=8),parameter:: sphveg = 3.5d3                            !average vegetation specific heat (J/kg*K)


! zero-out variables
      sn = 0
      sn2 = 0
      sn3 = 0
      d1i = 0
      d2i = 0
      d3i = 0
      d4i = 0
      is1 = 0
      if1 = 0
      is2 = 0
      if2 = 0
      fflag = 0
      iflag = 0
      cflag = 0
      kthu = 0d0
      kthl = 0d0
      kthuo = 0d0
      kthlo = 0d0
      bet = 0d0
      f1 = 0d0
      rhsv = 0d0
      sph = 0d0
      out = 0d0
      in  = 0d0
      lhe = 0d0
      rhow = 0d0
      rhowo = 0d0
      rhsi = 0d0
      rhsl = 0d0
      rhsu  = 0d0
      qb = 0d0
      q1rhs = 0d0
      q1lhsa = 0d0
      q1lhsb = 0d0
      q1lhsc = 0d0
      q2rhs = 0d0
      q2lhsa = 0d0
      q2lhsb = 0d0
      q2lhsc  = 0d0
      rhsDv = 0d0
      rhsDi = 0d0
      rhsDl = 0d0
      W = 0d0
      Wo = 0d0
      didT = 0d0
      rhsw = 0d0
      rhoi = 0d0
      ftempt = 0d0
      maxwat = 0d0
      maxice = 0d0
      rhv = 0d0
      rhovs = 0d0
      pres = 0d0
      dtop = 0d0
      maxdel = 0d0
      tsign = 0d0
      deep = 0d0
      mxd1 = 0d0
      d1 = 0d0
      d2 = 0d0
      t1 = 0d0
      t2 = 0d0
      t3 = 0d0
      dfwi = 0d0
      dfiw = 0d0
      y = 0d0
      fs = 0d0
      tf2 = 0d0
      minice = 0d0
      f2 = 0d0
      f3 = 0d0
      f3o = 0d0
      smm = 0d0
      icm = 0d0
      stc = 0d0
      dzl = 0d0
      dzu = 0d0
      dhdT = 0d0
      dhdTo = 0d0
      rhovo = 0d0
      rhovi = 0d0
      rho = 0d0
      rh = 0d0
      po = 0d0
      tu = 0d0
      tl = 0d0
      tuo = 0d0
      tlo = 0d0
      vpress = 0d0
      mixr = 0d0
      lhs1 = 0d0
      lhs2 = 0d0
      lhs3 = 0d0
!      qm = 0d0                                                           !constant heat source in profile

! zero out arrays
      do i=1,ntot
        delstt(i) = 0d0
        A(i) = 0d0
        B(i) = 0d0
        C(i) = 0d0
        D(i) = 0d0
        lhs4(i) = 0d0
        lhs5(i) = 0d0
        gam(i) = 0d0
        tcheck(i) = 0d0
        icheck(i) = 0d0
        wcheck(i) = 0d0
        scheck(i) = 0d0
        pcheck(i) = 0d0
      end do

      d0i = 0   !water vapor
      d1i = 1   !water
      d2i = 2   !ice
      d3i = 3   !air
      d4i = 4   !snow

      fs = 1d0/float(step)
      tf2 = lhfus
      if(meltfl == 's') tf2 = lhsub

      f3 = 1d0 - sigfl
      f3o = 1d0 - sigflo

! solve for parameters needed in the tridiagonal matrix
      rhs_errort = -dabs(mflag)
      errort = -dabs(mflag)

      do i=1,ntemp
        if(ntype(i) /= 27) then
          q1lhsa = 0d0
          q1lhsb = 0d0
          q1lhsc = 0d0
          q1rhs = 0d0
          q2lhsb = 0d0
          q2lhsc = 0d0
          q2rhs = 0d0
          qb = 0d0
          lhs1 = 0d0
          lhs2 = 0d0
          lhs3 = 0d0

          lhe = 2500775.6d0 - 2369.729d0*(stt(i) - Tref)                 !latent heat of evaporation, J/kg = (m/s)^2

          f1 = deltat/(2d0*delzt(i))                                     !s/m

          rhow = dense(stt(i),0d0,d1i)                                   !kg/m^3
          rhowo = dense(told(i),0d0,d1i)                                 !kg/m^3
          rhoi = dense(stt(i),0d0,d2i)

          out = 5d-1*(spheats(stt(i),d1i)*rhow*sink(i)                  &
                                 + spheats(told(i),d1i)*rhowo*sinko(i))  !(J/m^3*K)*m/s
          if(dabs(out) < eps) out = 0d0

          in = 5d-1*(spheats(stt(i),d1i)*rhow*source(i)                 &
                                + spheats(told(i),d1i)*rhowo*sourceo(i)) !(J/m^3*K)*m/s 
          if(dabs(in) < eps) in = 0d0

          rhsDi = 0d0
          rhsi = 0d0
          rhsDl = 0d0
          rhsw = 0d0
          rhsDv = 0d0
          rhsv = 0d0

! calculate phase change energies
          if(i <= nnodes.and.ntype(i) /= 26) then
            sph = 5d-1*(grspheat(i) + grspheato(i))                      !J/m^3*K
            rh = soilhumid(i,phead(i),soil_moist(i),stt(i))
            rho = soilhumid(i,ph_old(i),sm_old(i),told(i))
            pres = dmet(11) + 1d-2*((elev - nz(i)) + phead(i))          &
                                                             *rhow*grav  !mbar
            po = dmet(11) + 1d-2*((elev - nz(i)) + ph_old(i))           &
                                                            *rhowo*grav  !mbar
            dhdT = 0d0
            if(dabs(rh) > eps) dhdT = -(Rv/grav)*dlog(rh)
            dhdTo = 0d0
            if(dabs(rho) > eps) dhdTo = -(Rv/grav)*dlog(rho)
!            dhdT = (7d0/71.89d0)*(-0.1425d0 - 9.52d-4*(stt(i) - Tref))  &
!                                                              *phead(i)
!            dhdTo = (7d0/71.89d0)*(-0.1425d0 - 9.52d-4*(told(i) - Tref))&
!                                                             *ph_old(i)

! NOTE:
! freezing releases heat -> warms soil/water around it
! ice - iceo > 0 -> warming; ice - iceo < 0 -> cooling

!soil
! freezing/thawing energy
            if(stt(i) > Tref) then
              didT = -(rhow*lhfus/(rhoi*grav*Tref))*dsmdh(i)             !1/K
            else
              didT = -(rhow*lhfus/(rhoi*grav*stt(i)))*dsmdh(i)           !1/K
            end if
            didT = didT !*fs

            lhs5(i) = didT                                               !1/K, associated with freezing/thawing

            rhsDi = 0d0
            rhsi = 0d0
            if(ice(i)+iceo(i) > eps) then
!              if(stt(i) < Tref) then
                rhsDi = 5d-1*((lhfus + spheats(stt(i),d2i)              &
                                            *dabs(stt(i) - Tref))*rhoi  &
                      + (lhfus + spheats(told(i),d2i)                   &
                          *dabs(told(i) - Tref))*dense(told(i),0d0,d2i)) !J/m^3
!              else
!                rhsDi = 5d-1*lhfus*(rhoi + dense(told(i),0d0,d2i))
!              end if
              rhsDi = rhsDi
              rhsi = -rhsDi*(ice(i) - iceo(i))                           !J/m^3
              lhs1 = lhs5(i)
            end if

! water warming/cooling
            if(soil_moist(i)+sm_old(i) > eps) then
              t1 = 5d-1*(spheats(stt(i),d1i)*dabs(stt(i) - Tref)*rhow   &
                            + spheats(told(i),d1i)*dabs(told(i) - Tref) &
                                                                 *rhowo) !J/m^3
              t1 = anint(t1*1d15)*1d-15

              Wo = grav*(ph_old(i) - told(i)*dhdTo)
              W = grav*(phead(i) - stt(i)*dhdT)                          !(m/s)^2 = J/kg
              t2 = 5d-1*(W*rhow + Wo*rhowo)                              !J/m^3
              rhsDl = anint((t1 + t2)*1d15)*1d-15
              rhsw = rhsDl*(soil_moist(i) - sm_old(i))                   !J/m^3
              lhs3 = dhdT*dsmdh(i)
            end if

! evaporation/condensation energy
            rhsDv = 0d0
            rhsv = 0d0
            if(wvc(i)+wvco(i) > eps) then
              rhovo = 0.622d0*vap_press(i,rho,po)/(Rv*told(i))
              rhovi = 0.622d0*vap_press(i,rh,pres)/(Rv*stt(i))
              rhsDv = 5d-1*(lhe + spheats(stt(i),d0i)                   &
                                            *dabs(stt(i) - Tref))*rhovi &
                     + (lhe + spheats(told(i),d0i)                      &
                                           *dabs(told(i) - Tref))*rhovo  !J/m^3
!rhsDv = lhe*5d-1*(rhovi + rhovo)
              rhsDv = rhsDv
              rhsv = rhsDv*(wvc(i) - wvco(i))                            !J/m^3
            end if

            t1 = dthvdt(i)*(nsoilp(i,2) - (soil_moist(i) + ice(i)))
            t2 = -thvc(i)*(dhdT + lhs1)*dsmdh(i)

            lhs4(i) =  t1 + t2                                           !1/K, associated with evaporation/condensation
            if(wvc(i)+wvco(i) > eps) lhs2 = lhs4(i)

          else if(i <= nnodes.and.node_type(i) == 'WA') then            
!open water
            sph = 5d-1*(grspheat(i) + grspheato(i))                      !J/m^3*K

            t1 = 5d-1*(spheats(stt(i),d1i)*dabs(stt(i) - Tref)*rhow     &
                      + spheats(told(i),d1i)*dabs(told(i) - Tref)*rhowo) !J/m^3
            rhsDl = t1
            rhsw = rhsDl*(soil_moist(i) - sm_old(i))                     !J/m^3
            lhs3 = 0d0

            didT = -(rhow/rhoi)*lhfus/(timstep*3.6d3)
            lhs5(i) = didT                                               !1/K, associated with freezing/thawing
            if(ice(i) > eps.or.iceo(i) > eps) then
!              if(stt(i) <= Tref) then
                rhsDi = 5d-1*((lhfus + spheats(stt(i),d2i)              &
                                            *dabs(stt(i) - Tref))*rhoi  &
                                        + (lhfus + spheats(told(i),d2i) &
                          *dabs(told(i) - Tref))*dense(told(i),0d0,d2i)) !J/m^3
!              else
!                rhsDi = 5d-1*lhfus*(rhoi + dense(told(i),0d0,d2i))
!              end if
              rhsi = rhsDi*(ice(i) - iceo(i))                            !J/m^3
              lhs1 = lhs5(i)
            end if

          else if(i > nnodes) then
!snow
! snow freezing/thawing energy
            if(node_type(i) == 'HM') then
!              if(stt(i) < Tref) then
                rhsDi = (lhfus*rhotot + sphm*dabs(stt(i) - Tref))        !J/m^3
!              else
!                rhsDi = lhfus*rhotot
!              end if

              if(met(iw,ip_sd) > eps) then
                rhsi = rhsDi*(-refreeze/met(iw,ip_sd))                   !J/m^3
                didT = tf2*(-refreeze/met(iw,ip_sd))                    &
                                            /(met(iw,ip_sd)*grav*Tref)   !1/K
              else if(hm > eps) then
                rhsi = rhsDi*(-refreeze/hm)                              !J/m^3
                didT = tf2*(-refreeze/hm)/(hm*grav*Tref)                 !1/K
              end if

              rhsi = rhsi*fs
              didT = didT*fs
              lhs1 = didT                                                !1/K, associated with freezing/thawing

! snow evaporation/condensation energy
              if(wvc(i) > eps.or.wvco(i) > eps) then
                rhsDv = 0d0
                rhsv = 0d0
!                if(stt(i) < Tref) then
!                  rhsDv = (lhe*rhotot + sphm*(stt(i) - Tref))            !J/m^3
!                else
!                  rhsDv = lhe*rhotot
!                end if
                rhsDv = (lhe*rhotot + sphm*dabs(stt(i) - Tref))          !J/m^3
!                rhsv = rhsDv*(wvc(i) - wvco(i))*fs                       !J/m^3

                lhs2 = -(rhov(i)/rhow)*didT                              !1/K, associated with evaporation/condensation
              end if
            end if
          end if

          lhs1 = anint(lhs1*1d20)*1d-20
          lhs2 = anint(lhs2*1d20)*1d-20
          lhs3 = anint(lhs3*1d20)*1d-20
          lhs4(i) = anint(lhs4(i)*1d20)*1d-20
          lhs5(i) = anint(lhs5(i)*1d20)*1d-20

! interior and bottom nodes
          if(i <= nnodes-1) then
            if(i == 1) then                                             !bottom node
              kthl = 0d0                                                 !W/m*K
              kthlo = 0d0                                                !W/m*K

              if(node_type(1) /= 'SN') qb = 2d0*7.5d-2                   !W/m^2, constant deep-earth heat flux
            else                                                        !interior nodes
              kthl = (grthcond(i-1)*delzt(i-1) + grthcond(i)*delzt(i))  &
                                               /(delzt(i-1) + delzt(i))  !W/m*K
              kthlo = (grthcondo(i-1)*delzt(i-1)                        &
                       + grthcondo(i)*delzt(i))/(delzt(i-1) + delzt(i))  !W/m*K
            end if

            kthu = (grthcond(i+1)*delzt(i+1) + grthcond(i)*delzt(i))    &
                                               /(delzt(i+1) + delzt(i))  !W/m*K
            kthuo = (grthcondo(i+1)*delzt(i+1) + grthcondo(i)*delzt(i)) &
                                               /(delzt(i+1) + delzt(i))  !W/m*K

!            if(i == ?)
!              qm = ??
!            else
!              qm = 0d0
!            end if

! top soil node, veg & snow nodes if present
          else

! no veg, no snow
            if(icase == 0) then
              sn = i

              kthl = (grthcond(i-1)*delzt(i-1) + grthcond(i)*delzt(i))  &
                                               /(delzt(i-1) + delzt(i))  !W/m*K
              kthlo = (grthcondo(i-1)*delzt(i-1)                        &
                       + grthcondo(i)*delzt(i))/(delzt(i-1) + delzt(i))  !W/m*K

              kthu = 0d0
              kthuo = 0d0

              q1lhsb = disurfg                                            !W/m^2*K
              q1rhs = (isurfoldg + isurfg)                                !W/m^2

! snow/no veg .or. veg buried by snow
            else if(icase == 1.or.icase == 4) then
              if(i == ntemp) then                                       !snow layer
                sn = i

                kthl = (grthcond(i-1)*delzt(i-1)                        &
                                       + (f3*km + sigfl*kveg)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K
                kthlo = (grthcondo(i-1)*delzt(i-1)                      &
                                    + (f3o*kmo + sigflo*kveg)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K

                kthu = 0d0
                kthuo = 0d0

                out = 0d0
                in = 0d0
                sph = 5d-1*(sphm*f3 + sphmo*f3o)                          &
                              + 5d-1*sphveg*(sigfl*rhow + sigflo*rhowo)  !J/m^3*K

                q1lhsb = disurfg                                         !W/m^2*K
                q1rhs = (isurfoldg + isurfg)                             !W/m^2
              else                                                      !ground layer
                kthl = (grthcond(i-1)*delzt(i-1)                        &
                        + grthcond(i)*delzt(i))/(delzt(i-1) + delzt(i))  !W/m*K
                kthlo = (grthcondo(i-1)*delzt(i-1)                      &
                                               + grthcondo(i)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K

                kthu = (grthcond(i)*delzt(i)                            &
                                     + (f3*km + sigfl*kveg)*delzt(i+1)) &
                                               /(delzt(i) + delzt(i+1))  !W/m*K
                kthuo = (grthcondo(i)*delzt(i)                          &
                                  + (f3o*kmo + sigflo*kveg)*delzt(i+1)) &
                                               /(delzt(i) + delzt(i+1))  !W/m*K

                q1lhsb = 0d0
                q1rhs = 0d0
              end if

! vegetation, no snow
            else if(icase == 2) then
              if(i == ntemp) then                                       !veg layer
!              sn = i

                kthl = (grthcond(i-1)*delzt(i-1)                        &
                                                 + sigfl*kveg*delzt(i)) &
                                        / (delzt(i-1) + sigfl*delzt(i))  !W/m*K
                kthlo = (grthcondo(i-1)*delzt(i-1)                      &
                                                + sigflo*kveg*delzt(i)) &
                                        /(delzt(i-1) + sigflo*delzt(i))  !W/m*K

                kthu = 0d0
                kthu = 0d0

                out = 0d0
                in = 0d0

                sph = 5d-1*sphveg*(rhow + rhowo)                         !J/m^3*K
!                sph = 5d-1*sphveg*(sigfl*rhow + sigflo*rhowo)            !J/m^3*K

                q1lhsb = disurff                                         !W/m^2*K (veg)
                q2lhsb = -disurfgf                                        !W/m^2*K (ground)

                q1lhsa = disurffg                                        !W/m^2*K (veg)
                q2lhsa = -disurfg                                         !W/m^2*K (ground)

                q1rhs = (isurfoldf + isurff)                             !W/m^2
                q2rhs = -(isurfoldg + isurfg)                             !W/m^2
              else                                                      !ground layer
                sn = i

                kthl = (grthcond(i-1)*delzt(i-1)                        &
                                                + grthcond(i)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K
                kthlo = (grthcondo(i-1)*delzt(i-1)                      &
                                               + grthcondo(i)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K

                kthu = sigfl*(grthcond(i)*delzt(i) + kveg*delzt(i+1))   &
                                               /(delzt(i) + delzt(i+1))  !W/m*K
                kthuo = sigflo*(grthcondo(i)*delzt(i)                   &
                                                     + kveg*delzt(i+1)) &
                                               /(delzt(i) + delzt(i+1))  !W/m*K


                q1lhsb = disurfg                                         !W/m^2 (ground)
                q1lhsc = disurfgf                                        !W/m^2 (veg)
                q1rhs = (isurfoldg + isurfg)                             !W/m^2
              end if

! vegetation, snow (ntemp = nnodes + 2)
            else if(icase == 3) then
              if(i == ntemp) then                                       !veg layer
!              sn = i

                kthl = ((f3*km + sigfl*kveg)*delzt(i-1)                 &
                                                 + sigfl*kveg*delzt(i)) &
                                         /(delzt(i-1) + sigfl*delzt(i))  !W/m*K
                kthlo = ((f3o*kmo + sigflo*kveg)*delzt(i-1)             &
                                                + sigflo*kveg*delzt(i)) &
                                        /(delzt(i-1) + sigflo*delzt(i))  !W/m*K

                kthu = 0d0
                kthuo = 0d0

                out = 0d0
                in = 0d0

                sph = 5d-1*sphveg*(rhow + rhowo)                         !J/m^3*K
!                sph = 5d-1*sphveg*(sigfl*rhow + sigflo*rhowo)            !J/m^3*K

                q1lhsb = disurff                                         !W/m^2*K (veg)
                q2lhsb = -disurfgf                                        !W/m^2*K (ground)

                q1lhsa = disurffg                                        !W/m^2*K (veg)
                q2lhsa = -disurfg                                         !W/m^2*K (ground)

                q1rhs = (isurfoldf + isurff)                             !W/m^2
                q2rhs = -(isurfoldg + isurfg)                             !W/m^2
              else if(i == ntemp-1) then                                !snow layer
                sn = i

                kthl = (grthcond(i-1)*delzt(i-1)                        &
                                       + (f3*km + sigfl*kveg)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K
                kthlo = (grthcondo(i-1)*delzt(i-1)                      &
                                    + (f3o*kmo + sigflo*kveg)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K

                kthu = sigfl*((f3*km + sigfl*kveg)*delzt(i)             &
                                                     + kveg*delzt(i+1)) &
                                               /(delzt(i) + delzt(i+1))  !W/m*K
                kthuo = sigflo*((f3o*kmo + sigflo*kveg)*delzt(i)        &
                                                     + kveg*delzt(i+1)) &
                                               /(delzt(i) + delzt(i+1))  !W/m*K

                out = 0d0
                in = 0d0
                sph = 5d-1*(sphm*f3 + sphmo*f3o)                        &
                              + 5d-1*sphveg*(sigfl*rhow + sigflo*rhowo)  !J/m^3*K

                q1lhsb = disurfg                                         !W/m^2 (ground)
                q1lhsc = disurfgf                                        !W/m^2 (veg)
                q1rhs = (isurfoldg + isurfg)                             !W/m^2
              else                                                      !ground layer
                kthl = (grthcond(i-1)*delzt(i-1)                        &
                                                + grthcond(i)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K
                kthlo = (grthcondo(i-1)*delzt(i-1)                      &
                                               + grthcondo(i)*delzt(i)) &
                                               /(delzt(i-1) + delzt(i))  !W/m*K

                kthu = (grthcond(i)*delzt(i)                            &
                                     + (f3*km + sigfl*kveg)*delzt(i+1)) &
                                               /(delzt(i) + delzt(i+1))  !W/m*K
                kthuo = (grthcondo(i)*delzt(i)                          &
                                  + (f3o*kmo + sigflo*kveg)*delzt(i+1)) &
                                               /(delzt(i) + delzt(i+1))  !W/m*K
              end if
            end if  !interior nodes; veg vs no veg, snow vs no snow....
          end if  !partition between nodes

          if(i == 1) then
            dzl = 0d0
            dzu = 1d0/(zt(i+1) - zt(i))
            tl = stt(i)
            tlo = told(i)
            tu = 5d-1*(stt(i) + stt(i+1))
            tuo = 5d-1*(told(i) + told(i+1))
            rhsl = qb                                                    !W/m^2
            rhsu = kthu*(stt(i+1) - stt(i))*dzu                         &
                                      + kthuo*(told(i+1) - told(i))*dzu  !W/m^2
          else if(i <= nnodes-1.and.i > 1) then
            dzl = 1d0/(zt(i) - zt(i-1))
            dzu = 1d0/(zt(i+1) - zt(i))
            tl = 5d-1*(stt(i) + stt(i-1))
            tlo = 5d-1*(told(i) + told(i-1))
            tu = 5d-1*(stt(i) + stt(i+1))
            tuo = 5d-1*(told(i) + told(i+1))
            rhsl = kthl*(stt(i) - stt(i-1))*dzl                         &
                                      + kthlo*(told(i) - told(i-1))*dzl  !W/m^2  + 5d-1*qm
            rhsu = kthu*(stt(i+1) - stt(i))*dzu                         &
                                 + kthuo*(told(i+1) - told(i))*dzu       !W/m^2  + 5d-1*qm
          else
!          dzl = -1d0/dmax1(1d-2,dabs(zt(i-1) - zt(i)))
            dzl = 1d0/(zt(i) - zt(i-1))
            dzu = 0d0
            tl = 5d-1*(stt(i) + stt(i-1))
            tlo = 5d-1*(told(i) + told(i-1))
            rhsl = kthl*(stt(i) - stt(i-1))*dzl                         &
                              + kthlo*(told(i) - told(i-1))*dzl + q2rhs  !W/m^2
            if(i == ntemp) then
              tu = 5d-1*(stt(i) + dmet(4))
              tuo = 5d-1*(told(i) + airo)
              rhsu = q1rhs                                               !W/m^2
            else
              tu = 5d-1*(stt(i) + stt(i+1))
              tuo = 5d-1*(told(i) + told(i+1))
!            dzu = -1d0/dmax1(1d-2,dabs(zt(i) - zt(i+1)))
              dzu = 1d0/(zt(i+1) - zt(i))
              rhsu = kthu*(stt(i+1) - stt(i))*dzu                       &
                              + kthuo*(told(i+1) - told(i))*dzu + q1rhs  !W/m^2
            end if
          end if

! A(i) = coeff. for delstt(i-1)
! B(i) = coeff. for delstt(i)
! C(i) = coeff. for delstt(i+1)
! D(i) = rhs
! z is positive upwards from the bottom
          if(i == 1) then
            B(i) =  f1*kthu/(zt(i+1) - zt(i)) - f1*(in - out) + sph     &
                     + f1*5d-1*(flowu(i) - flowl(i))                    &
                     - rhsDi*lhs1 + rhsDv*lhs2 + rhsDl*lhs3              !J/m^3*K

            C(i) = -f1*kthu/(zt(i+1) - zt(i)) - f1*5d-1*flowu(i)         !J/m^3*K
          else if (i == ntemp) then
            A(i) = -f1*kthl/(zt(i) - zt(i-1)) - f1*5d-1*flowl(i)        &
                     + f1*(q2lhsa - q1lhsa)                              !J/m^3*K

            B(i) =  f1*kthl/(zt(i) - zt(i-1)) - f1*(in - out) + sph     &
                     + f1*5d-1*(flowu(i) - flowl(i))                    &
                     - rhsDi*lhs1 + rhsDv*lhs2 + rhsDl*lhs3             &
                     + f1*(q1lhsb - q2lhsb)                              !J/m^3*K
          else
            A(i) = -f1*kthl/(zt(i) - zt(i-1)) - f1*5d-1*flowl(i)         !J/m^3*K

            B(i) =  f1*(kthl/(zt(i) - zt(i-1)) + kthu/(zt(i+1) - zt(i)))&
                     + f1*5d-1*(flowu(i) - flowl(i)) - f1*(in - out)    &
                     - rhsDi*lhs1 + rhsDv*lhs2 + rhsDl*lhs3             &
                     - f1*q1lhsb + sph                                   !J/m^3*K

            C(i) = -f1*kthu/(zt(i+1) - zt(i)) - f1*5d-1*flowu(i)        &
                     + f1*(q1lhsc - q2lhsc)                              !J/m^3*K
          end if

          D(i) = -sph*(stt(i) - told(i)) + rhsi - rhsv - rhsw           &
                   + f1*(rhsu - rhsl) - f1*(fv1(i) + fv1o(i))           &
                   + f1*(in - out)*(stt(i) + told(i))                   &
                   - f1*(flowu(i)*dabs(tu - Tref)                       &
                   - flowl(i)*dabs(tl - Tref))                          &
                   - f1*(flowuo(i)*dabs(tuo - Tref)                     &
                   - flowlo(i)*dabs(tlo - Tref))                         !J/m^3

          A(i) = anint(A(i)*1d15)*1d-15
          B(i) = anint(B(i)*1d15)*1d-15
          C(i) = anint(C(i)*1d15)*1d-15
          D(i) = anint(D(i)*1d15)*1d-15

          if(dabs(D(i)) > rhs_errort) then
            rhs_errort = dabs(D(i))
            sert = D(i)
          end if

        else if(ntype(i) == 27) then
          D(i) = 0d0
          if(sn2 == 0) sn2 = i
          sn3 = i
        end if  !ntype(i) /= 27
      end do  !i=1,ntemp
      rhs_errort = anint(rhs_errort*1d10)*1d-10

! solve the matrix eqn. for delta sst(i)
      is2 = 0
      if2 = 0
      if(sn2 == 0.and.sn3 == 0) then                                     !no air nodes
        is1 = 1
        if1 = ntemp
      else if(sn2 == 1.and.sn3 == ntemp) then                            !all air nodes
        is1 = 0
        if1 = 0
      else if(sn2 > 1.and.sn3 == ntemp) then                             !air on top
        is1 = 1
        if1 = sn2 - 1
      else if(sn2 == 1.and.sn3 < ntemp) then                             !air on bottom
        is1 = sn3 + 1
        if1 = ntemp
      else if(sn2 > 1.and.sn3 < ntemp) then                              !air layer
        is1 = 1
        if1 = sn2 - 1
        is2 = sn3 + 1
        if2 = ntemp
      end if

      if(is1 > 0.and.if1 > 0) then
        bet = B(is1)
        if(dabs(bet) <= 1d-10) then
          tsign = 1d0
          if(bet < 0d0) tsign = -1d0                                     !tridiag fails
          delstt(is1) = tsign*dabs(D(is1))
        else
          delstt(is1) = D(is1)/bet                                       !K
        end if

        do j=is1+1,if1
          if(ntype(j) /= 27) then
            gam(j) = C(j-1)/bet                                          !unitless
            bet = B(j) - A(j)*gam(j)                                     !1/K
            if(dabs(bet) <= 1d-10.or.dabs(A(j)) <= eps) then             !tridiag fails
              tsign = 1d0
              if(bet < 0d0) tsign = -1d0
              delstt(j) = tsign*dabs(D(j))
            else
              delstt(j) = (D(j) - A(j)*delstt(j-1))/bet                  !K
            end if
          else
            delstt(j) = 0d0
          end if
        end do

        do j=if1-1,is1,-1
          delstt(j) = delstt(j) - gam(j+1)*delstt(j+1)                   !K
        end do
      end if

      if(is2 > 0.and.if2 > 0) then
        bet = B(is2)
        if(dabs(bet) <= 1d-10) then
          tsign = 1d0
          if(bet < 0d0) tsign = -1d0                                     !tridiag fails
          delstt(is2) = tsign*dabs(D(is2))
        else
          delstt(is2) = D(is2)/bet                                       !K
        end if

        do j=is2+1,if2
          if(ntype(j) /= 27) then
            gam(j) = C(j-1)/bet                                          !unitless
            bet = B(j) - A(j)*gam(j)                                     !1/K
            if(dabs(bet) <= 1d-10.or.dabs(A(j)) <= eps) then             !tridiag fails
              tsign = 1d0
              if(bet < 0d0) tsign = -1d0
              delstt(j) = tsign*dabs(D(j))
            else
              delstt(j) = (D(j) - A(j)*delstt(j-1))/bet                  !K
            end if
          else
            delstt(j) = 0d0
          end if
        end do

        do j=if2-1,is2,-1
          delstt(j) = delstt(j) - gam(j+1)*delstt(j+1)                   !K
        end do
      end if

! determine new temps, ice and vapor contents; adjust water contents accordingly
      errort = -99999.9d0
      dtop = 2.5d0
      if(iter == 0) dtop = 1d-1*dtop !2.5d-1
      mxd1 = 4d1*(nz(nnodes) - nz(1))
      do i=ntemp,1,-1
        delstt(i) = anint(delstt(i)*1d15)*1d-15

        tcheck(i) = stt(i)
        wcheck(i) = wvc(i)
        icheck(i) = ice(i)
        scheck(i) = soil_moist(i)
        pcheck(i) = phead(i)

        if(dabs(delstt(i)) > errort.and.i <= nnodes)                    &
                                               errort = dabs(delstt(i))

        if(iter > 0.and.(dabs(ice(i)-iceo(i)) > eps.and.                &
                                                     i <= nnodes)) then
          dtopi = 1d-1*dtop
        else
          dtopi = dtop
        end if

        if(dabs(delstt(i)) > dtopi) then
          tsign = 0d0
          if(delstt(i) /= 0d0) then
            tsign = dabs(delstt(i))/delstt(i)
          end if
          delstt(i) = tsign*dtopi
          delstt(i) = anint(delstt(i)*1d15)*1d-15
        end if

        stt(i) = stt(i) + delstt(i)

        if(i == nnodes.and.mstflag == 1) then                            !use measured soil temp
          if(aint(dabs(dmet(13)-mflag)*1d5)*1d-5 > eps) then
            stt(nnodes) = dmet(13)
            fflag = 1                                                    !must use this temp
          end if
        end if
!!      if(i == 1) stt(i) = Fuquan_value !!!

        if(ntype(i) == 27) then
          if(i > nnodes.or.(ntype(nnodes) == 27.or.ntype(1) == 27)) then
            stt(i) = dmet(4)                                             !air
          else
            if(i < ntemp.and.i > 1) then
              stt(i) = 5d-1*(stt(sn3+1) + stt(sn2-1))
            else
              stt(i) = dmet(4)
            end if
          end if
        end if

        maxdel = 1d1*timstep

        if(iw /= 1) then
          if(dabs(dmet(4)-airo) > eps) then
            maxdel = dmin1(mxd1,1d1*dabs(dmet(4)-airo))
          end if
        end if

        d1 = 0d0
        if(i < nnodes)  then
          d1 = dmin1(1d0,elev - nz(i))
        else if(i >= nnodes.and.ntemp /= nnodes) then
          d1 = dmax1(1d0,hfol,hm)
        end if

        deep = 1d0
        if(d1 > eps) then
          deep = 4.388d0*(d1**4d0) - 11.498d0*(d1**3d0)                 &
                              + 10.947d0*d1*d1 - 4.7845d0*d1 + 9.546d-1
          deep = dmin1(1d0,dmax1(1d-1,deep))
        end if

        maxdel = dmax1(maxdel*deep,mxd1*deep)
!        if(node_type(i) == 'VG') maxdel = maxdel*5d-1

        if(elev-nz(i) >= 1d0) then
          y = 1d0 + dmin1(1d0,1d-2*((elev - nz(i))-1d0))
        else if(i >=nnodes.and.hm >= 1d0) then
          y = 1d0 + dmin1(1d0,1d-2*(hm-1d0))
        else if(i >=nnodes.and.hfol >= 1d0) then
          y = 1d0 + dmin1(1d0,1d-2*(hfol-1d0))
        else
          y = 1d0 - deep
        end if

        y = 6d1 !dmin1(dmax1(maxdel,4d1*y),6d1)

! check for continuity in time
        if(dabs(stt(i)-too(i)) > maxdel) then
          if(stt(i) > too(i)) then
            stt(i) = too(i) + maxdel
          else
            stt(i) = too(i) - maxdel
          end if
        end if

! check for vertical continuity
        if(ntype(i) /= 27) then                                          !not air
          if(i == ntemp.or.(i /= ntemp.and.cflag == 0)) then
            cflag = 1
            if(dabs(dmet(4)-stt(i)) > 4d1) then
              if(stt(i) > dmet(4)) then
                j = 1
                do while(stt(i) > dmet(4)+4d1)
                  stt(i) = tcheck(i) - j*maxdel
                  j = j + 1
                end do
              else 
                j = 1
                do while(stt(i) < dmet(4)-4d1)
                  stt(i) = tcheck(i) + j*maxdel
                  j = j + 1
                end do
              end if
            end if

            if(i == nnodes.and.fflag == 1) stt(i) = dmet(13)             !must use measured temp
            tmelt(iw) = stt(i)
!!      if(i == 1) stt(i) = Fuquan_value !!!
          else
            k = sn
            if(sn2 > 0) k = sn2 - 1

            if(dabs(stt(k)-stt(i)) > y.and.k > 1) then
              if(stt(i) > stt(k)) then
                stt(i) = stt(k) + y
              else
                stt(i) = stt(k) - y
              end if
            end if

            if(i == nnodes.and.fflag == 1) stt(i) = dmet(13)             !must use measured temp
            if(icase == 3.and.i == sn) tmelt(iw) = stt(i)
!!      if(i == 1) stt(i) = Fuquan_value !!!
          end if
        end if

        stt(i) = anint(stt(i)*1d10)*1d-10

        if(node_type(i) == 'VG') then
          ftempt = stt(i)
        else
          ftempt = dmet(4)
        end if

        if(stt(i) <= Tref) then
          if(told(i) <= Tref) then
            rhow = dense(273.16d0,0d0,d1i)
          else
            rhow = dense(told(i),0d0,d1i)
          end if
          rhoi = dense(stt(i),0d0,d2i)
        else
          rhow = dense(stt(i),0d0,d1i)
          if(told(i) <= Tref) then
            rhoi = dense(told(i),0d0,d2i)
          else
            rhoi = dense(Tref,0d0,d2i)
          end if
        end if
        dfwi = rhow/rhoi                                                 !water/ice
        dfiw = rhoi/rhow                                                 !ice/water

! update ice content
        if(i <= nnodes.and.ntype(i) <= 25) then
          smm = soil_moist(i)*rhow                                       !kg
          icm = ice(i)*rhoi                                              !kg
          stc = anint((stt(i) - Tref)*1d15)*1d-15

          if(ii == 1.and.iter == 0) bftm(i) = smm + icm                   !kg
          bftm(i) = anint(bftm(i)*1d15)*1d-15

          if(stc < eps.or.ice(i) > 0d0) then
            minice = 0d0
            maxice = dmax1(0d0,dmin1(icm + smm,bftm(i)))                 !kg

            if(dabs(stc) > eps) then
              y = dabs(lhs5(i))*dmax1(dabs(delstt(i)),dabs(stc))         !m^3   ! *nsoilp(i,2)
            else
              y =  dabs(lhs5(i)*delstt(i))
            end if

            if((dabs(lhs5(i)) <= eps.and.ice(i) > eps)) then
              didT = -(rhow/rhoi)*lhfus/(grav*Tref*dabs(delzt(i)))
              if(dabs(stc) > eps) then
                y = dabs(didT)*dmax1(dabs(delstt(i)),dabs(stc))         !m^3  !*nsoilp(i,2)
              else
                y = dabs(didT*delstt(i))
              end if
            end if

            if(node_type(i) == 'WA') then
              if(dabs(stc) > eps) then
                y = dabs(lhs5(i))*dmin1(dabs(delstt(i)),dabs(stc))
              else
                y = dabs(lhs5(i)*delstt(i))
              end if
            end if
 
            if(y < eps) y = 0d0
            y = anint(y*1d15)*1d-15

            if(stc >= eps) then   !.or.((i == nnodes.and.hm <= eps).and.          &
!                                              dmet(4) > Tref+2d0)) then  !melting
              icm = (icheck(i) - y)*rhoi                                 !kg
              if(icm < minice) then
                y = icheck(i)*rhoi - minice
                icm = minice
              end if
              smm = scheck(i)*rhow + y*rhoi                              !kg

            else if(stc < eps.and.soil_moist(i) > eps) then              !freezing
              icm = (icheck(i) + y)*rhoi                                 !kg

              if(icm > maxice) then
                y = icm - maxice
                icm = maxice
              end if
              smm = dmax1(0d0,scheck(i)*rhow - y*rhoi)                   !kg

            else if(stc < eps.and.node_type(i) == 'WA') then
              icm = (icheck(i) + y)*rhoi                                 !kg
              if(icm > maxice) then
                y = icm - maxice
                icm = maxice
              end if
              smm = dmax1(0d0,scheck(i)*rhow - y*rhoi)                   !kg
            end if

            icm = anint(icm*1d15)*1d-15                                  !kg
            ice(i) = dmax1(0d0,icm/rhoi)
            ice(i) = anint(ice(i)*1d10)*1d-10
            smm = anint(smm*1d15)*1d-15                                  !kg
            soil_moist(i) = dmax1(0d0,smm/rhow)
            soil_moist(i) = anint(soil_moist(i)*1d10)*1d-10

!            if(dabs(ice(i)-icheck(i)) > 1d-8) then
!!            if(dabs(ice(i)-icheck(i)) > 1d-4*bftm(i)/rhow.and.soil_moist(i) > 0d0)   &
            if((icheck(i) <= eps.and.ice(i) > eps).or.                  &
                              (icheck(i) > eps.and.ice(i) <= eps)) then
              if(iw /= istart) stt(i) = Tref
            end if

            if(dabs((smm+icm)-bftm(i)) > 1d-5.and.                      &
                                             node_type(i) /= 'WA') then
              smm = bftm(i) - icm                                        !kg
              soil_moist(i) = anint((smm/rhow)*1d10)*1d-10

              if(ice(i) > eps) then
                if(soil_moist(i) <= 1d-10) soil_moist(i) = 0d0
              else
                if(soil_moist(i) < nsoilp(i,15)) then
                  soil_moist(i) = nsoilp(i,15)
                else if(soil_moist(i) > nsoilp(i,24)) then
                  soil_moist(i) = nsoilp(i,24)
                end if
              end if
              if(node_type(i) /= 'WA'.and.node_type(i) /= 'AI')         &
                                       phead(i) = head(i,soil_moist(i))

            else if(ice(i) <= eps) then
              if(soil_moist(i) < nsoilp(i,15)) then
                soil_moist(i) = nsoilp(i,15)
                if(node_type(i) /= 'WA'.and.node_type(i) /= 'AI')       &
                  phead(i) = head(i,soil_moist(i))
              else if(soil_moist(i) > nsoilp(i,24)) then
                soil_moist(i) = nsoilp(i,24)
                if(node_type(i) /= 'WA'.and.node_type(i) /= 'AI')       &
                  phead(i) = head(i,soil_moist(i))
              else if(node_type(i) == 'WA') then
                soil_moist(i) = 1d0
              end if
            end if

!            if((i == nnodes.and.fflag == 1).and.                        &
!                                      dabs(stt(i)-dmet(13)) > eps) then
!              stt(i) = dmet(13)                                          !must use measured temp
!              if(stc >= eps.and.ice(i) > eps) then
!                soil_moist(i) = soil_moist(i) + ice(i)*dfiw
!                ice(i) = 0d0
!              else if(stc < eps.and.soil_moist(i) > eps) then
!                ice(i) = ice(i) + soil_moist(i)*dfwi
!                soil_moist(i) = 0d0
!              end if
!              if(node_type(i) /= 'WA'.and.node_type(i) /= 'AI')         &
!                phead(i) = head(i,soil_moist(i))
!            end if
            if(ice(i) > eps) iflag = 1
          end if

! update water vapor content
          if((ntype(i) /= 26.or.ntype(i) /= 27)                         &
                                        .and.node_type(i) /= 'VG') then
            rh = soilhumid(i,phead(i),soil_moist(i),stt(i))
            pres = dmet(11) + 1d-2*((elev - nz(i)) + phead(i))          &
                                                             *rhow*grav  !mbar
          else
            rh = dmet(5)*1d-2
            if(ntype(i) == 26) rh = 1d0
            pres = dmet(11)
          end if

          vpress = vap_press(i,rh,pres)                                  !Pa
          mixr = 0d0
          if(pres*1d2-vpress > eps)                                     &
            mixr = 0.622d0*vpress/(pres*1d2 - vpress)                    !unitless (kg/kg)

          rhov(i) = 0.622d0*vpress/(Rv*stt(i))                           !kg/m^3
          rhov(i) = anint(rhov(i)*1d20)*1d-20
          rhoda(i) = dmax1(0.95d0,dmin1(2.8d0,(pres*1d2 - vpress)/      &
                                                          (Rd*stt(i))))  !kg/m^3
          rhoda(i) = anint(rhoda(i)*1d20)*1d-20

          t1 = dmin1(1d0,dmax1(mixr*rhoda(i)/(rhow + mixr*rhoda(i))     &
                                                                 ,0d0))
          t1 = anint(t1*1d20)*1d-20
          wvc(i) = t1*(nsoilp(i,2) - (soil_moist(i) + ice(i)))
          if(ntype(i) == 27) wvc(i) = t1

          wvc(i) = dmax1(0d0,dmin1(wvc(i),nsoilp(i,2) - (soil_moist(i)  &
                                                             + ice(i))))
          wvc(i) = anint(wvc(i)*1d15)*1d-15
        end if
      end do

      end subroutine soil_tmp

      subroutine soil_moisture(sm_old,wvco,iceo,sourceo,sinkro,vino,    &
                               told,thvc,dthvdh,runoff,rhs_errorm,      &
                               errorm,smerror,sert1)                  !,ii,iter

      use fasst_global

! this subroutine calcultes the soil moisture profile
! this version uses the Newton-Raphson technique

! no subroutines called

! uses the function: dense,head,soilhumid

      implicit none

      real(kind=8),intent(in):: sm_old(ntot),wvco(ntot),iceo(ntot)
      real(kind=8),intent(in):: sourceo(ntot),sinkro(ntot),vino(ntot)
      real(kind=8),intent(in):: told(ntot),thvc(ntot),dthvdh(ntot)
      real(kind=8),intent(inout):: runoff(ntot)
      real(kind=8),intent(out):: rhs_errorm,errorm,smerror,sert1

! local variables
      integer(kind=4):: i,j,ni,d0i,d1i,d2i,nb,iflag,ni1,nb1,sn2,sn3
      real(kind=8):: f1,f2,f3,rhow,rhoi,sd,sm_old1,t3,dwi,didh,dvdh
      real(kind=8):: rsinkl,rsinkh,ff,w1,c1,dense,head,bet,fts,t1
      real(kind=8):: smt,totmoisture,oldmoist,totsink,pho,smh,excess
      real(kind=8):: tsign,smv,sumvin,sms,diw
      real(kind=8):: gam(nnodes),D(nnodes)
      real(kind=8):: delphead(nnodes),A(nnodes),B(nnodes),C(nnodes)


! zero-out variables
      d0i = 0
      d1i = 0
      d2i = 0
      ni = 0
      nb = 0
      ni1 = 0
      nb1 = 0
      sn2 = 0
      sn3 = 0
      f1 = 0d0
      f2 = 0d0
      f3 = 0d0
      rsinkl = 0d0
      rsinkh = 0d0
      ff = 0d0
      w1 = 0d0
      c1 = 0d0
      bet = 0d0
      fts = 0d0
      smt = 0d0
      totmoisture = 0d0
      oldmoist = 0d0
      totsink = 0d0
      pho = 0d0
      smh = 0d0
      excess = 0d0
      t1  = 0d0
      t3 = 0d0
      sm_old1 = 0d0
      sd = 0d0
      rhow = 0d0
      rhoi = 0d0
      tsign = 0d0
      dwi = 0d0
      didh = 0d0
      dvdh = 0d0
      smv = 0d0
      sumvin = 0d0
      sms = 0d0
      sert1 = 0d0
      
! zero out arrays
      do i=1,nnodes
        delphead(i) = 0d0
        A(i) = 0d0
        B(i) = 0d0
        C(i) = 0d0
        D(i) = 0d0
        gam(i) = 0d0
        runoff(i) = 0d0
        sinkr(i) = 0d0
      end do

! initialize errors; constants
      rhs_errorm = -dabs(mflag)
      errorm = -dabs(mflag)
      smerror = -dabs(mflag)

      d1i = 1
      d2i = 2
      fts = 1d0/deltat

! determine sources and sinks
      do i=1,nnodes
        if(ntype(i) /= 27) then
          if(veg_flagl == 0.or.dabs(sigfl) <= eps) then                  !low vegetation
            rsinkl = 0d0                                                 !root uptake (unitless)
            frl(i) = 0d0
          else
            if(trmlm > eps) then
              ff = 1d0
              if(stt(i) <= Tref) ff = 0d0
              rsinkl = 1.5d-7*sigfl*ff*frl(i)                            !m/s
            else
              rsinkl = 0d0
            end if
          end if

          if(veg_flagh == 0) then                                        !high vegetation (trees)
            rsinkh = 0d0
            frh(i) = 0d0
          else
            if(trmhm > eps) then
              ff = 1d0
              if(stt(i) <= Tref) ff = 0d0
              rsinkh = 1.5d-7*sigfh*ff*frh(i)                            !m/s
            else
              rsinkh = 0d0
            end if
          end if

          sinkr(i) = (rsinkl + rsinkh)*deltat/delzs(i)                   !unitless
        else
          sinkr(i) = 0d0
        end if
      end do

! solve for parameters needed in the tridiagonal matrix
! A(i) = coeff. for delphead(i-1)
! B(i) = coeff. for delphead(i)
! C(i) = coeff. for delphead(i+1)
! D(i) = rhs
! z is positive upwards from the bottom

      do i=1,nnodes !ni+1,nb-1
        if(ntype(i) /= 27) then
          dwi = 0d0
          diw = 0d0
          didh = 0d0
          dvdh = 0d0
          if(ice(i)+iceo(i) > eps) then
            dwi = 5d-1*(dense(stt(i),0d0,d1i)/dense(stt(i),0d0,d2i)     &
                        + dense(told(i),0d0,d1i)/dense(told(i),0d0,d2i)) !unitless (water/ice)
            diw = 5d-1*(dense(stt(i),0d0,d2i)/dense(stt(i),0d0,d1i)     &
                        + dense(told(i),0d0,d2i)/dense(told(i),0d0,d1i)) !unitless (ice/water)
            didh = -dwi*dsmdh(i)
          end if

          if(wvc(i)+wvco(i) > eps)                                      &
            dvdh = dthvdh(i)*(nsoilp(i,2) - (soil_moist(i) + ice(i)))   &
                                            - thvc(i)*(1 - dwi)*dsmdh(i) !unitless (vapor/water)

          f1 = deltat/(2d0*delzs(i))                                     !s/m
          f2 = (soil_moist(i) - sm_old(i)) + (wvc(i) - wvco(i))         &
                                                + diw*(ice(i) - iceo(i)) !unitless
          f3 = f1*(source(i) + sourceo(i)) - 5d-1*(sinkr(i) + sinkro(i)) !unitless

          A(i) = -f1*khl(i)                                              !1/m
          B(i) = f1*(khl(i) + khu(i)) + dsmdh(i) + didh + dvdh           !1/m
          C(i) = -f1*khu(i)                                              !1/m
          D(i) = -f2 - f1*(vin(i) + vino(i)) + f3                        !unitless

          A(i) = anint(A(i)*1d15)*1d-15
          B(i) = anint(B(i)*1d15)*1d-15
          C(i) = anint(C(i)*1d15)*1d-15
          D(i) = anint(D(i)*1d15)*1d-15

          if(dabs(D(i)) > rhs_errorm) then
            rhs_errorm = dabs(D(i))
            sert1 = D(i)
          end if
        else if(ntype(i) == 27) then
          D(i) = 0d0
          if(sn2 == 0) sn2 = i
          sn3 = i
        end if
      end do

! solve the matrix eqn. for delta phead(i)
      if(sn2 == 0.and.sn3 == 0) then                                     !no air nodes
        ni = 1
        nb = nnodes
      else if(sn2 == 1.and.sn3 == nnodes) then                           !all air nodes
        ni = 0
        nb = 0
      else if(sn2 > 1.and.sn3 == nnodes) then                            !air on top
        ni = 1
        nb = sn2 - 1
      else if(sn2 == 1.and.sn3 < nnodes) then                            !air on bottom
        ni = sn3 + 1
        nb = nnodes
      else if(sn2 > 1.and.sn3 < nnodes) then                             !air layer
        ni = 1
        nb = sn2 - 1
        ni1 = sn3 + 1
        nb1 = nnodes
      end if

      if(ni > 0.and.nb > 0) then
        bet = B(ni)                                                      !1/m
        if(dabs(bet) <= 1d-10) then
          tsign = 1d0
          if(bet < 0d0) tsign = -1d0
          bet = tsign*1d-10
          delphead(ni) = tsign*dabs(D(ni))
        else
          delphead(ni) = D(ni)/bet                                       !m
        end if

        do j=ni+1,nb
          gam(j) = C(j-1)/bet                                            !unitless
          bet = B(j) - A(j)*gam(j)                                       !1/m
          if(dabs(bet) <= 1d-10) then  !.or.dabs(A(j)) <= eps
            tsign = 1d0
            if(bet < 0d0) tsign = -1d0
            bet = 1d-10*tsign                                            !tridiag fails, not enough moisture flow
            delphead(j) = tsign*dabs(D(j))
          else
            delphead(j) = (D(j) - A(j)*delphead(j-1))/bet                !m
          end if
        end do

        do j=nb-1,ni,-1
          delphead(j) = delphead(j) - gam(j+1)*delphead(j+1)             !m
        end do
      end if

      if(ni1 > 0.and.nb1 > 0) then
        bet = B(ni1)                                                     !1/m
        if(dabs(bet) <= 1d-10) then
          tsign = 1d0
          if(bet < 0d0) tsign = -1d0
          bet = tsign*1d-10
          delphead(ni1) = tsign*dabs(D(ni1))
        else
          delphead(ni1) = D(ni1)/bet                                     !m
        end if

        do j=ni1+1,nb1
          gam(j) = C(j-1)/bet                                            !unitless
          bet = B(j) - A(j)*gam(j)                                       !1/m
          if(dabs(bet) <= 1d-10) then  !.or.dabs(A(j)) <= eps
            tsign = 1d0
            if(bet < 0d0) tsign = -1d0
            bet = 1d-10*tsign                                            !tridiag fails, not enough moisture flow
            delphead(j) = tsign*dabs(D(j))
          else
            delphead(j) = (D(j) - A(j)*delphead(j-1))/bet                !m
          end if
        end do

        do j=nb1-1,ni1,-1
          delphead(j) = delphead(j) - gam(j+1)*delphead(j+1)             !m
        end do
      end if

! determine the maximum error, update phead(i)
! calculate the corresponding soil moisture
      do i=nnodes,1,-1
        delphead(i) = anint(delphead(i)*1d15)*1d-15

        if(ntype(i) /= 27) then
          rhow = dense(stt(i),0d0,d1i)
          rhoi = dense(stt(i),0d0,d2i)
          dwi = rhoi/rhow
          f1 = soil_moist(i) + ice(i)*dwi
          f3 = anint(delphead(i)*1d15)*1d-15

          tsign = 0d0
          if(dabs(vin(i)) > eps) then
            tsign = -dabs(vin(i))/vin(i)
          else if(dabs(vin(i)) <= eps.and.                              &
                                          dabs(delphead(i)) > eps) then
            tsign = dabs(delphead(i))/delphead(i)
          end if

          if(dabs(delphead(i)) > dabs(pheadmin(i))*2.5d-2)              &
            delphead(i) = dabs(pheadmin(i))*2.5d-2*tsign

          delphead(i) = anint(delphead(i)*1d15)*1d-15

          if(dabs(delphead(i)) > errorm) errorm = dabs(delphead(i))      !m

          sm_old1 = soil_moist(i)
          pho = phead(i)
          excess = 0d0

          if(i < ni) then
            phead(i) = phead(i)                                          !m
            smt = soil_moist(i)                                          !unitless
          else
            phead(i) = phead(i) + delphead(i) + sd*delzs(i)              !m

            if(phead(i) < 0d0.and.phead(i) > pheadmin(i)) then
              iflag = 0
              c1 = 1d0 + (nsoilp(i,10)*dabs(phead(i)*1d2))**nsoilp(i,11)
              w1 = dmax1(0d0,1d0/(c1**nsoilp(i,12)))
              smh = nsoilp(i,8) + w1*(nsoilp(i,9) - nsoilp(i,8))
              smh = dmax1(0d0,smh)
            else
              if(phead(i) >= 0d0) then
                iflag = 2
              else
                iflag = 1
              end if
              smh = nsoilp(i,24) - ice(i)*rhoi/rhow
            end if

            t3 = dsmdh(i)*delphead(i)
            sms = dmax1(0d0,soil_moist(i) + t3 + sd)                     !unitless

            smv = sm_old1 - vin(i)*(deltat/delzs(i)) + sd
            smv = dmax1(0d0,smv)

            if(iflag == 0.and.ice(i)+iceo(i) < eps) then
              smt = smh
              if(smt > nsoilp(i,24).or.dabs(f3-delphead(i)) > eps)      &
                smt = sms
            else
              if(iflag == 2.and.anint(met(iw,ip_pt)) == 2) then
                smt = smh
              else
                smt = sms
!                if(ice(i) > eps) smt = 5d-1*(sms + smv)
              end if
            end if

            if((node_type(i) == 'CO'.or.node_type(i) == 'AS')           &
                                    .or.node_type(i) == 'RO') smt = smv

            smt = anint(smt*1d10)*1d-10
            f2 = smt + ice(i)*dwi

! check for continuity in time
            t1 = dabs(vin(i)*deltat/delzs(i)) + sd
            if(t1 < eps) t1 = 1d-2

            if(dabs(f2-f1) > t1) then
              if(f2 > f1) then
                f2 = f1 + t1
              else
                f2 = f1 - t1
              end if
              smt = f2 - ice(i)*dwi
            end if
            smt = anint(smt*1d10)*1d-10

! check for over/under allowed water content
            excess = 0d0
            if(ice(i) <= eps) then
              if(smt-nsoilp(i,15) < eps) then
                smt = nsoilp(i,15)
              else if(smt-nsoilp(i,24) > eps) then
                excess = smt - nsoilp(i,24)                              !unitless
                smt = nsoilp(i,24)
              end if
            else if(ice(i) > eps) then
              if(smt < eps) then
                smt = 0d0
              else if(smt > nsoilp(i,24)) then
                excess = smt - nsoilp(i,24)
                smt = nsoilp(i,24)
              end if
            end if
            smt = anint(smt*1d10)*1d-10
            excess = anint(excess*1d10)*1d-10

            sd = 0d0
            if(excess > 0d0) then
!              if(i /= 1) then
!                sd = excess                                             !unitless
!              else if(i == 1) then
                runoff(i) = runoff(i) + excess*delzs(i)*fts              !m/s
!                end if
            end if
          end if
          sink(i) = sinkr(i)*delzs(i)*fts                                !m/s
          sink(i) = anint(sink(i)*1d10)*1d-10
          runoff(i) = anint(runoff(i)*1d10)*1d-10

          soil_moist(i) = smt
          phead(i) = head(i,soil_moist(i))

          totmoisture = totmoisture + soil_moist(i)*delzs(i)             !m
          totsink = totsink + sinkr(i)*delzs(i) + runoff(i)*deltat       !m
          oldmoist = oldmoist + sm_old1*delzs(i)                         !m
          sumvin = sumvin + vin(i)*deltat                                !m
          smerror = dabs((totmoisture - oldmoist) - (totsink + sumvin))  !m
        else
          if(dabs(delphead(i)) > errorm) errorm = dabs(delphead(i))      !m
          soil_moist(i) = 1d-2*met(iw,ip_rh)
        end if
      end do

      end subroutine soil_moisture


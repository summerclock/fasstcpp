      subroutine smooth(ntemp,ftempo)
! this subroutines smooths the temperature, soil moisture and ice variables

      use fasst_global

! no subroutines called
! uses the function: head

      implicit none

      integer(kind=4),intent(in):: ntemp
      real(kind=8),intent(in):: ftempo

! saved variables
      real(kind=8):: avestt(0:4,maxn+2),aveft(0:4),avesm(0:4,maxn)

      save:: avestt,aveft,avesm

      integer(kind=4):: kk,i,j,ii,ll,ij,tsign(maxn),iv
      real(kind=8):: temp,sumt,sumsm,mdel,deep1(maxn),tf,delsnow
      real(kind=8):: delst(maxn),delat,deep(maxn),pfac,sfac,tfac
      real(kind=8):: delsm(maxn),tempsm,maxdel !,head


! initialize variables
      ll = 0
      ii = 0
      temp = 0d0
      tf = 0d0
      delsnow = 0d0
      pfac = 0d0
      sfac = 0d0
      tfac = 0d0

      pfac = 1d0
      if(sigfl+sigfh > eps)                                             &
        pfac = dmin1(1d0,dmax1(7.5d-1,1d0-(sigfl + sigfh)))

      sfac = 1d0
!      if(hsaccum > eps) then
!        if(hsaccum <= 5d-2) then
!          sfac = 1d0  !0.75d0
!        else if(hsaccum >= 5d-1) then
!          sfac = 7.5d-1  !2.5d-1
!        else
!!          sfac = -(10d0/9d0)*hsaccum + 7.25d0/9d0
!          sfac = -(5d0/9d0)*hsaccum + 9.25d0/9d0
!        end if
!      end if

      if(dabs(pfac-1d0) <= eps) then
        tfac = sfac
      else if(dabs(sfac-1d0) <= eps) then
        tfac = pfac
      else
        tfac = 5d-1*(pfac + sfac) !dmax1(pfac,sfac)
      end if

      ll = 1
      if(timstep >= 2) ll = 1

      if(dabs(hm) <= eps) kk = 0

      if(iw == istart) then
        do i=1,ntot
          do j=0,ll
            avestt(j,i) = 0d0
            avesm(j,i) = 0d0
          end do
          avestt(0,i) = too(i)
          avesm(0,i) = smoo(i)
          if(infer_test == 1) then
            avestt(1,i) = too(i)
            avesm(1,i) = smoo(i)
          end if
        end do

        do j=0,ll
          aveft(j) = 0d0
        end do
        aveft(0) = ftempo
        if(infer_test == 1) aveft(1) = ftempo
      end if

! air
      j = 0
      delat = 0d0
      ii = iw + 1
      if(ii > ll) ii = ll
      if(iw > ll) then
        do while(j <= ii-1)
          delat = delat + dabs(met(iw-j,ip_tmp) - met(iw-(j+1),ip_tmp))
          j = j + 1
        end do
        delat = delat/real(ii)
      end if

! snow
      if(hm > eps) then
        i = nnodes + 1
        deep1(i) = 10d0*timstep !*5

        if(icaseo == 0.or.icaseo == 2) then
          kk = 1
          avestt(0,i) = stt(i)
        else
          kk = kk + 1
        end if

        if(kk > ll) then
          j = 0
          temp = avestt(j+1,i)
          do while (j < ll)
            avestt(j,i) = anint(temp*1d15)*1d-15
            j = j + 1
            if(j <= ll-1) temp = avestt(j+1,i)
          end do
          avestt(ll,i) = anint(stt(i)*1d15)*1d-15
          ii = ll
          j = ll
        else
          avestt(kk,i) = anint(stt(i)*1d15)*1d-15
          ii = kk
          j = kk
        end if
        iv = ii

        j = 0
        delst(i) = 0d0
        tsign(i) = 0
        do while(j <= ii-1)
          delst(i) = delst(i) + dabs(avestt(j+1,i) - avestt(j,i))            
          if(j == 0.and.delst(i) /= 0d0)                                &
            tsign(i) = dabs(avestt(j+1,i) - avestt(j,i))/               &
                                          (avestt(j+1,i) - avestt(j,i))
          j = j + 1
        end do

        if(ii /= 1) delst(i) = delst(i)/real(ii-1)

        if(dabs(avestt(ii,i)-avestt(ii-1,i)) > deep1(i)) then
          if(avestt(ii,i) > avestt(ii-1,i)) then
            avestt(ii,i) = avestt(ii-1,i) + deep1(i)
          else
            avestt(ii,i) = avestt(ii-1,i) - deep1(i)
          end if
          delsnow = deep1(i)
        else
          delsnow = delst(i)
        end if

        if(avestt(ii,i) > met(iw,ip_tmp)+273.15d0+2d1) then
          avestt(ii,i) = met(iw,ip_tmp) + 273.15d0 + 2d1
        else if(avestt(ii,i) < met(iw,ip_tmp)+273.15d0-2d1) then
          avestt(ii,i) = met(iw,ip_tmp) + 273.15d0 - 2d1
        end if

        sumt = 0d0
        do while(j > -1)
          sumt = sumt + avestt(j,i)
          j = j - 1
        end do

        if(iw /= 1) stt(i) = anint((sumt/real(iv + 1))*1d15)*1d-15

        tmelt(iw) = stt(i)
        if(stt(i)-273.15d0 > eps) then
          stt(i) = 273.15d0
          avestt(ii,i) = stt(i)
        end if
      end if

! ground
      do i=nnodes,1,-1
        if(node_type(i) /= 'AI') then
          if(iw > ll) then
            j = 0
            temp = avestt(j+1,i)
            tempsm = avesm(j+1,i)
            do while (j < ll)
              avestt(j,i) = temp
              avesm(j,i) = tempsm
              j = j + 1
              if(j <= ll-1) then
                temp = avestt(j+1,i)
                tempsm = avesm(j+1,i)
              end if
            end do
            avestt(ll,i) = stt(i)
            avesm(ll,i) = soil_moist(i)
            ii = ll
            j = ll
            ij = j
          else
            j = iw
            avestt(j,i) = stt(i)
            avesm(j,i) = soil_moist(i)
            ii = iw
            j = iw
            ij = j
          end if

          j = 0
          delst(i) = 0d0
          delsm(i) = 0d0
          tsign(i) = 0
          do while(j <= ii-1)
            delst(i) = delst(i) + dabs(avestt(j+1,i) - avestt(j,i))
            delsm(i) = delsm(i) + dabs(avesm(j+1,i) - avesm(j,i))           
            if(j == 0.and.delst(i) /= 0d0)                              &
              tsign(i) = dabs(avestt(j+1,i) - avestt(j,i))/             &
                                          (avestt(j+1,i) - avestt(j,i))
            j = j + 1
          end do

          if(ii /= 1) then
            delst(i) = delst(i)/real(ii-1)
            delsm(i) = delsm(i)/real(ii-1)
          end if
        end if  !if(node_type(i) /= 'AI')
      end do


      do i=nnodes,1,-1
        if(node_type(i) /= 'AI') then
          if(i == nnodes) then
            deep(i) = 1d1*timstep
            deep1(i) = deep(i)
          else
            deep(i) = 4d0 - dlog(elev - nz(i))
            deep1(i) = dmax1(1d-1,dabs(3d0 - dlog(elev - nz(i))))
          end if
          deep1(i) = deep1(i)*tfac

          j = ij
          mdel = dmax1(delat+5d0,1d1*timstep)
          if(iw == 1) mdel = 15d0
          mdel = mdel*tfac

          if(delst(i) > mdel) then
            mdel = delat
            if(iw == 1) mdel = delat + 5d0
            mdel = mdel*tfac

            if(i == nnodes) then
              avestt(ii,i) = avestt(ii-1,i) + mdel*tsign(i)
            else if(i == 1) then
              avestt(ii,i) = avestt(ii-1,i) + mdel*tsign(i)             &
                                           *dexp(2.5d-1*(nz(i) - elev))
            else
              if(delst(i+1) <= mdel.and.delst(i-1) <= mdel) then
                avestt(ii,i) = avestt(ii,i+1) - (avestt(ii,i+1)         &
                                    - avestt(ii,i-1))*(nz(i+1) - nz(i)) &
                                                   /(nz(i+1) - nz(i-1))
              else if(avestt(ii,i+1) <= mdel.and.                       &
                                            avestt(ii,i-1) > mdel) then
                avestt(ii,i) = avestt(ii,i+1) + mdel*tsign(i)           &
                                           *dexp(2.5d-1*(nz(i) - elev))
              else if(avestt(ii,i+1) > mdel.and.                        &
                                           avestt(ii,i-1) <= mdel) then
                avestt(ii,i) = avestt(ii,i-1) + mdel*tsign(i)           &
                                           *dexp(2.5d-1*(nz(i) - elev))
              else
                avestt(ii,i) = avestt(ii-1,i) + mdel*tsign(i)           &
                                           *dexp(2.5d-1*(nz(i) - elev))
              end if
            end if 
          end if

          if(dabs(avestt(ii,nnodes)-avestt(ii,i)) > deep(i)) then
            if(avestt(ii,nnodes) > avestt(ii,i)) then
              avestt(ii,i) = avestt(ii,nnodes) - deep(i)
            else
              avestt(ii,i) = avestt(ii,nnodes) + deep(i)
            end if
          end if

          if(ii > 1.and.dabs(avestt(ii,i)-avestt(ii-1,i))               &
                                                       > deep1(i)) then
            if(avestt(ii,i) > avestt(ii-1,i)) then
              avestt(ii,i) = avestt(ii-1,i) + deep1(i)
            else
              avestt(ii,i) = avestt(ii-1,i) - deep1(i)
            end if
          end if

          if(ntemp == nnodes.and.i.eq.nnodes) then
            if(avestt(ii,i) > met(iw,ip_tmp)+273.15d0+4d1) then
              avestt(ii,i) = met(iw,ip_tmp) + 273.15d0 + 4d1
            else if(avestt(ii,i) < met(iw,ip_tmp)+273.15d0-4d1) then
              avestt(ii,i) = met(iw,ip_tmp) + 273.15d0 - 4d1
            end if
          end if

          sumt = 0d0
          do while(j > -1)
            sumt = sumt + avestt(j,i)
            j = j - 1
          end do

          if(i /= nnodes) then
            stt(i) = sumt/real(ii+1)

            if(node_type(i) == 'SN'.and.stt(i)-273.15d0 > eps)          &
              stt(i) = 273.15d0
          else
            if(mstflag == 0) then
              stt(i) = sumt/real(ii+1)

              if(hm <= eps) tmelt(iw) = anint(stt(i)*1d15)*1d-15
              if(node_type(i) == 'SN'.and.stt(i)-273.15d0 > eps)        &
                stt(i) = 273.15d0
            end if  !if(mstflag == 0)
          end if  !if(i == nnodes)

          stt(i) = aint(stt(i)*1d15)*1d-15

! moisture check
! check for continuity in time
          j = ij
          if(iw /= 1) then
            if(aint(met(iw-1,ip_pt)) == 2) then                        !timestep after rain
              maxdel = 2.5d-1*nsoilp(i,24)   !2.5d-1
            else if(aint(met(iw,ip_pt)) == 2) then
              maxdel = 5d-1*nsoilp(i,24)    !5d-1
            else                                                       !other
              maxdel = 1.5d-1*nsoilp(i,24)          !3d-1
            end if
          else
            if(aint(met(iw,ip_pt)) == 2) then                          !rain
              maxdel = 5d-1*nsoilp(i,24)   !5d-1
            else                                                       !other
              maxdel = 1.5d-1*nsoilp(i,24)          !3d-1
            end if
          end if

          if(dabs(avesm(j,i)-smoo(i)) > maxdel) then
            if(avesm(j,i) > smoo(i)) then
              avesm(j,i) = smoo(i) + maxdel
            else
              avesm(j,i) = smoo(i) - maxdel
            end if
          end if

!check for vertical continuity
          if(i < nnodes.and.node_type(i) == node_type(i+1)) then
            maxdel = maxdel*dexp(nz(i) - elev)
            if(dabs(avesm(j,i+1)-avesm(j,i)) > maxdel) then
              if(avesm(j,i) > avesm(j,i+1)) then
                avesm(j,i) = avesm(j,i+1) + maxdel
              else
                avesm(j,i) = avesm(j,i+1) - maxdel
              end if
            end if
          end if

          if(ice(i) <= eps) then
            if(avesm(j,i)-nsoilp(i,15) < eps) then
              avesm(j,i) = nsoilp(i,15)
            else if(avesm(j,i)-nsoilp(i,24) > eps) then
              avesm(j,i) = nsoilp(i,24)
            end if
          else
            if(avesm(j,i)-(1d0 - ilim)*nsoilp(i,9) < eps) then
              avesm(j,i) = (1d0 - ilim)*nsoilp(i,9)
            else if(avesm(j,i)-ilim*nsoilp(i,9) > eps) then
              avesm(j,i) = ilim*nsoilp(i,9)
            end if
          end if
          avesm(j,i) = anint(avesm(j,i)*1d10)*1d-10

          sumsm = 0d0
          do while(j > -1)
            sumsm = sumsm + avesm(j,i)
            j = j - 1
          end do
!          soil_moist(i) = sumsm/dreal(ii+1)
!          soil_moist(i) = anint(soil_moist(i)*1d10)*1d-10
!          phead(i) =  head(i,soil_moist(i))
        end if  !(node_type(i) /= 'AI')
      end do

! veg
      if(icase == 2.or.icase == 3) then
        deep1(ntemp) = 5d0*timstep

        do i=nnodes+1,ntemp
          if(iw <= ll) then
            aveft(iw) = ftemp
            ii = iw
            j = iw
            ij = j
          else
            j = 0
            temp = aveft(j+1)
            do while (j < ll)
              aveft(j) = temp
              j = j + 1
              if(j <= ll-1) then
                temp = aveft(j+1)
              end if
            end do
            aveft(ll) = anint(ftemp*1d15)*1d-15
            ii = ll
            j = ll
            ij = j
          end if

        j = 0
        delst(i) = 0d0
        tsign(i) = 0
        do while(j <= ii-1)
          delst(i) = delst(i) + dabs(aveft(j+1) - aveft(j))            
          if(j == 0.and.delst(i) /= 0d0)                                &
            tsign(i) = dabs(aveft(j+1) - aveft(j))/                     &
                                                (aveft(j+1) - aveft(j))
          j = j + 1
        end do

        if(ii /= 1) then
          delst(i) = delst(i)/real(ii-1)
        end if

        mdel = dmax1(delat+5d0,10d0*timstep)
        if(iw == 1) mdel = 10d0
        if(delst(i) > mdel) then
          mdel = deltat
          if(iw == 1) mdel = deltat + 5d0
          aveft(ii) = aveft(ii-1) + mdel*tsign(i)
        end if

          if(dabs(aveft(ii)-aveft(ii-1)) > deep1(ntemp)) then
            if(aveft(ii) > aveft(ii-1)) then
              aveft(ii) = aveft(ii-1) + deep1(ntemp)
            else
              aveft(ii) = aveft(ii-1) - deep1(ntemp)
            end if
          end if

          tf = 4d1
!          if(met(iw,ip_prec)+met(iw,ip_prec2) > eps) tf = 2d1
          if(aveft(ii) > met(iw,ip_tmp)+273.15d0+tf) then
            aveft(ii) = met(iw,ip_tmp) + 273.15d0 + tf
          else if(aveft(ii) < met(iw,ip_tmp)+273.15d0-tf) then
            aveft(ii) = met(iw,ip_tmp) + 273.15d0 - tf
          end if

          sumsm = 0d0
          do while(ij > -1)
            sumsm = sumsm + aveft(ij)
            ij = ij - 1
          end do

          temp = ftemp
          ftemp = sumsm/real(ii+1)

          if(icase == 2) then
            stt(ntemp) = ftemp
          else if(icase == 3) then
            if(hfol_tot-hm <= 1d-3) ftemp = stt(ntemp-1)
            stt(ntemp) = ftemp
          end if

          ftemp = anint(ftemp*1d15)*1d-15
          stt(ntemp) = anint(stt(ntemp)*1d15)*1d-15
        end do
      else if(icase == 4) then
        ftemp = stt(nnodes) + (hfol_tot/hm)*(stt(nnodes+1)              &
                                                         - stt(nnodes))
        aveft(iv) = ftemp
        ftemp = anint(ftemp*1d15)*1d-15
      end if  !if(icase == 2.or.icase == 3)

      end subroutine smooth
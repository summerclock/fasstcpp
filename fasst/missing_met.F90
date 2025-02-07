      subroutine missing_met(wmode,ipos,ie,met1)

      use fasst_global
      use module_radiation

! calls the following subroutines:
!     sol_zen
!     Solflx
!     emisatm
!     dnirflx

! uses the function: met_date

!     case(1): infer_test = 1, continuation of previous run
!     case(2): fill in missing parameters
!     case(3): fill in missing time steps


      implicit none

      integer(kind=4),intent(in):: ipos,ie
      integer(kind=4),intent(in):: wmode
      real(kind=8),intent(in):: met1(1,maxcol)

! local variables
      integer(kind=4):: i,j,k,c3,msteps,tstps,icld(3),ib,ib0
      real(kind=8):: met2(maxlines,maxcol),prcp,met_date,rtstps,modyo
      real(kind=8):: ematm,cover(3),hgt(3),mcd,hcd,t1,t2o,mody,f1,daylim


      tstps = 0
      rtstps = 0d0
      mody = 0d0
      modyo = 0d0
      daylim = 0d0
      prcp = 0d0
      ematm = 0d0
      mcd = 0d0
      hcd = 0d0
      t1 = 0d0
      t2o = 0d0
      f1 = 0d0

      rtstps = 24d0/timstep
      tstps = int(rtstps)

      select case(wmode)
! ******************************************************************************
        case(1)
! initialize variables
          c3 = 0

          do i=1,maxlines
            do j=1,maxcol
              met2(i,j) = mflag
            end do
          end do

! determine initial times of new and old data
! old
          t1 = met_date(met1(1,ip_year),met1(1,ip_doy),                 &
                        met1(1,ip_hr),met1(1,ip_min))

! new
          t2o = met_date(met(istart,ip_year),met(istart,ip_doy),        &
                         met(istart,ip_hr),met(istart,ip_min))

          if(t2o <= t1) then                                             !files overlap
            do j=istart,iend
              c3 = c3 + 1
              do k=1,maxcol
                met2(c3,k) = met(j,k)
              end do
            end do
          else
            c3 = c3 + 1
            do k=1,maxcol
              met2(c3,k) = met1(1,k)
            end do

            mody = met(istart,ip_year) -                                 &
                                    aint(met(istart,ip_year)*2.5d-1)*4d0
            modyo = met1(1,ip_year) - aint(met1(1,ip_year)*2.5d-1)*4d0
            daylim = 1.042d0  !1.003d0
!            if(mody /= modyo) daylim = 1.07d0
            if(met(istart,ip_year) /= met(1,ip_year)) daylim = 25.07d0

            if(dabs(mody) <= eps) then
              msteps = int((t2o - t1)*367d0*tstps)
              if((t2o - t1)*367d0*tstps < daylim) msteps = 0
            else
              msteps = int((t2o - t1)*366d0*tstps)
              if((t2o - t1)*366d0*tstps < daylim) msteps = 0
            end if

            if(msteps > int(mgap*tstps)) then
              write(*,*) 'Met files too far apart, stopping'
              stop
            else if(msteps /= 0.and.msteps <= aint(2d0*tstps)) then
              do i=1,msteps
                c3 = c3 + 1

                met2(c3,1) = met1(1,1)
                met2(c3,2) = met1(1,2)

                mody = timstep*10d0 - aint(timstep)*10d0
                if(dabs(mody) <= eps) then
                  met2(c3,3) = met1(1,3) + i*timstep
                  met2(c3,4) = met1(1,4)
                else
                  met2(c3,3) = met1(1,3) + i*aint(timstep)
                  met2(c3,4) = met1(1,4) + i*(mody*1d-1)
                end if

                if(met2(c3,4) >= 6d1) then
                  met2(c3,3) = met1(1,3) + aint(met2(c3,4)/6d1)
                  met2(c3,4) = met2(c3,4) - aint(met2(c3,4)/6d1)*6d1
                end if

                if(met2(c3,3) >= 24d0) then
                  met2(c3,2) = met1(1,2) + aint(met2(c3,3)/24d0)
                  met2(c3,3) = met2(c3,3) - aint(met2(c3,3)/24d0)*24d0
                end if

                mody = met1(1,1) - aint(met1(1,1)*2.5d-1)*4d0
                if(met2(c3,2) > 366d0.and.dabs(mody) <= eps) then
                  met2(c3,2) = 1d0
                  met2(c3,1) = met2(c3,1) + 1d0
                else if(met2(c3,2) > 365d0.and.dabs(mody) > eps) then
                  met2(c3,2) = 1d0
                  met2(c3,1) = met2(c3,1) + 1d0
                end if

                call sol_zen(met2(c3,ip_year),met2(c3,ip_doy),          &
                             met2(c3,ip_hr),met2(c3,ip_min),            &
                             met2(c3,ip_zen),met2(c3,ip_az))

                f1 = 1d0/(msteps + 1d0)
                do j=5,ncols
                  if(j /= 29.or.j /= 30) then                            !ip_zen, ip_az
                    if(j == 11.or.j == 13) then
                      met2(c3,j) = met1(1,j)
                    else if(j == 16.or.j == 19) then
                      met2(c3,j) = met1(1,j)
                    else if(j == 22.or.j == 35) then
                      met2(c3,j) = met1(1,j)
                    else
                      met2(c3,j) = met1(1,j) + i*(met(1,j)              &
                                                        - met1(1,j))*f1
                    end if
                  end if
                end do

                met2(c3,ip_sd) = mflag

                if(dabs(met2(c3,ip_prec)) <= eps) then
                  met2(c3,ip_pt) = 1
                else if(met2(c3,ip_prec) > eps.and.                     &
                                            met2(c3,ip_tmp) > 0d0) then
                  met2(c3,ip_pt) = 2
                else if(met2(c3,ip_prec) > eps.and.                     &
                                            met2(c3,ip_tmp) <= 0d0) then
                  met2(c3,ip_pt) = 3
                end if

                if(dabs(met2(c3,ip_prec2)) <= eps) then
                  met2(c3,ip_pt2) = 1
                else
                  met2(c3,ip_pt2) = 3
                end if

                if(met2(c3,ip_zen) >= 90d0) then
                  met2(c3,ip_tsol) = 0d0
                  met2(c3,ip_dir) = 0d0
                  met2(c3,ip_dif) = 0d0
                  met2(c3,ip_upsol) = 0d0
                else
                  met2(c3,ip_tsol) = mflag
                  met2(c3,ip_dir) = mflag
                  met2(c3,ip_dif) = mflag
                  met2(c3,ip_upsol) = mflag

                  cover(1) = met2(c3,ip_lcd)
                  cover(2) = met2(c3,ip_mcd)
                  cover(3) = met2(c3,ip_hcd)
                  hgt(1) = met2(c3,ip_lhgt)
                  hgt(2) = met2(c3,ip_mhgt)
                  hgt(3) = met2(c3,ip_hhgt)
                  icld(1) = int(met2(c3,ip_lct))
                  icld(2) = int(met2(c3,ip_mct))
                  icld(3) = int(met2(c3,ip_hct))

                  prcp = met2(c3,ip_prec) + met2(c3,ip_prec2)

                  call Solflx(icld,met2(c3,ip_zen),met2(c3,ip_doy),prcp,&
                              met2(c3,ip_tsol),cover,hgt,               &
                              met2(c3,ip_dir),met2(c3,ip_dif))

                  if(dabs(met2(c3,ip_tsol)) <= eps)                     &
                    met2(c3,ip_upsol) = 0d0
                  met2(c3,ip_lcd) = cover(1)
                  met2(c3,ip_mcd) = cover(2)
                  met2(c3,ip_hcd) = cover(3)
                  met2(c3,ip_lhgt) = hgt(1)
                  met2(c3,ip_mhgt) = hgt(2)
                  met2(c3,ip_hhgt) = hgt(3)
                  met2(c3,ip_lct) = real(icld(1))
                  met2(c3,ip_mct) = real(icld(2))
                  met2(c3,ip_hct) = real(icld(3))     
                end if

                mcd = met2(c3,ip_mcd)
                hcd = met2(c3,ip_hcd)
                call emisatm(met2(c3,ip_tmp),met2(c3,ip_rh),ematm)
                call dnirflx(ematm,met2(c3,ip_tmp),met2(c3,ip_doy),     &
                             met2(c3,ip_lhgt),met2(c3,ip_mhgt),         &
                             met2(c3,ip_hhgt),met2(c3,ip_lcd),          &
                             mcd,hcd,met2(c3,ip_ir))
              end do   !i=1,msteps

              do i=istart,iend
                c3 = c3 + 1
                do j=1,maxcol
                  met2(c3,j) = met(i,j)
                end do
              end do
            else if(msteps == 0) then
              do i=istart,iend
                c3 = c3 + 1
                do j=1,maxcol
                  met2(c3,j) = met(i,j)
                end do
              end do
!              c3 = c3 - 1
            end if   !if(msteps > mgap)
          end if   !if(t2o > t1)

          do i=1,c3
            do j=1,maxcol
              met(i,j) = met2(i,j)
            end do
          end do

          iend = c3

! ******************************************************************************
        case(2)

          do i=ipos,ie
! time
            if(aint(dabs(met(i,3)-mflag)*1d5)*1d-5 <= eps) then          !hour
              if(i /= istart) then
                met(i,3) = met(ipos,3) + (i - ipos)*timstep
              else
                if(dabs(met(ipos+1,3)) <= eps) met(ipos+1,3) = 24d0
                met(i,3) = met(ipos+1,3) - (i - ipos + 1)*timstep
              end if
            end if

            if(aint(dabs(met(i,4)-mflag)*1d5)*1d-5 <= eps) then          !minute
              if(i /= istart) then
                met(i,4) = met(ipos,4)                                  & 
                                 + (i - ipos)*(timstep - aint(timstep))
              else
                met(i,4) = met(ipos+1,4)                                & 
                              - (i - ipos + 1)*(timstep - aint(timstep))
              end if
            end if

            if(met(i,4) >= 6d1) then
              met(i,3) = met(ipos,3) + aint(met(i,4)/6d1)
              met(i,4) = met(i,4) - aint(met(i,4)/6d1)*6d1
            end if

            if(met(i,3) >= 24d0) then
              met(i,2) = met(ipos,2) + aint(met(i,3)/24d0)
              met(i,3) = met(i,3) - aint(met(i,3)/24d0)*24d0
            end if

            mody = met(i,1) - aint(met(i,1)*2.5d-1)*4d0
            if(aint(met(i,2)) > aint(366d0).and.dabs(mody) <= eps) then
              met(i,2) = 1d0
              met(i,1) = met(i,1) + 1d0
            else if(aint(met(i,2)) > aint(365d0).and.dabs(mody) > eps) then
              met(i,2) = 1d0
              met(i,1) = met(i,1) + 1d0
            end if

! wind, pressure, temp, precipitation, rh
            do j=5,12
              if(aint(dabs(met(i,j)-mflag)*1d5)*1d-5 <= eps) then
                if(i == istart) then
                  met(i,j) = met(ipos+1,j)
                  if(aint(dabs(met(i,j)-mflag)*1d5)*1d-5 <= eps) then
                    if(j == 5) met(i,j) = 1d3
                    if(j == 7) met(i,j) = 6d1
                    if(j == 8) met(i,j) = 1.2d0
                    if(j == 9) met(i,j) = 0d0
                    if(j == 10.or.j == 12) met(i,j) = 0d0
                    if(j == 11) then                                     !precip type
                      if(dabs(met(i,ip_prec)) <= eps) then
                        met(i,ip_pt) = 1
                      else if(met(i,ip_prec) > eps.and.                 &
                                            met(i,ip_tmp) > 0d0) then
                        met(i,ip_pt) = 2
                      else if(met(i,ip_prec) > eps.and.                 &
                                            met(i,ip_tmp) <= 0d0) then
                        met(i,ip_pt) = 3
                      end if
                    end if
                    if(j == 13) then                                     !precip2 type
                      if(met(i,12) > eps) then
                        met(i,j) = 3d0
                      else
                        met(i,j) = 1d0
                      end if
                    end if
                  end if
                else if(i >= ie) then
                  met(i,j) = met(ie-1,j)
                else
                  if(aint(dabs(met(i+1,j)-mflag)*1d5)*1d-5 > eps) then
                    met(i,j) = 5d-1*(met(i-1,j) + met(i+1,j))
                  else
                    met(i,j) = met(i-1,j)
                  end if
                end if
              end if
            end do

! clouds
            do j=13,22
              if(aint(dabs(met(i,j)-mflag)*1d5)*1d-5 <= eps) then
                if(i == istart) then
                  met(i,j) = met(ipos+1,j)
                  if(aint(dabs(met(i,j)-mflag)*1d5)*1d-5 <= eps) then
                    if(met(i,10)+met(i,12) > eps) then
                      met(i,14) = 1d0
                      met(i,16) = 6d0
                      met(i,17) = 1d0
                      met(i,19) = 3d0
                      met(i,20) = 1d0
                      met(i,22) = 5d0
                    else
                      met(i,14) = 5d-1
                      met(i,16) = 6d0
                      met(i,17) = 0d0
                      met(i,19) = 0d0
                      met(i,20) = 0d0
                      met(i,22) = 0d0
                    end if
                  end if
                else if(i >= ie) then
                  met(i,j) = met(ie-1,j)
                else
                  if(aint(dabs(met(i+1,j)-mflag)*1d5)*1d-5 > eps) then
                    met(i,j) = 5d-1*(met(i-1,j) + met(i+1,j))
                  else
                    met(i,j) = met(i-1,j)
                  end if
                end if
              end if
            end do

            do j=29,30
              if(aint(dabs(met(i,j)-mflag)*1d5)*1d-5 <= eps) then
                call sol_zen(met(i,ip_year),met(i,ip_doy),met(i,ip_hr), &
                             met(i,ip_min),met(i,ip_zen),met(i,ip_az))
              end if
            end do

! solar
            do j=23,26
              if(aint(dabs(met(i,j)-mflag)*1d5)*1d-5 <= eps) then
                if(met(i,ip_zen) >= 90d0) then
                  met(i,ip_tsol) = 0d0
                  met(i,ip_dir) = 0d0
                  met(i,ip_dif) = 0d0
                  met(i,ip_upsol) = 0d0
                else
                  met(i,ip_tsol) = mflag
                  met(i,ip_dir) = mflag
                  met(i,ip_dif) = mflag
                  met(i,ip_upsol) = mflag

                  cover(1) = met(i,ip_lcd)
                  cover(2) = met(i,ip_mcd)
                  cover(3) = met(i,ip_hcd)
                  hgt(1) = met(i,ip_lhgt)
                  hgt(2) = met(i,ip_mhgt)
                  hgt(3) = met(i,ip_hhgt)
                  icld(1) = int(met(i,ip_lct))
                  icld(2) = int(met(i,ip_mct))
                  icld(3) = int(met(i,ip_hct))

                  prcp = met(i,ip_prec) + met(i,ip_prec2)

                  call Solflx(icld,met(i,ip_zen),met(i,ip_doy),prcp,    &
                              met(i,ip_tsol),cover,hgt,met(i,ip_dir),   &
                              met(i,ip_dif))

                  met(i,ip_lcd) = cover(1)
                  met(i,ip_mcd) = cover(2)
                  met(i,ip_hcd) = cover(3)
                  met(i,ip_lhgt) = hgt(1)
                  met(i,ip_mhgt) = hgt(2)
                  met(i,ip_hhgt) = hgt(3)
                  met(i,ip_lct) = real(icld(1))
                  met(i,ip_mct) = real(icld(2))
                  met(i,ip_hct) = real(icld(3))     
                end if
              end if
            end do

! ir
            do j=27,28
              if(aint(dabs(met(i,j)-mflag)*1d5)*1d-5 <= eps) then
                if(j == 27) then
                  mcd = met(i,ip_mcd)
                  hcd = met(i,ip_hcd)
                  call emisatm(met(i,ip_tmp),met(i,ip_rh),ematm)
                  call dnirflx(ematm,met(i,ip_tmp),met(i,ip_doy),       &
                               met(i,ip_lhgt),met(i,ip_mhgt),           &
                               met(i,ip_hhgt),met(i,ip_lcd),            &
                               mcd,hcd,met(i,ip_ir))
                end if
              end if
            end do

          end do   !i=

 ! *****************************************************************************
        case(3)

! initialize variables
          c3 = 0

          do i=1,maxlines
            do j=1,maxcol
              met2(i,j) = mflag
            end do
          end do

! determine gap in data
! beginning
          t1 = met_date(met1(1,ip_year),met1(1,ip_doy),met1(1,ip_hr),   &
                        met1(1,ip_min))

! end
          t2o = met_date(met(ie,ip_year),met(ie,ip_doy),                &
                         met(ie,ip_hr),met(ie,ip_min))

          mody = met(ie,ip_year) - aint(met(ie,ip_year)*2.5d-1)*4d0
          modyo = met1(1,ip_year) - aint(met1(1,ip_year)*2.5d-1)*4d0
          daylim = 1.042d0  !1.003d0
!          if(mody /= modyo) daylim = 1.07d0
          if(met1(1,ip_year) /= met(ie,ip_year)) daylim = 25.07d0

          if(dabs(mody) <= eps) then
            msteps = int((t2o - t1)*367d0*tstps)
            if((t2o - t1)*367d0*tstps < daylim) msteps = 0
          else
            msteps = int((t2o - t1)*366d0*tstps)
            if((t2o - t1)*366d0*tstps < daylim) msteps = 0
          end if
          
          if(msteps > int(mgap*tstps)) then                              ! >= 2 (user defined) days missing

            write(*,'(''Missing met data over 2 days apart. Filled in'',&
      &'' parameters except solar and IR are dubious.  Persistence '',  &
      &''assumed.'',/,'' Year '',i5,'', DOY '',i4,'', Hour '',i3,       &
      &'', Minute '',i3,/)')int(met(ie,ip_year)),int(met(ie,ip_doy)),   &
      int(met(ie,ip_hr)),int(met(ie,ip_min))
           write(10,'(''Missing met data over 2 days apart. Filled in'',&
      &'' parameters except solar and IR are dubious.  Persistence '',  &
      &''assumed.'',/,'' Year '',i5,'', DOY '',i4,'', Hour '',i3,       &
      &'', Minute '',i3,/)')int(met(ie,ip_year)),int(met(ie,ip_doy)),   &
      int(met(ie,ip_hr)),int(met(ie,ip_min))
          end if

          ib = ipos
          if(ib /= 1) ib = ib - 1
          do while (dabs((met(ib,ip_hr)+met(ib,ip_min)/6d1)             &
                       - (met(ipos,ip_hr)+met(ipos,ip_min)/6d1)) > eps)
            ib = ib - 1
            if(ib == 0) exit
          end do
          if(ib /= 0) ib = ib + 1
          ib0 = ib

          do i=1,ipos
            c3 = c3 + 1
            do j=1,maxcol
              met2(c3,j) = met(i,j)
            end do
          end do

          do i=1,msteps
            c3 = c3 + 1

            met2(c3,1) = met1(1,1)
            met2(c3,2) = met1(1,2)

            mody = timstep*1d1 - aint(timstep)*1d1
            if(dabs(mody) <= eps) then
              met2(c3,3) = met1(1,3) + i*timstep
              met2(c3,4) = met1(1,4)
            else
              met2(c3,3) = met1(1,3) + i*aint(timstep)
              met2(c3,4) = met1(1,4) + i*(mody/1d1)
            end if

            if(met2(c3,4) >= 6d1) then
              met2(c3,3) = met1(1,3) + aint(met2(c3,4)/6d1)
              met2(c3,4) = met2(c3,4) - aint(met2(c3,4)/6d1)*6d1
            end if

            if(met2(c3,3) >= 24d0) then
              met2(c3,2) = met1(1,2) + aint(met2(c3,3)/24d0)
              met2(c3,3) = met2(c3,3) - aint(met2(c3,3)/24d0)*24d0
            end if

            mody = met1(1,1) - aint(met1(1,1)*2.5d-1)*4d0
            if(met2(c3,2) > 366d0.and.dabs(mody) <= eps) then
              met2(c3,2) = 1d0
              met2(c3,1) = met2(c3,1) + 1d0
            else if(met2(c3,2) > 365d0.and.dabs(mody) > eps) then
              met2(c3,2) = 1d0
              met2(c3,1) = met2(c3,1) + 1d0
            end if

            call sol_zen(met2(c3,ip_year),met2(c3,ip_doy),              &
                         met2(c3,ip_hr),met2(c3,ip_min),                &
                         met2(c3,ip_zen),met2(c3,ip_az))

            f1 = 1d0/(msteps + 1d0)
            do j=5,ncols
              if(j /= 29.or.j /= 30) then                                !ip_zen, ip_az
                if(ib == 0) then
                  if(j == 11.or.j == 13) then
                    met2(c3,j) = met1(1,j)
                  else if(j == 16.or.j == 19) then
                    met2(c3,j) = met1(1,j)
                  else if(j == 22.or.j == 35) then
                    met2(c3,j) = met1(1,j)
                  else
                    met2(c3,j) = met1(1,j) + i*(met(ie,j)               &
                                                        - met1(1,j))*f1
                  end if
                else
                  met2(c3,j) = met(ib,j)
                end if
              end if
            end do

            met2(c3,ip_sd) = mflag

            if(dabs(met2(c3,ip_prec)) <= eps) then
              met2(c3,ip_pt) = 1
            else if(met2(c3,ip_prec) > eps.and.                         &
                                           met2(c3,ip_tmp) > 0d0) then
              met2(c3,ip_pt) = 2
            else if(met2(c3,ip_prec) > eps.and.                         &
                                           met2(c3,ip_tmp) <= 0d0) then
              met2(c3,ip_pt) = 3
            end if

            if(dabs(met2(c3,ip_prec2)) <= eps) then
              met2(c3,ip_pt2) = 1
            else
              met2(c3,ip_pt2) = 3
            end if

            if(met2(c3,ip_zen) >= 90d0) then
              met2(c3,ip_tsol) = 0d0
              met2(c3,ip_dir) = 0d0
              met2(c3,ip_dif) = 0d0
              met2(c3,ip_upsol) = 0d0
            else
              met2(c3,ip_tsol) = mflag
              met2(c3,ip_dir) = mflag
              met2(c3,ip_dif) = mflag
              met2(c3,ip_upsol) = mflag

              cover(1) = met2(c3,ip_lcd)
              cover(2) = met2(c3,ip_mcd)
              cover(3) = met2(c3,ip_hcd)
              hgt(1) = met2(c3,ip_lhgt)
              hgt(2) = met2(c3,ip_mhgt)
              hgt(3) = met2(c3,ip_hhgt)
              icld(1) = int(met2(c3,ip_lct))
              icld(2) = int(met2(c3,ip_mct))
              icld(3) = int(met2(c3,ip_hct))

              prcp = met2(c3,ip_prec) + met2(c3,ip_prec2)

              call Solflx(icld,met2(c3,ip_zen),met2(c3,ip_doy),prcp,    &
                          met2(c3,ip_tsol),cover,hgt,                   &
                          met2(c3,ip_dir),met2(c3,ip_dif))

              if(dabs(met2(c3,ip_tsol)) <= eps) met2(c3,ip_upsol) = 0d0
              met2(c3,ip_lcd) = cover(1)
              met2(c3,ip_mcd) = cover(2)
              met2(c3,ip_hcd) = cover(3)
              met2(c3,ip_lhgt) = hgt(1)
              met2(c3,ip_mhgt) = hgt(2)
              met2(c3,ip_hhgt) = hgt(3)
              met2(c3,ip_lct) = real(icld(1))
              met2(c3,ip_mct) = real(icld(2))
              met2(c3,ip_hct) = real(icld(3))     
            end if

            mcd = met2(c3,ip_mcd)
            hcd = met2(c3,ip_hcd)
            call emisatm(met2(c3,ip_tmp),met2(c3,ip_rh),ematm)

            call dnirflx(ematm,met2(c3,ip_tmp),met2(c3,ip_doy),         &
                         met2(c3,ip_lhgt),met2(c3,ip_mhgt),             &
                         met2(c3,ip_hhgt),met2(c3,ip_lcd),              &
                         mcd,hcd,met2(c3,ip_ir))

            if(ib /= 0) then
              ib = ib + 1
              if(ib > ipos) ib = ib0
            end if
          end do   !i=1,msteps

          do i=ie,iend
            c3 = c3 + 1
            do j=1,maxcol
              met2(c3,j) = met(i,j)
            end do
          end do

          do i=1,c3
            do j=1,maxcol
              met(i,j) = met2(i,j)
            end do
          end do

          iend = c3

      end select

      end subroutine missing_met

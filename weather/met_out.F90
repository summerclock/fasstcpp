      subroutine met_out(iend,met)

      use met_global
      use module_radiation

! calls the following subroutines:
!     sol_zen
!     Solflx
!     emisatm
!     dnirflx
!     met_date (appended to this subroutine)

      implicit none

      integer(kind=4),intent(in):: iend
      real(kind=8),intent(in):: met(nlines,maxcol)

! local variables
      integer(kind=4),parameter:: ndays = 365
      integer(kind=4):: i,j,k,c3,c2,keep,tstps,wstart,ncols1,nlines1
      integer(kind=4):: ntlines,msteps,vitd_index1,nlines2,icld(3)
      integer(kind=4):: pcols,mcount,err,ib,ib0,diy,io
      real(kind=4):: timstep1,timeoffset1
      real(kind=8):: ematm,cover(3),hgt(3),mcd,hcd,t1,t1o,t2,t2o,prcp
      real(kind=8):: mymod,rtstps
      real(kind=8),allocatable:: met1(:,:),met2(:,:)
      character(len=300) header

      save:: nlines1

! initialize variables
      pcols = 0
      pcols = 34
      c3 = 0
      c2 = 0
      keep = 0
      tstps = 0
      wstart = 0
      ncols1 = 0
      nlines1 = 0
      ntlines = 0
      msteps = 0
      vitd_index1 = 0
      nlines2 = 0
      mcount = 0
      err = 0
      ib = 0
      ib0 = 0
      diy = 0
      timstep1 = 0.0
      timeoffset1 = 0.0
      ematm = 0d0
      mcd = 0d0
      hcd = 0d0
      t1 = 0d0
      t1o = 0d0
      t2 = 0d0
      t2o = 0d0
      prcp = 0d0
      mymod = 0d0
      rtstps = 0d0

      do i=1,3
        icld(i) = 0
        cover(i) = 0d0
        hgt(i) = 0d0
      end do

! determine the time step
      rtstps = 24d0/timstep
      tstps = int(rtstps)
      mcount = ndatpos

      io = 0
      if(infer_test == 1.and.io /= -1) then
        read(3,1000,iostat=io) nlines1,ncols1,lat,mlong,elev,           &
                               vitd_index1,mcount,timeoffset1,timstep1, &
                               mflag,iheight
        read(3,'(a)',iostat=io) header
        read(3,'(a)',iostat=io) header

        if(nlines1 > 0) then
          allocate(met1(nlines1,ncols1),stat=err)                        !size met1 array

          do i=1,nlines1
            do j=1,ncols1
              met1(i,j) = mflag
            end do
          end do

          if(vitd_index == vitd_index1) then                             !find the correct position
            do i=1,nlines1
              read(3,*,iostat=io) (met1(i,j),j=1,ncols1)
            end do

! last time step of old file
            call met_date(diy,met1(nlines1,ip_year),                    &
                          met1(nlines1,ip_doy),met1(nlines1,ip_hr),     &
                          met1(nlines1,ip_min),t1o)

! new
            call met_date(diy,met(1,ip_year),met(1,ip_doy),met(1,ip_hr),&
                          met(1,ip_min),t2o)

            msteps = int((t2o - t1o)*diy*tstps)
            if(msteps < 0) msteps = 0

            ntlines = nlines1 + nlines + msteps
            pcols = max(ncols1,maxcol)
            allocate(met2(ntlines,pcols),stat=err)                       !size met2 array

            do i=1,ntlines
              do j=1,pcols
                met2(i,j) = mflag
              end do
            end do

! first timestep of old data
            call met_date(diy,met1(1,ip_year),met1(1,ip_doy),           &
                          met1(1,ip_hr),met1(1,ip_min),t1)

! use old data first
            c3 = 1
            i = 1
            do while(t1 < t2o.and.i <= nlines1)
              do k=1,ncols1
                met2(c3,k) = met1(i,k)
              end do
              c3 = c3 + 1

              i = i + 1
              if(i > nlines1) exit

              call met_date(diy,met1(i,ip_year),met1(i,ip_doy),         &
                            met1(i,ip_hr),met1(i,ip_min),t1)
            end do

! gap between data files
            c2 = i - 1
            i = 1
            do while(msteps /= 0.and.i <= msteps)
              if(msteps > (int(2.0*tstps))) then                         ! >= 2 days missing
                write(*,'(''Missing met data over 2 days apart. '',     &
     &''Filled in parameters except solar and IR are dubious. '',       &
     &''Persistance assumed.'',/,'' Time 1: Year '',i5,'', DOY '',i4,   &
     &'', Hour '',i3,'', Minute '',i3,/,'' Time 2: Year '',i5,'',       &
     &DOY '',i4,'', Hour '',i3,'', Minute '',i3,/)')                    &
                  int(met1(nlines1,ip_year)),int(met1(nlines1,ip_doy)), &
                  int(met1(nlines1,ip_hr)),int(met1(nlines1,ip_min)),   &
                  int(met(1,ip_year)),int(met(1,ip_doy)),               &
                  int(met(1,ip_hr)),int(met(1,ip_min))

                write(10,'(''Missing met data over 2 days apart. '',    &
     &''Filled in parameters except solar and IR are dubious. '',       &
     &''Persistance assumed.'',/,'' Time 1: Year '',i5,'', DOY '',i4,   &
     &'', Hour '',i3,'', Minute '',i3,/,'' Time 2: Year '',i5,          &
     &'', DOY '',i4,'', Hour '',i3,'', Minute '',i3,/)')                &
                  int(met1(nlines1,ip_year)),int(met1(nlines1,ip_doy)), &
                  int(met1(nlines1,ip_hr)),int(met1(nlines1,ip_min)),   &
                  int(met(1,ip_year)),int(met(1,ip_doy)),               &
                  int(met(1,ip_hr)),int(met(1,ip_min))
              end if

              ib = c2 - 1
              do while (dabs((met(ib,ip_hr)+met(ib,ip_min)/6d1)-        &
                             (met(c2,ip_hr)+met(c2,ip_min)/6d1)) > eps)
                ib = ib - 1
                if(ib == 0) exit
              end do
              if(ib /= 0) ib = ib + 1
              ib0 = ib

              do i=1,msteps
                met2(c3,1) = met1(c2,1)
                met2(c3,2) = met1(c2,2)

                mymod = timstep*1d1 - aint(timstep)*1d1
                if(dabs(mymod) <= eps) then
                  met2(c3,3) = met1(c2,3) + i*timstep
                  met2(c3,4) = met1(c2,4)
                else
                  met2(c3,3) = met1(c2,3) + i*int(timstep)
                  met2(c3,4) = met1(c2,4) + i*(mymod/1d1)
                end if

                if(met2(c3,4) >= 6d1) then
                  met2(c3,3) = met1(c2,3) + aint(met2(c3,4)/6d1)
                  met2(c3,4) = met2(c3,4) - aint(met2(c3,4)/6d1)*6d1
                end if

                if(met2(c3,3) >= 24d0) then
                  met2(c3,2) = met1(c2,2) + aint(met2(c3,3)/24d0)
                  met2(c3,3) = met2(c3,3) - aint(met2(c3,3)/24d0)*24d0
                end if

                mymod = met1(c2,1) - aint(met1(c2,1)/4d0)*4d0
                if(met2(c3,2) > 366d0.and.dabs(mymod) <= eps) then
                  met2(c3,2) = 1d0
                  met2(c3,1) = met2(c3,1) + 1d0
                else if(met2(c3,2) > 365d0.and.dabs(mymod) > eps) then
                  met2(c3,2) = 1d0
                  met2(c3,1) = met2(c3,1) + 1d0
                end if

                call sol_zen(met2(c3,ip_year),met2(c3,ip_doy),          &
                             met2(c3,ip_hr),met2(c3,ip_min),            &
                             met2(c3,ip_zen),met2(c3,ip_az))

                do j=5,ncols
                  if(j /= 29.or.j /= 30) then                            !ip_zen, ip_az
                    if(ib == 0) then
                      if(j == 11.or.j == 13.or.j == 16.or.j == 19.or.   &
                         j == 22.or.j == 35) then
                        met2(c3,j) = met1(c2,j)
                      else
                        met2(c3,j) = met1(c2,j) + i*(met(1,j)           &
                                     - met1(c2,j))/(msteps + 1d0)
                      end if
                    else
                      met2(c3,j) = met(ib,j)
                    end if
                  end if
                end do

                met2(c3,ip_sd) = mflag

                if(dabs(met2(c3,ip_prec)) <= eps) then
                  met2(c3,ip_pt) = 1
                else if(dabs(met2(c3,ip_prec)) > eps.and.               &
                                            met2(c3,ip_tmp) > 0d0) then
                  met2(c3,ip_pt) = 2
                else if(dabs(met2(c3,ip_prec)) > eps.and.               &
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

                  if(met2(c3,ip_tsol) == 0d0) met2(c3,ip_upsol) = 0d0
                  met2(c3,ip_lcd) = cover(1)
                  met2(c3,ip_mcd) = cover(2)
                  met2(c3,ip_hcd) = cover(3)
                  met2(c3,ip_lhgt) = hgt(1)
                  met2(c3,ip_mhgt) = hgt(2)
                  met2(c3,ip_hhgt) = hgt(3)
                  met2(c3,ip_lct) = float(icld(1))
                  met2(c3,ip_mct) = float(icld(2))
                  met2(c3,ip_hct) = float(icld(3))     
                end if

                mcd = met2(c3,ip_mcd)
                hcd = met2(c3,ip_hcd)
                call emisatm(met2(c3,ip_tmp),met2(c3,ip_rh),ematm)

                call dnirflx(ematm,met2(c3,ip_tmp),met2(c3,ip_doy),     &
                             met2(c3,ip_lhgt),met2(c3,ip_mhgt),         &
                             met2(c3,ip_hhgt),met2(c3,ip_lcd),          &
                             mcd,hcd,met2(c3,ip_ir))

                c3 = c3 + 1
                if(ib /= 0) then
                  ib = ib + 1
                  if(ib > c2) ib = ib0
                end if
              end do   !i=1,msteps
            end do

! end with new data
            i = 1
            do while(i <= iend)
              do k=1,maxcol
                met2(c3,k) = met(i,k)
              end do
              c3 = c3 + 1
              i = i + 1
            end do
            c3 = c3 - 1

! see if there are too many timesteps being kept
            keep = iend + nlines1 + msteps                              !int(ndays*tstps)
            if(c3 > keep) then
              wstart = c3 - keep
              nlines2 = keep
            else
              wstart = 1
              nlines2 = c3
            end if
          end if   !if(vitd_index == vitd_index1)
        end if   !if(nlines1 > 0) then
      else
        err = 0
        allocate(met2(nlines,pcols),stat=err)                            !size met2 array

        do i=1,iend
          do j=1,pcols
            met2(i,j) = mflag
          end do
        end do

        wstart = 1
        c3 = iend
        nlines2 = iend
        do i=1,iend
          do j=1,pcols
            met2(i,j) = met(i,j)
          end do
        end do
      end if   !if(infer_test == 1.and.io /= -1) then

      write(2,1000) nlines2,maxcol-2,lat,mlong,elev,vitd_index,mcount,  &
                    timeoffset,timstep,mflag,iheight

! write the column headers
      write(2,'(''Year  JD  Hr   M   APres   ATemp     RH     WnSp   '',&
     &'' WDir    Prec   PT   Prec2  PT2   LCAmt    LCHt  LCT   MCAmt '',&
     &''   MCHt  MCT   HCAmt    HCHt  HCT     STot     SDir     SDif '',&
     &''   UpSol      IR     IRUp     SZen      SAz      Snow     '',   &
     &''Grtemp       Vis       Aer'')')    
      write(2,'(''                    mbar     C       %       m/s   '',&
     &''       mm/stp       mm/stp                 km                '',&
     &''   km                   km           W/m^2    W/m^2    W/m^2'', &
     &''    W/m^2    W/m^2    W/m^2                         m       '', &
     &''   K         km'')')

      do i=wstart,c3
        write(2,1001) int(met2(i,ip_year)),int(met2(i,ip_doy)),         &
          int(met2(i,ip_hr)),int(met2(i,ip_min)),met2(i,ip_ap),         &
          met2(i,ip_tmp),met2(i,ip_rh),met2(i,ip_ws),met2(i,ip_wdir),   &
          met2(i,ip_prec),int(met2(i,ip_pt)),met2(i,ip_prec2),          &
          int(met2(i,ip_pt2)),met2(i,ip_lcd),met2(i,ip_lhgt),           &
          int(met2(i,ip_lct)),met2(i,ip_mcd),met2(i,ip_mhgt),           &
          int(met2(i,ip_mct)),met2(i,ip_hcd),met2(i,ip_hhgt),           &
          int(met2(i,ip_hct)),met2(i,ip_tsol),met2(i,ip_dir),           &
          met2(i,ip_dif),met2(i,ip_upsol),met2(i,ip_ir),met2(i,ip_irup),&
          met2(i,ip_zen),met2(i,ip_az),met2(i,ip_sd),met2(i,ip_tsoil),  &
          met2(i,ip_vis),int(met2(i,ip_aer))
      end do

 1000 format(i6,1x,i3,1x,f7.3,1x,f9.3,1x,f9.3,1x,i10,1x,i6,1x,f6.2,1x,  &
             f5.2,1x,f8.2,1x,f8.3)

 1001 format(4(i4),5(f8.2),2(f8.3,i5),3(2(f8.2),i5),4(1x,f8.2),         &
            2(1x,f8.2),2(1x,f8.2),1x,f11.3,1x,f8.2,1x,f9.2,5x,i5)

      if(infer_test == 1.and.nlines1 /= 0) then
        deallocate(met1,met2)
      else
        deallocate(met2)
      end if

      end subroutine met_out

! ******************************************************************************
      subroutine met_date(diy,year,doy,hr,minute,mdate)

!     calculate the decimal calendar date

      implicit none

      integer(kind=4),intent(out):: diy
      real(kind=8),intent(in):: year,doy,hr,minute
      real(kind=8),intent(out):: mdate

!internal variables
      real(kind=8):: mody

      mdate = 0d0
      mody = year - aint(year*2.5d-1)*4d0

      if(mody <= epsilon(1d0)) then
        mdate = year + doy/366d0 + hr/(24d0*366d0)                      &
                                           + minute/(6d1*24d0*366d0)
        diy = 366
      else
        mdate = year + doy/365d0 + hr/(24d0*365d0)                      & 
                                           + minute/(6d1*24d0*365d0)
        diy = 365
      end if

      end subroutine met_date

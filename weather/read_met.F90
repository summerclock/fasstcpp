      subroutine read_met(shdrlines,rst_cutoff,met)

      use met_global
      use module_radiation

! calls the following subroutines:
!     cloudbase (in radiation module)
! this subroutine uses the function rh (appended to this subroutine)

! Read  met data. The pointers to the columns in the met file are in the 
! include file. The basic approach is to read in all the observations, check
! for missing data, and parse out the data to the array met.  For certain missing 
! parameters we will use climatological values or other values that are appropriate.

      implicit none

      integer(kind=4),intent(in):: shdrlines
      real(kind=8),intent(in):: rst_cutoff
      real(kind=8),intent(inout):: met(nlines,maxcol)

! local variables
      integer(kind=4),parameter:: no_data = -901
      integer(kind=4):: i,j,rhflag,ftime,io,havep,linel,lai,laf,loi,lof
      integer(kind=4):: eli,elf,vii,vif
      real(kind=8):: in_data(mxcol),temp_prec(nlines),tmp(mxcol)
      real(kind=8):: dewpt(nlines),rh1,rh2,rh,zlcld,zmcld,zhcld
      character(len=300):: header
      

! initialize variables
      rhflag = 0
      ftime = 0
      havep = 0
      lai = 0
      laf = 0
      loi = 0
      lof = 0
      eli = 0
      elf = 0
      vii = 0
      vif = 0
      rh1 = 0d0
      rh2 = 0d0

      do i=1,mxcol
        in_data(i) = 0d0
        tmp(i) = 0d0
      end do

      do i=1,nlines
        temp_prec(i) = 0d0
        dewpt(i) = 0d0
      end do

! if multi-point file, read two lines of labeling info
      if(shdrlines > 0) then
        do i=1,shdrlines
          read(1,'(a)') header

          if(i == 1) then
            j = 1
            do while(j <= len_trim(header))
              if(header(j:j) == char(76).or.                            &
                                         header(j:j) == char(108)) then  !'L' or 'l'
                if(header(j+1:j+1) == char(65).or.                      &
                                      header(j+1:j+1) == char(97)) then  !'A' or 'a'
                  lai = j
                  do while(header(j:j) /= char(32))
                    j = j + 1 
                  end do
                  laf = j - 1
                else if(header(j+1:j+1) == char(79).or.                 &
                                     header(j+1:j+1) == char(111)) then  !'O' or 'o'
                  loi = j
                  do while(header(j:j) /= char(32))
                    j = j + 1 
                  end do
                  lof = j - 1
                end if
              else if(header(j:j) == char(69).or.                       &
                                         header(j:j) == char(101)) then  !'E' or 'e'
                eli = j
                do while(header(j:j) /= char(32))
                  j = j + 1 
                end do
                elf = j - 1
              else if(header(j:j) == char(77).or.                       &
                                         header(j:j) == char(109)) then  !'M' or 'm'
                vii = j
                do while(header(j:j) /= char(32))
                  j = j + 1 
                end do
                vif = j - 1
              else
                j = j + 1
              end if
            end do

            read(header(laf+1:loi-1),*) lat
            lat = aint(lat*1d10)*1d-10
            read(header(lof+1:eli-1),*) mlong
            mlong = aint(mlong*1d10)*1d-10
            read(header(elf+1:vii-1),*) elev
            elev = aint(elev*1d10)*1d-10
            read(header(vif+1:len_trim(header)),*) vitd_index
          end if
        end do
      end if

!      do icnt=1,nlines
      io = 0
      icnt = 1
      do while(icnt <= nlines.and.io /= -1)
        read(1,*)(in_data(i),i=1,mxcol)                                  !read in a line of met data

! Year
        if(y_col /= no_data) met(icnt,ip_year) = in_data(y_col)
        if(y_col == no_data) met(icnt,ip_year) = year0
        if(met(icnt,ip_year) < 100d0.and.year0 < 2000d0)                & 
           met(icnt,ip_year) = met(icnt,ip_year) + 1900d0
        if(met(icnt,ip_year) < 100d0.and.year0 >= 2000d0)               &
           met(icnt,ip_year) = met(icnt,ip_year) + 2000d0

! Day Of Year	
        if(jday_col == no_data) then
          write(*,'('' No DOY information, STOPPING'')')
!          call exit(-1)
          stop
        else
          met(icnt,ip_doy) = in_data(jday_col)

          if(icnt /= 1) then
            if(met(icnt,ip_doy) == 1.and.(met(icnt-1,ip_doy) == 365.or. &
                                        met(icnt-1,ip_doy) == 366)) then
              if(met(icnt,ip_year) == met(icnt-1,ip_year)) then
                if(y_col == no_data) then
                  if(ftime == 0) then
                    year0 = year0 + 1
                    ftime = 1
                  end if
                  met(icnt,ip_year) = year0
                end if
              end if
            end if
          end if
        end if

! Local time of simulation
        if(hr_col == no_data) then
          write(*,'('' No hour information, STOPPING'')')
!          call exit(-1)
          stop
        else if(hr_col /= no_data.and.m_col /= no_data) then
          met(icnt,ip_hr) = in_data(hr_col)
          met(icnt,ip_min)= in_data(m_col)
        else if(hr_col /= no_data.and.m_col == no_data) then
          if(in_data(hr_col) > 24d0) then
            met(icnt,ip_hr) = aint(in_data(hr_col)/100d0)
            met(icnt,ip_min) = dmod(in_data(hr_col),100d0)
          else
            met(icnt,ip_hr) = aint(in_data(hr_col)/1d0)
            met(icnt,ip_min) = dmod(in_data(hr_col),1d0)
          end if
        end if

        if(icnt /= 1) then
          if(met(icnt,ip_hr) == 24d0.and.dabs(met(icnt,ip_doy)-         &
                                       met(icnt-1,ip_doy)) <= eps) then
            met(icnt,ip_doy) = met(icnt,ip_doy) + 1
            met(icnt,ip_hr) = 0d0
          end if
          if(met(icnt,ip_hr) == 0d0.and.dabs(met(icnt,ip_doy)-          &
                                       met(icnt-1,ip_doy)) <= eps) then
            met(icnt,ip_doy) = met(icnt,ip_doy) + 1
          end if
        end if

! air pressure (mbar)
        if(ap_col /= no_data) then
          met(icnt,ip_ap) = in_data(ap_col)
          if(apflag == 'pa'.and.dabs(met(icnt,ip_ap)-mflag) > eps)      &
             met(icnt,ip_ap) = met(icnt,ip_ap)*1d-2
          if(dabs(met(icnt,ip_ap)-mflag) <= eps) then
            if(icnt == 1) then
              met(icnt,ip_ap) = 1d3
              write(10,'(''1     used default air pressure of 1000.0 '',&
                    &''mbar, met_id = '',i10)') vitd_index
            else
              met(icnt,ip_ap) = met(icnt-1,ip_ap)
            end if
          end if
        else
          met(icnt,ip_ap) = 1d3
          if(icnt == 1)write(10,'(''1     used default air pressure '', &
            &''of 1000.0 mbar, met_id = '',i17)') vitd_index
        end if

        if((met(icnt,ip_ap) < 0d0).or.(met(icnt,ip_ap) > 1500d0)) then
          write(10,'('' Air pressure out of range, met_id = '',i37)'    &
                &) vitd_index
          write(10,'(''   0 < air pressure < 1500 mbar, your value is'',&
                &'':                 '',f10.5)') met(icnt,ip_ap)
          write(10,'(''   On day '',i5,'', hour '',i3)')                &
                int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
          write(10,'(''   used default air pressure of 1000.0 '',       &
                &''mbar, met_id = '',i20)') vitd_index
          met(icnt,ip_ap) = 1d3
        end if

! Check for valid air temperature(C) information
        if(tmp_col /= no_data) then
          if(ttflag == '1') then             !regular air temp data
            met(icnt,ip_tmp) = in_data(tmp_col)
            if(tflag == 'tk'.and.dabs(met(icnt,ip_tmp)-mflag) > eps)    &
               met(icnt,ip_tmp) = met(icnt,ip_tmp) - 273.15d0
            if(tflag == 'tf'.and.dabs(met(icnt,ip_tmp)-mflag) > eps)    &
               met(icnt,ip_tmp) = (met(icnt,ip_tmp) - 32d0)*5d0/9d0
            if(dabs(met(icnt,ip_tmp)-mflag) <= eps) then
              if(icnt /= 1) met(icnt,ip_tmp) = met(icnt-1,ip_tmp)
              if(icnt == 1) then
                write(*,'('' No temperature information, STOPPING'')')
                stop
              end if
            end if
          else if(ttflag == '2') then        !mean air temp data
            met(icnt,ip_tmp) = in_data(tmp_col)
            if(tflag == 'tk'.and.dabs(met(icnt,ip_tmp)-mflag) > eps)    &
               met(icnt,ip_tmp) = met(icnt,ip_tmp) - 273.15d0
            if(tflag == 'tf'.and.dabs(met(icnt,ip_tmp)-mflag) > eps)    &
               met(icnt,ip_tmp) = (met(icnt,ip_tmp) - 32d0)*5d0/9d0
            if(dabs(met(icnt,ip_tmp)-mflag) <= eps) then
              if(icnt /= 1) met(icnt,ip_tmp)=met(icnt-1,ip_tmp)
              if(icnt == 1) then
                write(*,'('' No temperature information, STOPPING'')')
                stop
              end if
            end if
          end if
        else if(tmp1_col /= no_data.and.tmp2_col /= no_data) then        !max and min
          if(dabs(in_data(tmp1_col)-mflag) <= eps.or.                   &
             dabs(in_data(tmp2_col)-mflag) <= eps) then
            if(icnt /= 1) met(icnt,ip_tmp) = met(icnt-1,ip_tmp)
            if(icnt == 1) then
              write(*,'('' No temperature information, STOPPING'')')
              stop
            end if
          else
            met(icnt,ip_tmp)=(in_data(tmp1_col)+in_data(tmp2_col))/2d0
            if(tflag == 'tk'.and.dabs(met(icnt,ip_tmp)-mflag) > eps)    &
               met(icnt,ip_tmp) = met(icnt,ip_tmp) - 273.15d0
            if(tflag == 'tf'.and.dabs(met(icnt,ip_tmp)-mflag) > eps)    &
               met(icnt,ip_tmp) = (met(icnt,ip_tmp) - 32d0)*5d0/9d0
          end if
        else
          write(*,'('' No temperature information, STOPPING'')')
          stop
        end if

        if(dabs(met(icnt,ip_tmp)-mflag) <= eps)then
          write(*,'('' No temperature information, STOPPING'')')
          stop
        else if((met(icnt,ip_tmp) < -100d0).or.                         &
                                       (met(icnt,ip_tmp) > 100d0)) then
          write(*,'('' Air temperature out of range, STOPPING'')')
          write(*,'('' -100 < air temp < 100 C, your value is: '',      &
                &f10.3)')met(icnt,ip_tmp)
          write(*,'('' On day '',i5,'' hour '',i3)')                    &
                int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
          stop
        end if

!   Check for valid RH information otherwise set to 60% 
        if(rh_col /= no_data) then
          met(icnt,ip_rh) = in_data(rh_col)
          if((met(icnt,ip_rh) >= 0d0.and.met(icnt,ip_rh) <= 1d0).and.   &
                                     dabs(met(icnt,ip_rh)-mflag) > eps) &
             met(icnt,ip_rh) = met(icnt,ip_rh)*100d0
          if(dabs(met(icnt,ip_rh)-mflag) <= eps) then
            if(icnt == 1) then
              met(icnt,ip_rh) = 60d0
              write(10,'(''1     used default relative humidity of '',  &
                    &''60 %, met_id = '',i19)') vitd_index
            else
              met(icnt,ip_rh) = met(icnt-1,ip_rh)
            end if
          end if
        else if(dt_col /= no_data) then
          rhflag = 0
          dewpt(icnt) = in_data(dt_col)
          if(dtflag == 'tk'.and.dabs(dewpt(icnt)-mflag) > eps)          & 
             dewpt(icnt) = dewpt(icnt) - 273.15d0
          if(dtflag == 'tf'.and.dabs(dewpt(icnt)-mflag) > eps)          &
             dewpt(icnt) = (dewpt(icnt) - 32d0)*5d0/9d0
          if(dabs(dewpt(icnt)-mflag) <= eps) then
            if(icnt /= 1) dewpt(icnt) = dewpt(icnt-1)
            if(icnt == 1) then
              met(icnt,ip_rh) = 60d0
              rhflag = 1
            end if
          end if
          if(rhflag == 0) then
            rh1 = rh(dewpt(icnt))
            rh2 = rh(met(icnt,ip_tmp))
            met(icnt,ip_rh) = (rh1/rh2)*100d0
          end if
        else if(dd_col /= no_data) then
          rhflag = 0
          dewpt(icnt) = in_data(dd_col)
          if(ddflag == 'tk'.and.dabs(dewpt(icnt)-mflag) > eps)          &
             dewpt(icnt) = dewpt(icnt) - 273.15d0
          if(ddflag == 'tf'.and.dabs(dewpt(icnt)-mflag) > eps)          &
             dewpt(icnt) = (dewpt(icnt) - 32d0)*5d0/9d0
          if(dabs(dewpt(icnt)-mflag) <= eps) then
            if(icnt /= 1) dewpt(icnt) = dewpt(icnt-1)
            if(icnt == 1) then
              met(icnt,ip_rh) = 60d0
              rhflag = 1
            end if
          end if
          if(rhflag == 0) then
            dewpt(icnt) = met(icnt,ip_tmp) - dewpt(icnt)
            rh1 = rh(dewpt(icnt))
            rh2 = rh(met(icnt,ip_tmp))
            met(icnt,ip_rh) = (rh1/rh2)*1d2
          end if
        else
          met(icnt,ip_rh) = 60d0
          if(icnt == 1) write(10,'(''1     used default relative '',    &
            &''humidity of 60 %, met_id = '',i19)') vitd_index
        end if

        if(dabs(met(icnt,ip_rh)-mflag) > eps) then
          if((met(icnt,ip_rh) < 0d0).or.(met(icnt,ip_rh) > 1d2)) then
            write(10,'('' Relative humidity out of range, '',           &
                  &''met_id = '',i32)') vitd_index
            write(10,'(''   0 < RH < 100, your value is:            '', &
                  &''                     '',f10.5)') met(icnt,ip_rh)
            write(10,'(''   On day '',i5,'', hour '',i3)')              &
                  int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
            write(10,'(''   used default relative humidity of '',       &
                  &''60 %, met_id = '',i22)') vitd_index
            met(icnt,ip_rh) = 60d0
          end if
        end if

! If the wind information is missing set the value to 1.2 m/sec default for SNAP
        if(ws_col /= no_data) then
          met(icnt,ip_ws) = in_data(ws_col)
          if(wsflag == 'fps'.and.dabs(met(icnt,ip_ws)-mflag) > eps)     &
             met(icnt,ip_ws) = met(icnt,ip_ws)*0.3048d0
          if(dabs(met(icnt,ip_ws)-mflag) <= eps) then
            if(icnt == 1) then
              met(icnt,ip_ws) = 1.2d0
              write(10,'(''1     used default wind speed of 1.2 m/s, '',&
                    &''met_id = '',i23)') vitd_index
            else
              met(icnt,ip_ws) = met(icnt-1,ip_ws)
            end if
          end if
        else
          met(icnt,ip_ws) = 1.2d0
          if(icnt == 1) write(10,'(''1     used default wind speed '',  &
            &''of 1.2 m/s, met_id = '',i23)') vitd_index
        end if

        if((met(icnt,ip_ws) < 0d0).or.(met(icnt,ip_ws) > 110d0)) then
          write(10,'('' Wind speeds out of range, met_id = '',i38)')    &
                vitd_index
          write(10,'(''   0 < WS < 110 m/s, your value is:          '', &
                &''                   '',f10.5)') met(icnt,ip_ws)
          write(10,'(''   On day '',i5,'', hour '',i3)')                &
                             int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
          write(10,'(''   used default wind speed of 2.1 m/s, '',       &
                &''met_id = '',i26)') vitd_index
          met(icnt,ip_ws) = 2.1d0
        end if

! Calculate the 2D wind
        if(wd_col /= no_data) then
          met(icnt,ip_wdir) = in_data(wd_col)
          if(dabs(met(icnt,ip_wdir)-mflag) <= eps) then
            if(icnt /= 1) met(icnt,ip_wdir) = met(icnt-1,ip_wdir)
            if(icnt == 1) met(icnt,ip_wdir) = 0d0
          end if
          if((met(icnt,ip_wdir) < 0d0).or.                              &
                                      (met(icnt,ip_wdir) > 360d0)) then
            write(10,'('' Wind direction out of range, met_id = '',     &
                  &i35)') vitd_index
            write(10,'(''   0 < wind dir < 360, your value is:       '',&
                  &''                    '',f10.3)') met(icnt,ip_wdir)
            write(10,'(''   On day '',i5,'', hour '',i3)')              &
                             int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
            write(10,'(''   used default wind direction of 0.0, '',     &
                  &''met_id = '',i26)') vitd_index
            met(icnt,ip_wdir) = 0d0
          end if
        end if

! Check for valid prec (mm/h) information otherwise set to 0 
        if(prcp_col /= no_data.and.prcp_flag == '1') then
          met(icnt,ip_prec) = in_data(prcp_col)
          if(pflag == 'iph'.and.dabs(met(icnt,ip_prec)-mflag) > eps)    &
             met(icnt,ip_prec) = met(icnt,ip_prec)*25.4d0
          if(pflag == 'fph'.and.dabs(met(icnt,ip_prec)-mflag) > eps)    &
             met(icnt,ip_prec) = met(icnt,ip_prec)*304.8d0
          if(pflag == 'mph'.and.dabs(met(icnt,ip_prec)-mflag) > eps)    &
             met(icnt,ip_prec) = met(icnt,ip_prec)*1000d0
          if(dabs(met(icnt,ip_prec)-mflag) <= eps) then
            met(icnt,ip_prec) = 0d0
            if(icnt == 1) write(10,'(''1     used default precipita'',  &
              &''tion rate of 0 mm/step, met_id = '',i13)') vitd_index
          end if
        else if(pa_col /= no_data) then
          temp_prec(icnt) = in_data(pa_col)
          if(paflag == 'i'.and.dabs(temp_prec(icnt)-mflag) > eps)       &
             temp_prec(icnt) = temp_prec(icnt)*25.4d0
          if(paflag == 'f'.and.dabs(temp_prec(icnt)-mflag) > eps)       &
             temp_prec(icnt) = temp_prec(icnt)*304.8d0
          if(paflag == 'a'.and.dabs(temp_prec(icnt)-mflag) > eps)       &
             temp_prec(icnt) = temp_prec(icnt)*1000d0
          if(paflag == 'c'.and.dabs(temp_prec(icnt)-mflag) > eps)       &
             temp_prec(icnt) = temp_prec(icnt)*10d0

          if(dabs(temp_prec(icnt)-mflag) <= eps) then
            temp_prec(icnt) = 0d0
            if(icnt == 1) write(10,'(''1     used default precipita'',  &
              &''tion rate of 0 mm/step, met_id = '',i13)') vitd_index
          end if
!          temp_prec(icnt) = temp_prec(icnt)/timstep
          met(icnt,ip_prec) = temp_prec(icnt)
          if(icnt > 1) then
            if(temp_prec(icnt)-temp_prec(icnt-1) >= 0d0) then
              met(icnt,ip_prec) = temp_prec(icnt) - temp_prec(icnt-1)
            end if
          end if
        else
          met(icnt,ip_prec) = 0d0
          if(icnt == 1) write(10,'(''1     used default precipitation'',&
            &'' rate of 0 mm/step, met_id = '',i13)') vitd_index
        end if

! correct negative values from the above procedure
        if(met(icnt,ip_prec) < 0d0) met(icnt,ip_prec) = 0d0

! Precipitation Type
        if(pt_col /= no_data) then                                       !precip type
          met(icnt,ip_pt) = aint(in_data(pt_col))
!	    if(aint(met(icnt,ip_pt)) == 0) then         !note, these fixes are for grayling, yuma and soro data
!	      met(icnt,ip_pt) = 1d0   !none
!	    else if(aint(met(icnt,ip_pt)) == 1) then
!	     met(icnt,ip_pt) = 2d0   !rain
!	    else if(aint(met(icnt,ip_pt)) == 2.or.                        &
!                                 aint(met(icnt,ip_pt)) == -1) then   &
!            met(icnt,ip_pt) = 3d0   !snow
!	    else if(aint(met(icnt,ip_pt)) > 2.and.                        &
!                               dabs(met(icnt,ip_prec)) <= eps) then  &
!            met(icnt,ip_pt) = 1d0
          if(dabs(met(icnt,ip_pt)-mflag) <= eps.or.                     &
                             aint(met(icnt,ip_pt)) == aint(mflag)) then
            if(dabs(met(icnt,ip_prec)) <= eps) then
              met(icnt,ip_pt) = 1d0                                      !none
            else
              if(met(icnt,ip_tmp) > rst_cutoff) met(icnt,ip_pt) = 2d0    !rain
              if(met(icnt,ip_tmp) <= rst_cutoff) met(icnt,ip_pt) = 3d0   !snow
            end if
          end if
          if(aint(met(icnt,ip_pt)) == 0.and.                            &
            dabs(met(icnt,ip_prec)) <= eps) met(icnt,ip_pt) = 1
        else 
          if(dabs(met(icnt,ip_prec)) > eps.and.                         &
                   met(icnt,ip_tmp) > rst_cutoff) met(icnt,ip_pt) = 2d0  !rain
          if(dabs(met(icnt,ip_prec)) > eps.and.                         &
                   met(icnt,ip_tmp) <= rst_cutoff)met(icnt,ip_pt) = 3d0  !snow
          if(dabs(met(icnt,ip_prec)) <= eps) met(icnt,ip_pt) = 1d0       !none
        end if

! Check for valid snow prec (mm/h) information otherwise set to 0 
        if(prcp2_col /= no_data.and.prcp2_flag == '1') then
          met(icnt,ip_prec2) = in_data(prcp2_col)
          if(p2flag == 'iph'.and.dabs(met(icnt,ip_prec2)-mflag) > eps)  &
             met(icnt,ip_prec2) = met(icnt,ip_prec2)*25.4d0
          if(p2flag == 'fph'.and.dabs(met(icnt,ip_prec2)-mflag) > eps)  &
             met(icnt,ip_prec2) = met(icnt,ip_prec2)*304.8d0
          if(p2flag == 'mph'.and.dabs(met(icnt,ip_prec2)-mflag) > eps)  &
             met(icnt,ip_prec2) = met(icnt,ip_prec2)*1000d0
          if(dabs(met(icnt,ip_prec2)-mflag) <= eps) then
            met(icnt,ip_prec2) = 0d0
            if(icnt == 1) write(10,'(''1     used default snow precip'',&
              &''itation rate of 0 mm/step, met_id = '',i10)')          &
              vitd_index
            end if
        else if(pa2_col /= no_data) then
          temp_prec(icnt) = in_data(pa2_col)
          if(pa2flag == 'i'.and.dabs(temp_prec(icnt)-mflag) > eps)      &
             temp_prec(icnt) = temp_prec(icnt)*25.4d0
          if(pa2flag == 'f'.and.dabs(temp_prec(icnt)-mflag) > eps)      &
             temp_prec(icnt) = temp_prec(icnt)*304.8d0
          if(pa2flag == 'a'.and.dabs(temp_prec(icnt)-mflag) > eps)      &
             temp_prec(icnt) = temp_prec(icnt)*1000d0
          if(pa2flag == 'c'.and.dabs(temp_prec(icnt)-mflag) > eps)      &
             temp_prec(icnt) = temp_prec(icnt)*10d0
          if(dabs(temp_prec(icnt)-mflag) <= eps) then
            temp_prec(icnt) = 0d0
            if(icnt == 1) write(10,'(''1     used default snow precip'',&
              &''itation rate of 0 mm/step, met_id = '',i10)')          &
              vitd_index
          end if
!          temp_prec(icnt) = temp_prec(icnt)/timstep
          met(icnt,ip_prec2) = temp_prec(icnt)
          if(icnt > 1) then
            if(temp_prec(icnt)-temp_prec(icnt-1) >= 0d0) then
              met(icnt,ip_prec2) = temp_prec(icnt) - temp_prec(icnt-1)
            end if
          end if
        else
          met(icnt,ip_prec2) = 0d0
          if(icnt == 1) write(10,'(''1     used default snow precipit'',&
            &''ion rate of 0 mm/step, met_id = '',i10)') vitd_index
        end if

! correct negative values from the above procedure
        if(met(icnt,ip_prec2) < 0d0) met(icnt,ip_prec2) = 0d0

! Snow Precipitation Type
        if(pt2_col /= no_data) then                                      !precip type
          met(icnt,ip_pt2) = aint(in_data(pt2_col))
!	    if(aint(met(icnt,ip_pt2)) == 0) then         !note, these fixes are for grayling, yuma and soro data
!	      met(icnt,ip_pt2) = 1d0   !none
!	    else if(aint(met(icnt,ip_pt2)) == 2.or.                         &
!                                   aint(met(icnt,ip_pt2)) == -1) then  &
!            met(icnt,ip_pt) = 3d0   !snow
!	    else if(aint(met(icnt,ip_pt2)) > 2.and.                         &
!                                dabs(met(icnt,ip_prec2)) <= eps) then  &
!            met(icnt,ip_pt2) = 1d0
          if(dabs(met(icnt,ip_pt2)-mflag) <= eps.or.                    &
             dabs(aint(met(icnt,ip_pt2))-aint(mflag)) <= eps) then
            if(dabs(met(icnt,ip_prec2)) <= eps) then
              met(icnt,ip_pt2) = 1d0                                     !none
            else
              met(icnt,ip_pt2) = 3d0                                     !snow
            end if
          end if
          if(aint(met(icnt,ip_pt)) == 0.and.                            &
     &      dabs(met(icnt,ip_prec)) <= eps) met(icnt,ip_pt) = 1
        else 
          if(dabs(met(icnt,ip_prec2)) > eps) met(icnt,ip_pt2) = 3d0      !snow
          if(dabs(met(icnt,ip_prec2)) <= eps) met(icnt,ip_pt2) = 1d0     !none
        end if

        havep = 0
        if(dabs(met(icnt,ip_pt)-1d0) > eps.or.                          &
                                 dabs(met(icnt,ip_pt2)-1d0) > eps) then
          havep = 1
        end if
    
! Now check to see if this data set reported cloud information.  If not set to
! a climatological value for the low middle, and high cloud amount, hgt and type.
        if(lcld_col /= no_data) then                                     !low cloud amount
          met(icnt,ip_lcd) = in_data(lcld_col)
          if(dabs(met(icnt,ip_lcd)-mflag) <= eps) then
            if(icnt == 1) then
              if(havep == 0) then
                met(icnt,ip_lcd) = 0.5d0
              else if(havep == 1) then
                met(icnt,ip_lcd) = 1d0
              end if

              write(10,'(''1     used default low cloud amount of 0.5'',&
                    &'' (no precipitation)'',/,''        or 1.0 '',     &
                    &''(precipitation), met_id = '',i33)') vitd_index
            else
              if(dabs(tmp(lcld_col)-mflag) <= eps) then
                if(havep == 0) then
                  met(icnt,ip_lcd) = 0d0
                else if(havep == 1) then
                  met(icnt,ip_lcd) = 1d0
                end if
              else
                met(icnt,ip_lcd) = met(icnt-1,ip_lcd)
                if(havep == 1) met(icnt,ip_lcd) = 1d0
              end if
            end if
          end if
        else
          if(havep == 0) then
            met(icnt,ip_lcd) = 0.5d0
          else if(havep == 1) then
            met(icnt,ip_lcd) = 1d0
          end if

          if(icnt == 1) write(10,'(''1     used default low cloud '',   &
            &''amount of 0.5 (no precipitation)'',/,''        or 1.0 '',&
            &''(precipitation), met_id = '',i33)') vitd_index
        end if
        met(icnt,ip_lcd) = dmax1(0d0,dmin1(1d0,met(icnt,ip_lcd)))

        if(lhgt_col /= no_data)then                                      !low cloud height
          met(icnt,ip_lhgt) = in_data(lhgt_col)
          if(dabs(met(icnt,ip_lhgt)-mflag) > eps) then
            if(lcflag == 'm')                                           &
              met(icnt,ip_lhgt) = met(icnt,ip_lhgt)*1.609344d0
          else
            write(10,'(''2     model calculated low cloud height, '',   &
                  &''met_id = '',i24)') vitd_index
          end if
        else
          if(icnt == 1) write(10,'(''2     model calculated low cloud'',&
                              &'' height, met_id = '',i24)') vitd_index
        end if

! if no low cloud type information, assume stratus
        if(lct_col /= no_data) then                                      !low cloud type
          met(icnt,ip_lct) = in_data(lct_col)
          if(dabs(met(icnt,ip_lct)-mflag) <= eps.or.                    &
             dabs(aint(met(icnt,ip_lct))-aint(mflag)) <= eps) then
            if(icnt == 1) then
              if(met(icnt,ip_lcd) <= eps) met(icnt,ip_lct) = 0d0
              if(met(icnt,ip_lcd) > eps) met(icnt,ip_lct) = 6d0
              write(10,'(''1     used default low cloud type of 6, '',  &
                    &''stratus nebulosus'',/,''        and/or startus'',&
                    &'' fractus, if clouds present, met_id = '',i14)')  &
                                                             vitd_index
            else
              if(met(icnt,ip_lcd) > eps) then
                if(met(icnt-1,ip_lct) > eps) then
                  met(icnt,ip_lct) = met(icnt-1,ip_lct)
                else
                  met(icnt,ip_lct) = 6d0
                end if
              else
                met(icnt,ip_lct) = 0d0
              end if
            end if
          end if
        else
          if(met(icnt,ip_lcd) <= eps) met(icnt,ip_lct) = 0d0
          if(met(icnt,ip_lcd) > eps) met(icnt,ip_lct) = 6d0
          if(icnt == 1) write(10,'(''1     used default low cloud '',   &
            &''type of 6, stratus nebulosus'',/,''        and/or '',    &
            &''startus fractus, if clouds present, met_id = '',i14)')   &
                                                             vitd_index
        end if

! middle clouds
        if(mcld_col /= no_data) then                                     !middle cloud amount
          met(icnt,ip_mcd) = in_data(mcld_col)
          if(dabs(met(icnt,ip_mcd)-mflag) <= eps) then
            if(icnt == 1) then
              if(havep == 0) then
                met(icnt,ip_mcd) = 0d0
              else if(havep == 1) then
                met(icnt,ip_mcd) = 1d0
              end if

              write(10,'(''1     used default middle cloud amount of '',&
                    &''0.0 (no precipitation)'',/''        or 1.0 '',   &
                    &''(precipitation), met_id = '',i33)') vitd_index
            else
              if(dabs(tmp(mcld_col)-mflag) <= eps) then
                if(havep == 0) then
                  met(icnt,ip_mcd) = 0d0
                else if(havep == 1) then
                  met(icnt,ip_mcd) = 1d0
                end if
              else
                met(icnt,ip_mcd) = met(icnt-1,ip_mcd)
                if(havep == 1) met(icnt,ip_mcd) = 1d0
              end if
            end if
          end if
        else
          if(havep == 0) then
            met(icnt,ip_mcd) = 0d0
          else if(havep == 1) then
            met(icnt,ip_mcd) = 1d0
          end if

          if(icnt == 1) write(10,'(''1     used default middle cloud'', &
            &'' amount of 0.0 (no precipitation)'',/,''        or 1.0'',&
            &'' (precipitation), met_id = '',i33)') vitd_index
        end if
        met(icnt,ip_mcd) = dmax1(0d0,dmin1(1d0,met(icnt,ip_mcd)))

        if(mhgt_col /= no_data)then                                      !middle cloud height
          met(icnt,ip_mhgt) = in_data(mhgt_col)
          if(dabs(met(icnt,ip_mhgt)-mflag) > eps) then
            if(mcflag == 'm')                                           &
              met(icnt,ip_mhgt) = met(icnt,ip_mhgt)*1.609344d0
          else
            write(10,'(''2     model calculated middle cloud height, '',&
                  &''met_id = '',i21)') vitd_index
          end if
        else
          if(icnt == 1) write(10,'(''2     model calculated middle '',  &
                              &''cloud height, met_id = '',i21)')       &
                                                             vitd_index
        end if

        if(mct_col /= no_data) then                                      !middle cloud type
          met(icnt,ip_mct) = in_data(mct_col)
          if(dabs(met(icnt,ip_mct)-mflag) <= eps.or.                    &
             dabs(aint(met(icnt,ip_mct))-aint(mflag)) <= eps) then
            if(icnt == 1) then
              if(met(icnt,ip_mcd) <= eps) met(icnt,ip_mct) = 0d0
              if(met(icnt,ip_mcd) > 0d0) met(icnt,ip_mct) = 3d0
              write(10,'(''1     used default middle cloud type of 3,'',&
                    &'' altocumulus'',/,''        translucidus, 1 '',   &
                    &''level, if clouds present, met_id = '',i10)')     &
                                                             vitd_index
            else
              if(met(icnt,ip_mcd) > eps) then
                if(met(icnt-1,ip_mct) > eps) then
                  met(icnt,ip_mct) = met(icnt-1,ip_mct)
                else
                  met(icnt,ip_mct) = 3d0
                end if
              else
                met(icnt,ip_mct) = 0d0
              end if
            end if
          end if
        else
          if(met(icnt,ip_mcd) <= eps) met(icnt,ip_mct) = 0d0
          if(met(icnt,ip_mcd) > 0d0) met(icnt,ip_mct) = 3d0              !alto cumulus, 1 layer
          if(icnt == 1) write(10,'(''1     used default middle cloud'', &
            &'' type of 3, altocumulus'',/,''        translucidus, 1 '',&
            &''level, if clouds present, met_id = '',i15)') vitd_index
        end if

! high clouds
        if(hcld_col /= no_data) then                                     !high cloud amount
          met(icnt,ip_hcd) = in_data(hcld_col)
          if(dabs(met(icnt,ip_hcd)-mflag) <= eps) then
            if(icnt == 1) then
              if(havep == 0) then
                met(icnt,ip_hcd) = 0d0
              else if(havep == 1) then
                met(icnt,ip_hcd) = 1d0
              end if

              write(10,'(''1     used default high cloud amount of '',  &
                    &''0.0 (no precipitation)'',/,''        or 1.0 '',  &
                    &''(precipitation), met_id = '',i33)') vitd_index
            else
              if(dabs(tmp(hcld_col)-mflag) <= eps) then
                if(havep == 0) then
                  met(icnt,ip_hcd) = 0d0
                else if(havep == 1) then
                  met(icnt,ip_hcd) = 1d0
                end if
              else
                met(icnt,ip_hcd) = met(icnt-1,ip_hcd)
                if(havep == 1) met(icnt,ip_hcd) = 1d0
              end if
            end if
          end if
        else
          if(havep == 0) then
            met(icnt,ip_hcd) = 0d0
          else if(havep == 1) then
            met(icnt,ip_hcd) = 1d0
          end if

          if(icnt == 1) write(10,'(''1     used default high cloud '',  &
            &''amount of 0.0 (no precipitation)'',/,''        or '',    &
            &''1.0 (precipitation), met_id = '',i33)') vitd_index 
        end if
        met(icnt,ip_hcd) = dmax1(0d0,dmin1(1d0,met(icnt,ip_hcd)))

        if(hhgt_col /= no_data)then                                      !high cloud height
          met(icnt,ip_hhgt) = in_data(hhgt_col)
          if(dabs(met(icnt,ip_hhgt)-mflag) > eps) then
            if(hcflag == 'm')                                           &
              met(icnt,ip_hhgt) = met(icnt,ip_hhgt)*1.609344d0
          else
            write(10,'(''2     model calculated high cloud height, '',  &
                  &''met_id = '',i23)') vitd_index
          end if
        else
          if(icnt == 1) write(10,'(''2     model calculated high '',    &
                              &''cloud height, met_id = '',i23)')       &
                                                             vitd_index
        end if

        if(hct_col /= no_data) then                                      !high cloud type
          met(icnt,ip_hct) = in_data(hct_col)
          if(dabs(met(icnt,ip_hct)-mflag) <= eps.or.                    &
             dabs(aint(met(icnt,ip_lct))-aint(mflag)) <= eps) then
            if(icnt == 1) then
              if(met(icnt,ip_hcd) <= eps) met(icnt,ip_hct) = 0d0
              if(met(icnt,ip_hcd) > 0d0) met(icnt,ip_hct) = 5d0
              write(10,'(''1     used default high cloud type of 5, '', &
                    &''cirrus and/or '',/,''        cirrostratus < '',  &
                    &''45 above horiz., if clouds present, met_id = '', &
                    &                                 i6)') vitd_index
            else
              if(met(icnt,ip_hcd) > eps) then
                if(met(icnt-1,ip_hct) > eps) then
                  met(icnt,ip_hct) = met(icnt-1,ip_hct)
                else
                  met(icnt,ip_hct) = 5d0
                end if
              else
                met(icnt,ip_hct) = 0d0
              end if
            end if
          end if
        else
          if(met(icnt,ip_hcd) <= eps) met(icnt,ip_hct) = 0d0
          if(met(icnt,ip_hcd) > 0d0) met(icnt,ip_hct) = 5d0              !cirrus and/or cirrostratus
          if(icnt == 1) write(10,'(''1     used default high cloud '',  &
            &''type of 5, cirrus and/or '',/,''        cirrostratus <'',&
            &'' 45 above horiz., if clouds present, met_id = '',i6)')   &
                                                             vitd_index
        end if

! calculate the cloud base heights
      call cloudbase(met(icnt,ip_doy),lat,met(icnt,ip_lhgt),            &
                     met(icnt,ip_mhgt),met(icnt,ip_hhgt),               &
                     met(icnt,ip_lcd),met(icnt,ip_mcd),met(icnt,ip_hcd),&
                     zlcld,zmcld,zhcld)
      met(icnt,ip_lhgt) = zlcld
      met(icnt,ip_mhgt) = zmcld
      met(icnt,ip_hhgt) = zhcld

!   Check for valid total solar flux information 
        if(tsol_col /= no_data) then
          met(icnt,ip_tsol) = in_data(tsol_col)
          if(dabs(met(icnt,ip_tsol)-mflag) <= eps) then
             met(icnt,ip_tsol) = mflag
!            if(icnt /= 1) met(icnt,ip_tsol) = met(icnt-1,ip_tsol)
!            if(icnt == 1) met(icnt,ip_tsol) = mflag
          end if
        else
          met(icnt,ip_tsol) = mflag
        end if
  
        if(met(icnt,ip_tsol) < 0d0.and.                                 &
          dabs(met(icnt,ip_tsol)-mflag) > eps) met(icnt,ip_tsol) = mflag !0d0

!   Check for valid direct solar flux information 
        if(dirsol_col /= no_data) then
          met(icnt,ip_dir) = in_data(dirsol_col)
          if(dabs(met(icnt,ip_dir)-mflag) <= eps) then
             met(icnt,ip_dir) = mflag
!            if(icnt /= 1) met(icnt,ip_dir) = met(icnt-1,ip_dir)
!            if(icnt == 1) met(icnt,ip_dir) = mflag
          end if
        else
          met(icnt,ip_dir) = mflag
        end if

        if(met(icnt,ip_dir) < 0d0.and.                                  &
          dabs(met(icnt,ip_dir)-mflag) > eps) met(icnt,ip_dir) = mflag !0d0

!   Check for valid diffuse solar flux information 
        if(difsol_col /= no_data) then
          met(icnt,ip_dif) = in_data(difsol_col)
          if(dabs(met(icnt,ip_dif)-mflag) <= eps) then
             met(icnt,ip_dif) = mflag
!            if(icnt /= 1) met(icnt,ip_dif) = met(icnt-1,ip_dif)
!            if(icnt == 1) met(icnt,ip_dif) = mflag
          end if
        else
          met(icnt,ip_dif) = mflag
        end if

        if(met(icnt,ip_dif) < 0d0.and.                                  &
          dabs(met(icnt,ip_dif)-mflag) > eps) met(icnt,ip_dif) = mflag !0d0

! make sure that total = direct + diffuse
!        if((dabs(met(icnt,ip_dif)-mflag) > eps.and.                     &
!                              dabs(met(icnt,ip_dir)-mflag) > eps).and.  &
!                              dabs(met(icnt,ip_tsol)-mflag) > eps) then
!          if(dabs(met(icnt,ip_tsol)-                                    &
!                       (met(icnt,ip_dir)+met(icnt,ip_dif))) > eps) then
!           if(met(icnt,ip_dir)+met(icnt,ip_dif) >  
!           met(icnt,ip_dif) = met(icnt,ip_tsol) - met(icnt,ip_dir)
!           if(met(icnt,ip_dif) < eps) met(icnt,ip_dif) = mflag
!        end if
!        if(met(icnt,ip_dir) < 0d0.and.                                  &
!                               dabs(met(icnt,ip_dir)-mflag) > eps) then
!          met(icnt,ip_dir) = 0d0
!          met(icnt,ip_dif) = met(icnt,ip_tsol)
!        end if

! If we have direct and diffuse flux, but not total set total to the sum
        if(dabs(met(icnt,ip_tsol)-mflag) <= eps.and.                    &
           (dabs(met(icnt,ip_dir)-mflag) > eps.and.                     &
            dabs(met(icnt,ip_dif)-mflag) > eps))                        & 
           met(icnt,ip_tsol) = met(icnt,ip_dir) + met(icnt,ip_dif)

        if(icnt == 1.and.dabs(met(icnt,ip_tsol)-mflag) <= eps)          &
          write(10,'(''2     model calculated total solar radiation, '',&
                &''met_id = '',i19)')  vitd_index
        if(icnt == 1.and.dabs(met(icnt,ip_dif)-mflag) <= eps)           &
          write(10,'(''2     model calculated diffuse solar radiation'',&
                &'', met_id = '',i17)') vitd_index
        if(icnt == 1.and.dabs(met(icnt,ip_dir)-mflag) <= eps)           &
          write(10,'(''2     model calculated direct solar radiation,'',&
                &'' met_id = '',i18)') vitd_index

        if(dabs(met(icnt,ip_tsol)-mflag) > eps) then
          if(met(icnt,ip_tsol) < 0d0.or.met(icnt,ip_tsol) > 1355d0) then
            write(10,'('' Incoming shortwave radiation out of range, '',&
                  &''met_id = '',i21)') vitd_index
            write(10,'(''   0 < incoming < 1355 W/m^2, your value is:'',&
                  &''                    '',f10.3)') met(icnt,ip_tsol)
            write(10,'(''   On day '',i5,'' hour '',i3)')               &
                  int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
            write(10,'(''   model calculated total solar radiat'',      &
                  &''ion, met_id = '',i22)') vitd_index
            met(icnt,ip_tsol) = mflag
          end if
        end if

        if(dabs(met(icnt,ip_dir)-mflag) > eps) then
          if((met(icnt,ip_dir) < 0d0).or.                               &
                                      (met(icnt,ip_dir) > 1355d0)) then
            write(10,'('' Direct shortwave radiation out of range, '',  &
                  &''met_id = '',i23)') vitd_index
            write(10,'(''   0 < direct < 1355 W/m^2, your value is: '', &
                  &''                     '',f10.3)') met(icnt,ip_dir)
            write(10,'(''   On day '',i5,'' hour '',i3)')               &
                  int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
            write(10,'(''   model calculated direct solar radiat'',     &
                  &''ion, met_id = '',i21)') vitd_index
            met(icnt,ip_dir) = mflag
          end if
        end if

        if(dabs(met(icnt,ip_dif)-mflag) > eps) then
          if((met(icnt,ip_dif) < 0d0).or.                               &
                                      (met(icnt,ip_dif) > 1355d0)) then
            write(10,'('' Diffuse shortwave radiation out of range, '', &
                  &''met_id = '',i22)') vitd_index
            write(10,'(''   0 < diffuse < 1355 W/m^2, your value is: '',&
                  &''                    '',f10.3)') met(icnt,ip_dif)
            write(10,'(''   On day '',i5,'' hour '',i3)')               &
                  int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
            write(10,'(''   model calculated diffuse solar radia'',     &
                  &''tion, met_id = '',i20)') vitd_index
            met(icnt,ip_dif) = mflag
          end if
        end if

!   Check for valid reflected solar flux information 
        if(upsol_col /= no_data) then
          met(icnt,ip_upsol) = in_data(upsol_col)
          if(dabs(met(icnt,ip_upsol)-mflag) <= eps) then
            if(icnt /= 1) met(icnt,ip_upsol) = met(icnt-1,ip_upsol)
            if(icnt == 1) met(icnt,ip_upsol) = mflag
          end if
          if(dabs(met(icnt,ip_upsol)-mflag) > eps) then
            if((met(icnt,ip_upsol) < 0d0).or.                           &
                                    (met(icnt,ip_upsol) > 1355d0)) then
              write(10,'('' Reflected shortwave radiation out of '',    &
                    &''range, met_id = '',i20)') vitd_index
              write(10,'(''   0 < reflected < 1355 W/m^2, your value '',&
                    &''is: '',f10.3)') met(icnt,ip_upsol)
              write(10,'(''   On day '',i5,'' hour '',i3)')             &
                    int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
              write(10,'(''   FASST will calculate reflected solar'',   &
                    &'' radiation, met_id = '',i14)') vitd_index
              met(icnt,ip_upsol) = mflag
            end if
          end if
        else
          met(icnt,ip_upsol) = mflag
          if(icnt == 1) write(10,'(''2     FASST will calculate '',     &
            &''reflected solar radiation, met_id = '',i11)') vitd_index
        end if

!   Check for valid downwelling infrared flux information 
        if(ir_col /= no_data) then
          met(icnt,ip_ir) = in_data(ir_col)
          if(dabs(met(icnt,ip_ir)-mflag) <= eps) then
            if(icnt /= 1) met(icnt,ip_ir) = met(icnt-1,ip_ir)
            if(icnt == 1) met(icnt,ip_ir) = mflag
          end if
          if(icnt == 1.and.dabs(met(icnt,ip_ir)-mflag) <= eps)          &
            write(10,'(''2     model calculated downwelling IR, '',     &
                  &''met_id = '',i26)') vitd_index
          if(dabs(met(icnt,ip_ir)-mflag) > eps) then
            if((met(icnt,ip_ir) < 0d0).or.                              &
                                       (met(icnt,ip_ir) > 1355d0)) then
              write(10,'('' Downwelling IR radiation out of range, '',  &
                    &''met_id = '',i25)') vitd_index
              write(10,'(''   0 < Downwelling IR < 1355 W/m^2, your '', &
                    &''value is:              '',f10.3)')               &
                                                       met(icnt,ip_ir)
              write(10,'(''   On day '',i5,'' hour '',i3)')             &
                    int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
              write(10,'(''   model calculated downwelling IR, '',      &
                    &''met_id = '',i29)') vitd_index
              met(icnt,ip_ir) = mflag
            end if
          end if
        else
          met(icnt,ip_ir) = mflag
          if(icnt == 1) write(10,'(''2     model calculated downwell'', &
            &''ing IR, met_id = '',i26)') vitd_index
        end if

!   Check for valid upwelling infrared flux information 
        if(irup_col /= no_data) then             
          met(icnt,ip_irup) = in_data(irup_col)
          if(dabs(met(icnt,ip_irup)-mflag) <= eps) then
            if(icnt /= 1) met(icnt,ip_irup) = met(icnt-1,ip_irup)
            if(icnt == 1) met(icnt,ip_irup) = mflag
          end if
          if(icnt == 1.and.dabs(met(icnt,ip_irup)-mflag) <= eps)        &
            write(10,'(''2     FASST will calculate upwelling IR, '',   &
                  &''met_id = '',i24)') vitd_index
          if(dabs(met(icnt,ip_irup)-mflag) > eps) then
            if((met(icnt,ip_irup) < 0d0).or.                            &
                                     (met(icnt,ip_irup) > 1355d0)) then
              write(10,'('' Upwelling IR radiation out of range, '',    &
                    &''met_id = '',i27)') vitd_index
              write(10,'(''   0 < upwelling IR < 1355 W/m^2, your '',   &
                    &''value is:                '',f10.3)')             &
                                                      met(icnt,ip_irup)
              write(10,'(''   On day '',i5,'' hour '',i3)')             &
                    int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
              write(10,'(''   FASST will calculate upwelling IR, '',    &
                    &''met_id = '',i27)') vitd_index
              met(icnt,ip_irup) = mflag
            end if
          end if
        else
          met(icnt,ip_irup) = mflag
          if(icnt == 1) write(10,'(''2     FASST will calculate upwel'',&
            &''ling IR, met_id = '',i24)') vitd_index
        end if

!   Check for valid solar zenith angle information 
        if(zen_col /= no_data) then
          met(icnt,ip_zen) = in_data(zen_col)
          if(icnt == 1.and.dabs(met(icnt,ip_zen)-mflag) <= eps)         &
            write(10,'(''2     model calculated solar zenith angle, '', &
                  &''met_id = '',i22)') vitd_index
        else
           met(icnt,ip_zen) = mflag
           if(icnt == 1) write(10,'(''2     model calculated solar '',  &
             &''zenith angle, met_id = '',i22)') vitd_index
        end if

        if(dabs(met(icnt,ip_zen)-mflag) > eps)then
          if(met(icnt,ip_zen) < 0d0.or.met(icnt,ip_zen) > 180d0) then
            write(10,'('' Solar Zenith angle out of range, vitd_'',     &
                  &''index = '',i27)') vitd_index
            write(10,'(''   0 < solar zenith < 180, your value is:   '',&
                  &''                    '',f10.4)') met(icnt,ip_zen)
            write(10,'(''   On day '',i5,'' hour '',i3)')               &
                  int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
            write(10,'(''   model calculated solar zenith angle, '',    &
                  &''met_id = '',i25)') vitd_index
            met(icnt,ip_zen) = mflag
          end if
        end if

!   Check for valid solar azimuth angle information 
        if(az_col /= no_data) then
          met(icnt,ip_az) = in_data(az_col)
          if(icnt == 1.and.dabs(met(icnt,ip_az)-mflag) <= eps)          &
            write(10,'(''2     model calculated solar azimuth angle, '',&
                  &''met_id = '',i21)') vitd_index
        else
          met(icnt,ip_az) = mflag
          if(icnt == 1) write(10,'(''2     model calculated solar '',   &
            &''azimuth angle, met_id = '',i21)') vitd_index
        end if

        if(dabs(met(icnt,ip_az)-mflag) > eps) then
          if(met(icnt,ip_az) < 0d0.or.met(icnt,ip_az) > 360d0) then
            write(10,'('' Solar Azimuth angle out of range, vitd_'',    &
                  &''index = '',i26)') vitd_index
            write(10,'(''   0 < solar azimuth < 360, your value is: '', &
                  &''                     '',f10.4)') met(icnt,ip_az)
            write(10,'(''   On day '',i5,'' hour '',i3)')               &
                  int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
            write(10,'(''   model calculated solar azimuth angle, '',   &
                  &''met_id = '',i24)') vitd_index
            met(icnt,ip_az) = mflag
          end if
        end if

!   Check for valid soil surface temperature(K) information 
        if(tsoil_col /= no_data) then
          met(icnt,ip_tsoil) = in_data(tsoil_col)
          if(stflag == 'tc'.and.dabs(met(icnt,ip_tsoil)-mflag) > eps)   &
             met(icnt,ip_tsoil) = met(icnt,ip_tsoil) + 273.15d0
          if(stflag == 'tf'.and.dabs(met(icnt,ip_tsoil)-mflag) > eps)   &
             met(icnt,ip_tsoil) = (met(icnt,ip_tsoil) - 32d0)*5d0/9d0   & 
                                   + 273.15
          if(dabs(met(icnt,ip_tsoil)-mflag) <= eps) then
            if(icnt == 1) then
              met(icnt,ip_tsoil) = mflag
              write(10,'(''2     FASST will calculate soil surface '',  &
                    &''temperature, met_id = '',i12)') vitd_index
            else
              met(icnt,ip_tsoil) = met(icnt-1,ip_tsoil)
            end if
          end if

          if(dabs(met(icnt,ip_tsoil)-mflag) > eps) then
            if(met(icnt,ip_tsoil) < 173.15d0.or.                        &
                                    met(icnt,ip_tsoil) > 373.15d0) then
              write(10,'('' Surface temperautres out of range, vitd_'', &
                    &''index = '',i25)') vitd_index
              write(10,'(''   173.15 < surface temp < 373.15 K, your '',&
                    &''value  is:            '',f10.3)')                &
                                                     met(icnt,ip_tsoil)
              write(10,'(''   On day '',i5,'' hour '',i3)')             &
                   int(met(icnt,ip_doy)),int(met(icnt,ip_hr))
              write(10,'(''   FASST will calculate soil surface '',     &
                    &''temperature, met_id = '',i15)') vitd_index
              met(icnt,ip_tsoil) = mflag
            end if
          end if
        else
          met(icnt,ip_tsoil) = mflag
          if(icnt == 1) write(10,'(''2     FASST will calculate soil '',&
            &''surface temperature, met_id = '',i12)') vitd_index
        end if

! density material (snow)
        if(denmat_col /= no_data) then
          met(icnt,ip_dom) = in_data(denmat_col)
          if(dabs(met(icnt,ip_dom)-mflag) > eps) then
            if(pdflag == 'k') met(icnt,ip_dom) = met(icnt,ip_dom)*1d-3
            if(pdflag == 'p')met(icnt,ip_dom) =met(icnt,ip_dom)*1.602d-2
          end if
        else
          met(icnt,ip_dom) = mflag
        end if

! Check for valid snow depth(m) information
! Note: want actual snow depth, not swe
        if(sndepth_col /= no_data) then
          met(icnt,ip_sd) = in_data(sndepth_col)
          if(dabs(met(icnt,ip_sd)-mflag) > eps) then
            if(sdflag2 == 'i')met(icnt,ip_sd) = met(icnt,ip_sd)*3.048d-1
            if(sdflag2 == 'f') met(icnt,ip_sd) = met(icnt,ip_sd)*2.54d-2
            if(sdflag2 == 'c') met(icnt,ip_sd) = met(icnt,ip_sd)*1d-2
            if(sdflag2 == 'm') met(icnt,ip_sd) = met(icnt,ip_sd)*1d-3
            if(sdflag1 == '2') then
              if(dabs(met(icnt,ip_dom)-mflag) > eps) then
                met(icnt,ip_sd) = met(icnt,ip_sd)/met(icnt,ip_dom)
              else
                met(icnt,ip_sd) = met(icnt,ip_sd)*3d0
              end if
            end if
          else
            if(icnt /= 1) met(icnt,ip_sd) = met(icnt-1,ip_sd)
            if(icnt == 1) then
              met(icnt,ip_sd) = mflag
              write(10,'(''2     FASST will calculate snow depth, '',   &
                &''met_id = '',i26)') vitd_index
            end if
          end if

          if(icnt /= 1.and.(met(icnt,ip_sd) > met(icnt-1,ip_sd).and.    &
                                       dabs(met(icnt,ip_prec)) <= eps)) & 
            met(icnt,ip_sd) = met(icnt-1,ip_sd)
          if(met(icnt,ip_sd) < 0d0) met(icnt,ip_sd) = 0d0
        else
          met(icnt,ip_sd) = mflag
          if(icnt == 1) write(10,'(''2     FASST will calculate snow '',&
            &''depth, met_id = '',i26)') vitd_index
       end if

! Visibility distance (km)
        if(visdis_col /= no_data) then
          met(icnt,ip_vis) = in_data(visdis_col)
          if(dabs(met(icnt,ip_vis)-mflag) > eps.and.vdflag == 'm')      &
            met(icnt,ip_vis) = met(icnt,ip_vis)*1.6099344d0
        else
          met(icnt,ip_vis) = mflag
        end if

! aersol type
        if(vistype_col /= no_data) then
          met(icnt,ip_aer) = in_data(vistype_col)
          if(dabs(aint(met(icnt,ip_aer))-aint(mflag)) <= eps)           &
            met(icnt,ip_aer) = mflag
        else
          met(icnt,ip_aer) = mflag
        end if

! precipitation diameter (mm)
        if(precdiam_col /= no_data) then
          met(icnt,ip_pd) = in_data(precdiam_col)
          if(dabs(met(icnt,ip_pd)-mflag) > eps) then
            if(pdflag == 'm') met(icnt,ip_pd) = met(icnt,ip_pd)*1d3
            if(pdflag == 'c') met(icnt,ip_pd) = met(icnt,ip_pd)*1d1
            if(pdflag == 'f') met(icnt,ip_pd) = met(icnt,ip_pd)*3.048d2
            if(pdflag == 'i') met(icnt,ip_pd) = met(icnt,ip_pd)*2.54d1
          end if
        else
          if(dabs(met(icnt,ip_prec)) <= eps.or.                         &
                                  dabs(met(icnt,ip_prec2)) <= eps) then 
            met(icnt,ip_pd) = 0d0
          else
            met(icnt,ip_pd) = mflag
          end if
       end if

        do i=1,mxcol
          tmp(i) = in_data(i)
        end do

        icnt = icnt + 1 
      end do

      end subroutine read_met

! ******************************************************************************
      real(kind=8) function rh(atemp)

      implicit none

      real(kind=8),intent(in):: atemp

! local variables
      real(kind=8):: num,denom


      num = 22.452d0*atemp
      denom = 272.55d0 + atemp

      rh = dexp(num/denom)

      end 

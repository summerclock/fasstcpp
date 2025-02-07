      program met_reader

      use met_global
      use module_radiation

! This program translates any type of met data into the form needed by FASST.
! If no data exists it either calculates it or makes assumptions as to the values.
! This program is meant to be run prior to FASST.

! calls the following subroutines:
!     datatrans
!     read_met
!     sol_zen
!     Solflx
!     emisatm
!     dnirflx
!     met_out


      implicit none

      integer(kind=4):: io,i,j,k,iw,hdrlines,ntlines,itest,ltest,err
      integer(kind=4):: icld(3),shdrlines,slines,nmlen
      real(kind=8):: cover(3),hgt(3)
      real(kind=8):: mcd,hcd,ematm,iir,prcp,sum_tsol,sum_dir,rst_cutoff
      real(kind=8):: sum_dif,sum_ir,hr,saz,szen,stot,sdir,sdif
      real(kind=8),allocatable:: met(:,:)
      character(len=200):: hdrline
      character(len=500):: a
      character(len=512):: ctemp
      character(len=512):: weather_input,met_file_name,output1
      character(len=512):: initmet_file_name, meta_file_name
!      character(len=:),allocatable:: weather_input,met_file_name,output1
!      character(len=:),allocatable:: initmet_file_name, meta_file_name


! initialize variables
      hdrlines = 0
      ntlines = 0
      itest = 0
      ltest = 0
      err = 0
      shdrlines = 0
      slines = 0
      nmlen = 0

      do i=1,3
        icld(i) = 0
        cover(i) = 0d0
        hgt(i) = 0d0
      end do

      mcd = 0d0
      hcd = 0d0
      ematm = 0d0
      iir = 0d0
      prcp = 0d0
      sum_tsol = 0d0
      sum_dir = 0d0
      sum_dif = 0d0
      sum_ir = 0d0
      hr = 0d0
      saz = 0d0
      szen = 0d0
      stot = 0d0
      sdir = 0d0
      sdif = 0d0
      rst_cutoff = 0d0
      ctemp = ' '

!! if command line argument exists, use its I/O filenames
      call getarg(1,weather_input)
!      call get_command_argument(1,weather_input)
      weather_input = weather_input(1:len_trim(weather_input))

!!      call getarg(1,ctemp)
!      call get_command_argument(1,ctemp)
!      nmlen = len_trim(ctemp)
!      allocate(character(len=nmlen):: weather_input,stat=err)
!      weather_input = ' '
!      weather_input = ctemp(1:nmlen)
      
      if(weather_input == ' ') weather_input = 'weather.inp'

! open and read the weather input file
      io = 0
      open(unit=15,file=weather_input,iostat=io,status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening weather input file, error='',i3)') io
        write(*,*)' '
        write(*,'(a)') weather_input
        stop
      end if

      read(15,'(a)') meta_file_name
      meta_file_name = meta_file_name(1:len_trim(meta_file_name))
!      ctemp = ' '
!      nmlen = 0
!      read(15,'(a)') ctemp
!      nmlen = len_trim(ctemp)
!      allocate(character(len=nmlen):: meta_file_name,stat=err)
!      meta_file_name = ' '
!      meta_file_name = ctemp(1:nmlen)

      read(15,'(a)') met_file_name
      met_file_name = met_file_name(1:len_trim(met_file_name))
!      ctemp = ' '
!      nmlen = 0
!      read(15,'(a)') ctemp
!      nmlen = len_trim(ctemp)
!      allocate(character(len=nmlen):: met_file_name,stat=err)
!      met_file_name = ' '
!      met_file_name = ctemp(1:nmlen)

      read(15,'(a)') output1
      output1 = output1(1:len_trim(output1))
!      ctemp = ' '
!      nmlen = 0
!      read(15,'(a)') ctemp
!      nmlen = len_trim(ctemp)
!      allocate(character(len=nmlen):: output1,stat=err)
!      output1 = ' '
!      output1 = ctemp(1:nmlen)

      read(15,*) infer_test
      if(infer_test == 1) read(15,'(a)') initmet_file_name   !original met output file
!      if(infer_test == 1) then
!        ctemp = ' '
!        nmlen = 0
!        read(15,'(a)') ctemp
!        nmlen = len_trim(ctemp)
!        allocate(character(len=nmlen):: initmet_file_name,stat=err)
!        initmet_file_name = ' '
!        initmet_file_name = ctemp(1:nmlen)
!      end if
      read(15,*) albedo_fasst,rst_cutoff

      close(15)

! open and read the met translator file; get the number of headerlines in the met file
      io = 0
      open(unit=5,file=meta_file_name,iostat=io,status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening meta_file_name, error='',i3)') io
        write(*,*)' '
        write(*,'(a)') meta_file_name
        stop
      end if

      call datatrans(hdrlines,shdrlines)

      close(5)

! open the met file
      io = 0
      open(unit=1,file=met_file_name,iostat=io,status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening met input file, error='',i3)') io
        write(*,*)' '
        write(*,'(a)') met_file_name
        stop
      end if

! open the inferencing flag output file
      io = 0
      open(unit=10,file='met_inferred.out',iostat=io,status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening the log file met_inferred.out,'',  &
              &'' error='',i3)') io
        write(*,*)' '
        write(*,'(a)') 'met_inferred.out'
        stop
      end if

! Open the main output file
      io = 0
      open(unit=2,file=output1,iostat=io,status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening output1, error='',i3)') io
        write(*,*)' '
        write(*,'(a)') output1
        stop
      end if

! open the old met file if it exists
      if(infer_test == 1) then
        io = 0
        open(unit=3,file=initmet_file_name,iostat=io,status='unknown')
        if(io /= 0) then
          write(*,'('' Error opening initmet_file, error='',i3)') io
          write(*,*)' '
          write(*,'(a)') initmet_file_name
          stop
        end if
      end if

! check that number of data lines is correct
      if(hdrlines /= 0.and.hdrlines-shdrlines /= 0) then
        do j=1,hdrlines-shdrlines
          read(1,'(a)') hdrline
        end do
      end if

      ntlines = 0
      io = 0
      if(shdrlines == 0) then
        do while(io /= -1)
          read(1,*,iostat=io)

          ntlines = ntlines + 1
        end do
        if(ntlines > 1) ntlines = ntlines - 1
      else
        do j=1,shdrlines
          read(1,'(a)',iostat=io) hdrline
        end do

        ltest = 0
        do while(io /= -1.and.ltest == 0)
          read(1,'(a)',iostat=io) a
          i = 1
          itest = 0
          do while(i <= len_trim(a).and.itest == 0)
            if(ichar(a(i:i)) /= 32) then
              if(ichar(a(i:i)) >= 45.and.ichar(a(i:i)) <= 57) then
                ntlines = ntlines + 1
                itest = 1
              else
                exit
              end if
            end if
            i = i + 1
          end do
          if(itest == 0) ltest = 1
        end do
      end if

      nlines = min(nlines,ntlines)

      rewind(1)

! size the met array
      allocate(met(nlines,maxcol),stat=err)                              !size met array

! read any header information except for the column headings and units.
      if(hdrlines /= 0) then
        do j=1,hdrlines-shdrlines
          read(1,'(a)') hdrline
        end do
      end if

! beginning of loop
      do j=1,ndatpos

! The met data can have up to 45 parameters-zero out the met array
        do i=1,nlines
          do k=1,maxcol
            met(i,k) = 0d0
            met(i,k) = mflag
          end do
        end do

        call read_met(shdrlines,rst_cutoff,met)

        slines = 1
        if(ttflag == '2'.and.ave_period >= 16) then
          if(ave_period == 16) then
            slines = 4
          else if(ave_period == 17) then
            slines = 6
          else if(ave_period == 18) then
            slines = 8
          else if(ave_period == 19) then
            slines = 12
          else if(ave_period == 20) then
            slines = 24
          end if
        end if

        io = 0
        iw = 1
        do while(iw <= nlines.and.io /= -1)
          sum_tsol = 0d0
          sum_dir = 0d0
          sum_dif = 0d0
          sum_ir = 0d0

          do k=1,slines
            hr = met(iw,ip_hr) + k - 1
            szen = met(iw,ip_zen)
            saz = met(iw,ip_az)
            stot = met(iw,ip_tsol)
            sdir = met(iw,ip_dir)
            sdif = met(iw,ip_dif)
            iir = met(iw,ip_ir)

! calculate the solar zenith angle if missing
            if(dabs(szen-mflag) <= eps.or.dabs(saz-mflag) <= eps) then
              call sol_zen(met(iw,ip_year),met(iw,ip_doy),hr,           &
                           met(iw,ip_min),szen,saz)
            end if

! calculate missing total total, direct and diffuse solar; check non-missing values
            if(szen >= 90d0) then                                        !night
              stot = 0d0
              sdir = 0d0
              sdif = 0d0
            else                                                         !day
              if(dabs(stot-mflag) > eps.and.                            &
                (dabs(sdir-mflag) > eps.and.                            &
                                          dabs(sdif-mflag) > eps)) then
                if(dabs(stot-(sdir+sdif)) > 1d-3) then
                  sdif = stot - sdir
                  if(sdif < eps) then
                    sdif = 0d0
                    sdir = stot
                  end if
                end if
              else if(dabs(stot-mflag) > eps.and.                       &
                     (dabs(sdir-mflag) > eps.and.                       &
                                         dabs(sdif-mflag) <= eps)) then
                sdif = stot - sdir
                  if(sdif < eps) then
                    sdif = 0d0
                    sdir = stot
                  end if
              else if(dabs(stot-mflag) > eps.and.                       &
                     (dabs(sdir-mflag) <= eps.and.                      &
                                         dabs(sdif-mflag) > eps)) then
                sdir = stot - sdif
                  if(sdir < eps) then
                    sdif = stot
                    sdir = 0d0
                  end if
              else if(dabs(stot-mflag) <= eps.and.                      &
                     (dabs(sdir-mflag) > eps.and.                       &
                                           dabs(sdif-mflag) > eps)) then
                stot = sdir + sdif
              else

                cover(1) = met(iw,ip_lcd)
                cover(2) = met(iw,ip_mcd)
                cover(3) = met(iw,ip_hcd)
                hgt(1) = met(iw,ip_lhgt)
                hgt(2) = met(iw,ip_mhgt)
                hgt(3) = met(iw,ip_hhgt)
                icld(1) = int(met(iw,ip_lct))
                icld(2) = int(met(iw,ip_mct))
                icld(3) = int(met(iw,ip_hct))

                prcp = met(iw,ip_prec) + met(iw,ip_prec2)

                call Solflx(icld,szen,met(iw,ip_doy),prcp,stot,cover,   &
                            hgt,sdir,sdif)

                met(iw,ip_lcd) = cover(1)
                met(iw,ip_mcd) = cover(2)
                met(iw,ip_hcd) = cover(3)
                met(iw,ip_lhgt) = hgt(1)
                met(iw,ip_mhgt) = hgt(2)
                met(iw,ip_hhgt) = hgt(3)
                met(iw,ip_lct) = float(icld(1))
                met(iw,ip_mct) = float(icld(2))
                met(iw,ip_hct) = float(icld(3))
              end if

              sum_tsol = sum_tsol + stot
              sum_dir = sum_dir + sdir
              sum_dif = sum_dif + sdif
            end if

! correct cloud parameters if necessary
            if(dabs(met(iw,ip_lcd)) <= eps) then
              met(iw,ip_lhgt) = 0d0
              met(iw,ip_lct) = 0
            end if
            if(dabs(met(iw,ip_mcd)) <= eps) then
              met(iw,ip_mhgt) = 0d0
              met(iw,ip_mct) = 0
            end if
            if(dabs(met(iw,ip_hcd)) <= eps) then
              met(iw,ip_hhgt) = 0d0
              met(iw,ip_hct) = 0
            end if

!           If we are missing the downwelling irflux calculate a value, but first we must
!           calculate an effective atmospheric emissivity
            if(dabs(iir-mflag) <= eps) then
              mcd = met(iw,ip_mcd)
              hcd = met(iw,ip_hcd)
              call emisatm(met(iw,ip_tmp),met(iw,ip_rh),ematm)

              call dnirflx(ematm,met(iw,ip_tmp),met(iw,ip_doy),         &
                           met(iw,ip_lhgt),met(iw,ip_mhgt),             &
                           met(iw,ip_hhgt),met(iw,ip_lcd),              &
                           mcd,hcd,iir)
            end if
            sum_ir = sum_ir + iir

            if(slines == 1) then
              met(iw,ip_zen) = szen
              met(iw,ip_az) = saz
            end if
          end do   !k=1,slines

          met(iw,ip_tsol) = sum_tsol/float(slines)
          met(iw,ip_dir) = sum_dir/float(slines)
          met(iw,ip_dif) = sum_dif/float(slines)
          met(iw,ip_ir) = sum_ir/float(slines)

!         snow depth (m)  Note: want actual snow depth, not swe
!          met(iw,ip_sd) = mflag    !!!!!!! for wrf data only
          iw = iw + 1
        end do   !while  iw=1,nlines

!       write to the output file
        call met_out(iw-1,met)
      end do

! close all opened files
      close(1)
      close(2)
      if(infer_test == 1) close(3)
      close(10)

      deallocate(met)

      stop
      end program met_reader

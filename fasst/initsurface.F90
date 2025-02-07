      subroutine initsurface(nt0,nm0,sprint,sname,ttest,mtest,          &
                             nstr_flag,code1,phie,pdens,oldsd,oldhi,    &
                             sdensi)

      use fasst_global
      use module_canopy
      use module_lowveg

! calls the following subroutines:
!     upr_case (appended to US_soil_tools.f)
!     read_old_data (only in multi-run scenarios)
!     missing_met
!     veg_propl
!     veg_proph
!     get_user_soil_params
!     initprofile

! uses the function: met_date

      implicit none

      integer(kind=4),intent(in):: nt0,nm0,sprint
      character(len=4),intent(in):: sname(maxl)
      character(len=1),intent(in):: ttest,mtest
      integer(kind=4),intent(out):: nstr_flag(maxn)
      integer(kind=4),intent(inout):: code1
      real(kind=8),intent(inout):: phie,pdens,oldsd,oldhi,sdensi


! local variables
      integer(kind=4):: i,j,str_flag(maxl),d1i
      integer(kind=4):: wcheck,wstart,wend,tstps
      real(kind=8):: met1(1,maxcol),rtstps,modyo,daylim
      real(kind=8):: time1,time2,timeo(moverlap,4),met_date,mody


! STEP 1
! initialize certain check/test variables
      iseason = 0
      code1 = 0
      tstps = 0

      wcheck = 0
      wstart = 0
      wend = 0
      d1i = 0

      do i=1,maxl
        str_flag(i) = 0
      end do

      do i=1,maxn
        nstr_flag(i) = 0
      end do

      oldsd = 0d0
      oldhi = 0d0
      sdensi = 0d0
      zh = 0d0
      pdens = 0d0
      rtstps = 0d0
      mody = 0d0
      modyo = 0d0
      daylim = 0d0
     
      do i=1,nclayers
        laif(i) = 0d0
        dzveg(i) = 0d0
      end do

      iheightn = iheight
      rtstps = 24d0/timstep
      tstps = int(rtstps)


! STEP 2
! initialize surface conditions and soil profile, max/min values, etc
      if(infer_test == 1) then
        do j=1,maxcol
          met1(1,j) = mflag                                              !oldpos = 1
        end do

        d1i = 1
        wstart = max(moverlap,istart)
        call read_old_data(d1i,mpos,wstart,sdensi,phie,timeo)

        do i=1,maxcol
          met1(1,i) = met(istart-1,i)
        end do

        d1i = 1
        call missing_met(d1i,oldpos,oldpos,met1)

! double check met parameters for missing values
        do j=istart,iend
          do i=1,maxcol-5
            if(i /= 26.and.i /= 28) then
              if(aint(dabs(met(j,i)-mflag)*1d5)*1d-5 <= eps.or.         &
                 aint(dabs(met(j,i)-aint(mflag))*1d5)*1d-5 <= eps) then
                wcheck = 1
                if(wstart == 0) wstart = j
                wend = j
              end if
            end if
          end do
        end do

        if(wcheck == 1) then
          do j=1,maxcol
            met1(1,j) = mflag
          end do

          d1i = 2
          call missing_met(d1i,wstart,wend,met1)
        end if

        sdens(oldpos) = sdensi

        ft(oldpos) = ftemp
        tt(oldpos) = toptemp
        airt(1,oldpos) = met(oldpos,ip_tmp) + Tref
        airt(2,oldpos) = dmet1(oldpos,4) + Tref

        if(vegl_type /= 0) call veg_propl(biome_source,new_vtl,         &
                                          vegl_type)                     !get constant low veg params
        if(vegh_type /= 0) call veg_proph(biome_source,new_vth,         &
                                          vegh_type)                     !get constant high veg params
      else if(infer_test == 0) then
        do i=1,nlayers
          do j=1,maxp
            isoilp(i,j) = spflag
            soilp(i,j) = spflag
          end do
        end do

        do i=1,nlayers
          str_flag(i) = 0
          if(soiltype(i) == -1) then                                     !user supplied soil parameters
            str_flag(i) = 3
            d1i = len_trim(sname(i))
            call upr_case(d1i,sname(i))

            call get_user_soil_params(sname(i),sclass(i),soiltype(i),   &
                                      isoilp(i,1),isoilp(i,2),          &
                                      isoilp(i,3),isoilp(i,4),          &
                                      isoilp(i,5),isoilp(i,6),          &
                                      isoilp(i,7),isoilp(i,8),          &
                                      isoilp(i,9),isoilp(i,10),         &
                                      isoilp(i,11),isoilp(i,13),        &
                                      isoilp(i,12),isoilp(i,14),        &
                                      isoilp(i,18),isoilp(i,19),        &
                                      isoilp(i,20),isoilp(i,21),        &
                                      isoilp(i,22),isoilp(i,23))

            rewind(90)
          else 
            if(soiltype(i) == 0) then                                    !unknown soil type
              soiltype(i) = 7                                            !default value is dirty sand
              str_flag(i) = 1
            else if(soiltype(i) == -2) then                              !unknown fine grained
              soiltype(i) = 12
              str_flag(i) = 2
            else if(soiltype(i) == -3) then                              !unknown coarse grained
              soiltype(i) = 6
              str_flag(i) = 1
            else if(soiltype(i) == -4) then                              !unknown disturbed
              soiltype(i) = 7
              str_flag(i) = 0
              rho_fac(i) = 8d-1
            end if
            code1 = 1
            error_code = 1
            if(single_multi_flag == 0) then
              write(10,'('' Used default soil type. '',/)')
            end if
          end if
        end do

! double check met parameters for missing values
        do j=istart,iend
          do i=1,maxcol-5
            if(i > 4.and.(i /= 26.and.i /= 28)) then
              if(aint(dabs(met(j,i)-mflag)*1d5)*1d-5 <= eps.or.         &
                 aint(dabs(met(j,i)-aint(mflag))*1d5)*1d-5 <= eps) then
                wcheck = 1
                if(wstart == 0) wstart = j
                wend = j
              end if
            end if
          end do
        end do
        if(wcheck == 1) then
          do j=1,maxcol
            met1(1,j) = mflag
          end do

          d1i = 2
          call missing_met(d1i,wstart,wend,met1)
        end if

! initialize snow depth
        if(met(istart,ip_sd) > eps.and.                                 &
               aint(dabs(met(istart,ip_sd)-mflag)*1d5)*1d-5 > eps) then
          if(dabs(hsaccum) <= eps.or.                                   &
                       aint(dabs(hsaccum-spflag)*1d5)*1d-5 <= eps) then 
            hsaccum = met(istart,ip_sd)
          else if(dabs(hsaccum) >= eps.and.                             &
                        aint(dabs(hsaccum-spflag)*1d5)*1d-5 > eps) then
            hsaccum = dmax1(hsaccum,met(istart,ip_sd))
          end if
        end if
        hsaccum = anint(hsaccum*1d15)*1d-15
        hm = anint((hsaccum + newsd + hi)*1d15)*1d-15                    !m

! initialize soil profile
        call initprofile(ttest,nt0,mtest,nm0,str_flag,nstr_flag)

        if(veg_flagl == 1) call veg_propl(biome_source,new_vtl,         &
                                          vegl_type)                     !get constant low veg params

        if(veg_flagh == 1) call veg_proph(biome_source,new_vth,         &
                                          vegh_type)                     !get constant high veg params

        if(hsaccum > eps.or.hi > eps) toptemp = Tref
        sdensi = sdensw
        storll = 0d0
        storls = 0d0
      end if   !if(infer_test == 1) then

! initialize soil roughness length
      if(rough <= eps) then
        if(ntype(nnodes) >= 19) then                                     !concrete ,asphalt, bedrock, glaciers
          rough = 1d-3
        else if(ntype(nnodes) <= 4.or.ntype(nnodes) == 15) then          !gravels, peat
          rough = 5d-2
        else if(ntype(nnodes) > 4.and.ntype(nnodes) <= 6) then           !rough sands
          rough = 1d-2
        else                                                             !all else
          rough = 1d-3
        end if
      end if
      rough = anint(rough*1d15)*1d-15

! initialize snow/ice
      oldsd = hsaccum
      oldhi = hi
      if(dabs(pdens) <= eps) pdens = 5d2
      if(oldsd+oldhi > eps) phie = (1d0 - pdens*1d-3)*0.95d0

! write the header lines for the node input file
      if(single_multi_flag == 0) then
        write(3,'('' Total Number of Nodes: '',i4,'',  freq_id:'',i10)')& 
              &nnodes,freq_id
        write(3,'('' Year   JD   Hr    M    N  USCS    Depth   Grtemp   &  
              &Water         Ice       W + I       Vapor     F/T'')')
        write(3,'(''                                     m       K'')')
      end if

! write the snow output file header line
      if(sprint == 1) then
        write(55,1000) freq_id,iend,9,lat,mlong,elev,vitd_index,        &
                         met_count,timeoffset,timstep,mflag
 1000   format(i10,1x,i6,1x,i3,1x,f10.6,1x,f11.6,1x,f11.6,1x,i8,1x,i8,  &
               1x,f6.2,1x,f5.2,1x,f8.2)
        write(55,92)
  92    format(6x,'year',8x,'doy',8x,'sdold',10x,'add',12x,'meta',11x,  &
               'dwind',10x,'atop',11x,'abot',11x,'sd',11x,'delta sd',   &
               3x,'mode',7x,'diameter',12x,'swe',11x,'%wet')
      end if

! check met file for missing time steps
      do i=istart,iend
        if(i /= 1) then
          time1 = met_date(met(i-1,1),met(i-1,2),met(i-1,3),met(i-1,4))
          time2 = met_date(met(i,1),met(i,2),met(i,3),met(i,4))

          mody = met(i,1) - aint(met(i,1)*2.5d-1)*4d0
          modyo = met(i-1,1) - aint(met(i-1,1)*2.5d-1)*4d0
          daylim = 1.042d0  !1.003d0
!          if(mody /= modyo) daylim = 1.07d0
          if(met(i-1,1) /= met(i,1)) daylim = 25.07d0

          if(dabs(mody) <= eps) then
            if((time2 - time1)*367d0*tstps >= daylim) then
              do j=1,maxcol
                met1(1,j) = met(i-1,j)                                   !oldpos = 1
              end do

              d1i = 3
              j = i - 1
              call missing_met(d1i,j,i,met1)
            end if
          else
            if((time2 - time1)*366d0*tstps >= daylim) then
              do j=1,maxcol
                met1(1,j) = met(i-1,j)                                   !oldpos = 1
              end do

              d1i = 3
              j = i - 1
              call missing_met(d1i,j,i,met1)
            end if
          end if
        end if
      end do

      end subroutine initsurface

      subroutine initwater(sprint,phie,pdens,oldsd,oldhi,sdensi)

      use fasst_global
      use module_canopy
      use module_lowveg

! calls the following subroutines:
!     read_old_data (only in multi-run scenarios)
!     missing_met
!     veg_propl
!     veg_proph

! uses the function: met_date

      implicit none

      integer(kind=4),intent(in):: sprint
      real(kind=8),intent(inout):: phie,pdens,oldsd,oldhi,sdensi

! local variables
      integer(kind=4):: i,j,d1i,wcheck,wstart,wend
      real(kind=8):: met1(1,maxcol)
      real(kind=8):: time1,time2,timeo(moverlap,4),met_date


! STEP 1
! initialize certain check/test variables
      iseason = 0

      wcheck = 0
      wstart = 0
      wend = 0
      d1i = 0

      oldsd = 0d0
      oldhi = 0d0
      sdensi = 0d0
      zh = 0d0
      pdens = 0d0
     
      do i=1,nclayers
        laif(i) = 0d0
        dzveg(i) = 0d0
      end do

      iheightn = iheight


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

        sdens(oldpos) = sdensi

        ft(oldpos) = ftemp
        tt(oldpos) = toptemp
        airt(1,oldpos) = met(oldpos,ip_tmp) + Tref
        airt(2,oldpos) = dmet1(oldpos,4) + Tref

!        if(vegl_type /= 0) call veg_propl(biome_source,new_vtl,         &
!                                          vegl_type)                     !get constant low veg params
!        if(vegh_type /= 0) call veg_proph(biome_source,new_vth,         &
!                                          vegh_type)                     !get constant high veg params
      else if(infer_test == 0) then

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
        hsaccum = anint(hsaccum*1d5)*1d-5
        hm = anint((hsaccum + newsd + hi)*1d5)*1d-5                      !m

! initialize soil profile
!        if(veg_flagl == 1) call veg_propl(biome_source,new_vtl,         &
!                                          vegl_type)                     !get constant low veg params

!        if(veg_flagh == 1) call veg_proph(biome_source,new_vth,         &
!                                          vegh_type)                     !get constant high veg params

        if(hsaccum > eps.or.hi > eps) toptemp = Tref
        sdensi = sdensw
        storll = 0d0
        storls = 0d0
      end if   !if(infer_test == 1) then

! more snow depth initialization
      oldsd = hsaccum
      oldhi = hi
      if(dabs(pdens) <= eps) pdens = 5d2
      if(oldsd+oldhi > eps) phie = (1d0 - pdens*1d-3)*0.95d0

! write the snow output file header line
      if(sprint == 1) then
        write(55,1000) freq_id,iend,9,lat,mlong,elev,vitd_index,        &
                         met_count,timeoffset,timstep,mflag
 1000   format(i10,1x,i6,1x,i3,1x,f10.6,1x,f11.6,1x,f11.6,1x,i8,1x,i8,  &
               1x,f6.2,1x,f5.2,1x,f8.2)
        write(55,92)
  92    format(6x,'doy',8x,'sdold',10x,'add',12x,'meta',11x,'dwind',    &
               10x,'atop',11x,'abot',11x,'sd',13x,'delta sd',1x,'mode', &
               7x,'diameter')
      end if

! check met file for missing time steps
      do i=istart,iend
        if(i /= 1) then
          time1 = met_date(met(i-1,1),met(i-1,2),met(i-1,3),met(i-1,4))
          time2 = met_date(met(i,1),met(i,2),met(i,3),met(i,4))
          if(dabs(timstep-(time2-time1)*365d0*24d0) > 1d-1) then
            do j=1,maxcol
              met1(1,j) = met(i-1,j)                                     !oldpos = 1
            end do

            d1i = 3
            j = i - 1
            call missing_met(d1i,j,i,met1)
          end if
        end if
      end do

      end subroutine initwater
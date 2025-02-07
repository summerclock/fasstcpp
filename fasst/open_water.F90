      subroutine open_water(wtype,sprint,wvel,wdepth,sdensi,first_loop, &
                            code2,code4,oldsd)

      use fasst_global
      use module_canopy
      use module_lowveg
      use module_snow

! This subroutine calculates the open water temperautre based soley on latitude. 
! It is based on a routine received from Guy Seelye at AER
! NOTE: This subroutine is only a stop-gap measure

! calls the following subroutines:
!     read_old_data
!     missing_met
!     find_index (attached to this subroutine)
!     snow_model

! uses the function: met_date

      implicit none

      integer(kind=4),intent(in)::wtype,sprint
      real(kind=8),intent(in):: wvel,wdepth,sdensi
      integer(kind=4),intent(inout):: first_loop,code2,code4
      real(kind=8),intent(inout):: oldsd

! saved variables
      integer(kind=4):: first_snow

      save:: first_snow

! local variables
      integer(kind=4):: i,index1,lflag,tdate,wstart,d1i,f1i 
      real(kind=8):: scorr(12),watervel,waterdep,dh,phie,c1,c2,t1
      real(kind=8):: mixrgr,rhsurf,rhoaf,cf,pheatf,taf,stempt
      real(kind=8):: rpp,mixra,mixrf,dqdtf,wetbulba,pdens,smm
      real(kind=8):: dense,dates,xtemp,timeo(moverlap,4),walbedo,wemis


      data scorr/-12.5d0, -11.475d0, -10.45d0, -6.5d0, -2.55d0,         &
                  -0.05d0,  2.45d0,    1.225d0, 0d0,   -4.225d0,        &
                  -8.45d0,-10.475d0/ 
      
      data walbedo, wemis/0.9d0, 0.98d0/

! zero-out variables
      index1 = 0
      tdate = 0
      d1i = 0
      f1i = 0
      c1 = 0d0
      c2 = 0d0
      t1 = 0d0
      watervel = 0d0
      waterdep = 0d0
      dh = 0d0
      mixrgr = 0d0
      rhsurf = 0d0
      rhoaf = 0d0
      cf = 0d0
      pheatf = 0d0
      rpp = 0d0
      mixra = 0d0
      mixrf = 0d0
      wetbulba = 0d0
      pdens = 0d0
      smm = 0d0
      dates = 0d0
      xtemp = 0d0
      taf = 0d0
      stempt = 0d0

      if(first_loop == 0) then
        first_snow = 0                                                   !first time through snow routines
        vsmelt = 0d0                                                     !snow melt (m)
        vimelt = 0d0                                                     !ice melt (m)
        first_loop = 1
      end if

      lflag = 1                                                          !northern hemisphere
      if(lat < 0d0) lflag = -1                                           !southern hemisphere

!  !!!!!!!!!!!!!!!! Main Calculation Loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tdate = int(met(iw,ip_doy))

      call find_index(tdate,int(met(iw,ip_year)),lflag,index1)

! water velocity
      if(wtype == 1.or.wtype == aint(spflag)) then                       !open water
       watervel = 0d0                                                    !m/s
       waterdep = 10d0                                                   !m
      else                                                               !rivers, canals...
        if(dabs(wvel) <= eps.or.dabs(wvel-spflag) <= eps) then           !unknown water velocity
          watervel = 1d0
        else if(dabs(wvel-1) <= eps) then                                !>0 - <=1.5 m/s
          watervel = 0.75d0
        else if(dabs(wvel-2) <= eps) then                                !>1.5 m/s
          watervel = 3d0
        end if

        if(dabs(wdepth) <= eps.or.dabs(wdepth-spflag) <= eps) then       !unknown water depth
          waterdep = 0.5d0
        else if(dabs(wdepth-1) <= eps) then                              !>0 - <=0.8 m
          waterdep = 0.4d0
        else if(dabs(wdepth-2) <= eps) then                              !>0.8 m
          waterdep = 1d0
        end if
      end if

! water surface temperature
      met(iw,ip_tsoil) = Tref  + 56.182d0 - 0.95454d0*dabs(lat)         & 
                                      + 0.0021307d0*dabs(lat)*dabs(lat) &
                                                        + scorr(index1)  !K
      tt(iw) = met(iw,ip_tsoil)

! initialize the met calc variable dmet1(iw,)
      dmet1(iw,1) = met(iw,ip_tsol)
      if(dmet1(iw,1) < 1d-2) dmet1(iw,1) = 0d0
      dmet1(iw,2) = met(iw,ip_ir)
      dmet1(iw,3) = met(iw,ip_ws)
      dmet1(iw,4) = met(iw,ip_tmp) + 273.15d0                            !K
      dmet1(iw,5) = dmax1(1d-1,met(iw,ip_rh))                            !%
      if(met(iw,ip_prec) <= eps) then
        met(iw,ip_prec) = 0d0
        met(iw,ip_pt) = 1
      end if
      dmet1(iw,6) = met(iw,ip_prec)*1d-3                                 !m/timstep
      if(met(iw,ip_prec2) <= eps) then
        met(iw,ip_prec2) = 0d0
        met(iw,ip_pt2) = 1
      end if
      dmet1(iw,7) = met(iw,ip_prec2)*1d-3                                !m/timstep
      dmet1(iw,8) = met(iw,ip_upsol)
      if(dabs(dmet1(iw,1)) <= eps) dmet1(iw,8) = 0d0
      dmet1(iw,9) = met(iw,ip_dir)
      dmet1(iw,10) = met(iw,ip_dif)
      dmet1(iw,11) = met(iw,ip_ap)
      dmet1(iw,12) = met(iw,ip_irup)
      dmet1(iw,13) = met(iw,ip_tsoil)

! calculate the canopy temperature and moisture
      if(veg_flagh == 1) then
        call canopy_met(oldsd)                                           !correct met for canopy
        if(dmet1(iw,7) > 0d0) met(iw,ip_pt2) = float(3)
      else
        surfemisc(iw) = mflag
        do i=1,nclayers
          canopy_temp(i,iw) = mflag
         end do
      end if
      if(error_type == 2) code2 = 1

      airt(1,iw) = met(iw,ip_tmp) + Tref
      airt(2,iw) = dmet1(iw,4)

! get the low vegetation properties
      if(veg_flagl == 1) then
        call low_veg_prop(dmet1(iw,4),dzveg(nclayers))
      else
        sigfl = 0d0
      end if

! determine density and amount of new snow if no existing snow on ground
      newsd = 0d0
      pdens = 0d0
      smm = dmet1(iw,6) + dmet1(iw,7)
      if((aint(met(iw,ip_pt)) == 3.or.aint(met(iw,ip_pt2)) == 3)        &
                                                   .and.smm > 0d0) then
        d1i = 4
        pdens = dense(dmet1(iw,4),dmet1(iw,3),d1i)                       !kg/m^3
        if(aint(met(iw,ip_pt)) == 3) then
          newsd = dmet1(iw,6)*(pdens*1d-3)                               !m
        else
          newsd = dmet1(iw,7)*(pdens*1d-3)
        end if
        newsd = anint(newsd*1d10)*1d-10
        dsnow = 2d-2*((3d0/(4d0*pi))*5d-5)**(1d0/3d0)
        dsnow = anint(dsnow*1d10)*1d-10
        phie = (1d0 - pdens*1d-3)*0.95d0
        phie = anint(phie*1d10)*1d-10
        if(dabs(hsaccum) <= eps.and.newsd < 2d-3) then
          dsnow = 0d0
          phie = 0d0
          newsd = 0d0
          pdens = 0d0
        end if
      end if

      hm = hsaccum + newsd + hi                                          !m

! calculate the low vegetation temperature
      if(veg_flagl == 1) then
        if(hm > eps) f1i = 1
        c2 = 1d0
        call sp_humid(f1i,dmet1(iw,11),met(iw,ip_tsoil),c2,c1,mixrgr,  &
                      t1,t1,t1,t1,t1,t1,t1,t1,t1)

        d1i = 0
        rhsurf = 1d0
        if(hm <= eps) then
          stempt = Tref
        else
          stempt = dmet1(iw,13)
        end if
        taf = (1d0 - 7d-1*sigfl)*dmet1(iw,4) + sigfl*(6d-1*ftemp        &
                                                        + 1d-1*stempt)   !foliage temp at the ground surface (K)
        call lowveg_met(d1i,d1i,int(met(iw,ip_pt)),int(met(iw,ip_pt2)), &
                        mixrgr,rhsurf,taf,dmet1,wetbulba,rhoaf,cf,      &
                        pheatf,rpp,mixra,mixrf,dqdtf)
        ft(iw) = ftemp
        surfemisf(iw) = epf
      else
        ft(iw) = mflag
        surfemisf(iw) = mflag
      end if

! compute any ice growth
      if(met(iw,ip_tsoil) <= Tref) then
        if(iw == 1.or.dabs(met(iw,ip_hi)) <= eps) then
          met(iw,ip_hi) = dsqrt(2d0*ithcond/(idens*lhfus))              &
                            *dsqrt(dabs(met(iw,ip_tmp))*timstep*3600d0)
          if(met(iw,ip_tmp) > 0d0) met(iw,ip_hi) = 0d0
        else
          dh = (1d0/(idens*lhfus))*(-(met(iw,ip_tmp))/                  &
                            (met(iw,ip_hi)/ithcond  + hsaccum/sthcond))
          met(iw,ip_hi) = met(iw-1,ip_hi) + dh*timstep*3.6d3
          met(iw,ip_hi) = dmax1(0d0,met(iw,ip_hi))
        end if
      end if
      if(aint(dabs(met(iw,ip_hi)-mflag)*1d5)*1d-5 <= eps)               &
        met(iw,ip_hi) = 0d0

! compute snow depth
      if(hsaccum > eps.or.                                              &
          (aint(met(iw,ip_pt)) == 3.or.aint(met(iw,ip_pt2)) == 3)) then
        call snow_model(sprint,water_flag,first_snow,sdensi,hsaccum,    &
                        phie)

        vsmelt = vsmelt*1d-2                                             !m
        if(error_type == 4) code4 = 1
      else
        vsmelt = 0d0

        if(sprint == 1) then
          dates = met(iw,ip_doy) + met(iw,ip_hr)/24d0                   &
                                            + met(iw,ip_min)/(24d0*60d0)
          xtemp = 0d0
          meltfl = '0'
          write(55,91) dates,xtemp,xtemp,xtemp,xtemp,xtemp,xtemp,        &
                       hsaccum,xtemp,meltfl,xtemp
 91       format(f10.2,8(f15.8),a5,f15.8)
        end if
      end if
      met(iw,ip_sd) = hsaccum
      oldsd = hsaccum

! compute the met parameters
      if(dabs(met(iw,ip_hi)) <= eps) then
        met(iw,ip_upsol) = met(iw,ip_tsol)*walbedo
        if(met(iw,ip_tsoil) > eps)                                      &
          met(iw,ip_irup) = wemis*sigma*(met(iw,ip_tsoil)**4d0)
        surfemis(iw) = wemis
      else 
        met(iw,ip_upsol) = met(iw,ip_tsol)*ialbedo
        if(met(iw,ip_tsoil) > eps)                                      &
          met(iw,ip_irup) = iemis*sigma*(met(iw,ip_tsoil)**4d0)
        surfemis(iw) = iemis
      end if

      surfci(iw) = mflag
      surfrci(iw) = mflag
      surfcbr(iw) = mflag
      surfd(iw) = mflag
      if(dabs(met(iw,ip_hi)) <= eps) then
        surfmoist(iw) = 1d0
        surfice(iw) = 0d0
        surfd(iw) = 1d3
      else
        surfmoist(iw) = 0d0
        surfice(iw) = 1d0
        surfd(iw) = idens
      end if
      surfmoistp(iw) = surfmoist(iw)
      surficep(iw) = surfice(iw)

      if(met(iw,ip_hi) > eps) tt(iw) = Tref

! write the new input file for the next model run
! this option is used for the imets weather
!        if(dstart /= 2) then
!          wstart = max(iend-moverlap,istart)
!          if(iw == iend-(dstart-1).or.iw == iend-(dstart-1)+4.or.       &
!                                                       iw == iend) then
!            d1i = 2
!           call read_old_data(d1i,iw,wstart,sdensi,phie,timeo)
!          end if
!        else if(dstart == 2) then
          wstart = max(iend-(moverlap-1),istart)
          if(iw >= wstart) then
            d1i = 2
            call read_old_data(d1i,iw,wstart,sdensi,phie,timeo)
          end if
!        end if

      end subroutine open_water

! ******************************************************************************
      subroutine find_index(tdate,year,lflag,index1)


      implicit none

      integer(kind=4),intent(in):: tdate,lflag
      integer(kind=4),intent(in):: year
      integer(kind=4),intent(out):: index1

! local variables
      integer(kind=4):: dayimnt(12),dayemnt(12),indexm(12),indexp(12)
      integer(kind=4):: i,iday,eday,tflag
      real(kind=8):: mody

      data dayimnt / 1,32,60, 91,121,152,182,213,244,274,305,335/
      data dayemnt /31,59,90,120,151,181,212,243,273,304,334,365/

      data indexm /1,2,3,4,5,6,7,8,9,10,11,12/
      data indexp /7,8,9,10,11,12,1,2,3,4,5,6/


      iday = 0
      eday = 0
      tflag = 0
      i = 1
      do while(i <= 12.and.tflag == 0)
        mody = year - aint(year*2.5d-1)*4d0
        if(dabs(mody) <= epsilon(1d0).and.i > 2) iday = 1
        if(dabs(mody) <= epsilon(1d0).and.i > 1) eday = 1

        if(tdate >= dayimnt(i)+iday.and.tdate <= dayemnt(i)+eday) then
          if(lflag < 0) index1 = indexp(i)
          if(lflag > 0) index1 = indexm(i)
          tflag = 1
        end if

        i = i + 1
      end do

      end subroutine find_index

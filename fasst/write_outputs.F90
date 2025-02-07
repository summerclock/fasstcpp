      subroutine write_outputs(flprint,keep,vegint)

      use fasst_global

!     no subroutines called

      implicit none

      integer(kind=4),intent(in):: flprint,keep
      real(kind=8),intent(in):: vegint

! local variables
      integer(kind=4):: i,j,k,slippery,fw_flag,intval,t1(maxlines,4)
      integer(kind=4):: istartp,ncol
      real(kind=8):: max_rain,ptime,ttime,swater,swaterp,f1,swater1
      real(kind=8):: ttemp,mody,rintval,t2(maxlines,15),sd,f2,f3

! for stdmod output
      integer ice_code

! zero-out parameters
      slippery = 0
      fw_flag = 0
      intval = 0
      ice_code = 0
      ncol = 0
      max_rain = 0d0
      ptime = 0d0
      ttime = 0d0
      swater = 0d0
      swaterp = 0d0
      rintval = 0d0
      sd = 0d0
      f2 = 0d0
      f3 = 0d0

      ncol = 43

      if(keep < 0) then
        istartp = istart
      else
        istartp = iend - keep
        if(istartp < 1) istartp = 1
      end if

      do i=1,maxlines
        do j=1,4
          t1(i,j) = 0
        end do
        do j=1,15
          t2(i,j) = 0d0
        end do
      end do

! ******************************************************************************
! Main FASST output file
! ******************************************************************************
      write(2,1000) freq_id,(iend - istartp + 1),ncol,lat,mlong,elev,   &
                    vitd_index,met_count,timeoffset,timstep,mflag

! write the output file column headers
      write(2,'(''Year  JD  Hr   M   APres   ATemp     RH     WnSp   '',&  
      &'' WDir    Prec   PT   Prec2   PT2  LCAmt    LCHt   LCT  '',     &
      &''MCAmt    MCHt   MCT   HCAmt   HCHt   HCT   STot    SDir    '', &
      &''SDif   UpSol       IR     IRUp     SZen     SAz    Snow  '',   &
      &''   Ice  Grtemp  Grmoistv  Grmoistm    RCI6   RCI12    CBR '',  &
      &'' Frz Dp. T/U Dp. F/T/U    Sdens   Slpry FWType   SOdens     '',&
      &''   Emis       Vis     Aer'')')    

      write(2,'(''                    mbar     K               m/s   '',&
      &''        mm/stp       mm/stp                 km           '',   &
      &''        km                   km         W/m^2   W/m^2   '',    &
      &''W/m^2   W/m^2    W/m^2    W/m^2                     m      '', &
      &'' m       K        % vol    % mass                          '', &
      &''    m      m            kg/m^3                 kg/m^3      '', &
      &''             km'')')     

! print the met data
      max_rain = 0.25d0*25.4d0*timstep*2d0                               !(mm) needed for slippery

      do i=istartp,iend
        ttime = met(i,ip_doy) + met(i,ip_hr)/24d0
        if(met(i,ip_prec) >= max_rain) ptime = ttime

        if(met(i,ip_prec) > eps.and.aint(met(i,ip_pt)) == 1) then
          met(i,ip_pt) = met(i,ip_pt2)
          met(i,ip_pt2) = 1d0
        end if

        if(aint(dabs(met(i,ip_sd)-mflag)*1d5)*1d-5 <= eps)              &
          met(i,ip_sd) = 0d0

        if(water_flag == 0) then
          if(met(i,ip_sd) > eps) then
            slippery = 3
            if(met(i,ip_hi) > eps) then 
              fw_flag = 6
              if(iw /= 1) then
                if(met(i,ip_hi) < met(i-1,ip_hi)) fw_flag = 1
              end if
            else if (i /= 1) then
              if(met(i,ip_sd) < met(i-1,ip_sd)) fw_flag = 1
            else if (slushy(i) == 1) then
              fw_flag = 4
            else
              fw_flag = 5
            end if
          else if(met(i,ip_hi) > eps.and.met(i,ip_sd) <= eps) then
            slippery = 4    !2
            fw_flag = 0
          else if(ttime-ptime <= 0.25d0) then
            slippery = 2    !1
            fw_flag = 3
          else
            slippery = 1    !0
            fw_flag = 3
          end if
        else
          slippery = 1
          fw_flag = 3
          if(met(i,ip_hi) > eps) then
            slippery = 4
            fw_flag = 0
            if(iw /= 1) then
              if(met(i,ip_hi) < met(i-1,ip_hi)) fw_flag = 1
            end if
          else if(met(i,ip_sd) > eps) then
            slippery = 3
            fw_flag = 6
            if(iw /= 1) then
              if(met(i,ip_sd) < met(i-1,ip_sd)) fw_flag = 1
            end if
          end if
        end if

        if(dabs(met(i,ip_sd)) <= eps) sdens(i) = 0.0

        swater = surfmoist(i) + surfice(i)
        swaterp = surfmoistp(i) + surficep(i)

! determine surface freeze or thaw depth and state
        if(water_flag > 1) then
          frthick(i) = met(i,ip_hi)
          twthick(i) = 0d0
          sstate(i) = 'u'
          if(frthick(i) > eps.or.met(i,ip_tsoil) <= Tref)               &
            sstate(i) = 'f'
        end if

        write(2,1001) int(met(i,ip_year)),int(met(i,ip_doy)),           &
            int(met(i,ip_hr)),int(met(i,ip_min)),met(i,ip_ap),          &
            met(i,ip_tmp)+Tref,met(i,ip_rh),met(i,ip_ws),               &
            met(i,ip_wdir),met(i,ip_prec),int(met(i,ip_pt)),            &
            met(i,ip_prec2),int(met(i,ip_pt2)),                         &
            met(i,ip_lcd),met(i,ip_lhgt),int(met(i,ip_lct)),            &
            met(i,ip_mcd),met(i,ip_mhgt),int(met(i,ip_mct)),            &
            met(i,ip_hcd),met(i,ip_hhgt),int(met(i,ip_hct)),            &
            met(i,ip_tsol),met(i,ip_dir),met(i,ip_dif),met(i,ip_upsol), &
            met(i,ip_ir),met(i,ip_irup),met(i,ip_zen),met(i,ip_az),     &
            met(i,ip_sd),met(i,ip_hi),met(i,ip_tsoil),swater,swaterp,   &
            surfci(i),surfrci(i),surfcbr(i),frthick(i),twthick(i),      &
            sstate(i),sdens(i),slippery,fw_flag,surfd(i),surfemis(i),   &
            met(i,ip_vis),int(met(i,ip_aer))  !,ponding(i),tot_moist(i)
 
      end do

! ******************************************************************************
! vegetation and surface temperature output file
! ******************************************************************************
      if(vegint > eps) then
        rintval = timstep*60d0/vegint
        intval = int(rintval)
        if(intval < 1) intval = 1
        if(intval > maxlines) intval = maxlines

        write(28,1000) freq_id,int(intval*(iend - istartp + 1)),16,lat, &
                       mlong,elev,vitd_index,met_count,timeoffset,      &
                       timstep,mflag

! write the output file column headers
        write(28,'(''Year  JD  Hr   M   LVegTyp   LVegTmp     LVegEms'',&
       &''      LVegEvp  CTyp    CTopTmp    CMTmp    CLTmp   CAvEms  '',&   
       &''    CTotEvp   SurfTmp    Grtemp  Grmoistv         GrEvp  '',  &
       &''Snow+Ice     ATTmp    ASTmp'')')

        write(28,'(''                                 K              '',&   
       &''    kg/m^2*s              K        K        K              '',&
       &''  kg/m^2*s       K         K       % vol      kg/m^2*s     '',&
       &'' m          K        K'')')

        do i=istartp,iend
          swater = surfmoist(i) + surfice(i)

          if(i == istartp.or.intval == 1) then
            if(aint(dabs(met(i,ip_sd)-mflag)*1d5)*1d-5 > eps.and.       &
                    aint(dabs(met(i,ip_hi)-mflag)*1d5)*1d-5 > eps) then
              t2(i,8) = met(i,ip_sd) + met(i,ip_hi)
            else
              t2(i,8) = 0d0
            end if

            f2 = lheatf(i)/lhes(i)
            f3 = lheat(i)/lhes(i)
            write(28,1004) int(met(i,ip_year)),int(met(i,ip_doy)),      &
                 int(met(i,ip_hr)),int(met(i,ip_min)),vegl_type,ft(i),  &
                 surfemisf(i),f2,vegh_type,canopy_temp(1,i),            &
                 canopy_temp(2,i),canopy_temp(3,i),surfemisc(i),        &
                 cevap(i),tt(i),met(i,ip_tsoil),swater,f3,t2(i,8),      &
                airt(1,i),airt(2,i)
          else
            swater1 = surfmoist(i-1) + surfice(i-1)

            do j=1,intval
              if(j == intval) then
                if(aint(dabs(met(i,ip_sd)-mflag)*1d5)*1d-5 > eps.and.   &
                    aint(dabs(met(i,ip_hi)-mflag)*1d5)*1d-5 > eps) then
                  t2(i,8) = met(i,ip_sd) + met(i,ip_hi)
                else
                  t2(i,8) = 0d0
                end if

                f2 = lheatf(i)/lhes(i)
                f3 = lheat(i)/lhes(i)
                write(28,1004) int(met(i,ip_year)),int(met(i,ip_doy)),  &
                     int(met(i,ip_hr)),int(met(i,ip_min)),vegl_type,    &
                     ft(i),surfemisf(i),f2,vegh_type,canopy_temp(1,i),  &
                     canopy_temp(2,i),canopy_temp(3,i),surfemisc(i),    &
                     cevap(i),tt(i),met(i,ip_tsoil),swater,f3,t2(i,8),  &
                     airt(1,i),airt(2,i)
              else
                t1(j,1) = int(met(i-1,ip_year))
                t1(j,2) = int(met(i-1,ip_doy))
                t1(j,3) = int(met(i-1,ip_hr))
                t1(j,4) = int(met(i-1,ip_min) + j*vegint)

! corrections may be necessary due to inconsistant data recording
                k = 0
                do while (int(t1(j,4)) >= 60)
                  t1(j,4) = t1(j,4) - 60
                  k = k + 1
                end do
                t1(j,3) = t1(j,3) + k
                if(t1(j,3) >= 24) then
                  t1(j,3) = t1(j,3) - 24
                  t1(j,2) = t1(j,2) + 1
                end if

                mody = met(i,ip_year) - aint(met(i,ip_year)*2.5d-1)*4d0
                if(dabs(t1(j,2)-366d0) <= eps.and.dabs(mody) > eps) then
                  t1(j,2) = 1
                  t1(j,1) = t1(j,1) + 1
                else if(dabs(t1(j,2)-367d0) <= eps.and.                 &
                                               dabs(mody) <= eps) then
                  t1(j,2) = 1
                  t1(j,1) = t1(j,1) + 1
                end if

                f1 = (j*1d2/intval)*1d-2
                t2(j,1) = ft(i-1) + (ft(i) - ft(i-1))*f1                 !low veg temp

                do k=1,nclayers
                  t2(j,k+1) = canopy_temp(k,i-1)                        &
                              + (canopy_temp(k,i) - canopy_temp(k,i-1)) &
                                                                    *f1  !canopy temps
                end do

                t2(j,5) = tt(i-1) + (tt(i) - tt(i-1))*f1                 !surface temp
                t2(j,6) = met(i-1,ip_tsoil)                             &
                          + (met(i,ip_tsoil) - met(i-1,ip_tsoil))*f1     !ground temp 
                t2(j,7) = swater1 + (swater - swater1)*f1                !ground ice + water 
       
                if(aint(dabs(met(i,ip_sd)-mflag)*1d5)*1d-5 > eps.and.   &
                     aint(dabs(met(i,ip_hi)-mflag)*1d5)*1d-5 > eps) then
                  t2(j,8) = (met(i-1,ip_sd) + met(i-1,ip_hi))           & 
                             + ((met(i,ip_sd) + met(i,ip_hi))           &
                                - (met(i-1,ip_sd) + met(i-1,ip_hi)))*f1  !surface icing + snow depth
                else
                  t2(j,8) = 0d0
                end if

                t2(j,9) = airt(1,i-1) + (airt(1,i) - airt(1,i-1))*f1     !air temp above canopy
                t2(j,10) = airt(2,i-1) + (airt(2,i) - airt(2,i-1))*f1    !air temp below canopy
                t2(j,11) = surfemisf(i-1) + (surfemisf(i)               &
                                                   - surfemisf(i-1))*f1  !low veg emissivity
                t2(j,12) = surfemisc(i-1) + (surfemisc(i)               &
                                                   - surfemisc(i-1))*f1  !canopy emissivity

                f2 = lheatf(i-1)/lhes(i-1)
                t2(j,13) = f2 + (lheatf(i)/lhes(i) - f2)*f1              !low veg condensation(+)/evaporation(-)
                t2(j,14) = cevap(i-1) + (cevap(i) - cevap(i-1))*f1       !total canopy condensation(+)/evaporation(-)
                f3 = lheat(i-1)/lhes(i-1)
                t2(j,15) = f3 + (lheat(i)/lhes(i) - f3)*f1               !ground condensation(+)/evaporation(-)

                write(28,1004) t1(j,1),t1(j,2),t1(j,3),t1(j,4),         &
                     vegl_type,t2(j,1),t2(j,11),t2(j,13),vegh_type,     &
                     t2(j,2),t2(j,3),t2(j,4),t2(j,12),t2(j,14),t2(j,5), &
                     t2(j,6),t2(j,7),t2(j,15),t2(j,8),t2(j,9),t2(j,10)
              end if
            end do
          end if
        end do
      end if

! ******************************************************************************
! print the ground flux information -> Does not include vegetation component
! ******************************************************************************
      if(water_flag <= 1.and.flprint == 1) then
        write(24,1000) freq_id,(iend - istartp + 1),13,lat,mlong,elev,  &
                       vitd_index,met_count,timeoffset,timstep,mflag

        write(24,'(''Year  JD  Hr   M     SDown       SUp    IRdown  '',&
      &''    IRup     LHeat     SHeat     PHeat     CHeat    CHeat1'',  &      
      &''       Sum    SnDpth     MlTmp     SfTmp'')')
        write(24,'(''                     W/m^2     W/m^2     W/m^2  '',& 
      &''   W/m^2     W/m^2     W/m^2     W/m^2     W/m^2     W/m^2  '',&
      &''   W/m^2       m         K         K'')')

        do i=istartp,iend
          ttemp = tt(i)
          sd = met(i,ip_sd)
          if(i /= istartp) sd = met(i-1,ip_sd)

          write(24,1002) int(met(i,ip_year)),int(met(i,ip_doy)),        &
              int(met(i,ip_hr)),int(met(i,ip_min)),sdown(i),sup(i),     &
              irdown(i),irup(i),lheat(i),sheat(i),pheat1(i),cheat(i),   &
              cheat1(i),evap_heat(i),sd,tmelt(i),ttemp
        end do
      end if

! ******************************************************************************
! temporary for stdmod readin file
!     bedrock = 99.0
!     write(99,*) icnt,' 8 random / arbitrary terrain'
!     do i=1,icnt
!       write(99,3000) i,slippery,tan(slope*pi/180.0)*100.0,             &
!                node_type(1),kwi(i),t_strength(i),l_strength(i),bedrock
!       if(met(i,ip_sd) > eps.or.met(i,ip_hi) > eps.                   &
!           or.sstate(i) == 'f'.or.sstate(i) == 't') then
!         if(met(i,ip_hi) > eps) ice_code = 1
!         write(99,3001) met(i,ip_sd)*100.0/2.54,sdens,                  &
!              fr_depth(i)*100.0/2.54,tw_depth(i)*100.0/2.54,            &
!              t_moist(i),ice_code
!       end if
!     end do
! 3000 format(i8,'  SMRO  1.0000  0  ',i5,f8.2,'  0.1  300  ',a2,i5,     &
!             f8.3,f8.3,f8.3)
! 3001 format('   WIN ',f6.2,f6.2,f8.2,f8.2,f8.3,i3)
! ******************************************************************************

 1000 format(i10,1x,i6,1x,i3,1x,f10.6,1x,f11.6,1x,f11.6,1x,i8,1x,i8,1x, &
             f6.2,1x,f5.2,1x,f8.2)

 1001 format(4(i4),6(f8.2),i5,f8.2,i5,3(2(f8.2),i5),4(f8.2),1x,         &
             f8.2,1x,f8.2,1x,2(f8.2),2(f8.5),f8.2,2x,f8.5,2x,f8.5,      &
             3(f8.2),f8.4,f8.4,3x,a1,3x,f8.2,7x,i1,3x,i1,4x,f8.3,3x,    &
             f9.3,3x,f8.2,3x,i4)

 1002 format(4(i4),10(f10.2),f10.5,2(f10.2))

 1004 format(4(i4),3x,i4,5x,f8.2,3x,f9.3,3x,e10.4,1x,i4,3x,4(f9.2),3x,e10.4,1x,f9.2,2x,f8.2,2x,   &
             f8.5,4x,e10.4,1x,f8.4,2x,2(f9.2))

      end subroutine write_outputs

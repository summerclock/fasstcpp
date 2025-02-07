      subroutine datatrans(hdrlines,shdrlines)

      use met_global

! this subroutine calls the following subroutines:
!     word (appended to this subroutine)

! this subroutine reads and translates the meta file *.mta associated with the 
! meteorological file to get the column number of the various parameters
! and set unit flags

! NOTE: SEDRIS nomenclature is used where ever possible

! Parameter order in the file is:
!	number of header lines
!     Fractional Latitude (lat)(+ North) and longitude (long)(+ West)
!     location elevation (elev) (m)
!     timeoffset localtime-GMT or 9999 for missing (timeoffset) (hr)
!     met timestep in fractional hours (timstep) (hr)
!     missing flag (mflag)
!	number of observations
!	number of columns
!	number of positions
!	instrument height
!	start year
!     year column number (y_col)
!	day of the year column number (jday_col)
!     hour column number (hr_col)
!     minute column number (m_col)
!     atmospheric pressure column number (ap_col) (mbar)
!     temp column number (tmp_col) (C)
!     mean temp column number (tmp_col) (C): USED IF NO TEMP DATA
!     period for which mean temp is valid (SEDRIS nomenclature)
!         0 = unknown     7 = 1 min       13 = 1 hr
!         1 = 1 sec       8 = 2 min       14 = 2 hr
!         2 = 2 sec       9 = 5 min       15 = 3 hr
!         3 = 5 sec      10 = 10 imn      16 = 4 hr
!         4 = 10 sec     11 = 15 min      17 = 6 hr
!         5 = 15 sec     12 = 30 min      18 = 8 hr
!         6 = 30 sec                      19 = 12 hr
!         997 = withheld             20 = 24 hr
!         998 = not applicable
!     min/max temp column numbers (tmp1_col/tmp2_col) (C): USED IF NO TEMP DATA
!     RH column number (rh_col) (%)
!     dewpoint temp column number (dt_col) (C): USED IF NO RH DATA
!     dewpoint depression column number (dd_col) (C): USED IF NO RH 
!     wind speed column number (ws_col) (m/s)
!     wind direction column number (wd_col)
!     precipitation rate column number (prcp_col, prcp2_col) (mm/timestep)
!     precipitation type column number (pt_col,pt2_col) (SEDRIS nomenclature)
!         0    unknown         4   freezing rain      99   other
!         1    none            5   sleet
!         2    rain            6   hail
!         3    snow            7   graupel
!	NOTE: If there are 2 columns for precipitation, the second one is reserved for snow and the first for rain
!     precipitation accumulation column number (pa_col, pa2_col) (mm in a timestep)
!     low cloud amt column number (lcld_col) (fraction)
!     low cloud cover type column number (lct_col) (SEDRIS nomenclature)
!         0   none                                    6   stratus nebulosus and/or startus fractus
!         1   cumulus humulis or cumulus fractus      7   startus fractus and/or cumulus fractus of bad weather
!         2   cumulus mediocris or congestus          8   cumulus and stratocumulus
!         3   cumulonimbus calvus with/out cumulus, stratocumulus or stratus
!         4   stratocumulus cumulogentius             9   cumulonimbus capillatus
!         5   other stratocumulus types               99  no cloulds visible
!     low cloud hgt column number (lhgt_col) (km)
!     middle cloud amt column number (mcld_col) (fraction)
!     middle cloud type (mct_col) (SEDRIS nomenclature)
!         0   none                                   6    altocumulus cumulogentis or cumulonimbo genitus
!         1   altostratus translucidus               7    altocumulus translucidus or opacus, multi layer 
!         2   altostratus opacus or nimbostratus     8    altocumulus castellanus or floccus
!         3   altocumulus translucidus, 1 level      9    chaotic altocumulus
!         4   altocumulus translucidus, many levels, varying
!         5   altocumulus translucidus in bands      99   no clouds visible
!     middle cloud hgt column number (mhgt_col) (km)
!     high cloud amt column number (hcld_col) (fraction)
!     high cloud type (hct_col) (SEDRIS nomenclature)
!         0   none                                     5   cirrus and/or cirrostratus < 45 above horiz.
!         1   cirrus fibratus                          6   cirrus and/or cirrostratus > 45 above horiz.
!         2   cirrus spissatu, patchy                  7   cirrostratus full cover                  
!         3   cirrus spissatus cumulonimbo genitus     8   cirrostratus not full cover    
!         4   cirrus unicinus and/or fibratus          9   cirrocumulus
!         99  no clouds visible
!     high cloud hgt column number (hhgt_col) (km)
!     total solar (direct + diffuse) column number (tsol_col) (W/m^2)
!     Direct solar column number (dirsol_col) (W/m^2)
!     diffuse solar column number (difsol_col) (W/m^2)
!     reflected solar column number (upsol_col) (W/m^2)
!     IR dwn column number (ir_col) (W/m^2)
!     IR up column number (irup_col) (W/m^2)
!     solar zenith column number (zen_col) (degree, + Clockwise)
!     solar azimuth column number (az_col) (degree, + Clockwise)
!     soil surface temperature column number (tsoil_col) (C)
!     measured snow depth column number (sndepth_col) (m)
!	density of material column number (denmat_col)
!	precipitation diameter column number (precdiam_col)
!	visibility distance (km) column number (visdis_col) (km)
!	visibility type column number (vistype_col) (SEDRIS nomenclature)
!         0   unknown                      5   haze
!         1   none                         6   ocean spray
!         2   mist                         7   sand
!         3   dust                         8   volcanic ash
!         4   smoke                        9   other 

!	calculates the number of met columns based on the above (ncols)

      implicit none

      integer(kind=4),intent(out)::hdrlines,shdrlines

! local variables
      integer(kind=4):: jj,io,yr_found
      real(kind=8):: tmp_min_per,tmp_max_per
      character(len=120):: card
      character(len=100):: name
      character(len=35):: unit
      character(len=10):: var

! initialize global column values
      y_col = -901
      jday_col = -901
      hr_col = -901
      m_col = -901
      ap_col = -901
      tmp_col = -901
      tmp1_col = -901
      tmp2_col = -901
      rh_col = -901
      dt_col = -901
      dd_col = -901
      ws_col = -901
      wd_col = -901
      prcp_col = -901
      pa_col = -901
      pt_col = -901
      prcp2_col = -901
      pa2_col = -901
      pt2_col = -901
      lcld_col = -901
      lhgt_col = -901
      lct_col = -901
      mcld_col = -901
      mhgt_col = -901
      mct_col = -901
      hcld_col = -901
      hhgt_col = -901
      hct_col = -901
      tsol_col = -901
      dirsol_col = -901
      difsol_col = -901
      upsol_col = -901
      ir_col = -901
      irup_col = -901
      zen_col = -901
      az_col = -901
      tsoil_col = -901
      sndepth_col = -901
      visdis_col = -901
      vistype_col = -901
      precdiam_col = -901
      denmat_col = -901

      ave_period = 0
      iheight = 0
      year0 = 0
      yr_found = 0
      nlines = 0
      mxcol = 0
      ndatpos = 0

      lat = 0d0
      mlong = 0d0
      elev = 0d0
      timeoffset = 0d0
      timstep = 0d0
      mflag = 0d0

      ttflag = '0'
      prcp_flag = '0'
      prcp2_flag = '0'
      sdflag1 = '0'
      vdflag = '0'
      dflag = '0'
      pdflag = '0'

! initialize local and output variables
      hdrlines = 0
      shdrlines = 0
      jj = 0
      io = 0
      yr_found = 0
      tmp_min_per = 0d0
      tmp_max_per = 0d0

! determine which variable is in which column and the associated units
      do while(io /= -1)
        read(5,'(a)',iostat=io) card

        jj = len_trim(card)
        if(jj /= 0) then
          card = card(1:jj)
          call word(card,jj,name,var,unit)

          if(name == 'NHRD_LINES') then
            read(var,11) hdrlines

          else if(name == 'NSUBHRD_LINES') then
            read(var,11) shdrlines

          else if(name == 'EDCS_AC_SPATIAL_GEODETIC_LATITUDE') then
             read(var,12) lat

          else if(name == 'EDCS_AC_SPATIAL_GEODETIC_LONGITUDE') then
             read(var,12) mlong

          else if(name == 'EDCS_CC_TERRAIN_ELEVATION') then              !m
            read(var,12) elev
            if(unit == 'EDCS_UNITS_FOOT') elev = elev*3.048d-1

          else if(name == 'GMT-LOCAL_TIME') then
            read(var,*) timeoffset

          else if(name == 'TIME_STEP') then                              !hours
            read(var,*) timstep
            if(unit == 'EDCS_UNITS_MINUTE') timstep = timstep/60d0

          else if(name == 'MISSING_FLAG') then
            read(var,*) mflag

          else if(name == 'NO_OBSERVATIONS') then
            read(var,14) nlines

          else if(name == 'NO_DATA_COLUMNS') then
            read(var,16) mxcol

          else if(name == 'NO_POSITIONS') then
            read(var,18) ndatpos

          else if(name == 'INSTRUMENT_HEIGHT') then                      !m
            read(var,*)iheight
            if(unit == 'EDCE_UNITS_FOOT')iheight = iheight*3.048d-1

          else if(name == 'EDCS_AC_YEAR_COMMON_ERA'.and.yr_found == 0)  &
                                                                   then
            if(yr_found == 0) read(var,16) year0
            yr_found = yr_found + 1

          else if(name == 'EDCS_AC_YEAR_COMMON_ERA'.and.yr_found == 1)  &
                                                                   then 
            read(var,16) y_col

          else if(name == 'EDCS_AC_DAY_OF_YEAR') then
            read(var,16) jday_col

          else if(name == 'EDCS_AC_HOUR_OF_DAY') then
            read(var,16) hr_col

          else if(name == 'EDCS_AC_MINUTE_OF_HOUR') then
            read(var,16) m_col

          else if(name == 'EDCS_AC_ATMOSPHERIC_PRESSURE') then
            read(var,16) ap_col                                          !mbar
            if(ap_col /= -901) then
              apflag = 'mb'
              if(unit == 'EDCS_UNITS_PASCAL') apflag = 'pa'
            end if

          else if(name == 'EDCS_AC_AIR_TEMPERATURE') then
            read(var,16) tmp_col                                         !Celsius
            if(tmp_col /= -901.and.ttflag == '0') then
              ttflag = '1'
              if(unit == 'EDCS_UNITS_DEGREE_CELSIUS') tflag = 'tc'
              if(unit == 'EDCS_UNITS_DEGREE_FAHRENHEIT') tflag = 'tf'
              if(unit == 'EDCS_UNITS_KELVIN') tflag = 'tk'
            end if

          else if(name == 'EDCS_AC_AIR_TEMPERATURE_MEAN') then
            if(ttflag == '0') then
              read(var,16) tmp_col
              if(tmp_col /= -901) then
                ttflag = '2'
                if(unit == 'EDCS_UNITS_DEGREE_CELSIUS') tflag = 'tc'
                if(unit == 'EDCS_UNITS_DEGREE_FAHRENHEIT') tflag = 'tf'
                if(unit == 'EDCS_UNITS_KELVIN') tflag = 'tk'
              end if
            end if

          else if(name == 'EDCS_AC_AVERAGING_PERIOD') then
            read(var,16) ave_period

          else if(name == 'EDCS_AC_AIR_TEMPERATURE_MINIMUM') then
            if(ttflag == '0') then
              read(var,16) tmp1_col
              if(tmp1_col /= -901) then
                if(unit == 'EDCS_UNITS_DEGREE_CELSIUS') tflag = 'tc'
                if(unit == 'EDCS_UNITS_DEGREE_FAHRENHEIT') tflag = 'tf'
                if(unit == 'EDCS_UNITS_KELVIN') tflag = 'tk'
              end if
            end if

          else if(name == 'EDCS_AC_AIR_TEMPERATURE_MAXIMUM') then
            if(ttflag == '0') then
              read(var,16) tmp2_col
              if(tmp2_col /= -901) then
                if(unit == 'EDCS_UNITS_DEGREE_CELSIUS') t2flag = 'tc'
                if(unit == 'EDCS_UNITS_DEGREE_FAHRENHEIT') t2flag = 'tf'
                if(unit == 'EDCS_UNITS_KELVIN') t2flag = 'tk'
              end if
            end if

          else if(name == 'EDCS_AC_AIR_TEMPERATURE_MINIMUM_PERIOD') then
            read(var,12) tmp_min_per
            if(tmp1_col /= -901) then
              if(unit == 'EDCS_UNITS_HOUR') t3flag = 'h'
              if(unit == 'EDCS_UNITS_MINUTE') t3flag = 'm'
              if(unit == 'EDCS_UNITS_DAY') t3flag = 'd'
            end if

         else if(name == 'EDCS_AC_AIR_TEMPERATURE_MAXIMUM_PERIOD') then
            read(var,12) tmp_max_per
            if(tmp2_col /= -901) then
              if(unit == 'EDCS_UNITS_HOUR') t3flag = 'h'
              if(unit == 'EDCS_UNITS_MINUTE') t3flag = 'm'
              if(unit == 'EDCS_UNITS_DAY') t3flag = 'd'
            end if

          else if(name == 'EDCS_AC_RELATIVE_HUMIDITY') then
            read(var,16) rh_col                                          !%

          else if(name == 'EDCS_AC_DEWPOINT_TEMPERATURE') then
            read(var,16) dt_col                                          !C
            if(dt_col /= -901) then
              if(unit == 'EDCS_UNITS_DEGREE_CELSIUS') dtflag = 'tc'
              if(unit == 'EDCS_UNITS_DEGREE_FAHRENHEIT') dtflag = 'tf'
              if(unit == 'EDCS_UNITS_KELVIN') dtflag = 'tk'
            end if

          else if(name == 'EDCS_AC_DEWPOINT_DEPRESSION') then
            read(var,16) dd_col                                          !C
            if(dd_col /= -901) then
              if(unit == 'EDCS_UNITS_DEGREE_CELSIUS') ddflag = 'tc'
              if(unit == 'EDCS_UNITS_DEGREE_FAHRENHEIT') ddflag = 'tf'
              if(unit == 'EDCS_UNITS_KELVIN') ddflag = 'tk'
            end if

          else if(name == 'EDCS_AC_WIND_SPEED') then
            read(var,16) ws_col                                          !m/s
            if(ws_col /= -901) then
              if(unit == 'EDCS_UNITS_METERS_PER_SECOND') wsflag = 'mps'
              if(unit == 'EDCS_UNITS_FEET_PER_SECOND') wsflag = 'fps'
            end if

          else if(name == 'EDCS_AC_WIND_DIRECTION') then
            read(var,16) wd_col

          else if(name == 'EDCS_AC_PRECIPITATION_RATE') then
            read(var,16) prcp_col                                        !mm/hr
            if(prcp_col /= -901.and.prcp_flag /= '1') then
              prcp_flag = '1'
              if(unit == 'EDCS_UNITS_MILLIMETERS_PER_HOUR')             &
                                                         pflag = 'mmh'
              if(unit == 'EDCS_UNITS_METERS_PER_HOUR') pflag = 'mph'
              if(unit == 'EDCS_UNITS_FEET_PER_HOUR') pflag = 'fph'
              if(unit == 'EDCS_UNITS_INCHES_PER_HOUR') pflag = 'iph'
            end if

          else if(name == 'EDCS_AC_PRECIPITATION_TYPE') then
            read(var,16) pt_col

          else if(name == 'EDCS_AC_PRECIPITATION_ACCUMULATED') then
            if(prcp_flag == '0') then
              read(var,16) pa_col
              if(pa_col /= -901.and.prcp_flag /= '3') then
                prcp_flag = '3'
                if(unit == 'EDCS_UNITS_METER') paflag = 'a'
                if(unit == 'EDCS_UNITS_INCH') paflag = 'i'
                if(unit == 'EDCS_UNITS_FOOT') paflag = 'f'
                if(unit == 'EDCS_UNITS_MILLIMETER') paflag = 'm'
                if(unit == 'EDCS_UNITS_CENTIMETER') paflag = 'c'
              end if
            end if

          else if(name == 'EDCS_AC_PRECIPITATION_RATE_SNOW') then
            read(var,16) prcp2_col                                       !mm/hr
            if(prcp2_col /= -901.and.prcp2_flag /= '1') then
              prcp2_flag = '1'
              if(unit == 'EDCS_UNITS_MILLIMETERS_PER_HOUR')             &
                                                       p2flag = 'mmh'
              if(unit == 'EDCS_UNITS_METERS_PER_HOUR') p2flag = 'mph'
              if(unit == 'EDCS_UNITS_FEET_PER_HOUR') p2flag = 'fph'
              if(unit == 'EDCS_UNITS_INCHES_PER_HOUR') p2flag = 'iph'
            end if

          else if(name == 'EDCS_AC_PRECIPITATION_TYPE_SNOW') then
             read(var,16) pt2_col

          else if(name == 'EDCS_AC_PRECIPITATION_ACCUMULATED_SNOW') then
            if(prcp2_flag == '0') then
              read(var,16) pa2_col
              if(pa2_col /= -901) then
                prcp2_flag = '3'
                if(unit == 'EDCS_UNITS_METER') pa2flag = 'a'
                if(unit == 'EDCS_UNITS_INCH') pa2flag = 'i'
                if(unit == 'EDCS_UNITS_FOOT') pa2flag = 'f'
                if(unit == 'EDCS_UNITS_MILLIMETER') pa2flag = 'm'
                if(unit == 'EDCS_UNITS_CENTIMETER') pa2flag = 'c'
              end if
            end if

          else if(name == 'EDCS_AC_CLOUD_COVER_FRACTION_LOW') then
             read(var,16) lcld_col                                       !0 - 1

          else if(name == 'EDCS_AC_CLOUD_LOW_TYPE') then
            read(var,16) lct_col

          else if(name == 'EDCS_AC_CLOUD_BASE_HEIGHT_LOW_LEVEL') then
            read(var,16) lhgt_col                                        !km
            if(lhgt_col /= -901) then
              if(unit == 'EDCS_UNITS_KILOMETER') lcflag = 'k'
              if(unit == 'EDCS_UNITS_MILE') lcflag = 'm'
            end if

          else if(name == 'EDCS_AC_CLOUD_COVER_FRACTION_MIDDLE') then
            read(var,16) mcld_col

          else if(name == 'EDCS_AC_CLOUD_MIDDLE_TYPE') then
            read(var,16) mct_col

          else if(name == 'EDCS_AC_CLOUD_BASE_HEIGHT_MIDDLE_LEVEL') then
            read(var,16) mhgt_col
            if(mhgt_col /= -901) then
              if(unit == 'EDCS_UNITS_KILOMETER') mcflag = 'k'
              if(unit == 'EDCS_UNITS_MILE') mcflag = 'm'
            end if

          else if(name == 'EDCS_AC_CLOUD_COVER_FRACTION_HIGH') then
            read(var,16) hcld_col

          else if(name == 'EDCS_AC_CLOUD_HIGH_TYPE') then
            read(var,16) hct_col

          else if(name == 'EDCS_AC_CLOUD_BASE_HEIGHT_HIGH_LEVEL') then
            read(var,16) hhgt_col
            if(hhgt_col /= -901) then
              if(unit == 'EDCS_UNITS_KILOMETER') hcflag = 'k'
              if(unit == 'EDCS_UNITS_MILE') hcflag = 'm'
            end if

!! replace with SEDRIS nomenclature when new list comes out !!
          else if(name == 'RADIATION_SOLAR_GLOBAL') then
            read(var,16) tsol_col                                        !W/m^2

          else if(name == 'EDCS_AC_RADIATION_SOLAR_DIRECT') then
            read(var,16) dirsol_col

          else if(name == 'EDCS_AC_RADIATION_SOLAR_DIFFUSED') then
            read(var,16) difsol_col

!! replace with SEDRIS nomenclature when new list comes out !!
          else if(name == 'RADIATION_SOLAR_REFLECTED') then
            read(var,16) upsol_col

          else if(name == 'EDCS_AC_RADIATION_HEAT_FLUX_INFRARED') then
            read(var,16) ir_col

!! replace with SEDRIS nomenclature when new list comes out !!
          else if(name == 'RADIATION_HEAT_FLUX_INFRARED_REFLECTED') then
            read(var,16) irup_col
 
          else if(name == 'EDCS_AC_ANGLE_ELEVATION_SOLAR_HORIZONTAL')   &
                                                                   then
            read(var,16) zen_col

          else if(name == 'EDCS_AC_ANGLE_AZIMUTH_SOLAR_TRUE_NORTH') then
            read(var,16) az_col

          else if(name == 'EDCS_AC_SOIL_TEMPERATURE_SURFACE_LAYER') then
            read(var,16) tsoil_col                                       !K
            if(tsoil_col /= -901) then
              if(unit == 'EDCS_UNITS_KELVIN')stflag = 'tk'
              if(unit == 'EDCS_UNITS_DEGREE_FAHRENHEIT') stflag = 'tf'
              if(unit == 'EDCS_UNITS_DEGREE_CELSIUS') stflag = 'tc'
            end if

         else if(name == 'EDCS_AC_PRECIPITATION_SNOW_DEPTH') then
           read(var,16) sndepth_col
           if(sndepth_col /= -901) then
              sdflag1 = '1'
              if(unit == 'EDCS_UNITS_METER') sdflag2 = 'a'
              if(unit == 'EDCS_UNITS_FOOT') sdflag2 = 'f'
              if(unit == 'EDCS_UNITS_MILLIMETER') sdflag2 = 'm'
              if(unit == 'EDCS_UNITS_INCH') sdflag2 = 'i'
              if(unit == 'EDCS_UNITS_CENTIMETER') sdflag2 = 'c'
            end if

          else if(name == 'EDCS_AC_PRECIPITATION_SNOW_DEPTH_EQUIVALENT')   &
                                                                   then 
            if(sdflag1 == '0') then
              read(var,16) sndepth_col
              if(sndepth_col /= -901) then
                sdflag1 = '2'
                if(unit == 'EDCS_UNITS_METER') sdflag2 = 'a'
                if(unit == 'EDCS_UNITS_FOOT') sdflag2 = 'f'
                if(unit == 'EDCS_UNITS_MILLIMETER') sdflag2 = 'm'
                if(unit == 'EDCS_UNITS_INCH') sdflag2 = 'i'
                if(unit == 'EDCS_UNITS_CENTIMETER') sdflag2 = 'c'
              end if
            end if

          else if(name == 'EDCS_AC_DENSITY_MATERIAL') then
            read(var,16) denmat_col                                      !g/cm^3
            if(denmat_col /= -901) then
              if(unit == 'EDCS_UNITS_KILOGRAMS_PER_CUBIC_METER')        &
                dflag ='k'
              if(unit == 'EDCS_UNITS_POUNDS_PER_CUBIC_FOOT') dflag = 'p'
            end if

          else if(name == 'PRECIPITATION_DIAMETER') then
            read(var,16) precdiam_col                                    !mm
            if(precdiam_col /= -901) then
              if(unit == 'EDCS_UNITS_METERS') pdflag = 'm'
              if(unit == 'EDCS_UNITS_CENTIMETERS') pdflag = 'c'
              if(unit == 'EDCS_UNITS_FOOT') pdflag = 'f'
              if(unit == 'EDCS_UNITS_INCH') pdflag = 'i'
            end if

          else if(name == 'VISIBILITY_DISTANCE') then
            read(var,16) visdis_col
            if(visdis_col /= -901) then
              if(unit == 'EDCS_UNITS_MILE') vdflag = 'm'
            end if
 
          else if(name == 'VISIBILITY_TYPE') then
            read(var,16) vistype_col

          end if

        end if   !jj != 0
      end do   !while not eof

 11   format(i2)
 12   format(f8.2)
 14   format(i8)
 16   format(i4)
 18   format(i5)

! calculate the number of columns of met data
      ncols = max(y_col,jday_col,hr_col,m_col,ap_col,tmp_col,rh_col,    &
                  ws_col,wd_col,prcp_col,pa_col,pt_col,prcp2_col,       &
                  pa2_col,pt2_col,lcld_col,lhgt_col,lct_col,mcld_col,   &
                  mhgt_col,mct_col,hcld_col,hhgt_col,hct_col,tsol_col,  &
                  dirsol_col,difsol_col,upsol_col,ir_col,irup_col,      &
                  zen_col,az_col,sndepth_col,tsoil_col,vistype_col,     &
                  visdis_col,precdiam_col,denmat_col)

      if(ncols > mxcol) write(*,'(''Calculated number of columns does'',&
     &'' not agree with the maximum number of columns in the mta file.''&
     &)')

      if(ncols > mxcol)write(*,*)'calculated ',ncols,' mta file ',mxcol

! check to see if longitude and gmt-offset agree
      if(mlong*timeoffset > 0d0) then
        write(*,*) 'longitude does not agree with the gmt offset'
        stop
!        call exit(-1)
      end if

      end subroutine datatrans

! ******************************************************************
      SUBROUTINE WORD (A,N,NAME,VAR,UNIT)

!       GIVEN TEMP IMAGE (A),LENGTH OF A AND STARTING COLUMN (ISTR),
!       SCAN FOR NEXT WORD AND FIND IT'S END.  WORDS ARE
!       SEPARATED BY BLANKS.  ANY OTHER CHARACTER IS PART 
!       OF THE ALPHA STRING OR NUMBER.

!       RETURN ISTR* AS COLUMN OF BEGINNING OF WORD.
!              IEND* AS COLUMN OF END OF WORD.
!       CARD IMAGE (A) IS ASSUMMED TO BE 120 BYTES LONG.

      implicit none

      integer(kind=4),intent(in):: N
      character(len=N),intent(in):: A
      character,intent(out):: NAME*(*),VAR*(*),UNIT*(*)

! local variables
      INTEGER(kind=4):: I,ISTRA,IENDA,ISTRB,IENDB,ISTRC,IENDC
      CHARACTER(len=1):: A1

! initialize varaibles
      I = 0
      ISTRA = 0
      IENDA = 0
      ISTRB = 0
      IENDB = 0
      ISTRC = 0
      IENDC = 0

! FIND FIRST WORD
      I = 1
      ISTRA = I
      A1 = A(I:I)
      do while(I <= N)
        A1 = A(I:I)
        if(ichar(A1) == 9.or.ichar(A1) == 32) exit   !tab,blank
        I = I + 1
      end do

      IENDA = I - 1
      NAME = A(ISTRA:IENDA)

! FIND THE SECOND WORD
      do while(ichar(A1) == 9.or.ichar(A1) == 32)   !tab,blank
        A1 = A(I:I)
        I = I + 1
        IF (I > N) exit
      end do

      ISTRB = I - 1

      do while(I <= N)
        A1 = A(I:I)
        if(ichar(A1) == 9.or.ichar(A1) == 32) exit   !tab,blank
        I = I + 1
      end do

      IENDB = I - 1
      VAR = A(ISTRB:IENDB)

! FIND THE THIRD WORD
      do while(ichar(A1) == 9.or.ichar(A1) == 32)   !tab,blank
        A1 = A(I:I)
        I = I + 1
        IF (I > N) exit
      end do

      ISTRC = I - 1

      do while(I <= N)
        A1 = A(I:I)
        if(ichar(A1) == 9.or.ichar(A1) == 32) exit   !tab,blank
        I = I + 1
      end do

      IENDC = I - 1
      UNIT = A(ISTRC:IENDC)

      END SUBROUTINE WORD

      program bld_met_mta

c********************************************************************************
c This program builds the meta (.mta) file associated with a met data file
c********************************************************************************

      implicit none

	integer max_col
	parameter(max_col=99)

	integer ncol,funit,dtest,ivar,ihelp,icount,i,icol,npos
	integer col_array(max_col)
	real*4 rmflag,rvar
	character enter_file*1,file_name*160,munits*1

	parameter(rmflag=8888.0)
      data funit/90/


c does the user want to enter a new met file
      write(*,'(/,'' Do you want to enter a new met file, yes or no?
     &...'',$)')
	read(*,*) enter_file

	do while (enter_file.eq.'y'.or.enter_file.eq.'Y')
       do i=1,max_col            !initialize column array
          col_array(i) = 0
	  end do
	  icount = 0

        write(*,'(//,'' Enter the name of the meta file, up to 160 
     &characters long.'',/''   NOTE: It must have a .mta extension.'',/)
     &')
        read(*,*) file_name
        open(unit=funit,file=file_name,status='unknown')

	  write(*,'(//,'' To get help for any entry enter 8888. You will b
     &e asked about general data'',/,'' and site information first.'',//
     &)')

c number of header lines
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Number of header lines at the top of the data s
     &et? The maximum is 99...'',$)')
	    read(*,*) ivar
	    if(ivar.eq.int(rmflag).or.ivar.gt.99) then
	      call help(1)
	    else
	      ihelp = 1
	    end if
	  end do 
	  write(funit,'(''NHRD_LINES'',43x,i2,5x,''EDCS_UNITS_UNITLESS'')
     &)') ivar

c number of sub-header lines in a multi-point file
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Number of sub-header lines between the data set
     &s? The maximum is 99...'',$)')
	    read(*,*) ivar
	    if(ivar.eq.int(rmflag).or.ivar.gt.99) then
	      call help(42)
	    else
	      ihelp = 1
	    end if
	  end do 
	  write(funit,'(''NSUBHRD_LINES'',40x,i2,5x,''EDCS_UNITS_UNITLESS
     &''))') ivar

c site latitude
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the site latitude, Positive values North,
     & Negative values'',/,'' South....'',$)')
	    read(*,*) rvar
	    if(rvar.eq.rmflag) then
	      call help(2)
	    else
	      ihelp = 1
	    end if

	    if(abs(rvar).gt.90.0) then
		  ihelp = 0
	      write(*,*) 'latitude out of range, try again.'
	    end if
	  end do
	  write(funit,'(''EDCS_AC_SPATIAL_GEODETIC_LATITUDE'',19x,f6.2,2x,
     &''EDCS_UNITS_DEGREE_ARC'')') rvar

c site longitude
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the site longitude, Positive values East 
     &from Zulu, Negative values'',/,'' West....'',$)')
	    read(*,*) rvar
	    if(rvar.eq.rmflag) then
	      call help(3)
	    else
	      ihelp = 1
	    end if

	    if(abs(rvar).gt.180.0) then
	      ihelp = 0
	      write(*,*) 'longitude out of range, try again.'
	    end if
	  end do
	  write(funit,'(''EDCS_AC_SPATIAL_GEODETIC_LONGITUDE'',17x,f7.2,
     &2x,''EDCS_UNITS_DEGREE_ARC'')') rvar

c site elevation
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the site elevation.....'',$)')
	    read(*,*) rvar
          if(rvar.eq.rmflag) then
	      call help(4)
	    else
	      ihelp = 1
	    end if
	  end do

c elevation units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, feet or meters
     &...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_CC_TERRAIN_ELEVATION'',26x,f6.1,3x,''ED
     &CS_UNITS_METER'')') rvar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_CC_TERRAIN_ELEVATION'',25x,f7.1,3x,''ED
     &CS_UNITS_FOOT'')') rvar
          else
            write(*,*) 'Wrong units, try again.'
          end if
	  end do

c GMT offset
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the time offset between local time and
     & Greenwich Mean Time. Negative '',/,'' values West from
     & Zulu, Positive values East in hours...'',$)')
          read(*,*) rvar
	    if(rvar.eq.rmflag) then
	      call help(5)
	    else
	      ihelp = 1
	    end if
	  end do
	  write(funit,'(''GMT-LOCAL_TIME'',39x,f4.1,3x,''EDCS_UNITS_HOUR''
     &)') rvar
 
c data time step
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the data time step...'',$)')
          read(*,*) rvar
	    if(rvar.eq.rmflag) then
	      call help(6)
	    else
	      ihelp = 1
	    end if
	  end do

c data time step units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, minutes or 
     &hours...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''TIME_STEP'',44x,f4.1,3x,''EDCS_UNITS_MINUTE
     &'')') rvar
	    else if(munits.eq.'h'.or.munits.eq.'H') then
	      ihelp = 1
	      write(funit,'(''TIME_STEP'',44x,f4.1,3x,''EDCS_UNITS_HOUR'')
     &') rvar
          else
            write(*,*) 'Wrong units, try again.'
          end if
        end do

c missing data flag
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the missing data flag value...'',$)')
          read(*,*) rvar
	    if(rvar.eq.rmflag) then
	      call help(7)
	    else
	      ihelp = 1
	    end if
        end do
	  write(funit,'(''MISSING_FLAG'',38x,f7.1,3x,''EDCS_UNITS_UNITLESS
     &'')') rvar 

c number of data points
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the number of observations...'',$)')
          read(*,*) ivar
          if(ivar.eq.int(rmflag)) then
	      call help(8)
	    else
	      ihelp = 1
	    end if
	  end do
	  write(funit,'(''NO_OBSERVATIONS'',33x,i7,5x,''EDCS_UNITS_UNITLES
     &S'')') ivar

c number of data columns
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the number of data columns...'',$)')
          read(*,*) ncol
	    if(ncol.eq.int(rmflag)) then
	      call help(9)
          else
	      ihelp = 1
	    end if
	  end do
        write(funit,'(''NO_DATA_COLUMNS'',37x,i3,5x,''EDCS_UNITS_UNITLES
     &S'')') ncol

c number of meteorological data collection positions
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the number of meteorological data
     & collection positions...'',$)')
          read(*,*) npos
	    if(npos.eq.int(rmflag)) then
	      call help(40)
          else
	      ihelp = 1
	    end if
	  end do
        write(funit,'(''NO_POSITIONS'',40x,i3,5x,''EDCS_UNITS_UNITLES
     &S'')') npos

c instrument height
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the average instrument height...'',$)')
          read(*,*) rvar
	    if(rvar.eq.rmflag) then
	      call help(10)
	    else
	      ihelp = 1
	    end if
	  end do

c instrument height units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, feet or meters
     &...'',$)')
	    read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''INSTRUMENT_HEIGHT'',36x,f4.1,3x,''EDCS_UNITS
     &_METER'')') rvar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''INSTRUMENT_HEIGHT'',36x,f4.1,3x,''EDCS_UNITS
     &_FOOT'')') rvar
          else
            write(*,*) 'Wrong units, try again.'
          end if
	  end do

c year data was taken
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the year when the data set begins in YYYY
     & format...'',$)')
	    read(*,*) ivar
	    if(ivar.eq.int(rmflag)) then
	      call help(11)
	    else
	      ihelp = 1
	    end if
	  end do
	  write(funit,'(''EDCS_AC_YEAR_COMMON_ERA'',28x,i4,5x,''EDCS_UNITS
     &_YEAR'')') ivar

	  write(funit,'(''************************Meteorological_Parameter
     &_Column_Number************************'')')

        write(*,'(/,'' The remaining prompts will ask in which column
     & the  data is found, followed'',/,'' by the units of measurement.
     & If no data exists, enter -901.'',/,'' Remember to type 8888 for
     & help.'',/)')

c year column
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the YEAR data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(0)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_YEAR_COMMON_ERA'',28x,i4,5x,''EDCS_
     &UNITS_YEAR'')') ivar

c day column
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the DAY OF THE YEAR data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(0)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_DAY_OF_YEAR'',32x,i4,5x,''EDCS_UNITS_
     &DAY'')') ivar

c hour column
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the HOUR OF THE DAY data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(12)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_HOUR_OF_DAY'',32x,i4,5x,''EDCS_UNITS_
     &HOUR'')') ivar

c minute column
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the MINUTE OF THE HOUR data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(13)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_MINUTE_OF_HOUR'',29x,i4,5x,''EDCS_UNITS_
     &MINUTE'')') ivar

c air pressure
       icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the ATMOSPHERIC PRESSURE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(14)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_ATMOSPHERIC_PRESSURE'',23x,i4,5x,''
     &EDCS_UNITS_MILLIBAR'')') ivar

c air temperature
	  write(*,'(/,'' AIR TEMPERATURE data will be asked for next. You
     & can enter these several'',/,'' ways: 1- regular, 2- mean, 3- min
     & and max.'',/)')

c "typical" temperatures
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the AIR TEMPERATURE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(15)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c air temperature units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, Kelvin or Cels
     &ius or Fahrenheit...'',$)')
	    read(*,*) munits
	    if(munits.eq.'K'.or.munits.eq.'k') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE'',28x,i4,5x,''EDCS_
     &UNITS_KELVIN'')') ivar
          else if(munits.eq.'C'.or.munits.eq.'c') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE'',28x,i4,5x,''EDCS_
     &UNITS_DEGREE_CELSIUS'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE'',28x,i4,5x,''EDCS_
     &UNITS_DEGREE_FAHRENHEIT'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c mean air temperature
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the MEAN AIR TEMPERATURE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(16)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c units for mean air temperature
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, Kelvin or Celsi
     &us or Fahrenheit...'',$)')
	    read(*,*) munits
	    if(munits.eq.'K'.or.munits.eq.'k') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MEAN'',23x,i4,5x,''
     &EDCS_UNITS_KELVIN'')') ivar
          else if(munits.eq.'C'.or.munits.eq.'c') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MEAN'',23x,i4,5x,''
     &EDCS_UNITS_DEGREE_CELSIUS'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MEAN'',23x,i4,5x,''
     &EDCS_UNITS_DEGREE_FAHRENHEIT'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c period over which mean air temperature is valid
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the Mean Air Temperature AVERAGING PERIOD
     & CODE. See the associated'',/,'' help menu...'',$)')
	    read(*,*) ivar
	    if(ivar.eq.int(rmflag))then
	      call help(17)
	    else
	     ihelp = 1
	    end if
        end do
	  write(funit,'(''EDCS_AC_AVERAGING_PERIOD'',27x,i4,5x,''EDCS_
     &UNITS_ENUMERATION'')') ivar

c minimum air temperture
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the MINIMUM AIR TEMPERATURE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(18)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c units for minimum air temperature
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, Kelvin or Celsi
     &us or Fahrenheit...'',$)')
	    read(*,*) munits
	    if(munits.eq.'K'.or.munits.eq.'k') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM'',20x,i4,5x,
     &''EDCS_UNITS_KELVIN'')') ivar
          else if(munits.eq.'C'.or.munits.eq.'c') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM'',20x,i4,5x,
     &''EDCS_UNITS_DEGREE_CELSIUS'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM'',20x,i4,5x,
     &''EDCS_UNITS_DEGREE_FAHRENHEIT'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c period over which the minimum air temperature is valid
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the TIME PERIOD for the Minimum Air
     & Temperature. NOTE: Enter the'',/,'' actual time value...'',$)')
          read(*,*) rvar
	    if(rvar.eq.rmflag)then
	      call help(19)
          else
	      ihelp = 1
	    end if
	  end do

c units for the minimum air temperature period
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, hours or minut
     &es or days...'',$)')
	    read(*,*) munits
	    if(munits.eq.'h'.or.munits.eq.'H') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM_PERIOD'',14x
     &,f5.1,3x,''EDCS_UNITS_HOUR'')') rvar
	    else if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM_PERIOD'',14x
     &,f5.1,3x,''EDCS_UNITS_MINUTE'')') rvar
	    else if(munits.eq.'d'.or.munits.eq.'D') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM_PERIOD'',14x
     &,f5.1,3x,''EDCS_UNITS_DAY'')') rvar
	    else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do	

c maximum air temperature
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the MAXIMUM AIR TEMPERATURE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(18)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c units for maximum air temperature
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, Kelvin or Celsi
     &us or Fahrenheit...'',$)')
	    read(*,*) munits
	    if(munits.eq.'K'.or.munits.eq.'k') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM'',20x,i4,5x,
     &''EDCS_UNITS_KELVIN'')') ivar
          else if(munits.eq.'C'.or.munits.eq.'c') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM'',20x,i4,5x,
     &''EDCS_UNITS_DEGREE_CELSIUS'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM'',20x,i4,5x,
     &''EDCS_UNITS_DEGREE_FAHRENHEIT'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c period over which the maximum air temperature is valid
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the TIME PERIOD for the Maximum Air
     & Temperature. NOTE: Enter the'',/,'' actual time value...'',$)')
          read(*,*) rvar
	    if(rvar.eq.rmflag)then
	      call help(19)
          else
	      ihelp = 1
	    end if
	  end do

c units for the maximum air temperature period
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, hours or minut
     &es or days...'',$)')
	    read(*,*) munits
	    if(munits.eq.'h'.or.munits.eq.'H') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM_PERIOD'',14x
     &,f5.1,3x,''EDCS_UNITS_HOUR'')') rvar
	    else if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM_PERIOD'',14x
     &,f5.1,3x,''EDCS_UNITS_MINUTE'')') rvar
	    else if(munits.eq.'d'.or.munits.eq.'D') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_AIR_TEMPERATURE_MINIMUM_PERIOD'',14x
     &,f5.1,3x,''EDCS_UNITS_DAY'')') rvar
	    else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do	 

c relative humidity
        write(*,'(/,'' Moisture information is entered next. There are
     & three ways of'',/,'' describing it: 1- relative humidity,'',/,
     &'' 2- dewpoint temperature,'',/,'' or 3- dewpoint depression.'',/)
     &')

        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the RELATIVE HUMIDITY data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(20)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c relative humidity units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' If the Relative Humidity values are between 0
     & and 1: enter 1.'',/,'' If the Relative Humidity values are betwee
     &n 0 and 100: enter 2...'',$)')
	    read(*,*) dtest
	    if(dtest.eq.1) then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_RELATIVE_HUMIDITY'',26x,i4,5x,''EDCS
     &_UNITS_UNITLESS'')') ivar
	    else if(dtest.eq.2) then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_RELATIVE_HUMIDITY'',26x,i4,5x,''EDCS
     &_UNITS_PERCENT'')') ivar
	    else
	      write(*,*) 'Variable out of range, try again.'
	    end if
	  end do

c dewpoint temperature
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the DEWPOINT TEMPERATURE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(21)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c dewpoint temperature units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, Kelvin or Celsi
     &us or Fahrenheit...'',$)')
	    read(*,*) munits
	    if(munits.eq.'K'.or.munits.eq.'k') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DEWPOINT_TEMPERATURE'',23x,i4,5x,''
     &EDCS_UNITS_KELVIN'')') ivar
          else if(munits.eq.'C'.or.munits.eq.'c') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DEWPOINT_TEMPERATURE'',23x,i4,5x,''
     &EDCS_UNITS_DEGREE_CELSIUS'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DEWPOINT_TEMPERATURE'',23x,i4,5x,''
     &EDCS_UNITS_DEGREE_FAHRENHEIT'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c dewpoint depression
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)

	    !get a possible column input
          write(*,'(/,'' Enter the DEWPOINT DEPRESSION data column...'',
     &$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(21)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar

c dewpoint depression units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, Kelvin or Celsi
     &us or Fahrenheit...'',$)')
	    read(*,*) munits
	    if(munits.eq.'K'.or.munits.eq.'k') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DEWPOINT_DEPRESSION'',24x,i4,5x,''
     &EDCS_UNITS_KELVIN'')') ivar
          else if(munits.eq.'C'.or.munits.eq.'c') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DEWPOINT_DEPRESSION'',24x,i4,5x,''
     &EDCS_UNITS_DEGREE_CELSIUS'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DEWPOINT_DEPRESSION'',24x,i4,5x,''
     &EDCS_UNITS_DEGREE_FAHRENHEIT'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c wind speed
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the WIND SPEED data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(22)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c wind speed units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, meters per seco
     &nd or feet per second...'',$)')
	    read(*,*) munits
          if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_WIND_SPEED'',33x,i4,5x,''EDCS_UNITS_
     &METERS_PER_SECOND'')') ivar
          else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_WIND_SPEED'',33x,i4,5x,''EDCS_UNITS_
     &FEET_PER_SECOND'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c wind direction
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the WIND DIRECTION data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(22)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_WIND_DIRECTION'',29x,i4,5x,''EDCS_UNITS_
     &DECREE_ARC'')') ivar

c precipitation information
	  write(*,'(/,'' You are allowed 2 columns for precipitation. If
     & you do so, the first column is for rain/freezing rain and the
     & second is for snow.'',$)')

c precipitation rate
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the PRECIPITATION RATE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(23)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c precipitation rate units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, millimeters per
     & hour-m-'',/,'' or meters per hour-t-'',/,'' or inches per hour-i-
     &'',/,'' feet per hour-f-...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_RATE'',25x,i4,5x,''
     &EDCS_UNITS_MILLIMETERS_PER_HOUR'')') ivar
          else if(munits.eq.'t'.or.munits.eq.'T') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_RATE'',25x,i4,5x,''
     &EDCS_UNITS_METERS_PER_HOUR'')') ivar
          else if(munits.eq.'i'.or.munits.eq.'I') then
		  ihelp = 1		 
	      write(funit,'(''EDCS_AC_PRECIPITATION_RATE'',25x,i4,5x,''
     &EDCS_UNITS_INCHES_PER_HOUR'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_RATE'',25x,i4,5x,''
     &EDCS_UNITS_FEET_PER_HOUR'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c precipitation rate & precipitation accumulation types
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the PRECIPITATION TYPE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(24)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_PRECIPITATION_TYPE'',25x,i4,5x,''EDCS_
     &UNITS_ENUMERATION'')') ivar

c precipitation accumulation
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the PRECIPITATION ACCUMULATION data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(23)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c precipitation accumulation units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, millimeters-m-'
     &',/,'' or meters-t- or inches-i-'',/,'' or feet-f-'',/,''
     & or centimeters-c-...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED'',18x,i4,
     &5x,''EDCS_UNITS_MILLIMETER'')') ivar
          else if(munits.eq.'t'.or.munits.eq.'T') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED'',18x,i4,
     &5x,''EDCS_UNITS_METER'')') ivar
          else if(munits.eq.'i'.or.munits.eq.'I') then
		  ihelp = 1		 
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED'',18x,i4,
     &5x,''EDCS_UNITS_INCH'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED'',18x,i4,
     &5x,''EDCS_UNITS_FOOT'')') ivar
          else if(munits.eq.'c'.or.munits.eq.'C') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED'',18x,i4,
     &5x,''EDCS_UNITS_CENTIMETER'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c snow precipitation rate
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SNOW PRECIPITATION RATE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(23)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c snow precipitation rate units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, millimeters per
     & hour-m-'',/,'' or meters per hour-t-'',/,'' or inches per hour-i-
     &'',/,'' feet per hour-f-...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_RATE_SNOW'',20x,i4,5x,
     &''EDCS_UNITS_MILLIMETERS_PER_HOUR'')') ivar
          else if(munits.eq.'t'.or.munits.eq.'T') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_RATE_SNOW'',20x,i4,5x,
     &''EDCS_UNITS_METERS_PER_HOUR'')') ivar
          else if(munits.eq.'i'.or.munits.eq.'I') then
		  ihelp = 1		 
	      write(funit,'(''EDCS_AC_PRECIPITATION_RATE_SNOW'',20x,i4,5x,
     &''EDCS_UNITS_INCHES_PER_HOUR'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_RATE_SNOW'',20x,i4,5x,
     &''EDCS_UNITS_FEET_PER_HOUR'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c snow precipitation rate & snow precipitation accumulation types
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SNOW PRECIPITATION TYPE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(24)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
		  else
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          else
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_PRECIPITATION_TYPE_SNOW'',20x,i4,5x,''
     &EDCS_UNITS_ENUMERATION'')') ivar

c snow precipitation accumulation
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SNOW PRECIPITATION ACCUMULATION data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(23)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c snow accumulation units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, millimeters-m-'
     &',/,'' or meters-t- or inches-i-'',/,'' or feet-f-'',/,''
     & or centimeters-c-...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED_SNOW'',13x
     &,i4,5x,''EDCS_UNITS_MILLIMETER'')') ivar
          else if(munits.eq.'t'.or.munits.eq.'T') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED_SNOW'',13x
     &,i4,5x,''EDCS_UNITS_METER'')') ivar
          else if(munits.eq.'i'.or.munits.eq.'I') then
		  ihelp = 1		 
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED_SNOW'',13x
     &,i4,5x,''EDCS_UNITS_INCH'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED_SNOW'',13x
     &,i4,5x,''EDCS_UNITS_FOOT'')') ivar
          else if(munits.eq.'c'.or.munits.eq.'C') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_ACCUMULATED_SNOW'',13x
     &,i4,5x,''EDCS_UNITS_CENTIMETER'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c low cloud fraction
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the LOW CLOUD COVER FRACTION data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(25)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_CLOUD_COVER_FRACTION_LOW'',19x,i4,5x,
     &''EDCS_UNITS_UNITLESS'')') ivar

c low cloud type
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the LOW CLOUD COVER TYPE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(26)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_CLOUD_LOW_TYPE'',29x,i4,5x,''EDCS_UNITS_
     &ENUMERATION'')') ivar

c low cloud height
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the LOW CLOUD COVER HEIGHT data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(27)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c low cloud height units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, kilometers or
     & miles...'',$)')
	    read(*,*) munits
	    if(munits.eq.'k'.or.munits.eq.'K') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_CLOUD_BASE_HEIGHT_LOW_LEVEL'',16x,
     &i4,5x,''EDCS_UNITS_KILOMETER'')') ivar
	    else if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_CLOUD_BASE_HEIGHT_LOW_LEVEL'',16x,
     &i4,5x,''EDCS_UNITS_MILE'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c middle cloud fraction
      icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the MIDDLE CLOUD COVER FRACTION data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(25)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_CLOUD_COVER_FRACTION_MIDDLE'',16x,i4,5x,
     &''EDCS_UNITS_UNITLESS'')') ivar

c middle cloud type
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the MIDDLE CLOUD COVER TYPE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(28)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_CLOUD_MIDDLE_TYPE'',26x,i4,5x,''EDCS_
     &UNITS_ENUMERATION'')') ivar

c middle cloud height
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the MIDDLE CLOUD COVER HEIGHT data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(27)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c middle cloud height units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, kilometers or
     & miles...'',$)')
	    read(*,*) munits
	    if(munits.eq.'k'.or.munits.eq.'K') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_CLOUD_BASE_HEIGHT_MIDDLE_LEVEL'',
     &13x,i4,5x,''EDCS_UNITS_KILOMETER'')') ivar
	    else if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_CLOUD_BASE_HEIGHT_MIDDLE_LEVEL'',
     &13x,i4,5x,''EDCS_UNITS_MILE'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c high cloud fraction
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the HIGH CLOUD COVER FRACTION data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(25)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_CLOUD_COVER_FRACTION_HIGH'',18x,i4,5x,
     &''EDCS_UNITS_UNITLESS'')') ivar

c high cloud type
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the HIGH CLOUD COVER TYPE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(29)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_CLOUD_HIGH_TYPE'',28x,i4,5x,''EDCS_UNITS
     &_ENUMERATION'')') ivar

c high cloud height
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the HIGH CLOUD COVER HEIGHT data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(27)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c high cloud height units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, kilometers or
     & miles...'',$)')
	    read(*,*) munits
	    if(munits.eq.'k'.or.munits.eq.'K') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_CLOUD_BASE_HEIGHT_HIGH_LEVEL'',15x,
     &i4,5x,''EDCS_UNITS_KILOMETER'')') ivar
	    else if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_CLOUD_BASE_HEIGHT_HIGH_LEVEL'',15x,
     &i4,5x,''EDCS_UNITS_MILE'')')ivar
	    else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c radiation information
        write(*,'(/,'' Short and Long wave radiation information is
     & input next. The units for all'',/,'' radiation measurements are
     & Watts per square meter.'',/)')

c total incoming shortwave radiation
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the INCOMING GLOBAL SOLAR RADIATION
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(30)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''RADIATION_SOLAR_GLOBAL'',29x,i4,5x,''EDCS_UNITS_
     &WATTS_PER_SQUARE_METER'')') ivar

c direct shortwave solar radiation
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the INCOMING DIRECT SOLAR RADIATION
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(30)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_RADIATION_SOLAR_DIRECT'',21x,i4,5x,''
     &EDCS_UNITS_WATTS_PER_SQUARE_METER'')') ivar

c diffuse shortwave radiation
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the INCOMING DIFFUSE SOLAR RADIATION
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(30)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_RADIATION_SOLAR_DIFFUSED'',19x,i4,5x,
     &''EDCS_UNITS_WATTS_PER_SQUARE_METER'')') ivar

c reflected shortwave radiation
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the REFLECTED SOLAR RADIATION
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(31)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''RADIATION_SOLAR_REFLECTED'',26x,i4,5x,''EDCS_
     &UNITS_WATTS_PER_SQUARE_METER'')') ivar

c incoming IR
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the DOWNWELLING INFRARED RADIATION
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(32)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_RADIATION_HEAT_FLUX_INFRARED'',15x,i4,
     &5x,''EDCS_UNITS_WATTS_PER_SQUARE_METER'')') ivar

c emitted IR
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the EMITTED INFRARED RADIATION
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(33)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''RADIATION_HEAT_FLUX_INFRARED_REFLECTED'',13x,i4,
     &5x,''EDCS_UNITS_WATTS_PER_SQUARE_METER'')') ivar

c solar zenith angle
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SOLAR ELEVATION (ZENITH) ANGLE
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(34)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_ANGLE_ELEVATION_SOLAR_HORIZONTAL'',11x,
     &i4,5x,''EDCS_UNITS_DEGREE_ARC'')') ivar

c solar azimuth angle
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SOLAR ASPECT (AZIMUTH) ANGLE
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(35)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''EDCS_AC_ANGLE_AZIMUTH_SOLAR_TRUE_NORTH'',13x,i4,
     &5x,''EDCS_UNITS_DEGREE_ARC'')') ivar

c soil surface temperature
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SOIL SURFACE TEMPERATURE
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(36)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c soil surface temperature units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, Kelvin or Cels
     &ius or Fahrenheit...'',$)')
	    read(*,*) munits
	    if(munits.eq.'K'.or.munits.eq.'k') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_SOIL_TEMPERATURE_SURFACE_LAYER'',
     &13x,i4,5x,''EDCS_UNITS_KELVIN'')') ivar
          else if(munits.eq.'C'.or.munits.eq.'c') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_SOIL_TEMPERATURE_SURFACE_LAYER'',
     &13x,i4,5x,''EDCS_UNITS_DEGREE_CELSIUS'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_SOIL_TEMPERATURE_SURFACE_LAYER'',
     &13x,i4,5x,''EDCS_UNITS_DEGREE_FAHRENHEIT'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c snow depth
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SNOW DEPTH
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(37)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c snow depth units
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the units of measurement, millimeters-m-
     &'',/,'' or meters-t-'',/,'' or inches-i-'',/,'' or feet-f-'',/,''
     & or centimeters-c-...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH'',19x,i4,
     &5x,''EDCS_UNITS_MILLIMETER'')') ivar
          else if(munits.eq.'t'.or.munits.eq.'T') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH'',19x,i4,
     &5x,''EDCS_UNITS_METER'')') ivar
          else if(munits.eq.'i'.or.munits.eq.'I') then
		  ihelp = 1		 
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH'',19x,i4,
     &5x,''EDCS_UNITS_INCH'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH'',19x,i4,
     &5x,''EDCS_UNITS_FOOT'')') ivar
          else if(munits.eq.'c'.or.munits.eq.'C') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH'',19x,i4,
     &5x,''EDCS_UNITS_CENTIMETER'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c snow water equivalent (SWE)
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SNOW WATER EQUIVALENT DEPTH
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(37)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c SWE units
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the units of measurement, millimeters-m-
     &'',/,'' or meters-t-'',/,'' or inches-i-'',/,'' or feet-f-'',/,''
     & or centimeters-c-...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH_EQUIVALENT'
     &',8x,i4,5x,''EDCS_UNITS_MILLIMETER'')') ivar
          else if(munits.eq.'t'.or.munits.eq.'T') then
            ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH_EQUIVALENT'
     &',8x,i4,5x,''EDCS_UNITS_METER'')') ivar
          else if(munits.eq.'i'.or.munits.eq.'I') then
		  ihelp = 1		 
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH_EQUIVALENT'
     &',8x,i4,5x,''EDCS_UNITS_INCH'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH_EQUIVALENT'
     &',8x,i4,5x,''EDCS_UNITS_FOOT'')') ivar
          else if(munits.eq.'c'.or.munits.eq.'C') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_PRECIPITATION_SNOW_DEPTH_EQUIVALENT'
     &',8x,i4,5x,''EDCS_UNITS_CENTIMETER'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c snow density
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the SNOW DENSITY
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(38)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c snow density units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, grams per cm^3
     &'',/,'' or kilograms per m^3'',/,'' or pounds per ft^3...'',$)')
	    read(*,*) munits
	    if(munits.eq.'g'.or.munits.eq.'G') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DENSITY_MATERIAL'',27x,i4,5x,''EDCS_
     &UNITS_GRAMS_PER_CUBIC_CENTIMETER'')') ivar
	    else if(munits.eq.'k'.or.munits.eq.'K') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DENSITY_MATERIAL'',27x,i4,5x,''EDCS_
     &UNITS_KILOGRAMS_PER_CUBIC_METER'')') ivar
	    else if(munits.eq.'p'.or.munits.eq.'P') then
	      ihelp = 1
	      write(funit,'(''EDCS_AC_DENSITY_MATERIAL'',27x,i4,5x,''EDCS_
     &UNITS_POUNDS_PER_CUBIC_FOOT'')') ivar
	    else
	      write(*,*) 'Wrong units, try again.'
	    end if
	  end do

c snow precipitation diameter
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the PRECIPITATION DIAMETER
     & data column...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(39)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c snow precipitation diameter units
        ihelp = 0
	  do while (ihelp.eq.0)
          write(*,'(/,'' Enter the units of measurement, millimeters-m-'
     &',/,'' or meters-t-'',/,'' or inches-i-'',/,'' or feet-f-'',/,''
     & or centimeters-c-...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	      write(funit,'(''PRECIPITATION_DIAMETER'',
     &29x,i4,5x,''EDCS_UNITS_MILLIMETER'')') ivar
          else if(munits.eq.'t'.or.munits.eq.'T') then
	      ihelp = 1
	      write(funit,'(''PRECIPITATION_DIAMETER'',
     &29x,i4,5x,''EDCS_UNITS_METER'')') ivar
          else if(munits.eq.'i'.or.munits.eq.'I') then
	      ihelp = 1		 
	      write(funit,'(''PRECIPITATION_DIAMETER'',
     &29x,i4,5x,''EDCS_UNITS_INCH'')') ivar
	    else if(munits.eq.'f'.or.munits.eq.'F') then
	      ihelp = 1
	      write(funit,'(''PRECIPITATION_DIAMETER'',
     &29x,i4,5x,''EDCS_UNITS_FOOT'')') ivar
          else if(munits.eq.'c'.or.munits.eq.'C') then
	      ihelp = 1
	      write(funit,'(''PRECIPITATION_DIAMETER'',
     &29x,i4,5x,''EDCS_UNITS_CENTIMETER'')') ivar
          else
	      write(*,*) 'Wrong units, try again.'
          end if
	  end do

c visibility distance
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the VISIBILITY DISTANCE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(0)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns

c visibility units
        ihelp = 0
	  do while (ihelp.eq.0)
 	    write(*,'(/,'' Enter the units of measurement, miles or kilome
     &ters...'',$)')
          read(*,*) munits
	    if(munits.eq.'m'.or.munits.eq.'M') then
	      ihelp = 1
	  write(funit,'(''VISIBILITY_DISTANCE'',32x,i4,5x,''EDCS_UNITS_
     &MILE'')') ivar
	    else if(munits.eq.'k'.or.munits.eq.'K') then
	      ihelp = 1
	  write(funit,'(''VISIBILITY_DISTANCE'',32x,i4,5x,''EDCS_UNITS_
     &KILOMETER'')') ivar
          else
            write(*,*) 'Wrong units, try again.'
          end if
	  end do

c visibility type
        icount = icount + 1
        ihelp = 0
	  do while (ihelp.eq.0)
          
	    !get a possible column input
	    write(*,'(/,'' Enter the VISIBILITY TYPE data column
     &...'',$)')
          read(*,*) ivar

          !check if user asked for help instead
	    if(ivar.eq.int(rmflag)) then
	      call help(41)
		else
	      if(ivar.gt.ncol) then !column entry too large
	        write(*,'(/,'' Column number too large. try
     & again...'',$)')
	        ivar = rmflag !placeholder value, request new input from user
	      end if

	      if(ivar.ne.-901) then !don't check to see if -901 is used twice
		    do i=1,icount !check all other previous columns for duplicates
			  if(ivar.eq.col_array(i)) then
		        write(*,'(/,'' Column number already used. try
     & again...'',$)')
                  ivar = rmflag !placeholder value, request new input from user
	          end if
	        end do
	      end if
	    end if

		!if acceptable column was input, ivar is not equal rmflag
          if(ivar.ne.int(rmflag)) ihelp = 1
	  end do
	  col_array(icount) = ivar !add new column to list of used columns
	  write(funit,'(''VISIBILITY_TYPE'',36x,i4,5x,''EDCS_UNITS_
     &ENUMERATION'')') ivar


c subroutine lister, verifies correct columns used!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	call lister(col_array, icount, ncol)

c does the user want to enter the information for another met file
        write(*,'(/,'' Do you want to enter another met file, yes or 
     & no?...'',$)')
	  read(*,*) enter_file
	  if(enter_file.eq.'y'.or.enter_file.eq.'Y') funit = funit + 1

	end do  !while (enter_file.eq.'y'.or.enter_file.eq.'Y')

 1000 format(i4)

      stop
      end

c********************************************************************************
c********************************************************************************
      subroutine help(ic)

	implicit none

	integer,intent(in):: ic


	select case (ic)


	case(0)    !if no help
	  write(*, 
     &'(/,''No help available.'',/)')

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      case(1)   !number of header lines
        write(*,
     &'(/,'' Make sure that you include the column names and units and
     & all blank lines'',/, '' in the total number of header lines.'',/)
     &')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(2)   !latitude
        write(*,
     &'(/,'' The values range from +90 to -90 and must be in decimal
     & form, to 2 places.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(3)   !longitude
        write(*,
     &'(/,'' The values range from +180 to -180 and must be in decimal
     & form, to 2 places.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(4)   !elevation
        write(*,
     &'(/,'' The elevation may be above (+) or below (-) sea level in
     & either meters or feet. You will be asked for the units of 
     & measurement next.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(5)   !GMT offset
        write(*,
     &'(/,'' This should be a positive number between 0 and 24 hours.'',
     &/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(6)   !time step
        write(*,
     &'(/,'' The time step units should be in minutes or hours. You will
     & be asked for'',/,'' the units next.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(7)   !missing data
        write(*,
     &'(/,'' The missing data flag should be the same for the whole data
     & set.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(8)   !number of data lines
        write(*,
     &'(/,'' This is the number of actual data lines in the met file.'',
     &/,'' The maximum number of lines is 9,999,999.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(9)   !number of columns in input met file
        write(*,
     &'(/,'' This is the number of columns of data in the met file.'',/,
     &'' Maximum number of columns is 999.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(10)  !instrument height (assumes that all are the same)
        write(*,
     &'(/,'' Enter the instrument height at which the air temperature,
     & wind and dew point'',/,'' temperature or relative humidity were
     & measured.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(11)  !year
        write(*,
     &'(/,'' This should be a 4 digit value.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(12)  !hour and/or hour + minutes
        write(*,
     &'(/,'' This should be in 24 hour format (hhmm) if minutes are
     & included in the time stamp, otherwise it should go from 0 - 23
     & or 1 - 24. If minutes are included here, do not enter them as a
     & seperate data column'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(13)  !minutes only
        write(*,
     &'(/,'' Only enter a non -901 value if the minutes are in a 
     & seperate column from the hour data.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(14)  !air pressure
        write(*,
     &'(/,'' You will be asked for the units of air pressure next. The
     & acceptable units'',/,'' are millibars or Pascals.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(15)  !air temperature
        write(*,
     &'(/,'' This is the measured air temperature. You will be asked for
     & the units next.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(16)  !mean air temperature
        write(*,
     &'(/,'' This is the mean of the measured air temperature. You will
     & be asked for'',/,'' the units next followed by the averaging
     & period.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(17)  !period over which case(16) applies
        write(*,
     &'(/,'' Instead of entering the actual time over which the air
     & temperature is'',/,'' averaged, the SEDRIS codes are entered.
     & The codes are:'',/,
     &''   0 = unknown     7 = 1 min       13 = 1 hr'',/,
     &''   1 = 1 sec       8 = 2 min       14 = 2 hr'',/,
     &''   2 = 2 sec       9 = 5 min       15 = 3 hr'',/,
     &''   3 = 5 sec      10 = 10 imn      16 = 4 hr'',/,
     &''   4 = 10 sec     11 = 15 min      17 = 6 hr'',/,
     &''   5 = 15 sec     12 = 30 min      18 = 8 hr'',/,
     &''   6 = 30 sec                      19 = 12 hr'',/,
     &'' 997 = data withheld               20 = 24 hr'',/,
     &'' 998 = not applicable'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(18)  !min/max air temperature
        write(*,
     &'(/,'' This is the minimum or maximum of the measured air
     & temperature over a'',/,'' given time period. The units of
     & measurement will be entered next, followed'',/,'' by the time
     & period.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(19)  !period over which case(18) applies
        write(*,
     &'(/,'' Enter the actual time period for which the minimum and
     & maximum measured'',/,'' air temperature are given. The units of
     & measurement will be entered next.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(20)  !relative humidity
        write(*,
     &'(/,'' The relative humidity is the ratio of the actual vapor
     & pressure in the air'',/,'' to the maximum capacity of the vapor
     & pressure in the air. It can be measured'',/,'' from 0.0 - 1.0 or
     & from 0.0 - 100.0%.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(21)  !dewpoint temperature
        write(*,
     &'(/,'' The dewpoint (wet bulb) temperature is the temperature at
     & which saturation'',/,'' occurs. The difference between the air
     & temperature (dry bulb) and the'',/,'' dewpoint temperature is the
     & dewpoint depression.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(22)  !wind speed
        write(*,
     &'(/,'' The wind speed is often measured at a reference level of 2m
     & above the ground.'',/,'' The wind direction is measured +
     & clockwise from North in degrees.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(23)  !precipitation rate
        write(*,
     &'(/,'' The precipitation rate is a measure of how fast the
     & precipitation is'',/,'' accumulating versus the accumulated
     & precipitation which is a measure of how'',/,'' much has actually
     & fallen during a given time interval. The precipitation'',/,
     &'' data can be entered either way. The units of measurement for
     & both methods'',/,'' will be asked for after the column
     & information.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
      case(24)  !precipitation type
        write(*,
     &'(/,'' The precipitation type is a code describing the form of the
     & precipitation,'',/,'' i.e., snow, rain, none. SEDRIS nomenclature
     & is used. The SEDRIS codes are:'',/,
     &''        0 = unknown             5 = sleet'',/,
     &''        1 = none                6 = hail'',/,
     &''        2 = rain                7 = graupel'',/,
     &''        3 = snow               99 = other'',/,
     &''        4 = freezing rain'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

      case(25)  !cloud cover fraction
        write(*,
     &'(/,'' Cloud cover fraction is either measured in 8ths or 10ths.
     &'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(26)  !low cloud type
        write(*,
     &'(/,'' Low clouds are water clouds. The possible cloud types
     & following SEDRIS'',/,'' nomenclature are:'',/,'' 0  none'',/,
     &'' 1  cumulus humulis or cumulus fractus'',/,
     &'' 2  cumulus mediocris or congestus'',/,
     &'' 3  cumulonimbus calvus with/out cumulus, stratocumulus or
     & stratus'',/,'' 4  stratocumulus cumulogentius'',/,
     &'' 5  other stratocumulus types'',/,
     &'' 6  stratus nebulosus and/or startus fractus'',/,
     &'' 7  startus fractus and/or cumulus fractus of bad weather'',/,
     &'' 8  cumulus and stratocumulus'',/,
     &'' 9  cumulonimbus capillatus'',/,
     &'' 99 no cloulds visible'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(27)  !cloud height
        write(*,
     &'(/,'' Cloud height is measured from the ground to the bottom of
     & the cloud layer.'',/,'' Both low and middle clouds are water
     & clouds while high clouds are ice clouds.'',/,'' The range of base
     & heights at mid latitudes for the different cloud types is:'',/,
     &'' Low clouds:         0 - 6,500  ft or     1 - 1,980  m'',/,
     &'' Middle clouds:  6,500 - 23,000 ft or 1,980 - 7,000  m'',/,
     &'' High clouds:   16,500 - 45,000 ft or 5,000 - 13,700 m'',/ )')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(28)  !middle cloud type
        write(*,
     &'(/,'' Middle clouds are water clouds. The possible cloud types
     & using SEDRIS'',/,'' nomenclature are:'',/,'' 0  none'',/,
     &'' 1  altostratus translucidus'',/,
     &'' 2  altostratus opacus or nimbostratus'',/,
     &'' 3  altocumulus translucidus, 1 level'',/,
     &'' 4  altocumulus translucidus, many levels, varying'',/,
     &'' 5  altocumulus translucidus in bands'',/,
     &'' 6  altocumulus cumulogentis or cumulonimbo genitus'',/,
     &'' 7  altocumulus translucidus or opacus, multi layer'',/,
     &'' 8  altocumulus castellanus or floccus'',/,
     &'' 9  chaotic altocumulus'',/,
     &'' 99 no clouds visible'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(29)  !high cloud type
        write(*,
     &'(/,'' High clouds are ice clouds. The possible cloud types
     & using SEDRIS'',/,'' nomenclature are:'',/,'' 0  none'',/,
     &'' 1  cirrus fibratus'',/,
     &'' 2  cirrus spissatu, patchy'',/,                  
     &'' 3  cirrus spissatus cumulonimbo genitus'',/,
     &'' 4  cirrus unicinus and/or fibratus'',/,
     &'' 5  cirrus and/or cirrostratus < 45 above horiz.'',/,
     &'' 6  cirrus and/or cirrostratus > 45 above horiz.'',/,
     &'' 7  cirrostratus full cover'',/,
     &'' 8  cirrostratus not full cover'',/,
     &'' 9  cirrocumulus'',/,
     &'' 99 no clouds visible'',/)') 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(30)  !total incoming solar (shortwave) radiation
        write(*,
     &'(/,'' Incoming (downwelling) solar radiation consists of
     & direct and a diffuse'',/,'' components. The global solar
     & radiation is the sum of the two.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(31)  !reflected solar (shortwave) radiation
        write(*,
     &'(/,'' The reflected solar radiation is equal to the following:''
     &,/,''      (incoming solar radiation)*(1 - albedo)'',/,'' The net
     & solar radiation is equal to the incoming minus the reflected.'',/
     &)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(32)  !incoming longwave (IR) radiation
        write(*,
     &'(/,'' This is the infrared (longwave) radiation emmitted by the
     & sky.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(33)  !emitted longwave (IR) radiation
        write(*,
     &'(/,'' This is the infrared (longwave) radiation emmitted by the
     & ground. It is'',/,'' proportional to the temperature of the
     & emissivity*ground^4. The net infrared'',/,'' radiation is the
     & incoming (sky) minus the emitted (ground).'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(34)  !solar zenith angle
        write(*,
     &'(/,'' The zenith (elevation, altitude) angle is the angle that a
     & direct ray from'',/,'' the sun makes with the horizontal. It is
     & measured from the normal to the'',/,'' surface to the suns ray.''
     &,/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(35)  !solar azimuth angle
        write(*,
     &'(/,'' The solar azimuth (aspect) angle is the angle that a
     & horizontal'',/,'' projection of a direct ray from the sun makes
     & with true (geodetic) north.'',/,'' Positive angles are clockwise
     & from north.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(36)  !soil surface temperature
        write(*,
     &'(/,'' This is taken to be the temperature of the top of the soil.
     & It is not'',/,'' a depth averaged temperature.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(37)  !snow depth
        write(*,
     &'(/,'' Snow depth can be entered as the quantity directly measured
     &, or its water'',/,'' equivalent. The later is about a third of
     & the actual amount. The snow water'',/,'' equivalent depth data
     & column will be asked for next.'',/)') 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(38)  !snow density
        write(*,
     &'(/,'' This is the density of the snow as a function of time.'',/)
     &')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(39)  !precipitation diameter
        write(*,
     & '(/,'' This is intended for snow and is used by some models to
     & calculate the new'',/,'' snow density.'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(40)  !number of meteorological data collection positions
        write(*,
     & '(/,'' This is the number of meterological data collection
     & positions'',/)')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      case(41)  !visibility type
        write(*,
     & '(/,'' The possible visibility types using SEDRIS nomenclature ar
     &e:'',/,'' 0  unknown'',/,'' 1  none'',/,'' 2  mist'',/,     
     &'' 3  dust'',/,'' 4  smoke'',/,'' 5  haze'',/,'' 6  ocean spray''
     &,/,'' 7  sand'',/,'' 8  volcanic ash'',/,'' 9  other'',/)') 
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      case(42)  !number of sub-header lines
        write(*,
     &'(/,'' In a multi-point met file, this is the number of header lin
     &es lines between data groups.'',/)
     &')
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end select

	return
	end




c********************************************************************************
c********************************************************************************
	!This Subroutine will sort the columns used, and report back to the user
	!for verification (in order!) the columns that are reported. -901s are
	!not included. The program above has checked that no column is used twice.

	subroutine lister(col_array, count, ncol)

	implicit none
	integer,intent(in)::count !number of columns(including -901's)
	integer,intent(in)::ncol !number of columns of data user claimed above
	integer,intent(in)::col_array(count) !entered array of columns, unsorted
	
	integer I, J !loop variables for old, new arrays
	integer sorted(count) !new, sorted array, add as we go, initialized
							!to be big enough to fit as if no -901's
	integer alreadyadded, minsofar !already added is lowest already added
						!to new array. minthispass is the next smallest found
						!so far through the old array
	integer max !largest number in old array. says when we're done


	!!!!!!!!!!!!!!!used when searching for unused columns!!!!!!!!!!!!!!!!
	integer done !says when all actually used have been compared
					!0 is not finished, 1 is finished
	integer n !counts off columns said to be used
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
	J = 0 !get ready to add to sorted array
	alreadyadded = 0 !haven't added any
	max = col_array(1) !largest yet found is first item.

	!finds max in old array, so we know when to stop
	do 3037 I=1, count
		if (col_array(I).gt.max) then
			max = col_array(I)
		else
		end if
 3037	continue

	minsofar = max !start case for search. we should find lower things soon.

	do while (alreadyadded.lt.max)!until we've added up to max

		do 3050 I=1, count !look through old list
			!find the min not yet added
			if (col_array(I).lt.minsofar.and.
     &col_array(I).gt.alreadyadded) 
     &        then
				minsofar = col_array(I)
			else
			endif
 3050		continue 

		J = J+1 !go to next spot in sorted array
		sorted(J) = minsofar !add that min to the next sorted bit
		alreadyadded=minsofar !remember we've added it
		minsofar = max !reset min counter. we haven't found anything yet
	end do


	print*, 'For verification purposes:'
	print*, 'You have entered data from the following columns:'


	do 3065 I=1, J !re-using I, we're printing through sorted list up to last
					!meaningful value (J)
		write(*,*) sorted(I)
 3065	continue
	

	!compares columns used to columns claimed would be used. prints out columns unused.
	print*, 'you did not use these columns:'
	
	!Two staged. First searches for holes within/below sorted columns
	!then adds columns said would be used greater than max used

	I = 1 !counts off array locations of columns in sorted
	done = 0 !see declaration.
	do 3092 n=1, ncol !integer column numbers which should have been used

		if (done.eq.1) then
			print*, n !if we've passed all used columns, then
							!the remaining ones are unused. print them.
		else
			!check to see if we have skipped some columns
			if(sorted(I).gt.n) then
				print*, n !this number was skipped.
			else
				if (I.lt.J) then !there's still more in sorted array
					I = I + 1
				else
					done = 1 !used up all actually used slots in array
				end if
			end if
		end if
		
 3092	continue


	return
	end



      program make_user_soil_file

      implicit none

	     integer(kind=4):: ians,ihelp,icount,io,stype
	     real(kind=4):: rmflag,rvar
	     character(len=160):: usersoilfile
      character(len=4):: name
      character(len=120):: dump

	     data rmflag/8888.0/


	     write(*,'(/,'' Enter the file name where the data will be stored: &
            &'')')
	     read(*,'(a)') usersoilfile
      io = 0
      open(unit=90,file=usersoilfile,iostat=io,status='unknown')
      if(io /=0) then
        write(*,'('' Error opening user soil file, error ='',i3)') io
        stop
      end if

! find the end of the file
      icount = 0
      io = 0
      do while(io /= -1)
  	     read(90,'(a)',iostat=io) dump
	       icount = icount + 1
	     end do

      if(icount > 0) write(90,'(''  '')')

! start taking in the new data
      ians = 1
	     do while(ians == 1)

!       soil id
 	      write(*,'(/,'' Enter the new soil material name. The total name &
	             &'',/,''should not exceed four alpha-numeric characters ''&
	             &,/,'' all characters should be capatilized.....  '',$)')
	       read(*,*) name
	       write(90,'(a4,18x,''!new soil name'')') name

!       soil type
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' numerical soil classification ?'',/,'' For help&
                & enter 88. '',/, $)')
	         read(*,*) stype
	         if(stype == 88) then
	           call help(1)
	         else
	           ihelp = 1
	         end if
	       end do
	       write(90,'(4x,i3,15x,''!soil type and class  '')') stype

!       help information
        write(*,'(/,'' To get help for any entry enter 8888. If the ''  &
              &,/,''value is still unknown after getting help, enter 999&
              &.0 '',//)')

!       density
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Bulk density of dry material(g/cm^3)?....  '', &
                $)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(2)
	         else
	           ihelp = 1
	          end if
        end do
        write(90,'(3x,f8.3,11x,''!bulk density of dry material (g/cm^3) &
              &       '')') rvar

!       specific gravity, intrinsic density
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Intrinsic density of dry material(g/cm^3)?.... &
                & '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(3)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(3x,f8.3,11x,''!intrinsic density of dry material (g/c&
              &m^3)'')') rvar

!       solids fraction
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Volume fraction of solids (0.0 - 1.0)?....  '',&
                $)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(4)
	         else
	           ihelp = 1
	         end if
        end do
	       write(90,'(4x,f7.3,11x,''!volume fraction of solids (0.0 - 1.0) &
              &'')') rvar

!       porosity
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Porosity (0.0 - 1.0)?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(5)
	         else
	           ihelp = 1
	         end if
        end do
	       write(90,'(4x,f7.3,11x,''!porosity (0.0 - 1.0)'')') rvar

!       albedo
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Surface Albedo at normal incidence (0.0 - 1.0)?&
                &.....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(6)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(4x,f6.2,12x,''!Surface albedo-normal incidences'')') &
                                                                   rvar

!       emissivity
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Longwave surface emissivity (0.0 - 1.0)?.....  &
                &'',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(7)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(4x,f6.2,12x,''!Longwave surface emissivity'')') rvar

!       quartz content
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Quartz Content(0.0 - 1.0)?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(8)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(4x,f7.3,11x,''!quartz content(0.0 - 1.0)'')') rvar

!       organic fraction
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Organic fraction of solids (0.0 - 1.0)?....  ''&
                ,$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(9)
	         else
	           ihelp = 1
	         end if
        end do
	       write(90,'(4x,f7.3,11x,''!organic fraction (0.0 - 1.0)'')') rvar

!       thermal conductivity
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Thermal conductivity of the bulk dry material (&
                &W/m*K)?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(10)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(4x,f8.4,10x,''!thermal conductivity of the bulk dry m&
              &aterial(W/m*K)'')') rvar

!       specific heat
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Specific heat of dry material(J/kg*K)?....  '' &
                ,$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(11)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(3x,f7.2,12x,''!specific heat of dry material (J/kg*K)&
              &'')') rvar

!       saturated hydraulic conductivity
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Saturated hydraulic conductivity(cm/s)?....  ''&
                ,$)')
	        read(*,*) rvar
	        if(rvar == rmflag) then
	           call help(12)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(4x,f13.9,5x,''!saturated hydraulic conductivity (cm/&
              &s)'')') rvar

!       residual water content
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Residual water content (vol/vol) (0.0 - 1.0)?..&
                &..  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(13)
	         else
	           ihelp = 1
	         end if
        end do
	       write(90,'(4x,f7.3,11x,''!residual water content (vol/vol) (0.0 &
              &- 1.0)'')') rvar

!       maximum water content
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Maximum water content (vol/vol) (0.0 - 1.0)?...&
                &.  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(14)
	         else
	           ihelp = 1
	         end if
        end do
	       write(90,'(4x,f7.3,11x,''!maximum water content (vol/vol) (0.0 -&
              & 1.0)'')') rvar

!       van Genuchten pressure head
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' van Genuchten bubbling pressure head(cm)?....  &
                &'',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(15)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(2x,f7.1,13x,''!van Genuchten Bubbling pressure head (&
              &cm)'')') rvar

!       van Genuchten exponent
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' van Genuchten "n"?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(16)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(4x,f7.3,11x,''!van Genuchten exponent(n)'')') rvar

!       percent sand
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Percent sand?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(17)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(4x,f7.3,11x,''!percent sand'')') rvar

!       percent silt
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Percent silt?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(18)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(2x,f10.4,10x,''!percent silt'')') rvar

!       percent clay
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Percent clay?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(19)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(4x,f7.3,11x,''!percent clay'')') rvar

!       percent carbon content
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Percent carbon content?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(20)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(2x,f10.4,10x,''!carbon content'')') rvar

!       plastic limit
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Plastic limit?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(21)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(2x,f10.4,10x,''!plastic limit'')') rvar

!     percent passing #200 sieve
        ihelp = 0
        do while(ihelp == 0)
          write(*,'(/,'' Percent passing #200 sieve?....  '',$)')
	         read(*,*) rvar
	         if(rvar == rmflag) then
	           call help(22)
	         else
	           ihelp = 1
	         end if
        end do
        write(90,'(2x,f10.4,10x,''!percent passing #200 sieve'')') rvar

! write a blank line
        write(90,'(''  '')')

	       write(*,'(/,'' Would you like to added another material type? ye&
              &s=1, no=2.....  '')')
        read(*,*) ians
      end do  !while(ians == 1)

	close(90)

	stop
	end

! ********************************************************************************
	subroutine help(ic)

	implicit none
	
	integer,intent(in):: ic


	     select case (ic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	       case(0)    !if no help
	         write(*,'(/,'' No help available.'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(1)
	         write(*,'(/,'' The available USCS soil types are:  '',/,      &
                '' 1  GW Well graded gravels, gravel-sand mixtures, litt&
                         &le or no fines '',/,&
                '' 2  GP Poorly graded gravels, gravel-sand mixtures, li&
                         &ttle or no fines '',/,&
                '' 3  GM Silty gravels, gravel-sand-silt mixtures '',/, &
                '' 4  GC Clayey gravels, gravel-sand-clay mixtures '',/,&
                '' 5  SW Well graded sands, gravelly sands, little or no&
                         & fines '',/,&
                '' 6  SP Poorly graded sands, gravelly sands, little or &
                         &no fines'',/,&
                '' 7  SM Silty sands, sand-silt mixtures'',/,           &
                '' 8  SC Clayey sands, sand-clay mixtures '',/,         &
                '' 9  ML Inorganic silts and very fine sands, rock flour&
                         &, silty or clayey '',/,''        fine sands or&
                         & calyey silts with slight plasticity '',/,    &
                '' 10 CL Inorganic clays of low to medium plasticity, gr&
                         &avelly clays,'',/,''        sandy clays, silty&
                         & clays, lean clays '',/,                      &
                '' 11 OL Organic silts and organic silty clays of low pl&
                         &asticity '',/,&
                '' 12 CH Inorganic silts, micaceous or diatomaceous fine&
                         & sandy or silty soils, '',/,''        elastic &
                         &silts '',/,                                   &
                '' 13 MH Inorganic clays of high plasticity, fat clays  &
                         &'',/,&
                '' 14 OH Organic clays of medium to high plasticity, org&
                         &anic silts '',/,&
                '' 15 Pt Peat and other highly organic soils '',/,      &
                '' 16 MC SMSC Mixed SM and SC soils '',/,               &
                '' 17 CM CLML Mixed CL and ML soils '',/,               &
                '' 18 EV Evaporites (Salt Pans, similar properties to SW&
                         &) '',/,&
                '' 20 CO Concrete '',/,'' 21 AS Asphalt '',/,           &
                '' 25 RO Bed rock '',/,'' 26 WA Water '',/,             &
                '' 27 AI Air (Used for modeling bridges) '',/,          &
                '' 30 SN Permanent snow field or glacier'',/,/,         &
                '' The available USDA soil types are: '',/,             &
                '' 101 Sand '',/,'' 102 Loamy Sand '',/,                &
                '' 103 Sandy Loam '',/,'' 104 Silt Loam '',/,           &
                '' 105 Silt '',/,'' 106 Loam '',/,                      &
                '' 107 Sandy Clay Loam '',/,'' 108 Silty Clay Loam'',/, &
                '' 109 Clay Loam '',/,'' 110 Sandy Clay '',/,           &
                '' 111 Silty Clay '',/,'' 112 Clay '',/,                &
                '' 113 Peat '',/,'' 114 Water '',/,'' 115 Bedrock '',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(2)
	         write(*,'(/,'' The Bulk density of dry materials is the mass p&
                &er total volume (g/cm^3). '',/,'' (Example silty-sand =&
                & 1.49)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(3)
	         write(*,'(/,'' The intrinsic density of dry materials is the m&
                &ass per volume of'',/,'' solids (g/cm^3). (Example silt&
                &y-sand = 2.73)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(4)
          write(*,'(/,'' The volume fraction of solids equals the ratio &
                &bulk density/intrinsic'',/,'' density. (Example silty-s&
                &and = 0.546)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(5)
          write(*,'(/'' The porosity equals one minus the ratio of the b&
                &ulk density/intrinsic'',/,'' density. (Example silty-sa&
                &nd = 0.454)'',/)')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(6)
	         write(*,'(/,'' The albedo is the fraction of the total solar r&
                &adiation that is reflected'',/,'' when the surface is v&
                &iewed at normal incidence. (Example silty-sand = 0.4)''&
                ,/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(7)
	         write(*,'(/,'' The longwave surface emissivity is the fraction&
                & of thermal emitted radiation'',/,'' relative to that o&
                &f a blackbody at the same temperature.'',/,'' (Example &
                &silty-sand = 40.92)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(8)
	         write(*,'(/,'' The quartz is used to calculate the thermal con&
                &ductivity of the dry material.'',/,'' Suggested values &
                &for nonorganic silts is 0.35, nonorganic clays is 0.05,&
                &'',/,'' organic silts and clays is 0.20, gravels is 0.6&
                &5 and sands is 0.80.'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(9)
          write(*,'(/,'' The organic fraction is expressed as a fraction&
                & of the solids. '',/,'' It is roughly twice the carbon &
                &content. '',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(10)
	         write(*,'(/,'' The thermal conductivity of the dry material(W/&
           &m*K) is the proportionality'',/,'' factor relating heat flow&
           & to the temperature gradient.'',/,'' (Example silty-sand = 0&
           &.281)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(11)
	         write(*,'(/,'' The specific heat of dry materials (J/kg*K) is &
                &the number of Joules needed'',/,'' to raise one kilogra&
                &m of dry material one degree Kelvin.'',/,'' (Example sa&
                &nd = 730)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(12)
          write(*,'(/'' The saturated hydraulic conductivity(cm/s) is a &
                &proportionality factor'',/,'' relating fluid flow to it&
                &s gradient. (Example silty-sand = 0.0007987)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(13)
          write(*,'(/'' The residual water content (vol/vol) is the amou&
                &nt of water that is always'',/,'' there, even in arid c&
                &onditions. (Example sand = 0.026)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(14)
          write(*,'(/'' The maximum water content (vol/vol) is the water&
                & content under saturated'',/,'' conditions. (Example si&
                &lty-sand = 0.408)'',/)')
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(15)
          write(*,'(/'' The van Genuchten Bubbling pressure head(cm) is &
                &the pressure needed to'',/,'' force the fluid from the &
                &largest pores. (Example silty-sand = 2310)'',/)') 
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(16)
          write(*,'(/'' The van Genuchten exponent(n) is a measure of th&
                &e pore-size distribution.'',/,'' (Example silty-sand = &
                &1.5)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(17)
          write(*,'(/'' Percent sand (0.0 - 1.0). (Example silty-sand = &
                &0.82)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(18)
          write(*,'(/'' Percent silt (0.0 - 1.0). (Example silty-sand = &
                &0.12)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(19)
          write(*,'(/'' Percent clay (0.0 - 1.0). (Example silty-sand = &
                &0.06)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(20)
          write(*,'(/'' Percent carbon (0.0 - 1.0). (Example silty-sand &
                &= 0.01)'',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(21)
          write(*,'(/'' Plastic Limit (0.0 - 200.0). (Example silty-sand&
                & = 28 - 50) '',/)')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        case(22)
          write(*,'(/'' Percent passing #200 sieve (0.0 - 100.0). '',/, &
                '' (Example silty-sand is greater than 50%)'',/)')

	end select

	return
	end
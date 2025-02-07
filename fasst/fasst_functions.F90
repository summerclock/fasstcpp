! functions used by fasst
! ******************************************************************************
      real(kind=8) function dense(temp1,wind1,f1)

      use fasst_global

! density (kg/m^3) as a function of temperature
! sources:
!   air & water vapor: http://users.wpi.edu/~icardi/PDF (temperature range 100K - 1600K and standard pressure)
!   water:             Hillel (1998) "Environmental soil Physics", Noborio et al. (1996) Soil Dci Soc Am J 60: 1010-1021
!   ice                http://www.engineeringtoolbox.com


      implicit none

      integer(kind=4),intent(in):: f1
      real(kind=8),intent(in):: temp1,wind1

! local variables
      real(kind=8):: t1,t2

      t1 = 0d0
      t2 = 0d0
      dense = 0d0

      t1 = temp1 - 273.15d0                                              !Celsius

      select case (f1)

        case(0)  !water vapor
          dense = 4.192d-12*temp1*temp1*temp1*temp1                     & 
                 - 1.25128d-8*temp1*temp1*temp1                         &
                 + 1.45079d-5*temp1*temp1 - 8.12253d-3*temp1 + 2.17634d0
          dense = dmin1(1.14d0,dmax1(6d-1,dense))

        case(1)  !water
          dense = -3d-07*t1*t1*t1*t1 + 7d-05*t1*t1*t1 - 9.92d-3*t1*t1   & 
                   + 8.666d-2*t1 + 999.81d0
          dense = dmin1(1d3,dmax1(962d0,dense))

        case(2)  !ice
          dense = -2d-10*t1*t1*t1*t1*t1*t1 - 7d-8*t1*t1*t1*t1*t1        &
                  - 1d-5*t1*t1*t1*t1 - 7d-4*t1*t1*t1                    &
                  - 2.37d-2*t1*t1 - 4.36d-1*t1 + 9.1612d2
          dense = dmin1(9.257d2,dmax1(9.162d2,dense))

        case(3)  !air
          if(temp1 > 0d0) then
            dense = 3.6077819d2*(temp1**(-1.00336d0))
          else
            dense = 1d0
          end if
          dense = dmin1(2.05d0,dmax1(0.948d0,dense))
!          if(temp1 > eps) rhoa = 3.48d-3*(wind1*1d2)/temp1               !kg/m^3  ?dry

        case(4)  !snow
          if(temp1 > 260.15d0.and.temp1 <= 276.65d0) then
            if(wind1 > 0d0) then
              t2 = 1.4d0*((277.15d0 - temp1)**(-1.15d0))                &
                                                  + 8d-3*(wind1**1.7d0)
            else
              t2 = 1.4d0*((277.15d0 - temp1)**(-1.15d0))
            end if
            if(t2 > 5d1) then
              dense = 1d0
            else
              dense = 1d0 - 0.951d0*dexp(-t2)
            end if
          else if(temp1 <= 260.15d0) then
            if(wind1 > 0d0) then
              t2 = 8d-3*(wind1**1.7d0)
            else
              t2 = 8d-3*(1d-1**1.7d0)
            end if
            if(t2 > 5d1) then
              dense = 1d0
            else
              dense = 1d0 - 0.904d0*dexp(-t2)
            end if
          else
            dense = 1d0
          end if

          dense = dmax1(sdensd,dmin1(sdensw,dense*1d3))  !7d2  !5d2

      end select

      dense = anint(dense*1d20)*1d-20

      end function dense

! ******************************************************************************
      real(kind=8) function dddt(temp1,f1)

! derivative of density (kg/m^3) w.r.t. temperature
! sources:
!   air & water vapor: http://users.wpi.edu/~icardi/PDF
!   water:             Hillel (1998) "Environmental soil Physics", Noborio et al. (1996) Soil Sci Soc Am J 60: 1010-1021
!   ice                http://www.engineeringtoolbox.com


      implicit none

      integer(kind=4),intent(in):: f1
      real(kind=8),intent(in):: temp1

! local variables
      real(kind=8):: t1

      t1 = 0d0
      dddt = 0d0

      t1 = temp1 - 273.15d0                                              !Celsius

      select case (f1)

        case(0)  !water vapor
          dddt = 4d0*(4.192d-12*temp1*temp1*temp1)                      & 
                 - 3d0*(1.25128d-8*temp1*temp1)                         &
                 + 2d0*(1.45079d-5*temp1) - 8.12253d-3

        case(1)  !water
          dddt = 4d0*(-3d-07*t1*t1*t1) + 3d0*(7d-05*t1*t1)              &
                   - 2d0*(9.92d-3*t1) + 8.666d-2

        case(2)  !ice
          dddt = 6d0*(-2d-10*t1*t1*t1*t1*t1) - 5d0*(7d-8*t1*t1*t1*t1)   &
                 - 4d0*(1d-5*t1*t1*t1) - 3d0*(7d-4*t1*t1)               &
                 - 2d0*(2.37d-2*t1) - 4.36d-1

        case(3)  !air
          if(temp1 > 0d0) then
            dddt = -1.00336d0*(3.6077819d2*(temp1**(-2.00336d0)))
          else
            dddt = -1.00336d0
          end if

      end select

      dddt = anint(dddt*1d20)*1d-20

      end function dddt

! ******************************************************************************
      real(kind=8) function spheats(temp1,f1)

! specific heat (J/kg*K) as a function of temperature
! sources: 
!   air & water vapor: http://users.wpi.edu/~ierardi/PDF
!   water:             http://www.engineeringtoolbox.com
!   ice               Jordan (1991) CRREL SR 91-16, p.17


      implicit none

      integer(kind=4),intent(in):: f1
      real(kind=8),intent(in):: temp1

! local variables
      real(kind=8):: t1

      t1 = 0d0
      spheats = 0d0

      t1 = temp1 - 273.15d0                                              !Celsius

      select case(f1)

        case(0)  !water vapor
          spheats = 2.3888d-8*temp1*temp1*temp1*temp1                   &
                     - 6.5129d-5*temp1*temp1*temp1                      &
                     + 6.6178d-2*temp1*temp1 - 2.9086d1*temp1 + 6.6256d3
          spheats = dmin1(3.26d3,dmax1(2d3,spheats))

        case(1)  !water
          spheats = -1d-9*t1*t1*t1*t1*t1*t1 + 4d-7*t1*t1*t1*t1*t1       &
                     - 4d-5*t1*t1*t1*t1 +1.6d-3*t1*t1*t1 +4.5d-3*t1*t1  &
                     -1.8731d0*t1 + 4210d0
          spheats = dmin1(4219d0,dmax1(4178d0,spheats))

        case(2)  !ice
          spheats = -13.3d0 + 7.8d0*temp1
          if(spheats > 2050d0) spheats = 2050d0
          if(spheats < 1389d0) spheats = 1389d0

        case(3)  !dry air
          spheats = 1.9327d-10*temp1*temp1*temp1*temp1                  & 
                     - 7.9999d-7*temp1*temp1*temp1                      &
                     + 1.1407d-3*temp1*temp1 - 4.489d-1*temp1 + 1.0575d3
          spheats = dmin1(1.25d3,dmax1(1d3,spheats))

      end select

      spheats = anint(spheats*1d20)*1d-20

      end function spheats

! ******************************************************************************
      real(kind=8) function thconds(temp1,f1)

! thermal conductivity (W/m*K)
! sources:
!   air & water vapor: http://users.wpi.edu/~icardi/PDF
!   water:             Farouki (1981)
!   ice                http://www.engineeringtoolbox.com


      implicit none

      integer(kind=4),intent(in):: f1
      real(kind=8),intent(in):: temp1

! local variables
      real(kind=8):: t1


      t1 = 0d0
      thconds = 0d0

      t1 = temp1 - 273.15d0                                              !Celsius

      select case(f1)

        case(0)  !water vapor
          thconds = 8.3154d-5*temp1 - 7.4556d-3
          thconds = dmin1(2.36d-2,dmax1(6.95d-3,thconds))

        case(1)  !water
          thconds = 1.8d-3*temp1 + 0.0787d0
          thconds = dmin1(7.5d-1,dmax1(5.8d-1,thconds))

        case(2)  !ice
          thconds = 4d-7*t1*t1*t1 + 1d-4*t1*t1 - 6.9d-3*t1 + 2.2174d0
!         if(dabs(temp1) > eps) thconds = 488.19d0/temp1 + 0.4685d0      !? source, came from thcond.f
          thconds = dmin1(3.48d0,dmax1(2.2174d0,thconds))

        case(3)  !dry air
          thconds = 1.5207d-11*temp1*temp1*temp1                        &
                     - 4.8574d-8*temp1*temp1                            &
                     + 1.0184d-4*temp1 - 3.9333d-4
!         ka = 2.37d-2 + 6.41d-5*t1   !? source, came from thcond.f
          thconds = dmin1(1d-1,dmax1(1.59d-2,thconds))

      end select

      thconds = anint(thconds*1d20)*1d-20

      end function thconds

! ******************************************************************************
      real(kind=8) function met_date(year,doy,hr,minute)

! calculate the decimal calendar date


      implicit none

      real(kind=8),intent(in):: year,doy,hr,minute

! local variables
      real(kind=8):: mody


      mody = 0d0
      met_date = 0d0

      mody = year - aint(year*2.5d-1)*4d0

      if(dabs(mody) <= epsilon(1d0)) then
        met_date = year + doy/367d0 + hr/(24d0*367d0)                   &
                                              + minute/(6d1*24d0*367d0)
      else
        met_date = year + doy/366d0 + hr/(24d0*366d0)                   & 
                                              + minute/(6d1*24d0*366d0)
      end if
      met_date = anint(met_date*1d10)*1d-10

      end function met_date

! ******************************************************************************
      real(kind=8) function head(i,smt)

      use fasst_global

! calculate the pressure head based on soil water content using Van Genuchten (1980)


      implicit none

      integer(kind=4),intent(in):: i
      real(kind=8),intent(in):: smt

! local variables
      real(kind=8):: w1,e1,e2


      w1 = 0d0
      e1 = 0d0
      e2 = 0d0
      head = 0d0

      if((nsoilp(i,11) > eps.and.nsoilp(i,12) > eps).and.               &
                                               nsoilp(i,10) > eps) then
        if(smt > nsoilp(i,15).and.smt < nsoilp(i,24)) then
          w1 = (smt - nsoilp(i,8))/(nsoilp(i,9) - nsoilp(i,8))           !unitless
          e1 = -1d0/nsoilp(i,12)
          e2 = 1d0/nsoilp(i,11)
          if(w1 > eps.and.dabs((w1**e1) - 1d0) > eps)                   &
            head = (-(1d0/nsoilp(i,10))*((w1**e1) - 1d0)**e2)*1d-2
        else if(smt <= nsoilp(i,15)) then
          head = pheadmin(i)
        else if(smt >= nsoilp(i,24)) then
          head = 0d0
        end if
      else
        head = pheadmin(i)
      end if

      head = anint(head*1d20)*1d-20                                      !m

      end function head

! ******************************************************************************
      real(kind=8) function soilhumid(i,ph,sms,st)

      use fasst_global

! calculate the soil relative humidity using Campbell (1985)

      implicit none

      integer(kind=4),intent(in):: i
      real(kind=8):: ph,sms,st

! local variables
      real(kind=8):: t1


      t1 = 0d0
      soilhumid = 0d0

! note: shead <= 0d0 for unsaturated soil
      if(ntype(i) < 20) then
        if(ph >= eps) then
          soilhumid = 1d0
        else
          if(dabs(st) > eps) t1 = grav*ph/(Rv*st)

          if(dabs(t1) > 7d0.or.sms <= eps) then
            soilhumid = sms/nsoilp(i,9)
          else
!            soilhumid = dmax1(sm/nsoilp(i,9),dexp(t1))
            soilhumid = dexp(t1)
          end if
        end if
      else
        if(ice(i) <= eps) then
          soilhumid = 5d0*sms/nsoilp(i,9)
        else
          soilhumid = sms/nsoilp(i,9)
        end if
      end if

      soilhumid = dmin1(1d0,dmax1(0d0,soilhumid))

      soilhumid = anint(soilhumid*1d20)*1d-20

      end function soilhumid

! ******************************************************************************
      real(kind=8) function vap_press(i,rh,ap)

      use fasst_global

! calculates the vapor pressure (Pa)

      implicit none

      integer(kind=4):: i
      real(kind=8):: rh,ap

! local variables
      real(kind=8):: a,b,t1

      a = 0d0
      b = 0d0
      t1 = 0d0
      vap_press = 0d0

      if((ice(i) <= 0d0.and.i <= nnodes).or.                            &
                                       (i > nnodes.and.hm <= eps)) then  !over water; no snow
        a = 17.269d0
        b = 35.86d0
      else                                                               !over ice/snow
        a = 21.8745d0
        b = 7.66d0
      end if

      t1 = stt(i) - b                                                    !K
      if(dabs((a*(stt(i) - Tref))/t1) > 5d1.or.dabs(t1) <= eps) then
        vap_press = ap*1d2*rh                                            !Pa
      else
        vap_press = 610.78d0*rh*dexp(a*(stt(i) - Tref)/t1)               !Pa
      end if

      vap_press = anint(vap_press*1d20)*1d-20

      end function vap_press

! ******************************************************************************
      real(kind=8) function maxinfiltrate(i,smi)

      use fasst_global


      implicit none

      integer(kind=4),intent(in):: i
      real(kind=8),intent(in):: smi 

! local variables
      integer(kind=4):: j,n
      real(kind=8):: s2,w,theta,integral,e1,e2,dw,wold
      real(kind=8):: i1,i2,i3


      s2 = 0d0
      w = 0d0
      theta = 0d0
      integral = 0d0
      e1 = 0d0
      e2 = 0d0
      dw = 0d0
      wold = 0d0
      maxinfiltrate = 0d0

      if(nsoilp(i,12) > eps) e1 = 1d0/nsoilp(i,12)
      if(nsoilp(i,11) > eps) e2 = 1d0/nsoilp(i,11)

      j = 1
      n = 150
      do while(j < n)
        theta = smi + (nsoilp(i,9) - smi)*real(j)/real(n)
        w = (theta - nsoilp(i,8))/(nsoilp(i,9) - nsoilp(i,8))
        w = dmin1(dmax1(eps,w),1d0)

        dw = w - wold
        if(w > eps.and.e1 > eps) then
          if(dabs(-e1-5d-1) > eps) then
            i1 = w**(-e1 - 5d-1)
          else
            i1 = 0d0
          end if

          if((nsoilp(i,12) > eps.and.dabs(1d0 - w**e1) > eps).and.      &
                    dabs(1d0 - (1d0 - w**e1)**nsoilp(i,12)) > eps) then
            i2 = (1d0 - (1d0 - w**e1)**nsoilp(i,12))**2d0
          else
            i2 =  1d0 
          end if

          if(dabs(w**(-e1) - 1d0) > eps.and.dabs(e2-1d0) > eps) then
            i3 = (w**(-e1) - 1d0)**(e2 - 1d0)
          else
            i3 = 0d0
          end if

        else
          i1 = 1d0
          i2 = 1d0
          i3 = 0d0
        end if

        integral = i1*i2*i3*dw
        wold = w
        s2 = s2 + integral

        j = j + 1
      end do

      integral = s2
      if(nsoilp(i,10)*nsoilp(i,11)*nsoilp(i,12) > eps)                 &
      s2 = 2d0*nsoilp(i,7)*integral/                                   &
                               (nsoilp(i,10)*nsoilp(i,11)*nsoilp(i,12))  !cm^2/s

      maxinfiltrate = anint((s2*1d-4)*1d20)*1d-20                        !m^2/s

      end function maxinfiltrate

! ******************************************************************************
      integer(kind=4) function map_USDA_SoilType_to_USCS(lis)


      implicit none

      integer(kind=4),intent(in):: lis


! LIS/STATSGO soil types:
!     1 = sand                2 = loamy sand          3 = sandy loam
!     4 = silt loam           5 = silt                6 = loam
!     7 = sandy clay loam     8 = silty clay loam     9 = clay loam
!    10 = sandy clay         11 = silty clay         12 = clay
!    13 = peat               14 = open water         15 = bedrock

! FASST/USCS soil types:
! SEDRIS EDCS_AC_SOIL_TYPES (stype(maxl))
!   0 = unknown     5 = SW     9 = ML     12 = CH       15 = PT
!   1 = GW          6 = SP    10 = CL     13 = MH       16 = MC (SMSC) nonSEDRIS
!   2 = GP          7 = SM    11 = OL     14 = OH       17 = CM (CLML)
!   3 = GM          8 = SC                              18 = EVaporites (not used)
!   4 = GC
!  20 = COncrete   21 = ASphalt  Note: both of these are nonSEDRIS
!  25 = ROck       30 = SNow     Note: both of these are nonSEDRIS
!  26 = WAter                    Note: this is nonSEDRIS
!  27 = AIr                      Note: this is nonSEDRIS

! local variables
      integer(kind=4):: convert(15)

!                  1,2,3,4,5, 6,7, 8, 9,10,11,12,13,14,15
      data convert/6,7,7,9,9, 9,8,10,10, 8,12,12,15,26,25/  !Me
!      data convert/6,7,8,7,9,10,8,10,10, 8,10,12,14,26, 3/  !GSL

      if(lis <= 0) then
        map_USDA_SoilType_to_USCS = lis
      else
        map_USDA_SoilType_to_USCS = convert(lis)
      end if

      end function map_USDA_SoilType_to_USCS

! ******************************************************************************
      integer(kind=4) function USCS_sandsiltclay(sand,silt,clay,p200,   &
                                                 plimit,llimit,orgf)


      implicit none

      real(kind=8),intent(in):: sand,silt,clay,p200,plimit,llimit,orgf

! local variables
      real(kind=8):: gravel,aline,pli,eps


      eps = 0d0
      eps = epsilon(1d0)                                                 !tolerance limit for equality

      gravel = 0d0
      if(1d0-(sand+silt+clay) > eps) gravel = 1d0 - (sand + silt + clay)

      pli = 0d0
      aline = 0d0
      pli = plimit + llimit
      aline = 7.3d-1*(llimit - 2d1) 

      if(orgf < 3d1) then
        if(llimit < 5d1) then
          if(orgf < 2.5d1) then
            if(pli > 7d0.and.pli >= aline) then
              USCS_sandsiltclay = 10
            else if((pli >= 4d0.and.pli <= 7d0).and.pli >= aline) then
              USCS_sandsiltclay = 17
            else if(pli < 4d0.or.pli < aline) then
              USCS_sandsiltclay = 9
            end if
          else if(orgf >= 2.5d1) then
            USCS_sandsiltclay = 11
          end if
        else
          if(orgf < 2.5d1) then
            if(pli >= aline) then
              USCS_sandsiltclay = 12
            else if(pli < aline) then
              USCS_sandsiltclay = 13
            end if
          else
            USCS_sandsiltclay = 14
          end if
        end if
      else
        USCS_sandsiltclay = 15
      end if

      if(gravel > eps) then
        if(sand < gravel) then                                           !gravels
          if(p200 < 5d0) then
            USCS_sandsiltclay = 2
          else if(p200 >= 5d0.and.p200 <= 12d0) then
            if(USCS_sandsiltclay == 9.or.USCS_sandsiltclay == 13) then
              USCS_sandsiltclay = 3
            else
              USCS_sandsiltclay = 4
            end if
          else
            if(USCS_sandsiltclay == 9.or.USCS_sandsiltclay == 13) then
              USCS_sandsiltclay = 3
            else if(USCS_sandsiltclay == 10.or.                         &
                                          USCS_sandsiltclay == 12) then
              USCS_sandsiltclay = 4
            else if(USCS_sandsiltclay == 17) then
              USCS_sandsiltclay = 4
            end if
          end if
        else if(sand <= gravel) then                                     !sands
          if(p200 < 5d0) then
            USCS_sandsiltclay = 6
          else if(p200 >= 5d0.and.p200 <= 12d0) then
            if(USCS_sandsiltclay == 9.or.USCS_sandsiltclay == 13) then
              USCS_sandsiltclay = 7
            else
              USCS_sandsiltclay = 8
            end if
          else
            if(USCS_sandsiltclay == 9.or.USCS_sandsiltclay == 13) then
              USCS_sandsiltclay = 7
            else if(USCS_sandsiltclay == 10.or.                         &
                                          USCS_sandsiltclay == 12) then
              USCS_sandsiltclay = 8
            else if(USCS_sandsiltclay == 17) then
              USCS_sandsiltclay = 16
            end if
          end if
        end if
      end if

      end function USCS_sandsiltclay

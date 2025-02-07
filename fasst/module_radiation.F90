module module_radiation

      use fasst_global


   contains

      subroutine dnirflx(ematm,ta,doy,lcldbse,mcldbse,hcldbse,lcld,     &
                         mcld,hcld,flxir)

! calls the following subroutines:
!     cloudbase (appended to this subroutine)

!     The following routine calculates the downwelling IR flux.  In this
!     routine it is assumed that the major contribution from the clear
!     flux is emitted by the atmosphere below the cloud layers.  Thus,
!     the clear air flux is NOT weighted by (1-cloud cover).  If the
!     cloud amount is missing a global cloud mean of 0.53 is assumed.  In
!     this case the entire cloud amount is assigned to the low cloud type.
!     If low cloud is not missing but middle and high cloud amounts are
!     the middle and high cloud amount is set to 0.  The cloud base, if
!     missing is set using climatological information. The climatological cloud
!     base is a function of season and latitude and is based on the research work 
!     done by G Koenig in his investigation of albedo vs. greenhouse effects of
!     cloud on climate

!     parameter                 remarks
!     a,b,c,d            Parameterization coef for cloud base calculation
!     zlcld,zmcld,zhcld  Cloud base hgts

      implicit none

      real(kind=8),intent(in):: ematm,ta,doy,lcldbse,mcldbse,hcldbse
      real(kind=8),intent(in):: lcld
      real(kind=8),intent(inout):: mcld,hcld
      real(kind=8),intent(out):: flxir

! local variables
      real(kind=8):: zlcld,zmcld,zhcld,flxclr,flxcld


      zlcld = 0d0
      zmcld = 0d0
      zhcld = 0d0
      flxclr = 0d0
      flxcld = 0d0
      flxir = 0d0

! calculate the cloud base heights
      call cloudbase(doy,lat,lcldbse,mcldbse,hcldbse,lcld,mcld,hcld,    &
                     zlcld,zmcld,zhcld)

! calculate the effective middle and high cloud amounts assuming random overlap

      hcld = hcld*(1d0 - mcld)*(1d0 - lcld)
      mcld = mcld*(1d0 - lcld)

! Calculate the clear air flux and cloud flux
! Replace the cloud flux parameterization with stephen Boltzman approach using 
! a climatologically atmospheric profile based on season and latitude. The profile
! should be adjusted by comparing the climate temperature at the surface with the
! measured temperature. 
      if(ta+Tref > eps) flxclr = ematm*sigma*(ta + Tref)**4d0
      flxcld = lcld*(94d0 - 5.8d0*zlcld) + mcld*(94d0 - 5.8d0*zmcld)    &
                 + hcld*(94d0 - 5.8d0*zhcld)

!     Total downwelling flux
      flxir = flxclr + flxcld

      end subroutine dnirflx


! ******************************************************************************
      subroutine cloudbase(doy,plat,lcldbse,mcldbse,hcldbse,lcld,mcld,  &
                           hcld,zlcld,zmcld,zhcld)

      implicit none

! this subroutine calls no subroutines

      real(kind=8),intent(in):: doy,plat,lcldbse,mcldbse,hcldbse
      real(kind=8),intent(in):: lcld,mcld,hcld
      real(kind=8),intent(out):: zlcld,zmcld,zhcld

! local variables
      integer(kind=4):: i,j,isean,ilat
      real(kind=8):: coef(2,2,3,4),a,b,c,d

      data ((coef(1,1,i,j),j=1,4),i=1,3)                                 &
      /1.05d0,0.6d0,5d0,25d0,                                            &
       4.1d0 ,0.3d0,4d0,25d0,                                            & 
       7d0   ,1.5d0,3d0,30d0/
      data ((coef(1,2,i,j),j=1,4),i=1,3)                                 & 
      /1.05d0,0.6d0,1.5d0,25d0,                                          &
       4.1d0 ,2d0  ,1.7d0,25d0,                                          & 
       7d0   ,1.5d0,3d0  ,30d0/
      data ((coef(2,1,i,j),j=1,4),i=1,3)                                 &
      /1.15d0,0.45d0,5d0,25d0,                                           & 
       4.4d0 ,0.3d0 ,4d0,25d0,                                           &
       7d0   ,1.5d0 ,3d0,30d0/
      data ((coef(2,2,i,j),j=1,4),i=1,3)                                 &
      /1.15d0,0.6d0,1.5d0,25d0,                                          & 
       4.4d0 ,1.2d0,3d0  ,25d0,                                          &
       7d0   ,1.5d0,3d0  ,30d0/


      zlcld = 0d0
      zmcld = 0d0
      zhcld = 0d0

      ilat = 1
      if(plat >= dabs(25d0)) ilat = 2
      isean = 2                                                          !not winter
      if(plat > 0d0.and.(int(doy) > 330.or.int(doy) < 65))              &
        isean = 1                                                        !winter, northern hemisphere
      if(plat < 0d0.and.(int(doy) > 150.or.int(doy) < 250))             &
        isean = 1                                                        !winter, southern hemisphere

! Set cloud base altitude if missing

      if(aint(dabs(lcldbse-mflag)*1d5)*1d-5 <= eps.and.lcld > eps) then
        a = coef(isean,ilat,1,1)
        b = coef(isean,ilat,1,2)
        c = coef(isean,ilat,1,3)
        d = coef(isean,ilat,1,4)
        zlcld = a - b*(1d0 - dabs(dcos(c*(lat - d))))
      else
        zlcld = lcldbse
      end if
      if(zlcld <= eps) zlcld = 0d0

      if(aint(dabs(mcldbse-mflag)*1d5)*1d-5 <= eps.and.mcld > eps) then
        a = coef(isean,ilat,2,1)
        b = coef(isean,ilat,2,2)
        c = coef(isean,ilat,2,3)
        d = coef(isean,ilat,2,4)
        zmcld = a - b*(1d0 - dabs(dcos(c*(lat - d))))
      else
        zmcld = mcldbse
      end if
      if(zmcld <= eps) zmcld = 0d0

      if(aint(dabs(hcldbse-mflag)*1d5)*1d-5 <= eps.and.hcld > eps) then
        a = coef(isean,ilat,3,1)
        b = coef(isean,ilat,3,2)
        c = coef(isean,ilat,3,3)
        d = coef(isean,ilat,3,4)
        zhcld = a - b*(1d0 - dabs(dcos(c*(lat - d))))
      else
        zhcld = hcldbse
      end if
      if(zhcld <= eps) zhcld = 0d0


      end subroutine cloudbase

! ******************************************************************************
      subroutine emisatm(tmp,rh,ematm)

! this subroutine calls no subroutines

!     The longwave atmospheric emissivity is given as
!         ematm=-0.792+3.161*ematmp-1.573*ematmp^2
!     where
!         ematmp=-0.7+5.95*10**-5*vp*exp(1500/ta)
!     where vp is the vapor pressure in mbs and is obtained from
!     the relative humidity and ambient temperature using the Clausius
!     Clapeyrin equation.

! A new atmospheric emissivity formulation has been implemented based on 
! Todd M. Crawford & Claude E. Duchon, Am Improved parameterization for estimating
! effective atmospheric emissivity for use in calculating daytime downwelling longwave
! radiation, Jour of Applied meteorology, vol 18 april 1999, 474-480

!         ematm=1.24*[ea/Ta(K)]^(1/7)

!     parameter                 remarks
!         l                   latent heat of condensation (J/kg)
!         rv                  gas constant (j/(kg-K))
!        eso                  saturated vapor pressure (mb) at standard 
!                             conditions(1013 mb, 273.15K)
!        esa                  saturated vapor pressure (mb) at ta
!        ea                   vapor pressure
!      ematm                  atmospheric emissivity

      implicit none

      real(kind=8),intent(in):: tmp,rh
      real(kind=8),intent(out):: ematm
      
! local variables
      real(kind=8):: l,rv,eso,ea,ta


      l = 0d0
      rv = 0d0
      eso = 0d0
      ea = 0d0
      ta = 0d0
      ematm = 0d0

      rv = 461d0
      eso = 6.13d0

! retrieve the temp for this hour
      ta = tmp + Tref

! Calculate the vapor pressure based on the RH(in %)
      l = (-2.43d-3*ta + 3.166659d0)*10d0**6d0
      if(dabs(l/rv*(1d0/Tref - 1d0/ta)) < 5d1)                          &
        ea = (eso*dexp(l/rv*(1d0/Tref - 1d0/ta)))*rh*1d-2

!     Calculate the emissivity
!     ematmp=0.7+5.95E-05*ea*exp(1500/ta)                                 !idso 
!     ematm=-0.792+3.161*ematmp-1.573*ematmp*ematmp                       !idso with watchmann corr
      if(ea > 0d0) ematm = 1.24d0*(ea/ta)**(1d0/7d0)                      !special test of brutsaert

      end subroutine emisatm

! ******************************************************************************
      subroutine sol_zen(year,doy,hour,minute,szen,saz)

! this subroutine calls no subroutines

! calculate the solar zenith and azimuth angles. This routine is based on the NOAA
! web site (www.srrb.noaa.gov\highlights\sunrise\azel.html) which is based on the
! work of Meeus, Jean "Astronomical Algorithms"

      implicit none

      real(kind=8),intent(in):: year,doy,hour,minute
      real(kind=8),intent(out):: szen,saz

! local variables
      real(kind=8):: fyear,eqtime,decl,time_offset,tst,ha,zlat,zlong
      real(kind=8):: cosphi,phi,costheta,thour,f1,f2,mody


      fyear = 0d0
      eqtime = 0d0
      decl = 0d0
      time_offset = 0d0
      tst = 0d0
      ha = 0d0
      zlat = 0d0
      zlong = 0d0
      cosphi = 0d0
      phi = 0d0
      costheta = 0d0
      thour = 0d0
      f1 = 0d0
      f2 = 0d0
      mody = 0d0
      szen = 0d0
      saz = 0d0

      f1 = 1d0/180d0
      f2 = 1d0/pi
      zlat = pi*lat*f1                                                   !radians
      zlong = -pi*mlong*f1                                               !radians
      thour = hour
      if(dabs(thour-24d0) <= eps) thour = 0d0

! calculate the fractional year in radians
      mody = year - aint(year*2.5d-1)*4d0
      if(dabs(mody) > eps) then
        fyear = (2d0*pi/365d0)*(doy - 1d0 + (thour - 12d0)/24d0)
      else
        fyear = (2d0*pi/366d0)*(doy - 1d0 + (thour - 12d0)/24d0)
      end if

! calculate the equation of time in minutes
      eqtime = 229.18d0*(7.5d-5 + 1.868d-3*dcos(fyear)                  & 
                - 3.2077d-2*dsin(fyear) - 1.4615d-2*dcos(2d0*fyear)     &
                - 4.0849d-2*dsin(2d0*fyear))

! calculate the solar declination in radians
      decl = 6.918d-3 - 3.99912d-1*dcos(fyear) + 7.0257d-2*dsin(fyear)  &
                - 6.758d-3*dcos(2d0*fyear) + 9.07d-4*dsin(2d0*fyear)    &
                - 2.697d-3*dcos(3d0*fyear) + 1.48d-3*dsin(3d0*fyear)

! calculate the true solar time
      time_offset = eqtime + 4d0*mlong + 60d0*timeoffset                 !minutes

      tst = thour*60d0 + minute + time_offset                            !minutes

! calculate the solar hour angle
      ha = tst*2.5d-1 - 180d0                                            !degrees
      ha = ha*pi*f1                                                      !radians

! calculate the solar zenith and azimuth angles in degrees
      cosphi = dsin(zlat)*dsin(decl) + dcos(zlat)*dcos(decl)*dcos(ha)
      phi = acos(cosphi)
      szen = phi*180d0*f2

      costheta = (dsin(zlat)*cosphi-dsin(decl))/(dcos(zlat)*dsin(phi))
      if(zlat >= 0d0)then
        if(ha < 0d0) then
          saz = 180d0 - acos(costheta)*180d0*f2
        else
          saz = 180d0 + acos(costheta)*180d0*f2       
        end if
      else
        if(ha < 0d0) then
          saz = acos(costheta)*180d0*f2
        else
          saz = 360d0 - acos(costheta)*180d0*f2       
        end if
      end if
            
      end subroutine sol_zen

! ******************************************************************************
      subroutine Solflx(icld,szen,doy,prcp,tsol,cover,hgt,sdir,sdif)

! calls the following subroutines:
!     insol (appended to this subroutine)

!     This program implements Shapiro's parameterization of the
!     solar flux as coded in soil thermal model
!     jday: Julian day
!     cover: Fractional cloud cover for low, middle and high clouds
!     icl: Cloud type
!          Low: 0=none, 4=Sc or St, 5=Cu or Cb
!          Middle: 0=none, 3=As, Ac or any other type
!          High: 0=none, 1=Ci, 2=Cs
!     cover: fractional cloud amount for low, middle, and high
!     prcp: precipitation in meters/hour
!     albedo: surface albedo


      implicit none

      integer(kind=4),intent(inout):: icld(3)
      real(kind=8),intent(in):: szen,doy,prcp
      real(kind=8),intent(inout):: tsol,cover(3),hgt(3)
      real(kind=8),intent(out):: sdir,sdif

! local variables
      integer(kind=4):: ctest,d1
      real(kind=8):: cosz,frad,test,prcp1,sdowns,direct
      real(kind=8):: diffuse,f1


      ctest = 0
      d1 = 0
      cosz = 0d0
      frad = 0d0
      test = 0d0
      prcp1 = 0d0
      sdowns = 0d0
      direct = 0d0
      sdir = 0d0
      sdif = 0d0
      f1 = 0d0

      frad = 0.017453292d0                                                !pi/180


!      if(icld(1) /= 0) icld(1) = 4                                       !geoff's value was 4
!      if(icld(2) /= 0) icld(2) = 3                                       !geoff's value was 3

!      if(icld(3) /= 0) then                                              !geoff's value was 1 (0-0.4) or 2 (0.4-1)
!        if(cover(3) > 0d0.or.icld(3) == 5) icld(3) = 1
!        if(cover(3) > 0.4d0.or.icld(3) == 7) icld(3) = 2
!      end if

      if((cover(1) > 0.5d0).or.(cover(2) > 0.5d0).or.                   &
        (cover(3) > 0.5d0)) ctest = 1
      prcp1 = prcp*1d-3                                                  !m/timestep

      cosz = dcos(szen*frad)                                             !szen is the solar zenth angle

! Calculate solar flux and store in the met array
      d1 = int(doy)
      call insol(d1,icld,cosz,prcp1,hgt,cover,sdowns,direct,diffuse)

      if(aint(dabs(tsol-mflag)*1d5)*1d-5 <= eps) tsol = sdowns
      test = direct + diffuse
      if(sdowns > eps) f1 = tsol/sdowns

      if(dabs(test-tsol) > 1d-3) then
        if(ctest == 1) then
          direct = direct*f1
          diffuse = dmax1(0d0,tsol - direct)
          if(diffuse < eps) direct = tsol
        else
          diffuse = diffuse*f1
          direct = dmax1(0d0,tsol - diffuse)
          if(direct < 0d0) diffuse = tsol
        end if
      end if
      sdir = direct
      sdif = diffuse

      end subroutine Solflx

! ******************************************************************************
      subroutine insol(jday,icld,cosz,prcp1,hgt,cover,sdowns,direct,    &
                       diffuse)

!  Subroutine INSOL computes direct and diffuse solar radiation using
!  method of R. Shapiro 

! This program calls no subroutines

! Arguments
! jday : day from met data
! cover(3) : Fractional cloud cover
! icl(3) : Cloud type code
! prcp : Precipitation Value (m/hour)
! albedo : Albedo of surface node 
! cosz: cos of zenith angle.
! sdown : Incident solar radiative flux [W/m^2]
! sdowne: Solar insolation corrected for Earth's elliptical orbit.
! direct : Direct component of insolation.
! diffuse: diffuse component of solar radiation sdown.

! Local
! cosz: cos of zenith angle.
! coszcube: cube of cosz.
! coszsq: square of cosz.
! d1: d1=1d0-(rks(1)*rks(2)).
! d2: d2=1d0-(rks(2)*rks(3)).
! d3: d3=1d0-(rks(3)*rg).
! fr: fr=cover(i).
! i: looping index.
! j: j=icld(i).
! l: array element.
! pi: pi=3.14159265.
! rcld(4): Reflectivity for cloudy skies
!       =r2(j,1)+r2(j,2)*cosz+r2(j,3)*coszsq+r2(j,4)*coszcube.
! rclr(4): Reflectivity for clear skies
!       =r1(l,1)+r1(l,2)*cosz+r1(l,3)*coszsq+r1(l,4)*coszcube.
! rks(3): Combined reflectivity of clear and cloudy skies.
! rg: Albedo of ground (rg=albedo).
! r1(4,4): reflectance factor array for clear part of sky.
! r2(4,4): reflectance factor array for clouded part of sky.
! sdown0: sdown0=1369.2d0*sdowne*cosz.
! tcld(4): Transmissivity for cloudy skies
!          =t2(j,1)+t2(j,2)*cosz+t2(j,3)*coszsq+t2(j,4)*coszcube.
! tclr(4): Transmissivity for clear skies
!          =t1(l,1)+t1(l,2)*cosz+t1(l,3)*coszsq+t1(l,4)*coszcube.
! tdk(3): Direct componenet of tk(3).
! tk(3): Combined transmissivity of clear and cloudy skies.
! t1(4,4): transmittance factor array for clear part of sky.
! t2(4,4): transmittance factor array for clouded part of sky.
! wgt: Total weighted cloud fraction.
! wt(4,6): cloud cover factor array.
! icla(3): cloud type code (icla(1)=icl(3),icla(2)=icl(2),
!          icla(3)=icl(1))
! covera(3): cloud cover fraction (covera(1)=cover(3),covera(2)=cover(2),
!          covera(3)=cover(1))


      implicit none

      integer(kind=4),intent(in):: jday
      integer(kind=4),dimension(3),intent(inout):: icld
      real(kind=8),intent(in):: cosz,prcp1
      real(kind=8),dimension(3),intent(inout):: cover,hgt
      real(kind=8),intent(out):: sdowns,direct,diffuse

! local variables
      integer(kind=4):: i,j,l,icla(3)
      real(kind=8):: sdown0,wgt,rcld,tcld,fr,tclr
      real(kind=8):: coszsq,coszcube,d3,d2,rg,d1,rclr,sdowne
      real(kind=8):: rks(3),tk(3),tdk(3),covera(3)


! This next section was formerly blockdata dinsol
! ******************************************************************************
      real(kind=8):: r1(4,4),r2(4,4),t1(4,4),t2(4,4),wt(4,6)

! r1(4,4): = Cubic polynomial coefficients for Reflectivity
!            for clear sky
! r2(4,4): = Cubic relectivity coeffients for Reflectivity
!            for cloudy skies.
! t1(4,4): = Cubic polynomial coefficients for Transmissivty
!            for clear skies
! t2(4,4): = Cubic polynomial coefficients for Transmissivty
!            for cloudy skies.
! wt(4,6): = Biquadratic polynomial coefficients for the Weights.

! rclr=a0 + a1*cosz + a2*cosz^2 + a3*cosz^3, where in the
! following array the rows are a0, a1, a2 and a3, and
! the columns correspond to
!  1. Highest layer
!  2. Middle layer
!  3. Lowest layer
!  4. Smoke/fog (Not used in this model)

      data r1/0.12395d0,0.15325d0,0.15946d0,0.27436d0,                  &
             -0.34765d0,-0.3962d0,-0.42185d0,-0.43132d0,                &
              0.39478d0,0.42095d0,0.488d0,0.2692d0,                     &
             -0.14627d0,-0.142d0,-0.18492d0,-0.00447d0/

! *******************************
! rcld=a0 + a1*cosz + a2*cosz^2 + a3*cosz^3, where in the
! following array the rows are a0, a1, a2 and a3, and
! the columns correspond to 
!  1. Thin Ci/Cs (denoted as cloud type 1)
!  2. Thick Ci/Cs (denoted as cloud type 2)
!  3. As/Ac (denoted as cloud type 3)
!  4. Low cloud (denoted as cloud type 4 or 5)

      data r2/0.25674d0,0.42111d0,0.61394d0,0.69143d0,                  &
             -0.18077d0,-0.04002d0,-0.01469d0,-0.14419d0,               &
             -0.21961d0,-0.51833d0,-0.174d0,-0.051d0,                   &
              0.25272d0,0.4054d0,0.14215d0,0.06682d0/

! *******************************
! tclr=b0 + b1*cosz + b2*cosz^2 + b3*cosz^3, where in the
! following array the rows are b0, b1, b2 and b3, and
! the columns corespond to
!  1. Highest layer
!  2. Middle layer
!  3. Lowest layer
!  4. Smoke/fog (Not used in this model)

      data t1/0.76977d0,0.69318d0,0.68679d0,0.55336d0,                  &
              0.49407d0,0.68227d0,0.71012d0,0.61511d0,                  &
             -0.44647d0,-0.64289d0,-0.71463d0,-0.29816d0,               &
              0.11558d0,0.1791d0,0.22339d0,-0.06663d0/

! *******************************
! tcld=b0 + b1*cosz + b2*cosz^2 + b3*cosz^3, where in the
! following array the rows are b0, b1, b2 and b3, and
! the columns correspond to
!  1. Thin Ci/Cs (denoted as cloud type 1)
!  2. Thick Ci/Cs (denoted as cloud type 2)
!  3. As/Ac (denoted as cloud type 3)
!  4. Low cloud (denoted as cloud type 4 or 5)

      data  t2/0.63547d0,0.43562d0,0.23865d0,0.15785d0,                 &
               0.35229d0,0.26094d0,0.20143d0,0.3241d0,                  &
               0.08709d0,0.36428d0,-0.01183d0,-0.14458d0,               &
              -0.22902d0,-0.38556d0,-0.07892d0,0.01457d0/

! *******************************
! W=c0 + c1*cosz + c2*fk + c3*fk*cosz + c4*cosz^2 + c5*fk^2, 
! where in the following array the rows are co,c1,c2,c3,c4
! and c5, and the columns correspond to 
!  1. Thin Ci/Cs (denoted as cloud type 1)
!  2. Thick Ci/Cs (denoted as cloud type 2)
!  3. As/Ac (denoted as cloud type 3)
!  4. Low cloud (denoted as cloud type 4 or 5)

      data wt/0.675d0,1.552d0,1.429d0,1.512d0,                          &
             -3.432d0,-1.957d0,-1.207d0,-1.176d0,                       &
              1.929d0,-1.762d0,-2.008d0,-2.16d0,                        &
              0.842d0,2.067d0,0.853d0,1.42d0,                           &
              2.693d0,0.448d0,0.324d0,-0.032d0,                         &
             -1.354d0,0.932d0,1.582d0,1.422d0/


! initialize variables
      l = 0
      sdown0 = 0d0
      wgt = 0d0
      rcld = 0d0
      tcld = 0d0
      fr = 0d0
      tclr = 0d0
      coszsq = 0d0
      coszcube = 0d0
      d3 = 0d0
      d2 = 0d0
      rg = 0d0
      d1 = 0d0
      rclr = 0d0
      sdowne = 0d0
      sdowns = 0d0
      direct = 0d0
      diffuse = 0d0

      do i=1,3
        icla(i) = 0
        tk(i) = 0d0
        rks(i) = 0d0
        tdk(i) = 0d0
        covera(i) = 0d0
      end do

!  calculate cosine of zenith angle

! comput extra-terrestrial insolation
!  account for variation of sun/earth distance due to
!  eliptical orbit. Redone daily at 0500 hours.
      sdowne = 2d0*pi*float(jday - 2)/365.242d0
      sdowne = (1.0001399d0 + 1.67261d-2*dcos(sdowne))
      sdowne = sdowne*sdowne
      
      coszsq = cosz*cosz
      coszcube = coszsq*cosz
      sdown0 = 1369.2d0*sdowne*cosz

! Equate cloud array elements of Shapiro (high=1, midddle=2, low=3) with
! others (high=3, middle=2, low=1)
      icla(1) = icld(3)                                                  !low
      icla(2) = icld(2)                                                  !middle
      icla(3) = icld(1)                                                  !high
      covera(1) = cover(3)
      covera(2) = cover(2)
      covera(3) = cover(1)

      if(icla(3) /= 0) icla(3) = 4                                       !geoff's value was 4
      if(icla(2) /= 0) icla(2) = 3                                       !geoff's value was 3

      if(icla(1) /= 0) then                                              !geoff's value was 1 (0-0.4) or 2 (0.4-1)
        if(cover(1) > 0d0.or.icla(1) == 5) icla(1) = 1
        if(cover(1) > 0.4d0.or.icla(1) == 7) icla(1) = 2
      end if

!  if precipitation is reported, cover type and fraction
!  are completely determined
      if(prcp1 > eps) then                                               !prcpitation is reported
        covera(1) = 1d0
        covera(2) = 1d0
        covera(3) = 1d0
        icla(1) = 2
        icla(2) = 3
        icla(3) = 4

        if((dabs(cover(3)) <= eps).and.(dabs(cover(2)) <= eps).and.     &
             (dabs(cover(1)) <= eps)) then
          cover(3) = covera(1)
          hgt(3) = 8d0
          icld(3) = 5
          cover(2) = covera(2)
          hgt(2) = 4d0
          icld(2) = 3
          cover(1) = covera(3)
          hgt(1) = 1.5d0
          icld(1) = 6
        end if
!      else
!        cover(3) = covera(1) 
!        icld(3) = icla(1)
!        cover(2) = covera(2)
!        icld(2) = icla(2)
!        cover(1) = covera(3)
!        icld(1) = icla(3)
      end if

      do i=1,3
!  Initialize parameters
        wgt = 0d0
        rcld = 0d0
        tcld = 0d0
        rclr = 0d0
        tclr = 0d0
        j = icla(i)
        if(j > 4) j = 4
        fr = covera(i)
        l = i

        if(j /= 0) then
          wgt = wt(j,1) + wt(j,2)*cosz + wt(j,3)*fr + wt(j,4)*cosz*fr   &
                 + wt(j,5)*coszsq + wt(j,6)*fr*fr
          wgt = wgt*fr       

!          if(fr <  0.05d0) wgt = 0d0
!          if(fr >  0.95d0) wgt = 1d0

          if(wgt > 0d0) then
            rcld = r2(j,1) + r2(j,2)*cosz + r2(j,3)*coszsq              &
                     + r2(j,4)*coszcube
            tcld = t2(j,1) + t2(j,2)*cosz + t2(j,3)*coszsq              &
                     + t2(j,4)*coszcube
          end if
        end if
! Compute reflectivity and transmitivity for each layer
        if(wgt < 1d0) then
          rclr = r1(l,1) + r1(l,2)*cosz + r1(l,3)*coszsq                &
                   + r1(l,4)*coszcube
          tclr = t1(l,1) + t1(l,2)*cosz + t1(l,3)*coszsq                &
                   + t1(l,4)*coszcube
        end if

        rks(i) = wgt*rcld + (1d0 - wgt)*rclr
!        if(rks(i) < 5d-2) rks(i) = 0d0
        tk(i) = wgt*tcld + (1d0 - wgt)*tclr
!        if(tk(i) < 5d-2) tk(i) = 0d0

! Direct componenet of tk
        tdk(i) = dmax1(0d0,tk(i) - rks(i))
      end do

! calculation of insolation at the ground-sdown
      rg = albedo
      d1 = 1d0 - (rks(1)*rks(2))
      d2 = 1d0 - (rks(2)*rks(3))
      d3 = 1d0 - (rks(3)*rg)
      sdowns = d1*d2 - (rks(1)*rks(3)*tk(2)*tk(2))
      sdowns = d3*sdowns - (d1*rks(2)*rg*tk(3)*tk(3))
      sdowns = sdowns - (rks(1)*rg*tk(2)*tk(2)*tk(3)*tk(3))
      sdowns = (tk(1)*tk(2)*tk(3)*sdown0)/sdowns

      if(sdowns <= 0d0) then
        sdowns = 0d0
        direct = 0d0
        diffuse = 0d0
      else
        direct = dmax1(0d0,tdk(1)*tdk(2)*tdk(3)*sdown0)                  !Direct component of insolation
        diffuse = dmax1(0d0,sdowns - direct)                             !Diffuse component of insolation
!        if(diffuse < 1d0) then
!          diffuse = 0d0
!          direct = sdowns
!        end if
      end if

      end subroutine insol

end module module_radiation
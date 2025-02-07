module module_lowveg

      use fasst_global

   contains

      subroutine lowveg_met(ii,iter,pt1,pt2,mixrgr,wetness,taf,dmet, &
                            wetbulba,rhoaf,cf,pheatf,rpp,mixra,      & 
                            mixrf,dqdtf)

! This subroutine calculates the low vegetation root-uptake, precip interception
! and evaporation/transpiration

! calls the following subroutines:
!     sp_humid

! uses the functions: dense,spheats,vap_press   

      implicit none

      integer(kind=4),intent(in):: ii,iter,pt1,pt2
      real(kind=8),intent(in):: mixrgr,wetness,taf
      real(kind=8),intent(inout):: dmet(13),wetbulba
      real(kind=8),intent(out):: rhoaf,cf,pheatf,rpp,mixra,mixrf
      real(kind=8),intent(out):: dqdtf

! saved variables
      real(kind=8):: stlli,stlsi,interc1o,interco
      save:: stlli,stlsi,interc1o,interco

! local variables
      integer(kind=4):: i,f1i,d1i
      real(kind=8):: vpressa,vpressf,pdens,dqdta,wetbulbf,rpf,ra
      real(kind=8):: f1a,f2a,f2,f3,interc,max_wet,rptr,c1,c2,t0(6)
      real(kind=8):: ep,ef,etr,drip,ff,d1,qaf,pdens1,f1,spheat,drip1
      real(kind=8):: min_wat,max_wat,spheat1,interc1,meltlv,kths
      real(kind=8):: max_wetl,max_wets,dense,spheats,t1,vap_press


! initialize variables
      f1i = 0
      d1i = 0
      c1 = 0d0
      c2 = 0d0
      vpressa = 0d0
      vpressf = 0d0
      pdens = 0d0
      pdens1 = 0d0
      dqdta = 0d0
      dqdtf = 0d0
      wetbulbf = 0d0
      mixra = 0d0
      mixrf = 0d0
      qaf = 0d0
      rhoaf = 0d0
      rpf = 0d0
      rptr = 0d0
      ra = 0d0
      ep = 0d0
      ef = 0d0
      etr = 0d0
      ff = 0d0
      d1 = 0d0
      rpp = 0d0
      t1 = 0d0
      cf = 0d0
      spheat = 0d0
      spheat1 = 0d0
      min_wat = 0d0
      max_wat = 0d0
      kths = 0d0
      f1 = 0d0
      f1a = 0d0
      f2a = 0d0
      f2 = 0d0
      f3 = 0d0
      max_wet = 0d0
      max_wetl = 0d0
      max_wets = 0d0
      interc = 0d0
      interc1 = 0d0
      drip = 0d0
      drip1 = 0d0
      meltlv = 0d0
      pheatf = 0d0
      stll = 0d0
      stls = 0d0
      stlli = 0d0
      stlsi = 0d0

      do i=1,6
        t0(i) = 0d0
      end do

      stll = storll                                                      !m
      stls = storls                                                      !m
      if(icaseo /= 4.and.(ii == 0.and.iter == 0)) then
        stlli = stll
        stlsi = stls
      else if(icaseo == 4) then
        stll = 0d0
        stls = 0d0
      end if

      uaf = 8.3d-1*sigfl*dmet(3)*sqrt_chnf + (1d0 - sigfl)*dmet(3)       !foliage wind speed (m/s)

      if(uaf > 2d-1) then  !eps) then
        cf = (1d0 + 3d-1/uaf)*1d-2                                       !bulk transfer coefficient (unitless)
        ra = 1d0/(cf*uaf)                                                !atmospheric resistance to water vapor diffusion (s/m)
      else
        cf = 1d-2  !0d0
        ra = 0d0
      end if
      cf = anint(sigfl*cf*1d20)*1d-20

!      rhoaf = (3.48d-3*dmet(11)*1d2*5d-1)*(1d0/dmet(4) + 1d0/ftemp)      !foliage air density (kg/m^3)
      c1 = vap_press(nnodes+1,dmet(5)*1d-2,dmet(11))
!      rhoaf = (c1/Rv + (dmet(11)*1d2 - c1)/Rd)*5d-1*                    &
!                                              (1d0/dmet(4) + 1d0/ftemp)  !foliage air density (kg/m^3)
      rhoaf = (c1/Rv + (dmet(11)*1d2 - c1)/Rd)/taf                       !foliage air density (kg/m^3)
      rhoaf = anint(rhoaf*1d20)*1d-20

! solve for the water vapor densities
      if(hsaccum > eps.or.hi > eps.or.stls > eps) f1i = 1

      c2 = dmet(5)*1d-2
      call sp_humid(f1i,dmet(11),dmet(4),c2,c1,mixra,t1,vpressa,        &
                    wetbulba,t0(1),t0(2),t0(3),t0(4),t0(5),t0(6))

      c2 = 1d0
!      call sp_humid(f1i,dmet(11),ftemp,c2,c1,mixrf,t1,vpressf,wetbulbf, &
!                    t0(1),t0(2),t0(3),t0(4),t0(5),t0(6))
      call sp_humid(f1i,dmet(11),ftemp,c2,c1,mixrf,dqdtf,vpressf,       &
                    wetbulbf,t0(1),t0(2),t0(3),t0(4),t0(5),t0(6))

! calculate the stomatal resistance, rs, and max transpiration
! REF: http://www.ecmwf.int/research/ifsdocs/CY25r1/PHYSICS/
      f1a = (4d-3*dmet(1) + 5d-3)/(8.1d-1*(4d-3*dmet(1) + 1d0))
      f1 = 1d0/dmin1(1d0,dmax1(eps,f1a))                                 !unitless

! roots and transpiration
      trmlm = 0d0                                                        !maximum transpiration rate (m/s)
      do i=1,nnodes
        if(stt(i) > Tref) then
          frl(i) = rk(vegl_type,i)*(1d0 - (soil_moist(i) - nsoilp(i,16))&
                                        /(nsoilp(i,17) - nsoilp(i,16)))  !root water fraction
          frl(i) = dmax1(0d0,dmin1(1d0,anint(frl(i)*1d20)*1d-20))
         else
           frl(i) = 0d0
        end if

        trmlm = trmlm + frl(i)
        f2a = f2a + rk(vegl_type,i)*soil_moist(i)
        min_wat = min_wat + nsoilp(i,16)
        max_wat = max_wat + nsoilp(i,17)
      end do

      f2 = dmax1(0d0,dmin1(1d0,(f2a - min_wat)/(max_wat - min_wat)))
      if(f2 > eps) f2 = 1d0/f2   !1d-3

      ff = 1d0
!      if(stt(refn) <= Tref.or.hm > 5d-2) ff = 0d0
!      if(stt(nnodes) <= Tref.or.hm > 5d-2) ff = 0d0
!      if(ice(nnodes) > eps.or.hm > 5d-2) ff = 0d0
      trmlm = 1.5d-7*ff*trmlm  !*sigfl                                    !m/s

      t1 = 3d-4*(vpressf - vpressa)
      if(dabs(vpressf-vpressa) > eps.and.t1 < 5d1) then
        f3 = dexp(t1)                                                    !unitless
      else
        f3 = 1d0
      end if

      if(lail > eps.and.ftemp > Tref) then
        rsl = dmax1(0d0,(veg_prp(vegl_type,1)/lail)*f1*f2*f3)            !s/m, stomatal resistance
      else
        rsl = 0d0
      end if

!      if(rsl > eps.and.ra+rsl > eps) rpp = ra/(ra + rsl)                 !vegetation wetness
      if(ra+rsl > eps) rpp = ra/(ra + rsl)                               !vegetation wetness
      rpp = dmin1(dmax1(0d0,anint(rpp*1d20)*1d-20),1d0)

      d1 = 1d0 - sigfl*(6d-1*(1d0 - rpp) + 1d-1*(1d0 - wetness))
      qaf = ((1d0 - 7d-1*sigfl)*mixra + 6d-1*sigfl*rpp*mixrf            &
                                     + 1d-1*sigfl*wetness*mixrgr)/d1     !unitless
      d1i = 1
      ep = lail*cf*(rhoaf/dense(wetbulba,0d0,d1i))*uaf*                 &
                                                      (qaf - rpp*mixrf)  !potential evaporation (m/s)

! rain and snow density and specific heat
      d1i = 1
      pdens = dense(wetbulba,0d0,d1i)                                    !kg/m^3
      spheat = spheats(wetbulba,d1i)                                     !J/kg*K

      d1i = 4
      pdens1 = dense(wetbulba,dmet(3),d1i)                               !kg/m^3
      d1i = 2
      spheat1 = spheats(wetbulba,d1i)                                    !J/kg*K

! precip interception, maximum storage
!     max_wet = 0.935 + 0.498*lai - 0.00575*lai*lai   !H-H
      max_wet = veg_prp(vegl_type,7)*(lail + veg_prp(vegl_type,8))       !R & S (mm)
      max_wetl = anint((max_wet*1d-3)*1d20)*1d-20                        !m
      d1i = 1
      if(pdens1 > eps) then
        max_wets = (max_wet*1d-3)*dense(ftemp,0d0,d1i)/pdens1            !m
        max_wets = anint(max_wets*1d20)*1d-20
      end if

if(ii == 0.and.iter == 0) then
      t1 = 5d-1*(lail + veg_prp(vegl_type,8))
      if(t1 > 5d1) t1 = 5d1
      if(pt1 == 2.or.pt1 == 4)then           !rain
        interc = dmet(6)*(1d0 - dexp(-t1))                               !m
        interc = dmax1(0d0,dmin1(interc,max_wetl))
      else if(pt1 == 3) then                 !snow
        interc1 = dmet(6)*(1d0 - dexp(-t1))                              !m
        interc1 = dmax1(0d0,dmin1(interc1,max_wets))
      else if(pt2 == 3) then                 !snow
        interc1 = dmet(7)*(1d0 - dexp(-t1))                              !m
        interc1 = dmax1(0d0,dmin1(interc1,max_wets))
      end if
      interc = anint(interc*1d20)*1d-20
      interc1 = anint(interc1*1d20)*1d-20
      interco = interc
      interc1o = interc1

      stls = stls + interc1
      if(stls > max_wets) then
        drip1 = stls - max_wets
        stls = max_wets
      end if
      drip1 = anint(drip1*1d20)*1d-20
      stls = anint(stls*1d20)*1d-20
      if(stls < 1d-4) stls = 0d0

      meltlv = 0d0
      if(stls > 0d0) then
!        kths = 0.023d0 + (7.75d-5*pdens1                                &
!                       + 1.105d-6*pdens1*pdens1)*(2.29d0 - 0.023d0)      !W/m*K
        if(pdens1 < 1.56d-1) then                                        !W/m*K, th. cond. (Sturm et al)
          kths = 2.3d-2 + 2.34d-1*pdens1
          kths = dmin1(dmax1(2.3d-2,kths),1d0)
        else
          kths = 1.38d-1 - 1.01d0*pdens1 + 3.233d0*pdens1*pdens1
          kths = dmin1(dmax1(1.38d-1,kths),1d0)
        end if
        meltlv = dmax1(0d0,(kths/stls)*(dmet(4) - Tref)*                &
                                                (deltat/(lhfus*pdens1)))
        if(meltlv > stls*pdens/pdens1) meltlv = stls*pdens/pdens1
        meltlv = anint(meltlv*1d20)*1d-20
        stls = stls - meltlv

        stls = anint(stls*1d20)*1d-20
        if(stls < 1d-4) stls = 0d0
      end if

      stll = stll + interc + meltlv
      drip = 0d0
      if(stll > max_wetl) then
        drip = stll - max_wetl
        stll = max_wetl
      end if
      drip = anint(drip*1d20)*1d-20
      stll = anint(stll*1d20)*1d-20
      if(stll < 1d-6) stll = 0d0
else
  interc = interco
  interc1 = interc1o
end if

      if(rpp*mixrf >= qaf.or.trmlm <= eps) then
        if(f2 > min_wat.and.ep <= eps) ep = 0d0
      end if

      if(rsl+ra > eps) then
        if(stls > eps) then
          ff = (stls/max_wets)**(2d0/3d0)
        else if(stll > eps) then
          ff = (stll/max_wetl)**(2d0/3d0)
        end if
        rpf = 1d0 - (rsl/(rsl + ra))*(1d0 - ff)
        rptr = (ra/(rsl + ra))*(1d0 - ff)
      end if

      rptr = dmin1(dmax1(0d0,rptr),1d0)
      rptr = anint(rptr*1d20)*1d-20
      rpf = dmin1(rpf,1d0)
      rpf = anint(rpf*1d20)*1d-20

      ef = rpf*ep                                                        !leaf potential evaporation (m/s)
      etr = dmin1(trmlm,rptr*ep)                                         !leaf potential transpiration (m/s)

if(ii == 0.and.iter == 0) then
      if(stls > eps) then
        stls = stls - (ef - etr)*deltat                                  !m
        if(stls < 1d-4) stls = 0d0
        if(stls > max_wets) then
          drip1 = drip1 + (stls - max_wets)
          stls = max_wets
        end if
      else if(stll > eps) then
        stll = stll  - (ef - etr)*deltat                                 !m
        if(stll < 1d-6) stll = 0d0
        if(stll > max_wetl) then
          drip = drip + (stll - max_wetl)
          stll = max_wetl
        end if
      end if

      drip = dmax1(0d0,anint(drip*1d20)*1d-20)
      stll = dmax1(0d0,anint(stll*1d20)*1d-20)
      drip1 = dmax1(0d0,anint(drip1*1d20)*1d-20)
      stls = dmax1(0d0,anint(stls*1d20)*1d-20)

      dmet(6) = dmet(6) + sigfl*drip/(timstep*3.6d3)                           !m  !sigfl*
      dmet(7) = dmet(7) + sigfl*drip1/(timstep*3.6d3)                          !m  !sigfl*

      dmet(6) = anint(dmet(6)*1d15)*1d-15
      dmet(7) = anint(dmet(7)*1d15)*1d-15
end if

! HEAT LOSS DUE TO PRECIPITATION
! based on Jordan, R. CRREL Rep. 91-16, pp.29,15-17
! Note:  All freezing and thawing taken care of in snow_model.f and icethick.f

      pheatf = sigfl*(spheat*(pdens*interc) + spheat1*(pdens1*interc1))  !W/m^2*K
      pheatf = anint(pheatf*1d20)*1d-20

      end subroutine lowveg_met

! ******************************************************************************
      subroutine low_veg_prop(fatemp,c3_hgt)

! no subroutines called

      implicit none

      real(kind=8),intent(in):: fatemp,c3_hgt

! local variables
      integer(kind=4):: season
!      real(kind=8):: lowveg_min,lowveg_max,medveg_min,medveg_max,hiveg
!      real(kind=8):: epf_max,epf_min,fols_max,fols_min
      real(kind=8):: hgt1,ftg,hfol_old,sigfl_old,albf_old,epf_old
      real(kind=8):: hfoltot_old,fola_old,lail_old

!      data hiveg/185d0/                      !high veg height (cm)
!      data lowveg_min/5d0/                   !low veg height (cm)
!      data lowveg_max/50d0/
!      data medveg_min/5d0/                   !med veg height (cm)  !35d0
!      data medveg_max/90d0/                                        !60d0

!      data fols_min/0.70d0/                  !minimum fola(absorptivity = 1 - albedo)
!      data fols_max/0.85d0/                  !maximum fola(absorptivity = 1 - albedo)
!      data epf_max/0.96d0/                   !maximum foliage emissivity (unitless)
!      data epf_min/0.90d0/                   !minimum foliage emissivity (unitless)

! zero variables
      ftg = 0d0
      hgt1 = 0d0
      hfol_old = 0d0
      sigfl_old = 0d0
      albf_old = 0d0
      epf_old = 0d0
      lail_old = 0d0
      hfoltot_old = 0d0
      fola_old = 0d0

! determine season
      if(lat >= 0d0) then                                                !northern hemisphere
        if(aint(met(iw,ip_doy)) >= 335.or.                              &
           aint(met(iw,ip_doy)) <=  80) then
          season = 1                                                     !winter
        else if(aint(met(iw,ip_doy)) >=  81.and.                        &
                aint(met(iw,ip_doy)) <= 151) then
          season = 2                                                     !spring 
        else if(aint(met(iw,ip_doy)) >= 152.and.                        &
           aint(met(iw,ip_doy)) <= 243) then
          season = 3                                                     !summer
        else if(aint(met(iw,ip_doy)) >= 244.and.                        &
                aint(met(iw,ip_doy)) <= 334) then
          season = 4                                                     !fall
        end if
      else                                                               !southern hemisphere
        if(aint(met(iw,ip_doy)) >= 335.or.                              & 
           aint(met(iw,ip_doy)) <=  80) then
          season = 3                                                     !summer
        else if(aint(met(iw,ip_doy)) >=  81.and.                        &
                aint(met(iw,ip_doy)) <= 151) then
          season = 4                                                     !fall 
        else if(aint(met(iw,ip_doy)) >= 152.and.                        &
           aint(met(iw,ip_doy)) <= 243) then
          season = 1                                                     !winter
        else if(aint(met(iw,ip_doy)) >= 244.and.                        &
                aint(met(iw,ip_doy)) <= 334) then
          season = 2                                                     !spring
        end if
      end if

      if(iw == istart.and.infer_test == 0) iseason = season

!      ftg = 1d0 - 1.6d-4*(298d0 - toptemp)*(298d0 - toptemp)  !1.6d-3
      ftg = 1d0 - 1.6d-3*(298d0 - stt(refn))*(298d0 - stt(refn))
      ftg = anint(ftg*1d20)*1d-20

      sigfl_old = sigfl
      lail_old = lail
      hfol_old = hfol
      hfoltot_old = hfol_tot
      epf_old = epf
      fola_old = fola

      if(vegl_type == 2) then
        lail_old = lail
        lail = veg_prp(vegl_type,6) + ftg*(veg_prp(vegl_type,5)         &
                    - veg_prp(vegl_type,6))                              !leaf area index (% vol/vol)
        if(lail > veg_prp(vegl_type,5)) lail = veg_prp(vegl_type,5)
        if(lail < veg_prp(vegl_type,6).and.lail < ilail)                & 
          lail = veg_prp(vegl_type,6)
        if(iw == istart.and.infer_test == 0) ilail = lail
        if(season < iseason.and.lail > ilail) lail = ilail
        if(dabs(-0.75d0*lail) < 5d1) sigfl = 1d0 - dexp(-0.75d0*lail)
      else
        sigfl = veg_prp(vegl_type,2) - (1d0 - ftg)                      &
                   *(veg_prp(vegl_type,2) - veg_prp(vegl_type,3))
        if(sigfl < veg_prp(vegl_type,2)) sigfl = veg_prp(vegl_type,2)
        if(sigfl > veg_prp(vegl_type,3)) sigfl = veg_prp(vegl_type,3)
        sigfl = sigfl*1d-2
      end if

      if(season == iseason.and.dabs(isigfl-spflag) > eps) sigfl = isigfl
      if((season < iseason.and.sigfl > isigfl).and.                     &
         dabs(isigfl-spflag) > eps) sigfl = isigfl

      if(sigfl > eps) then
! low vegetation (short grass, desert)
        if(vegl_type == 2.or.vegl_type == 8) then
!          hfol = veg_prp(vegl_type,17) - (1d0 - 1d-1*ftg)*              &
          hfol = veg_prp(vegl_type,17) - (1d0 - ftg)*              &
                        (veg_prp(vegl_type,17) - veg_prp(vegl_type,16))
          if(hfol < veg_prp(vegl_type,16)) hfol = veg_prp(vegl_type,16)
          if(hfol > veg_prp(vegl_type,17)) hfol = veg_prp(vegl_type,17)
! high vegetation (evergreen & decid. shrubs)
        else if(vegl_type == 16.or.vegl_type == 17) then
          hfol = veg_prp(vegl_type,16)
! medium vegetation (crops, tall grass, tundra, irr. crops, semidesert)
        else
!          hfol = veg_prp(vegl_type,17) - (1d0 - 1d-1*ftg)*              &
          hfol = veg_prp(vegl_type,17) - (1d0 - ftg)*              &
                        (veg_prp(vegl_type,17) - veg_prp(vegl_type,16))
          if(hfol < veg_prp(vegl_type,16)) hfol = veg_prp(vegl_type,16)
          if(hfol > veg_prp(vegl_type,17)) hfol = veg_prp(vegl_type,17)
        end if

! adjust the foliage height if necessary
        if(season == iseason.and.dabs(ihfol-spflag) > eps) hfol = ihfol
        if(season < iseason.and.hfol > ihfol) hfol = ihfol

!        if(vegh_type /= 0.and.hfol >= c3_hgt*1d2)    &
!          hfol = (c3_hgt*1d2)*5d-1                                       !c3_hgt = canopy bottom layer height
        hfol_tot = hfol*1d-2
        if(iw /= 1) hfol_tot = 5d-1*(hfol_tot + hfoltot_old)
        hfol_tot = anint(hfol_tot*1d20)*1d-20

        hfol = dmax1(0d0,hfol*1d-2 - hm)

        if(iheightn <= hfol) hfol = iheightn*5d-1

        if(iw /= 1) then
          if(dabs(met(iw-1,ip_sd)+met(iw-1,ip_hi)) < eps)               &
            hfol = 5d-1*(hfol + hfol_old)
        end if

        hfol = anint(hfol*1d20)*1d-20

        hgt1 = iheightn - hm
        if(hgt1 <= hfol) hgt1 = hfol

! compute parameters used by latent and sensible heat terms
        if(hfol <= 0d0) then                                             !snow on ground
          z0l = 7.775d-3                                                 !snow surface roughness length (m)
          Zd = z0l                                                       !zero displacement height (m)
          sqrt_chnf = 0d0
        else
          z0l = 0.131d0*(hfol**0.997d0)                                  !low vegetation roughness length (m)
          z0l = dmin1(z0l,veg_prp(vegl_type,4))
          Zd = 0.701d0*(hfol**0.979d0)
          sqrt_chnf = vK/(dlog((hgt1 - Zd)/z0l))
        end if
        z0l = anint(z0l*1d20)*1d-20
        Zd = anint(Zd*1d20)*1d-20
        sqrt_chnf = anint(sqrt_chnf*1d20)*1d-20

! calculate leaf area index, foliage density
! from Ramirez and Senarath, J. Climate, 13(22), p.4050-
        if(vegl_type /= 2) then
          lail = veg_prp(vegl_type,6) + ftg*(veg_prp(vegl_type,5)       &
                    - veg_prp(vegl_type,6))                              !leaf area index (% vol/vol)
          if(lail > veg_prp(vegl_type,5)) lail = veg_prp(vegl_type,5)
          if(lail < veg_prp(vegl_type,6).and.lail < ilail)              & 
                                            lail = veg_prp(vegl_type,6)
          if(iw == istart.and.infer_test == 0) ilail = lail
          if(season < iseason.and.lail > ilail) lail = ilail
        end if

        epf = veg_prp(vegl_type,12) + ftg*(veg_prp(vegl_type,13) -      &
                                                veg_prp(vegl_type,12))   !emissivity
        if(epf < veg_prp(vegl_type,12)) epf = veg_prp(vegl_type,12)
        if(epf > veg_prp(vegl_type,13)) epf = veg_prp(vegl_type,13)
        if(season == iseason.and.dabs(iepf-spflag) > eps) epf = iepf
        if(season < iseason.and.epf > iepf) epf = iepf

        fola = veg_prp(vegl_type,15) - (1d0 - ftg)*                     &
                        (veg_prp(vegl_type,15) - veg_prp(vegl_type,14))  !absorption
        if(fola < veg_prp(vegl_type,14)) fola = veg_prp(vegl_type,14)
        if(fola > veg_prp(vegl_type,15)) fola = veg_prp(vegl_type,15)

        if(iw == istart.and.infer_test == 0) then
          if(dabs(ifola-spflag) > eps) fola = ifola
          ftemp = fatemp                                                !K
          uaf = met(iw,ip_ws)
        else
          lail = 5d-1*(lail_old + lail)
          sigfl = 5d-1*(sigfl_old + sigfl)
          fola = 5d-1*(fola_old + fola)
          epf = 5d-1*(epf_old + epf)
        end if

        if(hfol <= 0d0.and.hm <= eps) then
          sigfl = 0d0
          lail = 0d0
          fola = 0d0
          epf = 0d0
        end if

        albf = 1d0 - fola
      else
        epf = 0d0
        albf = 0d0
        lail = 0d0
        hfol = 0d0
        hfol_tot = 0d0
        hgt1 = iheightn - hm
        z0l = 0d0
        Zd = 0d0
        sqrt_chnf = 0d0
      end if

      sigfl = anint(sigfl*1d20)*1d-20
      lail = anint(lail*1d20)*1d-20
      fola = anint(fola*1d20)*1d-20
      albf = anint(albf*1d20)*1d-20
      epf = anint(epf*1d20)*1d-20
      hgt1 = anint(hgt1*1d20)*1d-20

      end subroutine low_veg_prop

! ******************************************************************************
      subroutine veg_propl(biome_source,new_vt,veg_type)

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: biome_source,new_vt,veg_type 

! local variables
      integer(kind=4):: i,j,jj,io,vid,hlines,fid,ntypes,err
      real(kind=8):: zup,zdw,t1,t2,t3,t4,ar,br,emissmin,emissmax
      real(kind=8):: srmax,srmin,cmax,cmin,rl,laimax,laimin,ddmax,sai
      real(kind=8):: folamin,folamax,heigmin,heigmax,defveg_prp(18,17)
      real(kind=8),allocatable:: newveg_prp(:,:)
      character(len=200):: header


! vegetation types (veg_type) based on BATS parameterizations
!     and numbering system
! LOW vegetation types: 1, 2, 7-11, 13, 16, 17
! HIGH vegetation types: 3-6, 18
! NON vegetation types: 12, 14, 15
!     1 = crop/mixed farming (l)          10 = irrigated crop (l)
!     2 = short grass (l)                 11 = semidesert (l)
!     3 = evergreen needle-leaf(h)        12 = ice cap/glacier (n)
!     4 = deciduous needle-leaf(h)        13 = bog/marsh (l)
!     5 = deciduous broadleaf(h)          14 = inland water (n)
!     6 = evergreen broadleaf(h)          15 = ocean (n)
!     7 = tall grass (l)                  16 = evergreen shrub (l)
!     8 = desert (l)                      17 = deciduous shrub (l)
!     9 = tundra (l)                      18 = mixed woodland(h)

! veg_prp(veg_type,i) where i = 
!    1 = minimum stomatal resistance (s/m)
!    2 = maximum coverage (%)
!    3 = minimum coverage (%)
!    4 = roughness length (m)
!    5 = maximum LAI
!    6 = minimum LAI
!    7 = maximum dew depth (mm)
!    8 = stem area index
!    9 = ar, root distribution calculation coeff.
!   10 = br, root distribution calculation coeff.
!   11 = maximum stomatal resistance (s/m)
!   12 = emissivity max (%)
!   13 = emissivity min (%)
!   14 = albedo min (%)
!   15 = albedo max (%)
!   16 = height min (cm)
!   17 = height max (cm)


! zero-out parameters
      zup = 0d0
      zdw = 0d0

      vid = 0
      hlines = 0
      fid = 0
      ntypes = 0
      srmax = 0d0
      srmin = 0d0
      cmax = 0d0
      cmin = 0d0
      rl = 0d0
      laimax = 0d0
      laimin = 0d0
      ddmax = 0d0
      sai = 0d0
      ar = 0d0
      br = 0d0
      emissmax = 0d0
      emissmin = 0d0
      folamin = 0d0
      folamax = 0d0
      heigmin = 0d0
      heigmax = 0d0

      do i=1,36
        read(31,'(a)') header
      end do

      jj = 0
      io = 0
      do while(io /= -1.and.jj == 0)
        read(31,*,iostat=io) vid,srmax,srmin,cmax,cmin,rl,laimax,laimin,&
                             ddmax,sai,ar,br,emissmin,emissmax,folamin, &
                             folamax,heigmin,heigmax

        if(veg_type == vid) then
          jj = 1
          defveg_prp(veg_type,1) = srmin
          defveg_prp(veg_type,2) = cmax
          defveg_prp(veg_type,3) = cmin
          defveg_prp(veg_type,4) = rl
          defveg_prp(veg_type,5) = laimax
          defveg_prp(veg_type,6) = laimin
          defveg_prp(veg_type,7) = ddmax
          defveg_prp(veg_type,8) = sai
          defveg_prp(veg_type,9) = ar
          defveg_prp(veg_type,10) = br
          defveg_prp(veg_type,11) = srmax
          defveg_prp(veg_type,12) = emissmin*1d-2
          defveg_prp(veg_type,13) = emissmax*1d-2
          defveg_prp(veg_type,14) = folamin*1d-2
          defveg_prp(veg_type,15) = folamax*1d-2
          defveg_prp(veg_type,16) = heigmin
          defveg_prp(veg_type,17) = heigmax
        end if
      end do

      rewind(31)
 
      if (biome_source > 0) then
        if(biome_source == 1000) then                                      !Modis_NOAH
          fid = 32
          hlines = 45
          ntypes = 20
        else if(biome_source == 2000) then                                 !UMD
          fid = 33
          hlines = 31
          ntypes = 14
        end if

        do i=1,hlines
          read(fid,'(a)') header
        end do

        allocate(newveg_prp(ntypes,17),stat=err)
        do i=1,ntypes
          do j=1,17
            newveg_prp(i,j) = 0d0
          end do
        end do
 
        jj = 0
        io = 0
        do while(io /= -1.and.jj == 0)
          read(fid,*,iostat=io) vid,srmax,srmin,cmax,cmin,rl,laimax,    &
                                laimin,ddmax,sai,ar,br,emissmin,        &
                                emissmax,folamin,folamax,heigmin,       &
                                heigmax

          if(new_vt == vid) then
            jj = 1
            newveg_prp(new_vt,1) = srmin
            newveg_prp(new_vt,2) = cmax
            newveg_prp(new_vt,3) = cmin
            newveg_prp(new_vt,4) = rl
            newveg_prp(new_vt,5) = laimax
            newveg_prp(new_vt,6) = laimin
            newveg_prp(new_vt,7) = ddmax
            newveg_prp(new_vt,8) = sai
            newveg_prp(new_vt,9) = ar
            newveg_prp(new_vt,10) = br
            newveg_prp(new_vt,11) = srmax
            newveg_prp(new_vt,12) = emissmin*1d-2
            newveg_prp(new_vt,13) = emissmax*1d-2
            newveg_prp(new_vt,14) = folamin*1d-2
            newveg_prp(new_vt,15) = folamax*1d-2
            newveg_prp(new_vt,16) = heigmin
            newveg_prp(new_vt,17) = heigmax
          end if
        end do
    
        rewind(fid)
      end if

      if(biome_source == 0) then
        do i=1,17
          veg_prp(veg_type,i) = defveg_prp(veg_type,i)
          veg_prp(veg_type,i) = anint(veg_prp(veg_type,i)*1d10)*1d-10
        end do
      else
        do i=1,17
          if(dabs(newveg_prp(new_vt,i)-spflag) <= eps) then 
            veg_prp(veg_type,i) = defveg_prp(veg_type,i)
          else
            veg_prp(veg_type,i) = newveg_prp(new_vt,i)
          end if
          veg_prp(veg_type,i) = anint(veg_prp(veg_type,i)*1d10)*1d-10
        end do
      end if

! root fraction in layer i for each veg type
      do i=1,nnodes
        if(i == 1) then
          zdw = elev - nzi(i)
          zup =  elev - (nzi(i) + nzi(i+1))*5d-1

        else if(i == nnodes) then
          zdw = elev - (nzi(i) + nzi(i-1))*5d-1
          zup = elev - nzi(i)
        else
          zdw = elev - (nzi(i) + nzi(i-1))*5d-1
          zup = elev - (nzi(i) + nzi(i+1))*5d-1
        end if

        t1 = 0d0
        t1 = veg_prp(veg_type,9)*zdw
        if(t1 > 5d1) t1 = 5d1

        t2 = 0d0
        t2 = veg_prp(veg_type,10)*zdw
        if(t2 > 5d1) t2 = 5d1

        t3 = 0d0
        t3 = veg_prp(veg_type,9)*zup
        if(t3 > 5d1) t3 = 5d1

        t4 = 0d0
        t4 = veg_prp(veg_type,10)*zup
        if(t4 > 5d1) t4 = 5d1

        rk(veg_type,i) = -5d-1*(dexp(-t1) + dexp(-t2)                   &
                                               - dexp(-t3) - dexp(-t4))
        rk(veg_type,i) = anint(rk(veg_type,i)*1d10)*1d-10
        if(dabs(rk(veg_type,i)) < eps) rk(veg_type,i) = 0d0
      end do

      deallocate(newveg_prp,stat=err)

      end subroutine veg_propl

end module module_lowveg

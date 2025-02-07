      subroutine fasst_main(nprint,sprint,frozen,nstr_flag,sdensi,      &
                            surfdepth,code2,code3,code4,code5,          &
                            first_loop,phie,pdens,oldsd,oldhi)


      use fasst_global
      use module_canopy
      use module_lowveg
      use module_snow

! The following program is the shell for integrating the solar and IR flux models,
! snow melt model, the soil temperature model, soil moisture and strength model, 
! and the freeze/thaw model.

! MAIN CALCULATION LOOP
!     STEP 1: Initialize met parameters; low vegetation properties

!     STEP 2: Calculate the canopy temperature and moisture if present
 
!     STEP 3: Calculate the soil temperature and soil moisture.

!     STEP 3a: Calculate the upward solar and IR fluxes.

!     STEP 4: Calculate soil strength.

!     STEP 5: If there is snow or ice on the ground, check to see if it is melting. If yes,
!             calculate by how much and also calculate the amount of run-off
!             generated.

!     STEP 6: Calculate ouput variables; print nodal information

! END MAIN CALCULATION LOOP

! ******************************************************************************
! UNIT CONVERSION NOTES: 
!     W = J/s; J/m = N; mbar = 100*Pa

!     % H2O by vol = vol H2O/tot vol; max vol H2O = porosity
!     % H2O by mass = mass water/mass soil [gravimetric soil water content]
!     % H2O by mass = (% H2O by vol)*(density water/density soil) assuming mass air neglig.
!     NOTE: maximum % H2O by vol = porosity 
! ******************************************************************************

! VARIABLES:

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

! SOIL PROPERTIES (isoilp(i,j):initial,nsoilp(k,j):nodal and soilp(i,j):layer
!                    where i=layers, k=nodes, j=soil property)
!  isoilp(i,1) = bulk dry density (g/cm^3)
!  isoilp(i,2) = porosity
!  isoilp(i,3) = albedo
!  isoilp(i,4) = emissivity
!  isoilp(i,5) = quartz content
!  isoilp(i,6) = dry thermal conductivity (W/m*K)
!  isoilp(i,7) = saturated hydraulic conductivity (cm/s)
!  isoilp(i,8) = minumum volumetric water content
!  isoilp(i,9) = maximum volumetric water content
!  isoilp(i,10) = van Genuchten's alpha (1/cm)
!  isoilp(i,11) = van Genuchten's n
!  isoilp(i,12) = van Genuchten's m (m = 1 - 1/n)
!  isoilp(i,13) = specific heat of dry soil (J/kg*K)
!  isoilp(i,14) = organic fraction (vol/vol)
!  isoilp(i,15) = minimum allowed volumetric water content
!  isoilp(i,16) = wilting point volumetric water content
!  isoilp(i,17) = field capacity volumetric water content
!  isoilp(i,18) = % sand
!  isoilp(i,19) = % silt
!  isoilp(i,20) = % clay
!  isoilp(i,21) = % carbon
!  isoilp(i,22) = plactic limit
!  isoilp(i,23) = percent fines passing #200 sieve
!  isoilp(i,24) = 0.999d0*soilp(i,9) => to prevent numerical instability, max water content
!  isoilp(i,25) = % gravel
! ******************************************************************************

! calls the following subroutines:
!     low_veg_prop
!     canopy_met
!     new_profile
!     soil_strength
!     snow_model
!     icethick

! uses the functions: dense,met_date

      implicit none

      integer(kind=4),intent(in):: nprint,sprint,frozen
      integer(kind=4),intent(in):: nstr_flag(maxn)
      real(kind=8),intent(in):: sdensi,surfdepth
      integer(kind=4),intent(inout):: code2,code3,code4,code5,first_loop
      real(kind=8),intent(inout):: phie,pdens,oldsd,oldhi

! saved variables
      integer(kind=4):: first_snow,first_time

      save:: first_snow,first_time

! local variables
      integer(kind=4):: i,j,k,fts,d1i,d2i,wstart
      real(kind=8):: surfm,surfi,surfmp,surfip,denom,met_date
      real(kind=8):: simelt,surfdd,smm,fsup,firup,num2,dense,dates
      real(kind=8):: xtemp,pfrozen,timeo(2,4),t1,sfdp,dwi,dwi1

      character(len=1):: nsstate


! initialize certain check/test variables
      d1i = 0
      d2i = 0
      wstart = 0
      simelt = 0d0
      num2 = 0d0
      sfdp = 0d0
      dwi = 0d0
      dwi1 = 0d0
      nsstate = ' '

      if(first_loop == 0) then
        first_time = 0                                                   !first time through new profile routines
        first_snow = 0                                                   !first time through snow routines
        vsmelt = 0d0                                                     !snow melt (m)
        vimelt = 0d0                                                     !ice melt (m)
        first_loop = 1
      end if

!  !!!!!!!!!!!!!!!! Main Calculation Loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(iw-anint(iw*1d-2)*1d2 == 0.and.single_multi_flag == 0) then
        write(*,*)' Start of run',iw,'of',iend,'runs'

      else if(iend > 50.and.(iw-anint(iw*2d-2)*5d1 == 0.and.            &
                                      single_multi_flag == 1)) then
        write(*,*)' Start of run',iw,'of',iend,'runs'
      end if

! initialize the met calc variable dmet1(iw,)
      dmet1(iw,1) = met(iw,ip_tsol)
      if(dmet1(iw,1) < 1d-2) dmet1(iw,1) = 0d0
      dmet1(iw,2) = met(iw,ip_ir)
      dmet1(iw,3) = met(iw,ip_ws)
      dmet1(iw,4) = met(iw,ip_tmp) + Tref                                !K
      dmet1(iw,5) = dmax1(1d-1,met(iw,ip_rh))                            !%
      if(met(iw,ip_prec) <= eps) then
        met(iw,ip_prec) = 0d0
        met(iw,ip_pt) = 1
      end if

      dmet1(iw,6) = met(iw,ip_prec)*1d-3                                 !m/timstep (rain and/or snow)
      if(dmet1(iw,6) <= eps) met(iw,ip_pt) = float(1)

      if(met(iw,ip_prec2) <= eps) then
        met(iw,ip_prec2) = 0d0
        met(iw,ip_pt2) = 1
      end if
      dmet1(iw,7) = met(iw,ip_prec2)*1d-3                                !m/timstep (snow only)
      if(dmet1(iw,7) <= eps) met(iw,ip_pt2) = float(1)

      dmet1(iw,8) = met(iw,ip_upsol)
      if(dabs(dmet1(iw,1)) <= eps) dmet1(iw,8) = 0d0
      dmet1(iw,9) = met(iw,ip_dir)
      dmet1(iw,10) = met(iw,ip_dif)
      dmet1(iw,11) = met(iw,ip_ap)
      dmet1(iw,12) = met(iw,ip_irup)
      dmet1(iw,13) = met(iw,ip_tsoil)

      do i=1,13
        dmet1(iw,i) = anint(dmet1(iw,i)*1d15)*1d-15
      end do

! get the low vegetation properties
      if(veg_flagl == 1) then
        call low_veg_prop(dmet1(iw,4),dzveg(nclayers))
      else
        sigfl = 0d0
      end if

! STEP 2
! calculate the canopy temperature and moisture
      if(veg_flagh == 1) then
        call canopy_met(oldsd)                                                  !correct met for canopy
        if(dmet1(iw,7) > 0d0) met(iw,ip_pt2) = float(3)
        if(dmet1(iw,6) <= eps) met(iw,ip_pt) = float(1)
        if(dmet1(iw,7) <= eps) met(iw,ip_pt2) = float(1)
      else
        surfemisc(iw) = mflag
        do i=1,nclayers
          canopy_temp(i,iw) = mflag
         end do
      end if
      if(error_type == 2) code2 = 1

      airt(1,iw) = met(iw,ip_tmp) + Tref
      airt(2,iw) = dmet1(iw,4)                          

! STEP 3
! determine density and amount of new snow if no existing snow on ground
      newsd = 0d0
      pdens = 0d0
      smm = dmet1(iw,6) + dmet1(iw,7)
      if((aint(met(iw,ip_pt)) == 3.or.aint(met(iw,ip_pt2)) == 3)        &
                                                   .and.smm > 0d0) then
        d1i = 4
        pdens = dense(dmet1(iw,4),dmet1(iw,3),d1i)                       !kg/m^3

        if(aint(met(iw,ip_pt)) == 3) then
          newsd = dmet1(iw,6)                                            !m snow
        else
          newsd = dmet1(iw,7)
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

! calculate the new moisture and temperature profiles
      simelt = vsmelt + vimelt
      hm = hsaccum + newsd + hi                                          !m

      call new_profile(code5,first_time,pdens,oldsd,simelt,phie,fsup,   &
                       firup)

      if(error_type == 3) code3 = 1

      met(iw,ip_tsoil) = stt(nnodes)
      surfemis(iw) = emis

      if(veg_flagl == 1) then
        ft(iw) = ftemp
        surfemisf(iw) = epf
      else
        ft(iw) = mflag
        surfemisf(iw) = mflag
      end if
      tt(iw) = toptemp

      vsmelt = 0d0
      vimelt = 0d0
      refreeze = 0d0
      refreezei = 0d0


! STEP 3a
! calculate upwelling/emitted solar/ir
      if(aint(dabs(met(iw,ip_upsol)-mflag)*1d5)*1d-5 <= eps)            &
        met(iw,ip_upsol) = sup(iw) + fsup                                !- to reflect away from surface
      if(dabs(met(iw,ip_tsol)) <= eps) met(iw,ip_upsol) = 0d0

      if(aint(dabs(met(iw,ip_irup)-mflag)*1d5)*1d-5 <= eps)             &
        met(iw,ip_irup) = irup(iw) + firup


! STEP 4
! calculate soil strength
      if(water_flag == 0)                                               &
        call soil_strength(nstr_flag,surfci(iw),surfrci(iw),surfcbr(iw))                 


! STEP 5
! snow_model calculates the snow depth
      if((hsaccum > eps.or.node_type(nnodes) == 'SN').or.               &
          (aint(met(iw,ip_pt)) == 3.or.aint(met(iw,ip_pt2)) == 3)) then
        if(node_type(nnodes) == 'WA'.and.ice(nnodes) > eps) then
          call snow_model(sprint,water_flag,first_snow,sdensi,hsaccum,  &
                          phie)

          vsmelt = vsmelt*1d-2                                             !m
          if(error_type == 4) code4 = 1
        else if((node_type(nnodes) == 'WA'.and.ice(nnodes) <= eps).or.  &
                                        node_type(nnodes) == 'AI') then
          vsmelt = 0d0
          hsaccum = 0d0

          if(sprint == 1) then
            dates = met(iw,ip_doy) + met(iw,ip_hr)/24d0                 &
                                           + met(iw,ip_min)/(24d0*60d0)
            xtemp = 0d0
            meltfl = '0'
            write(55,91) int(met(iw,ip_year)),dates,xtemp,xtemp,xtemp,  &
                         xtemp,xtemp,xtemp,hsaccum,xtemp,meltfl,xtemp,  &
                         xtemp,xtemp
          end if
        else if(node_type(nnodes) /= 'WA'.and.node_type(nnodes)         &
                                                         /= 'AI') then
          call snow_model(sprint,water_flag,first_snow,sdensi,hsaccum,  &
                          phie)

          if(error_type == 4) code4 = 1
        end if
      else
        vsmelt = 0d0
        hsaccum = 0d0

        if(sprint == 1) then
            dates = met(iw,ip_doy) + met(iw,ip_hr)/24d0                 &
                                           + met(iw,ip_min)/(24d0*60d0)
          xtemp = 0d0
          meltfl = '0'
          write(55,91) int(met(iw,ip_year)),dates,xtemp,xtemp,xtemp,    &
                       xtemp,xtemp,xtemp,hsaccum,xtemp,meltfl,xtemp,    &
                       xtemp,xtemp
        end if
      end if
 91   format(i10,2x,f10.3,8(f15.8),a5,3(f15.8))
      met(iw,ip_sd) = hsaccum
      oldsd = hsaccum

! surface ice thickness model
      if(node_type(nnodes) /= 'WA') then
        if(hi > eps.or.                                                 &
           (aint(met(iw,ip_pt)) == 2.or.aint(met(iw,ip_pt)) == 4)) then
          if(dabs(hsaccum) <= eps) call icethick(oldhi,hi)
          if(error_type == 2) code2 = 1
          oldhi = hi
        end if
      end if 
      met(iw,ip_hi) = hi


! STEP 6
! determine the average surface moisture and ice content
      d1i = 1
      d2i = 2
      surfm = 0d0
      surfi = 0d0
      surfmp = 0d0
      surfip = 0d0
      surfdd = 0d0
      denom = 0d0
      t1 = 0d0
      sfdp = surfdepth + 1d-4

      if(water_flag == 0) then
        j = nnodes
        do while(j >= 2.and.anint((elev-nzi(j))*1d5)*1d-5 <= sfdp)
          if(j == nnodes) then
            num2 = 5d-1*(elev - nzi(j-1))
          else if(anint((elev-nzi(j))*1d5)*1d-5 <= sfdp.and.            &
                          anint((elev-nzi(j-1))*1d5)*1d-5 <= sfdp) then
            num2 = 5d-1*(nz(j+1) - nzi(j-1))
          else if(aint((elev-nzi(j))*1d5)*1d-5 <= sfdp.and.             &
                           anint((elev-nzi(j-1))*1d5)*1d-5 > sfdp) then
            if(anint((elev-nzi(j))*1d5)*1d-5 - sfdp <= eps) then
              num2 = (elev - nzi(j)) - denom
            else
              num2 = 5d-1*((elev - nzi(j)) - surfdepth)
            end if
          end if
          denom = denom + num2
          t1 = nsoilp(j,1)  !/(1d0 - nsoilp(j,2))
          surfdd = surfdd + nsoilp(j,1)*num2

          surfm = surfm + soil_moist(j)*num2
          if(ntype(j) /= 27.and.dabs(t1) > eps) then
            surfmp = surfmp + (soil_moist(j)/t1)*num2
          else
            surfmp = surfm
          end if

          if(ice(j) > eps) then
            dwi = dense(stt(j),0d0,d2i)/dense(stt(j),0d0,d1i)
            surfi = surfi + ice(j)*num2*dwi
            if(ntype(j) /= 27.and.dabs(t1) > eps) then
              surfip = surfip + (ice(j)/t1)*num2*dwi
            else
              surfip = surfi
            end if
          end if

          j = j - 1
        end do
      else
        denom = (elev - nz(nnodes-1))
        surfdd = nsoilp(nnodes,1)*denom

        surfm = 5d-1*(soil_moist(nnodes) + soil_moist(nnodes-1))*denom
        if(ntype(nnodes) /= 27) then
          surfmp = 5d-1*(soil_moist(nnodes) + soil_moist(nnodes-1))     &
                                                *denom/nsoilp(nnodes,1)
        else
          surfmp = surfm
        end if

        if(ice(nnodes) > eps.or.ice(nnodes-1) > eps) then
          dwi = dense(stt(nnodes),0d0,d2i)/dense(stt(nnodes),0d0,d1i)
          dwi1 = dense(stt(nnodes-1),0d0,d2i)/                          &
                                           dense(stt(nnodes-1),0d0,d1i)
          surfi = 5d-1*(ice(nnodes)*dwi + ice(nnodes-1)*dwi1)*denom
          if(ntype(nnodes) /= 27) then
            surfip = 5d-1*(ice(nnodes)*dwi + ice(nnodes-1)*dwi1)*denom  &
                                                      /nsoilp(nnodes,1)
          else
            surfip = surfi
          end if
        end if
      end if

      if(denom > eps) then
        surfmoist(iw) = surfm/denom
        surfice(iw) = surfi/denom
        surfmoistp(iw) = surfmp/denom
        surficep(iw) = surfip/denom
        surfd(iw) = (surfdd/denom)*1d3                                   !kg/m^3
      else
        surfmoist(iw) = soil_moist(nnodes)
        surfice(iw) = ice(nnodes)
        surfmoistp(iw) = soil_moist(nnodes)/nsoilp(nnodes,2)
        surficep(iw) = ice(nnodes)/nsoilp(nnodes,2)
        surfd(iw) = nsoilp(nnodes,1)*1d3   !kg/m^3
      end if

! print the nodal information
      if(nprint == 1) then
        do j=1,nnodes
          k = nnodes - j + 1
          dwi = dense(stt(k),0d0,d2i)/dense(stt(k),0d0,d1i)
!            if(ice(k) > 0d0.or.(soil_moist(k) <= 1.01d0*nsoilp(k,8)     &
!                                          .and.stt(k) <= Tref)) then
          pfrozen = float(frozen)*1d-2*(ice(k)*dwi + soil_moist(k))
          if(ice(k)*dwi >= pfrozen.and.ice(k) > eps) then
            nsstate = 'f'
          else
           nsstate = 'u'
           if(ntype(k) > 19.and.stt(k) <= tref) nsstate = 'f'
          end if

          write(3,1001) int(met(iw,ip_year)),int(met(iw,ip_doy)),       &
                        int(met(iw,ip_hr)),int(met(iw,ip_min)),j,       &
                        node_type(k),elev-nz(k),stt(k),soil_moist(k),   &
                        ice(k),soil_moist(k)+ice(k),wvc(k),nsstate
        end do
      end if
 1001 format(5(i5),3x,a2,1x,2(f9.3),4(e12.4),3x,a1)

! determine surface freeze or thaw depth and state
      sstate(iw) = 'u'
      frthick(iw) = elev - nz(nnodes)
      twthick(iw) = elev - nz(nnodes)
      dwi = dense(stt(nnodes),0d0,d2i)/dense(stt(nnodes),0d0,d1i)
      pfrozen = float(frozen)*1d-2*(ice(nnodes)*dwi                     &
                                                  + soil_moist(nnodes))
!        pfrozen = float(frozen)*1d-2*ice(nnodes)                        &
!                                    /(ice(nnodes) + soil_moist(nnodes))
!      if(node_type(nnodes) == 'WA') pfrozen = 0.99d0

      fts = 0
      if(ice(nnodes)*dwi >= pfrozen.and.ice(nnodes) > eps) then
        sstate(iw) = 'f'
        do k=nnodes-1,1,-1
          dwi = dense(stt(k),0d0,d2i)/dense(stt(k),0d0,d1i)
          pfrozen = float(frozen)*1d-2*(ice(k)*dwi + soil_moist(k))
!            pfrozen = float(frozen)*1d-2*ice(k)                         &
!                                              /(ice(k) + soil_moist(k))
   
          if((ice(k)*dwi >= pfrozen.and.ice(k) > eps).and.fts == 0) then
            if(k >= 2) then
              frthick(iw) = (nz(k) + nz(k-1))*5d-1
              frthick(iw) = elev - frthick(iw)
            else if(k == 1) then
              twthick(iw) = 0d0
              frthick(iw) = elev - nz(k)
            end if
          else
            if(k == nnodes-1)                                           &
              frthick(iw) = elev - (nz(k) + nz(k+1))*5d-1
            fts = 1
            j = k
            do while (j >= 2.and.(ice(j) < pfrozen.or.                  &
                                                  dabs(ice(j)) <= eps))
              twthick(iw) = (nz(j) + nz(j-1))*5d-1
              twthick(iw) = elev - twthick(iw)
              j = j - 1
              dwi = dense(stt(j),0d0,d2i)/dense(stt(j),0d0,d1i)
              pfrozen = float(frozen)*1d-2*(ice(j)*dwi + soil_moist(j))
!                pfrozen = float(frozen)*1d-2*ice(j)                     &
!                                               /(ice(j) + soil_moist(j))

            end do

            if(j == 1.and.(ice(j)*dwi < pfrozen.or.ice(j) <= eps))      & 
              twthick(iw) = elev - nz(1)
            exit
          end if
        end do
      else
        if(iw /= 1) then
          if(sstate(iw-1) == 'f') sstate(iw) = 't'
          if(sstate(iw-1) == 't') sstate(iw) = 'u'
        end if

        twthick(iw) = (nz(nnodes) + nz(nnodes-1))*5d-1
        twthick(iw) = elev - twthick(iw)

        do k=nnodes-1,1,-1
          dwi = dense(stt(k),0d0,d2i)/dense(stt(k),0d0,d1i)
          pfrozen = float(frozen)*1d-2*(ice(k)*dwi + soil_moist(k))
!            pfrozen = float(frozen)*1d-2*ice(k)                         &
!                                              /(ice(k) + soil_moist(k))

          if((ice(k)*dwi < pfrozen.or.dabs(ice(k)) <= eps).and.         &
                                                         fts == 0) then
            if(k >= 2) then
              twthick(iw) = (nz(k) + nz(k-1))*5d-1
              twthick(iw) = elev - twthick(iw)
            else if(k <= 1) then
              frthick(iw) = 0d0
              twthick(iw) = elev - nz(k)
            end if
          else
            if(k == nnodes-1)                                           &
              twthick(iw) = elev - (nz(k) + nz(k+1))*5d-1
            fts = 1
            j = k
            do while (j >= 2.and.(ice(j) >= pfrozen.and.ice(j) > eps))
              frthick(iw) = (nz(j) + nz(j-1))*5d-1
              frthick(iw) = elev - frthick(iw)
              j = j - 1
              dwi = dense(stt(j),0d0,d2i)/dense(stt(j),0d0,d1i)
              pfrozen = float(frozen)*1d-2*(ice(j)*dwi + soil_moist(j))
!                pfrozen = float(frozen)*1d-2*ice(j)                     &
!                                               /(ice(j) + soil_moist(j))
            end do
            if(j == 1.and.ice(j)*dwi >= pfrozen)                        &
              frthick(iw) = elev - nz(1)
            exit
          end if
        end do
      end if
!      twthick(iw) = anint(twthick(iw)*1d2)*1d-2
!      frthick(iw) = anint(frthick(iw)*1d2)*1d-2
      if(node_type(nnodes) == 'WA') met(iw,ip_hi) = frthick(iw)

      error_type = 0

!      if(dstart /= 2) then
!        if((iw == iend-(dstart-1).or.iw == iend-(dstart-1)+4).or.     &
!                                                       iw == iend) then
!          d1i = 2
!          call read_old_data(d1i,iw,sdensi,phie,timeo)
!        end if
!      else if(dstart == 2) then

      wstart = max(iend-(moverlap-1),istart)
      if(iw >= wstart) then
        d1i = 2
        call read_old_data(d1i,iw,wstart,sdensi,phie,timeo)
      end if
!      end if
 
! !!!!!!!!!!!!!!!!!!!!!!!!end of main calculation loop !!!!!!!!!!!!!!!!!!!!!!!!!

      end subroutine fasst_main

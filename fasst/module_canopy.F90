module module_canopy

      use fasst_global

   contains

      subroutine canopy_met(oldsd)

!       @(#)aavctm.f  2.1   12/8/95   2.1

!         Dr. James A. Smith
!         Code 920
!         NASA Goddard Space Flight Center
!         Greenbelt, MD 20771

!         Telephone:  (301) 614-6020  
!         E-mail:     jasmith@hemlock.gsfc.nasa.gov

!     Modified by  11 Jan 2000
!         Dr. George Koenig
!         ERDC/CEERD-RF
!         72 Lyme Rd
!         Hanover, NH 03755
!
!         e-mail gkoenig@crrel.usace.army.mil
!         (603) 646-4556
         
!!          James R Jones  12 March 2001
!!          added year output to CRREL VIEW type data output for multi-year models
!!       JRJ  09/14/01 Accomodate file names in strings

!     modified by Dr. Susan Frankenstein Jan 2004 to incorporate into FASST
!       called from FASST subroutine 'fasst.f'

! ******************************************************************************
!  VERSION 2.1 (Beta)
!     This is the Beta 2.1 version of the updated thermal model.  
!     It is undergoing rapid modification, particularly with respect
!       to the sensible and latent heat formulations.

!     This version includes:

!      -Generalization to n-layers
!      -A wind profile calculation, appropriate for tree canopies 
!      -The computation of the Sij matrix necessary for long and short-wave
!          flux transfer calculations
!      -Incorporation of a canopy geometric clumping factor 
!      -The computation of average shortwave flux absorbed in the canopy a
!          using the radiosity approximation
!
!      -The general approach is to pre-calculate the geometrical dependences,
!       for example the imapct of the canopy structure on solar and infrared 
!       fluxes

!  Purpose:
!       Calculate temperature profiles and thermal infrared directional exitance 
!       for a general n-layer canopy
!
!       The model will also be linked with the CEERD soil model.  This model will be
!       used to calculate the radiative fluxes under the canopy.  The model will
!       run iteratively with the soil model. 

!  Input:

!     Input files required are:
!       Meteorological file:    
!     Output files required:
!       Arbitrary file for results:  e.g. vctm_output.
 
!  Output:
!       Energy Flux Terms, Canopy Layer Temperatures 
     
! ------------------------------------------------------------------------------
!    References: 

!    Smith, J. A., et al. 1981. "Thermal Vegetation Canopy Model Stuidies," 
!       Remote Sensing of Environment, 11:311-326.

!    Smith, J. A. and Goltz, S. M. (1994), "A Thermal Exitance and Energy
!      Balance Model for Forest Canopies," 
!       IEEE Trans. Geosciences and Rem. Sensing.
! ------------------------------------------------------------------------------
!  Parameter              Remarks
!  nclayers         number of canopy layers
!  max_thetar       number of view angles for integration over the hemisphere to
!                  obtain diffuse and IR flux
!  dzveg           canopy layer thickness
!  sabs            canopy layer shortwave absorption
!  epsc            canopy layer longwave emissivity
!  rho             canopy layer reflectance 
!  ta              air temperature (C)
!  tg              ground temperature (C)
! ------------------------------------------------------------------------------

! calls the following subroutines:
!     canopy_prop (attached to this subroutine)
!     scalc_2a (attached to this subroutine)
!     roughness (attached to this subroutine)
!     windprofile (attached to this subroutine)
!     radiosity (attached to this subroutine)
!     canopy_moist (attached to this subroutine)
!     feval (attached to this subroutine)
!     solve (attached to this subroutine)
!     output (attached to this subroutine)


      implicit none


! NOTE: nclayers is defined in fasst_global.F90. It is currently set to 3.

      real(kind=8),intent(in):: oldsd

! saved variables
      real(kind=8):: x(nclayers),epsc(nclayers)
      real(kind=8):: rho(nclayers),clump(nclayers),tau(nclayers)
      real(kind=8):: alp(nclayers),psi(nclayers+1)
      save:: x,epsc,rho,clump,tau,alp,psi

! local variables
! maximum number of viewing directions
      integer(kind=4),parameter:: max_thetar = 30
! n_thetar = number of zenith angles used to integrate over the sources of emitted and scattered radiation
      integer(kind=4),parameter:: n_thetar = 9
! maximum number of leaf angles
!      integer(kind=4),parameter:: max_thetak = 30
! n_thetak = number of leaf angles used to calculate leaf normal distribution of emitted and scattered radiation
      integer(kind=4),parameter:: n_thetak = 9
! maximum number of source azimuth angles
!      integer(kind=4),parameter:: max_phi = 30
! n_phi = number of azimuth angles used to integrate over the sources of emitted and scattered radiation
!      integer(kind=4),parameter:: n_phi = 18


      integer(kind=4):: i,iterations,ncsnow,season,icall,ll,j,ii,ij
      integer(kind=4):: dzflag,foliage_type(nclayers)                    !leaf angle type, max = 6
      integer(kind=4),parameter:: max_iter = 10

      real(kind=8):: sumdx,zo,zos,zdc,fmu,fmuh,glb,wscorr,fmuh1,aground
      real(kind=8):: zu,bta,btg,ta,tg,rh,epstg,sdepth,lw,htemp,ftg,laim
      real(kind=8):: bgr,sum_lai,sum_sai,tsign,rain,snow,tot_ep,dtop
      real(kind=8):: dx(nclayers),sabs(nclayers),fx(nclayers)
      real(kind=8):: a(nclayers),sigfhi(nclayers),bx(nclayers)
      real(kind=8):: rx(nclayers),qx(nclayers),rl(nclayers)
      real(kind=8):: wind_prof(nclayers),lai(nclayers),ful(nclayers)
      real(kind=8):: ra(nclayers),dzveg1(nclayers),sai(nclayers)
      real(kind=8):: ef(nclayers),def(nclayers),chf(nclayers)
      real(kind=8):: rh_prof(nclayers),tai(nclayers),pheat(nclayers)
      real(kind=8):: dpheat(nclayers)
      real(kind=8):: xco(nclayers),stcl(nclayers),stcs(nclayers)
      real(kind=8):: dfx(nclayers,nclayers)
      real(kind=8):: Sij(nclayers+2,nclayers+2)
      real(kind=8):: Wir(nclayers+2,max_thetar),ttemp,temp,sumt

      real(kind=8),parameter:: tol = 1d-2  !1d-4                         !error tolerance

! zero-out variables
      iterations = 0
      ncsnow = 0
      icall = 1
      ll = 0
      dzflag = 0
      zos = 0d0
      fmu = 0d0
      fmuh = 0d0
      glb = 0d0
      bta = 0d0
      btg = 0d0
      ta = 0d0
      tg = 0d0
      rh = 0d0
      epstg = 0d0
      sdepth = 0d0
      htemp = 0d0
      lw = 0d0
      ftg = 0d0
      laim = 0d0
      wscorr = 0d0
      fmuh1 = 0d0
      zu = 0d0
      rain = 0d0
      snow = 0d0
      temp = 0d0
      tot_ep = 0d0
      dtop = 0d0

      do i=1,nclayers
        sabs(i) = 0d0
        a(i) = 0d0
        lai(i) = 0d0
        ful(i) = 0d0
        dzveg1(i) = 0d0
        sai(i) = 0d0
      end do

      ttemp = toptemp
      toptemp = anint(toptemp*1d20)*1d-20

      call albedo_emis(oldsd,dmet1(iw,1),dmet1(iw,8))

      call canopy_prop(season,foliage_type,epsc,alp,rho,tau,clump,psi)

! Read met data (time, air tmp, grnd tmp, wind speed, rh, total solar, IR)
      fmu = dmax1(0.1d0,dmet1(iw,3))                                     !wind speed (m/s)
      ta = anint((dmet1(iw,4) -Tref)*1d20)*1d-20                         !air temp (C)
      tg = anint((toptemp - Tref)*1d20)*1d-20                            !surface temp (C)
      rh = dmet1(iw,5)                                                   !relative humidity (decimal)
      glb = dmet1(iw,1)                                                  !total incoming solar (W/m^2)
      lw = dmet1(iw,2)                                                   !incoming IR (W/m^2)
      bta = anint((lw/sigma)*1d20)*1d-20                                 !normalized by St.Bolz. constant
      if(aint(dabs(met(iw,ip_irup)-mflag)*1d5)*1d-5 <= eps) then
        if(ttemp > eps.and.ftemp > eps) then
          btg = (1d0 - sigfl)*((1d0 - emis)*dmet1(iw,2)                 & 
                                               + emis*sigma*ttemp**4d0) &
                                      +  sigfl*((1d0 - epf)*dmet1(iw,2) &
                                                + epf*sigma*ftemp**4d0)  !W/m^2 (-)
          btg = anint(btg*1d20)*1d-20
        end if                             
      else
        btg = met(iw,ip_irup)
      end if
      btg = anint((btg/sigma)*1d20)*1d-20

      fmuh = fmu

      if(iw == istart) then
        do i=1,nclayers
          x(i) = 0d0
          x(i) = 0.9d0*ta                                                !Initialize the layers to the air temperature (C)
          storcl(i) = 0d0                                                !initialize water storage on the leaves (m)
          storcs(i) = 0d0                                                !Initialize snow storage on the leaves
          if(infer_test == 0) avect(0,i) = x(i)
        end do

        if(infer_test == 1) then
          do i=1,nclayers
            x(i) = canopy_temp(i,oldpos) - Tref                          !C
            storcl(i) = stor(i)
            storcs(i) = stor(i+nclayers)
          end do
        end if
      end if

      if(iw >= 2) then
        do i=1,nclayers
          x(i) = canopy_temp(i,iw-1) - Tref                              !C
        end do
      end if

      do i=1,nclayers
        xco(i) = x(i)
      end do

      zu = zh + iheight                                                  !wind measurement height (m)
      sdepth = hsaccum + hi                                              !combined snow and ice thickness (m)

      htemp = 0d0
      if(hfol_tot >= sdepth) then
        htemp = hfol_tot
      else
        htemp = sdepth
      end if

      if(sdepth == 0d0) then                                             !no snow on ground
        if(vegl_type == 0) then                                          !no low vegetation
          zos = rough                                                    !ground roughness height (m)
        else                                                             !yes low vegetation
          zos = z0l                                                      !low veg roughness height (m) [low_veg_prop.f]
        end if
      else                                                               !yes snow on ground
        if(vegl_type == 0.or.(hfol_tot <= sdepth.or.                    &
                              dabs(hfol_tot-sdepth) < 1d-10)) then       !no low vegetation or snow depth more than foliage height
          zos = 7.775d-3                                                 !snow roughness height (m) (0.05 - 1.5mm)
        else                                                             !yes low vegetation, snow less than foliage height
          zos = z0l                                                      !low veg roughness height (m) [low_veg_prop.f]
        end if
      end if

!      ftg = 1d0 - 1.6d-4*(298d0 - toptemp)*(298d0 - toptemp)
      ftg = 1d0 - 1.6d-3*(298d0 - stt(refn))*(298d0 - stt(refn))
      ftg = anint(ftg*1d20)*1d-20

      laim = veg_prp(vegh_type,6) + ftg*(veg_prp(vegh_type,5)           & 
                                                - veg_prp(vegh_type,6))
      if(laim > veg_prp(vegh_type,5)) laim = veg_prp(vegh_type,5)
      if(laim < veg_prp(vegh_type,6)) laim = veg_prp(vegh_type,6)
      laim = anint(laim*1d20)*1d-20

! sigfh calculation from Ramirez and Senarath (2000), J. Climate, 13(22), p.4050-
      sigfh = veg_prp(vegh_type,2) - (1d0 - ftg)                        &
                  *(veg_prp(vegh_type,2) - veg_prp(vegh_type,3))
      if(sigfh < veg_prp(vegh_type,2)) sigfh = veg_prp(vegh_type,2)
      if(sigfh > veg_prp(vegh_type,3)) sigfh = veg_prp(vegh_type,3)
      sigfh = sigfh*1d-2
      if(season == iseason.and.dabs(isigfh-spflag) > eps)               &
        sigfh = isigfh
      if(dabs(isigfh-spflag) > eps.and.(sigfh > isigfh.and.             &
         iseason >= season)) sigfh = isigfh

      sigfh = anint(sigfh*1d20)*1d-20

      ncsnow = nclayers
      do i=1,nclayers
        lai(i) = laim*laif(i)                                            !partition between the layers
        if((iw == istart.and.infer_test == 0).and.ilai(i) /= spflag)    &
          lai(i) = ilai(i)
        lai(i) = anint(lai(i)*1d20)*1d-20

        dzveg1(i) = dzveg(i)
        if(i == nclayers) then
          dzveg1(i) = dzveg(i) - htemp
          if(dzveg1(i) <= eps) then
            ncsnow = ncsnow - 1
            dzveg1(i) = 0d0
            if(i /= 1) then
              dzveg1(i-1) = dzveg1(i-1) - (htemp - dzveg(i))
              if(dzveg1(i-1) < 0d0) then
                ncsnow = ncsnow - 1
                dzveg1(i-1) = 0d0
                dzveg1(i-2) = dzveg1(i-2) -                             &
                              (htemp - (dzveg(i) + dzveg(i-1)))
                if(dzveg1(i-2) < 0d0) then
                  ncsnow = ncsnow - 1
                  dzveg1(i-2) = 0d0
                end if
              end if   !if(dzveg1(i-1) < 0d0) then
            end if   !if(i /= 1) then
          end if
         end if
        dzveg1(i) = anint(dzveg1(i)*1d20)*1d-20
        if(dzveg1(i) <= eps) dzflag = 1
        sai(i) = veg_prp(vegh_type,7)*(1d0 - laif(i))
        sai(i) = anint(sai(i)*1d20)*1d-20
      end do

      if(ncsnow /= 0) then
! determine surface (ground) material, set parameters accordingly
!        psi(ncsnow+1) = (1d0 - sigfl)*albedo + sigfl*albf                !surface albedo, ground reflectance
        psi(ncsnow+1) = (1d0 - sigfl)*psi(nclayers+1) + sigfl*albf
        epstg = emis                                                     !surface emissivity

        call scalc_2a(ncsnow,max_thetar,n_thetar,n_thetak,foliage_type, &
                     lai,sai,clump,Sij,Wir)

! calculates displacement hgt roughness length in the canopy
        sum_lai = 0d0
        sum_sai = 0d0
        do i = 1,ncsnow
          sum_lai = sum_lai + lai(i)
          sum_sai = sum_sai + sai(i)
        end do

        call roughness(sum_lai,sum_sai,zh,htemp,zdc,zo)

! calculates canopy windprofile and rh profile
        call windprofile(ncsnow,zh,glb,sigfh,dzveg1,wind_prof,rh_prof,  &
                         tai)

! Use the wind scaling profile factor calculated in windprofile to obtain the variation
! of the wind speed within the canopy
        wscorr = dlog((zh - zdc)/zo)/dlog((zu - zdc)/zo)                 !log correction factor for wind
        fmuh = fmu*wscorr                                                !windspeed at top of canopy (cm/s)
        do i=1,ncsnow
          ful(i) = anint(fmuh*wind_prof(i)*1d20)*1d-20

          rh_prof(i) = 1d2 + (rh - 1d2)*rh_prof(i)
          if(rh_prof(i) > 1d2) rh_prof(i) = 1d2
          if(rh_prof(i) < 1d-1) rh_prof(i) = 1d-1
          rh_prof(i) = anint(rh_prof(i)*1d20)*1d-20

          tai(i) = anint(ta*tai(i)*1d20)*1d-20                             !C
        end do

        call radiosity(ncsnow,Sij,psi,sabs,aground)

        sumdx = 0d0
        do i=1,ncsnow
          a(i) = dmax1(0d0,glb*sabs(i))                                  !Short wave absorbed energy, glb shortwave at the top
                                                                         !of the canopy
          a(i) = anint(a(i)*1d20)*1d-20
          sumdx = sumdx + 1.1d0*tol
        end do

        do while(dabs(sumdx) > tol)
          call canopy_moist(ncsnow,rain,snow,tot_ep,lai,rl,ra,x,ful,ef, &
                            def,chf,sigfhi,rh_prof,tai,pheat,dpheat,    &
                            stcl,stcs)

          call feval(ncsnow,icall,Sij,x,epsc,alp,bta,btg,bx,a,qx,rx,ef, &
                     def,chf,sigfhi,fx,dfx,tai,pheat,dpheat)

          call solve(ncsnow,fx,dfx,dx)

          sumdx = 0d0
          dtop = 2.5d0 !1d1 !2.5d0 
!          if(iterations == 0) dtop = 1d-1*dtop
          do i=1,ncsnow
            tsign = 0d0
            if(dabs(dx(i)) > eps) tsign = anint(dabs(dx(i))/dx(i))
!            if(dabs(dx(i)) > dtop) dx(i) = tsign*dtop                    !C
            x(i) = x(i) + dx(i)                                          !C
            x(i) = anint(x(i)*1d15)*1d-15

            sumdx = sumdx + dabs(dx(i))
          end do

          iterations = iterations + 1
          if(iterations >= max_iter) exit
        end do  !do while(dabs(sumdx) > tol)

        if(iterations >= max_iter) then
          do i=1,ncsnow
            x(i) = tai(i)
          end do
          error_code = 1
          error_type = 2
!         if(single_multi_flag == 0) write(10,'(''freq_id'',i10,'' Too m&     
!           &any iterations in canopy_met; day'',i4,'' hour'',i4)')     &
!           freq_id,int(met(iw,ip_doy)),int(met(iw,ip_hr))
        end if

        call feval(ncsnow,icall,Sij,x,epsc,alp,bta,btg,bx,a,qx,rx,ef,   &
                   def,chf,sigfhi,fx,dfx,tai,pheat,dpheat)        

        icall = 0
!        call output(ncsnow,max_thetar,Sij,bx,btg,bta,Wir,bgr)
        call output(ncsnow,Sij,bx,btg,bta,bgr)
      end if  !if(ncsnow /= 0) then

      if(ncsnow < nclayers) then
        if(ncsnow == 0) then
          aground = 1d0
          bgr = lw
        end if
        do i=ncsnow+1,nclayers
          storcl(i) = 0d0
          storcs(i) = 0d0

          if(hfol_tot >= sdepth.and.hfol_tot > eps) then
            x(i) = ftemp - Tref                                          !C
          else if(sdepth > eps) then
            x(i) = stt(nnodes+1) - Tref   !0d0
          else
            x(i) = ta
          end if
        end do
      else if(dzflag == 1) then
        do i=1,ncsnow
          if(dzveg1(i) <= eps) then
            storcl(i) = 0d0
            storcs(i) = 0d0
            x(i) = tai(i) - Tref
          end if
        end do
      end if

      ll = 2
      if(timstep >= 2) ll = 1
      do i=1,nclayers
        storcl(i) = anint(stcl(i)*1d20)*1d-20
        storcs(i) = anint(stcs(i)*1d20)*1d-20

        if(iw > ll.or.(iw == ll.and.infer_test == 1)) then
          j = 0
          temp = avect(j+1,i)
          do while (j < ll)
            avect(j,i) = temp
            j = j + 1
            if(j <= ll-1) then
              temp = avect(j+1,i)
            end if
          end do
          avect(ll,i) = x(i)
          ii = ll
          j = ll
          ij = j
        else
          j = iw
          avect(j,i) = x(i)
          ii = iw
          j = iw
          ij = j
        end if

        sumt = 0d0
        do while(j > -1)
          sumt = sumt + avect(j,i)
          j = j - 1
        end do

        x(i) = sumt/float(ii+1)
        x(i) = anint(x(i)*1d10)*1d-10

        canopy_temp(i,iw) = anint((x(i) + Tref)*1d10)*1d-10              !K
        stor(i) = storcl(i)
        stor(i+nclayers) = storcs(i)
      end do

! Use canopy calculated effects to estimate air temperature, solar, IR, wind speed
      dmet1(iw,1) = aground*glb                                          !total solar above ground layer (W/m^2)
      if(dmet1(iw,1) < 1d-2) dmet1(iw,1) = 0d0
      dmet1(iw,2) = bgr                                                  !IR above ground layer (W/m^2)

      if(ncsnow /= 0) then
        dmet1(iw,3) = ful(ncsnow)                                        !wind speed above ground (m/s)        
        dmet1(iw,4) = tai(ncsnow) + Tref                                 !air temp above ground (K)
        dmet1(iw,5) = rh_prof(ncsnow)                                    !relative humidity above ground (%)
      end if

      dmet1(iw,6) = rain
      dmet1(iw,7) = snow
      if(dmet1(iw,6) <= eps.and.aint(met(iw,ip_pt)) /= 1)               &
                                               met(iw,ip_pt) = float(1)
      if(dmet1(iw,6) > eps.and.aint(met(iw,ip_pt)) /= 2)                &
                                               met(iw,ip_pt) = float(2)
      if(dmet1(iw,7) <= eps.and.aint(met(iw,ip_pt2)) /= 1)              &
                                              met(iw,ip_pt2) = float(1)
      if(dmet1(iw,7) > eps.and.aint(met(iw,ip_pt2)) /= 3)               &
                                              met(iw,ip_pt2) = float(3)

      if(met(iw,ip_tsol) > eps.and.dmet1(iw,1) > eps) then
        if(aint(dabs(dmet1(iw,8)-mflag)*1d5)*1d-5 > eps)                & 
          dmet1(iw,8) = dmet1(iw,8)*dmet1(iw,1)/met(iw,ip_tsol)          !upwelling solar (W/m^2)

        dmet1(iw,9) = dmet1(iw,9)*dmet1(iw,1)/met(iw,ip_tsol)            !direct solar (W/m^2)
        dmet1(iw,10) = dmet1(iw,10)*dmet1(iw,1)/met(iw,ip_tsol)          !diffuse solar (W/m^2)
      else
        if(aint(dabs(dmet1(iw,8)-mflag)*1d5)*1d-5 > eps)                &
          dmet1(iw,8) = 0d0

        dmet1(iw,9) = 0d0
        dmet1(iw,10) = 0d0
      end if

      do i=1,13
        dmet1(iw,i) = anint(dmet1(iw,i)*1d15)*1d-15
      end do

      sumdx = 0d0
      do i=1,ncsnow
        sumdx = sumdx + epsc(i)
      end do
      if(ncsnow /= 0) then
        surfemisc(iw) = sumdx/float(ncsnow)
      else
        surfemisc(iw) = sumdx
      end if

      toptemp = ttemp
      cevap(iw) = tot_ep                                                 !total canopy evaporation/condensation (kg/m^2*s)

      end subroutine canopy_met

! ******************************************************************************
      subroutine canopy_prop(season,foliage_type,epsc,alp,rho,tau,clump,&
                             psi)

!  Purpose: Set simulation variables
!     nclayers    !no. of vegetation layers (default = 3)
!     zh          !canopy height (m)
!     (foliage_type(j),j=1,nclayers)
!     (lai(j),j=1,nclayers)        !leaf area index
!     (clump(j),j=1,nclayers)      !Markov clumping factors (typically 1.0)
! NOTE: clump = 1: homogeneous, randomly placed
!       clump = 0: completly clumped 
!     (rho(j),j=1,nclayers)        !layer 1-3 reflectance
!     (tau(j),j=1,nclayers)        !layer 1-3 transmission
!     (alp(j),j=1,nclayers)        !layer 1-3 shortwave absorption
!     (epsc(j),j=1,nclayers)       !layer 1-3 longwave emissivity
! NOTE: layer 1 = top, 2 = middle, 3 = trunk

! no subroutines called

      implicit none

      integer(kind=4),intent(out):: season,foliage_type(nclayers)
      real(kind=8),intent(inout):: epsc(nclayers),alp(nclayers)
      real(kind=8),intent(inout):: rho(nclayers),tau(nclayers)
      real(kind=8),intent(inout):: clump(nclayers),psi(nclayers+1)

! local variables
      integer(kind=4),parameter:: nc1 = 3
      integer(kind=4):: i,j,foliagetype,vind,vind1

      real(kind=8):: ftg,rho_max,tau_min,m,df_alp(nc1,nc1)
      real(kind=8):: rho_min(nc1),tau_max(nc1),mclump(nc1),dzveg1o(5)
      real(kind=8):: eps_old(nclayers),alp_old(nclayers)
      real(kind=8):: rho_old(nclayers),tau_old(nclayers)
      real(kind=8):: clump_old(nclayers),psi_old(nclayers+1)


! transmittance and reflectance values from Dorman, J.L. & P.J. Sellers (1989)
!                                           J. Appl. Met., V28, pp.833-855.
! leaf type based on the above reference
! absorptions are based on SWOE report
      data foliagetype /6/
      data mclump /0.5d0, 0.5d0, 1d0/                                    !top, middle, bottom
      data df_alp /0.229d0, 0.214d0, 0.079d0,                           & !needle
                   0.255d0, 0.046d0, 0.038d0,                           & !broad
                   0.242d0, 0.130d0, 0.058d0/                             !mixed
      data rho_min /0.21d0, 0.26d0, 0.24d0/                              !needle, broad, mixed
      data rho_max /0.28d0/
      data tau_min /0.001d0/
      data tau_max /0.150d0, 0.150d0, 0.100d0/                           !needle, broad, mixed
      data dzveg1o /11.2d0,9d0,13.7d0,16.3d0,13.5d0/                     !vegh_type 3,4,5,6,18


! zero-out variables
      season = 0
      vind = 0
      vind1 = 0
      ftg = 0d0
      m = 0d0

      do i=1,nclayers
        eps_old(i) = 0d0
        alp_old(i) = 0d0
        rho_old(i) = 0d0
        tau_old(i) = 0d0
        clump_old(i) = 0d0

        eps_old(i) = epsc(i)
        alp_old(i) = alp(i)
        rho_old(i) = rho(i)
        tau_old(i) = tau(i)
        clump_old(i) = clump(i)

        foliage_type(i) = 0
        epsc(i) = 0d0
        alp(i) = 0d0
        rho(i) = 0d0
        tau(i) = 0d0
        clump(i) = 0d0
      end do

      do i=1,nclayers+1
        psi_old(i) = 0d0
        psi_old(i) = psi(i)
        psi(i) = 0d0
      end do

      if(albedo > eps) then
        psi(nclayers+1) = albedo
      else
        if(psi_old(nclayers+1) > eps) then
          psi(nclayers+1) = psi_old(nclayers+1)
        else
         psi(nclayers+1) = nsoilp(nnodes,3)
        end if
      end if

! determine season
      if(lat >= 0d0) then                                               !northern hemisphere
        if(met(iw,ip_doy) >= 335.or. met(iw,ip_doy) <=  80) season = 1   !winter
        if(met(iw,ip_doy) >=  81.and.met(iw,ip_doy) <= 151) season = 2   !spring 
        if(met(iw,ip_doy) >= 152.and.met(iw,ip_doy) <= 243) season = 3   !summer
        if(met(iw,ip_doy) >= 244.and.met(iw,ip_doy) <= 334) season = 4   !fall
      else                                                              !southern hemisphere
        if(met(iw,ip_doy) >= 335.or. met(iw,ip_doy) <=  80) season = 3   !summer
        if(met(iw,ip_doy) >=  81.and.met(iw,ip_doy) <= 151) season = 4   !fall 
        if(met(iw,ip_doy) >= 152.and.met(iw,ip_doy) <= 243) season = 1   !winter
        if(met(iw,ip_doy) >= 244.and.met(iw,ip_doy) <= 334) season = 2   !spring
      end if

      if(iw == istart.and.infer_test == 0) iseason = season

!      ftg = 1d0 - 1.6d-4*(298d0 - toptemp)*(298d0 - toptemp)
      ftg = 1d0 - 1.6d-3*(298d0 - stt(refn))*(298d0 - stt(refn))

      if(vegh_type == 3.or.vegh_type == 6) then                          !evergreen needle & broad leaf
        vind = 1
        vind1 = 1
        if(vegh_type == 6) then
          vind = 2
          vind1 = 4
        end if

        do j=1,nclayers
          foliage_type(j) = foliagetype
          clump(j) = mclump(j)
          if(j == nclayers) then
            rho(j) = rho_max
            tau(j) = tau_min
            epsc(j) = veg_prp(vegh_type,12)
            alp(j) = veg_prp(vegh_type,14)
          else
            rho(j) = rho_min(vind)
            tau(j) = tau_max(vind)
            epsc(j) = veg_prp(vegh_type,13)
            if(j == 1) then
              alp(j) = veg_prp(vegh_type,15)
            else
              m = (df_alp(1,vind) - df_alp(2,vind))*dzveg(1)/           &
                                                         dzveg1o(vind1)
!              alp(j) = anint((veg_prp(vegh_type,15)-m)*1d3 + 5d-1)*1d-3
              alp(j) = anint((veg_prp(vegh_type,15)-m)*1d20)*1d-20
            end if
          end if
        end do

      else if(vegh_type == 4.or.vegh_type == 5) then                     !deciduous needle & broad leaf
        vind = 1
        vind1 = 2
        if(vegh_type == 5) then
          vind = 2
          vind1 = 3
        end if

        do j=1,nclayers
          foliage_type(j) = foliagetype
          clump(j) = mclump(j)
          if(j == nclayers) then
            rho(j) = rho_max
            tau(j) = tau_min
            epsc(j) = veg_prp(vegh_type,12)
            alp(j) = veg_prp(vegh_type,14)
          else
            rho(j) = rho_min(vind) + (1d0 - ftg)*                       &
                                               (rho_max - rho_min(vind))
            if(rho(j) < rho_min(1)) rho(j) = rho_min(1)
            if(rho(j) > rho_max) rho(j) = rho_max

            tau(j) = tau_max(vind) - (1d0 - ftg)*                       &
                                               (tau_max(vind) - tau_min)
            if(tau(j) < tau_min) tau(j) = tau_min
            if(tau(j) > tau_max(vind)) tau(j) = tau_max(vind)

            epsc(j) = veg_prp(vegh_type,13) - (1d0 - ftg)*              &
                        (veg_prp(vegh_type,13) - veg_prp(vegh_type,12))
            if(epsc(j) < veg_prp(vegh_type,12))                         &
              epsc(j) = veg_prp(vegh_type,12)
            if(epsc(j) > veg_prp(vegh_type,13))                         &
              epsc(j) = veg_prp(vegh_type,13)

            if(j == 1) then
              alp(j) = veg_prp(vegh_type,15)
            else
              m = (df_alp(1,vind) - df_alp(2,vind))*dzveg(1)/           &
                                                         dzveg1o(vind1)
!              alp(j) = anint((veg_prp(vegh_type,15)-m)*1d3 + 5d-1)*1d-3
              alp(j) = anint((veg_prp(vegh_type,15)-m)*1d20)*1d-20
            end if
          end if
        end do

      else if(vegh_type == 18) then                                      !mixed
        vind = 3
        vind1 = 5

        do j=1,nclayers
          foliage_type(j) = foliagetype
          clump(j) = mclump(j)
          if(j == nclayers) then
            rho(j) = rho_max
            tau(j) = tau_min
            epsc(j) = veg_prp(vegh_type,12)
            alp(j) = veg_prp(vegh_type,14)
          else
            rho(j) = rho_min(3) + (1d0 - ftg)*(rho_max - rho_min(3))
            if(rho(j) < rho_min(3)) rho(j) = rho_min(3)
            if(rho(j) > rho_max) rho(j) = rho_max

            tau(j) = tau_max(3) - (1d0 - ftg)*(tau_max(3) - tau_min)
            if(tau(j) < tau_min) tau(j) = tau_min
            if(tau(j) > tau_max(3)) tau(j) = tau_max(3)

            epsc(j) = veg_prp(vegh_type,13) - (1d0 - ftg)*              &
                        (veg_prp(vegh_type,13) - veg_prp(vegh_type,12))
            if(epsc(j) < veg_prp(vegh_type,12))                         &
              epsc(j) = veg_prp(vegh_type,12)
            if(epsc(j) > veg_prp(vegh_type,13))                         &
              epsc(j) = veg_prp(vegh_type,13)
            if(j == 1) then
              alp(j) = veg_prp(vegh_type,15)
            else
              m = (df_alp(1,vind) - df_alp(2,vind))*dzveg(1)/           &
                                                         dzveg1o(vind1)
!              alp(j) = anint((veg_prp(vegh_type,15)-m)*1d3 + 5d-1)*1d-3
              alp(j) = anint((veg_prp(vegh_type,15)-m)*1d20)*1d-20
            end if
          end if
        end do
      end if

      if(season == iseason) then
        do i=1,nclayers
          if(ifoliage_type(i) /= spflag)                                & 
               foliage_type(i) = int(ifoliage_type(i))
          if(iclump(i) /= spflag) clump(i) = iclump(i)
          if(irho(i) /= spflag) rho(i) = irho(i)
          if(itau(i) /= spflag) tau(i) = itau(i)
          if(ialp(i) /= spflag) alp(i) = ialp(i)
          if(ieps(i) /= spflag) epsc(i) = ieps(i)
        end do
      end if

! The leaf scattering is approximated as 1/2 of the leaf reflectance 
! plus 1/2 of the leaf transmission 
      do i=1,nclayers
        if(rho(i) + tau(i) > 1d0) then
          if(single_multi_flag == 0)                                    &
            write(*,*) ' rho + tau > 1 in layer ',i
          tau(i) = 0.99d0 - rho(i)
        end if
        psi(i) = (rho(i) + tau(i))  !*5d-1
      end do

      if(iw /= istart) then
        do i=1,nclayers
          rho(i) = 5d-1*(rho(i) + rho_old(i))
          tau(i) = 5d-1*(tau(i) + tau_old(i))
          psi(i) = 5d-1*(psi(i) + psi_old(i))
          alp(i) = 5d-1*(alp(i) + alp_old(i))
          epsc(i) = 5d-1*(epsc(i) + eps_old(i))
          clump(i) = 5d-1*(clump(i) + clump_old(i))
        end do
        psi(nclayers+1) = 5d-1*(psi(nclayers+1) + psi_old(nclayers+1))
      end if

      do i=1,nclayers
        if(dzveg(i) > eps) then
          rho(i) = anint(rho(i)*1d20)*1d-20
          tau(i) = anint(tau(i)*1d20)*1d-20
          psi(i) = anint(psi(i)*1d20)*1d-20
          alp(i) = anint(alp(i)*1d20)*1d-20
          epsc(i) = anint(epsc(i)*1d20)*1d-20
          clump(i) = anint(clump(i)*1d20)*1d-20
        else
          rho(i) = 0d0
          tau(i) = 0d0
          if(albedo > eps) then
            psi(i) = albedo
          else
            if(psi_old(i) > eps) then
              psi(i) = psi_old(i)
            else
              psi(i) = nsoilp(nnodes,3)
            end if
          end if
          psi(i) = anint(psi(i)*1d20)*1d-20
          alp(i) = 1d0 - psi(i)
          alp(i) = anint(alp(i)*1d20)*1d-20
          epsc(i) = anint(emis*1d20)*1d-20
          clump(i) = 0d0
        end if
      end do
      psi(nclayers+1) = anint(psi(nclayers+1)*1d20)*1d-20

      end subroutine canopy_prop

! ******************************************************************************
      subroutine canopy_moist(ncsnow,rain,snow,tot_ep,lai,rl,ra,x,ful,  &
                              ef,def,chf,sigfhi,rh_prof,tai,pheat,      &
                              dpheat,stcl,stcs)

! Calculate canopy evaporation, precipitation interception

! calls the following subroutines:
!     sp_humid

! uses the function: dense,spheats,vap_press

      implicit none

      integer(kind=4),intent(in):: ncsnow
      real(kind=8),intent(in):: lai(nclayers),x(nclayers),ful(nclayers)
      real(kind=8),intent(in):: rh_prof(nclayers),tai(nclayers)
      real(kind=8),intent(out):: rain,snow,tot_ep
      real(kind=8),intent(out):: rl(nclayers),ra(nclayers)
      real(kind=8),intent(out):: ef(nclayers),def(nclayers)
      real(kind=8),intent(out):: chf(nclayers),sigfhi(nclayers)
      real(kind=8),intent(out):: pheat(nclayers),dpheat(nclayers)
      real(kind=8),intent(out):: stcl(nclayers),stcs(nclayers)

! saved variables
      real(kind=8):: stcli(nclayers),stcsi(nclayers)
      save:: stcli,stcsi

! local variables
      integer(kind=4):: i,f1i,d1i
      real(kind=8):: pdens1,f1a,f1,f2,f3,d1,qaf,f2a,mixra
      real(kind=8):: dqdta,vpressa,wetbulba,mixrf,dqdtf,vpressf
      real(kind=8):: cf,rpp,rpf,rptr,pdens,sheatcw,sheatcs,ff,fp
      real(kind=8):: taf,min_wat,max_wat,sheatc1,kths,c1,c2,t0(7)
      real(kind=8):: dense,rhoaf,etr,t1,ep,spheats,vap_press
      real(kind=8):: meltc(nclayers),interc(nclayers),max_wet(nclayers)
      real(kind=8):: max_wetl(nclayers),max_wets(nclayers),drip(nclayers)
      real(kind=8):: interc1(nclayers),drip1(nclayers),precip(nclayers)
      real(kind=8):: precip1(nclayers)

! zero-out variables
      f1i = 0
      d1i = 0
      pdens1 = 0d0
      f1a = 0d0
      f1 = 0d0
      f2 = 0d0
      f3 = 0d0
      c1 = 0d0
      c2 = 0d0
      t1 = 0d0
      d1 = 0d0
      qaf = 0d0
      f2a = 0d0
      mixra = 0d0
      dqdta = 0d0
      vpressa = 0d0
      wetbulba = 0d0
      mixrf = 0d0
      dqdtf = 0d0
      vpressf = 0d0
      ep = 0d0
      cf = 0d0
      rpp = 0d0
      rpf = 0d0
      rptr = 0d0
      pdens = 0d0
      sheatcw = 0d0
      sheatcs = 0d0
      ff = 0d0
      fp = 0d0
      etr = 0d0
      taf = 0d0
      min_wat = 0d0
      max_wat = 0d0
      sheatc1 = 0d0
      kths = 0d0
      rhoaf = 0d0
      tot_ep = 0d0

      do i=1,7
        t0(i) = 0d0
      end do

      do i=1,nclayers
        rl(i) = 0d0
        ra(i) = 0d0
        ef(i) = 0d0
        def(i) = 0d0
        chf(i) = 0d0
        sigfhi(i) = 0d0
        pheat(i) = 0d0
        dpheat(i) = 0d0

        precip(i) = 0d0
        precip1(i) = 0d0
        interc(i) = 0d0
        interc1(i) = 0d0
        drip(i) = 0d0
        drip1(i) = 0d0
        meltc(i) = 0d0
        max_wet(i) = 0d0
        max_wetl(i) = 0d0
        max_wets(i) = 0d0

        stcs(i) = 0d0
        stcl(i) = 0d0
        stcs(i) = storcs(i)
        stcl(i) = storcl(i)
        stcsi(i) = 0d0
        stcli(i) = 0d0
!        if(iterations == 0) then
          stcsi(i) = stcs(i)
          stcli(i) = stcl(i)
!        end if
      end do

! precip interception, maximum storage
      d1 = storcs(1) + storcs(2) + storcs(3)
      if((hsaccum > eps.or.hi > eps).or.d1 > eps) f1i = 1

      c1 = 0d0
      c2 = dmet1(iw,5)*1d-2
      call sp_humid(f1i,met(iw,ip_ap),dmet1(iw,4),c2,c1,mixra,t0(1),    &
                    vpressa,wetbulba,t0(2),t0(3),t0(4),t0(5),t0(6),     &
                    t0(7))

! rain and snow density and specific heat
      d1i = 1
      d1 = 0d0
      pdens = dense(wetbulba,d1,d1i)                                     !kg/m^3, water
      sheatcw = spheats(wetbulba,d1i)                                    !J/kg*K

      if(aint(met(iw,ip_pt)) == 3.or.aint(met(iw,ip_pt2)) == 3) then
        d1i = 4
        pdens1 = dense(wetbulba,dmet1(iw,3),d1i)                           !kg/m^3, snow
        d1i = 2
        sheatcs =spheats(wetbulba,d1i)                                     !J/kg*K
      else
        pdens1 = sdensw
        sheatcs = 2050d0
      end if

! storage of precip on foliage
! References:
!     Aston, A.R. (1979), J. Hydrology, 42, p.383-396
!     Hoyningen-Huene, J.v. (1983), Schriftenrieihe des Deutschen Verbandes fur
!         Wasserwirtshaft U. Kulterbau found at http://public.arcegmo.de/pdf-dateien/

      do i=1,ncsnow
        sigfhi(i) = anint(sigfh*laif(i)*1d20)*1d-20
!     max_wet = 0.935d0 + 0.498d0*lai(i) - 0.00575d0*lai(i)*lai(i)   !H-H
        max_wet(i) = veg_prp(vegh_type,7)*(lai(i)                       &
                                                + veg_prp(vegh_type,8))  !R & S (mm)

        max_wetl(i) = max_wet(i)*1d-3                                    !m
        max_wets(i) = (max_wet(i)*1d-3)*(pdens1/pdens)                   !m
        max_wetl(i) = anint(max_wetl(i)*1d20)*1d-20
        max_wets(i) = anint(max_wets(i)*1d20)*1d-20

        if(i == 1) then
          if(aint(met(iw,ip_pt)) == 2.or.aint(met(iw,ip_pt)) == 4) then !rain
            precip(i) = dmet1(iw,6)*timstep                                      !m/timestep(partial hour)
          else if(aint(met(iw,ip_pt)) == 3) then                        !snow
            precip1(i) = dmet1(iw,6)*timstep                                     !m/timestep(partial hour)
          else if(aint(met(iw,ip_pt2)) == 3)then                        !snow
            precip1(i) = dmet1(iw,7)*timstep                                     !m/timestep(partial hour)
          end if
        else
          precip(i) = dmax1(0d0,precip(i-1) - stcl(i-1))
          precip1(i) = dmax1(0d0,precip1(i-1) - stcs(i-1))
        end if
        precip(i) = anint(precip(i)*1d20)*1d-20
        precip1(i) = anint(precip1(i)*1d20)*1d-20

        t1 = 5d-1*(lai(i) + veg_prp(vegh_type,8))
        t1 = anint(t1*1d20)*1d-20

        if(t1 > 5d1) t1 = 5d1
        if(precip(i) > 0d0) then                                        !rain
          interc(i) = precip(i)*(1d0 - dexp(-t1))*sigfhi(i)              !m
          interc(i) = dmax1(0d0,dmin1(interc(i),max_wetl(i)))
        else if(precip1(i) > 0d0) then                                  !snow
          interc1(i) = precip1(i)*(1d0 - dexp(-t1))*sigfhi(i)            !m
          interc1(i) = dmax1(0d0,dmin1(interc1(i),max_wets(i)))
        end if

        interc(i) = anint(interc(i)*1d20)*1d-20
        interc1(i) = anint(interc1(i)*1d20)*1d-20

        stcs(i) = stcs(i) + interc1(i)
        if (stcs(i) > max_wets(i)) then
          drip1(i) = stcs(i) - max_wets(i)                               !m
          stcs(i) = max_wets(i)                                          !m
        end if
        drip1(i) = anint(drip1(i)*1d20)*1d-20
        stcs(i) = anint(stcs(i)*1d20)*1d-20
        if(stcs(i) < 1d-4) stcs(i) = 0d0

        if(stcs(i) > 0d0) then
!          kths = 2.3d-2 + (7.75d-5*pdens1                               & 
!                    + 1.105d-6*pdens1*pdens1)*(2.29d0 - 2.3d-2)          !W/m*K
          if(pdens1 < 1.56d-1) then                                      !W/m*K, th. cond. (Sturm et al)
            kths = 2.3d-2 + 2.34d-1*pdens1
            kths = dmin1(dmax1(2.3d-2,kths),1d0)
          else
            kths = 1.38d-1 - 1.01d0*pdens1 + 3.233d0*pdens1*pdens1
            kths = dmin1(dmax1(1.38d-1,kths),1d0)
          end if
          meltc(i) = (kths/stcs(i))*tai(i)*timstep*3.6d3/(pdens1*lhfus)
          meltc(i) = dmax1(0d0,dmin1(stcs(i)*pdens/pdens1,meltc(i)))
          meltc(i) = anint(meltc(i)*1d20)*1d-20
          stcs(i) = stcs(i) - meltc(i)
        end if
        stcs(i) = anint(stcs(i)*1d20)*1d-20
        if(stcs(i) < 1d-4) stcs(i) = 0d0

        stcl(i) = stcl(i) + interc(i) + meltc(i)
        if(stcl(i) > max_wetl(i)) then
          drip(i) = stcl(i) - max_wetl(i)
          stcl(i) = max_wetl(i)
        end if
        drip = anint(drip(i)*1d20)*1d-20
        stcl(i) = anint(stcl(i)*1d20)*1d-20
        if(stcl(i) < 1d-6) stcl(i) = 0d0
      end do

! roots
      trmhm = 0d0                                                        !maximum transpiration rate (m/s)
      do i=1,nnodes
        if(stt(i) > Tref) then
          frh(i) = rk(vegh_type,i)*(1d0 - (soil_moist(i) - nsoilp(i,16))&
                                        /(nsoilp(i,17) - nsoilp(i,16)))  !root water fraction
          frh(i) = dmax1(0d0,dmin1(1d0,anint(frh(i)*1d20)*1d-20))
        else
          frh(i) = 0d0
        end if

        trmhm = trmhm + frh(i)
        f2a = f2a + rk(vegh_type,i)*soil_moist(i)
        min_wat = min_wat + nsoilp(i,16)
        max_wat = max_wat + nsoilp(i,17)
      end do

      ff = 1d0
!      if(stt(nnodes) <= Tref.or.hm > eps) ff = 0d0
      trmhm = 1.5d-7*sigfh*ff*trmhm                                      !m/s
      trmhm = anint(trmhm*1d20)*1d-20

! calculate the stomatal resistance for each layer
! REF: http://www.ecmwf.int/research/ifsdocs/CY25r1/PHYSICS/
      f1a = (4d-3*dmet1(iw,1) + 5d-3)/(8.1d-1*(4d-3*dmet1(iw,1) + 1d0))
      f1 = 1d0/dmin1(1d0,f1a)                                            !unitless

      f2 = dmax1(0d0,dmin1(1d0,(f2a - min_wat)/(max_wat - min_wat)))
      if(f2 > eps) f2 = 1d0/f2   !1d-3

      do i=1,ncsnow
        c2 = rh_prof(i)*1d-2
        call sp_humid(f1i,met(iw,ip_ap),tai(i)+Tref,c2,c1,mixra,t0(1),  &
                      vpressa,wetbulba,t0(2),t0(3),t0(4),t0(5),t0(6),   &
                      t0(7))  

        c2 = 1d0
        call sp_humid(f1i,met(iw,ip_ap),x(i)+Tref,c2,c1,mixrf,dqdtf,    &  !t0(1),    &
                      vpressf,t0(2),t0(3),t0(4),t0(5),t0(6),t0(7),t1)

        t1 = 3d-4*(vpressf - vpressa)
        if(vpressf /= vpressa.and.dabs(t1) < 5d1) then
          f3 = dexp(t1)                                                  !unitless
        else
          f3 = 1d0
        end if

        if(lai(i) > eps) then
          rl(i) = dmax1(0d0,(veg_prp(vegh_type,1)/lai(i))*f1*f2*f3)      !s/m
        else
          rl(i) = 0d0
        end if
        rl(i) = anint(rl(i)*1d20)*1d-20

        c1 = vap_press(nnodes+1,met(iw,ip_rh)*1d-2,met(iw,ip_ap))
!        rhoaf = (c1/Rv + (met(iw,ip_ap)*1d2 - c1)/Rd)*5d-1*             &
!                              (1d0/(tai(i) + Tref) + 1d0/(x(i) + Tref))  !foliage air density (kg/m^3)
!        rhoaf = (3.48d-3*met(iw,ip_ap)*1d2*5d-1)                        &                         
!                             *(1d0/(tai(i) + Tref) + 1d0/(x(i) + Tref))  !foliage air density(kg/m^3)
        taf = (1d0 - 0.65d0*sigfh)*(tai(i) + Tref)                       &
                                            + 0.65d0*sigfh*(x(i) + Tref)
        rhoaf =  (c1/Rv + (met(iw,ip_ap)*1d2 - c1)/Rd)/taf               !foliage air density (kg/m^3)
        rhoaf = anint(rhoaf*1d20)*1d-20

        if(ful(i) > eps) then
          cf = (1d0 + 0.3d0/ful(i))  !*1d-2                                 !bulk transfer coefficient (unitless)
          ra(i) = 1d0/(cf*ful(i))                                        !atmospheric resistance to water vapor diffusion (s/m)
        else
          cf = 0d0
          ra(i) = 0d0
        end if
        cf = anint(sigfh*cf*1d20)*1d-20
        ra(i) = anint(ra(i)*1d20)*1d-20

!        if(rl(i) > eps.and.ra(i)+rl(i) > eps)                           &
        if(ra(i)+rl(i) > eps)                           &
          rpp = ra(i)/(ra(i) + rl(i))                                    !vegetation wetness
        rpp = dmin1(dmax1(0d0,anint(rpp*1d20)*1d-20),1d0)

        d1 = 1d0 - 6.5d-1*sigfhi(i)*(1d0 - rpp)
        qaf = ((1d0 - 6.5d-1*sigfhi(i))*mixra                           & 
                                   + 6.5d-1*sigfhi(i)*rpp*mixrf)/d1
        d1i = 1
        ep = lai(i)*cf*(rhoaf/dense(wetbulba,d1,d1i))*ful(i)*           &
                                                      (qaf - rpp*mixrf)  !potential evaporation (m/s)

        if(rpp*mixrf >= qaf.or.trmhm <= eps) then
          if(f2 > min_wat.and.ep <= eps) ep = 0d0
        end if
        ep = anint(ep*1d20)*1d-20

        if(rl(i)+ra(i) > eps) then
          if(stcs(i) > 0d0) then
            ff = (stcs(i)/max_wets(i))**(2d0/3d0)
          else if(stcl(i) > 0d0) then
            ff = (stcl(i)/max_wetl(i))**(2d0/3d0)     
          end if
          rpf = 1d0 - (rl(i)/(rl(i) + ra(i)))*(1d0 - ff)
          rptr = (ra(i)/(rl(i) + ra(i)))*(1d0 - ff)
        end if

        rptr = dmin1(dmax1(0d0,rptr),1d0)
        rptr = anint(rptr*1d20)*1d-20
        rpf = dmin1(dmax1(0d0,rpf),1d0)
        rpf = anint(rpf*1d20)*1d-20

        ef(i) = rpf*ep                                                   !leaf potential evaporation (m/s)
        ef(i) = anint(ef(i)*1d20)*1d-20
        etr = dmin1(trmhm,rptr*ep)                                       !leaf transpiration (m/s)
        etr = anint(etr*1d20)*1d-20

        if(stcs(i) > eps) then
          stcs(i) = stcs(i) - (ef(i) - etr)*(timstep*3.6d3)              !m
          if(stcs(i) < 1d-4) stcs(i) = 0d0
          if (stcs(i) > max_wets(i)) then
            drip1(i) = drip1(i) + (stcs(i) - max_wets(i))                !m
            stcs(i) = max_wets(i)                                        !m
          end if
        else if(stcl(i) > eps) then
          stcl(i) = stcl(i) - (ef(i) - etr)*(timstep*3.6d3)              !m
          if(stcl(i) < 1d-6) stcl(i) = 0d0
          if (stcl(i) > max_wetl(i)) then
            drip(i) = drip(i) + (stcl(i) - max_wetl(i))                  !m
            stcl(i) = max_wetl(i)                                        !m
          end if
        end if

        drip(i) = dmax1(0d0,anint(drip(i)*1d20)*1d-20)
        stcl(i) = dmax1(0d0,anint(stcl(i)*1d20)*1d-20)
        drip1(i) = dmax1(0d0,anint(drip1(i)*1d20)*1d-20)
        stcs(i) = dmax1(0d0,anint(stcs(i)*1d20)*1d-20)

        ef(i) = rhoaf*lai(i)*cf*ful(i)*(qaf - rpp*mixrf)                 !kg/m^2*s
        def(i) = rhoaf*lai(i)*cf*ful(i)*rpp*dqdtf                       &
                                      *(0.65d0*sigfhi(i)/d1 - 1d0)

        ef(i) = anint(ef(i)*1d20)*1d-20
        def(i) = anint(def(i)*1d20)*1d-20

        tot_ep = tot_ep + ef(i)*rpf

!        taf = (1d0 - 0.65d0*sigfh)*tai(i) + 0.65d0*sigfh*x(i)

        if(dabs(lai(i)) <= eps) then
          chf(i) = 0d0
        else
          chf(i) = dmax1(0d0,1.1d0*lai(i)*rhoaf*spheats(taf,3)          &
                                                            *cf*ful(i))  !foliage (W/m^2*K)
        end if
        chf(i) = anint(chf(i)*1d20)*1d-20

        pheat(i) = 1d-2*sigfhi(i)*(sheatcw*pdens*interc(i)              &
                                          + sheatcs*pdens1*interc1(i))  &
                        *(wetbulba - (x(i) + Tref))
        pheat(i) = anint(pheat(i)*1d20)*1d-20
        dpheat(i) = -1d-2*sigfhi(i)*(sheatcw*pdens*interc(i)            & 
                                           + sheatcs*pdens1*interc1(i))
        dpheat(i) = anint(dpheat(i)*1d20)*1d-20
      end do

      rain = 0d0
      rain = dmax1(0d0,(precip(ncsnow) + drip(ncsnow)))                  !m  !sigfhi(ncsnow)*
      rain = anint(rain*1d20)*1d-20

      snow = 0d0
      snow = dmax1(0d0,(precip1(ncsnow) + drip1(ncsnow)))                !m   !sigfhi(ncsnow)*
      snow = anint(snow*1d20)*1d-20

      tot_ep = anint(tot_ep*1d15)*1d-15

      end subroutine canopy_moist

! ******************************************************************************
      subroutine roughness(lamda,beta,h,htemp,d,z0)

!       @(#)roughness.f   2.1   12/8/95   2.1

!     This routine calculates the roughness and zero-displacement 
!     height for a vegetation canopy.

!     Reference:
!     Raupach,1994, Simplified expressions for vegetation
!         roughness length and zero-plane displacement as
!         functions of canopy height and area index.
!         Boundary layer meteorology, vol. 71, pp. 211-216
 
!  Input:
!     lamda - canopy leaf area index     m^2 /m^2 (canopy area index = 2*leaf area index)
!     beta - canopy stem area index      m^2/m^2
!     h - height of the canopy top       m
!  Output:
!     d - zero displacement height of the canopy   m
!     z0 - roughness height of the canopy          m

! no subroutines called

      implicit none

      real(kind=8),intent(in):: lamda,beta,h,htemp
      real(kind=8),intent(out):: d,z0

! local variables
      real(kind=8):: usuh,psih,temp1,temp2,h2

      real(kind=8),parameter:: usuh_max = 0.3d0                          !maximum u*/Uh (unitless)
      real(kind=8),parameter:: cr = 0.3d0                                !drag coeff. of isolated roughness element on substrate 
                                                                         !(unitless)
      real(kind=8),parameter:: cs = 0.003d0                              !drag coeff. of sustrate at height h (unitless)
      real(kind=8),parameter:: cw = 2.0d0                                !constant, roughness sublayer depth (unitless)
      real(kind=8),parameter:: cd1 = 7.5d0                               !constant (unitless)

! zero-out parameters
      usuh = 0d0
      psih = 0d0
      temp1 = 0d0
      temp2 = 0d0
      h2 = 0d0
      d = 0d0
      z0 = 0d0

      psih = dlog(cw) - 1d0 + 1d0/cw                                     !roughness sublayer influence function
      temp1 = dsqrt(cs + cr*(lamda + beta)*5d-1)                         !unitless
      usuh = dmin1(temp1,usuh_max)                                       !u*/Uh (unitless)
      temp2 = dsqrt(cd1*(lamda + beta))
      h2 = h - htemp                                                     !canopy height - (low foliage and/or snow + ice) (m)
      
      if(temp2 > 5d1) then
        d = h2
      else
        d = h2*(1d0 - (1d0 - dexp(-temp2))/temp2)
      end if
      d = anint(d*1d20)*1d-20

      if(dabs(-vK/usuh - psih) > 5d1) then
        z0 = h2
      else
        z0 = h2*(1d0 - d/h2)*dexp(-vK/usuh - psih)
      end if
      z0 = anint(z0*1d20)*1d-20

      end subroutine roughness

! ******************************************************************************
      subroutine windprofile(num_layers,zh,glb,sigfh,dz, wind_prof,     &
                            rh_prof,tai)

!       @(#)windprofile.f 2.1   12/8/95   2.1

!  Purpose:
!       Calculate normalized wind profiles for a general n-layer canopy

!  Input:
!       zh - canopy height (meters)
!       zo - ground roughness length 
!       num_layers - number of layers in the canopy
!       dz(num_layers) - thickness for each layer of the canopy
!              from top to bottom layer
!       dz(num_layers+1) - thickness of "ground layer"
!  Output:
!       wind_prof(num_layers+1) - normalized wind profile factors

! ------------------------------------------------------------------------------
!    Reference: 
!
!    Yamazaki, T, et al. 1992. "A Head-Balance Model with a Canopy of One
!       or Two Layers and its Application to Field Experiments," Journal
!       of Applied Meterology, Volume 31, Pages 86-103.
!
! ------------------------------------------------------------------------------

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers
      real(kind=8),intent(in):: zh,sigfh,glb,dz(num_layers)
      real(kind=8),intent(out):: wind_prof(num_layers)
      real(kind=8),intent(out):: rh_prof(num_layers)
      real(kind=8),intent(out):: tai(num_layers)

! local variables 
      integer(kind=4):: i
      real(kind=8):: t1,cstar,xc,cstar_prime,f,zhalf_layer,zdcist
      real(kind=8):: zpercent,zb


! zero-variables
      t1 = 0d0
      cstar = 0d0                                                        !canopy density (unitless)
      xc = 0d0
      cstar_prime = 0d0
      f = 0d0
      zhalf_layer = 0d0
      zdcist = 0d0
      zpercent = 0d0

!      iheightn = 0d0
      do i=1,num_layers
        wind_prof(i) = 0d0
        rh_prof(i) = 0d0
        tai(i) = 0d0
      end do

      cstar = 1d0
      xc = dlog10(cstar)

! This is from equation (32) in the reference. The exp is the NASA version, The
! 10^ is the original as it appears in the paper cited
!     cstar_prime = 10**((((xc+0.26)+dsqrt((xc+0.26)**2+0.16)))/2.0
!    $  -0.3)
      t1 = xc + 0.26d0
      if(dabs((t1 + dsqrt(t1*t1 + 0.16d0))*5d-1 - 0.3d0) < 5d1) then
        cstar_prime = dexp((t1 + dsqrt(t1*t1 + 0.16d0))*5d-1 - 0.3d0)
      end if

! weighting function for cstar, this is from equation (31) in the reference
      f = 0.37d0 + (0.494d0*(xc + 0.8d0))/dsqrt((xc + 0.8d0)*           &
                                                (xc - 0.5d0) + 1.1d0)
                                                
! loop though all layers and calculate wind speed factors
      zdcist = 0d0
      zb = 1d0/zh
      do i=1,num_layers
! determine the midpoint of the layer in percentage of the entire canopy height
        zhalf_layer = dz(i)*5d-1
        if(i == 1) then
          zdcist = zdcist
        else
          zdcist = zdcist + dz(i-1)                                      !calculate distance from the top
        end if
        zpercent = (zh - (zhalf_layer + zdcist))*zb

! This is from equation (30) in the reference
!        wind_prof(i)=f*dexp((-cstar_prime/(2.0*vK*vK))*(1.0 - zpercent))&
!                        + (1.0 - f)*dlog(zh*zpercent/zo)/dlog(zh/zo)
        if(dabs(1.1d0*(zpercent - 1d0))*(1d0 - sigfh*1d-1) < 5d1)       &
          wind_prof(i) = dexp(1.1d0*(zpercent - 1d0))*(1d0 - sigfh*1d-1)
        if(wind_prof(i) <= eps) wind_prof(i) = 0d0

        if(dabs(5.5d0*(zpercent - 1d0))*(1d0 - sigfh*1d-1) < 5d1)       &
          rh_prof(i) = dexp(5.5d0*(zpercent - 1d0))*(1d0 - sigfh*1d-1)
        if(rh_prof(i) <= eps) rh_prof(i) = 0d0

        if(glb <= eps) then
          if(dabs(-5d-2*(zpercent - 1d0))*(1d0 - sigfh*1d-1) < 5d1)     &
            tai(i) = dexp(-5d-2*(zpercent - 1d0))*(1d0 - sigfh*1d-1)
        else
          if(dabs(5d-2*(zpercent - 1d0))*(1d0 - sigfh*1d-1) < 5d1)      &
            tai(i) = dexp(5d-2*(zpercent - 1d0))*(1d0 - sigfh*1d-1)
        end if
        if(tai(i) <= eps) tai(i) = 0d0

        wind_prof(i) = anint(wind_prof(i)*1d20)*1d-20
        rh_prof(i) = anint(rh_prof(i)*1d20)*1d-20
        tai(i) = anint(tai(i)*1d20)*1d-20
      end do

      end subroutine windprofile

! ******************************************************************************
      subroutine scalc_2a(num_layers,max_thetar, n_thetar,n_thetak,&
                         foliage_type,lai,sai,clump,Sij,Wir)
 

! @(#)scalc_2a.f  2.1   12/8/95

! calls the following subroutines:
!     check (attached to this subroutine) - currently not used
!     kernal (attached to this subroutine)
!     leaf_slope (attached to this subroutine)
!     mean_proj (attached to this subroutine)
!     gap_prob (attached to this subroutine)
!     fmatrix (attached to this subroutine)
!     smatrix (attached to this subroutine)

      implicit none

      integer(kind=4),intent(in):: num_layers,max_thetar,n_thetar
      integer(kind=4),intent(in):: n_thetak,foliage_type(num_layers)
      real(kind=8),intent(in):: lai(num_layers),sai(num_layers)
      real(kind=8),intent(in):: clump(num_layers)
      real(kind=8),intent(out):: Sij(num_layers+2,num_layers+2)
      real(kind=8),intent(out)::Wir(num_layers+2,max_thetar)

! local variables
      integer(kind=4),parameter:: max_thetak = 30
      integer(kind=4):: i,j,r
      real(kind=8):: fk(max_thetar,max_thetak),f(num_layers,max_thetak)
      real(kind=8):: sinr(max_thetar),cosr(max_thetar)
      real(kind=8):: thetar(max_thetar),thetak(max_thetak)
      real(kind=8):: g_bar(num_layers,max_thetar)
      real(kind=8):: po(num_layers,max_thetar)
      real(kind=8):: Fijr(num_layers+2,num_layers+2,max_thetar)


! zero-out variables
      do i=1,num_layers+2
        do j=1,max_thetar
          Wir(i,j) = 0d0
        end do
      end do

! *** Compute kernal used for the mean layer projection of leaf distribution angle
       call kernal(max_thetar,max_thetak,n_thetar,n_thetak,thetar,      &
                   thetak,sinr,cosr,fk)

! *** Compute leaf slope distributions 
      call leaf_slope(num_layers,max_thetak,n_thetak,foliage_type,      &
                      thetak,f)

! *** Compute mean projection
      call mean_proj(num_layers,max_thetak,max_thetar,n_thetar,         &
                     n_thetak,fk,f,g_bar)

! *** Compute gap probabilities
      call gap_prob(num_layers,max_thetar,n_thetar,lai,sai,clump,g_bar, &
                    cosr,po)

! *** Compute Fijr
      call fmatrix(num_layers,max_thetar,n_thetar,po,Fijr)

! *** Extract Canopy View Factors, Wir, from first row of Fijr
      do j= 1,num_layers+2
        do r = 1,n_thetar+1
          Wir(j,r) = anint(Fijr(1,j,r)*1d20)*1d-20
        end do
      end do

! *** Compute Sij
      call smatrix(num_layers,max_thetar,n_thetar,sinr,cosr,Fijr,Sij)   

      end subroutine scalc_2a

! ******************************************************************************
      subroutine kernal(max_thetar,max_thetak,nthetar,nthetak,thetar,   &
                        thetak,sinr,cosr,fk)

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: max_thetar,max_thetak,nthetar,nthetak
      real(kind=8),intent(out):: thetar(max_thetar),sinr(max_thetar)
      real(kind=8),intent(out):: cosr(max_thetar),thetak(max_thetak)
      real(kind=8),intent(out):: fk(max_thetar,max_thetak)
      
! local variables 
      integer(kind=4):: i,j
      real(kind=8):: phi_prime,fi,fj,f1,f2

! zero-out variables
      phi_prime = 0d0
      fi = 0d0
      fj = 0d0
      f1 = 0d0
      f2 = 0d0

      do i=1,max_thetar
        sinr(i) = 0d0
        cosr(i) = 0d0
        thetar(i) = 0d0
      end do

      do i=1,max_thetak
        thetak(i) = 0d0
      end do

      do i=1,max_thetar
        do j=1,max_thetak
          fk(i,j) = 0d0
        end do
      end do

      f2 = pi/180d0

! Sin and cos of the angles of integration over elevation evenly distributes
      f1 = 1d0/float(nthetar)
      do i = 2,nthetar
        fi = i-1
        thetar(i) = (fi*f1)*pi*5d-1                                      !thetar not 90.0
        sinr(i) = dsin(thetar(i))
        cosr(i) = dcos(thetar(i))
      end do

! Handle end point (0 and 90 degrees)
      thetar(1) = 0.1d0*f2
      sinr(1) = dsin(thetar(1))
      cosr(1) = dcos(thetar(1))
      thetar(nthetar+1) = 89.9d0*f2                                      !last thetar set to 89.9
      sinr(nthetar+1) = dsin(thetar(nthetar+1))
      cosr(nthetar+1) = dcos(thetar(nthetar+1))

! same procedure for leaf distribtion angles
      f1 = 1d0/float(nthetak)
      do j=2,nthetak
        fj = j-1
        thetak(j) = (fj*f1) * pi*5d-1                                    !thetak not 90.0
      end do
      thetak(1) = 0.1d0*f2
      thetak(nthetak+1) = 89.9d0*f2

! Calculate the kernal used for the mean layer projection of leaf distribution angle 
! thetak over directions thetar
      f1 = 2d0/pi
      do j=1,nthetak+1
        do i=1,nthetar+1
          if(thetak(j) <= (pi*5d-1 - thetar(i))) then 
            fk(i,j) = cosr(i) * dcos(thetak(j))
          else
            phi_prime = dacos(-(1d0/dtan(thetak(j)))*                   &
                                                 (1d0/dtan(thetar(i))))
            fk(i,j) = f1*cosr(i)*dcos(thetak(j))*                       &
                                (phi_prime - pi*5d-1 - dtan(phi_prime))
          end if
          fk(i,j) = fk(i,j)*f1                                           !sf added to match swoe rep
          fk(i,j) = anint(fk(i,j)*1d20)*1d-20
        end do
      end do

      end subroutine kernal

! ******************************************************************************
      subroutine leaf_slope(num_layers,max_thetak,nthetak,foliage_type, &
                            thetak,f)

! *** Calculate probability density function (PDF) for
! *** different element slope  distributions.

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers,max_thetak,nthetak
      integer(kind=4),intent(in):: foliage_type(num_layers)
      real(kind=8),intent(in):: thetak(max_thetak)
      real(kind=8),intent(out):: f(num_layers,max_thetak)

! local variables
      integer(kind=4):: i,j
      real(kind=8):: f1


! zero-out variables
      do i=1,num_layers
        do j=1,max_thetak
          f(i,j) = 0d0
        end do
      end do

      f1 = 0d0
      f1 = 2d0/pi

      do i=1,num_layers
        if((foliage_type(i) == 1).or.(foliage_type(i) == 2)) then
          do j = 1,nthetak+1
            f(i,j) = f1*(1d0 + dcos(2d0*thetak(j)))                      !planophile, erectophile
          end do
        elseif((foliage_type(i) == 3).or.(foliage_type(i) == 4)) then
          do j = 1,nthetak+1
            f(i,j) = f1*(1d0 - dcos(4d0*thetak(j)))                      !plagiophile, extremophile
          end do
        elseif(foliage_type(i) == 5) then
          do j=1,nthetak+1
            f(i,j) = f1                                                  !uniform
          end do
        else
          do j=1,nthetak+1
            f(i,j) = dsin(thetak(j))                                     !spherical
          end do
        end if

        do j=1,nthetak+1
          f(i,j) = anint(f(i,j)*1d20)*1d-20
        end do
      end do

      end subroutine leaf_slope

! ******************************************************************************
      subroutine mean_proj(num_layers,max_thetak,max_thetar,nthetar,    &
                           nthetak,fk,f,g_bar)

!  Purpose:
!       Compute mean projection of all (leaves) in a layer in
!       directions thetar
 
!  Input:
!       fk   - kernal function      (dimensionless)
!       f    - leaf slope function  (dimensionless)
!  Output:
!       g_bar    (dimensionless)
  
! ------------------------------------------------------------------------------
!    References: 
  
!    J.A. Smith, R.E. Oliver and J.K. Berry
!       A comparison of two photographic techniques for 
!       estimating foliage angle distribution  (See Appendix 1)
!       Aust. J. Bot, 1977, vol 25, pp 545-553
!
!    Press, et al.  Numerical Recipes in Fortran (Cambridge University
!      Press)  Using Formula 4.1.11 on Page 127, Second Edition, to
!      (Extended Trapezoidal Rule)
!      compute numerical integral on an closed interval.. i
!      have computed the leaf slope and kernal functions to not include
!      theta = 90 
!
!    int (x1 to xn) { f(x)dx} = h [ 1/2 f1 + f2 + f3 + ... fn-1 + 1/2 fn]

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers,max_thetak,max_thetar
      integer(kind=4),intent(in):: nthetar,nthetak
      real(kind=8),intent(in):: fk(max_thetar,max_thetak)
      real(kind=8),intent(in):: f(num_layers,max_thetak)
      real(kind=8),intent(out):: g_bar(num_layers,max_thetar)

! local variables
      integer(kind=4),parameter:: maxtk = 30
      integer(kind=4):: i,j,k
      real(kind=8):: h,sumf,integrand(num_layers,max_thetar,maxtk)


! zero-out variables
      h = 0d0
      sumf = 0d0

      do i=1,num_layers
        do j=1,max_thetar
          g_bar(i,j) = 0d0
          do k=1,maxtk
            integrand(i,j,k) = 0d0
          end do
        end do
      end do

      h = (1d0/nthetak)

      do i=1,num_layers
        do j=1,nthetar+1
          do k=1,nthetak+1
            integrand(i,j,k) = f(i,k)*fk(j,k)
          end do
        end do
      end do

      do i=1,num_layers
        do j=1,nthetar+1
          sumf = 0d0
          do k=2,nthetak
            sumf = sumf + integrand(i,j,k)
          end do
          g_bar(i,j) = h*(0.5d0*(integrand(i,j,1) +                     &
                                    integrand(i,j,nthetak+1)) + sumf)
          g_bar(i,j) = anint(g_bar(i,j)*1d20)*1d-20
        end do
      end do

      end subroutine mean_proj

! ******************************************************************************
      subroutine gap_prob(num_layers,max_thetar,nthetar,lai,sai,clump,  &
                          g_bar,cosr,po)

! *** Calculate gap probability along ray direction.
!     clump is markov clumping factor between 0 and 1
!       clump = 1 means no clumping (should be default)
!       clump = 0 means total clumping, i.e. po = 1

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers,max_thetar,nthetar
      real(kind=8),intent(in):: lai(num_layers),sai(num_layers)
      real(kind=8),intent(in):: clump(num_layers),cosr(max_thetar)
      real(kind=8),intent(in):: g_bar(num_layers,max_thetar)
      real(kind=8),intent(out):: po(num_layers,max_thetar)

! local variables
      integer(kind=4):: i,j
      real(kind=8):: t1,clmp(num_layers)


! zero-out variables
      t1 = 0d0

      do i=1,num_layers
        clmp(i) = 0d0
        do j=1,max_thetar
          po(i,j) = 0d0
        end do
      end do

      do i=1,num_layers
        clmp(i) = clump(i)
        if(clump(i) < eps) clmp(i) = 1d-2
        do j=1,nthetar+1
          t1 = (lai(i) + sai(i))*clmp(i)*g_bar(i,j)/cosr(j)
          if(dabs(t1) > 5d1) then
            t1 = 1.1d0*(lai(i) + sai(i))*clmp(i)*g_bar(i,j)/cosr(j-1)
          end if
          t1 = dmin1(1d0,t1)
!          if(cosr(j) >= 1d-6.and.dabs(t1) < 5d1) then
            po(i,j)= dexp(-t1)
!          else
!            po(i,j) = 0d0      
!          end if
          po(i,j) = anint(po(i,j)*1d20)*1d-20
        end do
      end do

      end subroutine gap_prob

! ******************************************************************************
      subroutine fmatrix(num_layers,max_thetar,nthetar,po,Fijr)                 

! *** Calculate view factor coefficients Fijr for sink /source
! *** layer for direction thetar

! *** Sink and source index definitions:
! ***   1   - Sky               
! ***   2   - Leaf layer 1  
! ***   3   - Leaf layer 2
! ***   .
! ***   n   - Leaf layer n
! ***   n+1 - Ground            

! *** Remember! LAI indices start at one, not two

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers,max_thetar,nthetar
      real(kind=8),intent(in):: po(num_layers,max_thetar)
      real(kind=8),intent(out):: Fijr(num_layers+2,num_layers+2,max_thetar)

! local variables
      integer(kind=4):: i,j,r,k,m
      real(kind=8):: prod,f


! zero-out variables
      prod = 0d0
      f = 0d0

      do i=1,num_layers+2
        do j=1,num_layers+2
          do r = 1,max_thetar
            Fijr(i,j,r) = 0d0
          end do
        end do
      end do

      do r=1,nthetar+1
! CONTRIBUTION COEFFICIENTS TO SKY.
        Fijr(1,1,r) = 0d0                                                !From Sky to Sky
            
        do k=2,num_layers+1                                              !From leaf layers to sky
          prod = 1d0
          if(k > 2) then
            do m=2,k-1
              prod = prod*po(m-1,r)
            end do
          end if
          Fijr(1,k,r) = (1d0 - po(k-1,r))*prod
        end do
 
        prod = 1d0                                                       !From Ground to Sky
        do k=2,num_layers+1
          prod = prod*po(k-1,r)
        end do
        Fijr(1,num_layers+2,r) = prod

! CONTRIBUTION COEFFICIENTS TO LEAF LAYERS.
        do i=1,num_layers+2
          do j=2,num_layers+1
            if(j < i) then
              prod = 1d0
              if(j+1 < i) then
                do k=j+1,i-1
                  prod = prod*po(k-1,r)
                end do
              end if
              if(i /= num_layers+2) then
!                f = (1d0 - po(j-1,r))*(1d0 - po(i-1,r))*prod
                f = dsqrt(po(j-1,r))*(1d0 - po(i-1,r))*prod
              else
!                f = (1d0 - po(j-1,r))*prod                              !GROUND
                f = dsqrt(po(j-1,r))*prod
              end if
            else if(i == j) then
!             f = 2.0*po(j-1,r)
              f = 2d0*(1d0 - dsqrt(po(j-1,r)))
            else                                                         !i_sink  >  i_source
              prod = 1d0
              if(j > i+1) then
                do k=i+1,j-1
                  prod = prod*po(k-1,r)
                end do
              end if

              if(i /= 1) then    !i.e. if not equal to sky
!                f = (1d0 - po(j-1,r))*(1d0 - po(i-1,r))*prod
                f = dsqrt(po(j-1,r))*(1d0 - po(i-1,r))*prod
              else
!                f = (1d0 - po(j-1,r))*prod                              !SKY
                f = dsqrt(po(j-1,r))*prod
              end if

            end if
            Fijr(j,i,r) = f
          end do
        end do

! CONTRIBUTION COEFFICIENTS TO GROUND.
        prod = 1d0                                                       !From sky to ground.
        do k=2,num_layers+1
          prod = prod*po(k-1,r)
        end do

        Fijr(num_layers+2,1,r) = prod

        do k=2,num_layers+1                                              !From leaf layer to ground.
          prod = 1d0
          if(k+1 <= num_layers+1) then
            do m=k+1,num_layers+1
              prod = prod*po(m-1,r)
            end do
          end if
          Fijr(num_layers+2,k,r) = (1d0 - po(k-1,r))*prod
        end do

        Fijr(num_layers+2,num_layers+2,r) = 0d0                          !From ground to ground

        do i=1,num_layers+2
          do j=1,num_layers+2
            Fijr(i,j,r) = anint(Fijr(i,j,r)*1d20)*1d-20
          end do
        end do
      end do

      end subroutine fmatrix

! ******************************************************************************
      subroutine smatrix(num_layers,max_thetar,nthetar,sinr,cosr,Fijr,  &
                         Sij)

!    References: 
  
!    Press, et al.  Numerical Recipes in Fortran (Cambridge University
!      Press)  Using Formula 4.1.11 on Page 127, Second Edition, to
!      (Extended Trapezoidal Rule)
!      compute numerical integral on an closed interval.. i
!      have computed the leaf slope and kernal functions to not include
!      theta = 90  (set to 89.9)

!    int (x1 to xn) { f(x)dx} = h [ 1/2 f1 + f2 + f3 + ... fn-1 + 1/2 fn]

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers,max_thetar,nthetar
      real(kind=8),intent(in):: sinr(max_thetar),cosr(max_thetar)
      real(kind=8),intent(in):: Fijr(num_layers+2,num_layers+2,max_thetar)
      real(kind=8),intent(out):: Sij(num_layers+2,num_layers+2)

! local variables
      integer(kind=4):: i,j,r
      real(kind=8):: h,sumf
      real(kind=8):: integrand(num_layers+2,num_layers+2,max_thetar)


! zero-out variables
      h = 0d0
      sumf = 0d0

      do i=1,num_layers+2
        do j=1,num_layers+2
          Sij(i,j) = 0d0
        end do
      end do

      do i=1,num_layers+2
        do j=1,num_layers+2
          do r=1,max_thetar
            integrand(i,j,r) = 0d0
          end do
        end do
      end do

! Integrate the view factors over all view angles to get the exchange coefficients
! between the source/sinks.
! NOTE: integrand = Cijk in the text
      h = (1d0/nthetar)
      do i=1,num_layers+2
        do j=1,num_layers+2
          do r=1,nthetar+1
            integrand(i,j,r) = sinr(r)*cosr(r)*Fijr(i,j,r)               !angles of integration
          end do
        end do
      end do

      do i=1,num_layers+2
        do j=1,num_layers+2
          sumf = 0d0
          do r=2,nthetar
            sumf = sumf + integrand(i,j,r)                               !integration-no end points
          end do
          Sij(i,j) = pi*h*(0.5d0*(integrand(i,j,1) +                    &
                       integrand(i,j,nthetar+1)) + sumf)                 !final integration
          Sij(i,j) = anint(Sij(i,j)*1d20)*1d-20
        end do
      end do

      end subroutine smatrix

! ******************************************************************************
      subroutine radiosity(num_layers,Sij,psi,absorb,aground)

! radiosity.f 2.1   12/8/95

! Ref: Smith, J.A. and S.M. Goltz (1994)"A Thermal Exitance and Energy Balance Model 
!      for Forest Canopies", IEEE Trans. on Geosci. and Rem. Sens., V.32(5), pp.1060-1066.

! calls the following subroutines:
!     matinv (attached to this subroutine)
!     ludcmp,lubksb (attached to this subroutine)
!     bfunc,rfunc,qfunc (attached to this subroutine)

      implicit none

      integer(kind=4),intent(in):: num_layers
      real(kind=8),intent(in):: Sij(num_layers+2,num_layers+2)
      real(kind=8),intent(in):: psi(num_layers+1)
      real(kind=8),intent(out):: absorb(num_layers+1),aground

! local variables
      integer(kind=4):: i,j,d1
      real(kind=8):: albedoc,total,total_absorption
      real(kind=8):: b(num_layers+1),abs_coeff(num_layers+1)
      real(kind=8):: emission(num_layers+1),scatt(num_layers+1)
      real(kind=8):: matrix(num_layers+1,num_layers+1)
      real(kind=8):: matinvs(num_layers+1,num_layers+1)
      real(kind=8),parameter:: etotal = 1d0


! zero-out parameters
      d1 = 0
      albedoc = 0d0
      total = 0d0
      total_absorption = 0d0
      aground = 0d0

      do i=1,num_layers+1
        abs_coeff(i) = 0d0
        emission(i) = 0d0
        do j=1,num_layers+1
          matrix(i,j) = 0d0
        end do
      end do

      do i=1,num_layers+1
        absorb(i) = 0d0
      end do

      d1 = num_layers + 1

      do i=1,d1
        emission(i) = psi(i)*etotal*Sij(i+1,1)                           !eqn. (8) of ref.
      end do

      do j=1,d1
        do i=1,d1
          matrix(i,j) = -psi(i)*Sij(i+1,j+1)
          if(i == j) matrix(i,j) = 1d0 + matrix(i,j)
        end do
      end do

      call matinv(d1,matrix,matinvs)                                     !solve for radiosity

      do i = 1,d1
        b(i) = 0d0
        do j=1,d1
          b(i) = b(i) + matinvs(i,j)*emission(j)                         !radiosity, eqn. (7) of ref.
        end do
      end do

! *** Calculate "abs_coeff(i)", the absorption coefficient.
      do i=1,num_layers
        abs_coeff(i) = 1d0 - 2d0*psi(i)                                  !Leaf layers
      end do
      abs_coeff(d1) = 1d0 -  psi(num_layers+1)                               !Ground layer

! *** Compute multiple scattering contribution
      do i = 1,d1
        scatt(i) = 0d0
        do j=1,d1
          scatt(i) = scatt(i) + b(j)*Sij(i+1,j+1)
        end do
      end do

! The total solar flux reaching the ground can be obtained from 
!   Total solar flux at the ground = total solar flux at the canopy times
!   Sij(num_layers+2,SKY) + scatt(num_layers+1) 
      aground = dmax1(0d0,Sij(num_layers+2,1) + scatt(d1))
      aground = anint(aground*1d20)*1d-20

! Energy balance check
      do i=1,d1
        absorb(i) = abs_coeff(i)*(etotal*Sij(i+1,1) + scatt(i))          !eqn. (9) of ref.
        absorb(i) = anint(absorb(i)*1d20)*1d-20

        total_absorption  = total_absorption  + absorb(i)
        albedoc = albedoc + b(i)*Sij(1,i+1)
      end do

      total = total_absorption + albedoc

      end subroutine radiosity

! ******************************************************************************
      subroutine feval(num_layers,icall,Sij,x,epsc,alp,bta,btg,bx,a,qx, &
                       rx,ef,def,chf,sigfhi,fx,dfx,tai,pheat,dpheat)

!       @(#)feval.f   2.1   12/8/95   2.1

!     Sign Conventions:  Rn = LE + H
!       i.e. In = Out
!            B, and A (radiative) are positive  for flux into system (canopy)
!            LE, H are positive for flux out of system (canopy)
!
!        F = Rn - LE - H

! calls the following subroutines:
!     bfunc (attached to this subroutine)
!     rfunc (attached to this subroutine) - currently not used
!     qfunc (attached to this subroutine)

      implicit none

      integer(kind=4),intent(in):: num_layers,icall
      real(kind=8),intent(in):: x(num_layers),epsc(num_layers)
      real(kind=8),intent(in):: alp(num_layers),ef(num_layers)
      real(kind=8),intent(in):: def(num_layers),tai(num_layers)
      real(kind=8),intent(in):: a(num_layers),chf(num_layers)
      real(kind=8),intent(in):: sigfhi(num_layers),pheat(num_layers)
      real(kind=8),intent(in):: dpheat(num_layers),bta,btg
      real(kind=8),intent(in):: Sij(num_layers+2,num_layers+2)
      real(kind=8),intent(out):: bx(num_layers),qx(num_layers)
      real(kind=8),intent(out):: rx(num_layers),fx(num_layers)
      real(kind=8),intent(out):: dfx(num_layers,num_layers)

! local variables
      integer(kind=4):: i,j,il
      real(kind=8):: lw_abs,lwabs,lwemit,lecan,lhevap
      real(kind=8):: two_fac(num_layers),drx(num_layers),dbx(num_layers)
      real(kind=8):: dqx(num_layers)


! zero-out variables
      lw_abs = 0d0
      lwabs = 0d0
      lwemit = 0d0
      lecan = 0d0
      lhevap = 0d0

      do i=1,num_layers
        bx(i) = 0d0                                                      !long wave
        qx(i) = 0d0                                                      !sensible
        rx(i) = 0d0                                                      !latent
        fx(i) = 0d0                                                      !net
        two_fac(i) = 0d0
        do j=1,num_layers
          dfx(i,j) = 0d0                                                 !d(net)/dT
        end do
      end do

      do i=1,num_layers
! scale H and LE by lai 
        call bfunc(epsc(i),x(i),bx(i),dbx(i))

        call qfunc(x(i),chf(i),sigfhi(i),qx(i),dqx(i),tai(i))

        lhevap = 2500775.6d0 - 2369.729d0*tai(i)                         !J/kg = (m/s)^2
        rx(i) = anint(ef(i)*lhevap*1d20)*1d-20
        drx(i) = anint(def(i)*lhevap*1d20)*1d-20
      end do

      do il=1,num_layers
        lw_abs = 0d0                                                    
        lecan = lecan + rx(il)
        do j=1,num_layers
          lw_abs =  lw_abs +  bx(j)*Sij(il+1,j+1)
          two_fac(il) = two_fac(il) + Sij(j+1,il+1)
        end do
        two_fac(il) = two_fac(il) + Sij(1,il+1) + Sij(num_layers+2,il+1) !this should be exactly 2.0, however is sum of Sij

        fx(il) = alp(il)*sigma*(bta*Sij(il+1,1) + lw_abs                &
                  + btg*Sij(il+1,num_layers+2))                         &
                  - two_fac(il)*sigma*bx(il)                            &
                  + a(il) + qx(il) + rx(il) + pheat(il)
        fx(il) = anint(fx(il)*1d20)*1d-20

        if(icall == 1) then
          lwabs = lwabs + alp(il)*sigma*(bta*Sij(il+1,1) + lw_abs       & 
                  + btg*Sij(il+1,num_layers+2))
          lwemit = lwemit + two_fac(il)*sigma*bx(il)
        end if
      end do

      do il=1,num_layers
        do j=1,num_layers
          if(j == il) then
            dfx(il,j) = alp(il)*sigma*dbx(il)*Sij(il+1,j+1)             &
                        - two_fac(il)*sigma*dbx(il) + dqx(il) + drx(il) &
                        + dpheat(il)
          else
            dfx(il,j) = alp(il)*sigma*dbx(j)*Sij(il+1,j+1)
          end if
          dfx(il,j) = anint(dfx(il,j)*1d20)*1d-20
        end do
      end do

      end subroutine feval

! ******************************************************************************
      subroutine bfunc(epsi,xi,bxi,dbxi)

! Calculate the emitted flux normalized by Stephen-Boltzman constant
! Note no sigma

! no subroutines called
     
      implicit none

      real(kind=8),intent(in):: epsi,xi
      real(kind=8),intent(out):: bxi,dbxi

! local variables
      real(kind=8):: tc


      bxi = 0d0
      dbxi = 0d0
      tc = 0d0

      tc = anint((xi + Tref)*1d20)*1d-20
      if(tc > eps) then
        bxi =  epsi*tc*tc*tc*tc                                          !emissivity*T^4
        bxi = anint(bxi*1d20)*1d-20

        dbxi = 4d0*epsi*tc*tc*tc                                         !derivative wrt to T
        dbxi = anint(dbxi*1d20)*1d-20
      end if

      end subroutine bfunc

! ******************************************************************************
      subroutine qfunc(xi,chfi,sigfh,qxi,dqxi,taii)
      
!  Purpose:
!       Calculate Sensible Heat Flux
 
!  Input:
!       xi  - canopy temperature     (degress C)
!       tac - air temperature        (degress C)
!       chfi - sensible heat coeff.  (W/m^2*K)
!       sigfh - foliage density      (unitless)
!       ful - wind speed profile     (m sec-1)
!  Output:
!       qxi - Sensible Heat Flux (watts m-2)
!       dqxi- deriviative wrt xi  (Used for Newton Method)
 
! ------------------------------------------------------------------------------
!    Reference: 
!    Tibbals, E.C., Carr, E.K. Gates, D.M. and Kreith, F. (1964),
!     Radiation and convection in conifers. Am. J. Bot. 51:529-538
! ------------------------------------------------------------------------------

!         We are in the process of totally modifying and
!         generalizing this routine.  It will be replaced in its
!         entirety.

! no subroutines called

      implicit none

      real(kind=8),intent(in):: xi,chfi,sigfh,taii
      real(kind=8),intent(out):: qxi,dqxi

! local variables
      real(kind=8):: hci,tafi

! zero-out variables
      hci = 0d0
      tafi = 0d0
      qxi = 0d0
      dqxi = 0d0
      
      tafi = (1d0 - 0.65d0*sigfh)*taii + 0.65d0*sigfh*xi
      qxi = chfi*(tafi - xi)  !*1d1
      qxi = anint(qxi*1d20)*1d-20

      dqxi = chfi*(0.65d0*sigfh - 1d0)  !*1d1
      dqxi = anint(dqxi*1d20)*1d-20

      end subroutine qfunc

! ******************************************************************************
      subroutine output(num_layers,Sij,bx,btg,bta,bgr)
!      subroutine output(ncsnow,max_thetar,Sij,bx,btg,bta,Wir,bgr)

! output.f    1.5   12/7/95

!  Purpose:
!       Print output results at each time step

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers !,max_thetar
      real(kind=8),intent(in):: Sij(num_layers+2,num_layers+2)
      real(kind=8),intent(in):: bx(num_layers),btg
      real(kind=8),intent(in):: bta  !,Wir(num_layers+2,max_thetar)
      real(kind=8),intent(out):: bgr

! local variables
      integer(kind=4):: i,r,fi
!!      real(kind=8):: ert(max_thetar)
      real(kind=8):: bx_tosky(num_layers+1),bx_tognd(num_layers)
      real(kind=8):: bx_fromsky(num_layers+1),bx_fromgnd(num_layers)

! zero-out variables
      r = 0
      fi = 0
      bgr = 0d0

      do i=1,num_layers
        bx_tognd(i) = 0d0
        bx_fromgnd(i) = 0d0
      end do

      do i=1,num_layers+1
        bx_tosky(i) = 0d0
        bx_fromsky(i) = 0d0
      end do

! note in the Sij matrix nclayers+2 is the ground, layer 1 is the sky 
! and layer 2 the top canopy. For all other arrays nclayers+1 is the grnd

      do i=1,num_layers
        bx_tosky(i) = sigma*bx(i)*Sij(1,i+1)                             !to sky from each layer
        bx_tognd(i) = sigma*bx(i)*Sij(num_layers+2,i+1)                  !to grnd from each layer
        bgr = bgr + bx_tognd(i)
        bx_fromsky(i) = sigma*bta*Sij(i+1,1)                             !from sky to each layer
        bx_fromgnd(i) = sigma*btg*Sij(i+1,num_layers+2)                  !from grnd to each layer
      end do

      bx_tosky(num_layers+1) = sigma*btg*Sij(1,num_layers+2)             !from grnd to sky
      bx_fromsky(num_layers+1) = sigma*bta*Sij(num_layers+2,1)           !from sky to grnd
      bgr = bgr + bx_fromsky(num_layers+1)                               !total reaching the ground
      bgr = anint(bgr*1d20)*1d-20

! Compute an effective canopy+grnd temperature for each view direction         
!     do r=1,nthetar+1
!       ert(r) = 0d0
!       sum = 0d0
!       do i=2,num_layers+1
!          sum = sum + Wir(i,r)*sigma*bx(i-1)
!        end do
!       ert(r) = sum + Wir(num_layers+2,r)*sigma*btg
!       if(ert(i) > eps) ert(r) = ((ert(r)/sigma)**0.25d0) - Tref
!     end do

      end subroutine output

! ******************************************************************************
      subroutine solve(num_layers,y,a,dx)

! solve.f 2.1   12/8/95

!     nclayers is physical dimension of assumed
!     largest jacobian--matrix to be 
!     inverted.  It affects working array dimensions.  
!     The Jacobian matrix, a, is actually of dimension
!      n by m and enters through the calling argument

!  Purpose:
!       Solve general non-linear equations y(x) = 0 
!     Expand y in Taylor series (J = Jacobian)
!     and find that x which makes y(x) = 0, i.e.
!
!     First    y(x) = y(xo) + J(x-xo) + ...
!     Then     0 = y(xo) + J(x-xo) 
!
!     or  using generalized inverse (i.e. can take inverse of
!     non square matrix and even better for square matrices)
!
!     x-xo = (J'J)^-1 J [-y]

!  Input:
!       y - equations to be solved 
!       ncs - number of equations to be solved (dimension of y)
!       a - Jacobian matrix, i.e. partial of y wrt unknowns x
!           dimension of a is ncs x m
!       m - number of unknowns (dimension of x or dx)
!
!     If ncs > m   this formulation gives least squares solution
!        ncs = m   use of  generalized inverse in this case okay too
!        ncs < m   may be underdetermined
        
!  Output:
!       dx = x-xo     of dimension m

!  subroutines Called:
!     matinv.f    inverse of a matrix
! ------------------------------------------------------------------------------

! calls the following subroutines:
!     matinv (attached to this subroutine)

      implicit none

      integer(kind=4),intent(in):: num_layers
      real(kind=8),intent(in):: y(num_layers),a(num_layers,num_layers)
      real(kind=8),intent(out):: dx(num_layers)

! local variables
      integer(kind=4):: i,j,k
      real(kind=8):: ata(num_layers,num_layers),aty(num_layers)
      real(kind=8):: ata_inv(num_layers,num_layers)


      do i=1,num_layers
        do j=1,num_layers
          ata(i,j) = 0d0
          do k=1,num_layers
            ata(i,j) = ata(i,j) + a(k,i)*a(k,j)
          end do
        end do
      end do

      call matinv(num_layers,ata,ata_inv)

      do i=1,num_layers
        aty(i) = 0d0
        do j=1,num_layers
          aty(i) = aty(i) + a(j,i)*(-y(j))
        end do
      end do

      do i=1,num_layers
        dx(i) = 0d0
        do j=1,num_layers
          dx(i) = dx(i) + ata_inv(i,j)*aty(j)
          dx(i) = anint(dx(i)*1d15)*1d-15
        end do
      end do

      end subroutine solve

! ******************************************************************************
      subroutine matinv(num_layers,a,y)

!       @(#)matinv.f  2.1   12/8/95   2.1
      
!   Purpose:
!      Find inverse of matrix A with logical dimension n by n
!      and physcial dimension np x np
!      Original Matrix, a,  is destroyed

!   Input:
!       a - original matrix    (destroyed)
!       num_layers - physical dimensionof a

!   Output:
!       y - inverse matrix

!   Reference:
!      From Numerical Recipes, (Fortran Version), page 38:      
!      By William H. Pres, et al. (1986)

! calls the following subroutines:
!     ludcmp (attached to this subroutine)
!     lubksb (attached to this subroutine)

      implicit none

      integer(kind=4),intent(in):: num_layers
      real(kind=8),intent(inout):: a(num_layers,num_layers)
      real(kind=8),intent(out):: y(num_layers,num_layers)
      
! local variables 
      integer(kind=4):: i,j,indx(num_layers)
      real(kind=8):: d


      do i=1,num_layers
        do j=1,num_layers
          y(i,j) = 0d0
        end do
        y(i,i) = 1d0
      end do
 
      call ludcmp(num_layers,indx,a,d)

      do j=1,num_layers
        call lubksb(num_layers,indx,a,y(1,j))
      end do

      end subroutine matinv

! ******************************************************************************
      subroutine ludcmp(num_layers,indx,a,d)

!     LU Decomposition Algorithm for Matrix Inversion
!     From Numerical Recipes (Fortran Version), page 35
!     By William H. Pres, et. al. (1986)

!     Matrix A has physical dimension  num_layers x num_layers
!     Original Matrix A is destroyed

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers
      integer(kind=4),intent(inout):: indx(num_layers)
      real(kind=8),intent(inout):: a(num_layers,num_layers)
      real(kind=8),intent(out):: d

! local variables
      integer(kind=4):: i,j,k,imax
      real(kind=8):: aamax,sumf,dum,vv(num_layers)
      real(kind=8),parameter:: TINY1 = 1d-10

! zero-out variables
      imax = 0
      aamax = 0d0
      sumf = 0d0
      dum = 0d0

      do i=1,num_layers
        vv(i) = 0d0
      end do

      d = 1d0
      do i=1,num_layers
        aamax = 0d0
        do j=1,num_layers
if(dabs(a(i,j)) <= TINY1)then
write(*,*)'ludcmb',iw,i,j,a(i,j)
stop
end if
          if(dabs(a(i,j)) > aamax) aamax = dabs(a(i,j))
        end do
        if (aamax <= TINY1) aamax = TINY1                                !singular matrix
        vv(i) = 1d0/aamax
      end do

      do j=1,num_layers
        if (j > 1) then
          do i=1,j-1
            sumf = a(i,j)
            if (i > 1) then
              do k=1,i-1
                sumf = sumf - a(i,k)*a(k,j)
              end do
              a(i,j) = sumf
            end if
          end do
        end if

        aamax = 0d0
        do i=j,num_layers
          sumf = a(i,j)
          if (j > 1) then
            do k=1,j-1
              sumf = sumf - a(i,k)*a(k,j)
            end do
            a(i,j) = sumf
          end if
          dum = vv(i)*dabs(sumf)

          if((dum >= aamax).or.(dabs(dum-aamax) < 1d-10)) then
            imax = i
            aamax = dum
          end if
        end do

        if (j /= imax) then
          do k=1,num_layers
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          end do
          d = -d
          vv(imax) = vv(j)

        end if

        if(imax /= 0) then
          indx(j) = imax
        else
          indx(j) = 1
        end if
!if(dabs(a(j,j)) <= TINY1)write(*,*)'ludcmb2',j,a(j,j)
        if(dabs(a(j,j)) <= TINY1) a(j,j) = TINY1
        if(j /= num_layers) then
          dum = 1d0/a(j,j)
          do i=j+1,num_layers
            a(i,j) = a(i,j)*dum
          end do
        end if
      end do

      d = anint(d*1d20)*1d-20
!if(dabs(a(num_layers,num_layers)) <= TINY1)write(*,*)'ludcmb3',num_layers,a(num_layers,num_layers)
      if(a(num_layers,num_layers) <= TINY1)                             &
                                       a(num_layers,num_layers) = TINY1

      end subroutine ludcmp

! ******************************************************************************
      subroutine lubksb(num_layers,indx,a,b)

!     LU Back Substituion Algorithm for Matrix Inversion

!     From Numerical Recipes (Fortran Version), page 36
!     By William H. Pres, et. al. (1986)

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: num_layers,indx(num_layers)
      real(kind=8),intent(in):: a(num_layers,num_layers)
      real(kind=8),intent(inout):: b(num_layers)

! local variables
      integer(kind=4):: i,j,ii,ll
      real(kind=8):: sumf
      real(kind=8),parameter:: TINY1 = 1d-10


! zero-out variables
      ii = 0
      sumf = 0d0

      do i=1,num_layers
        ll = indx(i)
        sumf = b(ll)
        b(ll) = b(i)
        if(ii /= 0) then
          do j=ii,i-1
            sumf = sumf - a(i,j)*b(j)
          end do
        else if(sumf /= 0d0) then
          ii = i
        end if
        b(i) = sumf
      end do

      do i=num_layers,1,-1
        sumf = b(i)
        if(i < num_layers) then
          do j=i+1,num_layers
            sumf = sumf - a(i,j)*b(j)
          end do
        end if
        if(dabs(a(i,i)) > TINY1) then
          b(i) = sumf/a(i,i)
        else
!          write(*,*)'lubskb',i,a(i,i)
          b(i) = sumf*1d10
        end if
      end do

      end subroutine lubksb

! ******************************************************************************
      subroutine veg_proph(biome_source,new_vt,veg_type)

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: biome_source,new_vt,veg_type

! local variables
      integer(kind=4):: i,j,vind,jj,io,vid,hlines,fid,ntypes,err
      real(kind=8):: zup,zdw,cinfo(7,5),sumf,t1,t2,t3,t4,ar,br
      real(kind=8):: srmax,srmin,cmax,cmin,rl,laimax,laimin,ddmax,sai
      real(kind=8):: folamin,folamax,heigmin,heigmax,emissmin,emissmax
      real(kind=8):: defveg_prp(18,17)
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
!   12 = emissivity min (%)
!   13 = emissivity max (%)
!   14 = absorbtion min (%)
!   15 = absorbtion max (%)
!   16 = height min (cm)
!   17 = height max (cm)

! cinfo(row,column) 
! rows are: needle leaf evergreen, broad leaf evergreen, needle leaf decid.
!           broad leaf decid., mixed
! columns are: canopy height, layer thicknesses ratio to total thickness (2-4), 
!              lai factor[sum=1](5-7)
! heights are based on Dorman & Sellers (1989)
! NOTE: heights are now read in from FASST_defaultvegprop.inp file or something similar
      data cinfo/25.5d0, 0.439215686d0, 0.462745098d0, 0.098039216d0, 0.4d0, 0.5d0, 0.1d0,    & !needle leaf evergreen (3)
                 36.0d0, 0.452777778d0, 0.477777778d0, 0.069444444d0, 0.4d0, 0.5d0, 0.1d0,    & !broadleaf evergreen (6)
                 21.0d0, 0.428571429d0, 0.452380952d0, 0.119047619d0, 0.4d0, 0.5d0, 0.1d0,    & !needle leaf deciduous (4)
                 31.5d0, 0.434920635d0,  0.46984127d0, 0.095238095d0, 0.4d0, 0.5d0, 0.1d0,    & !broadleaf deciduous (5)
                 30.0d0,        0.45d0, 0.466666667d0, 0.083333333d0, 0.4d0, 0.5d0, 0.1d0/      !mixed (18)

! zero-out parameters
      vind = 0
      hlines = 0
      fid = 0
      ntypes = 0
      zup = 0d0
      zdw = 0d0
      sumf = 0d0
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

      if(infer_test == 0) then
        zh = 0d0

        do i=1,nclayers
          dzveg(i) = 0d0
          laif(i) = 0d0
        end do
      end if

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
          defveg_prp(veg_type,16) = heigmin*1d-2
          defveg_prp(veg_type,17) = heigmax*1d-2
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
                               laimin,ddmax,sai,ar,br,emissmin,emissmax,&
                               folamin,folamax,heigmin,heigmax

          if(new_vt == vid) then
            jj = 1
            veg_prp(new_vt,1) = srmin
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
            newveg_prp(new_vt,16) = heigmin*1d-2
            newveg_prp(new_vt,17) = heigmax*1d-2
          end if
        end do
     
        rewind(fid)
      end if

      if(biome_source == 0) then
        do i=1,17
          veg_prp(veg_type,i) = defveg_prp(veg_type,i)
          veg_prp(veg_type,i) = anint(veg_prp(veg_type,i)*1d20)*1d-20
        end do
      else
        do i=1,17
          if(dabs(newveg_prp(new_vt,i)-spflag) <= eps) then 
            veg_prp(veg_type,i) = defveg_prp(veg_type,i)
          else
            veg_prp(veg_type,i) = newveg_prp(new_vt,i)
          end if
          veg_prp(veg_type,i) = anint(veg_prp(veg_type,i)*1d20)*1d-20
        end do
      end if

! root fraction in layer i for each veg type
      do i=1,nnodes
        if(i == 1) then
          zdw = elev - nzi(i)
          zup = elev - (nzi(i) + nzi(i+1))*5d-1
        else if(i == nnodes) then
          zdw = elev - (nzi(i) + nzi(i-1))*5d-1
          zup = elev - nzi(i)
        else
          zdw = elev - (nzi(i) + nzi(i-1))*5d-1
          zup = elev - (nzi(i) + nzi(i+1))*5d-1
        end if
        
        t1 = veg_prp(veg_type,9)*zdw
        if(t1 > 5d1) t1 = 5d1
        t2 = veg_prp(veg_type,10)*zdw
        if(t2 > 5d1) t2 = 5d1
        t3 = veg_prp(veg_type,9)*zup
        if(t3 > 5d1) t3 = 5d1
        t4 = veg_prp(veg_type,10)*zup
        if(t4 > 5d1) t4 = 5d1

        rk(veg_type,i) = -0.5d0*(dexp(-t1) + dexp(-t2)                  &
                                               - dexp(-t3) - dexp(-t4))
        rk(veg_type,i) = anint(rk(veg_type,i)*1d20)*1d-20
        if(dabs(rk(veg_type,i)) < eps) rk(veg_type,i) = 0d0
      end do

! canopy only
      if(infer_test == 0) then
        if(vegh_type == 3.or.vegh_type == 4) then                        !needle leaf
          if(vegh_type == 3) then
            vind = 1
          else
            vind = 3
          end if
          zh = veg_prp(vegh_type,17) !cinfo(1,vind)
          do j=1,nclayers
            dzveg(j) = anint(cinfo(j+1,vind)*veg_prp(vegh_type,17)*1d2) &
                                                                  *1d-2
            laif(j) = cinfo(j+4,vind)
          end do
        else if(vegh_type == 5.or.vegh_type == 6) then                   !broad leaf
          if(vegh_type == 5) then
            vind = 4
          else
            vind = 2
          end if
          zh = veg_prp(vegh_type,17)
          do j=1,nclayers
            dzveg(j) = anint(cinfo(j+1,vind)*veg_prp(vegh_type,17)*1d2) &
                                                                  *1d-2
            laif(j) = cinfo(j+4,vind)
          end do
        else if(vegh_type == 18) then                                    !mixed
          vind = 5
          zh = veg_prp(vegh_type,17)
          do j=1,nclayers
            dzveg(j) = anint(cinfo(j+1,vind)*veg_prp(vegh_type,17)*1d2) &
                                                                  *1d-2
            laif(j) = cinfo(j+4,vind)
          end do
        end if

        if((vegh_type >= 3.and.vegh_type <= 6).or.vegh_type == 18) then
          do i=1,nclayers
            if(dabs(izh-spflag) > eps.and.                              &
                                    dabs(idzveg(i)-spflag) <= eps) then
              dzveg(i) = izh*dzveg(i)/zh
            else if(dabs(izh-spflag) > eps.and.                         &
                                     dabs(idzveg(i)-spflag) > eps) then
              dzveg(i) = idzveg(i)
            end if
            dzveg(i) = anint(dzveg(i)*1d5)*1d-5
          end do
          if(dabs(izh-spflag) > eps) zh = izh

! check for layer thickness closure
          do j=1,nclayers
            sumf = sumf + dzveg(j)
          end do

          if(dabs(sumf - zh) > 1d-3) then
            if(single_multi_flag == 0) write(*,'('' canopy + ground '', &
              &''thickness not equal to total height'')')
            stop
          end if

!          if(dzveg(nclayers) > eps) then
!            iheightn = zh - (dzveg(1) + dzveg(2) + dzveg(nclayers)*5d-1)
!          else
!            iheightn = dmin1(iheight,zh - (dzveg(1) + dzveg(2)*5d-1))
!          end if
        end if
      end if   !infer_test == 0

      end subroutine veg_proph

end module module_canopy

      subroutine surfenergy(ii,iter,sn,kave,sthick,pdens,dmet,evapcm,   &
                            rhsurf,pheatf,pheatg,mixrf,mixraf,mixrgrs,  &
                            mixra,lhtf,lheatg,sheatf,sheatg,dqdtf,      &
                            dqdtg,sfac,d1,rpp,cc1,isurfg,isurff,        &
                            disurfg,disurff,disurffg,disurfgf,lh)

      use fasst_global
      use module_lowveg

! calls the following subroutines:
!     sp_humid
!     lowveg_met
!     sflux

! uses the function: dense,spheats,thconds,soilhumid,drag (appended to this subroutine)

! NOTE: W = J/s

      implicit none

      integer(kind=4),intent(in):: ii,iter,sn
      real(kind=8),intent(in):: kave,sthick
      real(kind=8),intent(inout):: pdens,dmet(13)
      real(kind=8),intent(out):: evapcm,rhsurf,pheatf,pheatg,mixrf
      real(kind=8),intent(out):: mixraf,mixrgrs,mixra,lhtf
      real(kind=8),intent(out):: lheatg,sheatf,sheatg,dqdtf,dqdtg,sfac
      real(kind=8),intent(out):: d1,rpp,cc1,isurfg,isurff,disurfg
      real(kind=8),intent(out):: disurff,disurffg,disurfgf,lh

! local variables
      integer(kind=4):: f1i,d1i,d2i
      real(kind=8):: vpress,ldragcoeff,sdragcoeff,spheat,pdens1
      real(kind=8):: wetbulba,radsu,dqdta,rhoaf,rhoag
      real(kind=8):: cf,vpressa,shmax,spheatr,spheats,kd,kv,ka,dh,dv
      real(kind=8):: rhov,rhoda,dense,epf1,smr,net,t1,t2,t3,t4,t5
      real(kind=8):: eps1,thconds,vpsat,sheats,lai,taf,c1,c2,mixrgr
      real(kind=8):: soilhumid,radsd,radld,radlu,lhw,sh,pht1,stempt
      real(kind=8):: drag,newdrag


! zero-out parameters
      f1i = 0
      d1i = 0
      d2i = 0
      c1 = 0d0
      c2 = 0d0
      vpress = 0d0
      ldragcoeff = 0d0
      sdragcoeff = 0d0
      spheat = 0d0
      lh = 0d0
      pdens1 = 0d0
      lai = 0d0
      wetbulba = 0d0
      radsu = 0d0
      radsd = 0d0
      radlu = 0d0
      radld = 0d0
      dqdta = 0d0
      rhoaf = 0d0
      rhoag = 0d0
      eps1 = 0d0
      cf = 0d0
      vpressa = 0d0
      kd = 0d0
      kv = 0d0
      ka = 0d0
      dh = 0d0
      dv  = 0d0
      rhov = 0d0
      rhoda = 0d0
      epf1 = 0d0
      spheatr = 0d0
      sheats = 0d0
      pdens = 0d0
      rhsurf = 0d0
      evapcm = 0d0
      pheatf = 0d0
      pheatg = 0d0
      mixrf = 0d0
      mixraf = 0d0
      mixrgr = 0d0
      mixrgrs = 0d0
      mixra = 0d0
      lhtf = 0d0
      lheatg = 0d0
      sheatf = 0d0
      sheatg = 0d0
      dqdtf = 0d0
      dqdtg = 0d0
      sfac = 0d0
      taf = 0d0
      d1 = 0d0
      rpp = 0d0
      cc1 = 0d0
      isurfg = 0d0
      isurff = 0d0
      disurfg = 0d0
      disurff = 0d0
      disurffg = 0d0
      disurfgf = 0d0
      lhw = 0d0
      sh = 0d0
      pht1 = 0d0
      stempt = 0d0
      t1 = 0d0
      t2 = 0d0
      t3 = 0d0
      t4 = 0d0
      t5 = 0d0
      shmax = 0d0
      smr = 0d0

      stempt = stt(sn)

! CALCULATE THE SENSIBLE AND LATENT HEAT TERMS
! The latent and sensible heat formulations are based on pg.1389-1390 of
! Hughes, P.A.,et al.(1993),"A mathematical model for the prediction of 
! temperature of man-made and natural surfaces",Int. J. of Remote Sensing,
! V.14,No.7,pp.1383-1412.
      if(hm > eps.or.(node_type(nnodes) == 'SN'.or.                     &
                                       node_type(nnodes) == 'WA')) then
        rhsurf = 1d0
        smr = 1d0
      else if(hpond > eps.or.met(iw,ip_prec)+met(iw,ip_prec2) > eps) then
        rhsurf = 1d0
        smr = 1d0
      else if(stt(nnodes) <= Tref.or.ice(nnodes) > eps) then
        rhsurf = 1d0
        smr = 1d0
      else
        rhsurf = soilhumid(nnodes,phead(nnodes),soil_moist(nnodes),     &
                                                            stt(nnodes))
        smr = soil_moist(nnodes)/nsoilp(nnodes,9)
      end if
      if(vegl_type == 8) rhsurf = smr

      if(stempt < 35.9d0) then
        if(single_multi_flag == 0) write(10,'(''freq_id'',i10,          &
          &''toptemp = '',f8.3,'' airtemp = '',f8.3,'' Top temp too '', &
          &''low in surfenergy.f; day'',i4,'' hour'',i4)') freq_id,     &
          stempt,dmet(4),int(met(iw,ip_doy)),int(met(iw,ip_hr))
        stempt = dmet(4)
        error_code = 1
        error_type = 3
      end if

      if(hm > eps) f1i = 1
      c2 = 1d0
      call sp_humid(f1i,dmet(11),stempt,c2,phead(nnodes),mixrgrs,dqdtg, &
                    t1,t2,rhov,rhoda,vpsat,t3,t4,t5)

      rhoag = rhov + rhoda                                               !kg/m^3

      if(sigfl <= eps.or.icase == 4) then
        ftemp = dmet(4)                                                  !air temperature
        if(icase == 4) then
          ftemp = stt(nnodes) + (hfol_tot/hm)                           &
                                         *(stt(nnodes+1) - stt(nnodes))
          ftemp = anint(ftemp*1d20)*1d-20
        end if
        taf = dmet(4)
        uaf = dmet(3) !dmax1(1d-1,dmet(3))
        epf1 = 0d0
        chnf = 0d0
        rpp = 0d0
        sheatf = 0d0
        pheatf = 0d0
        lhtf = 0d0
        d1 = 1d0

        c2 = dmet(5)*1d-2
        call sp_humid(f1i,dmet(11),dmet(4),c2,c1,mixra,dqdta,vpressa,   &
                      wetbulba,rhov,rhoda,vpsat,t1,t2,t3)
        mixraf = mixra
      else
        if(icase == 2) then
          ftemp = stt(nnodes+1)
        else if(icase == 3) then
          ftemp = stt(nnodes+2)
        end if
        epf1 = epf
        taf = (1d0 - 7d-1*sigfl)*dmet(4) + sigfl*(6d-1*ftemp           &
                                                        + 1d-1*stempt)   !foliage temp at the ground surface (K)

        call lowveg_met(ii,iter,int(met(iw,ip_pt)),int(met(iw,ip_pt2)), &
                        mixrgrs,rhsurf,taf,dmet,wetbulba,rhoaf,cf,      &
                        pheatf,rpp,mixra,mixrf,dqdtf)

        eps1 = epf1 + emis - epf1*emis
        if(eps1 > eps) cc1 = epf1*emis*sigfl*sigma/eps1

        d1 = 1d0 - sigfl*(6d-1*(1d0 - rpp) + 1d-1*(1d0 - rhsurf))
        mixraf = (1d0 - 7d-1*sigfl)*mixra + 6d-1*sigfl*rpp*mixrf        &
                                          + 1d-1*sigfl*rhsurf*mixrgrs
        mixraf = mixraf/d1
      end if
      d1 = anint(d1*1d20)*1d-20
      cc1 = anint(cc1*1d20)*1d-20
      mixraf = anint(mixraf*1d20)*1d-20

      meltfl = 'm'
! latent heat of evaporation (lhevap) is from the WES formulation
      lh = 2500775.6d0 - 2369.729d0*(dmet(4) - Tref)                     !lhevap, J/kg = (m/s)^2
!      if(hm > eps.and.(vpsat-6d2 <= eps.and.dmet(4) < Tref)) then
      if(vpsat-6d2 <= eps.and.dmet(4) < Tref) then
        meltfl = 's'
        lh = lhsub
      end if

! solve for the drag coefficient
      if(iheightn <= hm.or.iheightn <= eps) then
        sdragcoeff = 2d-3 + 6d-3*(elev*2d-4)                             !unitless
        if(hsaccum > eps.or.hi > eps) sdragcoeff = 1.35d-3               !Rachel
        ldragcoeff = sdragcoeff
      else
! dv, etc. from www.stanford.edu/group/efmh/FAMbook2dEd (Fundamentals of Atmospheric Modeling by M. Jacobson)
        if(dmet(4) > eps) dv = 2.11d-5*((dmet(4)/Tref)**1.94d0)         &
                                                  *(1013.25d0/dmet(11))  !molecular diffusion coeff. water vapor

        if(sigfl > eps.and.(icase == 2.or.icase == 3)) then
!          chnf = sqrt_chnf*sqrt_chnf
          d1i = 1
          chnf = newdrag(d1i,dmet(4),ftemp,dv)
        end if
!       ldragcoeff = drag(dmet(4),stempt,dv)
        d1i = 0
        ldragcoeff = newdrag(d1i,dmet(4),stempt,dv)
        ldragcoeff = dmax1(0d0,(1d0 - sigfl)*ldragcoeff + sigfl*chnf)

        d1i = 0
        kv = thconds(dmet(4),d1i)                                        !W/m*K, water vapor
        d1i = 3
        kd = thconds(dmet(4),d1i)                                        !W/m*K, dry air
        ka = kd*(1d0 - (1.17d0 - 1.02d0*kv/kd)*dmet(5)*1d-2)             !thermal conductivity moist air
        dh = ka/(rhoag*((1d0 + 0.87d0*mixra)*spheats(dmet(4),d1i)))      !molecular thermal diffusion coeff.

        if(sigfl > eps.and.(icase == 2.or.icase == 3)) then
!          chnf = sqrt_chnf*sqrt_chnf
          d1i = 1
          chnf = newdrag(d1i,dmet(4),ftemp,dh)
        end if
!       sdragcoeff = drag(dmet(4),stempt,dh)
        d1i = 0
        sdragcoeff = newdrag(d1i,dmet(4),stempt,dh)
        sdragcoeff = dmax1(0d0,(1d0 - sigfl)*sdragcoeff + sigfl*chnf)
      end if

! LATENT HEAT LOSS/GAIN, SURFACE EVAPORATION/SUBLIMATION
! lheat > 0 -> qair > qgr => condensation; surface warms, air cools
! lheat < 0 -> qair < qgr => evaporation;  surface cools, air warms

      if(hm <= eps) ldragcoeff = ldragcoeff
      lheatg = ldragcoeff*lh*uaf*rhoag
      lheatg = anint(lheatg*1d10)*1d-10

      evapcm = lheatg*(mixraf - rhsurf*mixrgrs)                         &
                                            /(dense(stempt,0d0,d1i)*lh)  !m/s
      evapcm = aint(evapcm*1d10)*1d-10

      lhtf = lail*cf*lh*uaf*rhoaf                                        !W/m^2
      lhtf = anint(lhtf*1d20)*1d-20

! SENSIBLE HEAT LOSS/GAIN
! sheat > 0 -> Tair > Tgr => surface warms, air cools
! sheat < 0 -> Tair < Tgr => surface cools, air warms
      d1i = 3
      spheat = spheats(taf,d1i)

      shmax = 0d0
      if(hm > eps) shmax = 2d0
      sheatg = dmax1(shmax,sdragcoeff*spheat*uaf*rhoag)                  !W/m^2*K
      sheatg = anint(sheatg*1d20)*1d-20

      if(dabs(lail) <= eps) then
        sheatf = 0d0
      else
        shmax = 0d0 !2d0
        if(hm > eps) shmax = 2d0
        sheatf = dmax1(shmax,1.1d0*lail*cf*spheat*uaf*rhoaf) !*1d1)         !foliage (W/m^2*K)
      end if
      sheatf = anint(sheatf*1d20)*1d-20

! HEAT LOSS DUE TO PRECIPITATION
! based on Jordan, R. CRREL Rep. 91-16, pp.29,15-17
! Note:  All freezing and thawing taken care of in snow_model.f and icethick.f
      if(aint(met(iw,ip_pt)) == 2.or.aint(met(iw,ip_pt)) == 4) then      !rain, freezing rain
         d1i = 1
         pdens1 = dense(wetbulba,0d0,d1i)                                !kg/m^3
         spheatr = spheats(wetbulba,d1i)                                 !J/kg*K
      else if(aint(met(iw,ip_pt)) == 3.or.dmet(7) > 0d0) then            !snow
        d1i = 4
        pdens = dense(wetbulba,dmet(3),d1i)
        sheats = 2.09d3                                                  !J/kg*K
      end if
      ptemp = wetbulba

      if(pdens > eps.and.dmet(7) <= eps) then
        pheatg = sheats*pdens*dmet(6)
      else
        pheatg = spheatr*pdens1*dmet(6) + sheats*pdens*dmet(7)           !W/m^2*K
      end if

      pheatg = (1d0 - sigfl)*pheatg
      pheatg = anint(pheatg*1d20)*1d-20

! calculate solar penetration term
      t2 = 0d0
      if(hm <= eps) then                                                 !no snow or ice
        sfac = 1d0                                                       !unitless (shortwave attenuation)
      else if(hsaccum > eps.or.newsd > eps) then                         !snow
        if(newsd <= eps) then
          if(iw /= 1) then
            t2 = (hsaccum)*sdens(iw-1)*3.795d-3/dsqrt(dsnow)
            if(t2 > 5d1) t2 = 5d1
            sfac = 1d0 - dexp(-0.8d0)*dexp(-t2)                          !unitless 
          else
            t2 = (hsaccum)*sdensw*3.795d-3/dsqrt(1d-3)
            if(t2 > 5d1) t2 = 5d1
            sfac = 1d0 - dexp(-0.8d0)*dexp(-t2)                          !unitless 
          end if
        else
          t2 = (newsd)*sdens(iw-1)*3.795d-3/dsqrt(dsnow)
          if(t2 > 5d1) t2 = 5d1
          sfac = 1d0 - dexp(-0.8d0)*dexp(-t2)                            !unitless
        end if
      else if(hi > eps.and.(hsaccum <= eps.and.newsd <= eps)) then       !ice, no snow
        t2 = 0.8d0 + hi*idens*3.795d-3/dsqrt(1d-3)
        if(t2 > 5d1) t2 = 5d1
        sfac = 1d0 - dexp(-t2)
      else if(node_type(nnodes) == 'WA') then
        if(ice(nnodes) > eps) then
          t2 = 0.8d0 + met(iw,ip_hi)*idens*3.795d-3/dsqrt(1d-3)
          if(t2 > 5d1) t2 = 5d1
          sfac = 1d0 - dexp(-t2)
        else
          sfac = 1d0 - dexp(0.4d0)
        end if
      end if
      sfac = anint(sfac*1d20)*1d-20

! calculate the solar and incoming IR terms
      d1i = 0
      d2i = 1
      t1 = 0d0
      t2 = 1d0
      t3 = 0d0
      call sflux(d1i,d2i,sn,dmet(8),dmet(1),dmet(12),dmet(2),sfac,      &
                 stempt,cc1,taf,pheatg,sheatg,lheatg,mixraf,mixrgrs,    &
                 kave,stt(sn-1),sthick,rhsurf,rpp,ice(nnodes),          &
                 wvc(nnodes),soil_moist(nnodes),dqdtg,dqdtf,d1,radsd,   &
                 radsu,radld,radlu,pht1,sh,lhw,t1,t3,net,disurfg,       &
                 disurfgf)

      isurfg = radsd + radsu + radld + radlu + sh + lhw + pht1
      isurfg = anint(isurfg*1d20)*1d-20

      if(icase == 2.or.icase == 3) then 
        d1i = 1
        t1 = 0d0
        t3 = 0d0
        t4 = 0d0
        call sflux(d1i,d2i,sn,dmet(8),dmet(1),dmet(12),dmet(2),t2,ftemp,&
                   cc1,taf,pheatf,sheatf,lhtf,mixraf,mixrf,kveg,stempt, &
                   hfol_tot,rpp,rhsurf,ice(nnodes),wvc(nnodes),         &
                   soil_moist(nnodes),dqdtf,dqdtg,d1,radsd,radsu,radld, &
                   radlu,pht1,sh,lhw,t1,t3,t4,disurff,disurffg)

        isurff = radsd + radsu + radld + radlu + sh + lhw + pht1
        isurff = anint(isurff*1d18)*1d-18
      else
        isurff = 0d0
        disurff = 0d0
        disurffg = 0d0
      end if

      end subroutine surfenergy
  
! ******************************************************************************	  
! Mehod 1: Compute ustar using Louis (1979). Use this value to compute z0h and z0q 
! using Jacobson (2005). Use these values to find gammah and gammae using Mascart (1995).
! Then calculate the drag coefficients.
  
      real(kind=8) function drag(taf,stempt,dp)

      use fasst_global

      implicit none

      real(kind=8),intent(in):: taf,stempt,dp

! local variables
      real(kind=8):: ht,z0g,Rib,f1,Cdng0,Gammam,ustar,z0p,f2,Cpng0
      real(kind=8):: Gammap,mu,chs,ph,ch,c,ug


      ht = 0d0
      z0g = 0d0
      Rib = 0d0
      f1 = 0d0
      Cdng0 = 0d0
      Gammam = 0d0
      ustar = 0d0
      z0p = 0d0
      f2 = 0d0
      Cpng0 = 0d0
      Gammap = 0d0
      mu = 0d0
      chs = 0d0
      ph = 0d0
      ch = 0d0
      c = 0d0
      ug = 0d0

      ht = iheightn - hm
      ug = uaf

      if(hm <= eps) then
        z0g = rough
      else if(hm > eps) then
        z0g = 7.775d-3                                                   !(0.05 - 1.5mm)
      end if

! calculate the Richardson Number
! Rib < 0 - unstable; Rib = 0 - neutral stablilty; Rib > 0 - stable
      if(dabs(taf-stempt) <= eps.or.ug <= eps) then
        Rib = 0d0
      else
        Rib = 2d0*grav*ht*(taf - stempt)/(ug*ug*(taf + stempt))
      end if

! calculate z0h and z0q; use Louis(1979) to calculate ustar
      f1 = 1d0/z0g
      Cdng0 = vK*vK/(dlog(ht*f1)*dlog(ht*f1))

      if(Rib < 0d0) then                                                 !unstable
        c = 7.4d0*Cdng0*9.4d0*dsqrt(ht*f1)
        Gammam = 1d0 - 9.4d0*Rib/(1d0 + c*dsqrt(dabs(Rib)))       
      else if(dabs(Rib) < eps) then                                      !neutral
        Gammam = 1d0
      else if(Rib > eps) then                                            !stable
        Gammam = 1d0/((1d0 + 4.7d0*Rib)*(1d0 + 4.7d0*Rib))
      end if

      ustar = dsqrt(Cdng0*Gammam)*ug
      if(dabs(ustar) <= eps) then
        z0p = z0g
      else
        z0p = dp/(4d-1*ustar)                                            !FAM2dEd
        if(z0p >= z0g) z0p = 0.9d0*z0g
      end if
      f2 = 1d0/z0p

! use the method of Mascart et al. (1995)
      Cpng0 = vK*vK/(0.74d0*dlog(ht*f1)*dlog(ht*f2))

      if(hm <= eps) then
        mu = dlog(z0g*f2)
        chs = 3.2165d0 + 4.3431d0*mu + 5.36d-1*mu*mu                    &
                                                   - 7.81d-2*mu*mu*mu
        ph = 5.802d-1 - 1.571d-1*mu + 3.27d-2*mu*mu - 2.6d-3*mu*mu*mu
        if(ph > eps) ch = chs*Cpng0*9.4d0*(dlog(ht*f1)/dlog(ht*f2))     &
                                                         *((ht*f2)**ph)
       Gammap = 1d0 - 1d1*Rib/(1d0 + ch*dsqrt(dabs(Rib)))
       if(Gammap <= 0d0) Gammap = (1d0 + 4.7d0*Rib)**(-2d0)
      else
        if(Rib < 0d0) then                                               !unstable
          mu = dlog(z0g*f2)
          chs = 3.2165d0 + 4.3431d0*mu + 5.36d-1*mu*mu                  &
                                                     - 7.81d-2*mu*mu*mu
          ph = 5.802d-1 - 1.571d-1*mu + 3.27d-2*mu*mu - 2.6d-3*mu*mu*mu
          if(ph > eps) ch = chs*Cpng0*9.4d0*(dlog(ht*f1)/dlog(ht*f2))   &
                                                         *((ht*f2)**ph)
          Gammap = 1d0 - 1d1*Rib/(1d0 + ch*dsqrt(dabs(Rib)))
       if(Gammap <= 0d0) Gammap = (1d0 + 4.7d0*Rib)**(-2d0)
        else if(Rib >= 0d0) then                                         !stable, neutral
          Gammap = (1d0 + 4.7d0*Rib)**(-2d0)
        end if
      end if
      Gammap = dlog(ht*f1)/dlog(ht*f2)*Gammap

      drag = dmax1(0d0,Gammap*((1d0 - sigfl)*Cpng0 + sigfl*chnf))
      drag = anint(drag*1d20)*1d-20

      end

! ******************************************************************************
! Mehod 2: Compute ustar using Louis (1979). Use this value to compute z0h and z0q 
! using Jacobson (2005). Use these values to find L using Li (2010). Use Beljaars 
! & Holstag (1991) to solve for a new ustar and the drag coefficients.

      real(kind=8) function newdrag(vflag,taf,stempt,dp)

      use fasst_global

      implicit none

      integer(kind=4),intent(in):: vflag
      real(kind=8),intent(in):: taf,stempt,dp

! local variables
      integer(kind=4):: counter
      real(kind=8):: ht,f1,f2,cdng0,Gammam,ustar1,z0g,z0p,Rib,zeta
      real(kind=8):: zeta0,ustar,ug,alpha,beta,x1,x2,y1,y2,psim1
      real(kind=8):: psim2,psih1,psih2,Rib2,Rib5,zeta2,zeta5,a,b,c,d

      integer(kind=4),parameter:: cmax = 5                              !maximum iterations
      real(kind=8),parameter:: R = 9.5d-1                               !Prandtl number

      counter = 0
      ht = 0d0
      z0g = 0d0
      z0p = 0d0
      Rib = 0d0
      ug = 0d0
      Cdng0 = 0d0
      Gammam = 0d0
      ustar = 0d0
      ustar1 = 0d0
      alpha = 0d0
      beta = 0d0
      zeta = 0d0
      zeta0 = 0d0
      psim1 = 0d0
      psim2 = 0d0
      psih1 = 0d0
      psih2 = 0d0
      x1 = 0d0
      x2 = 0d0
      y1 = 0d0
      y2 = 0d0
      f1 = 0d0
      f2 = 0d0
      Rib2 = 0d0
      Rib5 = 0d0
      zeta2 = 0d0
      zeta5 = 0d0
      a = 0d0
      b = 0d0
      c = 0d0
      d = 0d0

      ht = iheightn - hm
      ug = uaf

      if(vflag == 0) then
        if (hm <= eps) then
          z0g = rough
        else if (hm > eps) then
          z0g = 7.775d-3
        end if
      else
        z0g = z0l
      end if

! Calculate the Richardson Number
! Rib < 0 - unstable; Rib = 0 - neutral stability; Rib > 0 - stable
      if(dabs(taf-stempt) <= eps.or.ug <= eps) then
        Rib = 0d0
      else
        Rib = 2d0*grav*ht*(taf - stempt)/(ug*ug*(taf + stempt))
      end if

! calculate z0h and z0q; use Louis(1979) to calculate ustar
      f1 = 1d0/z0g
      Cdng0 = vK*vK/(dlog(ht*f1)*dlog(ht*f1))

      if(Rib < 0d0) then                                                 !unstable
        c = 7.4d0*Cdng0*9.4d0*dsqrt(ht*f1)
        Gammam = 1d0 - 9.4d0*Rib/(1d0 + c*dsqrt(dabs(Rib)))       
      else if(dabs(Rib) < eps) then                                      !neutral
        Gammam = 1d0
      else if(Rib > eps) then                                            !stable
        Gammam = 1d0/((1d0 + 4.7d0*Rib)*(1d0 + 4.7d0*Rib))
      end if

      ustar1 = dsqrt(Cdng0*Gammam)*ug

      do while(counter <= cmax)

      if(dabs(ustar1) <= eps) then
        z0p = z0g
      else
        z0p = dp/(4d-1*ustar1)                                            !FAM2dEd
        if(z0p >= z0g) z0p = 0.9d0*z0g
      end if
      f2 = 1d0/z0p

! Calculate zeta = z/L using Li (2010)
      alpha = dlog(ht*f1)
      beta = dlog(ht*f2)

      if(Rib <= 0d0) then                                                !unstable and neutral
        if (Rib < -5d0) then                                             !unstable, Li eqn. doesn't work
          Rib2 = -2d0
          zeta2 = 4.5d-2*alpha*Rib2*Rib2 + Rib2*((3d-3*beta + 5.9d-3)   &
                    *alpha*alpha + (-8.28d-2*beta + 8.845d-1)*alpha +   &
                       (1.739d-1*beta*beta - 9.213d-1*beta - 1.057d-1))
          Rib5 = -5d0
          zeta5 = 4.5d-2*alpha*Rib5*Rib5 + Rib5*((3d-3*beta + 5.9d-3)   &
                    *alpha*alpha + (-8.28d-2*beta + 8.845d-1)*alpha +   &
                       (1.739d-1*beta*beta - 9.213d-1*beta - 1.057d-1))

          zeta = ((zeta5 - zeta2)/(Rib5 - Rib2))*(Rib + Rib2) + zeta2      

        else if(Rib < 0d0.and.Rib >= -5d0) then                          !unstable, Li eqn. works
          zeta = 4.5d-2*alpha*Rib*Rib + Rib*((3d-3*beta + 5.9d-3)*alpha &
                  *alpha + (-8.28d-2*beta + 8.845d-1)*alpha +           &
                       (1.739d-1*beta*beta - 9.213d-1*beta - 1.057d-1))
        else                                                             !neutral
          zeta = 0d0
        end if

! Hogstrom (1995)
        zeta0 = z0p*zeta/ht
        x1 = (1d0 - 1.9d1*zeta)**2.5d-1
        psim1 = dlog((1d0 + x1*x1)*(1d0 + x1)*(1d0 + x1)/8d0)           &
                                              - 2d0*datan(x1) + pi*5d-1
        x2 = (1d0 - 1.9d1*zeta0)**2.5d-1
        psim2 = dlog((1d0 + x2*x2)*(1d0 + x2)*(1d0 + x2)/8d0)           &
                                              - 2d0*datan(x2) + pi*5d-1

        y1 = (1d0 - 1.16d1*zeta)**5d-1
        psih1 = 2d0*dlog((1d0 + y1)*5d-1)
        y2 = (1d0 - 1.16d1*zeta0)**5d-1
        psih2 = dlog((1d0 + y2)*5d-1)

      else if(Rib > 0d0) then                                            !stable
        if(Rib <= 2d-1) then                                             !weekly stable
          zeta = Rib*Rib*((5.738d-1*beta - 4.399d-1)*alpha +            &
                  (-4.901d0*beta + 5.25d1)) + Rib*((-5.39d-2*beta       &
                           + 1.54d0)*alpha + (-6.69d-1*beta - 3.282d0))
          zeta = dmax1(0d0,zeta)
        else if(Rib > 2d-1) then                                         !strongly stable
          zeta = (7.529d-1*alpha + 1.494d1)*Rib + 1.569d-1*alpha        &
                                              - 3.091d-1*beta - 1.303d0
        end if

! Beljaars & Holstag (1991)
        zeta0 = z0p*zeta/ht
        a = 1d0
        b = 6.667d-1
        c = 5d0
        d = 3.5d-1
        psim1 = -a*zeta - b*(zeta - c/d)*dexp(-d*zeta) - b*c/d
        psim2 = -a*zeta0 - b*(zeta0 - c/d)*dexp(-d*zeta0) - b*c/d

        psih1 = -(1d0 + 2d0*a*zeta/3d0)**(1.5d0) - b*(zeta - c/d)*      &
                                            dexp(-d*zeta) - b*c/d + 1d0
        psih2 = -(1d0 + 2d0*a*zeta0/3d0)**(1.5d0) - b*(zeta0 - c/d)*    &
                                           dexp(-d*zeta0) - b*c/d + 1d0
      end if

! Calculate ustar and the coefficients using Li (2012)
      newdrag = 0d0
      if(ug > 1d-3.or.dabs(Rib) >= eps) then
        ustar = ug*vK/(dlog(ht*f2) - psim1 + psim2)
        newdrag = dmax1(0d0,(ustar/ug)*(vK/R)/(dlog(ht*f2) - psih1      &
                                                              + psih2))
      else
        ustar1 = ustar
        newdrag = (vK/R)/(dlog(ht*f2) - psih1 + psih2)
      end if
      newdrag = anint(newdrag*1d20)*1d-20

        if(dabs(ustar-ustar1) <= 1d-3) then
          counter = cmax + 1
        else
          counter = counter + 1
          ustar1 = ustar
        end if
      end do

      end

module module_snow

      use fasst_global

  contains

      subroutine snow_model(sprint,water_flag,first_snow,sdensi,sdf,    &
                            phie1)
    
!     *********************************************************************
!     The routine snmelt was developed by Mary Albert and orginally 
!     coded in basic by G. Krajeski.  The Fortran version was coded  by 
!     G. Koenig.  Snmelt determines the water flow in a snow pack under the
!     influence of gravity for a bulk, or homogenuous, snow pack.  The model
!     also predicts the depth of the freeze front in the snow, the snow depth,
!     and the water outflow at the bottom of the snow pack.  A technical manual
!     is available from Mary Albert

!     The model requires both initialization information and meteorological 
!     information.

!     Variables saved  are listed below

!     Variable               Remarks
!  n(unitless)        the power of the effective saturation is raised to in the relationship between
!                       the permeability and the saturation (range 2.7 to 3.3)
!  phie(unitless)     effective porosity
!  sd(cm)             snow depth
!  swe(cm)            snow water equivalent
!  swi(unitless)      irreducible water saturation of snow pack
!  xf(cm)             amount of the wave that freezes this time step
!  a(cm)              current melt depth
!  at(cm)             net surface balance in freeze/melt equivalent units
!  waterheld(cm)      amount of water held in snow pack
!  timestp(sec)       model time step
!  rain2(cm)          precipitation added to the snow pack
!  maxwet(cm)         depth of wet snow in snow pack
!  rho(g/cm**3)       ice density in snow pack
!  tottim(sec)        total model time
!  maxwetold(cm)      depth of wet snow pack from previous time step
!  vd(cm^3)           average volume of snow crystals in dry pack
!  vw(cm^3)           average volume of snow crystals in wet pack
!  sdold(cm)          depth of snow from previous time step
!  kappa(unitless)    numerical factor used in model
!  wave( ,1) [sec]    time the wave orginated-used to get exit time
!  wave( ,2) [cm]     snow depth at the time of wave creation
!  wave( ,3) [cm]     current water volume of the wave
!  wave( ,4) [sec]    time it will take the wave to exit
!  wave( ,5) [sec]    time required to freeze the wave (-1 implies infinite time)
!  smwaves(unitless)  number of snow melt waves in the snow pack
!  tims(unitless)     number of time steps in refreezing
!  d(cm)              average crystal diameter
!  k(cm^2)            intrinsic permeablility
!  temphist(C)        average temp over present period of refreezing

!     Snmelt History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem
!     Sept 2001 ...................Modified for inclusion in FASST

! calls the following subroutines:
!     inivarivals (appended to this subroutine)
!     getwave (appended to this subroutine)
!     shoveout (appended to this subroutine)
!     bottom (appended to this subroutine)
!     checkwaves (appended to this subroutine)

      IMPLICIT NONE

      integer(kind=4),intent(in):: water_flag,sprint
      integer(kind=4),intent(inout):: first_snow
      real(kind=8),intent(in):: sdensi
      real(kind=8),intent(inout):: sdf
      real(kind=8),intent(out):: phie1

! saved variables
      integer(kind=4):: smwaves,timestp,tims
      real(kind=8):: n,phie,sd,swe,swi,xf,a,at,waterheld,maxwet,rho
      real(kind=8):: maxwetold,vd,vw,sdold,kappa,wave(100,5),rain2
      real(kind=8):: d,k,tottim,temphist,dend,sph

      save:: n,phie,sd,swe,swi,xf,a,at,waterheld,maxwet,rho,maxwetold
      save:: vd,vw,sdold,kappa,wave,rain2,d,k,tottim,temphist
      save:: smwaves,timestp,tims,dend,sph

! local variables
      integer(kind=4):: i,j,year,nsnow
      real(kind=8):: sdoldi,v,abot,qb_f1,d1,atop,dates
      real(kind=8):: abot1,meta,add,dwind,wet


! zero-out local variables
      year = 0
      nsnow = 0
      sdoldi = 0d0
      v = 0d0
      abot = 0d0
      qb_f1 = 0d0
      d1 = 0d0
      atop  = 0d0
      dates = 0d0
      abot1 = 0d0
      meta = 0d0
      add = 0d0
      dwind  = 0d0
      wet = 0d0

! Initialize variables
      phie1 = 0d0
      vsmelt = 0d0

      if(first_snow == 0) then
        n = 3.3d0
        phie = 6.5d-1
        if(dabs(hsaccum-spflag) > eps.and.                              &
                                          dabs(iswe-spflag) > eps) then
          sd = hsaccum*1d2
          swe = iswe*1d2                                                 !snow water equivalent
          if(hsaccum > eps.and.iswe > eps) then
            rho = iswe/hsaccum
          else
            rho = sdensi*1d-3
          end if
          if(swe <= eps.and.sd > eps) swe = sd/3d0
        else if(dabs(hsaccum-spflag) > eps.and.                         &
                                         dabs(iswe-spflag) <= eps) then
          sd = hsaccum*1d2
          rho = sdensi*1d-3
!          swe = sd*rho
          swe = sd/3d0
        else if(dabs(hsaccum-spflag) <= eps.and.                        &
                                          dabs(iswe-spflag) > eps) then
          swe = iswe*1d2
          rho = sdensi*1d-3
!          sd = iswe/rho
          sd = iswe*3d0
        end if
        swi = 2.88d-2
        xf = 0d0
        a = 0d0
        at = 0d0
        waterheld = 0d0
        maxwet = 0d0
        tottim = 0d0
        maxwetold = 0d0
        vd = 0d0
        vw = 0d0
        sdold = sdf*1d2
        kappa = 0d0
        do i=1,50
          do j=1,5
            wave(i,j) = 0d0
          end do
        end do
        rain2 = 0d0
        smwaves = 0
        timestp = int(timstep*3.6d2)                                     !seconds

        temphist = met(iw,ip_tsoil) - Tref                               !C
        d = 0d0
        k = 0d0
        tims = 0
        dend = 0d0
        sph = 0d0

        if(infer_test == 1.and.sdf > eps) then
          sdold = sn_stat(1)
          sd = sn_stat(2)
          rho = sn_stat(3)
          maxwet = sn_stat(4)
          maxwetold = sn_stat(5)
          vd = sn_stat(6)
          vw = sn_stat(7)
          tottim = sn_stat(8)
          temphist = sn_stat(9)
          phie = sn_stat(10)
          swi = sn_stat(11)
          swe = sn_stat(12)
          rain2 = sn_stat(13)
          dend = sn_stat(14)
          sph = sn_stat(15)
          tims = sn_istat(1)
          smwaves = sn_istat(2)
        end if
      end if

      if(dabs(sd) <= eps) then
        smwaves = 0
        tottim = 0d0
      end if

      if(aint(met(iw,ip_pt2)) == 3.or.aint(met(iw,ip_pt)) == 3)         &
        nsnow = 1
      tottim = tottim + float(timestp)                                   !total model run time in seconds

!     Initialize variables used by the model
      call inivarivals(timestp,nsnow,n,phie,sd,swi,sdold,maxwet,        &
                       maxwetold,rho,vd,vw,dend,sph,kappa,d,k)

!     This is the start of the main physics model
      if(water_flag == 0) then
        qb_f1 = ((grthcond(nnodes) + grthcond(nnodes-1))*5d-1)/         &
                      ((nz(nnodes) - nz(nnodes-1))*5d-1)                 !W/m^2 K
!        qb_f1 = km*grthcond(nnodes)/(grthcond(nnodes)*sd                 &
!                                                       + km*nz(nnodes))

      else
        if(met(iw,ip_hi) > eps) qb_f1 = ithcond/met(iw,ip_hi)
        if(met(iw,ip_hi) > eps)                                         &
                    qb_f1 = ithcond*km/(km*met(iw,ip_hi) + ithcond*sd)
      end if

      call getwave(water_flag,timestp,smwaves,tims,tottim,kappa,n,      &
                   sdold,qb_f1,d,dend,sph,rho,swi,phie,swe,sd,rain2,    &
                   waterheld,xf,a,temphist,maxwet,at,wave,abot,meta,    &
                   add,dwind,atop)

      if(abot > eps) then
        abot1 = abot*rho                                                 !cm swe
      else
        abot1 = 0d0
      end if

!     For each time step the model generates a wave (volume of 
!     water) that will move through the snow pack.  Successive
!     waves can over take preceding waves and in do so the waves
!     are combined.  Once the entire snow pack has reach swi the
!     wave moves through the snow pack and appears as outflow(assuming
!     there is no freezing). The subroutine shoveout will determine 
!     what waves will exit the snow pack and the flux of the outflow

      call shoveout(timestp,smwaves,kappa,n,maxwet,d,v,tottim,rho,swi,  &
                    phie,swe,sd,rain2,wave)

! print out snow information
      if(sprint == 1) then
        dates = met(iw,ip_doy) + met(iw,ip_hr)/24d0                     &
                                           + met(iw,ip_min)/(24d0*60d0)

!        if(dabs(sd) <= eps) d = 0d0
        if(sd > eps) wet = maxwet/sd

        write(55,91) int(met(iw,ip_year)),dates,sdold,add,meta,dwind,   &
                     atop,abot,sd,sd-sdold,meltfl,d,swe,wet
   91   format(i10,2x,f10.3,8(f15.8),a5,3(f15.8))
      end if

      if(sd < 1d-1) then   !2d-1
        abot1 = abot1 + sd*rho
        sd = 0d0
        rho = 0d0
      end if

      n = anint(n*1d10)*1d-10
      phie = anint(phie*1d10)*1d-10
      sd = anint(sd*1d10)*1d-10
      swe = anint(swe*1d10)*1d-10
      swi = anint(swi*1d10)*1d-10
      xf = anint(xf*1d10)*1d-10
      a = anint(a*1d10)*1d-10
      at = anint(at*1d10)*1d-10
      waterheld = anint(waterheld*1d10)*1d-10
      maxwet = anint(maxwet*1d10)*1d-10
      rho = anint(rho*1d10)*1d-10
      maxwetold = anint(maxwetold*1d10)*1d-10
      vd = anint(vd*1d10)*1d-10
      vw = anint(vw*1d10)*1d-10
      kappa = anint(kappa*1d10)*1d-10
      rain2 = anint(rain2*1d10)*1d-10
      d = anint(d*1d10)*1d-10
      k = anint(k*1d10)*1d-10
      tottim = anint(tottim*1d10)*1d-10
      temphist = anint(temphist*1d10)*1d-10
      dend = anint(dend*1d10)*1d-10
      sph = anint(sph*1d10)*1d-10
      abot1 = anint(abot1*1d10)*1d-10

      do i=1,50
        do j=1,5
          wave(i,j) = anint(wave(i,j)*1d10)*1d-10
        end do
      end do

      sdf = sd*1d-2                                                      !m
      sdf = anint(sdf*1d10)*1d-10
      atopf = atop*1d-2                                                  !m
      atopf = anint(atopf*1d10)*1d-10
      sdens(iw) = rho*1d3                                                !kg/m^3

      first_snow = 1
      dsnow = d*1d-2                                                     !m
      phie1 = phie                                                       !unitless

! melting needs energy -> cools surface
      refreeze = 0d0                                                     !m
      if(sdold > eps) then
        refreeze = 1d-2*(atop + abot)*rho
      else if(sdold <= eps.and.sd > eps) then
        refreeze = 1d-2*(atop + abot)*rho
      end if
      refreeze = aint(refreeze*1d15)*1d-15

      vsmelt = (v + abot1)*1d-2                                          !m
      vsmelt = anint(vsmelt*1d10)*1d-10

     slushy(iw) = 0
     if(sd > eps.and.maxwet/sd > 0.15) slushy(iw) = 1

      sn_stat(1) = sdold
      sn_stat(2) = sd
      sn_stat(3) = rho
      sn_stat(4) = maxwet
      sn_stat(5) = maxwetold
      sn_stat(6) = vd
      sn_stat(7) = vw
      sn_stat(8) = tottim
      sn_stat(9) = temphist
      sn_stat(10) = phie
      sn_stat(11) = swi
      sn_stat(12) = swe
      sn_stat(13) = rain2
      sn_stat(14) = dend
      sn_stat(15) = sph
      sn_istat(1) = tims
      sn_istat(2) = smwaves

      end subroutine snow_model

! ******************************************************************************
      subroutine inivarivals(timestp,nsnow,n,phie,sd,swi,sdold,maxwet,  &
                             maxwetold,rho,vd,vw,dend,sph,kappa,d,k)

!     Initialize required snow parameters
!
!     Variable               Remarks
!     vo(cm**3)        fresh snow crystal volume
!     l(unitless)      fraction of pore volume associated with swi
!     v1(cm**3/sec)    speed of crystal volume growth
!     vav(cm**3)       average crystal volume in entire pack
!     ke(cm^2)         effective permeability

!     inivarivals History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! no subroutines called

      IMPLICIT NONE

      integer(kind=4),intent(in):: timestp,nsnow
      real(kind=8),intent(in):: n,phie,sd,swi
      real(kind=8),intent(inout):: sdold,maxwet,maxwetold,rho,vd,vw
      real(kind=8),intent(inout):: dend,sph
      real(kind=8),intent(out):: kappa,d,k

      real(kind=8):: dcold
      save:: dcold

! local variables
      real(kind=8):: l,v1,ke,vo,vav,f1,f2,tave,dend_old,sph_old
      real(kind=8):: tgrad,tempt,f,g,h,dc,dend_new,sph_new,pwet,dtime
      real(kind=8):: add

      real(kind=8),parameter:: pi = 3.141592654d0

! zero-out local variables
      l = 0d0
      v1 = 0d0
      ke = 0d0
      vo = 0d0
      vav  = 0d0
      f1 = 0d0
      f2 = 0d0
      kappa = 0d0
      d = 0d0
      k = 0d0
      tgrad = 0d0
      tempt = 0d0
      f = 0d0
      g = 0d0
      h = 0d0
      dc = 0d0
      dend_new = 0d0
      sph_new = 0d0
      pwet = 0d0
      dtime = 0d0
      tave = 0d0
      dend_old = 0d0
      sph_old = 0d0
      add = 0d0

      dtime = float(timestp)/(6d1*2.4d1)                                !time step in days

      if(sd > eps) tgrad = dabs(toptemp - stt(nnodes))/(sd*1d-2)        !temperature gradient (K/m)
      tave = 5d-1*(toptemp + stt(nnodes))                               !average snow temperature (K)
      if(sd > eps) pwet = 1d2*maxwet*swi/sd                             !liquid water content of snow (unitless)

      if(dabs(rho) <= eps) rho = 5d-1*(sdensd + sdensw)*1d-3            !initialize density to middle of allowed range (g/cm^3)
      vo = 5d-5

! Compute the volume of the snow crystals in the wet and dry
! portions of the snow pack based on changes in the snow 
! depth and depth of the wet portion of the snow pack from
! the previous time step. Thus, there is no real crystal growth but
! we shift the crystals from the wet to the dry or the dry to the wet
! portion of the pack based on snowdepth and wetness

      if(maxwet > eps) f1 = 1d0/maxwet
      if(sd-maxwet > eps) f2 = 1d0/(sd - maxwet)

      if(maxwet > maxwetold) then                                        !increase in the depth of wet snow
        vw = vw*(maxwetold*f1) + vd*((maxwet - maxwetold)*f1)
        if(sd > sdold) then                                              !increase in snow depth this time step
          if(sd > maxwet) then                                           !snow pack is not entirely wet
            vd = vd*((sdold - maxwet)*f2) + vo*((sd - sdold)*f2)
          end if
        end if
      else if(maxwet < maxwetold) then                                   !decrease in depth of wet snow
        if(sd > sdold) then                                              !increase in snow depth this time step
          if(sd > maxwet) then                                           !snow pack is not entirely wet
            vd = vo*((sd - sdold)*f2) + vd*((sdold - maxwetold)*f2)     &
                                         + vw*((maxwetold - maxwet)*f2)
          end if
        else
          if(sd > maxwet) then                                           !snow pack is not entirely wet
            vd = vd*((sd - maxwetold)*f2) + vw*((maxwetold - maxwet)*f2)
          end if
        end if
      else
        if(sd > sdold) then                                              !increase in snow depth this time step
          if(sd > maxwet) then                                           !snow pack is not entirely wet
            vd = vo*((sd - sdold)*f2) + vd*((sdold - maxwet)*f2)
          end if
        end if
      end if

      if(maxwet >= sd) vd = 0d0  !vo                                           !maxwet should not be > than sd
      if(maxwet <= eps) vw = 0d0
      sdold = sd
      maxwetold = maxwet

! The grain growth equations are from Brun: Wet-snow metamorphism and LWC, 
! Annals of Glaciology 13 1989, pp.22-
      v1 = 1.28d-11   !(1.25)
      vd = vd + v1*float(timestp)                                        !cm^3
      if(pwet > eps) then
!        l = (swi/rho)*1d2
        l = pwet
        v1 = 1.28d-11 + 4.22d-13*(l**3d0)                                !cm^3/s
        vw = vw + v1*float(timestp)                                      !cm^3
      end if

      if (sd > eps) then
        f1 = 1d0/sd
        vav = vw*(maxwet*f1) + vd*((sd - maxwet)*f1)
        if(vav < 0d0) vav = 0d0  !vo
      else
        vav = 0d0  !vo
      end if

      if(vav > eps) d = 2d0*(((3d0/(4d0*pi))*vav)**(1d0/3d0))

      if(7.8d0*rho > 5d1.or.d <= 0d0) then
        k = 7.7d-2
      else
        k = 7.7d-2*(d**2d0)*dexp(-7.8d0*rho)
      end if

      if(dabs(n) > 0d0) f1 = 1d0/n
      if(k > eps.and.dabs(f1) > eps) ke = (k**f1)/phie
      if(dabs(f1) > eps) kappa = n*ke*(5.47d4**f1)                       !constant used in the model

! grain shape from Crocus (Vionette et al. 2012, Geoscientific Model Development (5), pp.773-791)
      if(nsnow == 1) then
        dend_new = dmin1(dmax1(1.29d0 - 1.7d-1*met(iw,ip_ws),2d-1),1d0)  !dendricity of new snow
        sph_new = dmin1(dmax1(8d-2*met(iw,ip_ws) + 3.8d-1,5d-1),9d-1)    !sphericity of new snow
        add = (dmet1(iw,7) + dmet1(iw,6))*1d2*dcos(sloper)
      end if

      f1 = dexp(-6d3/tave)
      dend_old = dend
      sph_old = sph
      if(pwet > eps) then                                                !wet snow
        if(dend_old > eps) then
          if(sph < 1d0) then                                             !dendritic snow
            dend = dend_old - dtime*pwet*pwet*pwet/1.6d1
            sph = sph_old
          else
            dend = dend_old
            sph = sph_old + dtime*pwet*pwet*pwet/1.6d1
          end if
        else                                                             !non-dendritic snow
          if(sph < 1d0) then
            sph = sph_old + dtime*pwet*pwet*pwet/1.6d1
          else
            sph = sph_old
          end if
        end if
      else                                                               !dry snow
        if(dend_old > eps) then
          if(tgrad <= 5d0) then
            dend = dend_old - 2d8*f1*dtime
            sph = sph_old + 1d9*f1*dtime
          else
            dend = dend_old - 2d8*f1*(tgrad**4d-1)*dtime
            sph = sph_old + 2d8*f1*(tgrad**4d-1)*dtime
          end if
        else
          dend = 0d0
          if(tgrad <= 5d0) then
            sph = sph_old + 1d9*f1*dtime
          else
            if(sph_old > eps) then
              sph = sph_old + 2d8*f1*(tgrad**4d-1)*dtime
            else
              sph = 0d0
            end if
          end if
        end if
      end if
      dend = dmin1(dmax1(dend,0d0),1d0)                                  !dendricity
      sph = dmin1(dmax1(sph,0d0),1d0)                                    !sphericity

      if(sd+add > eps) then
        dend = (dend*sd*rho + dend_new*add)/(sd*rho + add)
        sph = (sph*sd + sph_new*add)/(sd*rho + add)
      else
        dend = 0d0
        sph = 0d0
      end if

      dend = dmin1(dmax1(dend,0d0),1d0)                                  !dendricity
      dend = anint(dend*1d15)*1d-15
      sph = dmin1(dmax1(sph,0d0),1d0)                                    !sphericity
      sph = anint(sph*1d15)*1d-15

! grain size chage due from CROCUS
      tempt = toptemp - Tref

      if(tempt < -4d1) then
        f = 0d0
      else if(tempt >= -4d1.and.tempt > -2.2d1) then
        f = 1.1d-2*(tempt + 4d1)
      else if(tempt >= -2.2d1.and.tempt > -6d0) then
        f = 2d-1 + 5d-2*(tempt + 2.2d1)
      else
        f = 1d0 - 5d-2*tempt
      end if

      if(rho <= 1.5d-1) then
        h = 1d0
      else if(rho > 1.5d-1.and.rho <= 4d-1) then
        h = 1d0 - 4d-3*(rho*1d3 - 1.5d2)
      else
        h = 0d0
      end if

      if(tgrad < 1.5d1) then
        g = 0d0
      else if(tgrad >= 1.5d1.and.tgrad < 2.5d1) then
        g = 1d-2*(tgrad - 1.5d1)
      else if(tgrad >= 2.5d1.and.tgrad < 4d1) then
        g = 1d-1 + 3.7d-2*(tgrad - 2.5d1)
      else if(tgrad >= 4d1.and.tgrad < 5d1) then
        g = 6.5d-1 + 2d-2*(tgrad - 4d1)
      else if(tgrad >= 5d1.and.tgrad < 7d1) then
        g = 8.5d-1 + 7.5d-3*(tgrad - 5d1)
      else
        g = 1d0
      end if

      if(sd <= eps) then
        dc = 0d0
      else
        if(dcold <= eps) dcold = d
        if(sph < 1d0 .and. dend <= eps) then
          dc = dcold + f*h*g*1.0417d-9*float(timestp)*1d2                !cm
          d = dc
        else
          dc = d                                                         !cm
        end if
      end if
      dc = anint(dc*1d15)*1d-15
      dcold = dc


      end subroutine inivarivals

! ******************************************************************************
      subroutine getwave(water_flag,timestp,smwaves,tims,tottim,kappa,  &
                         n,sdold,qb_f1,d,dend,sph,rho,swi,phie,swe,sd,  &
                         rain2,waterheld,xf,a,temphist,maxwet,at,wave,  &
                         abot,meta,add,dwind,atop)

!     This subroutine reads the met data one line at a time (modify
!     at a later date to improve the efficiency of the model) and computes  
!     the energy and mass buget associated with the melt wave for the met 
!     conditions. The impact of the present wave on the snow pact and 
!     the interaction with previous smwaves are calculated
!
!     Variable               Remarks
!     q(w/m**2)     sum of sensible and radiation heat input
!     sd(cm)        snow depth
!     swe(cm)       snow water equivalent
!     rho(g/cm**3)  ice density in snow pack
!     a(cm)         current melt depth
!     datum()       contains the met parameters that have been read
!                   in. The array index gives parameter type. For exampe,
!                   tacol is air temp, tscol is surface temp, wscol is
!                   the wind speed, preccol is the precipitation amt,
!                   ptypecol is the prec type.
!     tottim        total model time
!     wave( ,1)     time the wave orginated-used to get exit time
!     wave( ,2)     snow depth at time of wave creation
!     wave( ,3)     current water volume of the wave
!     wave( ,4)     time it will take the wave to exit
!     wave( ,5)     lifetime of the wave, set to finite value when
!                   there is refreezing. -1 used as an infinite time flag 

!     getwave History
!     April 10 1996................Fortran coding completed
!       Nov 19 1999 .................eleminated a divide by zero problem

! calls the following subroutines:
!     predsndp (appended to this subroutine)
!     satvolume (appended to this subroutine)
!     predictswi (appended to this subroutine)
 
      implicit none

      integer(kind=4),intent(in):: water_flag,timestp
      integer(kind=4),intent(inout):: smwaves,tims
      real(kind=8),intent(in):: tottim,kappa,n,sdold,qb_f1,d,dend,sph
      real(kind=8),intent(inout):: rho,swi,phie,swe,sd,rain2,waterheld
      real(kind=8),intent(inout):: xf,a,temphist,maxwet,at
      real(kind=8),intent(inout):: wave(100,5)
      real(kind=8),intent(out):: abot,meta,add,dwind,atop

! local variables
      integer(kind=4):: wval,d1i
      real(kind=8):: precamtr,precamts,vol,bot,f1,d1,sdo,rain,qbot,atopw
      real(kind=8):: ustart,rhoa,taut,tau,G,b4,ht,z0g,c,Rib,Gammam
      real(kind=8):: ustars,cold,tave,spheats,sumf,t1,f2,f3,n1,Cdng0
      real(kind=8):: ustarc,Mo,driftp

      real(kind=8),parameter:: vk = 0.35d0                               !von Karmen's constant (unitless)


! zero-out local variables
      wval = 0
      d1i = 0
      precamtr = 0d0
      precamts = 0d0
      vol = 0d0
      bot = 0d0
      atopw = 0d0
      f1 = 0d0
      rain = 0d0
      qbot = 0d0
      n1 = 0d0
      ustart  = 0d0
      rhoa = 0d0
      taut = 0d0
      tau = 0d0
      G = 0d0
      b4 = 0d0
      sumf = 0d0
      t1 = 0d0
      ustarc = 0d0
      Mo = 0d0
      driftp = 0d0

! initialize variables              
      atop = 0d0                                                         !cm, top melt depth
      abot = 0d0                                                         !cm, bottom melt depth
      dwind = 0d0                                                        !cm, wind drifted loss
      meta = 0d0                                                         !cm, loss due to metamorphism and pressure
      add = 0d0                                                          !cm, gain due to precip
      sdo = sd

      if(aint(met(iw,ip_pt2)) == 3.or.aint(met(iw,ip_pt)) == 3) then
        precamts = (dmet1(iw,7) + dmet1(iw,6))*1d2*dcos(sloper)
      else if(aint(met(iw,ip_pt)) == 2) then
        precamtr = dmet1(iw,6)*1d2*dcos(sloper)                          !(cm) swe
      end if

! add snow or rain
      if(precamtr > eps) then                                            !have precipitation
        d1i = 2
        call predsndp(d1i,timestp,precamtr,0d0,maxwet,d,rho,swi,phie,   &
                      swe,sd,rain2,sumf,d1)
        rain = anint(d1*1d15)*1d-15
      else if(precamts > eps) then                                       !add new snow to the snow pack
        d1i = 1
        call predsndp(d1i,timestp,precamts,0d0,maxwet,d,rho,swi,phie,   &
                      swe,sd,rain2,sumf,d1)
        add = anint(d1*1d15)*1d-15
      end if

! get reported snow depth or calculate metamorphism and settling
      if(aint(dabs(met(iw,ip_sd)-mflag)*1d5)*1d-5 > eps) then            !reported snow depth
        sd = met(iw,ip_sd)*1d2
        if(dabs(sd) <= eps.and.sdo > eps) sd = sdo                       !!added 3 march 2003
        swe = sd*rho
      else                                                               !snow depth missing make prediction
        d1i = 0
        call predsndp(d1i,timestp,precamts,0d0,maxwet,d,rho,swi,phie,   &
                      swe,sd,rain2,sumf,d1)                              !calculate snow depth due to metamorphism, settling
        meta = -anint(d1*1d15)*1d-15
      end if

! calculate the top melt depth
      f3 = 1d0/(lhfus*rho*1d3)
      tave = (toptemp + met(iw,ip_tsoil))*5d-1
      d1i = 2
      cold = dmax1(0d0,-spheats(tave,d1i)*1d3*(1d-2*sd*rho)             &
                                                        *(tave - Tref))  !J/m^2
      cold = 1d2*cold*f3                                                 !cm

      if(sd > 0d0.and.melt(iw) > 0d0) then  !for old settings: deltat = varied, maxiter = 30
!        if(maxwet/sd >= 6.25d-2) then !6.5 6.25, 7
!          f1 = dmax1(1d0,maxwet/sd) !0.8d0
!        else
!          f1 = 0.5d0    !0.55d0 0.45d0
!        end if
         f1 = 0.85d0 !0.9d0
f1 = 1d6
!!        melt(iw) = dmin1(melt(iw),(km/sd)*tmelt(iw)*float(stepi))

        if(meltfl == 's'.and.tmelt(iw) < Tref) then
          atop = 1d2*float(timestp)*melt(iw)/(lhsub*rho*1d3)             !cm top sublimation depth
          atop = dmin1(sd,atop,f1*float(timestp)/3.6d2)
          atopw = 0d0
        else if((meltfl == 'm').or.tmelt(iw) >= Tref) then  !.and.cold <= eps
          atop = 1d2*float(timestp)*melt(iw)*f3                          !cm top melt depth
          atop = dmin1(sd,atop,f1*float(timestp)/3.6d2)
          atopw = atop*rho                                               !cm swe
        end if
        atop = anint(atop*1d15)*1d-15

        sd = dmax1(0d0,sd - atop)
      end if
      a = atopw + precamtr
      at = atopw

! calculate the bottom melt depth
      if(sd > 0d0) then
        if(water_flag == 0.and.met(iw,ip_tsoil) > Tref) then
!!          qbot = dmax1(0d0,qb_f1*(met(iw,ip_tsoil) - Tref))              !W/m^2
!          qbot = dmax1(0d0,qb_f1*(met(iw,ip_tsoil) - toptemp))            !W/m^2
          qbot = dmax1(0d0,qb_f1*(stt(nnodes-1) - stt(nnodes)))            !W/m^2
        else
          qbot = 0d0
        end if

        f1 = 1d6
        abot = dmin1(sd,1d2*float(timestp)*qbot*f3,                     &
                                               f1*float(timestp)/3.6d2)  !cm
        abot = anint(abot*1d15)*1d-15                                    !cm
        sd = dmax1(0d0,sd - abot)
      end if

! wind transport, based on Jordan & Andreas, 1999 JGR 104(C4), pp.7785-7806
! parameters picked for snow on seaice
      if(sd > eps) then
        n1 = 12.91d0 - 24.52d0*phie + 11.88d0*phie*phie
        ustart = 0.1d0 + 4d0*dsqrt((1d0 - phie)*n1*(0.2d0))              !m/s
        ustart = ustart

        ht = iheightn - (hsaccum + hi + newsd)
        do while(ht <= 0d0)
          ht = ht + 5d-1
        end do
        z0g = 7.775d-3                                                   !(0.05 - 1.5mm)

! calculate the Richardson Number
        if(dabs(dmet1(iw,4)-toptemp) <= eps.or.dmet1(iw,3) <= eps) then
          Rib = 0d0 !0.21d0
        else
          Rib = 2d0*grav*ht*(dmet1(iw,4) - toptemp)/                    &
                                             (dmet1(iw,3)*dmet1(iw,3)*  &
                                              (dmet1(iw,4) + toptemp))
        end if

! calculate z0h and z0q; use Louis(1979) to calculate ustar
        f2 = ht/z0g
        Cdng0 = vk*vk/(dlog(f2)*dlog(f2))

        if(Rib < 0d0) then                                               !unstable
          c = 7.4d0*Cdng0*9.4d0*dsqrt(f2)
          Gammam = 1d0 - 9.4d0*Rib/(1d0 + c*dsqrt(dabs(Rib)))
        else                                                             !stable
          Gammam = 1d0/((1d0 + 4.7d0*Rib)*(1d0 + 4.7d0*Rib))
        end if
      
        ustars = dsqrt(Cdng0*Gammam)*dmet1(iw,3)

        rhoa = 3.48d-3*(met(iw,ip_ap)*1d2/dmet1(iw,4))                   !kg/m^3
        taut = ustart*ustart*rhoa                                        !kg/m*s^2
        tau = ustars*ustars*rhoa                                         !kg/m*s^2
        G = (5d0/(dsqrt(rhoa)*grav))*(dsqrt(tau) - dsqrt(taut))         &
                                                          *(tau - taut)  !kg/m*s

        if(rho*1d3 < 4d2.or.dabs(rho*1d3 - 4d2) > 5d1) then
          b4 = 1d0
        else
          b4 = dexp(-4.6d-2*(rho*1d3 - 4d2))                             !s/m
        end if

! from Crocus
        f1 = 1.25d0 - 4.2d-3*(dmax1(sdensd,rho*1d3) - sdensd)
        if(dend > eps) then
          Mo = 3.4d-1*(7.5d-1*dend - 5d-1*sph) + 6.6d-1*f1
        else
          Mo = 3.4d-1*(-5.83d-1*(d*1d1) - 8.33d-1*sph + 8.33d-1)        &
                                                            + 6.6d-1*f1
        end if
        Mo = dmax1(Mo,0d0)
        driftp = -2.868d0*dexp(-8.5d-2*met(iw,ip_ws)) + 1d0 + Mo
        if(driftp > eps) ustarc = -dlog((1d0 + Mo)/2.868d0)/8.5d-2

        t1 = 4d-2*dmet1(iw,4) + 8.84d-2*(rho*1d3)*grav*1d-2
!        if(ustars >= ustart.and.t1 < 5d1) then
!        if(ustars >= ustarc.and.t1 < 5d1) then
        if((met(iw,ip_ws) >= ustarc.and.driftp > eps).and.t1 < 5d1) then
!          dwind = sd*2.66d-3*b4*G*dexp(-t1)                              !m/s
          dwind = sd*2.66d-2*b4*G*dexp(-t1)                              !m/s
        else
          dwind = 0d0
        end if
        dwind = 1d2*dwind*float(timestp)                                 !cm
        dwind = anint(dwind*1d15)*1d-15

        sd = dmax1(0d0,sd - dwind)
      end if
      sd = anint(sd*1d15)*1d-15

! Calculate any freezing or refreezing
      call satvolume(timestp,smwaves,tims,rho,tottim,at,kappa,n,sdold,  &
                     xf,a,temphist,maxwet,sd,swi,waterheld,wave)

! calculate irreducible liquid water in the snow pack
      d1i = 1
      call predictswi(d1i,timestp,0d0,sd,xf,tottim,at,sdold,maxwet,swi,&
                      a,waterheld)

!     if(a > 0d0) then                                                    !water in snow pack
      if(a > 1d-4) then                                                  !water in snow pack
        smwaves = smwaves + 1
        wave(smwaves,1) = -tottim
        wave(smwaves,2) = sd
        wave(smwaves,3) = a

        if(sd > 0d0) then                                                !snow on the ground
          vol = a   !atop
          wval = smwaves
          call bottom(wval,timestp,vol,kappa,n,sd,wave,bot)
          wave(smwaves,4) = bot                                          !calculate time for wave to 
        else                                                             !get to the bottom of the pack
          wave(smwaves,4) = 0d0
        endif

        wave(smwaves,5) = -1d0                                           !flag indicating infinite duration to reach bottom of pack 
        if(smwaves > 1)                                                 &
          call checkwaves(timestp,smwaves,kappa,n,tottim,sd,wave)        !combine smwaves if necessary
      else if(a > 0d0.and.a <= 1d-4.and.smwaves >= 2) then
        wave(smwaves-1,3) = wave(smwaves-1,3) + a
      end if

      end subroutine getwave

! ******************************************************************************
      subroutine checkwaves(timestp,smwaves,kappa,n,tottim,sd,wave)

!     This routine checks to see of a new wave will overtake any 
!     of the preceding  waves.  If wave n overtakes wave n-1 the 
!     two waves are combined.  The combined wave is checked to see 
!     if it overtakes wave n-2 and so on.  Thus, at any moment in time
!     the numbering order will be desending from the top of the 
!     snow pack to the bottom
!
!     Variable               Remarks
!     wave( ,1)     time the wave orginated-used to get exit time
!     wave( ,2)     snow depth at time of wave creation
!     wave( ,3)     current water volume of the wave
!     wave( ,4)     time it will take the wave to exit
!     wave( ,5)     lifetime of the wave, set to infinite (-1)value when
!                   when there is refreezing
!     checkwave History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: timestp
      integer(kind=4),intent(inout):: smwaves
      real(kind=8),intent(in):: kappa,n,tottim,sd
      real(kind=8),intent(inout):: wave(100,5)

! local variables
      integer(kind=4):: w1,w2,d4

      w1 = 0
      w2 = 0
      d4 = 1

      if(smwaves > 1) then
        do while(wave(smwaves-1,4)-wave(smwaves-1,1)                    &
                                     > wave(smwaves,4)-wave(smwaves,1))
          w1 = smwaves - 1
          w2 = smwaves

          call collide(w1,w2,d4,timestp,smwaves,kappa,n,tottim,sd,wave)
          if(smwaves <= 1) exit
        end do
      end if

      end subroutine checkwaves

! ******************************************************************************
      subroutine shoveout(timestp,smwaves,kappa,n,maxwet,d,v,tottim,rho,&
                          swi,phie,swe,sd,rain2,wave)

!     Determine the flow out of the pack for this time step

!     Variable               Remarks
!     h1 & h2       lower & upper bounds to volume flux integration
!     flipper       set to -1 if smwaves should collide
!     P1-5          parts of the volume flux integration solution
!     tb            time shoveout begins
!     timend        if -1 current wave has ended 
!     smwaves       number of smwaves
!     tb            model time for this call to shoveout
!     wave( ,1)     time the wave orginated-used to get exit time
!     wave( ,2)     snow depth at the time of wave creation
!     wave( ,3)     current water volume of the wave
!     wave( ,4)     time it will take the wave to exit
!     wave( ,5)     time required to freeze the wave
!     wd            number of smwaves this routine has processed
!     t1            lower time boundary for flux integration
!     smwavesout    number of smwaves that will have some outflow this time step
!     v             volume of water coming out the pack this time step

! calls the following subroutines:
!     removewave (appended to this subroutine)
!     collide (appended to this subroutine)
!     predsndp (appended to this subroutine)

      IMPLICIT NONE

      integer(kind=4),intent(in):: timestp
      integer(kind=4),intent(inout):: smwaves
      real(kind=8),intent(in):: kappa,n,maxwet,d
      real(kind=8),intent(inout):: v,tottim,rho,swi,phie,swe,sd,rain2
      real(kind=8),intent(inout):: wave(100,5)

! local variables 
      integer(kind=4):: wd,ctr,smwavesout,timend,flipper,d2,d3,d4,d1i
      real(kind=8):: tb,t1,h1,h2,p1,p2,p3,p4,p5,d1,sumf

! zero-out local variables
      wd = 0
      ctr = 0
      d1i = 0
      smwavesout = 0
      timend = 0
      flipper =0
      d2 = 0
      d3 = 0
      d4 = 0
      d1i = 0
      tb = 0d0
      t1 = 0d0
      h1 = 0d0
      h2 = 0d0
      p1 = 0d0
      p2 = 0d0
      p3 = 0d0
      p4 = 0d0
      p5 = 0d0
      d1 = 0d0
      sumf = 0d0

      if((smwaves > 0).and.(sd > 0d0)) then                              !smwaves & snow present
        tb = tottim                                                      !start time

! Determine number of smwaves coming out of the pack(smwavesout)this time step
        do ctr=1,smwaves
          if(wave(ctr,4)-wave(ctr,1) < tottim+float(timestp)) then
            smwavesout = smwavesout + 1
          end if
        end do

! determine how much water will flow out of the pack for each wave
        do while(tottim < (tb+float(timestp)).and.smwavesout > wd)
          wd = wd + 1                                                    !number of wave being processed by this routine

! set lower time boundary for integration to get flux out this time step
          if(tottim >= (wave(1,4)-wave(1,1))) then    
            t1 = tottim
          else
            t1 = wave(1,4) - wave(1,1)     
!             set t1 time of integration based on location of the lowest 
!             wave relative to the bottom of the pack
          end if

          if(wd == smwavesout) then                                      !is this last wave out
            if(dabs(wave(1,5)+1d0) <= eps) then                          !if this wave is going to last for ever
              tottim = tb + float(timestp)
            else                                                         !wave is going last longer than this timestep
              if((wave(1,5)-wave(1,1)) > (tb+float(timestp))) then
                tottim = tb + float(timestp)
              else
                tottim = wave(1,5) - wave(1,1)
                timend = -1
              end if
            end if
            flipper = 0
          else                                                           !not last wave
            if(dabs(wave(1,5)+1d0) <= eps) then                          !no freezing involved
              tottim = wave(2,4) - wave(2,1)                             !time second wave comes out
            else                                                         !freezing
              if(wave(1,5)-wave(1,1) > wave(2,4)-wave(2,1)) then
                tottim = wave(2,4) - wave(2,1)
              else
                tottim = wave(1,5) - wave(1,1)
                timend = -1
              end if
            end if
            flipper = -1                                                 !when = -1 smwaves will collide
          end if

          if((dabs(kappa) > eps.and.wave(1,2) > 0d0).and.               &
              dabs(n-1d0) > eps) p1 = (wave(1,2)/kappa)**(n/(n - 1d0))
          p2 = 1d0 - n
          p3 = tottim + wave(1,1)
          if(p3 <= eps) p3 = 1d0                                         !added to stop the program from bombing.
          p4 = t1 + wave(1,1)
          if(p4 <= eps) p4 = 1d0                                         !added to stop the program from bombing. This is
                                                                         !a safety check since double precisioning tottim
                                                                         !has solved the problem
          if(dabs(n-1d0) > eps) p5 = 1d0/(1d0 - n)
          if(p3 > eps.and.dabs(p5) > eps) h1 = p1*p2*(p3**p5)
          if(p4 > eps.and.dabs(p5) > eps) h2 = p1*p2*(p4**p5)            !h2 = -wave(#,3)
          v = v + h1 - h2                                                !amount of water flowing out of pack

          d2 = 1
          if(timend /= 0) call removewave(d2,smwaves,wave)
          if(flipper /= 0) then
            if(timend == 0) then
              d3 = 2
              d4 = 2
              call collide(d2,d3,d4,timestp,smwaves,kappa,n,tottim,sd,  &
                           wave)
            end if
          end if
          flipper = 0    
          timend = 0
        end do
        tottim = tb
      else if(smwaves /= 0) then                                         !waves but no snow
        do ctr=1,smwaves
          if(wave(ctr,4)-wave(ctr,1) < tottim) then
          if((dabs(kappa) > eps.and.wave(ctr,2) > 0d0).and.             &
              dabs(n-1d0) > eps) p1 = (wave(ctr,2)/kappa)**(n/(n - 1d0))
            p2 = 1d0 - n
            p3 = tottim + wave(ctr,1)
            if(p3 <= eps) p3 = 1.0                                       !added to stop the program from bombing.
            if(dabs(n-1d0) > eps) p5 = 1d0/(1d0 - n)
            if(p3 > eps.and.dabs(p5) > eps) h1 = p1*p2*(p3**p5)
            v = v - h1
          else
            v = v + wave(ctr,3)
          end if
          call removewave(ctr,smwaves,wave)
        end do
      end if

      d1i = 3
      if(v > 0d0) call predsndp(d1i,timestp,0d0,v,maxwet,d,rho,swi,     &
                                phie,swe,sd,rain2,sumf,d1)

      end subroutine shoveout

! ******************************************************************************
      subroutine predsndp(smode,timestp,newsnow,meltswe,maxwet,d,rho,   &
                          swi,phie,swe,sd,rain2,sumf,d1)

!     Predict the snow depth and snow water equavalent
!     assumes that newsnow is the snow water equivalent depth

!     Variable               Remarks
!     newdense    density of new snow
!     critdense   density at which metamorphic densification stops
!     c1          constant for rate of metamorphic densification
!     c5          constant for rate of settlement for overburden
!     c6          constant for rate of settlement for overburden
!     eda         constant associated with viscosity
!     meta        rate of metamorphic densification
!     avgsnpress  average pressure of snow
         
!     predictsndp History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! no subroutines called

! uses function: dense


      IMPLICIT NONE

      integer(kind=4),intent(in):: smode,timestp
      real(kind=8),intent(in):: newsnow,meltswe,maxwet,d
      real(kind=8),intent(inout):: sumf,rho,swi,phie,swe,sd,rain2
      real(kind=8),intent(out):: d1

! local variables
      integer(kind=4):: d1i
      real(kind=8):: splus,critdense,c1,eda0,c5,c6,c3,c4,newdense,meta
      real(kind=8):: avgsnpress,eda,over,delta,ttemp,dense,f1,sloper
      real(kind=8):: edac,overc,f1c,f2c

! zero-out local variables
      d1i = 0
      splus = 0d0
      critdense = 0d0
      c1 = 0d0
      eda0 = 0d0
      c5 = 0d0
      c6 = 0d0
      c3 = 0d0
      c4 = 0d0
      newdense = 0d0
      meta = 0d0
      avgsnpress = 0d0
      eda  = 0d0
      over = 0d0
      delta = 0d0
      ttemp = 0d0
      d1 = 0d0
      f1 = 0d0
      sloper = 0d0
      sloper = dmax1(0d0,dmin1(1.57d0,slope*pi/1.8d2))

      select case (smode)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case(0)
!     predict snow metamorphic densification

        if(sd > eps) f1 = 1d0/sd

!        critdense = dmin1(1d-1,1.15d0*rho)                               !g/cm^3   [100 kg/m^3]
        critdense = 1d-1
        c1 = 2.778d-6                                                    !1/seconds
        eda0 =  3.6d7  !9d8 !3.6d7                                              !g/cm*s            
        c5 = 1d-1  !8d-2                                                        !1/C
        c6 = 2.3d1                                                       !cm^3/g  
      
        if(sd > eps) then                                                !snow present
          if(rho <= critdense.and.4.6d1*(rho - critdense) < 5d1) then
            c3 = 1d0                                                     !unitless
          else
            c3 = dexp(-4.6d1*(rho - critdense))
          end if
          
!          c4 = dmax1(1d0,dmin1(2d0,1d0 + maxwet*f1))
          c4 = 1d0
          if(maxwet > 0d0) c4 = 2d0

          ttemp = (toptemp + stt(nnodes))*5d-1                           !K

          if(dabs(4d-2*(Tref - ttemp)) < 5d1) then
            meta = dmin1(0d0,-c1*c3*c4*dexp(-4d-2*(Tref - ttemp)))       !(1/s) metamorphic densification
          else
            meta = 0d0
          end if
          meta = anint(meta*1d15)*1d-15

          avgsnpress = 981d0*swe*(2d0/3d0)*dcos(sloper)                  !(g/cm*s**2) snow pressure
!          avgsnpress = 981d0*(sdensw*1d-3 - rho)*sd*dcos(sloper)            !(g/cm*s**2) snow pressure  !*(2d0/3d0)

!          ttemp = 0.5d0*(ttemp + ptemp)

          if(dabs(c5*(Tref - ttemp) + c6*rho) < 5d1) then
            eda = eda0*dexp(c5*(Tref - ttemp) + c6*rho)                  !(g/cm*s) snow viscosity
            over = dmin1(0d0,-avgsnpress/eda)                            !(1/s) overburden settling
          else
            over = 0d0
          end if
          over = anint(over*1d15)*1d-15

! crocus
          c5 = 1d-1                                                    !1/C
          c6 = 2.3d1                                                   !cm^3/g
          if(dabs(c5*(Tref - ttemp) + c6*rho) < 5d1) then
            f1c = 1d0/(1d0 + 6d1*maxwet*swi*f1)
            f2c = dmin1(4d0,dexp(1d1*dmin1(4d-1,d*1d1 - 2d-1)))
            eda0 = 7.62237d6                                           !kg/m*s
            edac = 1d1*eda0*f1c*f2c*rho/2.5d-1                         !g/cm*s
            eda = edac*dexp(c5*(Tref - ttemp) + c6*rho)                !g/cm*s
            overc = dmin1(0d0,-avgsnpress/eda)                         !(1/s) overburden settling
          else
            overc = 0d0
          end if
          overc = anint(overc*1d15)*1d-15

          delta = sd*(over + meta)*float(timestp)                        !cm
!delta = sd*(overc + meta)*float(timestp)
          d1 = delta                                                     !cm

!     New snow depth due to densification and overburden
          sd = sd + delta                                                !cm

          splus = swi*swe                                                !water required to satisfy swi requirement
!          rho = (swe - splus)/sd
          rho = swe/sd
          if(rho > sdensw*1d-3)then
            rho = sdensw*1d-3
            swi = 0.05d0
!            swe = rho*sd/(1d0 - swi)
            swe = rho*sd
          end if
          phie = (1d0 - rho)*(1d0 - swi)
        end if

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case(1)
!     add new snow to snow pack and calculate new ice density, effective 
!     porosity, snow depth, and snow water equivalent values

! New Snow Density (g/cm^3)
! based on Jordan et al. (1999) JGR 104(C4),p.7785-7806
        d1i = 4
        newdense = dense(dmet1(iw,4),dmet1(iw,3),d1i)*1d-3               !g/cm^3

! crocus
!        newdense = 1.09d2 + 6d0*met(iw,ip_tmp)                          &
!                                           + 2.6d1*dsqrt(met(iw,ip_ws))
!        newdense = 1d-3*dmax1(sdensd,dmin1(sdensw,newdense))

! reference: Anderson, 1976
!      if(ptemp > 258.16d0) then
!       newdense = 50d0 + 1.7d0*(ptemp - 258.16)**1.5d0                   !kg/m^3
!     else
!       newdense = 5d1
!     end if
!     newdense = newdense*1d-3                                            !g/cm^3

! reference: Hedstrom & Pomeroy (1998) Hdr. Proc. 12, p.1611-1625
!     f1 = dmet1(iw,4) - Tref
!     if(f1 < 0d0.and.dabs(f1/2.6d0) < 5d1) then
!       newdense = (67.9d0 + 51.3d0*dexp(f1/2.6d0))*1d-3
!     else
!       newdense = (119.2d0 + 20d0*f1)*1d-3
!     end if

! reference: meted.ucar.edu/nwp/pcu2/etsnow2c.htm; Eta snow initial conditions
!     newdense = 2d-4*(ptemp - Tref)*(ptemp - Tref) + 
!     &             9.6d-3*(ptemp - Tref) + 0.1495

        f1 = 1d0/newdense
        sumf = newdense
        rho = (rho*sd + newsnow)/(sd + newsnow*f1)
        phie = (1d0 - rho)*(1d0 - swi)
        sd = sd + newsnow*f1
        d1 = newsnow*f1
        swe = swe + newsnow

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case(2)
!     increase in liquid water in pack due to rainfall

        rain2 = newsnow  !rain2 + newsnow

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case(3)
!     decrease in SWE or water in snow pack due to rain runoff

        if(rain2 > meltswe) then                                         !runoff
          rain2 = rain2 - meltswe
        else                                                             !rain held in snow pack-increasing swe
          if(swe > eps) then
            swe = swe + rain2 - meltswe
            rain2 = 0d0
          end if
        end if
        if(swe < 0d0) then
!       if((swe < 0d0).or.(sd < 0d0)) then
          d1 = -sd
          swe = 0d0
          sd = 0d0
          rho = 0d0
        end if

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case default

        error_code = 1
        error_type = 4

!        if(single_multi_flag == 0) write(10,'(''freq_id'',i10,          &
!          '' Invalidmode in predsndp.f; day'',i4,'' hour'',i4)')        &
!          freq_id,int(met(iw,ip_doy)),int(met(iw,ip_hr))

      end select

      if(sd < 0d0) then
        sd = 0d0
        swe = 0d0
        rho = 0d0
      end if

      end subroutine predsndp

! ******************************************************************************
      subroutine satvolume(timestp,smwaves,tims,rho,tottim,at,kappa,n,  &
                           sdold,xf,a,temphist,maxwet,sd,swi,waterheld, &
                           wave)

!     Determines if the snow has reached the irreducible water
!     saturation. If the irreducible water saturation level has 
!     not been reached the routine determines how much of the
!     current water wave  will be used in reaching this level.
!     The routine also calculates any refreezing

!     Variable               Remarks
!     l(J/g)            latent heat of fusion for water
!     kfs(J/(s-cm-K))   thermal conductivity of snow 
!     xfold(cm)         amount that froze the previous time step
!     volgone(cm)       volume of water loss to refreezing for current wave
!     vol(cm)           vol of melt or vol of water that can be frozen due to heat loss
!     heatneed(J/m**2)  heat needed raise refrozen snow to melt temp
!     heatgot(J/m**2)   amount of heat input this time step
!     upneed(cm)        amount of water needed to reach swi

!     satvolume History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! calls the following subroutines:
!     predictswi (appended to this subroutine)
!     removewave (appended to this subroutine)
!     depth (appended to this subroutine)

      IMPLICIT NONE

      integer(kind=4),intent(in):: timestp
      integer(kind=4),intent(inout):: smwaves,tims
      real(kind=8),intent(in):: rho,tottim,at,kappa,n,sdold
      real(kind=8),intent(inout):: xf,a,temphist,maxwet,sd,swi
      real(kind=8),intent(inout):: waterheld,wave(100,5)

! local variables
      integer(kind=4):: ctr,ctr1,d1i
      real(kind=8):: l,p1,p2,p3,p5,volgone,heatneed,heatgot,parm1
      real(kind=8):: xfold,kfs,exdelta,vol,upneed,depthval,t1

! zero=out local variables
      ctr =0
      ctr1 = 0
      d1i = 0
      l = 0d0
      p1 = 0d0
      p2 = 0d0
      p3 = 0d0
      p5 = 0d0
      volgone = 0d0
      heatneed = 0d0
      heatgot = 0d0
      parm1 = 0d0
      depthval = 0d0
      xfold = 0d0
      kfs = 0d0
      exdelta = 0d0
      vol = 0d0
      upneed = 0d0
      t1 = 0d0

      if(rho > eps) then
        kfs = 0.023d0 + (7.75d-5*rho*1d3                                &
                  + 1.105d-6*((rho*1d3)**2d0))*(2.29d0 - 0.023d0)        !W/m*K
        kfs = kfs*1d-2
      else
        kfs = 0.0045d0
      end if
      l = 333.05d0
      p1 = 0d0
      p2 = 0d0
      p3 = 0d0
      p5 = 0d0
      volgone = 0d0
      heatneed = 0d0
      heatgot = 0d0
      upneed = 0d0
      xfold = xf
      vol = a
      t1 = dmet1(iw,4) - Tref

      if(vol <= 0d0) then
        if(t1 < 0d0) then                                               !no melt and freezing conditions
          temphist = (temphist*tims - t1)/(tims + 1d0)
          tims = tims + 1

!     dsqrt(l*swi) is the latent heat necessary to freeze snow
!     at the irreducible saturation level. If xf > xfold
!     freezing can occur and predictswi is called to determine
!     the change in the amount of unfrozen water in the snow pack
          xf = dsqrt((2d0*kfs*temphist*timestp*tims)/(l*swi))
          if(xf > xfold)then
            exdelta = xf - xfold
            d1i = 3
            call predictswi(d1i,timestp,exdelta,sd,xf,tottim,at,sdold,  &
                            maxwet,swi,a,waterheld)

          end if
          xfold = xf

!     If freezing to bottom of the pack remove smwaves that would
!     have exit the bottom of the pack
          if(xf >= sd) then
            do ctr=1,smwaves
              call removewave(ctr,smwaves,wave)
            end do

!     Freezing not to bottom of the pack, but some of the water has
!     frozen-remove smwaves accordingly
          else
            do ctr=smwaves,1,-1
              parm1 = dmax1(0d0,tottim + wave(ctr,1))
              ctr1 = ctr
              call depth(ctr1,parm1,kappa,n,sd,wave,depthval)            !find depth of wave in pack

!     If the depth of the freezing front is greater than the depth of the  
!     wave remove the wave
              if((xf-depthval-sd+wave(ctr,2)) > 0d0) then
                call removewave(ctr,smwaves,wave)

!     calculate the time required to freeze the wave.  If it takes 
!     less time to freeze the wave than it takes the wave to reach
!     the bottom of the pack remove the wave
              else
                if((wave(ctr,2)+xf-sd > 0d0.and.dabs(kappa) > eps.and.  &
                                              dabs(n-1d0) > eps) ) then
                  p1 = ((wave(ctr,2) + xf - sd)/kappa)**(n/(n - 1d0))
                  p2 = n - 1d0
                  p3 = tottim + wave(ctr,1)
                  if(p3 <= eps) p3 = 1d0
                  p5 = 1d0/(1d0 - n)
                  if(dabs(p5) > eps) volgone = p1*p2*(p3**p5)
                  if(volgone > eps)                                     &
                    wave(ctr,5) = (((n - 1d0)/volgone)**(n - 1d0))*     &
                                               ((wave(ctr,2)/kappa)**n)
                  if(wave(ctr,5) < wave(ctr,4)) then
                    call removewave(ctr,smwaves,wave)
                  end if
                end if
                exit
              end if
            end do
          end if
        end if

        vol = 0d0
      else if(vol > 0d0) then
        if(maxwet-sd > eps.and.sd > 0d0) then                            !snow pack wet with no frozen snow
          maxwet = sd
          if((tims /= 0).or.(temphist > eps)) then
            temphist = 0d0
            tims = 0
            xf = 0d0
            xfold = 0d0
          end if
        else
! determine heat needed to raise the temp of the frozen pack to melt level
          if((tims /= 0).or.(temphist > eps)) then
            heatneed = 1.0405d0*rho*xf*temphist    
            heatgot = vol*3330464d0
            if(heatgot > heatneed) then                                  !bring snow to melt temp
              vol = (1d0 - (heatneed/heatgot))*vol
              xf = 0d0
              tims = 0
              temphist = 0d0
            else
              vol = 0d0
              temphist = 0.9611d0*((heatneed - heatgot)/(rho*xf))
            end if
          end if

          if(maxwet < sd*swi) then                                       !pack not completely at swi level
            upneed = (sd - maxwet)*swi
          else
            upneed = 0d0
          end if

          if(vol >= upneed) then                                         !determine how much of the pack can be raised to swi level
            d1i = 2
            call predictswi(d1i,timestp,upneed,sd,xf,tottim,at,sdold,   &
                            maxwet,swi,a,waterheld)

            vol = vol - upneed                                           !remaining volume will leave snow pack
            if((tims /= 0).or.(temphist > eps)) then
              temphist = 0d0
              tims = 0
              xf = 0d0
              xfold = 0d0
            end if
          else
            d1i = 2
            call predictswi(d1i,timestp,vol,sd,xf,tottim,at,sdold,      &
                            maxwet,swi,a,waterheld)
            vol = 0d0
            if((tims /= 0).or.(temphist > eps)) then
              temphist = (temphist*tims - t1)/(tims + 1d0)
              tims = tims + 1
              if(temphist < eps) then
                temphist = 0d0
                tims = 0
                xf = 0d0
                xfold = 0d0
              end if
            end if                                   
          end if
        end if
      end if

      a = vol

      end subroutine satvolume

! ******************************************************************************
      subroutine predictswi(wmode,timestp,exdelta,sd,xf,tottim,at,sdold,&
                            maxwet,swi,a,waterheld)

!     Predicts the irreducible water saturation and volume of 
!     wet snow in the pack

!     Variable               Remarks
!     rainuse(cm)   amount of precipitation from input file
!     exdelta(cm)   is the amount of water that freezes this time step
         
!     predictswi History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! no subroutines called

      IMPLICIT NONE

      integer(kind=4),intent(in):: wmode,timestp
      real(kind=8),intent(in):: exdelta,sd,xf,tottim,at,sdold
      real(kind=8),intent(inout):: maxwet,swi,a,waterheld

! localc variables
      real(kind=8):: rainuse,freezeuse,wetuse,swiold,swidelta,tempconst
      real(kind=8):: rainconst,refreezeconst,heatbalconst,wetconst
      real(kind=8):: runoffconst,f1,omw,swiconst

! zero-out local variables
      rainuse = 0d0
      freezeuse = 0d0
      wetuse = 0d0
      swiold = 0d0
      swidelta = 0d0
      tempconst = 0d0
      omw = 0d0
      rainconst = 0d0
      refreezeconst = 0d0
      heatbalconst = 0d0
      wetconst = 0d0
      swiconst = 0d0
      runoffconst = 0d0

      select case (wmode)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case(1)                                                            !predict swi and adjust maxwet

        if(sd > eps)then
          if(aint(met(iw,ip_pt)) == 2) rainuse = dmet1(iw,6)*1d2

          if(xf > eps.and.sd > eps) then
            freezeuse = dmin1(1d0,dmax1(0d0,xf/sd))                      !fraction of snow with frozen water
          end if

          if(maxwet > eps.and.sd > eps) then
            wetuse = dmin1(1d0,dmax1(0d0,maxwet/sd))                     !fraction of snow at irreducible saturation
          end if

!     Calculate new irreducible water saturation value
          f1 = 1d0/3.6d2
          tempconst = -5d-6*f1
          rainconst = -4d-2*f1
          runoffconst = -1.5d-3*f1
          refreezeconst = -1.5d-4*f1

          swidelta = tempconst*(ptemp - Tref) + rainconst*rainuse +     &
                     refreezeconst*freezeuse + runoffconst*a +          &
                     heatbalconst*at + wetconst*wetuse + swiconst*swi
          swidelta = swidelta*timestp
          swiold = swi
          swi = dmax1(5d-2,dmin1(1d-1,swi + swidelta))

          if((sd > sdold).and.(tottim > (1.5d0*float(timestp)))) then
            swi = dmax1(5d-2,dmin1(1d-1,                                &
                                   (swi*sdold + 5d-2*(sd - sdold))/sd))
          end if

          maxwet = dmax1(0d0,dmin1(waterheld/swi,sd))

! Calculate current melt depth
          if((rainuse > 0d0.or.maxwet > sd-swi).or.                     &
                                         (a < at.or.swi < swiold)) then
            if((maxwet-sd)*swi > at-a) then
              a = at
            else
!              a = dmax1(0d0,dmin1(a + (maxwet - sd)*swi,sd))
            a = dmax1(0d0,a + maxwet - (sd - swi))
            end if
          end if

          if(maxwet-(sd-swi) > eps) then                                 !cal depth of wet snow & water held
            maxwet = dmax1(0d0,sd - swi)
            waterheld = maxwet*swi
          end if
        else
          swi = 5d-2
          waterheld = 0d0
          maxwet = 0d0
        end if

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case(2)                                                            !adding water to the snow pack
       
        waterheld = waterheld + exdelta
        maxwet = dmax1(0d0,dmin1(waterheld/swi,sd))
        if(maxwet-(sd-swi) > eps) maxwet = dmax1(0d0,sd - swi)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case(3) 
!     refreezing occurs or freezing beyond the last time freezing
!     depth occurs-remove water

        omw = maxwet
        maxwet = dmax1(0d0,maxwet - exdelta)                             !remove water that freezes
        if((maxwet+xf)-sd > eps) maxwet = sd - xf                        !check for physical consistency

        waterheld = dmax1(0d0,waterheld - (omw - maxwet)*swi)            !still may be some water in snow

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case default

        error_code = 1
        error_type = 4

!       if(single_multi_flag == 0) write(10,'(''freq_id'',i10,'' Invalid
!     & mode in predictswi.f; day'',i4,'' hour'',i4)') freq_id,
!     &int(met(iw,ip_doy)),int(met(iw,ip_hr))

      end select

      end subroutine predictswi

! ******************************************************************************
      subroutine bottom(wavenumber,timestp,volume,kappa,n,sd,wave,bot)

!     Computes the time for the wave to reach the bottom
!     of the snow pack using newton's method of approaximation
!     if there is more than one wave

!     Variable               Remarks
!     wavenumber    number of smwaves in pack
!     bot(sec)      time for wave to reach bottom of snow pack
!     volume(cm)    size of wave (depth of water enetering with wave)
!     tau(sec)      time to bottom for the wave
!     dta(sec)      time between the present wave and the preceding wave
!     diff(cm)      volume(depth per unit area) of water removed from
!                   previous wave due to freezing
!     ticker        number of overestimates for tau in iteration scheme
!     bwn(s)        time for wave to reach bottom of snow pack
!     bwn1(cm)      depth of back edge of preceding wave
!     d(cm)         depth of wave at time tau
                 
!     bottom History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! calls the following subroutines:
!     depth (appended to this subroutine)

      IMPLICIT NONE

      integer(kind=4),intent(in):: wavenumber,timestp
      real(kind=8),intent(in):: volume,kappa,n,sd,wave(100,5)
      real(kind=8),intent(out):: bot

! local variables
      integer(kind=4):: asign,ticker,iter,wnt
      real(kind=8):: tau,d,d1,d2,d3,d4,p1,p2,p5,diff,bwn,bwn1,dta,f1,f2

! zero-out local variables
      asign = 0
      ticker = 0
      iter = 0
      wnt = 0
      wnt = wavenumber
      tau = 0d0
      d = 0d0
      d1 = 0d0
      d2 = 0d0
      d3 = 0d0
      d4 = 0d0
      p1 = 0d0
      p2 = 0d0
      p5 = 0d0
      diff = 0d0
      bwn = 0d0
      bwn1 = 0d0
      dta = 0d0
      bot = 0d0
      f1 = 0d0
      f2 = 0d0
      if(dabs(1d0-n) > eps) f1 = 1d0/(1d0 - n)
      if(n > 0d0) f2 = 1d0/n

      select case (wavenumber)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case(1)
!     Calculate the time for a single wave to reach the bottom

        if((dabs(n-1d0) > eps.and.volume > eps).and.(wave(wnt,2) > eps  &
                                                     .and.kappa > eps)) &
          bot = (((n - 1d0)/volume)**(n - 1d0))*((wave(wnt,2)/kappa)**n)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case default
!     More than one wave in the pack and the previous wave may have 
!     some of the water removed due to freezing

        asign = -1

        if(wave(wnt,2) > 0d0) then                                       !snow present
          if(dabs(wave(wnt-1,5)+1d0) <= eps) then                        !no frozen water
            tau = timestp*5d-1                                           !use half a time step
            dta = wave(wnt-1,1) - wave(wnt,1)

            iter = 0
            call depth(wnt,tau,kappa,n,sd,wave,d)
            do while(dabs(d-wave(wnt,2)) >= 5d-3*wave(wnt,2)            &
                                                       .and.iter < 100)
! use newton iteration scheme to calculate tau & depth
              d1 = d
              d2 = 0d0
              d3 = 0d0
              if((tau > eps.and.dabs(f1) > eps).and.                    &
                                               dabs(tau+dta) > eps) then 
                d2 = (f2 - 1d0)*(((tau**f1) - ((tau + dta)**f1))**asign)       
                d3 = f1*(tau**(n*f1) - (tau + dta)**(n*f1))
              end if
              d4 = d1*d2*d3             
              if(dabs(d4) > eps) tau = tau - (d - wave(wnt,2))/d4
              if(tau < 0d0.and.abs(ticker-1) > 0) then                   !tau can not be negative-try again
                tau = timestp*(4d0**(ticker - 1d0))
                ticker = ticker - 1
              end if
              call depth(wnt,tau,kappa,n,sd,wave,d)
              iter = iter + 1
            end do   !if (dabs(d-wave(wnt,2)) >= 
            bot = tau

            if(iter >= 100) then
              error_code = 1
              error_type = 4

!             if(single_multi_flag == 0) write(10,'(''freq_id'',i10,     &
!              &'' Too many iterations in bottom.f 1; day'',i4,          &
!              &'' hour'',i4)')freq_id,int(met(iw,ip_doy)),              &
!               int(met(iw,ip_hr))
            end if

            if((volume > eps.and.dabs(n-1d0) > eps).and.                &
                       (wave(wnt,2) > eps.and.dabs(kappa) > eps)) then
              if(bot > (((n - 1d0)/volume)**(n - 1d0))                  &
                                        *((wave(wnt,2)/kappa)**n)) then
                error_code = 1
                error_type = 4

!               if(single_multi_flag == 0) write(10,'(''freq_id'',i10,   &
!                &'' Error in time to bottom for the wave, 1; day'',i4,  &
!                &'' hour'',i4)') freq_id,int(met(iw,ip_doy)),           &
!               int(met(iw,ip_hr))
              end if   !if(bot > (((n - 1d0)/volume)**(n - 1d0))
            end if

!    some of water in previous wave is frozen.  Use same appraoch as above
!     but take the frozen liquid into consideration
          else if(dabs(wave(wnt-1,5)+1d0) > eps) then                    !yes frozen water
            if((volume > eps.and.dabs(n-1d0) > eps).and.                &
                             (wave(wnt,2) > eps.and.dabs(kappa) > eps)) &
              bot = (((n - 1d0)/volume)**(n - 1d0))*                    &
                                               ((wave(wnt,2)/kappa)**n)
            if(bot-wave(wnt,1) < wave(wnt-1,5)-wave(wnt-1,1)) then
              tau = timestp*5d-1
              dta = wave(wnt-1,1) - wave(wnt,1)
              if((wave(wnt-1,2) > eps.and.dabs(kappa) > eps).and.       &
                  dabs(f1) > eps) p1 = (wave(wnt-1,2)/kappa)**(n*f1)
              p2 = n - 1d0
              p5 = f1

!     water missing in previous wave due to freezing
              if(dabs(wave(wnt-1,5)) > eps.and.dabs(p5) > eps)          &
                diff = p1*p2*(wave(wnt-1,5)**p5)

              iter = 0
              if(((dabs(1d0-f2) > eps.and.f2 > eps).and.                &
                  (dabs(tau) > eps.and.dabs(tau+dta) > eps)).and.       &
                               (volume > eps.and.dabs(diff) > eps)) then
                bwn = kappa*((volume*(-f1))**(1d0 - f2))*(tau**f2)
                bwn1 = kappa*((diff*(-f1))**(1d0 - f2))*                &
                                                       ((tau + dta)**f2)
                if(bwn > bwn1) then                                      !overtake previous wave
                  if(wave(wnt,3) > diff) then
                    d1 = kappa*(((wave(wnt,3) - diff)*(-f1))**          &
                                                  (1d0 - f2))*(tau**f2)
                    d2 = (1d0 - ((tau/(tau + dta))**(-f1)))**(f2 - 1d0)
                    d = d1*d2
                  else
                    d = bwn
                  end if   !if(wave(wavenumber,3) > diff) then
                else
                  d = bwn
                end if   !if(bwn > bwn1) then
              end if

              do while(dabs(d-wave(wnt,2)) >= 5d-3*wave(wnt,2)          &
                                                       .and.iter < 100)
                d1 = d
                d2 = 0d0
                d3 = 0d0
                if((dabs(tau) > eps.and.dabs(tau+dta) > eps).and.       &
                                                   dabs(f1) > eps) then
                  d2 = (f2 - 1d0)*(((tau**f1) - ((tau + dta)**f1))      &
                                                                **asign)   
                  d3 = f1*(tau**(n*f1) - (tau + dta)**(n*f1))
                end if
                d4 = d1*d2*d3
                if(dabs(d4) > eps) tau = tau - (d - wave(wnt,2))/d4
                if(tau < 0d0.and.abs(ticker-1) > 0) then
                  tau = timestp*(2d0**float((ticker - 1)))
                  ticker = ticker - 1
                end if

                if(((dabs(1d0-f2) > eps.and.f2 > eps).and.              &
                  (dabs(tau) > eps.and.dabs(tau+dta) > eps)).and.       &
                               (volume > eps.and.dabs(diff) > eps)) then
                  bwn = kappa*((volume*(-f1))**(1d0 - f2))*(tau**f2)
                  bwn1 = kappa*((diff*(-f1))**(1d0 - f2))               &
                                                    *((tau + dta)**f2)

                  if(bwn > bwn1) then                                    !overtake previous wave
                    if(wave(wnt,3) > diff) then
                      d1 = kappa*(((wave(wnt,3) - diff)*(-f1))**        &
                                                  (1d0 - f2))*(tau**f2)
                      d2 = (1d0 - ((tau/(tau + dta))**(-f1)))**(f2 - 1d0)
                      d = d1*d2
                    else
                      d = bwn
                    end if   !if(wave(wnt,3) > diff) then
                  else
                    d = bwn
                  end if   !if(bwn > bwn1) then
                end if

                iter = iter + 1
              end do   !if(dabs(d-wave(wnt,2)) >= 
              bot = tau

              if(iter > 100) then
                error_code = 1
                error_type = 4

!                if(single_multi_flag == 0) write(10,'(/,''freq_id'',i10,&
!                 &'' Too many iterations in bottom.f 2;        day'',i4,&
!                 &'' hour'',i4)') freq_id,int(met(iw,ip_doy)),          &
!                  int(met(iw,ip_hr))            
              end if   !if(iter > 100) then

            if((volume > eps.and.dabs(n-1d0) > eps).and.                &
                             (wave(wnt,2) > eps.and.dabs(kappa) > eps)) then
                if(bot > (((n - 1d0)/volume)**(n - 1d0))                  &
                                        *((wave(wnt,2)/kappa)**n)) then
                  error_code = 1
                  error_type = 4

!                if(single_multi_flag == 0) write(10,'(''freq_id'',i10,  &
!                 &'' Error in time to bottom for the wave, 2; day'',i4, &
!                 &'' hour'',i4)') freq_id,int(met(iw,ip_doy)),          &
!                  int(met(iw,ip_hr))
                end if
              end if   
            end if   !if((bot-wave(wnt,1)) < (wave(wnt-1,5)
          end if   !if(wave(wnt-1,5) == -1d0) then
        else
          bot = 0d0
        end if   !if(wave(wnt,2) > 0d0) then

      end select

      end subroutine bottom

! ******************************************************************************
      subroutine collide(w1,w2,mode,timestp,smwaves,kappa,n,tottim,sd,  &
                         wave)
      
!     This routine combines two smwaves when wave n overtakes
!     wave n-1.  In mode 1 wave n overtakes wave n-1 in the 
!     snow pack, while in mode 2 it overtakes it at the bottom
!     of the pack

!     Computes the time for the wave to reach the bottom
!     of the snow pack using newton's method of approaximation
!     if there is more than one wave

!     Variable               Remarks
!     w1            wave 1 earlier wave
!     w2            wave2  later wave
!     diff(cm)      amount of refrozen water
!     delta(sec)    time difference between formation of smwaves
         
!     collide History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! calls the following subroutines:
!     removewave (appended to this subroutine)
!     bottom (appended to this subroutine)

      IMPLICIT NONE

      integer(kind=4),intent(in):: w1,w2,mode,timestp
      integer(kind=4),intent(inout):: smwaves
      real(kind=8),intent(in):: kappa,n,tottim,sd
      real(kind=8),intent(inout):: wave(100,5)

! local variables
      integer(kind=4):: wval
      real(kind=8):: delta,p1,p2,p5,diff,vol,bot,x1,vc1,f1

! zero-out local variables
      wval = 0
      delta = 0d0
      p1 = 0d0
      p2 = 0d0
      p5 = 0d0
      diff = 0d0
      vol = 0d0
      bot = 0d0
      x1 = 0d0
      vc1 = 0d0
      f1 = 0d0
      if(dabs(n-1d0) > eps) f1 = 1d0/(n - 1d0)

      delta = wave(w1,1) - wave(w2,1)                                    !dif in orgination time of smwaves

      select case (mode)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      case(1)                                                            !both smwaves collide in snow pack

        wave(w1,1) = wave(w2,1)
        wave(w1,2) = wave(w2,2)
      
        if(dabs(wave(w1,5)+1d0) <= eps) then                             !no freezing 
          wave(w1,3) = wave(w1,3) + wave(w2,3)                           !combine water for smwaves
        else
          if((wave(w1,2) > eps.and.dabs(kappa) > eps).and.              &
                        dabs(f1) > eps) p1 = (wave(w1,2)/kappa)**(n*f1)
          p2 = n - 1d0
          p5 = -f1
          if(wave(w1,5) > eps.and.dabs(p5) > eps) then
            diff = p1*p2*(wave(w1,5)**p5)
          else
            diff = 0d0
          end if
          wave(w1,3) = dmax1(0d0,wave(w1,3) + wave(w2,3) - diff)         !combine smwaves minus
        end if                                                           !frozen water
        wave(w1,5) = wave(w2,5)
        call removewave(w2,smwaves,wave)                                 !remove wave w2

        if(wave(w1,2) > 0d0.and.wave(w1,3) > eps) then
          vol = wave(w1,3)
          wval = w1
          call bottom(wval,timestp,vol,kappa,n,sd,wave,bot)
          wave(w1,4) = bot                                               !time for wave to get to bottom of pack
        else
          wave(w1,4) = 0d0
        end if

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      case(2)                                                            !smwaves collide at bottom of pack

        if((wave(w1,2) > eps.and.dabs(kappa) > eps).and.                &
                                                   dabs(f1) > eps) then
          x1 = (wave(w1,2)/kappa)**(n*f1)
          vc1 = x1*(n - 1d0)*((tottim + wave(w1,1))**(-f1))
          if(wave(w1,5) /= -1d0) then                                    !test for frozen water
            p1 = (wave(w1,2)/kappa)**(n*f1)
            p2 = n - 1d0
            p5 = -f1
            if(dabs(p5) > eps) diff = p1*p2*(wave(w1,5)**p5)
            vc1 = dmax1(0d0,vc1 - diff)                                  !remove frozen water from volume
          end if
        end if

        wave(w1,1) = wave(w2,1)
        wave(w1,2) = wave(w2,2)
        wave(w1,3) = vc1 + wave(w2,3)
        wave(w1,4) = tottim + wave(w1,1)
        wave(w1,5) = wave(w2,5)
        call removewave(w2,smwaves,wave)

      end select

      end subroutine collide

! ******************************************************************************
      subroutine removewave(wn,smwaves,wave)

!     Remove smwaves that have exit the pack or have been combined
!     with other smwaves
!
!     Variable               Remarks
!     wave( ,1)     time the wave orginated-used to get exit time
!     wave( ,2)     snow depth at the time of wave creation
!     wave( ,3)     current water volume of the wave
!     wave( ,4)     time it will take the wave to exit
!     wave( ,5)     time required to freeze the wave
!     wn            present wavenumber
!
!     removewave History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! no subroutines called

      IMPLICIT NONE

      integer(kind=4),intent(in):: wn
      integer(kind=4),intent(inout):: smwaves
      real(kind=8),intent(inout):: wave(100,5)

! local variables
      integer(kind=4):: ic,icc


      do ic=wn,smwaves
        do icc=1,5                                                       !reset wave parameters
          wave(ic,icc) = wave(ic+1,icc)
        end do
      end do

      smwaves = smwaves - 1

      end subroutine removewave      

! ******************************************************************************
      subroutine depth(wn,t,kappa,n,sd,wave,d)

!       Compute the depth of the wave at a specific time.  There
!     can be several scanerios: a single wave in the pack, multiple
!     smwaves in the pack, and partial freezing of one of the smwaves

!     Variable               Remarks
!     d(cm)         depth of wave
!     diff(cm)      volume of wave lost to freezing
!     wn            wavenumber
!     t(sec)        time to depth of wave
         
!     depth History
!     April 10 1996................Fortran coding completed
!     Nov 19 1999 .................eleminated a divide by zero problem

! no subroutines called

      IMPLICIT NONE

      integer(kind=4),intent(in):: wn
      real(kind=8),intent(in):: t,kappa,n,sd,wave(100,5)
      real(kind=8),intent(out):: d

! local variables
      integer(kind=4):: aa
      real(kind=8):: dta,d1,d2,c,p1,p2,p5,diff,f1,f2

! zero-out local variables
      aa = 0
      dta = 0d0
      d1 = 0d0
      d2 = 0d0
      c = 0d0
      p1 = 0d0
      p2 = 0d0
      p5 = 0d0
      diff = 0d0
      d = 0d0
      f1 = 0d0
      f2 = 0d0
      if(dabs(n-1d0) > eps) f1 = 1d0/(n - 1d0)
      if(n > eps) f2 = 1d0/n

! Compute the depth of the wave at a specific time
      select case (wn)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      case(1)                                                            !single wave in the pack

          if((dabs(t) > eps.and.dabs(1d0-f2) > eps).and.                &
                                              dabs(wave(1,3)*f1) > eps) &
            d = kappa*((wave(1,3)*f1)**(1d0 - f2))*(t**f2)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      case default    !more than one wave in pack

        aa = -1

        if(dabs(wave(wn-1,5)+1d0) <= eps) then                           !no freezing
          dta = wave(wn-1,1) - wave(wn,1)

          if((dabs(t) > eps.and.dabs(1d0-f2) > eps).and.                &
                                        dabs(wave(wn,3)*f1) > eps) then
            d1 = kappa*((wave(wn,3)*f1)**(1d0 - f2))*(t**f2)
            if(dta <= 0d0) then
              d = sd
            else
              d2 = (1d0 - ((t/(t + dta))**f1))**(f2 - 1d0)
              d = d1*d2
            end if
            c = kappa*((wave(wn,3)*f1)**(1d0 - f2))*(t**f2)
          end if

!     c is the depth a single wave can have. When there are multiple
!     smwaves the depth d must be greater than the value c
          if(d < (c*0.99d0)) then
            error_code = 1
            error_type = 4

!            if(single_multi_flag == 0) write(10,'(''freq_id'',i10,      &
!             &'' Error in wave depth in function depth.f; day'',i4,     &
!             &'' hour'',i4)') freq_id,int(met(iw,ip_doy)),              &
!              int(met(iw,ip_hr))
          end if
        else                                                             !freezing associated with previous wave
          if((dabs(kappa) > eps.and.wave(wn-1,2) > eps).and.            &
              dabs(f1) > eps) p1 = (wave(wn-1,2)/kappa)**(n*f1)
          p2 = n - 1d0
          p5 = f1
          if(dabs(p5) > eps) diff = p1*p2*(wave(wn-1,5)**p5)
          dta = wave(wn-1,1) - wave(wn,1)

          if((wave(wn,3) > eps.and.dabs(diff) > eps).and.               &
             (dabs(n-1d0) > eps.and.dabs((wave(wn,3)/diff)**            &
                                              (n-1d0)-1d0) > eps)) then
            if(t < (((wave(wn,3)/diff)**(n-1d0)-1d0)**aa)*dta.or.       &
                                                diff > wave(wn,3)) then
              if((dabs(t) > eps.and.dabs(1d0-f2) > eps).and.            &
                                             dabs(wave(wn,3)*f1) > eps) &
                d = kappa*((wave(wn,3)*f1)**(1d0 - f2))*(t**f2)
            else
              if((dabs(t) > eps.and.dabs(1d0-f2) > eps).and.            &
                                        dabs(wave(wn,3)*f1) > eps) then
                d1 = kappa*(((wave(wn,3) - diff)*f1)**(1d0 - f2))       &
                                                               *(t**f2)
                if(dabs(dta) <= eps) then
                  d = sd
                else
                  d2 = (1d0 - ((t/(t + dta))**f1))**(f2 - 1d0)
                  d = d1*d2
                end if
              end if
            end if
          end if
        end if

      end select

      end subroutine depth

end module module_snow

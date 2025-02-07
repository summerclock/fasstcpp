      subroutine new_profile(code5,first_time,pdens,oldsd,simelt,phie,  &
                             fsup,firup)

      use fasst_global

! this program calculates the temperature, moisture and ice profiles of the soil

! calls the following subroutines:
!     albedo_emis
!     surfenergy
!     flow_param
!     sub_divide
!     th_param
!     soil_tmp
!     soil_moisture
!     smooth
!     sflux

! uses the functions: dense,soilhumid,maxinfiltrate,head,vap_press

      implicit none

      integer(kind=4),intent(inout):: code5,first_time
      real(kind=8),intent(in):: pdens,oldsd,simelt
      real(kind=8),intent(inout):: phie,fsup,firup

! saved variables
      integer(kind=4) ntemp,iwm,iwt,iwmr,iwtr,iwsm
      real(kind=8):: infl_cum,tcum,klhtop,s2i,thetai,klhtopi,airo,totice
      real(kind=8):: maxrt,maxt,maxrm,maxm
      real(kind=8):: x1(6),x(14,maxn),stcalc(maxn)

      save:: ntemp,iwm,iwt,iwmr,iwtr,iwsm,infl_cum,maxrt,maxt,maxrm,maxm
      save:: tcum,klhtop,s2i,thetai,klhtopi,airo,x1,x,totice,stcalc

! local variables
      integer(kind=4):: i,ii,iter,ll,sn,icode,d0i,d1i,d2i,d3i,lflag
      integer(kind=4):: iflag,iceflag,full,isn
      real(kind=8):: errort,rhs_errort,errorm,rhs_errorm,smerror
      real(kind=8):: fcheat1,t2,swd,swu,ird,maxinfiltrate
      real(kind=8):: evaprate,infl_max,a,b,sigflo,sphmo,kmo,rpp,d1
      real(kind=8):: t1,vp,stempt,precip,taf,sheatf,sheatg,dnet1,dnet2
      real(kind=8):: lhtf,lheatg,pheatf,pheatg,rhsurf,cc1,in_rate
      real(kind=8):: dqdtf,dqdtg,sfac,sumrunoff,iru,pht1,sh,lh,ch,ch1
      real(kind=8):: mixrf,mixrg,mixra,mixraf,dense,rh,pres
      real(kind=8):: kave,temp2,sthick,fsdown,firdown,fvir
      real(kind=8):: fpheat,pdensnew,excess,eh,rain,soilhumid,rstep
      real(kind=8):: isf,isg,fcheat,mixrgs,snow,dedt,cheatw,qtopi
      real(kind=8):: cheati,cheatv,fheat,rhotot,fcheck,f1,f2,ftempo
      real(kind=8):: disg,disf,disfg,disgf,plimit,lh1
      real(kind=8):: hpondi,ph,mat,sert,summoist,rhow,rhoi,qtopv
      real(kind=8):: sumin,overland,head,mat1,sumthick,sumo,extra,maxsm
      real(kind=8):: albedoi,emisi,Z,r,sert1,mat2,mat3,mat4,toticeo,hatm
      real(kind=8):: qtop,qbot,simeltsm,vap_press,mixrgo,qv_max,D
      real(kind=8):: kvhtop,kvttop,tsign,temp(ntot),t0(4)
      real(kind=8):: shs,lhs,mixas,mixgs

      real(kind=8):: delvar(maxcol),delvar1(12),dmet(13),sinko(ntot)
      real(kind=8):: zt(ntot),delzt(ntot),iceo(ntot)
      real(kind=8):: wvco(ntot),told(ntot),fv1o(ntot),ph_old(ntot)
      real(kind=8):: sourceo(ntot),vino(ntot),flowuo(ntot),flowlo(ntot)
      real(kind=8):: rhov(ntot),rhoda(ntot),runoff(ntot),sm_old(ntot)
      real(kind=8):: grspheato(ntot),grthcondo(ntot)
      real(kind=8):: sinkro(ntot),y1(6),y(14,ntot)
      real(kind=8):: thvc(ntot),dthvdt(ntot),dthvdh(ntot)


      integer(kind=4),parameter:: maxiter = 10                           !maximum allowed iterations
      real(kind=8),parameter:: allowerrort = 2d-2                        !K, acceptable error for temperature iterations (5d-2)
      real(kind=8),parameter:: allowerrorm = 2d-3 !1d0                   !m, acceptable error for pressure head interations (3d-4)
!      real(kind=8),parameter:: maxerrorm = 3d0                           !unitless, maximum allowed error for water balance (3d-1)
!      real(kind=8),parameter:: maxerrort = 50d0                          !J/m^3, maximum allowed error for surface heat flux (1d6)
!      real(kind=8),parameter:: maxsmerror = 1d-3                         !m, maximum allowed error in total soil moisture (1d-3)

! zero-out variables
      ii = 0
      iter = 0
      sn = 0
      d0i = 0
      d1i = 0
      d2i = 0
      d3i = 0
      iflag = 0
      iceflag = 0
      full = 0
      lflag = 0
      isn = 0
      errorm = 0d0
      rhs_errorm = 0d0
      smerror = 0d0
      errort = 0d0
      rhs_errort = 0d0
      evaprate = 0d0
      infl_max = 0d0
      a = 0d0
      b = 0d0
      t1 = 0d0
      vp = 0d0
      stempt = 0d0
      precip = 0d0
      taf = 0d0
      simeltsm  = 0d0
      sheatf  = 0d0
      sheatg   = 0d0
      lhtf = 0d0
      lheatg  = 0d0
      pheatf  = 0d0
      pheatg  = 0d0
      rhsurf = 0d0
      cc1 = 0d0
      d1 = 0d0
      rpp = 0d0
      dqdtf = 0d0
      dqdtg = 0d0
      sfac = 0d0
      sumrunoff = 0d0
      mixrf = 0d0
      mixrg  = 0d0
      mixra = 0d0
      mixraf = 0d0
      kave = 0d0
      temp2 = 0d0
      sthick = 0d0
      fsdown = 0d0
      fsup = 0d0
      firdown = 0d0
      fvir = 0d0
      firup = 0d0
      fpheat = 0d0
      pdensnew = 0d0
      excess = 0d0
      isf = 0d0
      isg = 0d0
      fcheat = 0d0
      fcheat1 = 0d0
      cheati = 0d0
      cheatv = 0d0
      cheatw = 0d0
      fheat = 0d0
      rhotot = 0d0
      iru = 0d0
      pht1 =0d0
      sh = 0d0
      lh = 0d0
      ch = 0d0
      ch1 = 0d0
      eh = 0d0
      mixrgs = 0d0
      fcheck = 0d0
      ftempo = 0d0
      dedt = 0d0
      t2 = 0d0
      swd = 0d0
      swu = 0d0
      ird = 0d0
      dnet1 = 0d0
      dnet2 = 0d0
      disf = 0d0
      disg = 0d0
      disfg = 0d0
      disgf = 0d0
      in_rate = 0d0
      plimit = 0d0
      qtopi = 0d0
      hpondi = 0d0
      ph = 0d0
      sert = 0d0
      qtopv = 0d0
      sumthick = 0d0
      melt(iw) = 0d0
      albedoi = 0d0
      emisi = 0d0
      Z = 0d0
      r = 0d0
      sert1 = 0d0
      mat = 0d0
      mat1 = 0d0
      mat2 = 0d0
      mat3 = 0d0
      mat4 = 0d0
      toticeo = 0d0
      rstep = 0d0
      qv_max = 0d0
      D = 0d0
      kvhtop = 0d0
      kvttop = 0d0
      tsign = 0d0
      lh1 = 0d0
      shs = 0d0
      lhs = 0d0
      mixas = 0d0
      mixgs = 0d0

      do i=1,4
        t0(i) = 0d0
      end do

      do i=1,ntot
        delzt(i) = 0d0
        zt(i) = 0d0
        flowu(i) = 0d0
        flowl(i) = 0d0
        fv1(i) = 0d0
        told(i) = stt(i)
        ph_old(i) = phead(i)
        vino(i) = vin(i)
        rhov(i) = 0d0
        rhoda(i) = 0d0
        runoff(i) = 0d0
      end do

      do i=1,6
        y1(i) = 0d0
      end do

      do i=1,13
        do ii=1,ntot
          y(i,ii) = 0d0
        end do
      end do

      d1i = 1
      d2i = 2
      d3i = 3

! determine the new timestep
      deltat = deltati
      step = stepi

      lflag = 0
      if(iw /= 1) then
        if(aint(met(iw-1,ip_pt)) == 2.or.aint(met(iw,ip_pt)) == 2) lflag = 1
      else
        if(aint(met(iw,ip_pt)) == 2) lflag = 1
      end if
      if(hm > eps) lflag = 2

!      if(aint(met(iw,ip_pt)) == 2) then !.or.(simelt > eps.or.hpond > eps)) then
      if(lflag == 1) then
        deltat = deltati*0.25d-1
        rstep = timstep*3.6d3/deltat
        step = aint(rstep)
!      else if(totice > eps.and.hm < eps) then
!        deltat = deltati*1d-1
!        rstep = timstep*3.6d3/deltat
!        step = aint(rstep)
       else if(lflag == 2) then
        deltat = deltati*2.5d-1
        rstep = timstep*3.6d3/deltat
        step = aint(rstep)
      end if

      f1 = 1d0/float(step)
      f1 = anint(f1*1d20)*1d-20
      f2 = 1d0/deltat                                                    !1/s
      f2 = anint(f2*1d20)*1d-20

! maximum fluid infiltration
      albedoi = albedo
      emisi = emis

      if(dabs(nsoilp(nnodes,7)-spflag) > eps) then
        infl_max = 1d-2*nsoilp(nnodes,7)                                 !m/s
      else
        infl_max = 0d0
      end if
      overland = 0d0                                                     !m
      plimit = 1d-1*rough*dmax1(0d0,dmin1(1d0,dcos(sloper)))             !m (maximim ponding depth due to soil roughness)

! inflow due to snow and ice melting
      simeltsm = simelt*f1
      simeltsm = anint(simeltsm*1d15)*1d-15

! initialize variables if first time through
      if(first_time == 0) then
        ntemp = 0
        ntemp = nnodes
        icase = 0
        icaseo = 0
        hpond = 0d0
        infl_cum = 0d0
        tcum = 0d0
        klhtop = 0d0
        klhtop = nsoilp(nnodes,7)*1d-2                                   !m/s
        s2i = 0d0
        thetai = 0d0
        klhtopi = 0d0
        airo = 0d0
        totice = 0d0
        phie = 6.5d-1

        maxrt = -dabs(mflag)
        maxt = -dabs(mflag)
        maxrm = -dabs(mflag)
        maxm = -dabs(mflag)
        maxsm = -dabs(mflag)
        iwm = 0
        iwt = 0
        iwmr = 0
        iwtr = 0
        iwsm = 0

        if(infer_test == 1) then
          do i=1,ntot
            told(i) = 0d0
            wvco(i) = 0d0
            iceo(i) = 0d0
            ph_old(i) = 0d0
            sm_old(i) = 0d0
            sourceo(i) = 0d0
            sinko(i) = 0d0
            sinkro(i) = 0d0
            flowuo(i) = 0d0
            flowlo(i) = 0d0
            fv1o(i) = 0d0
            vino(i) = 0d0
!            too(i) = 0d0
!            smoo(i) = 0d0
!            woo(i) = 0d0
!            ioo(i) = 0d0
!            phoo(i) = 0d0
            stcalc(i) = 0d0

            told(i) = stt(i) !too(i)
            wvco(i) = wvc(i) !woo(i)
            iceo(i) = ice(i) !ioo(i)
            ph_old(i) = phead(i) !phoo(i)
            sm_old(i) = soil_moist(i) !smoo(i)
            sourceo(i) = source(i)
            sinko(i) = sink(i)
            sinkro(i) = sinkr(i)
            flowuo(i) = flowu(i)
            flowlo(i) = flowl(i)
            fv1o(i) = fv1(i)
            vino(i) = vin(i)
            stcalc(i) = stt(i)
          end do
          airo = met(iw-1,ip_tmp) + Tref

! calculate the number of nodes; surface configuration
          ntemp = nnodes
          icaseo = icase
          icase = 0                                                      !no vegetation, no snow
          sn = nnodes
          if(veg_flagl == 0.or.sigfl <= eps) then
            if(hm > eps) then
              ntemp = nnodes + 1
              icase = 1                                                  !no vegetation, yes snow
              sn = ntemp
            end if
          else
            if(hm <= eps) then
              ntemp = nnodes + 1
              icase = 2                                                  !yes vegetation, no snow
            else
              if(hfol_tot-hm <= eps) then
                ntemp = nnodes + 1
                icase = 4                                                !yes vegetation buried by snow
                sn = ntemp
              else if(hfol_tot-hm > eps) then
                ntemp = nnodes + 2
                icase = 3                                                !yes vegetation, yes snow
                sn = ntot - 1
              end if
            end if
          end if

          if(veg_flagl == 1) sn = nnodes
          stempt = stt(sn)
          kave = 5d-1*(grthcond(sn) + grthcond(sn - 1))
          sthick = zt(sn) - zt(sn-1)
          if(icase == 1.or.icase == 4) then
            kave = km*(1d0 - sigfl) + sigfl*kveg
            sthick = hm
          end if

          dmet(1) = dmet1(iw-1,1)                                        !total incoming solar
          dmet(2) = dmet1(iw-1,2)                                        !incoming IR
          dmet(3) = dmet1(iw-1,3)                                        !wind speed
          dmet(4) = dmet1(iw-1,4)                                        !air temperature
          dmet(5) = dmet1(iw-1,5)                                        !relative humidity
          dmet(6) = (dmet1(iw-1,6)/(timstep*3.6d3))*dcos(sloper)         !m/timestep (rain or snow)
          dmet(6) = anint(dmet(6)*1d15)*1d-15
          if(met(iw,ip_pt) == 2)                                        &
            runoff = (dmet1(iw-1,6)/(timstep*3.6d3))*dsin(sloper)
          dmet(7) = (dmet1(iw-1,7)/(timstep*3.6d3))*dcos(sloper)         !m/timestep (snow only)
          dmet(7) = anint(dmet(7)*1d15)*1d-15
          dmet(8) = dmet1(iw-1,8)                                        !reflected solar
          dmet(9) = dmet1(iw-1,9)                                        !direct incoming solar
          dmet(10) = dmet1(iw-1,10)                                      !diffuse incoming solar
          dmet(11) = met(iw-1,ip_ap)                                     !air pressure
          dmet(12) = met(iw-1,ip_irup)                                   !reflected & emitted IR
          dmet(13) = met(iw-1,ip_tsoil)                                  !measured soil surface temp

          ii = 0
          iter = 0
          call surfenergy(ii,iter,sn,kave,sthick,pdens,dmet,evaprate,   &
                          rhsurf,pheatf, pheatg,mixrf,mixraf,mixrg,     &
                          mixra,lhtf,lheatg,sheatf,sheatg,dqdtf,dqdtg,  &
                          sfac,d1, rpp,cc1,isg,isf,disg,disf,disfg,     &
                          disgf,lh1)
!          isurfoldg = isg
!          isurfoldf = isf
        else
! zero out the arrays
          do i=1,ntot
            wvco(i) = 0d0
            iceo(i) = 0d0
            ph_old(i) = 0d0
            sourceo(i) = 0d0
            sinko(i) = 0d0
            sinkro(i) = 0d0
            flowuo(i) = 0d0
            flowlo(i) = 0d0
            fv1o(i) = 0d0
            vino(i) = 0d0

            too(i) = 0d0
            woo(i) = 0d0
            ioo(i) = 0d0
            smoo(i) = 0d0
            phoo(i) = 0d0

            if(iw == istart) then
              wvco(i) = wvc(i)
              iceo(i) = ice(i)
              ph_old(i) = phead(i)
              sourceo(i) = source(i)
              sinko(i) = sink(i)                                         !m/s
              sinkro(i) = sinkr(i)                                       !unitless
              flowuo(i) = flowu(i)
              flowlo(i) = flowl(i)
              fv1o(i) = fv1(i)
              vino(i) = vin(i)
              too(i) = stt(i)
              smoo(i) = soil_moist(i)
              woo(i) = wvc(i)
              ioo(i) = ice(i)
              phoo(i) = phead(i)
              stcalc(i) = 0d0
              stcalc(i) = stt(i)
            end if
          end do

          if(iw == istart) then
            do i=1,6
              x1(i) = 0d0
            end do

            do i=1,13
              do ii = 1,maxn
                x(i,ii) = 0d0
              end do
            end do
          end if

          qtop = 0d0
          qbot = 0d0
          sphmo = 0d0
          kmo = 0d0
          sigflo = 0d0
        end if
        pdensnew = pdens
      end if

! initialize "old" variables
      sphmo = sphm
      kmo = km
      sigflo = sigfl
      ftempo = ftemp

      icode = 0
      do i=1,ntot
        flowuo(i) = flowu(i)
        flowlo(i) = flowl(i)
        fv1o(i) = fv1(i)
        if(ice(i) > 0d0) icode = 1
      end do

! determine the met interpolation variables
      do ll=1,maxcol
        delvar(ll) = 0d0
        if(iw /= 1) then
          if(aint(dabs(met(iw,ll)-mflag)*1d5)*1d-5 > eps.and.           &
                    aint(dabs(met(iw-1,ll)-mflag)*1d5)*1d-5 > eps) then
            delvar(ll) = (met(iw,ll) - met(iw-1,ll))*f1
            if(ll == 10.or.ll == 12)                                    &
              delvar(ll) = dmax1(0d0,met(iw,ll)/(timstep*3.6d3))
!              delvar(ll) = dmax1(0d0,met(iw,ll)*f2)
          else
            delvar(ll) = mflag
          end if
        end if
        delvar(ll) = anint(delvar(ll)*1d20)*1d-20
      end do

      do ll=1,12
        delvar1(ll) = 0d0
        dmet(ll) = 0d0
        if(iw /= 1) then
          if(aint(dabs(dmet1(iw,ll)-mflag)*1d5)*1d-5 > eps.and.         &
                  aint(dabs(dmet1(iw-1,ll)-mflag)*1d5)*1d-5 > eps) then
            delvar1(ll) = (dmet1(iw,ll) - dmet1(iw-1,ll))*f1
            if(ll == 6.or.ll == 7)                                      &
             delvar1(ll) = dmax1(0d0,dmet1(iw,ll)/(timstep*3.6d3))       !m/s
!             delvar1(ll) = dmax1(0d0,dmet1(iw,ll)*f2)
          else
            delvar1(ll) = mflag
          end if
        end if
        delvar(ll) = anint(delvar(ll)*1d20)*1d-20
      end do

      isn = iw - 1
      if(iw == 1) isn = iw
      dmet(1) = dmet1(isn,1)                                            !total incoming solar
      dmet(2) = dmet1(isn,2)                                            !incoming IR
      dmet(3) = dmet1(isn,3)                                            !wind speed
      dmet(4) = dmet1(isn,4)                                            !air temperature
      dmet(5) = dmet1(isn,5)                                            !relative humidity
      dmet(6) = delvar1(6)*dcos(sloper)                                 !m/timestep (rain or snow)
      dmet(6) = anint(dmet(6)*1d15)*1d-15
      if(met(iw,ip_pt) == 2) runoff = delvar1(6)*dsin(sloper)
      dmet(7) = delvar1(7)*dcos(sloper)                                 !m/timestep (snow only)
      dmet(7) = anint(dmet(7)*1d15)*1d-15
      dmet(8) = dmet1(isn,8)                                            !reflected solar
      dmet(9) = dmet1(isn,9)                                            !direct incoming solar
      dmet(10) = dmet1(isn,10)                                          !diffuse incoming solar
      dmet(11) = dmet1(isn,11)                                          !air pressure
      dmet(12) = dmet1(isn,12)                                          !reflected & emitted IR
      dmet(13) = dmet1(isn,13)                                          !measured soil surface temp

      stt(ntot) = dmet(4)

! adjust number of nodes for presence of snow and/or low veg
      ntemp = nnodes
      icaseo = icase
      icase = 0                                                          !no vegetation, no snow
      if(veg_flagl == 0.or.sigfl <= eps) then
        if(hm > eps) then
          ntemp = nnodes + 1
          icase = 1                                                      !no vegetation, yes snow
        end if
      else
        if(hm <= eps) then
          ntemp = nnodes + 1
          icase = 2                                                      !yes vegetation, no snow
        else
          if(hfol_tot-hm <= eps) then
            ntemp = nnodes + 1
            icase = 4                                                    !yes vegetation buried by snow
          else if(hfol_tot-hm > eps) then
            ntemp = nnodes + 2
            icase = 3                                                    !yes vegetation, yes snow
          end if
        end if
      end if

      do i=1,ntot
        told(i) = too(i)
        iceo(i) = ioo(i)
        wvco(i) = woo(i)
        ph_old(i) = phoo(i)
        if(node_type(i) == 'HM'.or.node_type(i) == 'SN')                &
          stt(i) = stcalc(i)
        if(dabs(told(i)) <= eps) told(i) = 1d0                           !safety catch

        if(i <= nnodes) then                                             !soil/pavement/rock layers
          sm_old(i) = smoo(i)
          grthcondo(i) = grthcond(i)
          grspheato(i) = grspheat(i)
          zt(i) = nz(i)
          sourceo(i) = source(i)                                         !m/s; place holder at this point
          sinko(i) = sink(i)                                             !m/s; total water loss from a node
          sinkro(i) = sinkr(i)                                           !unitless; root uptake

          if(dabs(stt(i)) <= eps) stt(i) = 1d0                           !safety catch to stop program bombing
!          pres = dmet(11) - 1d-2*(elev - zt(i) + phead(i))              &
!                                            *dense(stt(i),0d0,d1i)*grav  !mbar
          pres = dmet(11) - 1d-2*phead(i)*dense(stt(i),0d0,d1i)*grav     !mbar
          if(node_type(i) /= 'WA'.and.node_type(i) /= 'AI')then
            rh = soilhumid(i,phead(i),soil_moist(i),stt(i))
!            if(ice(i) > eps) rh = 1d0
          else if(node_type(i) == 'WA') then
            rh = 1d0
            pres = dmet(11)
          else if(node_type(i) == 'AI') then
            rh = 1d-2*dmet(5)
            pres = dmet(11)
          end if

          call sp_humid(d0i,pres,stt(i),rh,phead(i),t0(1),t0(2),t0(3),  &
                        t0(4),rhov(i),rhoda(i),t1,thvc(i),dthvdt(i),    &
                        dthvdh(i))
          if(ntype(i) == 27) wvc(i) = rhov(i)/rhoda(i)
        else                                                             !snow/lowveg/air layers
          pres = dmet(11)
          sourceo(i) = 0d0
          sinko(i) = 0d0
          sinkro(i) = 0d0

          rh = 1d-2*dmet(5)
          vp = vap_press(i,rh,pres)                                      !Pa (vapor pressure)

          phead(i) = -vp/(Rv*stt(i)*grav) + qtopv*deltat                 !m (head due to vapor pressure)
          if(hpond > eps) phead(i) = 0d0
          if(phead(i) > 0d0) phead(i) = 0d0
          phead(i) = anint(phead(i)*1d20)*1d-20
          if(i == ntot) hatm = phead(i)

          if(icase == 0) then                                            !no foliage, no snow/ice
            sn = nnodes
            ftemp = dmet(4)
            zt(i) = elev + 50d0

            km = 0d0
            ice(i) = 0d0
            wvc(i) = 0d0
            stt(i) = dmet(4)
            node_type(i) = 'AI'

          else if(icase == 1) then                                       !no foliage, yes snow/ice
            sn = ntemp
            ftemp = dmet(4)
            if(i == ntot) then
              zt(i) = elev + 50d0
              node_type(i) = 'AI'
            else if(i == ntot-1) then
              zt(i) = elev + hm
              nz(i) = zt(i)

              node_type(i) = 'HM'
              if(icaseo == 0) then                                       !bare ground the last time step
                told(i) = dmin1(dmet(4),Tref)
                too(i) = dmin1(dmet(4),Tref)
                stt(i) = dmin1(dmet(4),Tref)
                y(1,i) = stt(i)
                x(1,i) = stt(i)
                toptemp = stt(i)
                iceo(i) = 0d0
                wvco(i) = phie !6.5d-1
                kmo = km
              end if

              rh = 1d0
              vp = vap_press(i,rh,pres)
              rhoda(i) =(pres - vp)/(Rd*stt(i)) 
              rhov(i) = vp/(Rv*stt(i))

              wvc(i) = phie
              ice(i) = (refreeze + refreezei)*f1
              ice(i) = anint(ice(i)*1d20)*1d-20

              phead(i) = 0d0
            end if

          else if(icase == 4) then                                       !foliage buried by snow/ice
            sn = ntemp
            ftemp = stt(nnodes) + (hfol_tot/hm)*(stt(nnodes+1)          &
                                                         - stt(nnodes))
            ftemp = anint(ftemp*1d20)*1d-20

            if(i == ntot) then
              zt(i) = elev + 50d0
              node_type(i) = 'AI'
            else if(i == ntot-1) then
              zt(i) = elev + hm
              nz(i) = zt(i)

              node_type(i) = 'HM'

              if(icaseo == 0) then                                       !bare ground the last time step
                told(i) = dmin1(dmet(4),Tref)
                too(i) = dmin1(dmet(4),Tref)
                stt(i) = dmin1(dmet(4),Tref)
                toptemp = stt(i)
                iceo(i) = 0d0
                wvco(i) = phie !6.5d-1
                kmo = km
              end if
              if(kmo <= eps) kmo = km

              rh = 1d0
              vp = vap_press(i,rh,pres)
              rhoda(i) =(pres - vp)/(Rd*stt(i)) 
              rhov(i) = vp/(Rv*stt(i))

              wvc(i) = phie
              ice(i) = (refreeze + refreezei)*f1 !f1  !0d0
              ice(i) = anint(ice(i)*1d10)*1d-10

              phead(i) = 0d0
            end if
          else if(icase == 2) then                                       !yes foliage, no snow/ice
            sn = nnodes

            km = 0d0
            kmo = 0d0
            if(i == ntot) then
              zt(i) = elev + 50d0
              node_type(i) = 'AI'
            else if(i == ntemp) then
              zt(i) = elev + hfol_tot
              nz(i) = zt(i)
              node_type(i) = 'VG'

              wvc(i) = 0d0
              wvco(i) = 0d0
              ice(i) = 0d0
              iceo(i) = 0d0
              if((icaseo == 4.or.icaseo == 3).or.                       &
                           (first_time == 0.and.infer_test == 0)) then   !snow melted in last time step
                stt(i) = ftemp
                told(i) = ftemp
                too(i) = ftemp
              end if
            end if
          else if(icase == 3) then                                       !yes foliage, yes snow/ice
            sn = ntot - 1
            if(i == ntot) then
              zt(i) = elev + hfol_tot
              nz(i) = zt(i)

              stt(i) = ftemp
              if(icaseo == 4.or.icaseo == 2) too(i) = ftemp
              node_type(i) = 'VG'

              if(icaseo == 2) then
                told(i) = told(i-1)
                too(i) = too(i-1)
              end if
              wvc(i) = 0d0
              wvco(i) = 0d0
              ice(i) = 0d0
              iceo(i) = 0d0
            else if(i == ntot-1) then
              zt(i) = elev + hm
              nz(i) = zt(i)

              node_type(i) = 'HM'
              if(icaseo == 2) then
                told(i) = dmin1(dmet(4),Tref)
                too(i) = dmin1(dmet(4),Tref)
                stt(i) = dmin1(dmet(4),Tref)
                wvco(i) = phie !6.5d-1
                toptemp = stt(i)

                iceo(i) = 0d0
                if(kmo <= eps) kmo = km
              else
                told(i) = toptemp
                too(i) = toptemp
              end if

              rh = 1d0
              vp = vap_press(i,rh,pres)
              rhoda(i) =(pres - vp)/(Rd*stt(i)) 
              rhov(i) = vp/(Rv*stt(i))

              wvc(i) = phie
              ice(i) = (refreeze + refreezei)*f1 !f1  !0d0
              ice(i) = anint(ice(i)*1d10)*1d-10

              phead(i) = 0d0
            end if
          end if   !if(ntemp == nnodes+1) then
        end if   !if(i <= nnodes) then
        source(i) = 0d0
        sink(i) = 0d0
        sinkr(i) = 0d0
      end do

      if(hm > eps) then
        call th_param(d2i,ntot,pdens,pdensnew,rhov,rhoda,rhotot)
      end if

! determine node thicknesses
      do i=1,ntot
        if(elev < 1d3) then
          zt(i) = anint(zt(i)*1d20)*1d-20
        else
          zt(i) = anint(zt(i)*1d20)*1d-20
        end if

        if(i > ntemp) then
          delzt(i) = 0d0
        else if(i >= nnodes) then
          if(i == ntemp) then
            delzt(i) = (zt(i) - zt(i-1)) !5d-1*                          !m
            if(i > nnodes) delzt(i) = delzt(i)*5d-1
          else 
            delzt(i) = 5d-1*(zt(i) - zt(i-1)) + 5d-1*(zt(i+1) - zt(i))   !m
          end if
          delzt(i) = anint(delzt(i)*1d18)*1d-18
        else
          delzt(i) = anint(delzs(i)*1d18)*1d-18                          !m
        end if
      end do

! calculate the surface energy flux
      stempt = stt(sn)
      kave = 5d-1*(grthcond(sn) + grthcond(sn - 1))
      sthick = zt(sn) - zt(sn-1)
      if(icase == 1.or.icase == 4) then
        kave = km*(1d0 - sigfl) + sigfl*kveg
        sthick = hm
      end if

      if(hpond > eps) then
        if(met(iw,ip_zen) >= 9d1) then
          sgralbedo = 0d0
        else
          Z = met(iw,ip_zen)*pi/180d0                                    !radians
          r = dasin(dsin(Z)/1.33d0)
          albedo = 5d-1*(dsin(Z-r)*dsin(Z-r)/(dsin(Z+r)*dsin(Z+r))   &
                           + dtan(Z-r)*dtan(Z-r)/(dtan(Z+r)*dtan(Z+r)))
        end if
        emis = nsoilp(nnodes,4)
      else
        albedo = albedoi
        emis = emisi
      end if

      ii = 0
      iter = 0
      call surfenergy(ii,iter,sn,kave,sthick,pdens,dmet,evaprate,rhsurf,&
                      pheatf, pheatg,mixrf,mixraf,mixrg,mixra,lhtf,     &
                      lheatg,sheatf,sheatg,dqdtf,dqdtg,sfac,d1,rpp,cc1, &
                      isg,isf,disg,disf,disfg,disgf,lh1)

      rain = dmet(6)
      snow = dmet(7)

      if(iw.eq.istart.and.infer_test == 0) then
        isurfoldg = isg
        isurfoldf = isf
      end if

!maximum infiltration rate
      if(tcum <= eps.or.infl_cum <= eps) then
        thetai = soil_moist(nnodes)                                      !unitless
        s2i = maxinfiltrate(nnodes,thetai)                               !m^2/s
        klhtopi = klhtop                                                 !m/s
        in_rate = infl_max
      else if((tcum > eps.and.infl_cum > eps).and.s2i > eps) then
        t1 = 2d0*(infl_cum - klhtopi*tcum)*infl_max/s2i                  !unitless
        if(dabs(t1) < 5d1) then
          t1 = dexp(2d0*(infl_cum - klhtopi*tcum)*infl_max/s2i)          !unitless
          if(dcos(sloper)-1d0/t1 >= eps) then
            in_rate = infl_max*dmin1(1d0,dcos(sloper))                   !m/s
          else
            in_rate = infl_max*dmin1(1d0,dcos(sloper) - 1d0/t1 + 1d0)    !m/s
!            in_rate = infl_max*dmin1(1d0,dcos(sloper) + 1d0/(t1 - 1d0))  !m/s
          end if
          in_rate = dmax1(0d0,in_rate)
        end if
      end if

! calculate maximum allowed evaporation rate
      isn = nnodes
      call flow_param(isn,ntemp,qtop,qtopv,simeltsm,zt,delzt,qbot,      &
                      klhtop,kvhtop,kvttop)
      if(hm <= eps) then
        qv_max = kvhtop*dmax1(0d0,(hatm - phead(nnodes) + 1d0)          &
                                      + kvttop*(dmet(4) - stt(nnodes)))
      else
        if(phie > eps.and.stt(sn) > eps)                                &
          D = dmax1(0d0,(phie**(5d0/3d0))*(2.12d-5*(stt(sn)/Tref)**2d0))
        kvhtop = dmax1(0d0,D*(vap_press(sn,1d0,dmet(11))/(Rv*stt(sn)))  &
                             *grav/(dense(stt(sn),0d0,d1i)*Rv*stt(sn)))
        qv_max = kvhtop*dmax1(0d0,(hatm - phead(sn) + 1d0))
      end if

      tsign = 0d0
      if(evaprate /= 0d0) then
        tsign = dabs(evaprate)/evaprate
      else if(qv_max /= 0d0) then
        tsign = dabs(qv_max)/qv_max
      end if
      qtopv = 0d0
      qtopv = -tsign*dmin1(dabs(evaprate),dabs(qv_max))                  !m/s; evaprate < 0 =>evaporation
      if(dabs(evaprate) <= eps) qtopv = qv_max

! qtop > 0 : moisture moving into soil
      extra = 0d0
      hpondi = hpond
      if(hm <= eps.and.node_type(nnodes) /= 'WA') then
        precip = rain
        if(aint(met(iw,ip_pt)) /= 2) precip = 0d0                        !only rain infiltrates

        if(precip > eps.and.qtopv > eps) qtopv = 0d0

        if(hpond > eps) then
          hpond = hpondi + (qtopv + precip)*deltat + simeltsm            !m
          extra = qtopv
          qtopv = 0d0
          if(hpond < eps) then
            extra = 0d0
            qtopv = hpond*f2                                             !m/s
            hpond = 0d0
          end if

          qtopi = hpond*f2
          if(hpond < eps) hpond = 0d0
        else
          qtopi = precip + hpond*f2 + simeltsm*f2                        !m/s
        end if
      else
        qtopi = hpond*f2 + simeltsm*f2                                   !m/s
        if(hm > eps) then
          extra = qtopv
          qtopv = 0d0
        end if
      end if
      qtopi = anint(qtopi*1d20)*1d-20  !(1d0 - sigfl)*
      qtopv = anint(qtopv*1d20)*1d-20  !(1d0 - sigfl)*

!!      qtop = dmin1(in_rate,qtopi)                                        !m/s
!!      if(nsoilp(nnodes,24)-(soil_moist(nnodes)+ice(nnodes)) < eps) then
!!        qtop = 0d0
!!      else
        qtop = qtopi
!!      end if

      if(node_type(nnodes) /= 'WA') then
!!        if((qtop > eps.and.dabs(qtop-qtopi) > eps).or.hpond > eps) then !.or.                  &
!                                    dabs(qtop-hpond*f2) <= 1d-10) then
        if(hpond > eps.or.nsoilp(nnodes,24)-(soil_moist(nnodes)          &
                                              +ice(nnodes)) < eps) then
          qtop = dmax1(0d0,(nsoilp(nnodes,24) - (soil_moist(nnodes)     &
                                      + ice(nnodes)))*delzs(nnodes)*f2)  !m/s
!          hpond = dmax1(0d0,(qtopi - qtop)*deltat)                       !m
          hpond = dmax1(0d0,hpond - qtop*deltat)                         !m
          if(hpond > plimit) then
            overland = overland + (hpond - plimit)                       !m
            hpond = plimit
          end if
        end if
      else if(node_type(nnodes) == 'WA') then
        qtop = 0d0
        hpond = 0d0
        overland = 0d0
      end if
      hpond = anint(hpond*1d20)*1d-20

      if(hpond > eps) then
        if(met(iw,ip_zen) >= 9d1) then
          sgralbedo = 0d0
        else
          Z = met(iw,ip_zen)*pi/180d0                                    !radians
          r = dasin(dsin(Z)/1.33d0)
          albedo = 5d-1*(dsin(Z-r)*dsin(Z-r)/(dsin(Z+r)*dsin(Z+r))      &
                           + dtan(Z-r)*dtan(Z-r)/(dtan(Z+r)*dtan(Z+r)))
        end if
        albedo = anint(albedo*1d20)*1d-20
        emis = nsoilp(nnodes,4)
      else
        albedo = albedoi
        emis = emisi
      end if

      hatm = hatm + qtopv*deltat

      isn = 1
      call flow_param(isn,ntemp,qtop,qtopv,simeltsm,zt,delzt,qbot,      &
                      klhtop,kvhtop,kvttop)

      if(first_time == 0) first_time = 1

      if(iw == 1) then
        sumo = 0d0
        do i=1,ntot
          temp(i) = 0d0
          if(i <= nnodes) then
          rhow = dense(stt(i),0d0,d1i)
          if(ice(i) > eps) then
            rhoi = dense(stt(i),0d0,d2i)
            temp(i) = soil_moist(i) + ice(i)*rhoi/rhow
          else
            temp(i) = soil_moist(i) 
          end if
          if(wvc(i) > eps) temp(i) = temp(i) + wvc(i)*rhov(i)/rhow

          sumo = sumo + temp(i)*delzs(i)
          end if
        end do
      else
        sumo = tot_moist(iw-1)
      end if

! ******************************************************************************
! MAIN CALCULATION LOOP
      do ii=1,step

! partition the met variables to the temp time step
        call sub_divide(ii,delvar,delvar1,dmet)

        call albedo_emis(oldsd,dmet(1),dmet(8))

        iter = 0
        iflag = 0
        errorm = 1.1d0*allowerrorm
!        rhs_errorm = 1.1d0*maxerrorm
        errort = 1.1d0*allowerrort
!        rhs_errort = 1.1d0*maxerrort
!        smerror = 1.1d0*maxsmerror
        sert = rhs_errort
        sert1 = rhs_errorm
! ------------------------------------------------------------------------------
! START INNER ITERATION LOOP
!        do while(((errorm > allowerrorm.or.rhs_errorm > maxerrorm).or.  &
!                  smerror > maxsmerror).or.                             &
!                  (errort > allowerrort.or.rhs_errort > maxerrort))
        do while(errorm > allowerrorm.or.errort > allowerrort)
          lflag = 0
          if(iter > 0) then
            mat1 = mat
            mat4 = mat3
            mat3 = mat2
          end if
          mat = sert
          mat2 = sert1

          if(ii > 1) then !.and.icase == icaseo) then
            if(iter == 0) then
              do i=1,ntot
                if((y(3,i) > eps.and.x(3,i) < eps).or.(y(3,i) < eps.and.&
                                                    x(3,i) > eps)) then
                  told(i) = y(1,i)
                  wvco(i) = y(2,i)
                  iceo(i) = y(3,i)
                  sm_old(i) = y(4,i)
                  ph_old(i) = y(5,i)
                else
                  told(i) = x(1,i)
                  wvco(i) = x(2,i)
                  iceo(i) = x(3,i)
                  sm_old(i) = x(4,i)
                  ph_old(i) = x(5,i)
                end if

                sinko(i) = x(6,i)
                sinkro(i) = x(7,i)
                sourceo(i) = x(8,i)

                if(i <= nnodes) then
                  flowuo(i) = y(9,i)
                  flowlo(i) = y(14,i)
                  fv1o(i) = y(10,i)
                  vino(i) = y(11,i)
                end if

                grthcondo(i) = y(12,i)
                grspheato(i) = y(13,i)
              end do
              isurfoldg = x1(1)
              isurfoldf = x1(2)
              sigflo = x1(3)
              kmo = x1(4)
              sphmo = x1(5)
              mixrgo = x1(6)
            else
              do i=1,ntot
                told(i) = y(1,i)
                wvco(i) = y(2,i)
                iceo(i) = y(3,i)
                sm_old(i) = y(4,i)
                ph_old(i) = y(5,i)

                sinko(i) = y(6,i)
                sinkro(i) = y(7,i)
                sourceo(i) = y(8,i)

                if(i <= nnodes) then
                  flowuo(i) = y(9,i)
                  flowlo(i) = y(14,i)
                  fv1o(i) = y(10,i)
                  vino(i) = y(11,i)
                end if

                grthcondo(i) = y(12,i)
                grspheato(i) = y(13,i)
             end do
              isurfoldg = y1(1)
              isurfoldf = y1(2)
              sigflo = y1(3)
              kmo = y1(4)
              sphmo = y1(5)
              mixrgo = y1(6)
            end if
          end if

          full = 0
          do i=1,ntot
            x(1,i) = stt(i)
            x(2,i) = wvc(i)
            x(3,i) = ice(i)
            x(4,i) = soil_moist(i)
            x(5,i) = phead(i)
            x(6,i) = sink(i)
            x(7,i) = sinkr(i)
            x(8,i) = source(i)
            if(i <= nnodes) then
              x(9,i) = flowu(i)
              x(14,i) = flowl(i)
              x(10,i) = fv1(i)
              x(11,i) = vin(i)
            end if
            x(12,i) = grthcond(i)
            x(13,i) = grspheat(i)
            if(nsoilp(i,24)-soil_moist(i) < eps) full = full + 1
          end do
          x1(1) = isg
          x1(2) = isf
          x1(3) = sigfl
          x1(4) = km
          x1(5) = sphm
          x1(6) = mixrg

!          if(full < nnodes.or.(full >= nnodes.and.qtop+qtopv+hpond      &
!                                                          <= 0d0)) then
!            if(((errorm > allowerrorm.or.rhs_errorm > maxerrorm).or.    &
!                                            smerror > maxsmerror)) then
            if(errorm > allowerrorm) then
              call soil_moisture(sm_old,wvco,iceo,sourceo,sinkro,vino,  &
                                 told,thvc,dthvdh,runoff,rhs_errorm,    &
                                 errorm,smerror,sert1)
            else
              lflag = 2
            end if
!          else
!            lflag = 1
!            errorm = allowerrorm
!           rhs_errorm = maxerrorm
!           smerror = maxsmerror
!          end if

          Call soil_tmp(ntemp,iceflag,kmo,sphmo,isg,isf,disg,disf,      &
                        disfg,disgf,sigflo,rhotot,airo,dmet,grthcondo,  &
                        grspheato,zt,delzt,told,iceo,wvco,flowuo,flowlo,&
                        fv1o,sinko,sourceo,ph_old,sm_old,thvc,dthvdt,   &
                        rhov,rhoda,rhs_errort,errort,sert,ii,iter)

          if(iw == istart.and.(ii == 1.and.iter == 0)) then
            do i=1,nnodes
              if(ice(i) > eps.and.iceo(i) <= eps) then
                iceo(i) = ice(i)
                ioo(i) = ice(i)
                sm_old(i) = soil_moist(i)
                smoo(i) = soil_moist(i)
                wvco(i) = wvc(i)
                woo(i) = wvc(i)
                ph_old(i) = phead(i)
                phoo(i) = phead(i)
              end if
            end do
          end if

! calculate the surface energy flux
          stempt = toptemp

          call surfenergy(ii,iter,sn,kave,sthick,pdens,dmet,evaprate,   &
                          rhsurf,pheatf,pheatg,mixrf,mixraf,mixrg,mixra,&
                          lhtf,lheatg,sheatf,sheatg,dqdtf,dqdtg,sfac,   &
                          d1,rpp,cc1,isg,isf,disg,disf,disfg,disgf,lh1)

! get the hydraulic parameters
          isn = 1
          call flow_param(isn,ntemp,qtop,qtopv,simeltsm,zt,delzt,qbot,  &
                          klhtop,kvhtop,kvttop)

          toptemp = stt(nnodes)
          ftemp = dmet(4)

          toticeo = totice
          totice = 0d0
          do i=1,ntot
            totice = totice + ice(i)
            if(dabs(stt(i)) <= eps) stt(i) = 1d0                         !safety catch to stop program bombing
            if(i <= nnodes) then                                         !soil/pavement/rock layers
!              pres = dmet(11) - 1d-2*(elev - zt(i) + phead(i))          &
!                                            *dense(stt(i),0d0,d1i)*grav  !mbar
              pres = dmet(11) - 1d-2*phead(i)*dense(stt(i),0d0,d1i)*grav !mbar
              if(node_type(i) /= 'WA'.and.node_type(i) /= 'AI')then
                rh = soilhumid(i,phead(i),soil_moist(i),stt(i))
!                if(ice(i) > eps) rh = 1d0
              else if(node_type(i) == 'WA') then
                rh = 1d0
                pres = dmet(11)
              else if(node_type(i) == 'AI') then
                rh = 1d-2*dmet(5)
                pres = dmet(11)
              end if

              call sp_humid(d0i,pres,stt(i),rh,phead(i),t0(1),t0(2),    &
                            t0(3),t0(4),rhov(i),rhoda(i),t1,thvc(i),    &
                            dthvdt(i),dthvdh(i))

              if(ntype(i) == 27) wvc(i) = rhov(i)/rhoda(i)
            else                                                         !snow/lowveg/air/surface layers
              pres = dmet(11)
              rh = 1d-2*dmet(5)
              vp = vap_press(i,rh,pres)                                   !Pa (vapor pressure)

              phead(i) = -vp/(Rv*stt(i)*grav) + qtopv*deltat              !head due to vapor pressure (m)
              if(hpond > eps) phead(i) = 0d0
              if(phead(i) > 0d0) phead(i) = 0d0
              phead(i) = anint(phead(i)*1d20)*1d-20
              if(i == ntot) hatm = phead(i)

              if(icase == 0) then                                        !bare ground
                stt(i) = dmet(4)
                ice(i) = 0d0
                wvc(i) = 0d0
              else if(icase == 1) then                                   !no foliage, yes snow/ice
                stt(ntot) = dmet(4)
                toptemp = stt(ntemp)
                if(i == ntot-1) then
                  rh = 1d0
                  vp = vap_press(i,rh,pres)
                  rhoda(i) =(pres - vp)/(Rd*stt(i)) 
                  rhov(i) = vp/(Rv*stt(i))
                  phead(i) = 0d0
                end if
              else if(icase == 4) then                                   !foliage buried by snow/ice
                stt(ntot) = dmet(4)
                toptemp = stt(ntemp)
                ftemp = stt(nnodes) + (hfol_tot/hm)*(stt(nnodes+1)       &
                                                         - stt(nnodes))
                ftemp = anint(ftemp*1d20)*1d-20

                if(i == ntot-1) then
                  rh = 1d0
                  vp = vap_press(i,rh,pres)
                  rhoda(i) =(pres - vp)/(Rd*stt(i)) 
                  rhov(i) = vp/(Rv*stt(i))
                  phead(i) = 0d0
                end if
              else if(icase == 2) then                                   !yes foliage, no snow/ice
                stt(ntot) = dmet(4)
                ftemp = stt(ntemp)
              else if(icase == 3) then                                   !yes foliage, yes snow/ice
                toptemp = stt(ntot-1)
                ftemp = stt(ntot)
                if(i == ntot-1) then
                  rh = 1d0
                  vp = vap_press(i,rh,pres)
                  rhoda(i) =(pres - vp)/(Rd*stt(i)) 
                  rhov(i) = vp/(Rv*stt(i))
                  phead(i) = 0d0
                end if
              end if
            end if
          end do

! get the thermal parameters
          call th_param(d1i,ntot,pdens,pdensnew,rhov,rhoda,rhotot)
          if(hm > eps) then
            call th_param(d2i,ntot,pdens,pdensnew,rhov,rhoda,rhotot)
          end if  

          kave = 5d-1*(grthcond(sn) + grthcond(sn - 1))
          sthick = zt(sn) - zt(sn-1)
          if(icase == 1.or.icase == 4) then
            kave = km
            sthick = hm
          end if

          iter = iter + 1
          if(iter >= maxiter) then
            iflag = 1
            exit
          else if(iter >= 3) then
            if(dabs(mat+sert) < 1d0) then
              iflag = 2
              exit
            else if(dabs(mat1-sert) < 5d0.and.(mat1*mat < eps.and.      &
                                                  mat*sert < eps)) then
              iflag = 3
              exit
            else if(dabs(mat1+sert) < 5d0.and.(mat1*mat < eps.and.      &
                                                  mat*sert < eps)) then
              iflag = 4
              exit
            else if(dabs(mat-sert) < 1d0) then
              iflag = 5
              exit
            else if(dabs(sert) >= dabs(mat).and.((totice <= 0d0.and.    &
                                toticeo > eps).or.(totice > 0d0.and.    &
                                                 toticeo <= eps))) then
              iflag = 6

              do i=1,ntot
                stt(i) = x(1,i)
               wvc(i) = x(2,i)
               ice(i) = x(3,i)
                soil_moist(i) = x(4,i)
                phead(i) = x(5,i)
                sink(i) = x(6,i)
                sinkr(i) = x(7,i)
                source(i) = x(8,i)
                flowu(i) = x(9,i)
                flowl(i) = x(14,i)
                fv1(i) = x(10,i)
                vin(i) = x(11,i)
                grthcond(i) = x(12,i)
                grspheat(i) = x(13,i)
              end do
              isg = x1(1)
              isf = x1(2)
              sigfl = x1(3)
              km = x1(4)
              sphm = x1(5)
              mixrg = x1(6)

              exit
!            else if((dabs(mat3-sert1) < 3d0.and.dabs(mat4-mat2) < 3d0)  &
!                                     .and.(mat4*mat3*mat2*sert1 > eps   &
!                                          .and.sert1 > maxerrorm)) then
            else if((dabs(mat3-sert1) < 3d0.and.dabs(mat4-mat2) < 3d0)  &
                                  .and.mat4*mat3*mat2*sert1 > eps) then
              if(lflag == 0) then
                iflag = 7
                exit
              end if
            end if
          end if
        end do  !while.....
! END OF INNER ITERATION LOOP
! ------------------------------------------------------------------------------

! check freezing
!        do i=1,nnodes
!          if(ntype(i) < 20) then
!            if(soil_moist(i) > eps.and.stt(i)-Tref <= eps) then
!             if(stt(i) > Tref-5d-2.or.ice(i) <= eps) then
!               stt(i) = Tref
!              else
!                ice(i) = bftm(i)/dense(stt(i),0d0,d2i)
!                soil_moist(i) = 0d0
!                phead(i) = head(i,soil_moist(i))
!              end if
!            else if(ice(i) > eps.and.stt(i)-Tref > eps) then
!             if(stt(i) < Tref+5d-2) then
!               stt(i) = Tref
!             else
!               soil_moist(i) = bftm(i)/dense(stt(i),0d0,d1i)
!               phead(i) = head(i,soil_moist(i))
!               ice(i) = 0d0
!              end if
!           end if
!          end if
!        end do

!        do i=1,nnodes
!          if(ntype(i) < 20) then
!            if(soil_moist(i) > 1d-10.and.stt(i)-Tref < 0d0) then
!              if(ice(i) > eps) then
!                stt(i) = Tref
!                ice(i) = bftm(i)/dense(stt(i),0d0,d2i)
!                soil_moist(i) = 0d0
!                phead(i) = head(i,soil_moist(i))
!              end if
!            else if(ice(i) > 1d-10.and.stt(i)-Tref >= 0d0) then
!              stt(i) = Tref
!              soil_moist(i) = bftm(i)/dense(stt(i),0d0,d1i)
!              phead(i) = head(i,soil_moist(i))
!              ice(i) = 0d0
!            end if
!          end if
!        end do

!!        do i=nnodes+1,ntemp
!!          if((node_type(i) == 'HM'.and.hm > eps).and.stt(i) > Tref)     &
!!            stt(i) = Tref                                                !snow temp
!!        end do

!!        do i=1,nnodes
!!          if(node_type(i) == 'SN'.and.stt(i) > Tref) stt(i) = Tref
!!        end do

!        if(rhs_errort > maxrt) then
!          maxrt = rhs_errort
!          iwtr = iw
!        end if
!        if(rhs_errorm > maxrm) then
!          maxrm = rhs_errorm
!          iwmr = iw
!        end if
!!        if(errort > maxt) then
!!          maxt = errort
!!          iwt = iw
!!        end if
!!        if(errorm > maxm) then
!!          maxm = errorm
!!          iwm = iw
!!        end if
!        if(smerror > maxsm) then
!          maxsm = smerror
!          iwsm = iw
!        end if

!!        if(iflag >= 1) then
!          if(ii <= step.and.((errorm > allowerrorm.or.                  &
!                              rhs_errorm > maxerrorm).or.               &
!                             (errort > allowerrort.or.                  &
!                              rhs_errort > maxerrort))) then
!!          if(ii <= step.and.(errorm > allowerrorm.or.                    &
!!                                             errort > allowerrort)) then
!!            error_code = 1
!!            code5 = 1
!!            ecount = ecount + 1
!!            if(single_multi_flag == 0) then
!!              write(10,'(''freq_id'',i10,'' No convergence in new_'',   &
!!                       &''profile; day'',i4,'', hour'',i3,'', minute'', &
!!                       &i3)') freq_id,int(met(iw,ip_doy)),              &
!!                               int(met(iw,ip_hr)),int(met(iw,ip_min))

!!              write(10,'(''iw = '',i10,'', ii = '',i4,'', iter = '',i3,''&
!!                          &, iflag = '',i2)') iw,ii,iter,iflag

!!              if(errorm > allowerrorm)                                  &
!!                             write(10,*) 'errordelm',errorm,allowerrorm
!              if(rhs_errorm > maxerrorm)                                &
!                          write(10,*) 'rhs_errorm',rhs_errorm,maxerrorm
!              if(smerror > maxsmerror)                                  &
!                             write(10,*) 'colerrorm',smerror,maxsmerror
!!              if(errort > allowerrort)                                  &
!!                             write(10,*) 'errordelt',errort,allowerrort
!              if(rhs_errort > maxerrort)                                &
!                          write(10,*) 'rhs_errort',rhs_errort,maxerrort
!!              write(10,*)' '
!!            end if
!!          end if
!!        end if

! procede with the next time (deltat) step
! update parameters

! soil profile final check
        if(ii == step) then
!          call smooth(ntemp,ftempo)

          do i=nnodes+1,ntemp
            if(node_type(i) == 'HM'.and.hm > eps) then
              if(stt(i) > Tref) stt(i) = Tref                            !snow temp
              shs = sheatg
              lhs = lheatg
              mixas = mixraf
              mixgs = mixrg
            end if
          end do

          do i=1,nnodes
            if(node_type(i) == 'SN'.and.stt(i) > Tref) stt(i) = Tref
          end do

          toptemp = stt(sn)
          temp2 = stt(sn-1)
          stempt = toptemp

          call surfenergy(ii,iter,sn,kave,sthick,pdens,dmet,evaprate,   &
                          rhsurf,pheatf,pheatg,mixrf,mixraf,mixrg,mixra,&
                          lhtf,lheatg,sheatf,sheatg,dqdtf,dqdtg,sfac,   &
                          d1,rpp,cc1,isg,isf,disg,disf,disfg,disgf,     &
                          lhes(iw))

          if(node_type(nnodes) == 'WA') then
            hpond = 0d0
            overland = 0d0
          end if
          hpond = anint(hpond*1d20)*1d-20

          call th_param(d1i,ntot,pdens,pdensnew,rhov,rhoda,rhotot)
          if(hm > eps) then
            call th_param(d2i,ntot,pdens,pdensnew,rhov,rhoda,rhotot)
          end if

          kave = 5d-1*(grthcond(sn) + grthcond(sn - 1))
          sthick = zt(sn) - zt(sn-1)
          if(icase == 1.or.icase == 4) then
            kave = km
            sthick = hm
          end if

          isn = 1
          call flow_param(isn,ntemp,qtop,qtopv,simeltsm,zt,delzt,qbot,  &
                          klhtop,kvhtop,kvttop)
        end if !ii == step

        do i=1,ntot
          y(1,i) = stt(i)
          y(2,i) = wvc(i)
          y(3,i) = ice(i)
          y(4,i) = soil_moist(i)
          y(5,i) = phead(i)
          y(6,i) = sink(i)
          y(7,i) = sinkr(i)
          y(8,i) = source(i)
          y(9,i) = flowu(i)
          y(14,i) = flowl(i)
          y(10,i) = fv1(i)
          y(11,i) = vin(i)
          y(12,i) = grthcond(i)
          y(13,i) = grspheat(i)

!          told(i) = y(1,i)
!          wvco(i) = y(2,i)
!          iceo(i) = y(3,i)
!          sm_old(i) = y(4,i)
!          ph_old(i) = y(5,i)

!          sinko(i) = y(6,i)
!          sinkro(i) = y(7,i)
!          sourceo(i) = y(8,i)
!          sink(i) = 0d0
!          sinkr(i) = 0d0
!          source(i) = 0d0

!          if(i <= nnodes) then
!            flowuo(i) = y(9,i)
!            flowlo(i) = y(14,i)
!            fv1o(i) = y(10,i)
!            vino(i) = y(11,i)
!          end if

!          grthcondo(i) = y(12,i)
!          grspheato(i) = y(13,i)
        end do

        y1(1) = isg
        y1(2) = isf
        y1(3) = sigfl
        y1(4) = km
        y1(5) = sphm
        y1(6) = mixrg

!        isurfoldg = y1(1)
!        isurfoldf = y1(2)
!        sigflo = y1(3)
!       kmo = y1(4)
!        sphmo = y1(5)
      end do   !ii=1,step

! END OF MAIN LOOP
! ****************************************************************************

! final freezing check
      iflag = 0
      do i=1,ntemp
        stcalc(i) = stt(i)
        if((node_type(i) == 'HM'.and.hm > eps).and.stt(i) > Tref)       &
          stt(i) = Tref                                                  !snow temp
        if(node_type(i) == 'SN'.and.stt(i) > Tref) stt(i) = Tref

        if(ntype(i) < 20) then
          if(soil_moist(i) > 1d-7.and.stt(i)-Tref < 0d0) then
!            if(ice(i) > eps) then
              ice(i) = bftm(i)/dense(stt(i),0d0,d2i)
              soil_moist(i) = 0d0
              phead(i) = head(i,soil_moist(i))
!            end if
          else if(ice(i) > 1d-7.and.stt(i)-Tref >= 0d0) then
            soil_moist(i) = bftm(i)/dense(stt(i),0d0,d1i)
            phead(i) = head(i,soil_moist(i))
            ice(i) = 0d0
          end if

!          if(dabs(stt(i)-Tref) <= eps.and.(ice(i) > 1d-7.and.soil_moist(i) > 1d-7)) then
!            if(ice(i) >= soil_moist(i).and.soil_moist(i) > 1d-7) then
!              ice(i) = bftm(i)/dense(stt(i),0d0,d2i)
!              soil_moist(i) = 0d0
!              phead(i) = head(i,soil_moist(i))
!            else
!              soil_moist(i) = bftm(i)/dense(stt(i),0d0,d1i)
!              phead(i) = head(i,soil_moist(i))
!              ice(i) = 0d0
!            end if
!          else if(stt(i) < Tref.and.soil_moist(i) > 1d-7) then
!            ice(i) = bftm(i)/dense(stt(i),0d0,d2i)
!            soil_moist(i) = 0d0
!            phead(i) = head(i,soil_moist(i))
!          else if(stt(i) > Tref.and.ice(i) > 1d-7) then
!            soil_moist(i) = bftm(i)/dense(stt(i),0d0,d1i)
!            phead(i) = head(i,soil_moist(i))
!            ice(i) = 0d0
!          end if
        end if

        if((i <= nnodes.and.(ice(i) > eps.or.ioo(i) > eps)).and.        &
                       (node_type(i) /= 'WA'.and.node_type(i) /= 'AI')) &
          iflag = 1
      end do

      toptemp = stt(sn)
      temp2 = stt(sn-1)
      stempt = toptemp

      if(iflag == 1) then
        iceflag = 0
        do i=nnodes,1,-1
          if((ice(i) > eps.or.ioo(i) > eps).and.(node_type(i) /= 'WA'   &
                                       .and.node_type(i) /= 'AI')) then
            d1 = (dense(Tref,0d0,d1i) - dense(stt(i),0d0,d2i))          &
                                                 /dense(stt(i),0d0,d2i)
            if(i == nnodes) then
              delzs(i) = delzs(i) + d1*(ice(i) - ioo(i))*delzs(i)/      &
                                                            nsoilp(i,2)
              if(delzs(i) < delzsi(i)) delzs(i) = delzsi(i)
              if(ice(i) <= eps.and.delzs(i) > delzsi(i))                &
                                                   delzs(i) = delzsi(i)
              delzs(i) = anint(delzs(i)*1d20)*1d-20

              sumthick = delzs(i)
            else
              delzs(i) = delzs(i) + d1*(ice(i) - ioo(i))*delzs(i)/      &
                                                            nsoilp(i,2)
              if(delzs(i) < delzsi(i)) delzs(i) = delzsi(i)
              if(ice(i) <= eps.and.delzs(i) > delzsi(i))                &
                                                   delzs(i) = delzsi(i)
              delzs(i) = anint(delzs(i)*1d20)*1d-20

              nz(i) = sumthick + 5d-1*(nz(i+1) - nz(i))
              nz(i) = anint(nz(i)*1d20)*1d-20
              nz(i) = elev - nz(i)
              nz(i) = anint(nz(i)*1d20)*1d-20

              sumthick = sumthick + delzs(i)
              sumthick = anint(sumthick*1d20)*1d-20
            end if
          else
            if(i == nnodes) then
              sumthick = delzs(i)
            else
              nz(i) = sumthick + 5d-1*(nz(i+1) - nz(i))
              nz(i) = anint(nz(i)*1d20)*1d-20
              nz(i) = elev - nz(i)
              nz(i) = anint(nz(i)*1d20)*1d-20

              sumthick = sumthick + delzs(i)
              sumthick = anint(sumthick*1d20)*1d-20
            end if
          end if
          if(dabs(delzs(i)-delzsi(i)) > eps) iceflag = 1
        end do

        if(iceflag == 0) then
          do i=1,nnodes
            delzs(i) = delzsi(i)
            nz(i) = nzi(i)
          end do
        end if
      end if

      do i=1,ntot
        too(i) = told(i)
        woo(i) = wvco(i)
        ioo(i) = iceo(i)
        phoo(i) = ph_old(i)
        smoo(i) = sm_old(i)
      end do

      storll = stll
      storls = stls
      airo = dmet(4)

      if(newsd > eps) then
        pdensnew = pdens
      else
        pdensnew = 0d0
      end if

! compute the energy terms for output purposes
! if lowveg present, check for closer at veg surface; fcheck should be 0.0
      taf = (1d0 - 7d-1*sigfl)*dmet1(iw,4) + 6d-1*sigfl*ftemp           & 
                                                   + 1d-1*sigfl*toptemp

      if(veg_flagl /= 0.and.(sigfl > eps.and.hfol_tot > 0d0)) then
        t1 = 0d0
        t2 = 1d0
        t0(1) = 0d0
        t0(2) = 0d0
        lh = lhtf
        sh = sheatf
        ph = pheatf

        call sflux(d1i,step,sn,dmet1(iw,8),dmet1(iw,1),met(iw,ip_irup), &
                   dmet1(iw,2),t2,ftemp,cc1,taf,ph,sh,lh,mixraf,mixrf,  &
                   kveg,toptemp,hfol_tot,rpp,rhsurf,t1,t0(1),t0(2),     &
                   dqdtf,dqdtg,d1,fsdown,fsup,firdown,firup,pheatf,     &
                   sheatf,lheatf(iw),fcheat,fcheat1,fcheck,dnet1,dnet2)
      end if
!write(*,*)iw,fcheck,fsdown+fsup,firdown+firup
!write(*,*)iw,sheatf,lheatf(iw),fcheat+fcheat1

! check for closer; evap_heat should be 0.0 unless there is snow
      t2 = 1d0
      call sflux(d0i,step,sn,dmet1(iw,8),dmet1(iw,1),met(iw,ip_irup),   &
                   dmet1(iw,2),t2,stempt,cc1,taf,pheatg,sheatg,lheatg,  &
                   mixraf,mixrg,kave,temp2,sthick,rhsurf,rpp,ioo(sn),   &
                   woo(sn),smoo(sn),dqdtg,dqdtf,d1,sdown(iw),sup(iw),   &
                   irdown(iw),irup(iw),pheat1(iw),sheat(iw),lheat(iw),  &
                   cheat(iw),cheat1(iw),evap_heat(iw),dnet1,dnet2)

!sup(iw) = sup(iw) + fsup
      if(iw+1 < iend) then
        if(aint(dabs(dmet1(iw,8)-mflag)*1d5)*1d-5 <= eps.and.           &
                        aint(dabs(dmet1(iw+1,8)-mflag)*1d5)*1d-5 > eps) &
        dmet1(iw,8) = albedo*dmet1(iw,1)
      end if

! melt energy
      if(hm > eps.or.node_type(nnodes) == 'SN') then
        if(tmelt(iw) > Tref) then
          taf = (1d0 - 7d-1*sigfl)*dmet1(iw,4) + 6d-1*sigfl*ftemp       & 
                                                 + 1d-1*sigfl*tmelt(iw)

          call sflux(d0i,step,sn,dmet1(iw,8),dmet1(iw,1),               &
                     met(iw,ip_irup),dmet1(iw,2),sfac,tmelt(iw),cc1,    &
                     taf,pheatg,shs,lhs,mixas,mixgs,kave,temp2,sthick,  &
                     rhsurf,rpp,ioo(ntemp),woo(ntemp),smoo(ntemp),      &
                     dqdtg,dqdtf,d1,swd,swu,ird,iru,pht1,sh,lh,ch,ch1,  &
                     eh,dnet1,dnet2)
          melt(iw) = dmax1(0d0,eh)
        else
          melt(iw) = dmax1(0d0,evap_heat(iw))
        end if
      else
        melt(iw) = 0d0
      end if
      melt(iw) = melt(iw)*float(stepi)                                   !W/m^2
      melt(iw) = anint(melt(iw)*1d10)*1d-10

      overland = anint(overland*1d10)*1d-10
      hpond = anint(hpond*1d10)*1d-10
      ponding(iw) = hpond

      sumrunoff = 0d0                                                    !m/s
      do i=ntemp,1,-1
        sumrunoff = sumrunoff + runoff(i)
      end do
      sumrunoff = anint(sumrunoff*1d20)*1d-20

      if(precip > eps.and.hm <= eps) then
        infl_cum = infl_cum + qtop*deltat                                !m
        tcum = tcum + timstep*3.6d3                                      !s
      else if(qtop > eps.and.hm > eps) then
        infl_cum = infl_cum + qtop*deltat
        tcum = tcum + timstep*3.6d3
      else
        infl_cum = 0d0
        tcum = 0d0
      end if
      if(infl_cum < eps) infl_cum = 0d0

! moisture self check; debug information
!      if(iw.eq.iend) then
!        write(*,*) iwmr,maxrm,iwtr,maxrt,'rhs'
!        write(*,*) iwm,maxm,iwt,maxt,'del'
!        write(*,*) iwsm,maxsm,'sm'
!     end if

      summoist = 0d0
      sumin = 0d0
      do i=1,ntot
        temp(i) = 0d0
        if(i <= nnodes) then
        rhow = dense(stt(i),0d0,d1i)
        if(ice(i) > eps) then
          rhoi = dense(stt(i),0d0,d2i)
          temp(i) = soil_moist(i) + ice(i)*rhoi/rhow
        else
          temp(i) = soil_moist(i) 
        end if

        if(wvc(i) > eps) temp(i) = temp(i) + wvc(i)*rhov(i)/rhow

        summoist = summoist + temp(i)*delzs(i)                           !m
        sumin = sumin + vin(i)
        end if
      end do
      tot_moist(iw) = anint(summoist*1d-10)*1d10                         !m

      if((hpond > eps.or.hpondi > eps).or.overland > eps) then
        extra = extra + (dmax1(0d0,hpond - hpondi) + overland)*f2        !m/s
      end if
      extra = anint(extra*1d-10)*1d10

!if(iw.eq.1)then
!write(88,*)' '
!write(88,*)'mbalance'
!write(88,*)'m'
!end if
!write(88,*)(qtop+qtopv+qbot-sumrunoff)-(summoist-sumo)*f2*f1
!      write(99,101) met(iw,ip_doy)+met(iw,ip_hr)/24.0+                  &
!                    met(iw,ip_min)/(24.0*60.0),(temp(i),i=1,13),        &
!                    precip*1d3,(summoist-sumo)*f2*f1,                   &
!                    (qtop+qtopv+qbot-sumrunoff)+extra,-sumin+extra,     &
!                    extra,hm
! 101  format(f8.3,1x,13(f11.8,1x),f10.6,1x,4(e11.4,1x),f8.5)

!      write(88,102) met(iw,ip_doy)+met(iw,ip_hr)/24.0+                  &
!                    met(iw,ip_min)/(24.0*60.0),(stt(i),i=1,13),         &
!                    met(iw,ip_prec),stt(nnodes+1),                      &
!                    met(iw,ip_tmp)+Tref,ftemp,hm
! 102  format(f8.3,1x,13(f10.5,1x),4(f10.5,1x),f8.5)

      end subroutine new_profile

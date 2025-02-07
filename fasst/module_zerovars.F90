module module_zerovars

      use fasst_global

   contains

      subroutine zero_parameters(tflag,nt0,nm0,stflag,lcount,wtype,wvel,&
                                 wdepth,ttest,mtest,sname)

! This subroutine zeros out all arrays

!     no subroutines called

      implicit none

      integer(kind=4),intent(out):: tflag,nt0,nm0,stflag,lcount,wtype
      real(kind=8),intent(out):: wvel,wdepth
      character(len=1),intent(out):: ttest,mtest
      character(len=4),intent(out):: sname(maxl)

! local variables
      integer(kind=4) i,j


      freq_id = 1                                                        !complex id

! various flags
      tflag = 0
      nt0 = 0
      nm0 = 0
      stflag = 0
      lcount = 0
      ecount = 0
      rough = 0d0
      ttest = 'n'
      mtest = 'n'

      do i=1,maxl
        sname(i) = '    '
      end do

! initial soil conditions
      nlayers = 0
      nnodes = 0
      sgralbedo = 0d0
      sgremis = 0d0
      slope = 0d0
      aspect = 0d0
      albedo = 0d0
      emis = 0d0
      isurfoldg = 0d0
      isurfoldf = 0d0

! snow & ice parameters
      sn_istat(1) = 0
      sn_istat(2) = 0
      hsaccum = 0d0
      hi = 0d0
      dsnow = 0d0
      refreeze = 0d0
      newsd = 0d0
      atopf = 0d0
      km = 0d0
      sphm = 0d0
      hm = 0d0
      refreezei = 0d0
      iswe = 0d0
      iswe = spflag
      do i=1,13
        sn_stat(i) = 0d0
      end do

! vegetation parameters
      isigfl = spflag                !initial low vegetation density
      iepf = spflag                  !initial low veg emissivity
      ifola = spflag                 !initial low veg absorptivity
      ihfol = spflag                 !initial low veg height (cm)
      vegl_type = 0
      vegh_type = 0
      veg_flagl = 0                  !low vegetation flag (0 = none; 1 = yes)
      veg_flagh = 0                  !high vegetation flag (0 = no trees; 1 = yes)
      izh = 0d0                      !initial canopy height (m)
      isigfh = 0d0                   !initial canopy density
      hfol = 0d0                     !vegetation height (m)
      fola = 0d0                     !vegetation absorptivity
      epf = 0d0                      !vegetation emissivity
      albf = 1d0 - fola              !vegetation albedo
      lail = 0d0                     !vegetation leaf area index (m^2/m^2)
      sigfl = 0d0
      state = 0d0
      ftemp = spflag
      z0l = 0d0
      Zd = 0d0
      chnf = 0d0
      sqrt_chnf = 0d0
      trmlm = 0d0
      trmhm = 0d0
      rsl = 0d0
      storll = 0d0
      storls = 0d0
      sigfh = 0d0
      uaf = 0d0
      hfol_tot = 0d0

      water_flag = 0
      wtype = int(spflag)            !water type: 0 = rivers, 1 = lakes
      wdepth = spflag                !water depth (m)
      wvel = spflag                  !water velocity (m/s)

      do i=0,4
        do j=1,nclayers
          avect(i,j) = 0d0
        end do
      end do

! calculation variables
      step = 0
      stepi = 0
      error_code = 0
      error_type = 0
      vitd_index = 0
      deltat = 0d0
      deltati = 0d0
      hpond = 0d0
      toptemp = 0d0
      ptemp = 0d0

      do i=1,maxn
        ntype(i) = 0
        icourse(i) = 0
        zti(i) = 0d0
        tm(i) = 0d0
        zm(i) = 0d0
        sm(i) = 0d0
        soil_moist(i) = 0d0
        grthcond(i) = 0d0
        grspheat(i) = 0d0
        nz(i) = 0d0
        nzi(i) = 0d0
        pheadmin(i) = 0d0
        frl(i) = 0d0
        frh(i) = 0d0
        sinkr(i) = 0d0
        delzs(i) = 0d0
        bftm(i) = 0d0

        node_type(i) = ' '

        do j=1,maxp
          nsoilp(i,j) = 0d0
        end do
      end do

      do i=1,maxn+2
        stt(i) = 0d0
        ice(i) = 0d0
        phead(i) = 0d0
        wvc(i) = 0d0
        khl(i) = 0d0
        khu(i) = 0d0
        vin(i) = 0d0
        flowu(i) = 0d0
        flowl(i) = 0d0
        fv1(i) = 0d0
        source(i) = 0d0
        sink(i) = 0d0
      end do

      do i=1,nclayers
        ifoliage_type(i) = spflag
        ilai(i) = spflag
        iclump(i) = spflag
        irho(i) = spflag
        itau(i) = spflag
        ialp(i) = spflag
        ieps(i) = spflag
        idzveg(i) = spflag
      end do

      do i=1,maxl
        soiltype(i) = 0
        lthick(i) = 0d0
        rho_fac(i) = 0d0
        rho_fac(i) = 1d0
        stype(i) = '  '
        do j=1,maxp
          isoilp(i,j) = spflag
          soilp(i,j) = spflag
        end do
      end do

      do i=1,nclayers+nclayers
        stor(i) = 0d0
      end do

      do i=1,18
        do j=1,17
          veg_prp(i,j) = 0d0
        end do

        do j=1,maxn
          rk(i,j) = 0d0
        end do
      end do

      end subroutine zero_parameters

! ******************************************************************************
      subroutine zero_mxl_params

      implicit none

! local variables
      integer(kind=4):: i,j

! initialize all non-allocatable arrays depending on maxlines
      do i=1,maxlines
        sdens(i) = 0d0
        surfci(i) = 0d0
        surfrci(i) = 0d0
        surfcbr(i) = 0d0
        surfice(i) = 0d0
        surfmoist(i) = 0d0
        surficep(i) = 0d0
        surfmoistp(i) = 0d0
        surfemis(i) = 0d0
        surfemisf(i) = 0d0
        surfemisc(i) = 0d0
        surfd(i) = 0d0
        airt(1,i) = 0d0
        airt(2,i) = 0d0
        ft(i) = 0d0
        tt(i) = 0d0
        frthick(i) = 0d0
        twthick(i) = 0d0

        cheat(i) = 0d0
        cheat1(i) = 0d0
        sdown(i) = 0d0
        sup(i) = 0d0
        irdown(i) = 0d0
        irup(i) = 0d0
        pheat1(i) = 0d0
        lheat(i) = 0d0
        evap_heat(i) = 0d0
        sheat(i) = 0d0
        melt(i) = 0d0
        tmelt(i) = 0d0
        lhes(i) = 0d0
        lheatf(i) = 0d0
        cevap(i) = 0d0
        ponding(i) = 0d0
        tot_moist(i) = 0d0

        sstate(i) = ' '

        do j=1,maxcol
          met(i,j) = 0d0
        end do

        do j=1,10
          dmet1(i,j) = 0d0
        end do
      end do

      do i=1,nclayers
        do j=1,maxlines
          canopy_temp(i,j) = 0d0
        end do
      end do

      end subroutine zero_mxl_params

end module module_zerovars

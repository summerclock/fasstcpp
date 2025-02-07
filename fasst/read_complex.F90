      subroutine read_complex(wtype,io,lcount,nm0,wvel,wdepth,mtest,    &
                              sname)

      use fasst_global

      implicit none

! no subroutines called     

      integer(kind=4),intent(inout):: wtype,io
      integer(kind=4),intent(out):: lcount,nm0
      real(kind=8),intent(inout):: wvel,wdepth
      character(len=1),intent(out):: mtest
      character(len=4),intent(out):: sname(maxl)

! local variables
      integer(kind=4):: i,j,jj,k,kk,islope,iaspect,ielev,iwdepth,iwvel
      integer(kind=4):: ivegldens,iveghdens,iveghht,r_type,iswater
      integer(kind=4):: irough,stemp,iveglht,igwl,soilt(maxl)
      integer(kind=4):: wflag2,nlayers1,ithick(maxl)
      integer(kind=4):: soiltype1(3),soiltype2(maxn+3),soiltype3(maxn+3)
      integer(kind=4):: soiltype4(maxn+3)
      real(kind=8):: swater,d1,d2,lthick1(3)
      real(kind=8):: lthick2(maxn+3),lthick3(maxn+3),lthick4(maxn+3)
!      character(len=2):: vstype(30)
      character(len=4):: isname(maxl)

!      data vstype/ 'GW','GP','GM','GC','SW','SP','SM','SC','ML','CL',   &
!                   'OL','CH','MH','OH','PT','MC','CM','EV','  ','CO',   &
!                   'AS','  ','  ','  ','RO','WA','AI','  ','  ','SN'/


! initialize variables
      lcount = 0

      freq_id = 0          !terrain id
      vitd_index = 0       !met id
      islope = 0           !slope
      iaspect = 0          !aspect
      ielev = 0            !polygon elevation (m)
      wtype = 0            !0 = none, 1 = lake, 2 = river
      iwdepth = 0          !water depth (m)
      iwvel = 0            !water velocity (m/s)
      vegl_type = 0        !low vegetation type
      ivegldens = 0        !initial low veg density (%)
      iveglht = 0          !initial low veg height (cm)
      vegh_type = 0        !high vegetation type
      iveghdens = 0        !initial high veg density (%)
      iveghht = 0          !initial high veg height (m)
      r_type = 0           !road type
      iswater = 0          !surface water content or bottom layer if paved (vol/vol)*100
      irough = 0           !surface roughness
      stemp = 0            !surface temperature (C)
      igwl = 0             !ground water level (m)
      nlayers = 0          !number of ground layers
      soilt = 0            !soil layer type or bottom layer if paved
      ithick = 0           !soil layer thickness or bottom layer if paved
      wflag2 = 0           !used to place layer thicknesses for open water
      nlayers1 = 0         !used to place layer thicknesses for open water
      isname = '    '      !user surface soil name

      read(5,*,iostat=io) freq_id,vitd_index,islope,iaspect,ielev,wtype,&
                          iwdepth,iwvel,vegl_type,ivegldens,iveglht,    &
                          vegh_type,iveghdens,iveghht,r_type,iswater,   & 
                          irough,stemp,igwl,nlayers,                    &
                          (soilt(i),ithick(i),isname(i),i=1,nlayers)

! terrain orientation
      slope = 0d0                                                        !degrees from flat
      slope = dmax1(0d0,dmin1(9d1,float(islope)))
      sloper = dmax1(0d0,dmin1(1.57d0,slope*pi/1.8d2))
      sloper = anint(sloper*1d20)*1d-20

      aspect = 0d0                                                       !degrees from N; 180 = south
      aspect = dmax1(0d0,dmin1(3.6d2,float(iaspect)))

      elev = 0d0                                                         !m above/below sea level
      elev = float(ielev)

! open water
      if(wtype < 0.or.wtype == 998) wtype = 0                            !1 = lake/pond; 2 = river

      water_flag = 0
      if(wtype > 0) then
        water_flag = 1
        if(wtype == 2) water_flag = 2

        wdepth = float(iwdepth)                                          !m
        if(wdepth < 0d0) wdepth = spflag

        wvel = float(iwvel)                                              !m/s
        if(wvel < 0d0) wvel = spflag
      else
        water_flag = 0
      end if

      if(water_flag <= 2) then
        nlayers1 = 3
        soiltype1(1) = 26                                                !'WA'
        lthick1(1) = 2d-2
        soiltype1(2) = 26
        if(dabs(aint((wdepth-spflag)*1d5)*1d-5) > eps) then
          lthick1(2) = wdepth - lthick(1)
        else
          if(wtype == 1) then
            lthick1(2) = 9.98d0
          else
            lthick1(2) = 2.98d0
          end if
         end if
        soiltype1(3) = 7                                                 !'SM'
        lthick1(3) = 1d0
      end if

! low vegetation
      if(vegl_type <= 0.or.vegl_type == 998) vegl_type = 0               !see veg_prop for allowed values

      isigfl = spflag                                                    !% ground covered
      ihfol = spflag                                                     !cm
      if(vegl_type > 0) then
        veg_flagl = 1
        isigfl = float(ivegldens)*1d-2
        if(isigfl > 0.98d0) isigfl = 0.98d0
        if(isigfl < 0d0) isigfl = spflag

        ihfol = float(iveglht)
        if(ihfol < 0d0) ihfol = spflag
      end if

      if(vegl_type == 13.or.(wtype > 0.and.vegl_type /= 0)) then         !bog/marsh
        if(vegl_type == 13) then
          water_flag = 1
          if(nlayers1 == 0) then
            nlayers1 = 3
            soiltype1(1) = 26                                            !'WA'
            lthick1(1) = 2d-2
            soiltype1(2) = 26
            if(dabs(aint((wdepth-spflag)*1d5)*1d-5) > eps) then
              lthick1(2) = wdepth - lthick1(1)
            else
              lthick1(2) = 1.5d0 - lthick1(1)
            end if
            wdepth = lthick1(1) + lthick1(2)
            soiltype1(3) = 7                                             !'SM'
            lthick1(3) = 1d0
          end if
        else if(vegl_type == 14) then                                    !inland water
          water_flag = 1
          vegl_type = 0
          nlayers = 3
          soiltype1(1) = 26                                              !'WA'
          lthick1(1) = 2d-2
          soiltype1(2) = 26
          if(dabs(aint((wdepth-spflag)*1d5)*1d-5) > eps) then
            lthick1(2) = wdepth - lthick1(1)
          else
            lthick1(2) = 10d0 - lthick1(1)
          end if
          wdepth = lthick(1) + lthick1(2)
          soiltype1(3) = 7                                               !'SM'
          lthick1(3) = 1d0
        else if(vegl_type == 15) then
          water_flag = 3
          vegl_type = 0
        end if
      end if

! trees, i.e., high vegetation
      if(vegh_type <= 0.or.vegh_type == 998) vegh_type = 0               !see veg_prop for allowed values

      isigfh = spflag                                                    !% ground covered
      izh = spflag                                                       !m
      if(vegh_type > 0) then
        veg_flagh = 1
        isigfh = float(iveghdens)*1d-2
        if(isigfh > 0.98d0) isigfh = 0.98d0
        if(isigfh < 0d0) isigfh = spflag

        izh = float(iveghht)
        if(izh < 0d0) izh = spflag
      end if

      if(wtype > 0.and.vegh_type > 0) then                               !swamp
        water_flag = 1
        if(nlayers1 == 0) then
          nlayers1 = 3
          soiltype1(1) = 26
          lthick1(1) = 2d-2
          soiltype1(2) = 26
          if(dabs(aint((wdepth-spflag)*1d5)*1d-5) > eps) then
            lthick1(2) = wdepth - lthick1(1)
          else
            lthick1(2) = 1.5d0 - lthick1(1)
          end if
          wdepth = lthick1(1) + lthick1(2)
          soiltype1(3) = 7
          lthick1(3) = 1d0
        end if
      end if

! determine the initial surface moisture factor
      swater = 0d0
      if(iswater == 0.or.iswater == -998) then
        swater = spflag
      else
        swater = float(iswater)*1d-2                                     !vol/vol
      end if

      mtest = 'n'
      nm0 = 0
      if(swater > 0d0.and.dabs(swater-spflag) > eps) then
        mtest = 'y'
        nm0 = 1
        sm(1) = swater
        zm(1) = 0d0
      end if

! surface roughness
      rough = 0d0                                                        !m
      if(irough == 1) then                                               !smooth
        rough = 1d-3
      else if(irough == 2) then                                          !medium
        rough = 5d-3
      else if(irough == 3) then                                          !rough
        rough = 1d-2
      else if(irough == 4) then                                          !very rough
        rough = 5d-2
      end if

! surface temperature
      toptemp = mflag
      if(stemp > -101.and.stemp < 101) then
        toptemp = float(stemp) + Tref                                    !K
      else
        toptemp = mflag
      end if

! groundwater level
      gwl = -1d0
      if(igwl > 0) gwl = float(igwl)

! determine the soil profile - dirt, pavement, rock, permanent snow
      if(r_type == 998) r_type = 0

      if(nlayers > 0) then
        do i=1,nlayers
          lthick(i) = float(ithick(i))*1d-2                              !m
          if(lthick(i) <= 0d0) then
            if(nlayers == 1) then
              lthick(i) = 1d0
            else
              lthick(i) = 0.25d0
            end if
          end if

          if(soilt(i) > 0.and.soilt(i) /= 998) then
!            stype(i) = vstype(soilt(i))
            soiltype(i) = soilt(i)

            if(soiltype(i) == 26) then
              water_flag = 1
              wflag2 = 1
              wtype = 1
              wvel = 0d0
            end if

            if(r_type > 1) rho_fac(i) = 1.2d0
          else if(soilt(i) == -1) then
            sname(i) = isname(i)
            if(sname(i) == '    '.or.sname(i) == 'UNK_'.or.             &
               sname(i) == 'NA__') soiltype(i) = 0
          else if(soiltype(i) == -2) then
            if(r_type > 1.and.nlayers > 1) rho_fac(i) = 1.2d0
          else if(soiltype(i) == -3) then
            if(r_type > 1.and.nlayers > 1) rho_fac(i) = 1.2d0
          else if(soiltype(i) == -4) then
            rho_fac(i) = 0.8d0
          else if(soilt(i) == 0.or.soilt(i) == 998) then
            soiltype(i) = 0
            if(r_type > 1.and.nlayers > 1) rho_fac(i) = 1.2d0
          end if
        end do
        
        if(nlayers == 1.and.r_type >= 1) then
          if(r_type == 1) then                                           !dirt
            nlayers = 2
            soiltype(2) = soiltype(1)
            lthick(2) = 0.5d0*lthick(1)
            rho_fac(2) = 1d0

            lthick(1) = lthick(2)
          else if(r_type == 2) then                                      !asphalt
            nlayers = 3
            soiltype(3) = soiltype(1)
            lthick(3) = lthick(1)
            rho_fac(3) = 1d0

            soiltype(1) = 21                                             !'AS'
            lthick(1) = 0.08d0
            rho_fac(1) = 1d0

            soiltype(2) = 1                                              !'GW'
            lthick(2) = 0.6d0
            rho_fac(2) = 1.2d0
          else if(r_type == 3) then                                      !concrete
            nlayers = 3
            soiltype(3) = soiltype(1)
            lthick(3) = lthick(1)
            rho_fac(3) = 1d0

            soiltype(1) = 20                                             !'CO'
            lthick(1) = 0.08d0
            rho_fac(1) = 1d0

            soiltype(2) = 1                                              !'GW'
            lthick(2) = 0.6d0
            rho_fac(2) = 1.2d0
          else if(r_type == 4) then                                      !grass
            nlayers = 2
            soiltype(2) = soiltype(1)
            lthick(2) = 0.5d0*lthick(1)
            rho_fac(2) = 1d0

            lthick(1) = lthick(2)

            if(vegl_type == 0) vegl_type = 2                             !short grass
          end if
        end if
      else if(nlayers <= 0.and.water_flag == 0) then
        if(r_type == 1) then                                             !dirt
          nlayers = 2
          soiltype(1) = 7                                                !'SM'
          soiltype(2) = 7

          lthick(1) = 5d-1
          lthick(2) = 1d0

          rho_fac(1) = 1.2d0
          rho_fac(2) = 1d0    
        else if(r_type == 2) then                                        !asphalt
          nlayers = 3
          soiltype(1) = 21                                               !'AS'
          lthick(1) = 0.08d0
          rho_fac(1) = 1d0

          soiltype(2) = 1                                                !'GW'
          lthick(2) = 0.6d0
          rho_fac(2) = 1.2d0

          soiltype(3) = 7                                                !'SM'
          lthick(3) = 0.4d0
          rho_fac(3) = 1d0
        else if(r_type == 3) then                                        !concrete
          nlayers = 3
          soiltype(1) = 20                                                  !'CO'
          lthick(1) = 0.08d0
          rho_fac(1) = 1d0

          soiltype(2) = 7                                                !'GW'
          lthick(2) = 0.6d0
          rho_fac(2) = 1.2d0

          soiltype(3) = 7                                                !'SM'
          lthick(3) = 0.4d0
          rho_fac(3) = 1d0
        else if(r_type == 4) then                                        !grass
          if(vegl_type == 0) vegl_type = 2                               !short grass
          nlayers = 2

          soiltype(1) = 7                                                !'SM'
          soiltype(2) = 7

          lthick(1) = 5d-1
          lthick(2) = 1d0

          rho_fac(1) = 1.2d0
          rho_fac(2) = 1d0
        else                                                             !no road
          nlayers = 1
          lthick(1) = 1d0
          soiltype(1) = 7
          rho_fac(1) = 1d0
        end if

      else if(wflag2 == 1.and.nlayers <= 1) then
        nlayers = 3
        soiltype(2) = soiltype(1)
        lthick(2) = lthick(1) - 2d-2
        soiltype(1) = 26                                                 !'WA'
        lthick(1) = 2d-2
        soiltype(3) = 7
        lthick(3) = 1d0
      else if(nlayers1 > 0) then
        kk = nlayers1 + nlayers
        j = 1
        d1 = 0d0
        d2 = 0d0
        do i=1,kk
          if(i <= nlayers1.and.soiltype1(i) == 26) then
            lthick2(i) = d1 + lthick1(i)
            soiltype2(i) = soiltype1(i)
          else if(i > nlayers1.and.soiltype(i-nlayers) == 26) then
            lthick2(i) = d2 + lthick(i-nlayers1)
            soiltype2(i) = soiltype(i-nlayers1)
          end if
          if(lthick2(i) > eps) j = i
        end do

        call sort2(j,lthick2,soiltype2,k,lthick3,soiltype3)

        do i=1,k
          lthick3(i) = lthick3(i) - d1
          d1 = d1 + lthick3(i)
        end do
        d2 = d1

        do i=1,kk
          k = i - j
          if(i <= nlayers1.and.soiltype1(i) /= 26) then
            lthick2(k) = d1 + lthick1(i)
            soiltype2(k) = soiltype1(i)
          else if(i > nlayers1.and.soiltype(i-nlayers) /= 26) then
            lthick2(k) = d2 + lthick(i-nlayers1)
            soiltype2(k) = soiltype(i-nlayers1)
          end if
        end do
        call sort2(k,lthick2,soiltype2,jj,lthick4,soiltype4)

        do i=1,jj
          lthick4(i) = lthick4(i) - d1
          d1 = d1 + lthick3(i)
        end do

        nlayers = k + jj
        do i=1,nlayers
          if(i <= k) then
            lthick(i) = lthick3(i)
            soiltype(i) = soiltype3(i)
          else
            lthick(i) = lthick4(i-k)
            soiltype(i) = soiltype4(i-k)
          end if
        end do
      end if   !if(nlayers == 0.and.water_flag == 0)

      end subroutine read_complex
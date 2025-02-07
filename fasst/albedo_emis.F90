      subroutine albedo_emis(oldsd,sd,su)

      use fasst_global

! no subroutines called

      implicit none

      real(kind=8),intent(in):: oldsd,sd,su

! internal variables
      integer(kind=4):: aflag
      real(kind=8):: sgremis1,albedoo,salbedoo,salbedon,t1,albedo1,a,Z,r

      aflag = 0
      t1 = 0d0
      Z = 0d0
      r = 0d0
      albedo = 0d0
      salbedoo = 0d0
      salbedon = 0d0
      sgremis1 = 0d0
      sgremis1 = nsoilp(nnodes,4)
      albedoo = albedo
      if(iw == 1.and.hsaccum+newsd+hi > 0d0) albedoo = soalbedo

      albedo = 1d0
      if(sd > eps.and.aint(dabs(sd-mflag)*1d5)*1d-5 > eps) then
        if(su > eps.and.aint(dabs(su-mflag)*1d5)*1d-5 > eps) then
          albedo1 = (1d0 - sigfl)*su/sd
          if(albedo1 > 1d0) albedo1 = 0.99d0
          aflag = 1
        end if
      else if(sd < eps.and.aint(dabs(sd-mflag)*1d5)*1d-5 > eps) then
        albedo1 = 0d0
        aflag = 1
      end if

      if(ntype(nnodes) <= 19) then
!        sgralbedo = -0.246d0*dlog(soil_moist(nnodes) + 0.249d0)
        a = 0d0
        a = dexp(-1d0) + 0.23865d0*soil_moist(nnodes)/nsoilp(nnodes,9)
        sgralbedo = -nsoilp(nnodes,3)*dlog(a)
        sgralbedo = dmax1(5d-1*nsoilp(nnodes,3),                        &
                                     dmin1(nsoilp(nnodes,3),sgralbedo))

!        sgremis1 = 0.9d0 + 0.18d0*soil_moist(nnodes)                     !Chung and Horton (1987), Water R.R.23(12),p.2179
        a = 0d0
        a = 0.99d0 - nsoilp(nnodes,4)
        sgremis1 = nsoilp(nnodes,4) + a*soil_moist(nnodes)/             &
                                                       nsoilp(nnodes,9)
        sgremis1 = dmin1(0.99d0,dmax1(nsoilp(nnodes,4),sgremis1))
      else if(ntype(nnodes) == 26) then                                  !water
        Z = met(iw,ip_zen)*pi/180d0
        r = dasin(dsin(Z)/1.33d0)
        if(Z+r > 0d0) then
          sgralbedo = 5d1*(dsin(Z-r)*dsin(Z-r)/(dsin(Z+r)*dsin(Z+r))    &
                           + dtan(Z-r)*dtan(Z-r)/(dtan(Z+r)*dtan(Z+r)))
        else
          sgralbedo = nsoilp(nnodes,3)
        end if
        sgremis1 = nsoilp(nnodes,4)
      else
        sgralbedo = nsoilp(nnodes,3)
        sgremis1 = nsoilp(nnodes,4)
      end if

      if(dabs(hsaccum+hi+newsd) <= eps) then                             !no snow or ice
        emis = dmax1(sgremis1,sgremis)
        albedo = dmin1(albedo,sgralbedo)
      else if(hi > 0d0.and.dabs(hsaccum+newsd) <= eps) then              !ice
        emis = iemis
        albedo = ialbedo
      else if(hsaccum > oldsd.or.newsd > eps) then                       !new snow
        emis = semis
        albedo = dmin1(albedo,snalbedo)
      else                                                               !old snow
        emis = semis

! Douville et al (1995) Climate Dynamics 12(1) p.21-35
        if(dabs(atopf) <= eps) then                                      !no melting
          salbedoo = albedoo - 8d-3*timstep/2.4d1
        else                                                             !melting
          if(timstep > 5d3) then
            salbedoo = snalbedo
          else
            salbedoo = snalbedo + (albedoo - snalbedo)*                 &
                        dexp(-0.24d0*timstep/2.4d1)
          end if
        end if
        salbedoo = dmax1(soalbedo,dmin1(snalbedo,salbedoo))

! Roesch (2000)
        t1 = toptemp - Tref
        salbedon = 0.5d0 - 0.07582627d0*t1 - 5.5360168d-3*t1*t1         & 
                    - 5.2966269d-5*t1*t1*t1 + 4.2372742d-6*t1*t1*t1*t1
        salbedon = dmax1(soalbedo,dmin1(snalbedo,salbedon))

        albedo = dmin1(albedo,dmax1(salbedon,salbedoo))
      end if

      if(aflag == 1) albedo = albedo1

      albedo = anint(albedo*1d20)*1d-20
      emis = anint(emis*1d20)*1d-20

      end subroutine albedo_emis

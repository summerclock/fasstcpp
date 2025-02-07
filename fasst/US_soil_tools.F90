      subroutine get_user_soil_params(name,usclass,ustype,rhod,poros,   &
                                      alb,em,qtz,kdry,sathdrc,minwc,    &
                                      maxwc,vgbpres,vgexp,spheat,vgm,   &
                                      orgf,psand,psilt,pclay,pcarbon,   &
                                      pl,p200)

      use fasst_global

! ******************************************************************************
! This routine will get the soil parameters for new soil types. If the soil name
! from the usmc.inp is not one of the default values,
! this routine will search the appropriate file to find the required soil 
! information.
! NAME is the four alpha-numeric characters (alpha characters should be
! in caps).
! ******************************************************************************

! calls the following subroutines:
!     upr_case (appended to this subroutine)
!     map_USDA_SoilType_to_USCS

      implicit none

      character(len=4),intent(inout):: name
      character(len=4),intent(out):: usclass
      integer(kind=4),intent(out):: ustype
      real(kind=8),intent(out):: rhod,poros,alb,em,qtz,kdry,sathdrc
      real(kind=8),intent(out):: minwc,maxwc,vgbpres,vgexp,spheat,vgm
      real(kind=8),intent(out):: orgf,psand,psilt,pclay,pcarbon,pl,p200

! local variables
      integer(kind=4):: ir,ptest,nslines,io,d1,tflag,lis
      integer(kind=4):: map_USDA_SoilType_to_USCS
      real(kind=8):: rhoid,thetas,lnn,lnalpha
      character(len=4):: sname
      character(len=120):: dump

! zero-out variables
      ir = 0
      io = 0
      ptest = 0
      nslines = 0
      d1 = 0
      tflag = 0
      lis = 0
      rhoid = 0d0
      thetas = 0d0
      lnn = 0d0
      lnalpha  = 0d0
      sname = '    '

      ustype = 0
      rhod = spflag
      poros = spflag
      alb = spflag
      em = spflag
      qtz = spflag
      kdry = spflag
      sathdrc = spflag
      minwc = spflag
      maxwc = spflag
      vgbpres = spflag
      vgexp = spflag
      spheat = spflag
      vgm = spflag
      orgf = spflag
      psand = spflag
      psilt = spflag
      pclay = spflag
      pcarbon = spflag
      pl = spflag
      p200 = spflag
      usclass = '  '
      
      nslines = 23     !number of data lines + blank

! Make sure the name is in caps
      d1 = len_trim(name)
      call upr_case(d1,name)
      write(*,*)'User soil name  ',name

      do while(io /= -1)
        read(90,'(a)',iostat=io) sname

        d1 = len_trim(sname)
        call upr_case(d1,sname)
        if(name == sname) then
          read(90,*,iostat=io) ustype
          if(abs(ustype) < 100) then                                     !USCS classification
            usclass = 'USCS'
          else
            usclass = 'USDA'

            if(ustype > 0) then
              lis = ustype - 100
              ustype = map_USDA_SoilType_to_USCS(lis)
            else
              ustype = -(abs(ustype) - 100)
            end if
          end if

          if((ustype == 20.or.ustype == 21).or.(ustype == 25            &
                                           .or.ustype == 26)) tflag = 1

          read(90,*,iostat=io) rhod                                      !bulk dry density (g/cm^3)
          if(rhod <= 0d0.and.dabs(rhod-spflag) > eps) then
            write(*,*) name,' bulk dry density out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) rhoid                                     !intrinsic dry density (g/cm^3)
          if(rhoid <= 0d0.and.dabs(rhoid-spflag) > eps) then
            write(*,*) name,' intrinsic dry density out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) thetas                                    !solids fraction (0.0 - 1.0)
          if((thetas <= 0d0.or.thetas > 1d0).and.                       &
                                        dabs(thetas-spflag) > eps) then 
            write(*,*) name,' solids fraction out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) poros                                     !porosity (0.0 - 1.0)
          if((poros <= 0d0.or.poros > 1d0).and.                         &
                                         dabs(poros-spflag) > eps) then
            write(*,*) name,' porosity out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) alb                                       !albedo (0.0 - 1.0)
          if((alb <= 0d0.or.alb > 1d0).and.                             &
                                           dabs(alb-spflag) > eps) then 
            write(*,*) name,' albedo out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) em                                        !emissivity (0.0 - 1.0)
          if((em <= 0d0.or.em > 1d0).and.dabs(em-spflag) > eps) then 
            write(*,*) name,' emissivity out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) qtz                                       !quarzt content (0.0 - 1.0)
          if((qtz < 0d0.or.qtz > 1d0).and.dabs(qtz-spflag) > eps) then
            write(*,*) name,' quartz content out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) orgf                                      !organic fraction of solids (0.0 - 1.0)
          if((orgf < 0d0.or.orgf > 1d0).and.                            &
                                          dabs(orgf-spflag) > eps) then 
            write(*,*) name,' organic fraction out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) kdry                                      !dry thermal conductivity (W/m*K)
          if(kdry <= 0d0.and.dabs(kdry-spflag) > eps) then
            write(*,*) name,' dry thermal conductivity out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) spheat                                    !specific heat (J/kg*K)
          if((spheat <= 0d0.or.spheat > 2200d0).and.                    &
                                         dabs(spheat-spflag) > eps)then
            write(*,*) name,' specific heat out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) sathdrc                                   !saturated hydraulic conductivity (cm/s)
          if(sathdrc <= 0d0.and.dabs(sathdrc-spflag) > eps) then
            if(tflag /= 1) then
              write(*,*) name,' saturated hydraulic conductivity out ', & 
                              &'of range'
              ptest = 1
            end if
          end if

          read(90,*,iostat=io) minwc                                     !residual (minimum) water content (vol/vol)
          if((minwc <= 0d0.or.minwc > poros).and.                       &
                                         dabs(minwc-spflag) > eps) then
            if(tflag /= 1) then
              write(*,*) name,' min water content out of range'
              ptest = 1
            end if
          end if

          read(90,*,iostat=io) maxwc                                     !maximum water content (vol/vol)
          if((maxwc <= 0d0.or.maxwc > poros).and.                       &
                                         dabs(maxwc-spflag) > eps) then
            if(tflag /= 1) then 
              write(*,*) name,' max water content out of range'
              ptest = 1
            end if
          end if

          read(90,*,iostat=io) vgbpres                                   !van Genuchten bubbling pressure head (cm)
          if(vgbpres <= 0d0.and.dabs(vgbpres-spflag) > eps) then
            if(tflag /= 1) then
              write(*,*) name,' van Genuchten pressure head out of ',   &
                              &'range'
              ptest = 1
            end if
          end if
          if(dabs(vgbpres-spflag) > eps.and.vgbpres > eps)              &
                                                 vgbpres = 1d0/vgbpres   !1/cm

          read(90,*,iostat=io) vgexp                                     !van Genuchten exponent, n
          if(vgexp <= 0d0.and.dabs(vgexp-spflag) > eps) then
            if(tflag /= 1) then
              write(*,*) name,' van Genuchten exponent out of range'
              ptest = 1
            end if
          else if(dabs(vgexp-1d0) <= eps) then
            write(*,*) name,' van Genuchten exponent can not equal 1'    !causes a divide by 0 error
            ptest = 1
          end if    
          if(dabs(vgexp-spflag) > eps.and.vgexp > eps)                  &
                                                  vgm = 1d0 - 1d0/vgexp  !unitless

          read(90,*,iostat=io) psand                                     !percent sand
          if((psand < 0d0.or.psand > 1d0).and.                          &
                                         dabs(psand-spflag) > eps) then
            write(*,*) name,' percent sand out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) psilt                                     !percent silt
          if((psilt < 0d0.or.psilt > 1d0).and.                          &
                                         dabs(psilt-spflag) > eps) then
            write(*,*) name,' percent silt out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) pclay                                     !percent clay
          if((pclay < 0d0.or.pclay > 1d0).and.                          &
                                         dabs(pclay-spflag) > eps) then
            write(*,*) name,' percent clay out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) pcarbon                                   !carbon content
          if((pcarbon < 0d0.or.pcarbon > 1d0).and.                      &
                                        dabs(pcarbon-spflag) > eps)then
            write(*,*) name,' percent carbon out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) pl                                        !plactic limit
            if((pl < 0d0.or.pl > 2d2).and.dabs(pl-spflag) > eps) then
            write(*,*) name,' plastic limit out of range'
            ptest = 1
          end if

          read(90,*,iostat=io) p200                                      !percent fines passing #200 sieve
          if((p200 < 0d0.or.p200 > 1d2).and.                            &
                                          dabs(p200-spflag) > eps) then
            write(*,*) name,' percent passing #200 sieve out of range'
            ptest = 1
          end if

! See if any parameters can be calculated from others
          if(dabs(rhod-spflag) <= eps.and.(dabs(rhoid-spflag) > eps     &
             .and.dabs(thetas-spflag) > eps)) rhod = rhoid*thetas

          if(dabs(rhoid-spflag) <= eps.and.(dabs(rhod-spflag) > eps     &
             .and.dabs(thetas-spflag) > eps)) rhoid = rhod/thetas

          if(dabs(thetas-spflag) <= eps.and.(dabs(rhod-spflag) > eps    &
             .and.dabs(rhoid-spflag) > eps)) thetas = rhod/rhoid

          if(dabs(poros-spflag) <= eps.and.                             &
                                       dabs(thetas-spflag) > eps) then
            poros = 1d0 - thetas
            if(poros <= 0d0.or.poros > 1d0) then
              write(*,*) name,' porosity out of range'
              ptest = 1
            end if
          end if

          if(dabs(rhod-spflag) <= eps.and.(dabs(rhoid-spflag) > eps     &
             .and.dabs(poros-spflag) > eps)) rhod = rhoid*(1d0 - poros)

          if(dabs(poros-spflag) > eps.and.dabs(maxwc-spflag) <= eps)    &
            maxwc = poros

          if(dabs(orgf-spflag) <= eps.and.                              &
                                       dabs(pcarbon-spflag) > eps) then  
            orgf = 2d0*pcarbon                                           !SSSA reccomendation
            if(orgf > 1d0) orgf = 1d0
          end if

          if(dabs(qtz-spflag) <= eps.and.dabs(psand-spflag) > eps)      &
            qtz = psand                                                  !Peters-Lidard et al., J. Atmos. Sci., V55 (1998)

! calculate max water content, min water content, van Genuchten parameters using
! PTF of Vereecken et al. (1989)"Estimating the soil moisture retention characteristic from
! texture, bulk density and carbon content", Soil Science, V.148, No.6, p.389-403
          if(tflag /= 1.and.ustype /= 30) then
            if(dabs(psand-spflag) > eps.and.dabs(pclay-spflag) > eps    &
               .and.dabs(pcarbon-spflag) > eps) then
              if(dabs(psilt-spflag) <= eps)                             &
                psilt = 1d0 - (psand + pclay)

              if(dabs(minwc-spflag) <= eps)                             &
                minwc = 1.5d-2 + 5d-3*pclay + 1.4d-2*pcarbon

              if(dabs(maxwc-spflag) <= eps)                             &
               maxwc = 0.81d0 - 0.283d0*rhod + 1d-3*pclay

              if(maxwc > poros) maxwc = poros

              if(dabs(vgbpres-spflag) <= eps) then
                lnalpha = -2.486d0 + 2.5d-2*psand                       &
                           - 0.351d0*pcarbon - 2.617d0*rhod             &
                           - 2.3d-2*pclay
                if(dabs(lnalpha) < 5d1) vgbpres = dexp(lnalpha)          !1/cm
              end if

              if(dabs(vgexp-spflag) <= eps) then
                lnn = 5.3d-2 - 9d-3*psand - 1.3d-2*pclay                &
                       + 1.5d-4*psand*psand
                if(dabs(lnn) < 5d1.and.dabs(lnn) > eps) then
                  vgexp = dexp(lnn)
                  vgm = 1d0 - 1d0/vgexp
                end if
              end if
            end if
          end if

          if(ptest == 1) then
             write(*,*) ' Soil Input Parameters out of range, STOPPING'
             stop
          end if
          return
        else
          do ir=1,nslines                                                !skip to the next soil type
            read(90,'(a)',iostat=io) dump
          end do
        end  if
      end do

! Unable to find the user defined soil type silty-sand will be used as a default
      write(*,'('' !! Can not find your soil, defaulting to sand !! '') &
            &')                                                          !Graves silty-sand, MA
      usclass = 'USCS'
      ustype = 7                                                         !silty sand
      rhod = 1.49d0                                                      !Guyman et al. (1993), p.55
      rhoid = 2.73d0                                                     !Guyman et al. (1993), p.55
      thetas = rhod/rhoid
      poros = 1d0 - thetas
      alb = 0.35d0                                                       !Sullivan et al. (1997), p.38
      em = 0.92d0                                                        !Sullivan et al. (1997), p.51
      qtz = 0.8d0                                                        !Tarnawski (1997), p.96
      orgf = 0.0d0
      kdry = 0.831d0   !0.281d0                                          !Jordan (2000)
      spheat = 830d0
      sathdrc = 0.000533d0                                               !Guyman et al. (1993), p.55
      minwc = 0.001d0
      maxwc = poros
      vgbpres = 1d0/23.5474d0                                            !Jordan (2000)
      vgexp = 1.50d0                                                     !Jordan (2000)
      vgm = 1d0 - 1d0/vgexp                                              !unitless
      psand = spflag                                                     !Loamy Sand, Cospy et al. (1984), p.683
      psilt = spflag                                                     !Loamy Sand, Cospy et al. (1984), p.683
      pclay = spflag                                                     !Loamy Sand, Cospy et al. (1984), p.683
      pcarbon = spflag                                                   !Vereecken et al. (1989), p.391
      pl = spflag
      p200 = spflag

      end subroutine get_user_soil_params

! ******************************************************************************
! ******************************************************************************
      SUBROUTINE upr_case(d1,AL)
                                                     
!     --------------------------------------------------------------------------
!     THIS ROUTINE TAKES AN INCOMING CHARACTER STRING VARIABLE 'AL' AND
!     CONVERTS ALL ELEMENTS TO UPPER CASE. ONLY ELEMENTS CORRESPONDING
                                                                       
!     FOR EXAMPLE, IF AL IS DEFINED BY                                  
!     DATA AL/'PATH'/                                                   
                                                                       
!     THEN AL IS RETURNED AS PATH                                                         
!     --------------------------------------------------------------------------

      implicit none

      integer(kind=4),intent(in):: d1
      character(len=d1),intent(inout):: AL

! loacal variables
      integer(kind=4):: IAL,I,INEW,ILZ                                         
      CHARACTER(len=1):: CAL(d1)
      CHARACTER(len=d1):: DUMMY 
                                                   
      DATA IAL/97/ILZ/122/ 
      
                                                   
      DUMMY=AL
      DO I=1,d1
      CAL(I) = DUMMY(I:I)                                                       
        INEW=ICHAR(CAL(I))                                                
        IF(INEW >= IAL) THEN
          IF(INEW <= ILZ) then
            CAL(I)=CHAR(INEW-32)
            DUMMY(I:I) = CAL(I)
          ENDIF
        ENDIF
      END DO
                                                          
      AL=DUMMY
                                                           
      END SUBROUTINE upr_case

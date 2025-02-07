      subroutine icethick(hinit,hfinal)

      use fasst_global

! this subroutine calculates ice accretion/depletion due to freezing rain
! if snow is present, then no ice accumulates

! no subroutines called

      implicit none

      real(kind=8),intent(in):: hinit
      real(kind=8),intent(out):: hfinal

! local variables
      real(kind=8):: hnew1,qtop,vimelt1,freeze_frac
      real(kind=8):: qb_f1,qbot,hnew2,vimelt2,lheat1

      real(kind=8),parameter:: pdens = 1d3                               !kg/m^3 (precip density)


! initialize variables
      hnew1 = 0d0
      qtop = 0d0
      vimelt1 = 0d0
      freeze_frac = 0d0
      qb_f1 = 0d0
      qbot = 0d0
      hnew2 = 0d0
      vimelt2 = 0d0
      hfinal = 0d0
      vimelt = 0d0

      lheat1 = melt(iw)                           

      if((aint(met(iw,ip_pt)) == 2.or.aint(met(iw,ip_pt)) == 4)         &
            .and.dabs(hsaccum) <= eps)then                              !rain or freezing rain
        if(dabs(hinit) <= eps) then                                     !no ice to start
          if(stt(nnodes) <= Tref.and.lheat1 < 0d0) then
            freeze_frac = dmax1(0d0,dmin1(1d0,lheat(iw)/lheat1))         !based on crrel rep 96-2
            hnew1 = freeze_frac*met(iw,ip_prec)*1d-3*timstep*           &
                                                          (pdens/idens)  !m
            hnew2 = 0d0
          else
            hnew1 = 0d0
            hnew2 = 0d0
          end if
          vimelt1 = 0d0                                                  !m
          vimelt2 = 0d0                                                  !m
        else if(hinit > eps) then                                        !ice to start
          if(toptemp <= Tref.and.lheat1 < 0d0) then                     !freezing
            freeze_frac = dmax1(0d0,dmin1(1d0,lheat(iw)/lheat1))         !based on crrel rep 96-2
            hnew1 = freeze_frac*met(iw,ip_prec)*1d-3*timstep*           &
                                                          (pdens/idens)  !m
            vimelt1 = 0d0
          else
            qtop = -lheat1                                               !W/m^2
            hnew1 = (timstep*3.6d3)*qtop/(idens*lhfus)                   !m
            vimelt1 = dabs(hnew1)*idens*1d-3                             !m (from the top)
          end if

          if(stt(nnodes) > Tref) then                                   !bottom temp > 0
            qb_f1 = ((grthcond(nnodes) + grthcond(nnodes-1))*5d-1)/     &
                                     ((nz(nnodes) - nz(nnodes-1))*5d-1)  !W/m^2 K
            qbot = qb_f1*(stt(nnodes) - Tref)                            !W/m^2
            hnew2 = dmax1(0d0,dmin1(hinit-dabs(hnew1),                  &
                                    (timstep*3.6d3)*qbot/(idens*lhfus))) !m, bottom melt depth
          else
            hnew2 = 0d0
          endif
          vimelt2 = hnew2*idens*1d-3
        end if
      else if(int(met(iw,ip_pt)) == 1.and.dabs(hsaccum) <= eps) then     !no rain
        if(hinit > eps) then                                             !ice to start
          if(toptemp > Tref) then
            qtop = -lheat1                                               !W/m^2
            hnew1 = dmax1(0d0,(timstep*3.6d3)*qtop/(idens*lhfus))
          else
            hnew1 = 0d0
          end if
          vimelt1 = dabs(hnew1)*idens*1d-3                               !m (from the top)

          if(stt(nnodes) > Tref) then                                   !bottom temp > 0
            qb_f1 = ((grthcond(nnodes)+grthcond(nnodes-1))*5d-1)/       &
                                    ((nz(nnodes) - nz(nnodes-1))*5d-1)   !W/m^2 K
            qbot = qb_f1*(stt(nnodes) - Tref)                            !W/m^2
            hnew2 = dmax1(0d0,dmin1(hinit-dabs(hnew1),                  &
                                   (timstep*3600d0)*qbot/(idens*lhfus))) !melt depth bottom temp >0
          else
            hnew2 = 0d0
          endif
          vimelt2 = hnew2*idens*1d-3
        else                                                             !no ice, no melting
          hnew1 = 0d0
          hnew2 = 0d0
          vimelt1 = 0d0
          vimelt2 = 0d0
        end if
      end if

      hfinal = hinit + hnew1 - hnew2
      if(hfinal < 1d-5) hfinal = 0d0
      hfinal = anint(hfinal*1d10)*1d-10

      vimelt = vimelt1 + vimelt2
      vimelt = anint(vimelt*1d15)*1d-15

      if(hinit > eps) then
        refreezei = (hnew1 - hnew2)/hinit
      else if(dabs(hinit) <= eps.and.hfinal > eps) then
        refreezei = (hnew1 - hnew2)/hfinal
      else
        refreezei = 0d0
      end if
      refreezei = anint(refreezei*1d15)*1d-15

      end subroutine icethick
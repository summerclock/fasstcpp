      subroutine read_old_data(iomode,ii,wstart,sdensi,phie,timeo)
 
      use fasst_global

! this subroutine reads in the previous runs' data, writes the final runs output

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: iomode,ii,wstart
      real(kind=8),intent(inout):: sdensi,phie
      real(kind=8),intent(out):: timeo(moverlap,4)

! local variable
      integer(kind=4):: io,i,j,ijunk,kk,junki
      real(kind=8):: junk,jjunk(maxcol)
      character(len=4):: jline
      character(len=2):: jjline


      select case (iomode)

! ******************************************************************************
      case(1)                                                            ! read in the old data

        ijunk = -1
        junki = 0
        io = 0
        junk = 0d0
        do j=1,maxcol
          jjunk(j) = 0d0
        end do
        jline = '    '
        jjline = '  '

        do while(io.ne.-1.and.ijunk.ne.freq_id)
          read(11,iostat=io) ijunk

          if(ijunk.eq.freq_id) then
            kk = 1

            do while(jline /= 'endr'.and.kk <= moverlap)
              if(kk.eq.ii) then
! time variables
                read(11,iostat=io) deltati,stepi,deltat,step,nnodes,    &
                                   ntot,refn

! site variables, surface temperature
                read(11,iostat=io) albedo,emis,sgralbedo,sgremis,       &
                                    hpond,gwl,isurfoldg,isurfoldf

! snow and ice variables
                read(11,iostat=io) hi,hsaccum,sdensi,toptemp,ptemp,     &
                                   iheightn,vsmelt,vimelt,km,sphm,      &
                                   atopf,refreeze,phie
                read(11,iostat=io) (sn_stat(i),i=1,15)
                read(11,iostat=io) sn_istat(1),sn_istat(2)

! low vegetation variables
                read(11,iostat=io) veg_flagl,vegl_type,iseason,isigfl,  &
                                   ihfol,iepf,ifola,ilail,ftemp,storll, &
                                   storls,hfol,hfol_tot,lail,epf,sigfl, &
                                   albf,fola

! canopy variables
                read(11,iostat=io) veg_flagh,vegh_type,zh,isigfh
                do j=1,nclayers
                  read(11,iostat=io) ifoliage_type(j),iclump(j),        &
                                     irho(j),itau(j),ialp(j),ieps(j),   &
                                     ilai(j),laif(j),dzveg(j),          &
                                     canopy_temp(j,oldpos)
                  read(11,iostat=io) (avect(i,j),i=0,4)
                end do 
                read(11,iostat=io) (stor(j),j=1,nclayers+nclayers)

! nodal properties
                do i=1,nnodes
                  read(11,iostat=io) node_type(i)
                  read(11,iostat=io) ntype(i),soil_moist(i),nzi(i),     &
                                     nz(i),delzsi(i),delzs(i),          &
                                     grthcond(i),grspheat(i),frl(i),    &
                                     frh(i),sinkr(i),pheadmin(i)
                  read(11,iostat=io) (nsoilp(i,j),j=1,maxp)
                end do

                do i=1,ntot
                  read(11,iostat=io) stt(i),phead(i),ice(i),wvc(i),     &
                                     vin(i),flowu(i),flowl(i),fv1(i),   &
                                     source(i),sink(i),khl(i),khu(i)
                  read(11,iostat=io) too(i),phoo(i),ioo(i),woo(i),smoo(i)
                end do

! initial met information
                read(11,iostat=io) (dmet1(oldpos,j),j=1,13)
                read(11,iostat=io) (met(oldpos,j),j=1,maxcol)

                read(11,iostat=io) jline

              else !not the correct location                    
                read(11,iostat=io) junk,junki,junk,junki,nnodes,ntot,   &
                                   junki
                read(11,iostat=io) (jjunk(j),j=1,8)

                read(11,iostat=io) (jjunk(j),j=1,13)
                read(11,iostat=io) (jjunk(j),j=1,13)
                read(11,iostat=io) junki,junki

                read(11,iostat=io) junki,junki,junki,(jjunk(j),j=1,15)

                read(11,iostat=io) junki,junki,junk,junk
                do i=1,nclayers
                  read(11,iostat=io) (jjunk(j),j=1,10)
                  read(11,iostat=io) (jjunk(j),j=1,5)
                end do
                read(11,iostat=io) (jjunk(j),j=1,nclayers+nclayers)

                do i=1,nnodes
                  read(11,iostat=io) jjline
                  read(11,iostat=io) junki,(jjunk(j),j=1,11)
                  read(11,iostat=io) (jjunk(j),j=1,maxp)
                end do

                do i=1,ntot
                  read(11,iostat=io) (jjunk(j),j=1,12)
                  read(11,iostat=io) (jjunk(j),j=1,5)
                end do

                read(11,iostat=io) (jjunk(j),j=1,13)
                read(11,iostat=io) (jjunk(j),j=1,maxcol)

                read(11,iostat=io) jline
              end if   !if(kk.eq.ii)
              kk = kk + 1
            end do   !while(jline /= 'endr')

          else
            do while(jline(1:4) /= 'endr')
              read(11,iostat=io) junk,junki,junk,junki,nnodes,ntot
              read(11,iostat=io) (jjunk(j),j=1,8)

              read(11,iostat=io) (jjunk(j),j=1,13)
              read(11,iostat=io) (jjunk(j),j=1,13)
              read(11,iostat=io) junki,junki

              read(11,iostat=io) junki,junki,junki,(jjunk(j),j=1,14)

              read(11,iostat=io) junki,junki,junk,junk
              do i=1,nclayers
                read(11,iostat=io) (jjunk(j),j=1,10)
                read(11,iostat=io) (jjunk(j),j=1,5)
              end do
              read(11,iostat=io) (jjunk(j),j=1,nclayers+nclayers)

              do i=1,nnodes
                read(11,iostat=io) jjline
                read(11,iostat=io) junki,(jjunk(j),j=1,11)
                read(11,iostat=io) (jjunk(j),j=1,maxp)
              end do

              do i=1,ntot
                read(11,iostat=io) (jjunk(j),j=1,12)
                read(11,iostat=io) (jjunk(j),j=1,5)
              end do

              read(11,iostat=io) (jjunk(j),j=1,13)
              read(11,iostat=io) (jjunk(j),j=1,maxcol)

              read(11,iostat=io) jline
            end do
          end if  !if(ijunk.eq.freq_id)
        end do

! ******************************************************************************
      case(2)                                                            ! write data to file used in "hot-start" runs

! write soil properties for hot run
        if(ii == wstart) write(17) freq_id

! time variables
        write(17) deltati,stepi,deltat,step,nnodes,ntot,refn

! site variables, surface temperature
        write(17) albedo,emis,sgralbedo,sgremis,hpond,gwl,isurfoldg,    &
                  isurfoldf

! snow and ice variables
        write(17) hi,hsaccum,sdensi,toptemp,ptemp,iheightn,vsmelt,      &
                  vimelt,km,sphm,atopf,refreeze,phie
        write(17) (sn_stat(i),i=1,15)
        write(17) sn_istat(1),sn_istat(2)

! low vegetation variables
        write(17) veg_flagl,vegl_type,iseason,isigfl,ihfol,iepf,ifola,  &
                  ilail,ftemp,storll,storls,hfol,hfol_tot,lail,epf,     &
                  sigfl,albf,fola

! canopy variables
        write(17) veg_flagh,vegh_type,zh,isigfh
        do j=1,nclayers
          write(17) ifoliage_type(j),iclump(j),irho(j),itau(j),         &
                    ialp(j),ieps(j),ilai(j),laif(j),dzveg(j),           &
                    canopy_temp(j,ii)
          write(17) (avect(i,j),i=0,4)
        end do
        write(17) (stor(j),j=1,nclayers+nclayers)

! node properties
        do i=1,nnodes
          write(17) node_type(i)

          write(17) ntype(i),soil_moist(i),nzi(i),nz(i),delzsi(i),      &
                    delzs(i),grthcond(i),grspheat(i),frl(i),frh(i),     &
                    sinkr(i),pheadmin(i)
          write(17) (nsoilp(i,j),j=1,maxp)
        end do

        do i=1,ntot
          write(17) stt(i),phead(i),ice(i),wvc(i),vin(i),flowu(i),      &
                    flowl(i),fv1(i),source(i),sink(i),khl(i),khu(i)
         write(17) too(i),phoo(i),ioo(i),woo(i),smoo(i)
        end do

! final met properties
        write(17) (dmet1(ii,j),j=1,13)
        write(17) (met(ii,j),j=1,maxcol)

        if(ii /= iend) then
          write(17) 'loop'
        else
          write(17) 'endr'
        end if

! write to the time stamp file
        if(ii == wstart) write(18) freq_id
        write(18) met(ii,ip_year),met(ii,ip_doy),met(ii,ip_hr),         &
                  met(ii,ip_min)
        if(ii /= iend) then
          write(18) 'loop'
        else
          write(18) 'endr'
        end if

! ******************************************************************************
      case(3)                                                            ! read in time stamp of the old data

        ijunk = -1
        io = 0
        junk = 0d0
        jline = '    '

        do while(io.ne.-1.and.ijunk.ne.freq_id)
          read(12,iostat=io) ijunk
          if(ijunk.eq.freq_id) then
            i = 1
            do while(jline(1:4) /= 'endr'.and.i <= moverlap)
              read(12,iostat=io) timeo(i,1),timeo(i,2),timeo(i,3),      &
                                 timeo(i,4),jline
              i = i + 1
            end do
          else
            do while(jline(1:4).ne.'endr')
              read(12,iostat=io) junk,junk,junk,junk,jline
            end do
          end if
        end do

! ******************************************************************************
      end select

      end subroutine read_old_data

      program fasst_driver

      use fasst_global
      use module_zerovars

! The following program is the shell for integrating the single-point and multi-point (btra)
! versions of operation. 

! STEP 0
! initialize all non-allocatable arrays

! STEP 1: Read in the control file that indicates the column location of the 
!           meteorological variables in the met file. Initialize pointers,
!           open the output files.
 
! STEP 2: Read in the meteorological information.

! STEP 3
! Read in old information if it exists, open files & write headers.

! STEP 4: Initialize the soil profile and state of the ground


! STEP 5: Output the data and calculated values

!*******************************************************************************
! UNIT CONVERSION NOTES: W = J/s; J/m = N
!*******************************************************************************

! calls the following subroutines:
!     upr_case (appended to US_soil_tools.f)
!     zero_parameters
!     read_complex (only in multi-point mode)
!     read_old_data(infer_test=1 only)
!     initsurface
!     fasst_main
!     initwater
!     open_water
!     write_outputs
!     makemtagroundout - currently not used
!     makemtafasstout - currently not used

! uses the functions: met_date, map_USDA_SoilType_to_USCS, sort2 (appended to this program)


      implicit none

      integer(kind=4):: io,system                                        !may need to change this line for other compilers
      integer(kind=4):: i,j,jj,k,kk,icnti,ijunk,tflag,err,nt0,d1i
      integer(kind=4):: nm0,stflag,pindex,met_col,tempcount,lis,mlen
      integer(kind=4):: lcount,wtype,npos,num_frq,ntlines,iendi,wstart
      integer(kind=4):: nprint,flprint,sprint,frozen,ltest,itest
      integer(kind=4):: mtsteps,msteps1,tstps,keep,dirpos,fdi
      integer(kind=4):: code1,code2,code3,code4,code5,wflag2,nlayers1
      integer(kind=4):: first_loop,map_USDA_SoilType_to_USCS,fnlen
      integer(kind=4):: nstr_flag(maxn),soiltype1(3)
      integer(kind=4):: soiltype2(maxn+3),soiltype3(maxn+3)
      integer(kind=4):: soiltype4(maxn+3),modis_type(20),umd_type(14)
      integer(kind=4),allocatable:: vitd_indexm(:)

      real(kind=8):: junk,sdensi,phie,met_date,it,et,wvel,wdepth,rtstps
      real(kind=8):: vegint,t1oo,t1o(moverlap),t1n,timeo(moverlap,4)
      real(kind=8):: oldsd,oldhi,pdens,mody,surfdepth,lthick1(3),modyo
      real(kind=8):: d1,d2,lthick2(maxn+3),lthick3(maxn+3),daylim
      real(kind=8):: lthick4(maxn+3) !,pest(maxlines,6)
      real(kind=8),allocatable:: latm(:),longm(:),metm(:,:,:)

      character(len=1):: ttest,mtest,fp 
      character(len=4):: sname(maxl)
      character(len=280):: header
      character(len=239):: vitd_att
      character(len=500):: a
      character(len=500):: fasst_input,met_file_name,fasstusersoil
      character(len=500):: fasstoutmta,groundoutmta,output1,output2
      character(len=500):: output3,output4,vitd_input,output5,ro1,ro2
!      character(len=512):: ctemp
!      character(len=:),allocatable:: fasst_input,met_file_name
!      character(len=:),allocatable:: fasstusersoil,fasstoutmta
!      character(len=:),allocatable:: groundoutmta,vitd_input,output1
!      character(len=:),allocatable:: output2,output3,output4,output5
!      character(len=:),allocatable:: ro1,ro2

      data met_col/32/                                                   !number of columns in the met data file

! convert modis and umd vegetation types to fasst defaults
! biome_sourece: FASST default = 0, MODIS/IGBP = 1000, UMD/LDAS = 2000
!                      1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,
      data modis_type/ 3, 6, 4, 5,18,16,16,18, 7, 2,13,10, 0, 1,12, 8,  &
!                     17,18,19,20
                      14, 9, 9, 9/
!                    1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14
      data umd_type/ 3, 6, 5, 4,18, 3, 2,11,11, 2, 1, 8, 0,15/


      it = 0d0
      et = 0d0
      call cpu_time(it)


! STEP 0
! initialize some variables
      single_multi_flag = 0
      infer_test = 0
      vegint = 0d0
      nprint = 0
      flprint = 0
      sprint = 0
      frozen = 100
      istart = 0
      iend = 0
      oldpos = 0
      mpos = 0
      iendi = 0
      wstart = 0
      wflag2 = 0
      nlayers1 = 0
      fasstusersoil = 'fasstusersoil.inp'
      fasstoutmta = 'fasstout.mta'
      groundoutmta = 'groundout.mta'

      biome_source = 0
      new_vtl = 0
      new_vth = 0
      icnti = 0
      ncols = 0
      pindex = 0
      met_count = 0
      mstflag = 0
      npos = 0
      ijunk = 0
      lis = 0
      fnlen = 0
      mlen = 0
      elev = 0d0
      lat = 0d0
      mlong = 0d0
      mflag = 0d0
      timeoffset = 0d0
      timstep = 0d0
      iheight = 0d0
      junk = 0d0
      rtstps = 0d0
      surfdepth = 0d0
      d1 = 0d0
      d2 = 0d0
      mody = 0d0
      modyo = 0d0
      daylim = 0d0
!      ctemp = ' '


! STEP 1
! This file contains information needed to run the model, for example the soil type,
! number of layers, any initial soil profile data and snow depth information

! if command line argument exists, use its I/O filenames
      call getarg(1,fasst_input)                                         !comment out this or next line depending on OS
!!      call get_command_argument(1,fasst_input)
      fasst_input = fasst_input(1:len_trim(fasst_input))
!      call getarg(1,ctemp)
!      call get_command_argument(1,ctemp)
!      fnlen = len_trim(ctemp)
!      allocate(character(len=fnlen):: fasst_input,stat=err)
!      fasst_input = ' '
!      fasst_input = ctemp(1:fnlen)

!      if(fasst_input == ' ') fasst_input = 'fassti.inp'

!      fasst_input = 'c:\fasst\inp_files\gr1_zip.inp'
!      fasst_input = 'c:\fasst\inp_files\test_frqcytable.inp'

!      fasst_input = 'c:\fasst\inp_files\yuma.inp'
!      fasst_input = 'c:\ols\frd3_fasst_cl.inp'
!      fasst_input = 'c:\fasst\inp_files\rw2003.inp'

!      fasst_input = 'c:\fasst\pest_runs\gr1\gr1.inp'
!      fasst_input = 'c:\aer\chris\ExecutionParameters.txt'
!      fasst_input = 'c:\btra\met_cell_test\test_frqcytable.inp'
!      fasst_input = 'c:\btra\met_cell_test\Lebanon\&
!ExecutionParametersRMS.txt'
!      fasst_input = 'c:\rae\esc2010\hqsc520_veg.inp'
!      fasst_input = 'c:\btra\met_cell_test\Lebanon\FasstRMS.inp'  !2 = no precip
!      fasst_input = 'c:\snow_61\sasp\SASPWinter0910_FASST2c.inp'
!       fasst_input = 'c:\Travers_2008\susan\gr1.inp'
!       fasst_input = 'c:\fasst_support\field_data\crrel_pavement\pvw.inp'
!      fasst_input = 'c:\snow\Tyler\DR_1.inp'

      dirpos = len_trim(fasst_input) + 1

      fp = fasst_input(dirpos:dirpos)
      do while((fp /= '\'.or.fp /= '\').and.dirpos > 1)
        dirpos = dirpos - 1
        fp = fasst_input(dirpos:dirpos)
      end do
      if(dirpos == 0) dirpos = 1

! open the input file
      io = 0
      open(unit=4,file=fasst_input,iostat=io,status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening FASST input file = '',a,             &
              &'' error='',i3)') fasst_input,io
        stop
      end if

! read the input file
      read(4,'(a)') met_file_name
      met_file_name = met_file_name(1:len_trim(met_file_name            &
      (1:index(met_file_name,'!')-1)))
!      ctemp = ' '
!      fnlen = 0
!      read(4,'(a)') ctemp
!      fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!      mlen = fnlen
!      allocate(character(len=fnlen):: met_file_name,stat=err)
!      met_file_name = ' '
!      met_file_name = ctemp(1:fnlen)

      read(4,*) infer_test,single_multi_flag,mgap,keep
      if(infer_test /= 1.and.infer_test /= 0) then
        write(*,'(/,''!!! The infer test flag must be 0 or 1,           & 
              &stopping. !!!'',/)')
        stop
      end if

      oldpos = 1
      istart = 1
      if(infer_test == 1) istart = 2

      read(4,*) gwl,vegint,nprint,flprint,sprint,frozen,mstflag          !mstflag: 0 = AFWA product, don't use measured soil temp
                                                                         !         1 = other
      if(single_multi_flag == 0) then                                    !single point calculation
        call zero_parameters(tflag,nt0,nm0,stflag,lcount,wtype,wvel,    &
                             wdepth,ttest,mtest,sname)

        read(4,*) vitd_index

        read(4,'(a)') output1                                            !met and surface data
        output1 = output1(1:len_trim(output1(1:index(output1,'!')-1)))
!        ctemp = ' '
!        fnlen = 0
!        read(4,'(a)') ctemp                                              !met and surface data
!        fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!        allocate(character(len=fnlen):: output1,stat=err)
!        output1 = ' '
!        output1 = ctemp(1:fnlen)
        read(4,'(a)') output2                                            !vertical profile data
        output2 = output2(1:len_trim(output2(1:index(output2,'!')-1)))
!        ctemp = ' '
!        fnlen = 0
!        read(4,'(a)') ctemp                                              !vertical profile data
!        fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!        allocate(character(len=fnlen):: output2,stat=err)
!        output2 = ' '
!        output2 = ctemp(1:fnlen)
        read(4,'(a)') output3                                            !surface flux data
        output3 = output3(1:len_trim(output3(1:index(output3,'!')-1)))
!        ctemp = ' '
!        fnlen = 0
!        read(4,'(a)') ctemp                                              !surface flux data
!        fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!        allocate(character(len=fnlen):: output3,stat=err)
!        output3 = ' '
!        output3 = ctemp(1:fnlen)
        read(4,'(a)') output4                                            !vegetation temperature output data
        output4 = output4(1:len_trim(output4(1:index(output4,'!')-1)))
!        ctemp = ' '
!        fnlen = 0
!        read(4,'(a)') ctemp                                              !vegetation temperature output data
!        fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!        allocate(character(len=fnlen):: output4,stat=err)
!        output4 = ' '
!        output4 = ctemp(1:fnlen)
        read(4,'(a)') output5                                            !snow processes output data
        output5 = output5(1:len_trim(output5(1:index(output5,'!')-1)))
!        ctemp = ' '
!        fnlen = 0
!        read(4,'(a)') ctemp                                              !snow processes output data
!        fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!        allocate(character(len=fnlen):: output5,stat=err)
!        output5 = ' '
!        output5 = ctemp(1:fnlen)

        read(4,*) slope,aspect                                           !degrees from horizontal, degrees from N
        sloper = dmax1(0d0,dmin1(1.57d0,slope*pi/1.8d2))
        sloper = anint(sloper*1d20)*1d-20
        if(slope > 90d0) then
          write(*,'(/,''!!! Slope can not be greater than 90, '',       & 
                &''stopping. !!!'',/)')
          stop
        end if

        if(abs(aspect) > 360d0) then
          write(*,'(/,''!!! Aspect can not be greater than +/-360, '',  &
                &''stopping. !!!'',/)')
          stop
        end if

        read(4,*) rough                                                  !surface roughness (m)

!       determine if the polygon is open water (1) or land (0)
!       add seperate subroutine for water. Don't go to rest of program.
        read(4,*) water_flag
        if(water_flag /= 1.and.water_flag /= 0) then
          write(*,'(/,''!!! Open water flag must be 0 or 1, '',         &
                &''stopping. !!!'',/)')
          stop
        end if

        if(water_flag == 1) then
          read(4,*) wtype                                                !1 = lakes, ponds;  2 = rivers
          if(wtype == 2) water_flag = 2
          read(4,*) wdepth
          read(4,*) wvel

          soiltype(1) = 26                                               !WAter
          lthick(1) = 2d-2
          soiltype(2) = 26
          soiltype(3) = 7                                                !SM
          if(dabs(aint((wdepth-spflag)*1d5)*1d-5) > eps) then
            lthick(2) = wdepth - lthick(1)
          else
            if(wtype == 1) then
              lthick(2) = 1d1 - lthick(1)
            else
              lthick(2) = 3d0 - lthick(1)
            end if
          end if
          lthick(3) = 1d0
          nlayers = 3
        end if

        if(infer_test == 0) then
          read(4,*) hsaccum,iswe                                         !initial snow depth (m), initial snow water equivalent (m)
          if(hsaccum < 0d0) then
            write(*,'(/,''!!! Initial snow depth can not be less than'',& 
                  &'' 0, stopping. !!!'',/)')
            stop
          end if

          if(iswe < 0d0) then
            write(*,'(/,''!!! Initial snow water equivalent can not '', &
                  &''be less than 0, stopping. !!!'',/)')
            stop
          end if

          read(4,*) hi                                                   !initial ice thickness (m)
          if(hi < 0d0) then
            write(*,'(/,''!!! Initial ice thickness can not be less '', &
                  &''than 0, stopping. !!!'',/)')
            stop
          end if

          if(water_flag == 0) then
            read(4,*) vegl_type                                          !low vegetation type on the surface

            if(vegl_type /= 0) then
              veg_flagl = 1

              if(vegl_type >= 1000.and.vegl_type < 2000) then
                biome_source = 1000                                      !MODIS/IGBP
              else if(vegl_type >= 2000.and.vegl_type < 3000) then 
                biome_source = 2000                                      !UMD/LDAS
              end if
              new_vtl = vegl_type - biome_source

              if(new_vtl > 0) then
                if(biome_source == 1000) then                            !MODIS
                  if(new_vtl /= 13.and.new_vtl /= 16) then
                    vegl_type = modis_type(new_vtl)
                  else
                    vegl_type = 0
                    veg_flagl = 0
                  end if
                else if(biome_source == 2000) then                       !UMD
                  if(new_vtl /= 13) then
                    vegl_type = umd_type(new_vtl)
                  else
                    vegl_type = 0
                    veg_flagl = 0
                  end if
                end if
              else
                vegl_type = 0
                veg_flagl = 0
              end if

              if(vegl_type >= 13.and.vegl_type <= 15) then
                wtype = 1
                wvel = 0d0
                if(vegl_type == 13) then                                 !swamp/bog
                  water_flag = 1
                  wdepth = 1.5d0
                  nlayers1 = 3
                  soiltype1(1) = 26
                  lthick1(1) = 2d-2
                  soiltype1(2) = 26
                  lthick1(2) = wdepth - lthick1(1)
                  soiltype1(3) = 7
                  lthick1(3) = 1d0
                else if(vegl_type == 14) then                            !inland water
                  water_flag = 1
                  veg_flagl = 0
                  vegl_type = 0
                  wdepth = 10d0
                  nlayers1 = 3
                  soiltype1(1) = 26
                  lthick1(1) = 2d-2
                  soiltype1(2) = 26
                  lthick1(2) = wdepth - lthick1(1)
                  soiltype1(3) = 7
                  lthick1(3) = 1d0
                else                                                     !ocean
                  water_flag = 3
                  vegl_type = 0
                  veg_flagl = 0
                end if
              end if

              read(4,*) isigfl                                           !initial foliage density
              if(isigfl > 0.98d0.and.dabs(aint((isigfl-spflag)          &
                                     *1d5)*1d-5) > eps) isigfl = 0.98d0
              read(4,*) iepf                                             !initial foliage emissivity
              read(4,*) ifola                                            !initial foliage absorptivity (1 - albedo)
              read(4,*) ihfol                                            !initial foliage height (cm)

              if(dabs(isigfl) <= eps) then
                vegl_type = 0
                veg_flagl = 0
              end if

              if(vegl_type == 12) veg_flagl = 0                          !ice cap/permanent snow
            end if

            read(4,*) vegh_type                                          !high vegetation type (trees)

            if(vegh_type /= 0) then
              veg_flagh = 1

              if(vegh_type >= 1000.and.vegh_type < 2000) then
                biome_source = 1000
              else if(vegh_type >= 2000.and.vegh_type < 3000) then 
                biome_source = 2000 
              end if
              new_vth = max(0,vegh_type - biome_source)

              if(new_vth > 0) then
                if(biome_source == 1000) then
                  vegh_type = modis_type(new_vth)
                else if(biome_source == 2000) then
                  vegh_type = umd_type(new_vth)
                end if
              else
                vegh_type = 0
                veg_flagh = 0
              end if

              if(vegh_type == 14.or.vegh_type == 15) then                !inland water, ocean
                wtype = 1
                vegh_type = 0
                veg_flagh = 0

                if(vegh_type == 14) then                                 !inland water
                  water_flag = 1
                  wdepth = 10d0
                  nlayers1 = 3
                  soiltype1(1) = 26
                  lthick1(1) = 2d-2
                  soiltype1(2) = 26
                  lthick1(2) = wdepth - lthick1(1)
                  soiltype1(3) = 7
                  lthick1(3) = 1d0
                else                                                     !ocean
                  water_flag = 3
                end if
                wdepth = 3d0
                wvel = 0d0
              end if

              read(4,*) izh,isigfh                                       !canopy height (m), canopy density (0 - 1)

              if(isigfh > 0.98d0.and.dabs(aint((isigfh-spflag)          &
                                     *1d5)*1d-5) > eps) isigfh = 0.98d0
!             for each layer read in leaf area index, Markox clumping factors (typically 1.0), 
!             thickness, reflectance, transmission, shortwave absorption (1 - albedo), 
!             longwave emissivity
              do j=1,nclayers  
                read(4,*) ifoliage_type(j),ilai(j),iclump(j),idzveg(j), & 
                      irho(j),itau(j),ialp(j),ieps(j)
                if(dabs(ilai(j)) <= eps) ilai(j) = 0.5d0
                if(dabs(iclump(j)) <= eps) iclump(j) = 1d-2
              end do

              if(dabs(isigfh) <= eps) then
                vegh_type = 0
                veg_flagh = 0
              end if

              if(vegh_type == 12) veg_flagh = 0                          !ice cap/permanent snow
            end if

!           soil type, number of layers, layer thickness
            read(4,*) nlayers                                            !number of soil layers (max = 10)
            do i=1,nlayers
              read(4,*) soiltype(i),lthick(i)                            !layer type, layer thickness (m)
              if(abs(soiltype(i)) < 100) then                            !USCS soil type
                sclass(i) = 'USCS'
              else if(abs(soiltype(i)) >= 100) then                      !USDA soil type
                sclass(i) = 'USDA'
                if(soiltype(i) > 0) then
                  lis = soiltype(i) - 100
                  soiltype(i) = map_USDA_SoilType_to_USCS(lis)
                else
                  soiltype(i) = -(abs(soiltype(i)) - 100)
                end if
              end if

              if(soiltype(i) == -1) then
                lcount = i
                read(4,'(a)') sname(i)
              end if

              if(soiltype(i) == 26) then
                water_flag = 1
                wflag2 = 1
                wtype = 1
                wvel = 0d0
              end if
            end do

            if(wflag2 == 1.and.nlayers <= 1) then
              nlayers = 3
              soiltype(2) = soiltype(1)
              lthick(2) = lthick(1) - 2d-2
              soiltype(1) = 26
              lthick(1) = 2d-2
              soiltype(3) = 7
              lthick(3) = 1d0
            else if(nlayers1 > 0) then
              kk = 0
              kk = nlayers1 + nlayers
              j = 1
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

              d1 = 0
              do i=1,k
                lthick3(i) = lthick3(i) - d1
                d1 = d1 + lthick3(i)
              end do
              d2 = d1

              do i=1,kk
                k = i - j
                if(k > 0) then
                  if(i <= nlayers1.and.soiltype1(i) /= 26) then
                    lthick2(k) = d1 + lthick1(i)
                    soiltype2(k) = soiltype1(i)
                  else if(i > nlayers1.and.soiltype(i-nlayers) /= 26) then
                    lthick2(k) = d2 + lthick(i-nlayers1)
                    soiltype2(k) = soiltype(i-nlayers1)
                  end if
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
            end if

            if(lcount /= 0) read(4,'(a)') fasstusersoil
            fasstusersoil = fasstusersoil(1:len_trim(fasstusersoil      &
                                        (1:index(fasstusersoil,'!')-1)))
!            if(lcount /= 0) then
!              ctemp = ' '
!              fnlen = 0
!              read(4,'(a)') ctemp
!              fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!              allocate(character(len=fnlen)::fasstusersoil,stat=err)
!              fasstusersoil = ' '
!              fasstusersoil = ctemp(1:fnlen)
!            end if

!           soil temp data info
            read(4,*) ttest                                              !"y" = have initial soil temps in C
            if(ttest == 'y'.or.ttest == 'Y') then
              read(4,*) nt0                                              !number of initial temps
              do j=1,nt0
                read(4,*) zti(j),tm(j)                                   !depth (m), temp (C)
                if(abs(tm(j)) > 100d0) then
                  write(*,'(/,''!!! Measured soil temperatures out of'',&
                        &'' range, stopping. !!!'',/)')
                  stop
                end if
                tm(j) = tm(j) + Tref                                     !convert to Kelvins
              end do
            end if

!           soil moisture data info
            read(4,*) mtest                                              !"y" = have initial soil moisture (fraction)
            if(mtest == 'y'.or.mtest == 'Y') then
              read(4,*) nm0                                              !number of initial temps
              do j=1,nm0
                read(4,*) zm(j),sm(j)                                    !depth (m), moisture (0.0 - 1.0)
                if(sm(j) > 1d0) then
                  write(*,'(/,''!!! Measured soil moistures out of '',  &
                        &''range, stopping. !!!'',/)')
                  stop
                end if
              end do
            end if
          end if   !water_flag = 0
        end if   !infer_test = 0

      else if(single_multi_flag == 1) then                               !multi point calculation
        read(4,'(a)') vitd_input
        vitd_input = vitd_input(1:len_trim(vitd_input(1:              &
                                             index(vitd_input,'!')-1)))

!        ctemp = ' '
!        fnlen = 0
!        read(4,'(a)') ctemp
!        fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!        allocate(character(len=fnlen):: vitd_input,stat=err)
!        vitd_input = ' '
!        vitd_input = ctemp(1:fnlen)
        read(4,'(a)') fasstusersoil
        fasstusersoil = fasstusersoil(1:len_trim(fasstusersoil(1:       &
                                          index(fasstusersoil,'!')-1)))
!        ctemp = ' '
!        fnlen = 0
!        read(4,'(a)') ctemp
!        fnlen = len_trim(ctemp(1:index(ctemp,'!')-1))
!        allocate(character(len=fnlen):: fasstusersoil,stat=err)
!        fasstusersoil = ' '
!        fasstusersoil = ctemp(1:fnlen)
      else
        write(*,'('' Model does not recognize single_multi_flag input'',&
              &'' of = '',i3)') single_multi_flag
      end if
      surfdepth = 5d-2

      close(4)                                                           !close the input file


! STEP 2
!     Read met from appropriate met file. All the met values are read-in. The model
!     loops over each value. The met conditions are the boundary conditions for 
!     solving a series of equations. 

! Open the met file
      io = 0
      open(unit=1,file=met_file_name,iostat=io,status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening met file, error='',i3)') io
        write(*,*) ' '
        write(*,'(a)') met_file_name
        stop
      end if

      read(1,*) icnti,ncols,lat,mlong,elev,pindex,                      & !#_lines, #_col, lat, long, elevation (m), pointer, 
                met_count,timeoffset,timstep,                           & !#_met_sites, offset_from_zulu,data_timestep, 
                mflag,iheight                                             !missing_data_flag, instrumentation_height

      lat = anint(lat*1d10)*1d-10
      mlong = anint(mlong*1d10)*1d-10
      elev = anint(elev*1d10)*1d-10
      timstep = anint(timstep*1d10)*1d-10
      mflag = anint(mflag*1d10)*1d-10

      if(iheight < eps) iheight = 2d0
      iheight = anint(iheight*1d10)*1d-10

! check that the number of met lines is correct
      ntlines = 0
      read(1,'(a)') header
      read(1,'(a)') header

      ltest = 0
      io = 0

      do while(io /= -1.and.ltest == 0)
        read(1,'(a)',iostat=io) a
        i = 1
        itest = 0
        do while(i <= len_trim(a).and.itest == 0)
          if(ichar(a(i:i)) /= 32) then
            if(ichar(a(i:i)) >= 45.and.ichar(a(i:i)) <= 57) then
              ntlines = ntlines + 1
              itest = 1
            else
              exit
            end if
          end if
          i = i + 1
        end do
        if(itest == 0) ltest = 1
      end do

!      icnti = ntlines
!      if(met_count > 1) icnti = icnti - 1
!      icnti = icnti - 1
      if(icnti /= ntlines) icnti = ntlines - 1

      rewind(1)

      dstart = 0
      if(timstep == 3d0.and.met_count > 1) then
        dstart = min(14,icnti)
      else
        dstart = 2
      end if 

      allocate(vitd_indexm(met_count),stat=err)                          !size arrays
      allocate(latm(met_count),longm(met_count),stat=err)                !size arrays
      allocate(metm(met_count,icnti,maxcol),stat=err)                    !size arrays

! zero out the met array
      do i=1,met_count
        do j=1,icnti
          do k=1,maxcol
            metm(i,j,k) = mflag
          end do
        end do
      end do

! read in the met arrray and lat/long/index information
      i = 1
      io = 0
      do while(io /= -1.and.i <= met_count)
        read(1,*,iostat=io) ijunk,ijunk,latm(i),longm(i),junk,          &
                            vitd_indexm(i),ijunk,junk,junk,junk
        latm(i) = anint(latm(i)*1d10)*1d-10
        longm(i) = anint(longm(i)*1d10)*1d-10

        read(1,'(a)',iostat=io) header
        read(1,'(a)',iostat=io) header

        do j=1,icnti
          if(ncols > met_col) then
            read(1,*,iostat=io) (metm(i,j,k),k=1,met_col),metm(i,j,34), &
                                metm(i,j,35)
            do k=1,met_col+2
              metm(i,j,k) = anint(metm(i,j,k)*1d15)*1d-15
            end do
          else
            read(1,*,iostat=io) (metm(i,j,k),k=1,met_col)
            do k=1,met_col
              metm(i,j,k) = anint(metm(i,j,k)*1d15)*1d-15
            end do
          end if

          if(metm(i,j,ip_prec) <= 1d-4) metm(i,j,ip_pt) = 1
          if(metm(i,j,ip_prec2) <= 1d-4) metm(i,j,ip_pt2) = 1
        end do
        i = i + 1
      end do
      if(metm(1,icnti,ip_year)-metm(1,icnti-1,ip_year) > 1) icnti = icnti - 1

      iend = istart + icnti - 1
      iendi = iend
      close(1)                                                           !close the met file


! STEP 3
! read in old information if it exists
! open files, write headers

! default layer properties file
      io = 0
      open(unit=30,file='FASST_defaultlayerprop.inp',iostat=io,         &
                                                      status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening default layerprop file, error = '',  &
              &                                                i3)') io
        write(*,*) ' '
        write(*,'(a)') 'FASST_defaultlayerprop.inp'
        stop
      end if

! default vegetation properties file
      io = 0
      open(unit=31,file='FASST_defaultvegprop.inp',iostat=io,           &
                                                      status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening default vegprop file, error='',i3)') &
                                                                     io
        write(*,*) ' '
        write(*,'(a)') 'FASST_defaultvegprop.inp'
        stop
      end if

! Modis vegetation properties file
      io = 0
      open(unit=32,file='MODIS_IGBP_vegprop.inp',iostat=io,             &
                                                      status='unknown')
      if(io/=0) then
        write(*,'(''Error opening MODIS vegprop file, error = '',i3)')  &
                                                                    io
        write(*,*) ' '
        write(*,'(a)') 'MODIS_IGBP_vegprop.inp'
        stop
      end if

! UMD vegetation properties file
      io = 0
      open(unit=33,file='UMD_LDAS_vegprop.inp',iostat=io,               &
                                                      status='unknown')
      if(io/=0)then
        write(*,'(''Error opening UMD vegprop file, error = '',i3)') io
        write(*,*) ' '
        write(*,'(a)') 'UMD_LDAS_vegprop.inp'
        stop
      end if

! User soil property file; default is 'fasstusersoil.inp' 
      io = 0
      if(lcount == 0.and.single_multi_flag == 0)                       &
        fasstusersoil = 'fasstusersoil.inp'
      open(unit=90,file=fasstusersoil,iostat=io,status='unknown')
      if(io /= 0) then
        write(*,'('' Error opening user soil property file, error = '', &
                   &                                           i3)') io
        write(*,*) ' '
        write(*,'(a)') fasstusersoil
        stop
        end if

! name the outputfiles for the multi-point case
      if(single_multi_flag == 1) then
        io = 0
        open(unit=5,file=vitd_input,iostat=io,status='unknown')
        if(io /= 0) then
          write(*,'('' Error opening complex file, error = '',i3)') io
          write(*,*) ' '
          write(*,'(a)') vitd_input
          stop
        end if

        npos = 500
!        npos = mlen
        do while(met_file_name(npos:npos) /= '.')
          npos = npos - 1
        end do
!        fnlen = 0
!        fnlen = npos + 3

!        allocate(character(len=fnlen)::output1,output2,output3,stat=err)
!        allocate(character(len=fnlen)::output4,output5,stat=err)

        output1 = ' '
        output1 = met_file_name(1:npos)//'fst'                           !met and surface data
        output2 = ' '
        output2 = met_file_name(1:npos)//'grd'                           !nodal information
        output3 = ' '
        output3 = met_file_name(1:npos)//'flx'                           !surface flux information
        output4 = ' '
        output4 = met_file_name(1:npos)//'vgt'                           !vegetation data
        output5 = ' '
        output5 = met_file_name(1:npos)//'snw'                           !snow processes data
      end if

      io = 0
      open(unit=2,file=output1,iostat=io,status='unknown')               !met, surface info output
      if(io /= 0) then
        write(*,'('' Error opening main output, error = '',i3)') io
        write(*,*) ' '
        write(*,'(a)') output1
        stop
      end if

      if(nprint == 1) then
        io = 0
        open(unit=3,file=output2,iostat=io,status='unknown')             !nodal output
        if(io /= 0) then
          write(*,'('' Error opening node output, error = '',i3)') io
        write(*,*) ' '
        write(*,'(a)') output2
          stop
        end if
      end if

      if(flprint == 1) then
        io = 0
        open(unit=24,file=output3,iostat=io,status='unknown')            !flux output
        if(io /= 0) then
          write(*,'('' Error opening flux output, error='',i3)') io
          write(*,*) ' '
          write(*,'(a)') output3
          stop
        end if

        write(24,'(''Surface Flux File'')')
        write(24,'(''Conventions used: '')')
        write(24,'(''     + => flux towards surface'')')
        write(24,'(''     sheat > 0 -> Tair > Tgr => surface warms, air &
              &cools'')')
        write(24,'(''     sheat < 0 -> Tair < Tgr => surface cools, air &
              &warms'')')
        write(24,'(''     lheat > 0 -> qair > qgr => condensation; surfa&
              &ce warms, air cools'')')
        write(24,'(''     lheat < 0 -> qair < qgr => evaporation; surfac&
              &e cools, air warms'')')
        write(24,'('' '')')
      end if

      if(vegint > eps) then
        io = 0
        open(unit=28,file=output4,iostat=io,status='unknown')            !veg output
        if(io /= 0) then
          write(*,'('' Error opening veg output, error = '',i3)') io
          write(*,*) ' '
          write(*,'(a)') output4
          stop
        end if
      end if

      if(sprint == 1) then
        io = 0
        open(unit=55,file=output5,iostat=io,status='unknown')            !snow output
        if(io /= 0) then
          write(*,'('' Error opening snow output, error = '',i3)') io
          write(*,*) ' '
          write(*,'(a)') output5
          stop
        end if
      end if

! open the inferencing flag output file, erase it first if it already exists
      open(unit=10,file='inferred.out',status='unknown')

!      fnlen = dirpos + 17
!      allocate(character(len=fnlen)::ro1,stat=err)
      ro1 = ' '
      ro1 = fasst_input(1:dirpos)//'fasst_ro_soil.dat'

!      fnlen = dirpos + 18
!      allocate(character(len=fnlen)::ro2,stat=err)
      ro2 = ' '
      ro2 = fasst_input(1:dirpos)//'fasst_ro_times.dat'

      open(unit=11,file='fasst_ro_soil.dat',status='unknown',           &
          !  form='binary')                                                !comment out this or next line depending on OS
          form='unformatted',access='stream')
      open(unit=12,file='fasst_ro_times.dat',status='unknown',          &
          !  form='binary')                                                !comment out this or next line depending on OS
          form='unformatted',access='stream')
      open(unit=17,file='temp_soil.dat',status='unknown',               &
          !  form='binary')                                                !comment out this or next line depending on OS
          form='unformatted',access='stream')
      open(unit=18,file='temp_times.dat',status='unknown',              &
          !  form='binary')                                                !comment out this or next line depending on OS
          form='unformatted',access='stream')

! determine number of missing lines in new met file and size arrays accordingly 
      if(single_multi_flag == 0) then
        i = 1
        tflag = 0
        do while(i <= met_count.and.tflag == 0)
          if(vitd_index == vitd_indexm(i)) then
            fdi = freq_id
            tflag = 1
            kk = i
          end if
          i = i + 1
        end do
        if(tflag == 0) then
          write(*,'('' Something wrong with input meteorological file, &
                       &stopping.'')')
          stop
        end if
      else if(single_multi_flag == 1) then
        io = 0
        i = 0
        vitd_att = ' '
        do while(io /= -1)
          read(5,'(a)',iostat=io) vitd_att
          i = i + 1
        end do
        i = i - 2
        num_frq = i
        rewind(5)

        io = 0
        vitd_att = ' '
        read(5,'(a)',iostat=io) vitd_att                                 !header line

        j = 1
        do while(j <= num_frq)
          call read_complex(wtype,io,lcount,nm0,wvel,wdepth,mtest,sname)
          i = 1
          tflag = 0
          do while(i <= met_count.and.tflag == 0)
            if(vitd_index == vitd_indexm(i)) then
              if(j == 1) fdi = freq_id
              tflag = 1
              kk = i
            end if
            i = i + 1
          end do
          j = j + 1
        end do
        rewind(5)
      end if
      i = i - 1

      if(infer_test == 1) then
        d1i = 3
        freq_id = fdi
        wstart = max(moverlap,istart)
        call read_old_data(d1i,dstart,wstart,sdensi,phie,timeo)

        do k=1,moverlap
          t1o(k) = 0d0
          t1o(k) = met_date(timeo(k,1),timeo(k,2),timeo(k,3),timeo(k,4))
        end do
        t1n = 0d0
        t1n = met_date(metm(kk,1,1),metm(kk,1,2),metm(kk,1,3),          &
                       metm(kk,1,4))

        mpos = 0
        do k=1,moverlap-1
          if(t1o(k) < t1n.and.t1o(k+1) >= t1n) mpos = k
        end do
        if(mpos == 0) mpos = moverlap
      end if

! check met file for missing time steps
      mtsteps = 0
      msteps1 = 0
      t1oo = 0d0
      t1n = 0d0

      if(timstep > 0d0) rtstps = 24d0/timstep
      tstps = int(rtstps)

      do k=2,icnti
        t1oo = met_date(metm(i,k-1,1),metm(i,k-1,2),metm(i,k-1,3),      &
                         metm(i,k-1,4))
        t1n = met_date(metm(i,k,1),metm(i,k,2),metm(i,k,3),             &
                         metm(i,k,4))

        mody = metm(i,k,1) - aint(metm(i,k,1)*2.5d-1)*4d0
        modyo = metm(i,k-1,1) - aint(metm(i,k-1,1)*2.5d-1)*4d0

        daylim = 1.042d0  !1.003d0
!        if(mody /= modyo) daylim = 1.07d0
        if(metm(i,k-1,1) /= metm(i,k,1)) daylim = 25.07d0

        msteps1 = 0
        if(dabs(mody) <= eps) then
          msteps1 = int((t1n - t1oo)*367d0*tstps)
          if((t1n - t1oo)*367d0*tstps < daylim) msteps1 = 0
        else
          msteps1 = int((t1n - t1oo)*366d0*tstps)
          if((t1n - t1oo)*366d0*tstps < daylim) msteps1 = 0
        end if
        mtsteps = mtsteps + msteps1
      end do

      if(infer_test == 0) then
        maxlines = iend + mtsteps
      else
        maxlines = iend + mtsteps + moverlap + 1
!        maxlines = iend + (dstart - oldpos) + mtsteps
      end if

      allocate(met(maxlines,maxcol),dmet1(maxlines,13),stat=err)
      allocate(sdens(maxlines),sstate(maxlines),stat=err)
      allocate(canopy_temp(nclayers,maxlines),stat=err)
      allocate(ft(maxlines),tt(maxlines),airt(2,maxlines),stat=err)
      allocate(surfice(maxlines),surfmoist(maxlines),stat=err)
      allocate(surficep(maxlines),surfmoistp(maxlines),stat=err)
      allocate(surfci(maxlines),surfrci(maxlines),stat=err)
      allocate(surfcbr(maxlines),surfd(maxlines),stat=err)
      allocate(frthick(maxlines),twthick(maxlines),stat=err)
      allocate(surfemis(maxlines),surfemisf(maxlines),stat=err)
      allocate(surfemisc(maxlines),slushy(maxlines),stat=err)

! surface energy fluxes
      allocate(cheat(maxlines),cheat1(maxlines),stat=err)
      allocate(sdown(maxlines),sup(maxlines),stat=err)
      allocate(irdown(maxlines),irup(maxlines),stat=err)
      allocate(pheat1(maxlines),lheat(maxlines),stat=err)
      allocate(evap_heat(maxlines),sheat(maxlines),stat=err)
      allocate(melt(maxlines),tmelt(maxlines),stat=err)
      allocate(lheatf(maxlines),lhes(maxlines),stat=err)
      allocate(cevap(maxlines),ponding(maxlines),stat=err)
      allocate(tot_moist(maxlines),stat=err)

! for stdmod output
!      integer(kind=4) kwi(maxlines)
!      real(kind=8) l_strength(maxlines),t_strength(maxlines),t_moist(maxlines)
!      real(kind=8) bedrock


! STEP 4
! main calculation loops
      if(single_multi_flag == 0) then                                    !only one point
        i = 1
        tflag = 0
        do while(i <= met_count.and.tflag == 0)
          if(vitd_index == vitd_indexm(i)) then
            tflag = 1
            lat = latm(i)
            mlong = longm(i)

            call zero_mxl_params

            do k=istart,iend
              kk = k
              if(infer_test == 1) kk = k - 1
              do j=1,maxcol
                met(k,j) = metm(i,kk,j)
              end do

              if(aint(dabs(met(k,ip_zen)-mflag)*1d5)*1d-5 > eps         &
                             .and.dabs(met(k,ip_zen)) > dabs(9d1)) then
                met(k,ip_tsol) = 0d0
                met(k,ip_dir) = 0d0
                met(k,ip_dif) = 0d0
                met(k,ip_upsol) = 0d0
              end if

              if(aint(dabs(met(k,ip_upsol)-mflag)*1d5)*1d-5 > eps) then
                if(met(k,ip_upsol)-met(k,ip_tsol) >  eps.or.            &
                  met(k,ip_upsol) < eps) met(k,ip_upsol) = mflag
              end if
            end do

            code1 = 0
            code2 = 0
            code3 = 0
            code4 = 0
            code5 = 0
            first_loop = 0                                               !first time through fasst_main and/or open_water
            sdensi = 0d0
            phie = 0d0
            oldsd = 0d0
            oldhi = 0d0

            if(water_flag <= 2) then
              call initsurface(nt0,nm0,sprint,sname,ttest,mtest,        &
                               nstr_flag,code1,phie,pdens,oldsd,oldhi,  &
                               sdensi)

              do iw=istart,iend
                call fasst_main(nprint,sprint,frozen,nstr_flag,sdensi,  &
                                surfdepth,code2,code3,code4,code5,      &
                                first_loop,phie,pdens,oldsd,oldhi)
!      pest(iw,1) = soil_moist(nnodes-1) + ice(nnodes)
!      pest(iw,2) = soil_moist(nnodes-3) + ice(nnodes-3)
!      pest(iw,3) = soil_moist(nnodes-6) + ice(nnodes-6)
!      pest(iw,4) = stt(nnodes-1)
!      pest(iw,5) = stt(nnodes-3)
!      pest(iw,6) = stt(nnodes-6)
              end do

! alert the user to any errors that might have occurred during the run
              if(error_code == 1.and.single_multi_flag == 0) then
                write(*,'(/,'' !!!!!! Look in inferred.out for '',      &
                  &''run-time messages, etc. for case # '',i8)') freq_id

                if(code2 == 1) write(10,'(/, '' Maximum interations '', &
                                 &''in canopy profile calculation.'')')
                if(code3 == 1) write(10,'(/, '' Error in surface '',    &
                                 &''energy calculation.'')')
                if(code4 == 1) write(10,'(/, '' Error in snow melt '',  &
                                 &''calculation.'')')
!                if(code5 == 1) write(10,'(/, '' Maximum iterations '',  &
!                                 &''in new profile calculation. '',     &
!                                 &''Occured '',i10,'' times.'')') ecount
              end if
 
              stflag = 1

! ******************************************************************************
! print the met data for PEST
 !     open(unit=60,file='surface_info.dat',status='unknown')
 !     write(60,'(''Year  JD  Hr   M             Grtemp           '',&
 !           &''Grmoist'')')
 !     write(60,'(''                              K'')')

 !      do i=istart,iend
 !        write(60,1005) int(met(i,ip_year)),int(met(i,ip_doy)),&
 !                      int(met(i,ip_hr)),int(met(i,ip_min)),&
 !                      met(i,ip_tsoil),surfmoistp(i)+surficep(i)
!!      write(60,1005) int(met(i,ip_year)),int(met(i,ip_doy)),
!!     &      int(met(i,ip_hr)),int(met(i,ip_min)),met(i,ip_sd),
!!     &      pest(i,1),pest(i,4)
 !     end do

 !     close(60)

! 1005 format(4(i4),2x,f20.16,2x,f18.16)
!! 1005 format(4(i4),2x,f18.16,2x,f18.16,2x,f20.16)

! ******************************************************************************

            else
              call initwater(sprint,phie,pdens,oldsd,oldhi,sdensi)

              do iw=istart,iend
                call open_water(wtype,sprint,wvel,wdepth,sdensi,        &
                                first_loop,code2,code4,oldsd)

              end do
            end if

            call write_outputs(flprint,keep,vegint)
          end if
          i = i + 1
        end do

        if(tflag == 0.and.i <= met_count) then
          write(10,'('' Not finding met data'',i10)') vitd_index
          write(*,'('' Not finding met data'',i10)') vitd_index
        end if

      else if(single_multi_flag == 1) then                               !multi point calculation
        tempcount = 0
        io = 0
        vitd_att = ' '
        read(5,'(a)',iostat=io) vitd_att                                 !header line

        do while(io /= -1)
          tempcount = tempcount + 1
          iend = iendi
          call zero_parameters(tflag,nt0,nm0,stflag,lcount,wtype,wvel,  &
                               wdepth,ttest,mtest,sname)

          call zero_mxl_params

          call read_complex(wtype,io,lcount,nm0,wvel,wdepth,mtest,sname)

          i = 1
          tflag = 0
          do while(i <= met_count.and.tflag == 0)
            if(vitd_index == vitd_indexm(i)) then
              write(*,*)' Running case',tempcount,' of',num_frq,        &
                        ' freq_id = ',freq_id

              tflag = 1
              stflag = 1
              lat = latm(i)
              mlong = longm(i)
              do k=istart,iend
                kk = k
                if(infer_test == 1) kk = k - 1
                do j=1,maxcol
                  met(k,j) = metm(i,kk,j)
                end do

                if(aint(dabs(met(k,ip_zen)-mflag)*1d5)*1d-5 > eps       &
                             .and.dabs(met(k,ip_zen)) > dabs(9d1)) then
                  met(k,ip_tsol) = 0d0
                  met(k,ip_dir) = 0d0
                  met(k,ip_dif) = 0d0
                  met(k,ip_upsol) = 0d0
                end if

                if(aint(dabs(met(k,ip_upsol)-mflag)*1d5)*1d-5           &
                                                            > eps) then
                  if(met(k,ip_upsol) > met(k,ip_tsol).or.               &
                  met(k,ip_upsol) < eps) met(k,ip_upsol) = mflag
                end if
              end do

              code1 = 0
              code2 = 0
              code3 = 0
              code4 = 0
              code5 = 0
              first_loop = 0                                             !first time through fasst_main and/or open_water
              sdensi = 0d0
              phie = 0d0
              oldsd = 0d0
              oldhi = 0d0
              if(water_flag <= 2) then
                call initsurface(nt0,nm0,sprint,sname,ttest,mtest,      &
                                 nstr_flag,code1,phie,pdens,oldsd,      &
                                 oldhi,sdensi)

                do iw=istart,iend
                  call fasst_main(nprint,sprint,frozen,nstr_flag,sdensi,&
                                  surfdepth,code2,code3,code4,code5,    &
                                  first_loop,phie,pdens,oldsd,oldhi)
!      pest(iw,1) = soil_moist(nnodes-1) + ice(nnodes)
!      pest(iw,2) = soil_moist(nnodes-3) + ice(nnodes-3)
!      pest(iw,3) = soil_moist(nnodes-6) + ice(nnodes-6)
!      pest(iw,4) = stt(nnodes-1)
!      pest(iw,5) = stt(nnodes-3)
!      pest(iw,6) = stt(nnodes-6)

                end do

! alert the user to any errors that might have occurred during the run
                if(error_code == 1.and.single_multi_flag == 0) then
                  write(*,'(/,'' !!!!!! Look in inferred.out for '',    &
                    &''run-time messages, etc. for case # '',i8)')      &
                                                                freq_id

                  if(code2 == 1) write(10,'(/, '' Maximum interations'',&
                                 &'' in canopy profile calculation.'')')
                  if(code3 == 1) write(10,'(/, '' Error in surface '',  &
                                 &''energy calculation.'')')
                  if(code4 == 1) write(10,'(/, '' Error in snow melt '',&
                                 &''calculation.'')')
!                  if(code5 == 1) write(10,'(/, '' Maximum iterations '',&
!                                 &''in new profile calculation. '',     &
!                                 &''Occured '',i10,'' times.'')') ecount
!                  if(code5 == 1) write(*,'(/, '' Maximum iterations '', &
!                                 &''in new profile calculation. ''.     &
!                                 &''Occured '',i10,'' times.'')') ecount

                end if
 
                stflag = 1

! ******************************************************************************
! print the met data for PEST
 !     open(unit=60,file='surface_info.dat',status='unknown')
 !     write(60,'(''Year  JD  Hr   M             Grtemp           '',&
 !           &''Grmoist'')')
 !     write(60,'(''                              K'')')

 !      do i=istart,iend
 !        write(60,1005) int(met(i,ip_year)),int(met(i,ip_doy)),&
 !                      int(met(i,ip_hr)),int(met(i,ip_min)),&
 !                      met(i,ip_tsoil),surfmoistp(i)+surficep(i)
!!      write(60,1005) int(met(i,ip_year)),int(met(i,ip_doy)),
!!     &      int(met(i,ip_hr)),int(met(i,ip_min)),met(i,ip_sd),
!!     &      pest(i,1),pest(i,4)
 !     end do

 !     close(60)

! 1005 format(4(i4),2x,f20.16,2x,f18.16)
!! 1005 format(4(i4),2x,f18.16,2x,f18.16,2x,f20.16)

! ******************************************************************************

              else
                call initwater(sprint,phie,pdens,oldsd,oldhi,sdensi)

                do iw=istart,iend
                  call open_water(wtype,sprint,wvel,wdepth,sdensi,      &
                                  first_loop,code2,code4,oldsd)

                end do
              end if

              call write_outputs(flprint,keep,vegint)
            end if
            i = i + 1
          end do

          if((tflag == 0.and.i >= met_count).and.vitd_index /= 0) then
            write(10,'(''freq_id'',i10, '' Not finding met data'',i10)')&
                  freq_id,vitd_index
            write(*,'(''freq_id'',i10, '' Not finding met data'',i10)') &
                  freq_id,vitd_index
          end if
        end do
      end if

! call the meta file creator subroutines
!     if (single_multi_flag == 0)                                       &
!        call makemtagroundout(groundoutmta,3,keep,met)
!     call makemtafasstout(fasstoutmta,2,keep,met)


! STEP 6
! close files, clean shop
      close(2)
      if(nprint == 1) close(3)
      if(single_multi_flag == 1) close(5)
      close(10)
      close(11)
      close(12)
      close(17)
      close(18)
      if(flprint == 1) close(24)
      if(vegint > eps) close(28)
      close(30)
      close(31)
      close(32)
      close(33)
      if(sprint == 1) close(55)
      close(90)

! get rid of "temporary" files
      if(stflag == 1) then
        io = system("mv temp_soil.dat fasst_ro_soil.dat")              !need to change this line for unix/linux/macs (move - mv)
        io = system("mv temp_times.dat fasst_ro_times.dat")            !need to change this line for unix/linux/macs (move - mv)
      else
        io = system("rm temp_soil.dat")                               !need to change this line for unix/linux/macs (del - rm)
        io = system("rm temp_times.dat")                              !need to change this line for unix/linux/macs (del - rm)
      end if

      deallocate(vitd_indexm,latm,longm,metm,stat=err)
      deallocate(met,dmet1,sdens,sstate,canopy_temp,stat=err)
      deallocate(ft,tt,airt,stat=err)
      deallocate(surfice,surfmoist,surficep,surfmoistp,surfd,stat=err)
      deallocate(surfci,surfrci,surfcbr,frthick,twthick,stat=err)
      deallocate(cheat,cheat1,sdown,sup,irdown,irup,stat=err)
      deallocate(pheat1,lheat,evap_heat,sheat,melt,tmelt,stat=err)
      deallocate(surfemis,surfemisf,surfemisc,slushy,stat=err)
      deallocate(lheatf,lhes,cevap,ponding,tot_moist,stat=err)
!      deallocate(fasst_input,met_file_name,fasstusersoil,stat=err)       !only for f2003 version
!      deallocate(vitd_input,output1,output2,output3,stat=err)
!      deallocate(output4,output5,ro1,ro2,stat=err)

      call cpu_time(et)
      write(*,*)' '
      write(*,*) 'total run time: ',et-it

      stop
      end program fasst_driver

! ******************************************************************************
      subroutine sort2(n,arr,arr1,ncount,b,b1)

      use fasst_global
! this subroutine sorts the nodes in ascending order

      implicit none

      integer(kind=4),intent(in):: n
      integer(kind=4),intent(out):: ncount,b1(maxn)
      integer(kind=4),intent(inout):: arr1(maxn+3)
      real(kind=8),intent(out):: b(maxn)
      real(kind=8),intent(inout):: arr(maxn+3)

 
! local variables
      integer(kind=4):: i,j,nm
      real(kind=8):: a

      a = 0d0
      do i=1,maxn
        b(i) = 0d0
        b1(i) = 0
      end do

      do j=2,n
        a = arr(j)
        nm = arr1(j)
        i = j - 1
        do while(i >= 1.and.arr(i) > a)
          arr(i+1) = arr(i)
          arr1(i+1) = arr1(i)
          i = i - 1
        end do
        if(i == 1.and.arr(i) > a) i = 0
        arr(i+1) = a
        arr1(i+1) = nm
      end do

      ncount = 1
      b(ncount) = arr(1)
      b1(ncount) = arr1(1)
      do j=2,n
        a = arr(j)*1000d0
        if(dabs(a-arr(j-1)*1d3) > 1d-4) then !eps) then
          ncount = ncount + 1
          b(ncount) = arr(j)
          b1(ncount) = arr1(j)
        end if
      end do

      end subroutine sort2


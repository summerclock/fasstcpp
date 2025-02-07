      subroutine soil_strength(nstr_flag,ss_strengtho,ll_strengtho,     &    !scio,srcio
                               scbro)

      use fasst_global

! use numerical_libraries   !may need to change or comment out this line for other compilers

! The equations for determining ci and rci come from Steve Grant
! NOTE: need to get the maximum values from him. Add this information as
!       extra columns to coefci and coefrci

! NOTES on CBR (California Bearing Ratio) obtained from the following web sites:
!     www.trb.org/mepdg/2appendices_cc.pdf
!  CBR = 28.09*(d60)^0.358  for coarse materials, i.e., weighted plasticity index (WPI) = 0
!  WPI = (% passing #200)*(Plasticity Index)
!  CBR = 75/(1 + 0.728*WPI), other materials, WPI > 0
!     www.globalsecurity.org/military/library/policy/army/fm/5-430-00-1/APPH.html
!  (CI,CBR) -> (10,0.25), (20,0.5), (60,1.5) => CBR = 0.025CI

! no subroutines called

      implicit none

      integer(kind=4),intent(in):: nstr_flag(maxn)
!     real(kind=8),intent(out):: sci,srci,ss_strength,ll_strength,scbr
      real(kind=8),intent(out):: ss_strengtho,ll_strengtho,scbro  !scio,srcio

! local variables
!      integer(kind=4),parameter:: nvar = 3
!      integer(kind=4),parameter:: nrow = 453
!      integer(kind=4),parameter:: ncol = 4
!      integer(kind=4),parameter:: k = 10                                 !number of nearest neighbors to search for
      integer(kind=4) i,j,scount,lcount,scourse(30)

      integer(kind=4):: ols_flag  !,io,ind(nvar),idiscr(nrow),iknn(k)

!      real(kind=8):: knn(nrow,ncol),part(nrow),xkey(nvar),knn_dist(k)
      real(kind=8):: strci,rci,mc,tot_ci,tot_rci,smsp(6,16),cbr_min(30)
      real(kind=8):: sstrengtho,tot_strengtho,tot_cio,tot_rcio,cbr
      real(kind=8):: lstrengtho,tot_cbro,scio,srcio  !,ss_strengtho,ll_strengtho

      real(kind=8):: coefci(4,14),coefrci(4,14),E,D
      real(kind=8):: cisg,rcisg,sigma0,lambda,a,b,c,smt,pht,w,dct
      real(kind=8):: sigma0t,cisgt,rcisgt,cisgmax,rcisgmax,cbrsg
      real(kind=8):: sstrength,ss_strength,tot_strength,sci,srci
      real(kind=8):: lstrength,ll_strength,scbr,tot_cbr,t1,t2

      real(kind=8):: ols_cbr,tot_olscbr,olscbr,ols_cbrsg,tot_olscbrsg
      real(kind=8):: olscbrsg,plt,f1,cbr_r1,cbr_r2,cbr_r3
      real(kind=8):: ols(2,20),PL(20) !,PLI(20),p4(20),p200(20)
      real(kind=8):: R1(3,20),R2(3,20),R3(3,20),st(11,nnodes)
      real(kind=8):: knn_cbr !,dd,knn_stats(2,nvar),sum_dist
!      character(len=200):: header

      real(kind=8),parameter:: convert = 0.0142230d0                     !convert psi to cm
      real(kind=8),parameter:: m = 1d0                                   !friction factor
      real(kind=8),parameter:: sgbeta = 165d0*pi/180d0                   !angle of the cone penetrometer

!                   CIa       CIb       RCIa      RCIb    CImax RCImax
      data smsp/  8.749d0, -1.1949d0, 12.542d0, -2.955d0, 1.5d2, 1.5d2,&  !GM
                  9.056d0, -1.3566d0, 12.542d0, -2.955d0, 1.5d2, 1.5d2,&  !GC
                  3.987d0,  0.815d0 ,  3.987d0,  0.815d0, 1.5d2, 1.5d2,&  !SW
                  3.987d0,  0.815d0 ,  3.987d0,  0.815d0, 1.5d2, 1.5d2,&  !SP
                  8.749d0, -1.1949d0, 12.542d0, -2.955d0, 1.5d2, 1.5d2,&  !SM
                  9.056d0, -1.3566d0, 12.542d0, -2.955d0, 1.5d2, 1.5d2,&  !SC
                 10.225d0, -1.565d0 , 11.936d0, -2.407d0,   3d2,   3d2,&  !ML
                 10.998d0, -1.848d0 , 15.506d0,  -3.53d0,   3d2,   3d2,&  !CL 
                 10.977d0, -1.754d0 , 17.399d0, -3.584d0,   3d2,   3d2,&  !OL
                 13.816d0, -2.583d0 , 13.686d0, -2.705d0,   3d2,   3d2,&  !CH
                 12.321d0, -2.044d0 , 23.641d0, -5.191d0,   3d2,   3d2,&  !MH
                 13.046d0, -2.172d0 , 12.189d0, -1.942d0,   3d2,   3d2,&  !OH
                      0d0,      0d0 ,      0d0,      0d0, 1.5d1, 1.5d1,&  !Pt 
                  9.056d0, -1.3566d0, 12.542d0, -2.955d0, 1.5d2, 1.5d2,&  !MC (SMSC)
                  9.454d0,  -1.385d0, 14.236d0, -3.137d0,   3d2,   3d2,&  !CM (CLML)
                  3.987d0,   0.815d0,  3.987d0,  0.815d0, 1.5d2, 1.5d2/   !EV

      data scourse/1,1,2,2,1,1,2,2,2,2,2,2,2,2,4,2,2,1,0,3,3,0,0,0,3,0, &
                   0,0,0,0/
     
      data cbr_min/ 60d0, 35d0, 30d0, 15d0, 20d0, 15d0, 20d0, 10d0,&  !GW,GP,GM,GC,SW,SP,SM,SC
                     5d0,  5d0,  1d0,  2d0,  1d0,  0d0, 15d0,  5d0,&  !ML,CL,OL,CH,MH,OH,Pt,MC
                    20d0,  0d0,100d0,100d0,  0d0,  0d0,  0d0,  0d0,&  !CM,EV,XX,CO,AS,XX,XX,XX
                   100d0,  0d0,  0d0,  0d0,  0d0,  1d0/               !RO,WA,AI,XX,XX,SN
       
!                       a         b           c      lambda
      data coefci/      0d0,        0d0,       0d0,         0d0,&  !SW
                        0d0,        0d0,       0d0,         0d0,&  !SP
                  90.2617d0,  57.9645d0,    -248d0,    1.1812d0,&  !SM
                   13.862d0,    175.1d0,  3080.4d0,    0.4803d0,&  !SC
                  -5.9508d0,  74.1596d0,    -100d0,    0.4877d0,&  !ML
                  42.7166d0,    2.935d0, 11.8003d0,    1.0988d0,&  !CL
                    118.1d0,   2.4527d0,  -506.3d0,    0.9892d0,&  !OL
                  33.9469d0,   2.0355d0,  -129.5d0,    1.3644d0,&  !CH
                  75.6996d0,   1.5039d0,  -206.9d0,    1.0707d0,&  !MH
                  87.8898d0,   0.0117d0,  -210.2d0,    1.2245d0,&  !OH
                        0d0,        0d0,       0d0,         0d0,&  !PT
                   13.862d0,    175.1d0,  3080.4d0,    0.4803d0,&  !MC - SMSC
                  -4.3898d0,  75.6072d0,    -500d0,    0.4880d0,&  !CM - CLML
                        0d0,        0d0,       0d0,         0d0/   !EV - evaporites

!                        a         b          c      lambda
      data coefrci/      0d0,       0d0,        0d0,      0d0,&    !SW
                         0d0,       0d0,        0d0,      0d0,&    !SP
                     3.47d-3,    2.5d-4,  -0.0216d0, 2.8795d0,&    !SM
                   59.9976d0,  1.9953d0, -54.8093d0, 0.6354d0,&    !SC
                   -2.3143d0, 42.0819d0,       -1d2, 0.3557d0,&    !ML
                         0d0,       0d0,        0d0,      0d0,&    !CL
                         0d0,       0d0,        0d0,      0d0,&    !OL
                   13.8991d0,     105d0, -17.6119d0, 0.4388d0,&    !CH
                         0d0,       0d0,        0d0,      0d0,&    !MH
                         0d0,       0d0,        0d0,      0d0,&    !OH
                         0d0,       0d0,        0d0,      0d0,&    !PT
                   59.9976d0,  1.9953d0,  -54.809d0, 0.6354d0,&    !MC - SMSC
                   -8.8264d0, 44.0451d0,       -1d2, 0.3553d0,&    !CM - CLML
                         0d0,       0d0,        0d0,      0d0/     !EV - evaporites

! CBR = a*CI**b
!                  a       b         a       b           a      b
      data ols/1.1392d0,0.4896d0,1.1392d0,0.4896d0,1.1392d0,0.4896d0,&      !GW,GP,GM
               1.1392d0,0.4896d0,1.1392d0,0.4896d0,1.1392d0,0.4896d0,&      !GC,SW,SP
               1.1392d0,0.4896d0,1.1392d0,0.4896d0,0.1111d0,0.7390d0,&      !SM,SC,ML
               0.1266d0,0.6986d0,0.1305d0,0.6776d0,0.1264d0,0.6979d0,&      !CL,OL,CH
               0.0820d0,0.7174d0,0.1305d0,0.6776d0,0.2985d0,0.5358d0,&      !MH,OH,Pt
               1.1392d0,0.4896d0,0.1281d0,0.6984d0,1.1392d0,0.4896d0,&      !MC,CM,EV
               1.1392d0,0.4896d0,0.1305d0,0.6776d0/                         !UK,UF

! Plastic Limit
      data PL/999d0,999d0, 27d0,   14d0,999d0,999d0, 27.5d0, 16d0, 27d0,&   !GW,GP,GM,GC,SW,SP,SM,SC,ML
             19.5d0, 42d0, 24d0, 42.5d0, 71d0,166d0, 18.5d0, 17d0,999d0,&   !CL,OL,CH,MH,OH,Pt,MC,CM,EV
             27.5d0, 24d0/                                                  !UK,UF

! Plasticity Index
!      data PLI/999d0,999d0,   9d0,  13d0,999d0,999d0,  5d0,11.5d0,  9d0,&   !GW,GP,GM,GC,SW,SP,SM,SC,ML
!              15.5d0,  3d0,37.5d0,28.5d0, 34d0, 41d0,5.5d0,   6d0,999d0,&   !CL,OL,CH,MH,OH,Pt,MC,CM,EV
!                 5d0,37.5d0/                                                !UK,UF

! percent passing #4 sieve (4.75mm diameter)
!      data p4/ 39d0, 33d0, 50d0, 44d0,  85d0, 73d0, 86d0, 80d0, 97d0,&   !GW,GP,GM,GC,SW,SP,SM,SC,ML
!               93d0, 95d0, 98d0, 95d0,96.5d0,999d0, 86d0, 96d0, 85d0,&   !CL,OL,CH,MH,OH,Pt,MC,CM,EV
!               86d0, 98d0/                                               !UK,UF

! percent passing #200 sieve (0.074mm diameter)
!      data p200/  4d0,   3d0, 20d0, 17d0,   3d0,  3d0, 22d0, 30d0, 70d0,& !GW,GP,GM,GC,SW,SP,SM,SC,ML
!                 65d0,67.5d0, 73d0, 80d0,76.5d0,999d0, 34d0, 65d0,  3d0,& !CL,OL,CH,MH,OH,Pt,MC,CM,EV
!                 22d0,  73d0/                                             !UK,UF

! CBR = a + b*MC**c
! NOTE: This one worked best with the OLS data (TR-08-13, pg 58, Fig 25)
!                a            b            c
      data R1/   999d0,      999d0,     999d0,&   !GW
              7.7814d0, -0.04242d0,   1.649d0,&   !GP
                 999d0,      999d0,     999d0,&   !GM
                 999d0,      999d0,     999d0,&   !GC
                 999d0,      999d0,     999d0,&   !SW
                 999d0,      999d0,     999d0,&   !SP
              7.7814d0,-0.042428d0,   1.649d0,&   !SM
                 999d0,      999d0,     999d0,&   !SC
              -1.815d0,     5480d0,  -2.336d0,&   !ML
              -5.671d0,      135d0, -0.9031d0,&   !CL
                 999d0,      999d0,     999d0,&   !OL
               -4.52d0,     77.7d0,   -0.66d0,&   !CH
              -3.184d0,      733d0,  -1.287d0,&   !MH
                 999d0,      999d0,     999d0,&   !OH
                 999d0,      999d0,     999d0,&   !Pt
              7.7814d0,-0.042428d0,   1.649d0,&   !MC
               -6.09d0,   139.08d0, -0.8953d0,&   !CM  
                 999d0,      999d0,     999d0,&   !EV
            -119.625d0, 130.8733d0,-0.01994d0,&   !UK
              0.5062d0,   405.96d9,  -1.591d0/    !UF

! CBR = a + b*exp(-(MC/PL)/c)
!                   a        b            c
      data R2/     999d0,  999d0,      999d0,&   !GW
                   999d0,  999d0,      999d0,&   !GP
                   999d0,  999d0,      999d0,&   !GM
                   999d0,  999d0,      999d0,&   !GC
                   999d0,  999d0,      999d0,&   !SW
                   999d0,  999d0,      999d0,&   !SP
                -0.418d0, 3771d0,   0.1129d0,&   !SM
                   999d0,  999d0,      999d0,&   !SC
                -0.169d0,165.5d0,   0.2117d0,&   !ML
               -0.5509d0,52.32d0,   0.3828d0,&   !CL
                   999d0,  999d0,      999d0,&   !OL
               -0.5957d0,26.97d0,   0.6256d0,&   !CH
                -2.356d0,24.86d0,   0.5914d0,&   !MH
                   999d0,  999d0,      999d0,&   !OH
                   999d0,  999d0,      999d0,&   !Pt
                -0.418d0, 3771d0,   0.1129d0,&   !MC
                     0d0,72.12d0,   0.3076d0,&   !CM
                   999d0,  999d0,      999d0,&   !EV
               0.02801d0,32.46d0,    0.426d0,&   !UK
              -0.03409d0,31.94d0,   0.5042d0/    !UF

! CBR = 10**(a + b*exp(-MC/c))
!                    a            b            c
      data R3/     999d0,        999d0,      999d0,&   !GW
                   999d0,        999d0,      999d0,&   !GP
                   999d0,        999d0,      999d0,&   !GM
                   999d0,        999d0,      999d0,&   !GC
                   999d0,        999d0,      999d0,&   !SW
                   999d0,        999d0,      999d0,&   !SP
              -1008791d0, 1008792.22d0, 17866100d0,&   !SM
                   999d0,        999d0,      999d0,&   !SC
                -1.273d0,      8.091d0,    12.22d0,&   !ML
                -3.257d0,      5.495d0,    47.16d0,&   !CL
                   999d0,        999d0,      999d0,&   !OL
               -0.8406d0,      3.108d0,    34.82d0,&   !CH
                -1.955d0,      5.031d0,     50.8d0,&   !MH
                   999d0,        999d0,      999d0,&   !OH
                   999d0,        999d0,      999d0,&   !Pt
              -1008791d0, 1008792.22d0, 17866100d0,&   !MC
                -2.453d0,      5.037d0,    33.53d0,&   !CM
                   999d0,        999d0,      999d0,&   !EV
                -1.779d0,      2.628d0,    101.9d0,&   !UK
               -0.9324d0,      2.033d0,    49.99d0/    !UF

! needed for k-nearest neighbor
!      data ind/1,2,3/

! k-nearest neighbor mean and standard deviation
!      data knn_stats/ 11.373d0,   4.41d0,&  !gravimetric moisture content (%)
!                     118.997d0,  9.279d0,&  !dry density (lb/ft^3)
!                      36.874d0,  27.43d0/   !percent passing #200 sieve (%)

! zero-out variables
      ols_flag = 0
      scount = 0
      lcount = 0
      strci = 0d0
      rci = 0d0
      mc = 0d0
      tot_ci = 0d0
      tot_rci = 0d0
      sstrengtho = 0d0
      tot_strengtho = 0d0
      tot_cio = 0d0
      tot_rcio = 0d0
      cbr = 0d0
      lstrengtho = 0d0
      tot_cbro = 0d0
      ss_strengtho  = 0d0
      ll_strengtho = 0d0
      E = 0d0
      D = 0d0
      cisg = 0d0
      rcisg = 0d0
      sigma0 = 0d0
      lambda = 0d0
      a = 0d0
      b = 0d0
      c = 0d0
      smt = 0d0
      pht = 0d0
      w = 0d0
      sigma0t = 0d0
      cisgt = 0d0
      rcisgt = 0d0
      cisgmax = 0d0
      rcisgmax = 0d0
      cbrsg = 0d0
      sstrength = 0d0
      ss_strength = 0d0
      tot_strength = 0d0
      sci = 0d0
      srci = 0d0
      lstrength = 0d0
      ll_strength = 0d0
      scbr = 0d0
      tot_cbr = 0d0
      ols_cbr = 0d0
      tot_olscbr = 0d0
      olscbr = 0d0
      ols_cbrsg = 0d0
      tot_olscbrsg = 0d0
      olscbrsg = 0d0
      plt = 0d0
      cbr_r1 = 0d0
      cbr_r2 = 0d0
      cbr_r3 = 0d0
!      dd = 0d0
!      knn_cbr = 0d0
!      sum_dist = 0d0
      t1 = 0d0
      t2 = 0d0

      scio = 0d0
      srcio = 0d0
      scbro = 0d0

! to normalize values: norm_val = (val - mean_val)/stdev_val
! to convert from g/cm^3 to lb/ft^3: g/cm^3 * lb/454g * (30.48cm/ft)^3 = lb/ft^3
! tangent sigmoid = 2/(1 + exp(-2x)) - 1 ~ tanh   !used for neural net work

! open the output file; write headers
! open the k-nn look-up table; form a k-d tree
      if(iw == istart.and.ols_flag == 1)then
! open a soil strenght information file (for OLS)
        open(unit=54,file='OLS_soil_strength.dat',status='unknown')

        write(54,'(''USCS = USCS soil type'')')
        write(54,'(''Z(m) = depth in meters'')')  
        write(54,'(''THETA = actaul volumetric water content'')')
        write(54,'(''WC = 100 * THETA/porosity (vol/vol)'')')
        write(54,'(''MC = 100 * THETA/dry density (mass/mass)'')')
        write(54,'(''SMSP_CI = cone index based on SMSPII curves'')')
        write(54,'(''SG_CI = cone index based on Steve Grant curves'')')
        write(54,'(''CBR_CI = Debrah CBR vs CI curves using SMSPII '',  & 
        &''CI values'')')
        write(54,'(''CBR_CI_SG = Debrah CBR vs CI curves using Steve '',&
        &''Grant CI values'')')
        write(54,'(''CBR_MC_R1 = Rosa CBR vs MC curve 1'')')
        write(54,'(''CBR_MC_R2 = Rosa CBR vs MC curve 2'')')
        write(54,'(''CBR_MC_R3 = Rosa CBR vs MC curve 3'')')
        write(54,'(''CBR_KNN = k-nearest neighbor CBR'')')
        write(54,'(''CBR_JRAC = JRAC CBR using SMSPII CI curves'')')
        write(54,'(''CBR_JRACSG = JRAC CBR using Steve Grant CI curves''&
        &)')
        write(54,'(''Dry Density (g/cm^3)'')')
        write(54,'(''NOTE: No SMSPII CI vs THETA curves exist for '',   &
        &''gravels or peat. Therefore,'',/,''      a default value of'',&
        &'' 150.0 is used for gravels, 15.0 for unfrozen peat'',/,      &
        &''      and 300.0 for frozen peat.'')')    
        write(54,'('' '')')
        write(54,1000) freq_id,(iend - istart + 1),14,lat,mlong,elev,   &
                       vitd_index,met_count,timeoffset,timstep,mflag
        write(54,'('' USCS     Z(m)      THETA      WC      MC      '', &
        &'' SMSP_CI        SG_CI       CBR_CI    CBR_CI_SG    '',       &
        &''CBR_MC_R1    CBR_MC_R2    CBR_MC_R3      CBR_KNN      '',    &
        &''CBR_JRAC    CBR_JRACSG     DENSITY'')')

!        io = 0
!        open(unit=25,file='c:\OLS\Lookup_Table_Lab_MC_DD_#200_stats_&
!             &normalized.txt',status='unknown',iostat=io)

!        if(io == 0) then
!          read(25,'(a)') header
!          do i=1,nrow
!            read(25,*) (knn(i,j),j=1,ncol)
!          end do
!        else 
!          write(*,'('' Error opening soil strenght lookup file, '',     &
!          &''error='',i3)') io
!          stop
!        end if

!        call quadt(nrow,nvar,ncol,knn,nrow,ind,3,idiscr,part)   !may need to change or comment out this line for other compilers
      end if

      do i=1,nnodes
        if(scourse(ntype(i)) == 3) then                                  !pavements, rock, snow
          strci = 3d2
          rci = 3d2
          cbr = cbr_min(ntype(i))
          ols_cbr = cbr
          ols_cbrsg = cbr
        else if(scourse(ntype(i)) == 0) then                             !air, water
          strci = 0d0
          rci = 0d0
          cbr = cbr_min(ntype(i))
          ols_cbr = cbr
          ols_cbrsg = cbr
        else if(ntype(i) <= 2) then                                      !GP, GW
          strci = 1.5d2
          rci = 1.5d2
          cbr = cbr_min(ntype(i))
          ols_cbr = ols(1,ntype(i))*strci**ols(2,ntype(i))
          ols_cbrsg = ols_cbr
        else if(ntype(i) == 15) then                                     !peat
          strci = smsp(5,ntype(i)-2)
          rci = smsp(6,ntype(i)-2)
          if(ice(i) > 0d0) then                                          !frozen peat
            strci = 3d2
            rci = 3d2
          end if
          cbr = cbr_min(ntype(i))
          ols_cbr = cbr
          ols_cbrsg = cbr
        else                                                             !all other soils
! SMSPII method
          mc = 1d2*soil_moist(i)/nsoilp(i,1)                             !% by weight
          if(mc <= 0d0) mc = 1d-3                                        !saftey net!!!!!!!!

          t1 = smsp(1,ntype(i)-2) + smsp(2,ntype(i)-2)*dlog(mc)
          if(dabs(t1)  < 5d1) strci = dexp(t1)
          if(ntype(i) <= 4) strci = strci*1.1d0
          if(strci > smsp(5,ntype(i)-2)) strci = smsp(5,ntype(i)-2)

          t1 = smsp(3,ntype(i)-2) + smsp(4,ntype(i)-2)*dlog(mc)
          if(dabs(t1)  < 5d1) rci = dexp(t1)
          if(ntype(i) <= 4) rci = rci*1.1d0
          if(rci > smsp(6,ntype(i)-2)) rci = smsp(6,ntype(i)-2)

          cbr = 2d-5*strci*strci + 6d-3*strci + 0.129d0
          if(cbr > 1d2) cbr = 1d2
          cbr = (cbr + cbr_min(ntype(i)))*5d-1

          if(nstr_flag(i) == 0) then
            ols_cbr = ols(1,ntype(i))*(strci**ols(2,ntype(i)))
          else if(nstr_flag(i) == 1) then                                !unknown, coarse-grained
            ols_cbr = ols(1,19)*(strci**ols(2,19))
          else                                                           !unknown, fine-grained
            ols_cbr = ols(1,20)*(strci**ols(2,20))
          end if

! Steve Grant, George Mason Method
! need another head and soil moisture value for steve's method
          smt = 1.01d0*soil_moist(i)
          dct = dcos(pi - sgbeta)/dsin(pi - sgbeta)

          if(smt >= nsoilp(i,9)) then                                    !saturated
             pht = 0d0
          else if(dabs(smt-nsoilp(i,8)) < 1d-3.or.                      &
                                              smt <= nsoilp(i,15)) then  !totally unsaturated
             pht = pheadmin(i)                                           !cm 
          else                                                           !general conditions
             w = (smt - nsoilp(i,8))/(nsoilp(i,9) - nsoilp(i,8))         !unitless
             if(w > eps.and.dabs(w**(-1d0/nsoilp(i,12)) - 1d0) > eps)   &
               pht = (w**(-1d0/nsoilp(i,12)) - 1d0)**(1d0/nsoilp(i,11)) 
             pht = -(1d0/nsoilp(i,10))*pht                               !cm
          end if

! CI
          if(ntype(i) > 4) then
            lambda = coefci(4,ntype(i)-4)
            E = -(dsin(lambda)*dcos(sgbeta) +                           &
                                         m*dcos(lambda)*dsin(sgbeta))/  &
                                          (dcos(lambda) - dcos(sgbeta))
            D = (dsin(lambda) + m*dsin(sgbeta))                         &
                                          /(dcos(lambda) - dcos(sgbeta))

            a = coefci(1,ntype(i)-4)
            if(soil_moist(i) < nsoilp(i,9)) then
              c = 0d0
              b = coefci(2,ntype(i)-4)
            else
              c = coefci(3,ntype(i)-4)
              b = 0d0
            end if

            sigma0 = a + b*(dabs(phead(i))*convert)                     &
                       + c*(soil_moist(i) - nsoilp(i,9))
            sigma0t = a + b*(dabs(pht)*convert) + c*(smt - nsoilp(i,9))

            if(dabs(a) > eps.and.E/sigma0 > 0d0) then
              cisg = (sigma0/dsqrt(3d0))*(3d0*pi - 2d0*sgbeta           &
                      + dasin(m) + dlog(E/(dsqrt(3d0)*sigma0))          &
                      + D*5d-1 + m*dct - dsqrt(1d0 - m*m))
              cisg = dmax1(0d0,dmin1(999d0,cisg))
            else
              cisg = 0d0
            end if
            if(dabs(a) > eps.and.E/sigma0t > 0d0) then
              cisgt = (sigma0t/dsqrt(3d0))*(3d0*pi - 2d0*sgbeta         & 
                       + dasin(m) + dlog(E/(dsqrt(3d0)*sigma0t))        &
                       + D*5d-1 + m*dct - dsqrt(1d0 - m*m))
              cisgt = dmax1(0d0,dmin1(999d0,cisgt))
            else
              cisgt = 0d0
            end if

!         cisgmax = smsp(5,ntype(i)-2)   !temporary
!         if(cisgt-cisg > 0d0.and.cisg < 0.8d0*cisgmax)                  & 
!               cisg = dmax1(0d0,dmin1(999d0,0.8d0*cisgmax))

            cbrsg = 2d-5*cisg*cisg + 6d-3*cisg + 0.129d0
            if(cbrsg > 100d0) cbrsg = 100d0
            cbrsg = (cbrsg + cbr_min(ntype(i)))*5d-1

            if(nstr_flag(i) == 0) then
              if(cisg > eps) ols_cbrsg = ols(1,ntype(i))*               &
                                                  cisg**ols(2,ntype(i))
            else
              if(cisg > eps) ols_cbr = ols(1,19)*cisg**ols(2,19)
            end if

! RCI
            lambda = coefrci(4,ntype(i)-4)
            E = -(dsin(lambda)*dcos(sgbeta) +m*dcos(lambda)             &
                            *dsin(sgbeta))/(dcos(lambda) - dcos(sgbeta))
            D = (dsin(lambda) + m*dsin(sgbeta))                         &
                                          /(dcos(lambda) - dcos(sgbeta))

            a = coefrci(1,ntype(i)-4)
            b = coefrci(2,ntype(i)-4)
            if(soil_moist(i) < nsoilp(i,9)) then
              c = 0d0
            else
              c = coefrci(3,ntype(i)-4)
            end if

            sigma0 = a + b*(-phead(i)*convert)                          & 
                       + c*(soil_moist(i) - nsoilp(i,9))
            sigma0t = a + b*(-pht*convert) + c*(smt - nsoilp(i,9))

            if(dabs(a) > eps.and.E/sigma0 > 0d0) then
              rcisg = (sigma0/dsqrt(3d0))*(3d0*pi - 2d0*sgbeta          & 
                       + dasin(m) + dlog(E/(dsqrt(3d0)*sigma0))         &
                       + D*5d-1 + m*dct - dsqrt(1d0 - m*m))
              rcisg = dmax1(0d0,dmin1(999d0,rcisg))
            else
              rcisg = 0d0
            end if
            if(dabs(a) > eps.and.E/sigma0t > 0d0) then
              rcisgt = (sigma0t/dsqrt(3d0))*(3d0*pi - 2d0*sgbeta        & 
                        + dasin(m) + dlog(E/(dsqrt(3d0)*sigma0t))       &
                        + D*5d-1+ m*dct - dsqrt(1d0 - m*m))
              rcisgt = dmax1(0d0,dmin1(999d0,rcisgt))
            else
              rcisgt = 0d0
            end if

            rcisgmax = smsp(6,ntype(i)-2)                                !temporary
            if(rcisgt-rcisg > 0d0.and.rcisg < 0.8d0*rcisgmax)           & 
               rcisg = dmax1(0d0,dmin1(999d0,0.8d0*rcisgmax))
          end if !steve grant loop, ntype(i) > 4
        end if   !ntype(i).....

        cisg = dmax1(0d0,cisg)
        strci = dmax1(0d0,strci)
        rcisg = dmax1(0d0,rcisg)
        rci = dmax1(0d0,rci)
        cbr = dmax1(0d0,cbr)

! determine strength based on fine/coarse
        if(scourse(ntype(i)) == 1) then                                  !coarse
          sstrength = cisg
          sstrengtho = strci
        else if(scourse(ntype(i)) == 2) then                             !fine
          sstrength = rcisg                                              !rci 
          sstrengtho = rci
        else                                                             !peat, pavements, air, snow
          sstrength = strci
          sstrengtho = strci
        end if

        if(elev-nz(i) <= 0.1524d0) then                                  !6"
!        if(elev-nz(i) <= 0.4d0) then                                     !40 cm
          scount = scount + 1
          tot_strength = tot_strength + sstrength
          tot_ci = tot_ci + cisg
          tot_rci = tot_rci +  rcisg
          tot_cbr = tot_cbr + cbrsg
          tot_strengtho = tot_strengtho + sstrengtho
          tot_cio = tot_cio + strci
          tot_rcio = tot_rcio + rci
          tot_cbro = tot_cbro + cbr 
          tot_olscbr = tot_olscbr + ols_cbr
          tot_olscbrsg = tot_olscbrsg + ols_cbrsg
        else if(elev-nz(i) > 0.1524d0.and.elev-nz(i) <= 0.3048d0) then   !6" - 12"
          lstrength = lstrength + sstrength
          lstrengtho = lstrengtho + sstrengtho
          lcount = lcount + 1
        else if(elev-nz(i) > 0.3048d0.and.lcount == 0) then
          lstrength = lstrength + sstrength
          lstrengtho = lstrengtho + sstrengtho
          lcount = lcount + 1
        end if
 
        if(ols_flag == 1) then
          cbr_r1 = spflag
          cbr_r2 = spflag
          cbr_r3 = spflag

          if(nstr_flag(i) == 0.or.nstr_flag(i) == 3) then
            if(dabs(r1(1,ntype(i))-spflag) > eps) then
              if(mc > eps)cbr_r1 = r1(1,ntype(i)) + r1(2,ntype(i))      &
                                           *(mc**r1(3,ntype(i)))
            else
              if(mc > eps) cbr_r1 = r1(1,19) + r1(2,19)*(mc**r1(3,19))
            end if

            if(dabs(nsoilp(i,22)-spflag) > eps) then
              plt = nsoilp(i,22)
            else
              plt = PL(ntype(i))
            end if

            t1 = (mc/plt)/r2(3,ntype(i))
            t2 = mc/r3(3,ntype(i))
            if(dabs(plt-spflag) > eps) then
              if(dabs(t1) < 5d1)                                        &
                cbr_r2 = r2(1,ntype(i)) + r2(2,ntype(i))*dexp(-t1)
                
              if(dabs(t2) < 5d1) then
                cbr_r3 = 1d1**(r3(1,ntype(i)) + r3(2,ntype(i))          &
                                                            *dexp(-t2))
              end if
            end if
          else if(nstr_flag(i) == 1) then                                !unknown, coarse-grained
            if(mc > eps) cbr_r1 = r1(1,19) + r1(2,19)*(mc**r1(3,19))

            t1 = (mc/plt)/r2(3,19)
            t2 = (mc/PL(19))/r2(3,19)
            if(dabs(nsoilp(i,22)-spflag) > eps) then
              if(dabs(t1) < 5d1) cbr_r2 = r2(1,19) + r2(2,19)*dexp(-t1)
            else
              if(dabs(t1) < 5d1) cbr_r2 = r2(1,19) + r2(2,19)*dexp(-t2)
            end if

            t1 = mc/r3(3,19)
            if(dabs(t1) < 5d1)                                          &
             cbr_r3 = 1d1**(r3(1,19) + r3(2,19)*dexp(-t1))
          else if(nstr_flag(i) == 2) then                                !unknown, fine-grained
            if(mc > eps) cbr_r1 = r1(1,20) + r1(2,20)*(mc**r1(3,20))

            t1 = (mc/plt)/r2(3,20)
            t2 = (mc/PL(20))/r2(3,20)
            if(dabs(nsoilp(i,22)-spflag) > eps) then
              if(dabs(t1) < 5d1) cbr_r2 = r2(1,20) + r2(2,20)*dexp(-t1)
            else
              if(dabs(t2) < 5d1) cbr_r2 = r2(1,20) + r2(2,20)*dexp(-t2)
            end if

            t1 = mc/r3(3,20)
            if(dabs(t1) < 5d1)                                          &
              cbr_r3 = 1d1**(r3(1,20) + r3(2,20)*dexp(-t1))
          end if

! k-nn method
!          if(ntype(i) < 19) then
!            xkey(1) = (mc - knn_stats(1,1))/knn_stats(2,1)
!            dd = nsoilp(i,1)*30.48*30.48*30.48/454.0
!            xkey(2) = (dd - knn_stats(1,2))/knn_stats(2,2)
!            if(dabs(nsoilp(i,23)-spflag) > eps) then
!              xkey(3) = (nsoilp(i,23) - knn_stats(1,3))/knn_stats(2,3)
!            else
!              xkey(3) = (p200(ntype(i)) - knn_stats(1,3))/knn_stats(2,3)
!            end if

!            call nghbr(nvar,xkey,k,nrow,ncol,knn,nrow,ind,3,idiscr,part,&
!                       1,iknn,knn_dist)   !may need to change or comment out this line for other compilers

!            sum_dist = 0d0
!            do j=1,k
!              if(abs(knn_dist(j)) < 5d1) then
!                sum_dist = sum_dist + exp(-knn_dist(j))
!              end if
!            end do

            knn_cbr = 0d0
!            do j=1,k
!              if(dabs(sum_dist) > eps) then
!                 if(abs(knn_dist(j)) < 5d1) then
!                  knn_cbr = knn_cbr + knn(iknn(j),4)*exp(-knn_dist(j)) &
!                                                             /sum_dist
!                  end if
!              else
!                knn_cbr = knn_cbr
!              end if
!            end do
!          else
            knn_cbr = spflag
!          end if

! write OLS info to temporary array to reverse output
          st(1,i) = mc
          st(2,i) = strci
          st(3,i) = cisg
          st(4,i) = ols_cbr
          st(5,i) = ols_cbrsg
          st(6,i) = cbr_r1
          st(7,i) = cbr_r2
          st(8,i) = cbr_r3
          st(9,i) = knn_cbr
          st(10,i) = cbr
          st(11,i) = cbrsg
        end if  !ols_flag == 1
      end do

! output the data
      if(ols_flag == 1) then
        do j=1,nnodes
          i = nnodes - j + 1
          write(54,555) node_type(i),nz(i),soil_moist(i),               &
                      1d2*soil_moist(i)/nsoilp(i,2),st(1,i),st(2,i),    &
                      st(3,i),st(4,i),st(5,i),st(6,i),st(7,i),st(8,i),  &
                      st(9,i),st(10,i),st(11,i),nsoilp(i,1)
        end do

        write(54,*)' '
  555 format(2x,a2,2x,f10.4,1x,f8.4,1x,2(f8.2),1x,10(f12.2,1x),f12.3)

 1000 format(i10,1x,i6,1x,i3,1x,f10.6,1x,f11.6,1x,f11.6,1x,i8,1x,i8,1x, &
             f6.2,1x,f5.2,1x,f8.2)
      end if

! Steve Grant, George Mason method
      if(scount /= 0) then
        f1 = 1d0/float(scount)
        ss_strength = tot_strength*f1
        sci = tot_ci*f1
        srci = tot_rci*f1
        scbr = tot_cbr*f1
        olscbrsg = tot_olscbrsg*f1
      else
        ss_strength = mflag
        sci = mflag
        srci = mflag
        scbr = mflag
        olscbrsg = mflag
      end if

      if(lcount /= 0) then
        ll_strength = lstrength/float(lcount)
      else
        ll_strength = mflag
      end if

! SMSPII method
      if(scount /= 0) then
        ss_strengtho = tot_strengtho/float(scount)
        scio = tot_cio/float(scount)
        srcio = tot_rcio/float(scount)
        scbro = tot_cbro/float(scount)
        olscbr = tot_olscbr/float(scount)
      else
        ss_strengtho = mflag
        scio = mflag
        srcio = mflag
        scbro = mflag
        olscbr = mflag
      end if

      if(lcount /= 0) then
        ll_strengtho = lstrengtho/float(lcount)
      else
        ll_strengtho = mflag
      end if

! close output files
      if(iw == iend.and.single_multi_flag == 0) then
        close(25)
        close(54)
      end if

      end subroutine soil_strength


module fasst_global

! met file pointer information
! The array met(timestep counter,ip_*) contains the met information.
      integer(kind=4),parameter:: ip_year = 1                            !year
      integer(kind=4),parameter:: ip_doy = 2                             !day of year
      integer(kind=4),parameter:: ip_hr = 3                              !hour (local time)
      integer(kind=4),parameter:: ip_min = 4                             !minute (local time)
      integer(kind=4),parameter:: ip_ap = 5                              !air pressure (mbar)
      integer(kind=4),parameter:: ip_tmp = 6                             !air temperature (C)
      integer(kind=4),parameter:: ip_rh = 7                              !relative humidity (%)
      integer(kind=4),parameter:: ip_ws = 8                              !wind speed (m/s)
      integer(kind=4),parameter:: ip_wdir = 9                            !wind direction (+ clockwise from N)
      integer(kind=4),parameter:: ip_prec = 10                           !precipitation rate (mm/hr)
      integer(kind=4),parameter:: ip_pt = 11                             !precipitation type
      integer(kind=4),parameter:: ip_prec2 = 12                          !precipitation snow rate (mm/hr)
      integer(kind=4),parameter:: ip_pt2 = 13                            !precipitation snow type
      integer(kind=4),parameter:: ip_lcd = 14                            !low cloud amount (fractional amount)
      integer(kind=4),parameter:: ip_lhgt = 15                           !low cloud hgt (km)
      integer(kind=4),parameter:: ip_lct = 16                            !low cloud type
      integer(kind=4),parameter:: ip_mcd = 17                            !middle cloud amount (fractional amount)
      integer(kind=4),parameter:: ip_mhgt = 18                           !middle cloud hgt (km)
      integer(kind=4),parameter:: ip_mct = 19                            !middle cloud type
      integer(kind=4),parameter:: ip_hcd = 20                            !high cloud amount (fractional amount)
      integer(kind=4),parameter:: ip_hhgt = 21                           !high cloud hgt (km)
      integer(kind=4),parameter:: ip_hct = 22                            !high cloud type
      integer(kind=4),parameter:: ip_tsol = 23                           !total solar flux (W/m^2)
      integer(kind=4),parameter:: ip_dir = 24                            !direct solar flux (W/m^2)
      integer(kind=4),parameter:: ip_dif = 25                            !diffuse solar flux (W/m^2)
      integer(kind=4),parameter:: ip_upsol = 26                          !reflected solar flux (W/m^2)
      integer(kind=4),parameter:: ip_ir = 27                             !downwelling IR flux (W/m^2)
      integer(kind=4),parameter:: ip_irup = 28                           !upwelling IR flux (W/m^2)
      integer(kind=4),parameter:: ip_zen = 29                            !solar zenith angle (degrees)
      integer(kind=4),parameter:: ip_az = 30                             !solar azimuth angle from North clockwise
      integer(kind=4),parameter:: ip_sd = 31                             !snow depth (m)
      integer(kind=4),parameter:: ip_tsoil = 32                          !soil temperature (C)
      integer(kind=4),parameter:: ip_hi = 33                             !surface ice thickness (m). This is for roads, etc.
      integer(kind=4),parameter:: ip_vis = 34                            !visibility (km). This is a pass-through param.
      integer(kind=4),parameter:: ip_aer = 35                            !aerosol type. This is a pass-through param.

! met data, calculated and raw, and met file info
      integer(kind=4),parameter:: maxcol = 35                            !maximum number of columns of met data
      integer(kind=4),parameter:: moverlap = 15                          !maximum number of lines of metfile overlap
      integer(kind=4):: maxlines,ncols,istart,iend,met_count
      integer(kind=4):: single_multi_flag,water_flag
      real(kind=8):: mflag,timeoffset,elev,lat,mlong,timstep,slope
      real(kind=8):: albedo,emis,iheight,aspect,sloper

      save:: maxlines,ncols,istart,iend,met_count,single_multi_flag
      save:: water_flag,mflag,timeoffset,elev,lat,mlong,timstep,slope
      save:: aspect,albedo,emis,iheight,sloper

! initial conditions, soil structure
      integer(kind=4),parameter:: maxl = 10                              !maximum number of layers           
      integer(kind=4),parameter:: maxn = 100                             !maximum number of nodes
      integer(kind=4),parameter:: maxp = 25                              !maximum number of soil parameters
      integer(kind=4):: nlayers,nnodes,refn,soiltype(maxl)
      integer(kind=4):: ntype(maxn),icourse(maxn)
      real(kind=8):: tm(maxn),zti(maxn),zm(maxn),sm(maxn)
      real(kind=8):: soil_moist(maxn)
      real(kind=8):: grthcond(maxn),grspheat(maxn),nz(maxn),delzs(maxn)
      real(kind=8):: stt(maxn+2),ice(maxn+2),phead(maxn+2),wvc(maxn+2)
      real(kind=8):: sgralbedo,sgremis,gwl,lthick(maxl),rho_fac(maxl)
      real(kind=8):: rough,isoilp(maxl,maxp),soilp(maxl,maxp)
      real(kind=8):: nsoilp(maxn,maxp),nzi(maxn),delzsi(maxn)
      character(len=2):: stype(maxl),node_type(maxn+2)
      character(len=4):: sclass(maxl),nclass(maxn+2)

      save:: nlayers,nnodes,refn,soiltype,ntype,icourse,lthick
      save:: tm,zti,zm,sm,sgralbedo,sgremis,gwl,stt,soil_moist,grthcond
      save:: grspheat,ice,nz,nzi,delzs,delzsi,isoilp,soilp,nsoilp,phead
      save:: wvc,rho_fac,rough,stype,node_type,sclass,nclass

! snow and ice parameters
      integer(kind=4):: sn_istat(2)
      real(kind=8):: hsaccum,hi,dsnow,refreeze,newsd,atopf,km,sphm
      real(kind=8):: refreezei,sn_stat(15),vsmelt,vimelt,iswe,hm

      real(kind=8),parameter:: snalbedo = 0.8d0 !0.7d0                   !new snow albedo (unitless)
      real(kind=8),parameter:: soalbedo = 0.5d0                          !old snow albedo (unitless)
      real(kind=8),parameter:: semis = 0.92d0                            !snow emissivity (unitless)
      real(kind=8),parameter:: sthdiff = 2d-07                           !snow thermal diffusivity (m^2/s)
      real(kind=8),parameter:: sthcond = 0.3492d0                        !snow thermal conductivity (W/m*K)
!      real(kind=8),parameter:: sdensw = 750d0                            !wet snow density (kg/m^3)
      real(kind=8),parameter:: sdensw = 550d0                            !wet snow density (kg/m^3)
      real(kind=8),parameter:: sdensd = 50d0                             !dry snow density (kg/m^3)
      real(kind=8),parameter:: ialbedo = 0.7d0                           !ice albedo (unitless)
      real(kind=8),parameter:: iemis = 0.9d0                             !ice emissivity (unitless)
      real(kind=8),parameter:: ithdiff = 1.167d-06                       !ice thermal diffusivity (m^2/s)              
      real(kind=8),parameter:: ithcond = 2.1648d0                        !ice thermal conductivity (W/m*K)
      real(kind=8),parameter:: idens = 916.5d0                           !ice density (kg/m^3)
      real(kind=8),parameter:: lhfus = 3.335d05                          !latent heat of fusion (J/kg)
      real(kind=8),parameter:: lhsub = 2.838d06                          !latent heat of sublimation (J/kg = (m/s)^2)

      save:: sn_istat,hsaccum,hi,dsnow,refreeze,newsd,atopf,km,sphm
      save:: refreezei,vsmelt,vimelt,iswe,sn_stat,hm

! foliage variables
      integer(kind=4),parameter:: nclayers = 3                           !number of canopy layers
      real(kind=8),parameter:: kveg = 0.38d0                             !average vegetation thermal conductivity (W/m*K)

      integer(kind=4):: vegl_type,vegh_type,veg_flagl,veg_flagh,iseason
      integer(kind=4):: biome_source,new_vtl,new_vth
      real(kind=8):: isigfl,iepf,ifola,ihfol,isigfh,izh,iheightn,zh,Zd
      real(kind=8):: albf,ftemp,z0l,chnf,trmlm,trmhm,rsl,storll,storls
      real(kind=8):: sqrt_chnf,sigfh,uaf,hfol_tot,lail,sigfl,state,fola
      real(kind=8):: ilail,epf,hfol,stll,stls,veg_prp(18,17)
      real(kind=8):: rk(18,maxn),frl(maxn),frh(maxn),sinkr(maxn)
      real(kind=8):: ifoliage_type(nclayers),ilai(nclayers)
      real(kind=8):: iclump(nclayers),irho(nclayers),laif(nclayers)
      real(kind=8):: itau(nclayers),ialp(nclayers),ieps(nclayers)
      real(kind=8):: dzveg(nclayers),storcl(nclayers),storcs(nclayers)
      real(kind=8):: idzveg(nclayers),stor(nclayers+nclayers)
      real(kind=8):: avect(0:4,nclayers)

      save:: vegl_type,vegh_type,veg_flagl,veg_flagh,iseason
      save:: biome_source,new_vtl,new_vth,veg_prp
      save:: rk,frl,frh,ifoliage_type,ilai,iclump,irho,itau,ialp,ieps
      save:: idzveg,laif,dzveg,storcl,storcs,stor,isigfl,iepf,ifola
      save:: ihfol,isigfh,izh,iheightn,zh,Zd,ilail,albf,ftemp,z0l,chnf
      save:: trmlm,trmhm,rsl,storll,storls,sqrt_chnf,sigfh,uaf,hfol_tot
      save:: lail,sigfl,state,fola,epf,hfol,stll,stls,sinkr,avect

! calculation variables
      integer(kind=4),parameter:: mt = 30                                !maximum number of soil types
!      integer(kind=4),parameter:: dstart = 14                            !position in met file where overlap begins
      real(kind=8),parameter:: pi = 3.141592654d0
      real(kind=8),parameter:: sigma = 5.669d-8                          !Stephan-Boltzman const.(W/m^2*K^4)
      real(kind=8),parameter:: spflag = 999d0                            !missing soil, veg parameter flag
      real(kind=8),parameter:: ilim = 0.998d0                            !closeness of ice to maximum water content
      real(kind=8),parameter:: wlim = 0.998d0                            !closeness of watervapor to porosity
      real(kind=8),parameter:: eps = epsilon(1d0)                        !tolerance limit for equality
      real(kind=8),parameter:: vK = 4d-1 !0.35d0                               !von Karmen's constant (unitless)
      real(kind=8),parameter:: Tref = 273.15d0                           !reference temperature (K)
      real(kind=8),parameter:: grav = 9.81d0                             !gravitational acceleration (m/s^2)
      real(kind=8),parameter:: Rv = 8.3143d0/1.8015d-2                   !specific gas constant for water vapor (J/kg*K)
      real(kind=8),parameter:: Rd = 8.3143d0/2.8964d-2                   !specific gas constant for dry air (J/kg*K)

      integer(kind=4):: iw,stepi,error_code,error_type,infer_test,oldpos
      integer(kind=4):: freq_id,vitd_index,ecount,mstflag,step
      integer(kind=4):: icase,icaseo,dstart,mpos
      real(kind=8):: deltat,deltati,pheadmin(maxn),hpond,toptemp,ptemp
      real(kind=8):: mgap,isurfoldg,isurfoldf
      character(len=1):: meltfl

      save:: step,stepi,error_code,error_type,infer_test,freq_id
      save:: oldpos,vitd_index,ecount,mstflag,icase,icaseo,dstart,mpos
      save:: deltat,deltati,pheadmin,hpond,toptemp,ptemp,mgap,meltfl
      save:: isurfoldg,isurfoldf

! profile parameters
      integer(kind=4):: ntot
      real(kind=8):: khu(maxn+2),khl(maxn+2)
      real(kind=8):: vin(maxn+2),flowu(maxn+2),flowl(maxn+2)
      real(kind=8):: fv1(maxn+2),source(maxn+2),sink(maxn+2)
      real(kind=8):: too(maxn+2),smoo(maxn+2),woo(maxn+2),ioo(maxn+2)
      real(kind=8):: phoo(maxn+2),bftm(maxn),dsmdh(maxn+2)

      save:: ntot,khl,khu,vin,flowu,flowl,fv1,source,sink,dsmdh
      save:: too,smoo,woo,ioo,phoo,bftm

! allocatable arrays - dependent on number of lines in the met file
      integer(kind=4),allocatable:: slushy(:)
      real(kind=8),allocatable:: dmet1(:,:),sdens(:),canopy_temp(:,:)
      real(kind=8),allocatable:: airt(:,:),ft(:),tt(:),surfice(:)
      real(kind=8),allocatable:: surfmoist(:),surficep(:),surfmoistp(:)
      real(kind=8),allocatable:: surfci(:),surfrci(:),surfcbr(:),lhes(:)
      real(kind=8),allocatable:: frthick(:),twthick(:),cheat(:),met(:,:)
      real(kind=8),allocatable:: cheat1(:),sdown(:),sup(:),irdown(:)
      real(kind=8),allocatable:: irup(:),pheat1(:),lheat(:),evap_heat(:)
      real(kind=8),allocatable:: sheat(:),melt(:),tmelt(:),surfd(:)
      real(kind=8),allocatable:: surfemis(:),surfemisf(:),surfemisc(:)
      real(kind=8),allocatable:: lheatf(:),cevap(:),ponding(:)
      real(kind=8),allocatable:: tot_moist(:)
      character(len=1),allocatable:: sstate(:)

end module fasst_global

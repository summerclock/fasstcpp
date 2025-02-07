module met_global


! met file column and unit flag information
      integer(kind=4):: y_col,jday_col,hr_col,m_col,ap_col,tmp_col
      integer(kind=4):: rh_col,ws_col,wd_col,prcp_col,pa_col,pt_col
      integer(kind=4):: prcp2_col,pa2_col,pt2_col,lcld_col,lhgt_col
      integer(kind=4):: lct_col,mcld_col,mhgt_col,mct_col,hcld_col
      integer(kind=4):: hhgt_col,hct_col,tsol_col,dirsol_col,difsol_col
      integer(kind=4):: upsol_col,ir_col,irup_col,zen_col,az_col
      integer(kind=4):: sndepth_col,tsoil_col,tmp1_col,tmp2_col,dt_col
      integer(kind=4):: dd_col,precdiam_col,visdis_col,vistype_col
      integer(kind=4):: denmat_col

      save:: y_col,jday_col,hr_col,m_col,ap_col,tmp_col,rh_col,ws_col
      save:: wd_col,prcp_col,pa_col,pt_col,prcp2_col,pa2_col,pt2_col
      save:: lcld_col,lhgt_col,lct_col,mcld_col,mhgt_col,mct_col
      save:: hcld_col,hhgt_col,hct_col,tsol_col,dirsol_col,difsol_col
      save:: upsol_col,ir_col,irup_col,zen_col,az_col,sndepth_col
      save:: tsoil_col,tmp1_col,tmp2_col,dt_col,dd_col,precdiam_col
      save:: visdis_col,vistype_col,denmat_col

! met flags
      character(len=1):: paflag,ttflag,t3flag,lcflag,mcflag,hcflag
      character(len=1):: sdflag1,sdflag2,prcp_flag,prcp2_flag,pa2flag
      character(len=1):: vdflag,pdflag,dflag
      character(len=2):: apflag,tflag,stflag,t2flag,dtflag,ddflag
      character(len=3):: wsflag,pflag,p2flag

      save:: paflag,ttflag,t3flag,lcflag,mcflag,hcflag,sdflag1,sdflag2
      save:: prcp_flag,prcp2_flag,pa2flag,vdflag,pdflag,dflag,apflag
      save:: tflag,stflag,t2flag,dtflag,ddflag,wsflag,pflag,p2flag

! met file pointer information
! The array met(timestep counter,ip_*) contains the met information.
      integer(kind=4),parameter:: ip_year = 1          !year
      integer(kind=4),parameter:: ip_doy = 2           !day of year
      integer(kind=4),parameter:: ip_hr = 3            !hour (local time)
      integer(kind=4),parameter:: ip_min = 4           !minute (local time)
      integer(kind=4),parameter:: ip_ap = 5            !air pressure (mbar)
      integer(kind=4),parameter:: ip_tmp = 6           !air temperature (C)
      integer(kind=4),parameter:: ip_rh = 7            !relative humidity (%)
      integer(kind=4),parameter:: ip_ws = 8            !wind speed (m/s)
      integer(kind=4),parameter:: ip_wdir = 9          !wind direction (+ clockwise from N)
      integer(kind=4),parameter:: ip_prec = 10         !precipitation rate (mm/hr)
      integer(kind=4),parameter:: ip_pt = 11           !precipitation type
      integer(kind=4),parameter:: ip_prec2 = 12        !precipitation snow rate (mm/hr)
      integer(kind=4),parameter:: ip_pt2 = 13          !precipitation snow type
      integer(kind=4),parameter:: ip_lcd = 14          !low cloud amount (fractional amount)
      integer(kind=4),parameter:: ip_lhgt = 15         !low cloud hgt (km)
      integer(kind=4),parameter:: ip_lct = 16          !low cloud type
      integer(kind=4),parameter:: ip_mcd = 17          !middle cloud amount (fractional amount)
      integer(kind=4),parameter:: ip_mhgt = 18         !middle cloud hgt (km)
      integer(kind=4),parameter:: ip_mct = 19          !middle cloud type
      integer(kind=4),parameter:: ip_hcd = 20          !high cloud amount (fractional amount)
      integer(kind=4),parameter:: ip_hhgt = 21         !high cloud hgt (km)
      integer(kind=4),parameter:: ip_hct = 22          !high cloud type
      integer(kind=4),parameter:: ip_tsol = 23         !total solar flux (W/m^2)
      integer(kind=4),parameter:: ip_dir = 24          !direct solar flux (W/m^2)
      integer(kind=4),parameter:: ip_dif = 25          !diffuse solar flux (W/m^2)
      integer(kind=4),parameter:: ip_upsol = 26        !reflected solar flux (W/m^2)
      integer(kind=4),parameter:: ip_ir = 27           !downwelling IR flux (W/m^2)
      integer(kind=4),parameter:: ip_irup = 28         !upwelling IR flux (W/m^2)
      integer(kind=4),parameter:: ip_zen = 29          !solar zenith angle (degrees)
      integer(kind=4),parameter:: ip_az = 30           !solar azimuth angle from North clockwise
      integer(kind=4),parameter:: ip_sd = 31           !snow depth (m)
      integer(kind=4),parameter:: ip_tsoil = 32        !soil temperature (C)
      integer(kind=4),parameter:: ip_vis = 33          !visibility (km). This is a pass-through param.
      integer(kind=4),parameter:: ip_aer = 34          !aerosol type. This is a pass-through param.
      integer(kind=4),parameter:: ip_dom = 35          !density of material (g/cm^3)
      integer(kind=4),parameter:: ip_pd = 36           !precipitation diameter (mm)

! met data, calculated and raw, and met file info
      integer(kind=4),parameter:: maxcol = 36
      integer(kind=4):: icnt,year0,ncols,nlines,mxcol,ave_period
      integer(kind=4):: ndatpos,vitd_index,infer_test
      real(kind=8):: mflag,timeoffset,elev,lat,mlong,timstep,slope
      real(kind=8):: iheight,albedo_fasst

      save:: icnt,year0,ncols,nlines,mxcol,ave_period,ndatpos
      save:: vitd_index,infer_test,mflag,timeoffset,elev,lat,mlong
      save:: timstep,slope,iheight,albedo_fasst

! general parameters
      real(kind=8),parameter:: pi = 3.141592654
      real(kind=8),parameter:: sigma = 5.669e-8     !Stephan-Boltzman const.(W/m^2*K^4)
      real(kind=8),parameter:: eps = epsilon(1d0)   !tolerance limit for equality
      real(kind=8),parameter:: Tref = 273.15d0      !melt temperature

end module met_global

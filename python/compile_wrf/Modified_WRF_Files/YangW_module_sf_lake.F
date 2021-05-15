MODULE module_sf_lake

! The lake scheme was retrieved from the Community Land Model version 4.5 
! (Oleson et al. 2013) with some modifications by Gu et al. (2013). It is a 
! one-dimensional mass and energy balance scheme with 20-25 model layers, 
! including up to 5 snow layers on the lake ice, 10 water layers, and 10 soil 
! layers on the lake bottom. The lake scheme is used with actual lake points and 
! lake depth derived from the WPS, and it also can be used with user defined 
! lake points and lake depth in WRF (lake_min_elev and lakedepth_default). 
! The lake scheme is independent of a land surface scheme and therefore 
! can be used with any land surface scheme embedded in WRF. The lake scheme 
! developments and evaluations were included in Subin et al. (2012) and Gu et al. (2013) 
!
!   Subin et al. 2012: Improved lake model for climate simulations, J. Adv. Model. 
!   Earth Syst., 4, M02001. DOI:10.1029/2011MS000072; 
!   Gu et al. 2013: Calibration and validation of lake surface temperature simulations 
!   with the coupled WRF-Lake model. Climatic Change, 1-13, 10.1007/s10584-013-0978-y. 
!======================================================================================
!MODULE module_sf_lake

! |- SUBROUTINE Lake (Call LakeMain inside)        ! For Coupling with WRF within timeloop
! |- SUBROUTINE LakeMain                           ! Integrated Lake Physical Solutions
! | |- CALL ShalLakeFluxes
! | | |- CALL QSat
! | | |- CALL MoninObukIni     (Or use FrictionVeloc_lake without dust instead)
! | | |- CALL FrictionVelocity (Or use FrictionVeloc_lake without dust instead)
! | | |- CALL QSat
! | |- CALL ShalLakeTemperature
! | | |- CALL SoilThermProp_Lake
! | | |- CALL Tridiagonal
! | | |- CALL PhaseChange_Lake
! | | |- CALL SoilThermProp_Lake
! | |- CALL ShalLakeHydrology
! | |- CALL BuildSnowFilter (use SnowHydro_lake)
! | |- CALL SnowWater
! | |- CALL SnowCompaction
! | |- CALL CombineSnowLayers
! | |- CALL DivideSnowLayers
! | |- CALL BuildSnowFilter

! |- SUBROUTINE ShalLakeFluxes
! |- SUBROUTINE ShalLakeTemperature
! |- SUBROUTINE ShalLakeHydrology
! |- SUBROUTINE SoilThermProp_Lake
! |- SUBROUTINE PhaseChange_Lake
! |- SUBROUTINE QSat
! |- SUBROUTINE Tridiagonal

! |- SUBROUTINE lakeini                           ! Initialization For LakeModel Run
! |- SUBROUTINE lakeini_offline                   ! Initialization For Offline LakeModel Run

! Call lakeini         in module_physics_init.F   <== Initialize LakeModel before integration
! Call lakeini_offline in YangW_Lake_Offline.f90  <== Initialize LakeModel before integration
! Call Lake            in module_surface_driver.F <== Excuting   LakeModel within timeloop
!=====================================================================================
 USE module_wrf_error
 USE module_model_constants, ONLY : rcp  

    implicit none 
    integer, parameter          ::    r8 = selected_real_kind(12) 
    real(r8), public, parameter :: spval = 1.e36  !special value for missing data (ocean)
  

  !--- Some Variables not defined in Gu_Version ---!
    real, parameter :: depth_c = 50.              ! below the level t_lake3d will be 277.0  !mchen

  !--- These are lake physical constants ---!
    integer, parameter  :: nlevsoil             = 10      ! number of soil layers
    integer, parameter  :: nlevlake             = 25      ! number of lake layers
    integer, parameter  :: nlevsnow             = 5       ! maximum number of snow layers

    integer,parameter  ::     lbp               = 1 ! pft-index bounds
    integer,parameter  ::     ubp               = 1
    integer,parameter  ::     lbc               = 1 ! column-index bounds
    integer,parameter  ::     ubc               = 1
    integer,parameter  ::     num_shlakec       = 1 ! number of columns in lake filter
    integer,parameter  ::     filter_shlakec(1) = 1 ! lake filter (columns)
    integer,parameter  ::     num_shlakep       = 1 ! number of pfts in lake filter
    integer,parameter  ::     filter_shlakep(1) = 1 ! lake filter (pfts)
    integer,parameter  ::     pcolumn(1)        = 1  
    integer,parameter  ::     pgridcell(1)      = 1  
    integer,parameter  ::     cgridcell(1)      = 1 ! gridcell index of column
    integer,parameter  ::     clandunit(1)      = 1 ! landunit index of column
    integer,parameter  ::     begg              = 1
    integer,parameter  ::     endg              = 1
    integer,parameter  ::     begl              = 1
    integer,parameter  ::     endl              = 1
    integer,parameter  ::     begc              = 1
    integer,parameter  ::     endc              = 1
    integer,parameter  ::     begp              = 1
    integer,parameter  ::     endp              = 1
    integer,parameter  ::     column            = 1
    logical,parameter  ::     lakpoi(1)         = .true.

    real(r8), parameter :: vkc    = 0.4_r8       !von Karman constant [-]
    real(r8), parameter :: pie    = 3.141592653589793_r8 ! pi
    real(r8), parameter :: grav   = 9.80616_r8   !gravity constant [m/s2]
    real(r8), parameter :: sb     = 5.67e-8_r8   !stefan-boltzmann constant  [W/m2/K4]
    real(r8), parameter :: tfrz   = 273.16_r8    !freezing temperature [K]
    real(r8), parameter :: denh2o = 1.000e3_r8   !density of liquid water [kg/m3]
    real(r8), parameter :: denice = 0.917e3_r8   !density of ice [kg/m3]
    real(r8), parameter :: cpice  = 2.11727e3_r8 !Specific heat of ice [J/kg-K]
    real(r8), parameter :: cpliq  = 4.188e3_r8   !Specific heat of water [J/kg-K]
    real(r8), parameter :: hfus   = 3.337e5_r8   !Latent heat of fusion for ice [J/kg]
    real(r8), parameter :: hvap   = 2.501e6_r8   !Latent heat of evap for water [J/kg]
    real(r8), parameter :: hsub   = 2.501e6_r8+3.337e5_r8 !Latent heat of sublimation    [J/kg]
    real(r8), parameter :: rair   = 287.0423_r8  !gas constant for dry air [J/kg/K]
    real(r8), parameter :: cpair  = 1.00464e3_r8 !specific heat of dry air [J/kg/K]
    real(r8), parameter :: tcrit  = 2.5          !critical temperature to determine rain or snow
    real(r8), parameter :: tkwat  = 0.6          !thermal conductivity of water [W/m/k]
    real(r8), parameter :: tkice  = 2.290        !thermal conductivity of ice   [W/m/k]
    real(r8), parameter :: tkairc = 0.023        !thermal conductivity of air   [W/m/k]
    real(r8), parameter :: bdsno = 250.          !bulk density snow (kg/m**3)
    
  !--- These are physical constants for SnowHydro ---!
    real(r8), parameter :: wimp   = 0.05    !Water impremeable if porosity less than wimp
    real(r8), parameter :: ssi    = 0.033   !Irreducible water saturation of snow
    real(r8), parameter :: cnfac  = 0.5     !Crank Nicholson factor between 0 and 1

  !--- Initialize water type constants ---!
    integer,parameter :: istsoil = 1        !soil         "water" type
    integer, private  :: i                  ! loop index 
    real(r8)          :: dtime              ! land model time step (sec)
    real(r8)          :: zlak(1:nlevlake)   !lake z  (layers)
    real(r8)          :: dzlak(1:nlevlake)  !lake dz (thickness)
    real(r8)          :: zsoi(1:nlevsoil)   !soil z  (layers)
    real(r8)          :: dzsoi(1:nlevsoil)  !soil dz (thickness)
    real(r8)          :: zisoi(0:nlevsoil)  !soil zi (interfaces)  
    real(r8)          :: sand(19)           ! percent sand
    real(r8)          :: clay(19)           ! percent clay
    data(sand(i), i=1,19)/92.,80.,66.,20.,5.,43.,60.,&
                          10.,32.,51., 6.,22.,39.7,0.,100.,54.,17.,100.,92./
    data(clay(i), i=1,19)/ 3., 5.,10.,15.,5.,18.,27.,&
                          33.,33.,41.,47.,58.,14.7,0., 0., 8.5,54.,  0., 3./
    real(r8)          :: watsat(1,nlevsoil) ! volumetric soil water at saturation (porosity)
    real(r8)          :: tksatu(1,nlevsoil) ! thermal conductivity, saturated soil [W/m-K]
    real(r8)          :: tkmg(1,nlevsoil)   ! thermal conductivity, soil minerals  [W/m-K]
    real(r8)          :: tkdry(1,nlevsoil)  ! thermal conductivity, dry soil (W/m/Kelvin)
    real(r8)          :: csol(1,nlevsoil)   ! heat capacity, soil solids (J/m**3/Kelvin)

    CONTAINS

SUBROUTINE Lake( t_phy        ,p8w            ,dz8w         ,qvcurr          ,&  !i
                 u_phy        ,v_phy          , glw         ,emiss           ,&
                 rainbl       ,dtbl           ,swdown       ,albedo          ,&
                 xlat_urb2d   ,z_lake3d       ,dz_lake3d    ,lakedepth2d     ,&
                 watsat3d     ,csol3d         ,tkmg3d       ,tkdry3d         ,&
                 tksatu3d     ,ivgtyp         ,ht           ,xland           ,& 
                 iswater, xice, xice_threshold, lake_min_elev                ,&
                 ids          ,ide            ,jds          ,jde             ,&
                 kds          ,kde            ,ims          ,ime             ,&
                 jms          ,jme            ,kms          ,kme             ,&
                 its          ,ite            ,jts          ,jte             ,&
                 kts          ,kte                                           ,&
                 h2osno2d     ,snowdp2d       ,snl2d        ,z3d             ,&  !h
                 dz3d         ,zi3d           ,h2osoi_vol3d ,h2osoi_liq3d    ,&
                 h2osoi_ice3d ,t_grnd2d       ,t_soisno3d   ,t_lake3d        ,&
                 savedtke12d  ,lake_icefrac3d                                ,& 
                 lakemask                                                    ,& 
                 hfx          ,lh             ,grdflx       ,tsk             ,&  !o
                 qfx          ,t2             ,th2          ,q2              ,&
                !--- YangW: Updated External Variables From the namelist ---!
                !--- beta & extinction coefficient ---!
                 lakedebug_yw,   beta_yw      ,eta_flag     , eta_scale_yw,  eta_yw,&
                !--- roughness & mixing processes ---!
                roughness_flag, diffusivity_flag, diffusivity_index_yw, tdmax, mixing_factor, mixing_factor_ked,&
                !--- YangW: albedo processes ---!
                !--- for calculating coszen ---!
                xlong_lsm, start_yr_yw, start_day_yw, start_sec_yw, coszen_lsm ,&  ! the xlong_lsm & coszen_lsm is input from the module_surface_driver.F itself
                !--- for calculating albedo ---!
                coszen_flag, albedo_flag,                                       &
                sabg_lakeout, cosz_lakeout, albedo_sols_yw, albedo_soll_yw, albedo_solsd_yw, albedo_solld_yw, sabg_factor, alb_wat, alb_ice, & ! YangW: Output 4-kind albedo 
                z0mg_yw, z0hg_yw, z0qg_yw, km_yw, ke_yw, ked_yw, kw_yw)    
!==============================================================================
! This subroutine was first edited by Hongping Gu and Jiming Jin for coupling
! 07/20/2010
!==============================================================================
    IMPLICIT NONE
  !--- input variables ---!
    INTEGER, INTENT(IN)                                 :: ids,ide, jds,jde, kds,kde,  &
                                                           ims,ime, jms,jme, kms,kme,  &
                                                           its,ite, jts,jte, kts,kte
    INTEGER, INTENT (IN)                                :: iswater
    REAL, INTENT(IN)                                    :: xice_threshold
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: LAKEMASK
    REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: XICE
    REAL,           DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)    :: t_phy  
    REAL,           DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)    :: p8w    
    REAL,           DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)    :: dz8w
    REAL,           DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)    :: qvcurr
    REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)    :: U_PHY
    REAL, OPTIONAL, DIMENSION( ims:ime, kms:kme, jms:jme ),INTENT(IN)    :: V_PHY
    REAL,           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)    :: glw
    REAL,           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)    :: emiss
    REAL,           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)    :: rainbl
    REAL,           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)    :: swdown
    REAL,           DIMENSION( ims:ime, jms:jme )         ,INTENT(INOUT) :: albedo
    REAL,           DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)    :: XLAND
    REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme )         ,INTENT(IN)    :: XLAT_URB2D
    INTEGER,        DIMENSION( ims:ime, jms:jme )         ,INTENT(INOUT) :: IVGTYP
    REAL,                                                   INTENT(IN)   :: dtbl
    REAL,           DIMENSION( ims:ime,1:nlevlake,jms:jme ),INTENT(IN)   :: z_lake3d
    REAL,           DIMENSION( ims:ime,1:nlevlake,jms:jme ),INTENT(IN)   :: dz_lake3d
    REAL,           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)   :: watsat3d
    REAL,           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)   :: csol3d
    REAL,           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)   :: tkmg3d
    REAL,           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)   :: tkdry3d
    REAL,           DIMENSION( ims:ime,1:nlevsoil,jms:jme ),INTENT(IN)   :: tksatu3d
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(IN)   :: lakedepth2d    
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(IN)   :: ht
    REAL                                                   ,INTENT(IN)   :: lake_min_elev
  !--- output variables ---!
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(OUT)  :: HFX
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(OUT)  :: LH
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(OUT)  :: GRDFLX
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(OUT)  :: TSK
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(OUT)  :: QFX   
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(OUT)  :: T2
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(OUT)  :: TH2
    REAL,           DIMENSION( ims:ime, jms:jme )          ,INTENT(OUT)  :: Q2
  !--- input & output variables ---!:
    real,           dimension(ims:ime,jms:jme )                 ,intent(inout)  :: savedtke12d 
    real,           dimension(ims:ime,jms:jme )                 ,intent(inout)  :: snowdp2d,       &    
                                                                                   h2osno2d,       &    
                                                                                   snl2d,          &    
                                                                                   t_grnd2d
    
    real,    dimension( ims:ime,1:nlevlake, jms:jme )           ,INTENT(inout)  :: t_lake3d,       &    
                                                                                   lake_icefrac3d
    real,    dimension( ims:ime,-nlevsnow+1:nlevsoil, jms:jme ) ,INTENT(inout)  :: t_soisno3d,     &    
                                                                                   h2osoi_ice3d,   &    
                                                                                   h2osoi_liq3d,   &    
                                                                                   h2osoi_vol3d,   &    
                                                                                   z3d,            &    
                                                                                   dz3d 
    real,    dimension( ims:ime,-nlevsnow+0:nlevsoil, jms:jme ) ,INTENT(inout)  :: zi3d    

  !--- local variable ---!
    REAL     :: SFCTMP,PBOT,PSFC,ZLVL,Q2K,EMISSI,LWDN,PRCP,SOLNET
    INTEGER  :: C,i,j,k
    real(r8) :: SOLDN

  !--- tempory varibles to be initialized ---!
    real(r8)  :: forc_t(1)           ! atmospheric temperature (Kelvin)
    real(r8)  :: forc_pbot(1)        ! atm bottom level pressure (Pa) 
    real(r8)  :: forc_psrf(1)        ! atmospheric surface pressure (Pa)
    real(r8)  :: forc_hgt(1)         ! atmospheric reference height (m)
    real(r8)  :: forc_hgt_q(1)       ! observational height of humidity [m]
    real(r8)  :: forc_hgt_t(1)       ! observational height of temperature [m]
    real(r8)  :: forc_hgt_u(1)       ! observational height of wind [m]
    real(r8)  :: forc_q(1)           ! atmospheric specific humidity (kg/kg)
    real(r8)  :: forc_u(1)           ! atmospheric wind speed in east direction (m/s)
    real(r8)  :: forc_v(1)           ! atmospheric wind speed in north direction (m/s)
   !real(r8)  :: forc_rho(1)         ! density (kg/m**3)
    real(r8)  :: forc_lwrad(1)       ! downward infrared (longwave) radiation (W/m**2)
    real(r8)  :: prec(1)             ! snow or rain rate [mm/s]
    real(r8)  :: sabg(1)             ! solar radiation absorbed by ground (W/m**2)
    real(r8)  :: lat(1)              ! latitude (radians)
    real(r8)  :: z_lake(1,nlevlake)  ! layer depth for lake (m)
    real(r8)  :: dz_lake(1,nlevlake) ! layer thickness for lake (m)
    real(r8)  :: lakedepth(1)        ! column lake depth (m)
    logical   :: do_capsnow(1)       ! true => do snow capping

   !--- Input & output variables ---!
    real(r8)  :: h2osoi_vol(1,-nlevsnow+1:nlevsoil) ! volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
    real(r8)  :: t_grnd(1)                          ! ground temperature (Kelvin)
    real(r8)  :: h2osno(1)                          ! snow water (mm H2O)
    real(r8)  :: snowdp(1)                          ! snow height (m)
    real(r8)  :: z(1,-nlevsnow+1:nlevsoil)          ! layer depth for snow & soil (m)
    real(r8)  :: dz(1,-nlevsnow+1:nlevsoil)         ! layer thickness for soil or snow (m)
    real(r8)  :: t_soisno(1,-nlevsnow+1:nlevsoil)   ! soil (or snow) temperature (Kelvin)
    real(r8)  :: t_lake(1,nlevlake)                 ! lake temperature (Kelvin)
    integer   :: snl(1)                             ! number of snow layers
    real(r8)  :: h2osoi_liq(1,-nlevsnow+1:nlevsoil) ! liquid water (kg/m2)
    real(r8)  :: h2osoi_ice(1,-nlevsnow+1:nlevsoil) ! ice lens (kg/m2)
    real(r8)  :: savedtke1(1)                       ! top level eddy conductivity from previous timestep (W/m.K)
    real(r8)  :: zi(1,-nlevsnow+0:nlevsoil)         ! interface level below a "z" level (m)
    real(r8)  :: lake_icefrac(1,nlevlake)           ! mass fraction of lake layer that is frozen
   !--- output variables ---!
    real(r8)  :: eflx_gnet(1)       !net heat flux into ground (W/m**2)
    real(r8)  :: eflx_lwrad_net(1)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8)  :: eflx_sh_tot(1)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8)  :: eflx_lh_tot(1)     ! total latent heat flux (W/m8*2)  [+ to atm]
    real(r8)  :: t_ref2m(1)         ! 2 m height surface air temperature (Kelvin)
    real(r8)  :: q_ref2m(1)         ! 2 m height surface specific humidity (kg/kg)
    real(r8)  :: taux(1)            ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8)  :: tauy(1)            ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8)  :: ram1(1)            ! aerodynamical resistance (s/m)
                                    ! for calculation of decay of eddy diffusivity with depth
                                    ! Change the type variable to pass back to WRF.

!=== Variables For Updated physical Part ===! 
 !--- YangW: Variables For Updated Roughness & Eddy Diffusivity Calculation ---!
 !--- Added in Subroutine LakeMain -> ShalLakeFluxes      -> FrictionVelocity ---!
 !---                              -> ShalLakeTemperature ->                  ---!
  !--- Local Declaration Variables 
    real(r8) :: z0mg(1)                           ! roughness length over ground, momentum      (m)
    real(r8) :: z0hg(1)                           ! roughness length over ground, momentum      (m)
    real(r8) :: z0qg(1)                           ! roughness length over ground, momentum      (m)
    real(r8) :: ws(1)                             ! surface friction velocity (m/s)
    real(r8) :: ks(1)                             ! coefficient passed to ShalLakeTemperature
    real(r8) :: u2m                               ! 2 m wind speed (m/s)

    real(r8) :: drhodz(nlevlake-1)                ! d [rhow] /dz (kg/m**4)
    !--- YangW: must set n2 >= 7.5e-5, or the ked will be infinity ---! 
    real(r8) :: n2(nlevlake-1)                    ! brunt-vaisala frequency (/s**2)
    real(r8) :: ri(nlevlake-1)                    ! richardson number
    real(r8) :: km                                ! molecular diffusion coefficient (m**2/s)
    real(r8) :: ke(nlevlake-1)                    ! eddy diffusion coefficient (m**2/s)
    real(r8) :: ked(nlevlake-1)                   ! enhanced diffusion coefficient (m**2/s)
    real(r8) :: kw(lbc:ubc,nlevlake)              ! km+ke+ked or km+ke coefficient (m**2/s)
 !-----------------------------------------------------------------!
!===========================================!

  !=== Some Variables Prepared For Coupled WRF-Lake namelist -> Subroutine Lake ===!
    logical, intent(in) :: lakedebug_yw
    real   , intent(in) :: beta_yw   

    !--- For Subroutine Lake -> LakeMain -> ShalLakeFluxes ---!
    logical, intent(in)    :: roughness_flag    

    !--- For Subroutine Lake -> LakeMain -> ShalLakeTemperatures ---!
    logical, intent(in) :: diffusivity_flag       ! True For using the New diffusivity scheme
    integer, intent(in) :: diffusivity_index_yw   ! 1 for kw=mixing_factor*(km+ke+ked)
                                                  ! 2 for kw=km+mixing_factor*ke+ked
    real, intent(in)    :: tdmax
    real, intent(in)    :: mixing_factor    
    real, intent(in)    :: mixing_factor_ked    
    logical, intent(in) :: eta_flag
    real, intent(in)    :: eta_scale_yw
    real, intent(in)    :: eta_yw
    real, intent(in)    :: sabg_factor
    real, intent(in)    :: alb_wat
    real, intent(in)    :: alb_ice

    !=== Albedo processes For Subroutine Lake ===!
      !*** Some InOut Vaiables ***!
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: sabg_lakeout       ! solar radiation absorbed by ground
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: cosz_lakeout       ! the final coszen output: 
                                                                                ! cosz_lakeout=coszen_yw  (when coszen_flag=.true.)
                                                                                ! cosz_lakeout=coszen_lsm (when coszen_flag=.true.)
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: albedo_sols_yw     ! albedo for vis direct  solar radiation
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: albedo_soll_yw     ! albedo for nir direct  solar radiation
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: albedo_solsd_yw    ! albedo for vis diffuse solar radiation
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: albedo_solld_yw    ! albedo for nir diffuse solar radiation
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: z0mg_yw            ! roughness length over ground, momentum      (m)
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: z0hg_yw            ! roughness length over ground, sensible heat (m)
        real, dimension( ims:ime,jms:jme ) ,intent(inout) :: z0qg_yw            ! roughness length over ground, latent heat   (m)

        real, dimension( ims:ime,jms:jme )              , intent(inout)  :: km_yw   ! molecular diffusion coefficient (m**2/s)
        real, dimension( ims:ime,1:nlevlake-1, jms:jme ), intent(inout)  :: ke_yw   ! eddy diffusion coefficient (m**2/s) 
        real, dimension( ims:ime,1:nlevlake-1, jms:jme ), intent(inout)  :: ked_yw  ! enhanced diffusion coefficient (m**2/s) 
        real, dimension( ims:ime,1:nlevlake, jms:jme ),   intent(inout)  :: kw_yw   ! km+ke+ked or km+ke coefficient (m**2/s) 
     
      !--- Some namelist variables ---!
        real, dimension( ims:ime, jms:jme ), intent(in)  :: xlong_lsm 
        real, dimension( ims:ime, jms:jme ), intent(in)  :: coszen_lsm
        real(r8)                                         :: coszen_lsm_temp(1)
        real(r8)                                         :: lon(1)              ! longitude (radians) read from xlong_lsm
        logical                            , intent(in)  :: coszen_flag         ! .true. for using the calculated cozsen, false for using the coszen_lsm
        integer                            , intent(in)  :: start_yr_yw         ! Julian Year For decreasing with the Integration
        integer                            , intent(in)  :: start_day_yw        ! Julian day  For decreasing with the Integration
        integer                            , intent(in)  :: start_sec_yw        ! Julian sec  For decreasing with the Integration
        logical                            , intent(in)  :: albedo_flag         ! true for using the new albedo scheme
      !--- Some Local Variables ---!
        integer               :: idate_indice = 1      ! 1 for Initilizing the idate(3) only once
        integer               :: idate(3)              ! year, julian day, seconds of the starting time	   
        real(r8)              :: calday                ! Julian cal day (1.xx to 365.xx)
        real(r8)              :: coszen_yw             ! YW calculated coszen
        real(r8)              :: coszen                ! cosine of the solar zenith angle defined by coszen_flag
        real(r8)              :: forc_sols             ! atm visible direct beam solar rad onto srf             [W/m2]	 
        real(r8)              :: forc_soll             ! atm Near-infrared (NIR) direct beam solar rad onto srf [W/m2]
        real(r8)              :: forc_solsd            ! atm vis diffuse solar rad onto srf [W/m2]	   
        real(r8)              :: forc_solld            ! atm nir diffuse solar rad onto srf [W/m2] 
        !--- Following sag_yw is calculated based on BATS, different from the snowage in ShalLakeHydrology ---!
        real(r8)              :: h2osno_old(1)         ! snow water (mm H2O) at the initilized or last integrated time
        integer               :: albedo_indice = 1     ! 1 for Giving the initialized sag_yw and albedo for initilizing sabg(radiation absorded by ground) 
        real(r8)              :: sag_yw                ! non-dimensional snow age
        real(r8), parameter   :: zlnd   = 0.01         ! Roughness length for soil [m]
        real(r8)              :: fsno_yw               ! snow fraction defined by YW
        real(r8)              :: alb_yw(2,2)           ! averaged albedo according to different landtype and wave lendths of solar radiation [-]
                                                       ! (1,1) for vis direct  radiation
                                                       ! (1,2) for vis diffuse radiation
                                                       ! (2,1) for nir direct  radiation
                                                       ! (2,2) for nir diffuse radiation
    !==============================================!
  !================================================================================!

  !=== YangW: For Debug ====!
   if (lakedebug_yw) then
       print *, "****** YangW Following are Lake Debug Processes ******"
       print *, "beta_yw = ",beta_yw
       print *, "roughness_flag = ",roughness_flag
       print *, "diffusivity_flag = ",diffusivity_flag
       print *, "tdmax = ",tdmax
       print *, "mixing_factor = ",mixing_factor
       print *, "eta_flag = ",eta_flag
       print *, "eta_scale_yw = ",eta_scale_yw
       print *, "eta_yw = ",eta_yw
       print *, "******************************************************"
   endif
  !=========================!

  dtime = dtbl  ! model Integrated time step, Used to get precipitation [mm/s]

  !=== YangW Updated: Get current time of model run ===!
  !--- Attention: The calday must be calculated at each time loop ---!
  !               And calday shouldn't be calculated in the i,j lakegrid loop        Pay Attention !!!
  !               Then the coszen can be calculated at the specified lat(c),lon(c)
    if (idate_indice < 2) then
        idate(1)     = start_yr_yw
        idate(2)     = start_day_yw
        idate(3)     = start_sec_yw
        idate_indice = idate_indice+1
    end if
    !=== YangW: For Debug ====!
     if (lakedebug_yw) then
         print *, "****** YangW Following are Lake Debug Processes ******"
         print *, "*** Check Whether the idate is updated at each time step ***"
         print *, "idate_indice=",idate_indice
         print *, "idate(1) = ",idate(1),"idate(2) = ",idate(2),"idate(3) = ",idate(3)
         print *, "******************************************************"
     endif
    !=========================!
  !====================================!

    DO J = jts,jte
       DO I = its,ite

          SFCTMP  = t_phy(i,1,j)
          PBOT    = p8w(i,2,j)
          PSFC    = P8w(i,1,j) 
          ZLVL    = 0.5 * dz8w(i,1,j) 
          Q2K     = qvcurr(i,1,j)/(1.0 + qvcurr(i,1,j))
          EMISSI  = EMISS(I,J) 
          LWDN    = GLW(I,J)*EMISSI 
          PRCP    = RAINBL(i,j)/dtbl
          SOLDN   = SWDOWN(I,J)                        ! SOLDN is total incoming solar
         !SOLNET  = SOLDN*(1.-ALBEDO(I,J))             ! use mid-day albedo to determine net downward solar ! Annotated by YW
                                                       ! (no solar zenith angle correction) 
     if (lakemask(i,j).eq.1) THEN
         do c = 1,column
            forc_t(c)          = SFCTMP           ! [K]
            forc_pbot(c)       = PBOT 
            forc_psrf(c)       = PSFC
            forc_hgt(c)        = ZLVL             ! [m]
            forc_hgt_q(c)      = ZLVL             ! [m]
            forc_hgt_t(c)      = ZLVL             ! [m]
            forc_hgt_u(c)      = ZLVL             ! [m]
            forc_q(c)          = Q2K              ! [kg/kg]
            forc_u(c)          = U_PHY(I,1,J)
            forc_v(c)          = V_PHY(I,1,J)
           !forc_rho(c)        = SFCPRS / (287.04 * SFCTMP * (1.0+ 0.61 * Q2K)) ![kg/m/m/m] 
            forc_lwrad(c)      = LWDN             ! [W/m/m]
            prec(c)            = PRCP             ! [mm/s]
           !sabg(c)            = SOLNET           ! Annotated by YW
            lat(c)             = XLAT_URB2D(I,J)*pie/180  ! [radian] 
           !--- YangW Updated: Read in xlong for calculating coszen ---!
            lon(c)             = xlong_lsm(I,J)*pie/180
            coszen_lsm_temp(c) = coszen_lsm(I,J)
            do_capsnow(c)      = .false.
            lakedepth(c)       = lakedepth2d(i,j)
            savedtke1(c)       = savedtke12d(i,j)
            snowdp(c)          = snowdp2d(i,j)
            h2osno(c)          = h2osno2d(i,j)
            snl(c)             = snl2d(i,j)
            t_grnd(c)          = t_grnd2d(i,j)
            do k = 1,nlevlake
               t_lake(c,k)        = t_lake3d(i,k,j)
               lake_icefrac(c,k)  = lake_icefrac3d(i,k,j)
               z_lake(c,k)        = z_lake3d(i,k,j)
               dz_lake(c,k)       = dz_lake3d(i,k,j)
            enddo
            do k = -nlevsnow+1,nlevsoil
               t_soisno(c,k)      = t_soisno3d(i,k,j)
               h2osoi_ice(c,k)    = h2osoi_ice3d(i,k,j)
               h2osoi_liq(c,k)    = h2osoi_liq3d(i,k,j)
               h2osoi_vol(c,k)    = h2osoi_vol3d(i,k,j)
               z(c,k)             = z3d(i,k,j)
               dz(c,k)            = dz3d(i,k,j)
            enddo   
            do k = -nlevsnow+0,nlevsoil
               zi(c,k)            = zi3d(i,k,j)
            enddo
            do k = 1,nlevsoil
               watsat(c,k)        = watsat3d(i,k,j)
               csol(c,k)          = csol3d(i,k,j)
               tkmg(c,k)          = tkmg3d(i,k,j)
               tkdry(c,k)         = tkdry3d(i,k,j)
               tksatu(c,k)        = tksatu3d(i,k,j)
            enddo
             
            !=== YangW Updated: Recalculate albedo and sabg ===!
                if (albedo_flag) then
                  !--- First: Calculate the cosine of solar zenith angle ---
                    calday    = float(idate(2))+float(idate(3))/86400.        ! Calculate the cosine of zenith angle
                    CALL orb_coszen(coszen_yw,calday,lon(c),lat(c)) 
                    if (coszen_flag) then
                        coszen    = max(0.000, coszen_yw) 
                    else
                        coszen    = max(0.000, coszen_lsm_temp(c)) ! YangW: I think this judgement is also essential For coszen_lsm
                    end if
                    if (lakedebug_yw) then
                        print *, "****** YangW Following are Lake Debug Processes ******"
                        print *, "*** Check Whether the calculated coszen_yw is equal to the input coszen_lsm_temp***"
                        print *, "calculated coszen_yw= ",coszen_yw,"; Input coszen_lsm_temp= ",coszen_lsm_temp(c)
                        print *, "******************************************************"
                    endif
!print 100, it,ir2,coszen,idate  
!100  format(2i8,f10.2,3i8)
                  !--- Second: Divide the Four solar radiation ---!
                    if(coszen <= 0) SOLDN=0   
                    !--- YangW: We may enlarge the SOLDN in Offline but not Coupled Model ---!
                    !SOLDN    = SOLDN*1.05         ! by Anning Huang
                    !--- Seperate the solar radiation into four kinds ---!
                    CALL rd_forcing(SOLDN,coszen,forc_sols,forc_soll,forc_solsd,forc_solld)    

                    h2osno_old(c) = h2osno(c)
                    !--- Give a initialized albedo & coresponding sabg ---! 
                    !--- Attention: After Initialization, all the albedo & sabg should be calculated after the updated h2osno by LakeMain ---!
                      if (albedo_indice < 2) then
                          fsno_yw = 0.         ! snow fraction
                          if (snowdp(c) > 0.) then
                              fsno_yw  = tanh(snowdp(c)/(2.5 * zlnd *h2osno(c)/snowdp(c)/100))
                          end if
                          if (snl(c)==0) then  ! No snow layer cover -> No snowage
                              sag_yw = 0.
                          end if 
                          !--- YangW: Calculate snowage when there exits h2osno although no snow layer ??---!
                          CALL snowage_bats(dtime,t_grnd(c),h2osno(c),h2osno_old(c),sag_yw)                      ! Get the sag_yw(snowage)
                          CALL alblake_yw(fsno_yw,sag_yw,coszen,h2osno(c),t_grnd(c),snl(c),alb_yw,lakedebug_yw)  ! Get the alb_yw(albedo) based on fsno_yw & sag_yw 
                          albedo_indice = albedo_indice+1
                      end if
                      if (lakedebug_yw) then
                          print *, "****** YangW Following are Lake Debug Processes ******"
                          print *, "*** Check Whether the albedo_indice will be re-declared during the Integration ***"
                          print *, "albedo_indice = ",albedo_indice
                          print *, "******************************************************"
                      endif

                      !--- Renew the sabg with 4-kind incoming solar radiation and corresponding albedo ---!
                      !--- YangW: Adopt the CoLM-Lake action: Renew sabg(c) only when forc_sols+forc_soll+forc_solsd+forc_solld>0. ---!
                      sabg(c)    = 0.  
                      if (forc_sols+forc_soll+forc_solsd+forc_solld>0.) then
                         sabg(c) =   forc_sols *(1.-alb_yw(1,1)) + forc_soll *(1.-alb_yw(2,1)) &
                                   + forc_solsd*(1.-alb_yw(1,2)) + forc_solld*(1.-alb_yw(2,2))
                         sabg(c) = sabg(c)*sabg_factor
                      endif
                 else
                      sabg(c) = SOLDN*(1.-ALBEDO(I,J)) ! Original albedo and sabg Scheme
                      sabg(c) = sabg(c)*sabg_factor
                 end if                                ! For albedo_flag
             !==================================================!
  
             !=== YangW: Output sabg_lakeout & albedo_lakeout ===!
              if (albedo_flag) then
                  albedo_sols_yw(I,J)    = alb_yw(1,1)
                  albedo_soll_yw(I,J)    = alb_yw(2,1)
                  albedo_solsd_yw(I,J)   = alb_yw(1,2)
                  albedo_solld_yw(I,J)   = alb_yw(2,2)
              end if
              sabg_lakeout(I,J)      = sabg(c)
              cosz_lakeout(I,J)      = coszen
             !===================================================!

          enddo                                  ! For c the Column Initialized Loop

          !--- Start Lake Physical Solutions ---!
           CALL LakeMain(forc_t,forc_pbot,forc_psrf,forc_hgt,forc_hgt_q, & !I  
                         forc_hgt_t,forc_hgt_u,forc_q, forc_u,           &
                         forc_v,forc_lwrad,prec, sabg,lat,               &
                         z_lake,dz_lake,lakedepth,do_capsnow,            &
                         h2osno,snowdp,snl,z,dz,zi,                      & !H
                         h2osoi_vol,h2osoi_liq,h2osoi_ice,               &
                         t_grnd,t_soisno,t_lake,                         &
                         savedtke1,lake_icefrac,                         &
                         eflx_lwrad_net,eflx_gnet,                       & !O 
                         eflx_sh_tot,eflx_lh_tot,                        &
                         t_ref2m,q_ref2m,                                &
                         taux,tauy,ram1,                                 &
                         !--- YangW: Variables For Updated Roughness & Eddy Diffusivity Calculation ---!
                          !--- Local Declaration Variables 
                          z0mg,z0hg,z0qg,ws,ks,u2m,drhodz,n2,ri,km,ke,ked,kw,           & 
                         !--- YangW : Updated External Variables From the namelist ---!
                          lakedebug_yw, beta_yw      ,eta_flag         , eta_scale_yw, eta_yw,        &
                          roughness_flag ,diffusivity_flag , diffusivity_index_yw, tdmax, mixing_factor, mixing_factor_ked)
         
          !--- Exchange In&Out Variables with land model ---!        
           do c = 1,column
              HFX(I,J)    = eflx_sh_tot(c)            ![W/m/m]
              LH(I,J)     = eflx_lh_tot(c)            !W/m/m]
              GRDFLX(I,J) = eflx_gnet(c)              !W/m/m]
              TSK(I,J)    = t_grnd(c)                 ![K]
              T2(I,J)     = t_ref2m(c)
              TH2(I,J)    = T2(I,J)*(1.E5/PSFC)**RCP
              Q2(I,J)     = q_ref2m(c) 
              albedo(i,j) = ( alb_ice * lake_icefrac(c,1) ) + ( (1.0-lake_icefrac(c,1)) * alb_wat)  
              if (tsk(i,j) >= tfrz) then
                 qfx(i,j) = eflx_lh_tot(c)/hvap
              else
                 qfx(i,j) = eflx_lh_tot(c)/hsub       ! heat flux (W/m^2)=>mass flux(kg/(sm^2))
              endif

              !--- YangW : Output roughness ---!
                z0mg_yw(i,j) = z0mg(1)
                z0hg_yw(i,j) = z0hg(1)
                z0qg_yw(i,j) = z0qg(1)

           enddo
          !--- Renew Lake State Varialbes:(14) For Next Loop ---!
           do c = 1,column
              savedtke12d(i,j)         = savedtke1(c)
              snowdp2d(i,j)            = snowdp(c)
              h2osno2d(i,j)            = h2osno(c)
              snl2d(i,j)               = snl(c)
              t_grnd2d(i,j)            = t_grnd(c)
              do k = 1,nlevlake
                 t_lake3d(i,k,j)       = t_lake(c,k)
                 lake_icefrac3d(i,k,j) = lake_icefrac(c,k)
                 kw_yw(i,k,j)          = kw(1,k)        ! YangW: Output Total diffusivity
              enddo
              do k = -nlevsnow+1,nlevsoil
                 z3d(i,k,j)            = z(c,k)
                 dz3d(i,k,j)           = dz(c,k) 
                 t_soisno3d(i,k,j)     = t_soisno(c,k)
                 h2osoi_liq3d(i,k,j)   = h2osoi_liq(c,k)
                 h2osoi_ice3d(i,k,j)   = h2osoi_ice(c,k)
                 h2osoi_vol3d(i,k,j)   = h2osoi_vol(c,k)
              enddo
              do k = -nlevsnow+0,nlevsoil
                 zi3d(i,k,j)           = zi(c,k)
              enddo

              !--- YangW : Output extra diffusivity ---!
                km_yw(i,j)               = km
                do k = 1,nlevlake-1
                   ke_yw(i,k,j)          = ke(k)
                   ked_yw(i,k,j)         = ked(k)
                enddo
              !----------------------------------------!
                
           enddo

          !=== As the h2osno updated, Renew the alb_yw for next integration time step ===!
          !=== YangW Attention : Without the column loop, the initialized c=2 and all the input var(c) is zero, must be wrong ===!
           if (albedo_flag) then
               do c = 1, column
                  fsno_yw = 0.         ! snow fraction
                  if (snowdp(c) > 0.) then
                      fsno_yw  = tanh(snowdp(c)/(2.5 * zlnd *h2osno(c)/snowdp(c)/100))
                  end if
                  if (snl(c)==0) then  ! No snow layer cover -> No snowage
                      sag_yw = 0.
                  end if
                 !--- YangW: Calculate snowage when there exits h2osno although no snow layer ??---!
                  CALL snowage_bats(dtime,t_grnd(c),h2osno(c),h2osno_old(c),sag_yw)                      ! Get the sag_yw(snowage)
                  CALL alblake_yw(fsno_yw,sag_yw,coszen,h2osno(c),t_grnd(c),snl(c),alb_yw,lakedebug_yw)  ! Get the alb_yw(albedo) based on fsno_yw & sag_yw 
               enddo
            end if
          !===============================================================================!

        endif ! For lakemask(i,j) if: Run LakeModel on LakeGrids 

       ENDDO  ! For i,j loop
    ENDDO

    !=== YangW Updated: Prepare idate For the Next integration step ===!
      CALL ticktime (dtime, idate)                       ! Get the idate(3) at next integration step
    !==================================================================!

END SUBROUTINE Lake


    SUBROUTINE LakeMain(forc_t,forc_pbot,forc_psrf,forc_hgt,forc_hgt_q, & !I  
                          forc_hgt_t,forc_hgt_u,forc_q, forc_u,         &   
                          forc_v,forc_lwrad,prec, sabg,lat,             &   
                          z_lake,dz_lake,lakedepth,do_capsnow,          &
                          h2osno,snowdp,snl,z,dz,zi,                    & !H
                          h2osoi_vol,h2osoi_liq,h2osoi_ice,             &
                          t_grnd,t_soisno,t_lake,                       &  
                          savedtke1,lake_icefrac,                       &
                          eflx_lwrad_net,eflx_gnet,                     & !O 
                          eflx_sh_tot,eflx_lh_tot,                      &
                          t_ref2m,q_ref2m,                              &
                          taux,tauy,ram1,                               &
                        !--- YangW: Updated Roughness & Eddy Diffusivity Calculation ---!
                        !--- Following are Local Variabls ---!
                         z0mg,z0hg,z0qg,ws,ks,u2m,drhodz,n2,ri,km,ke,ked,kw,          & 
                        !--- YangW : Updated External Variables From the namelist ---!
                          lakedebug_yw, beta_yw, eta_flag,  eta_scale_yw, eta_yw,       &
                          roughness_flag, diffusivity_flag, diffusivity_index_yw, tdmax,  mixing_factor, mixing_factor_ked)
    implicit none
!in: 
    real(r8),intent(in) :: forc_t(1)           ! atmospheric temperature (Kelvin)
    real(r8),intent(in) :: forc_pbot(1)        ! atm bottom level pressure (Pa) 
    real(r8),intent(in) :: forc_psrf(1)        ! atmospheric surface pressure (Pa)
    real(r8),intent(in) :: forc_hgt(1)         ! atmospheric reference height (m)
    real(r8),intent(in) :: forc_hgt_q(1)       ! observational height of humidity [m]
    real(r8),intent(in) :: forc_hgt_t(1)       ! observational height of temperature [m]
    real(r8),intent(in) :: forc_hgt_u(1)       ! observational height of wind [m]
    real(r8),intent(in) :: forc_q(1)           ! atmospheric specific humidity (kg/kg)
    real(r8),intent(in) :: forc_u(1)           ! atmospheric wind speed in east direction (m/s)
    real(r8),intent(in) :: forc_v(1)           ! atmospheric wind speed in north direction (m/s)
   !real(r8),intent(in) :: forc_rho(1)         ! density (kg/m**3)
    real(r8),intent(in) :: forc_lwrad(1)       ! downward infrared (longwave) radiation (W/m**2)
    real(r8),intent(in) :: prec(1)             ! snow or rain rate [mm/s]
    real(r8),intent(in) :: sabg(1)             ! solar radiation absorbed by ground (W/m**2)
    real(r8),intent(in) :: lat(1)              ! latitude (radians)
    real(r8),intent(in) :: z_lake(1,nlevlake)  ! layer depth for lake (m)
    real(r8),intent(in) :: dz_lake(1,nlevlake) ! layer thickness for lake (m)
    real(r8),intent(in) :: lakedepth(1)        ! column lake depth (m)
    logical,intent(in)  :: do_capsnow(1)       ! true => do snow capping
    !!!!!!!!!!!!!!!!tep(in),hydro(in)   
    !real(r8), intent(in) :: watsat(1,1:nlevsoil)      ! volumetric soil water at saturation (porosity)
    !!!!!!!!!!!!!!!!hydro

!in&out
    real(r8),intent(inout) :: h2osoi_vol(1,-nlevsnow+1:nlevsoil)  ! volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
    real(r8),intent(inout) :: t_grnd(1)                           ! ground temperature (Kelvin)
    real(r8),intent(inout) :: h2osno(1)                           ! snow water (mm H2O)
    real(r8),intent(inout) :: snowdp(1)                           ! snow height (m)
    real(r8),intent(inout) :: z(1,-nlevsnow+1:nlevsoil)           ! layer depth for snow & soil (m)
    real(r8),intent(inout) :: dz(1,-nlevsnow+1:nlevsoil)          ! layer thickness for soil or snow (m)
    real(r8),intent(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)    ! soil (or snow) temperature (Kelvin)
    real(r8),intent(inout) :: t_lake(1,nlevlake)                  ! lake temperature (Kelvin)
    integer ,intent(inout) :: snl(1)                              ! number of snow layers
    real(r8),intent(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)  ! liquid water (kg/m2)
    real(r8),intent(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)  ! ice lens (kg/m2)
    real(r8),intent(inout) :: savedtke1(1)                        ! top level eddy conductivity from previous timestep (W/m.K)
    real(r8),intent(inout) :: zi(1,-nlevsnow+0:nlevsoil)          ! interface level below a "z" level (m)
    real(r8),intent(inout) :: lake_icefrac(1,nlevlake)            ! mass fraction of lake layer that is frozen

!out:
    real(r8),intent(out) :: eflx_gnet(1)       ! net heat flux into ground (W/m**2)
    real(r8),intent(out) :: eflx_lwrad_net(1)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8),intent(out) :: eflx_sh_tot(1)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8),intent(out) :: eflx_lh_tot(1)     ! total latent heat flux (W/m8*2)  [+ to atm]
    real(r8),intent(out) :: t_ref2m(1)         ! 2 m height surface air temperature (Kelvin)
    real(r8),intent(out) :: q_ref2m(1)         ! 2 m height surface specific humidity (kg/kg)
    real(r8),intent(out) :: taux(1)            ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8),intent(out) :: tauy(1)            ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8),intent(out) :: ram1(1)            ! aerodynamical resistance (s/m)
                                               ! for calculation of decay of eddy diffusivity with depth
                                               ! Change the type variable to pass back to WRF.

!local output
    real(r8) :: begwb(1)              ! water mass begining of the time step
    real(r8) :: t_veg(1)              ! vegetation temperature (Kelvin)
    real(r8) :: eflx_soil_grnd(1)     ! soil heat flux (W/m**2) [+ = into soil]
    real(r8) :: eflx_lh_grnd(1)       ! ground evaporation heat flux (W/m**2) [+ to atm]
    real(r8) :: eflx_sh_grnd(1)       ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8) :: eflx_lwrad_out(1)     ! emitted infrared (longwave) radiation (W/m**2)
    real(r8) :: qflx_evap_tot(1)      ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8) :: qflx_evap_soi(1)      ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8) :: qflx_prec_grnd(1)     ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8) :: forc_snow(1)          ! snow rate [mm/s]
    real(r8) :: forc_rain(1)          ! rain rate [mm/s]
    real(r8) :: qflx_snomelt(1)       !snow melt (mm H2O /s) tem(out),snowwater(in)
    integer  :: imelt(1,-nlevsnow+1:nlevsoil)      ! flag for melting (=1), freezing (=2), Not=0 (new)
    real(r8) :: endwb(1)              ! water mass end of the time step
    real(r8) :: snowage(1)            ! non dimensional snow age [-]
    real(r8) :: snowice(1)            ! average snow ice lens
    real(r8) :: snowliq(1)            ! average snow liquid water
    real(r8) :: t_snow(1)             ! vertically averaged snow temperature
    real(r8) :: qflx_drain(1)         ! sub-surface runoff (mm H2O /s)
    real(r8) :: qflx_surf(1)          ! surface runoff (mm H2O /s)
    real(r8) :: qflx_infl(1)          ! infiltration (mm H2O /s)
    real(r8) :: qflx_qrgwl(1)         ! qflx_surf at glaciers, wetlands, lakes
    real(r8) :: qcharge(1)            ! aquifer recharge rate (mm/s)
    real(r8) :: qflx_snowcap(1)       ! excess precipitation due to snow capping (mm H2O /s) [+]
    real(r8) :: qflx_snowcap_col(1)   ! excess precipitation due to snow capping (mm H2O /s) [+]
    real(r8) :: qflx_snow_grnd_pft(1) ! snow on ground after interception (mm H2O/s) [+]
    real(r8) :: qflx_snow_grnd_col(1) ! snow on ground after interception (mm H2O/s) [+]
    real(r8) :: qflx_rain_grnd(1)     ! rain on ground after interception (mm H2O/s) [+]
    real(r8) :: frac_iceold(1,-nlevsnow+1:nlevsoil)  ! fraction of ice relative to the tot water
    real(r8) :: qflx_evap_tot_col(1)  ! pft quantity averaged to the column (assuming one pft)
    real(r8) :: soilalpha(1)          ! factor that reduces ground saturated specific humidity (-)
    real(r8) :: zwt(1)                ! water table depth
    real(r8) :: fcov(1)               ! fractional area with water table at surface
    real(r8) :: rootr_column(1,1:nlevsoil)           ! effective fraction of roots in each soil layer
    real(r8) :: qflx_evap_grnd(1)     ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8) :: qflx_sub_snow(1)      ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8) :: qflx_dew_snow(1)      ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8) :: qflx_dew_grnd(1)      ! ground surface dew formation (mm H2O /s) [+]
    real(r8) :: qflx_rain_grnd_col(1) ! rain on ground after interception (mm H2O/s) [+]

    !=== Variables For Updated Physical Part ===!
     !--- YangW: Updated Roughness & Eddy Diffusivity Calculation ---!
     !--- Following are Local Variabls Delivered into Subroutine ShaLakeFluxes & ShaLakeTemperature---!
       real(r8),intent(out) :: z0mg(1)            ! roughness length over ground, momentum      (m)
       real(r8),intent(out) :: z0hg(1)            ! roughness length over ground, sensible heat (m)
       real(r8),intent(out) :: z0qg(1)            ! roughness length over ground, latent heat   (m)
       real(r8),intent(out) :: ws(1)              ! surface friction velocity (m/s)
       real(r8),intent(out) :: ks(1)              ! coefficient passed to ShalLakeTemperature
       real(r8),intent(out) :: u2m                ! 2 m wind speed (m/s)

       real(r8) :: drhodz(nlevlake-1)             ! d [rhow] /dz (kg/m**4)
       real(r8) :: n2(nlevlake-1)                 ! brunt-vaisala frequency (/s**2)
       real(r8) :: ri(nlevlake-1)                 ! richardson number
       real(r8),intent(out) :: km                   ! molecular diffusion coefficient (m**2/s)
       real(r8),intent(out) :: ke(nlevlake-1)       ! eddy diffusion coefficient (m**2/s)
       real(r8),intent(out) :: ked(nlevlake-1)      ! enhanced diffusion coefficient (m**2/s)
       real(r8),intent(out) :: kw(lbc:ubc,nlevlake) ! molecular + eddy diffusion coefficient (m**2/s)
     !---------------------------------------------------!
    !===========================================!

  !=== Some Variables Prepared For Coupled WRF-Lake namelist -> Subroutine Lake ===!
    logical, intent(in) :: lakedebug_yw
    real   , intent(in) :: beta_yw   

    !--- For Subroutine Lake -> LakeMain -> ShalLakeFluxes ---!
    logical, intent(in) :: roughness_flag    

    !--- For Subroutine Lake -> LakeMain -> ShalLakeTemperatures ---!
    logical, intent(in) :: diffusivity_flag       ! True For using the New diffusivity scheme
    integer, intent(in) :: diffusivity_index_yw   ! 1 for kw=mixing_factor*(km+ke+ked)
                                                  ! 2 for kw=km+mixing_factor*ke+ked
    real, intent(in)    :: tdmax
    real, intent(in)    :: mixing_factor    
    real, intent(in)    :: mixing_factor_ked    
    logical, intent(in) :: eta_flag
    real, intent(in)    :: eta_scale_yw
    real, intent(in)    :: eta_yw
  !================================================================================!


    !--- Judge Rainfall-kind of Snow-kind precipitation with the critical temperature tfrz+tcrit(2.5K) ---!
    if (prec(1)> 0.) then
        if ( forc_t(1) > (tfrz + tcrit)) then     
            forc_rain(1) = prec(1)
            forc_snow(1) = 0.
           !flfall(1) = 1.
         else
            forc_rain(1) = 0.
            forc_snow(1) = prec(1)
           !if ( forc_t(1) <= tfrz) then      ! Not Implemented ?? YangW
           !    flfall(1) = 0.
           !else if ( forc_t(1) <= tfrz+2.) then
           !         flfall(1) = -54.632 + 0.2 *  forc_t(1)
           !     else
           !         flfall(1) = 0.4
         endif
    else
         forc_rain(1) = 0.
         forc_snow(1) = 0.
        !flfall(1) = 1.
    endif

    CALL ShalLakeFluxes(forc_t,forc_pbot,forc_psrf,forc_hgt,forc_hgt_q,   & !i
                          forc_hgt_t,forc_hgt_u,forc_q,                   &
                          forc_u,forc_v,forc_lwrad,forc_snow,             &
                          forc_rain,t_grnd,h2osno,snowdp,sabg,lat,        &
                          dz,dz_lake,t_soisno,t_lake,snl,h2osoi_liq,      &
                          h2osoi_ice,savedtke1,                           &
                          qflx_prec_grnd,qflx_evap_soi,qflx_evap_tot,     & !o
                          eflx_sh_grnd,eflx_lwrad_out,eflx_lwrad_net,     &
                          eflx_soil_grnd,eflx_sh_tot,eflx_lh_tot,         &
                          eflx_lh_grnd,t_veg,t_ref2m,q_ref2m,taux,tauy,   &
                          ram1,eflx_gnet,                                 &
                        !--- YangW: Updated Roughness Calculation---!
                        !--- Following are Local Variabls ---!
                          u2m,ws,ks,z0mg,z0hg,z0qg,z_lake,lake_icefrac,lakedepth,   &
                        !--- YangW : Updated External Variables From the namelist ---!
                          lakedebug_yw, beta_yw      ,eta_flag         , eta_scale_yw, eta_yw,         &
                          roughness_flag , diffusivity_flag, diffusivity_index_yw, tdmax, mixing_factor, mixing_factor_ked)

    CALL ShalLakeTemperature(t_grnd,h2osno,sabg,dz,dz_lake,z,zi,         & !i
                             z_lake,ws,ks,snl,eflx_gnet,lakedepth,       &
                             lake_icefrac,snowdp,                        & !i&o
                             eflx_sh_grnd,eflx_sh_tot,eflx_soil_grnd,    & !o
                             t_lake,t_soisno,h2osoi_liq,                 &
                             h2osoi_ice,savedtke1,                       &
                             frac_iceold,qflx_snomelt,imelt,             & 
                        !--- YangW: Updated Eddy Diffusivity Calculation---!
                        !--- Following are Local Variabls ---!
                             drhodz,n2,ri,km,ke,ked,kw,                      &
                        !--- YangW : Updated External Variables From the namelist ---!
                          lakedebug_yw, beta_yw      ,eta_flag         , eta_scale_yw, eta_yw,        &
                          roughness_flag , diffusivity_flag, diffusivity_index_yw, tdmax, mixing_factor, mixing_factor_ked)

    CALL ShalLakeHydrology(dz_lake,forc_rain,forc_snow,                      & ! i
                           begwb,qflx_evap_tot,forc_t,do_capsnow,            &
                           t_grnd,qflx_evap_soi,                             &
                           qflx_snomelt,imelt,frac_iceold,                   & ! i add by guhp
                           z,dz,zi,snl,h2osno,snowdp,lake_icefrac,t_lake,    & ! i&o
                           endwb,snowage,snowice,snowliq,t_snow,             & ! o
                           t_soisno,h2osoi_ice,h2osoi_liq,h2osoi_vol,        &
                           qflx_drain,qflx_surf,qflx_infl,qflx_qrgwl,        &
                           qcharge,qflx_prec_grnd,qflx_snowcap,              &
                           qflx_snowcap_col,qflx_snow_grnd_pft,              &
                           qflx_snow_grnd_col,qflx_rain_grnd,                &
                           qflx_evap_tot_col,soilalpha,zwt,fcov,             &
                           rootr_column,qflx_evap_grnd,qflx_sub_snow,        &
                           qflx_dew_snow,qflx_dew_grnd,qflx_rain_grnd_col)
!==================================================================================
! !DESCRIPTION:
! Calculation of Shallow Lake Hydrology. Full hydrology of snow layers is
! done. However, there is no infiltration, and the water budget is balanced with 
                       
   END SUBROUTINE LakeMain


SUBROUTINE ShalLakeFluxes(forc_t,forc_pbot,forc_psrf,forc_hgt,forc_hgt_q, &  !i
                          forc_hgt_t,forc_hgt_u,forc_q,                   &
                          forc_u,forc_v,forc_lwrad,forc_snow,             &
                          forc_rain,t_grnd,h2osno,snowdp,sabg,lat,        &
                          dz,dz_lake,t_soisno,t_lake,snl,h2osoi_liq,      &
                          h2osoi_ice,savedtke1,                           &
                          qflx_prec_grnd,qflx_evap_soi,qflx_evap_tot,     &  !o
                          eflx_sh_grnd,eflx_lwrad_out,eflx_lwrad_net,     &
                          eflx_soil_grnd,eflx_sh_tot,eflx_lh_tot,         &
                          eflx_lh_grnd,t_veg,t_ref2m,q_ref2m,taux,tauy,   &
                          ram1,eflx_gnet,                                 &           
                        !--- YangW: Updated Roughness Calculation---!
                        !--- Following are Local Variabls ---!
                          u2m,ws,ks,z0mg,z0hg,z0qg,z_lake,lake_icefrac,lakedepth,   &
                        !--- YangW : Updated External Variables From the namelist ---!
                          lakedebug_yw, beta_yw      ,eta_flag         , eta_scale_yw, eta_yw,         &
                          roughness_flag , diffusivity_flag , diffusivity_index_yw, tdmax, mixing_factor, mixing_factor_ked)
!==============================================================================
! DESCRIPTION:
! Calculates lake temperatures and surface fluxes for shallow lakes.
!
! Shallow lakes have variable depth, possible snow layers above, freezing & thawing of lake water,
! and soil layers with active temperature and gas diffusion below.
!
! WARNING: This subroutine assumes lake columns have one and only one pft.
!
! REVISION HISTORY:
! Created by Zack Subin, 2009
! Reedited by Hongping Gu, 2010 
!==============================================================================
    implicit none
!in: 
    real(r8),intent(in) :: forc_t(1)          ! atmospheric temperature (Kelvin)
    real(r8),intent(in) :: forc_pbot(1)       ! atmospheric pressure (Pa)
    real(r8),intent(in) :: forc_psrf(1)       ! atmospheric surface pressure (Pa)
    real(r8),intent(in) :: forc_hgt(1)        ! atmospheric reference height (m)
    real(r8),intent(in) :: forc_hgt_q(1)      ! observational height of humidity [m]
    real(r8),intent(in) :: forc_hgt_t(1)      ! observational height of temperature [m]
    real(r8),intent(in) :: forc_hgt_u(1)      ! observational height of wind [m]
    real(r8),intent(in) :: forc_q(1)          ! atmospheric specific humidity (kg/kg)
    real(r8),intent(in) :: forc_u(1)          ! atmospheric wind speed in east direction (m/s)
    real(r8),intent(in) :: forc_v(1)          ! atmospheric wind speed in north direction (m/s)
    real(r8),intent(in) :: forc_lwrad(1)      ! downward infrared (longwave) radiation (W/m**2)
   !real(r8),intent(in) :: forc_rho(1)        ! density (kg/m**3)
    real(r8),intent(in) :: forc_snow(1)       ! snow rate [mm/s]
    real(r8),intent(in) :: forc_rain(1)       ! rain rate [mm/s]
    real(r8),intent(in) :: h2osno(1)          ! snow water (mm H2O)
    real(r8),intent(in) :: snowdp(1)          ! snow height (m)
    real(r8),intent(in) :: sabg(1)            ! solar radiation absorbed by ground (W/m**2)
    real(r8),intent(in) :: lat(1)             ! latitude (radians)
    real(r8),intent(in) :: dz(1,-nlevsnow+1:nlevsoil)         ! layer thickness for soil or snow (m)
    real(r8),intent(in) :: dz_lake(1,nlevlake)                ! layer thickness for lake (m)
    real(r8),intent(in) :: t_soisno(1,-nlevsnow+1:nlevsoil)   ! soil (or snow) temperature (Kelvin)
    real(r8),intent(in) :: t_lake(1,nlevlake)                 ! lake temperature (Kelvin)
    integer ,intent(in) :: snl(1)                             ! number of snow layers
    real(r8),intent(in) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil) ! liquid water (kg/m2)
    real(r8),intent(in) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil) ! ice lens (kg/m2)
    real(r8),intent(in) :: savedtke1(1)       ! top level eddy conductivity from previous timestep (W/m.K)
!inout:
    real(r8),intent(inout) :: t_grnd(1)       ! ground temperature (Kelvin)
!out:
    real(r8),intent(out):: qflx_prec_grnd(1)  ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8),intent(out):: qflx_evap_soi(1)   ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8),intent(out):: qflx_evap_tot(1)   ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8),intent(out):: eflx_sh_grnd(1)    ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8),intent(out):: eflx_lwrad_out(1)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8),intent(out):: eflx_lwrad_net(1)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8),intent(out):: eflx_soil_grnd(1)  ! soil heat flux (W/m**2) [+ = into soil]
    real(r8),intent(out):: eflx_sh_tot(1)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8),intent(out):: eflx_lh_tot(1)     ! total latent heat flux (W/m8*2)  [+ to atm]
    real(r8),intent(out):: eflx_lh_grnd(1)    ! ground evaporation heat flux (W/m**2) [+ to atm]
    real(r8),intent(out):: t_veg(1)           ! vegetation temperature (Kelvin)
    real(r8),intent(out):: t_ref2m(1)         ! 2 m height surface air temperature (Kelvin)
    real(r8),intent(out):: q_ref2m(1)         ! 2 m height surface specific humidity (kg/kg)
    real(r8),intent(out):: taux(1)            ! wind (shear) stress: e-w (kg/m/s**2)
    real(r8),intent(out):: tauy(1)            ! wind (shear) stress: n-s (kg/m/s**2)
    real(r8),intent(out):: ram1(1)            ! aerodynamical resistance (s/m)
    real(r8),intent(out):: eflx_gnet(1)       !net heat flux into ground (W/m**2)
                                              ! Change the type variable to pass back to WRF.
!OTHER LOCAL VARIABLES:
    integer , parameter :: islak  = 2         ! index of lake, 1 = deep lake, 2 = shallow lake
    integer , parameter :: niters = 3         ! maximum number of iterations for surface temperature
    real(r8), parameter :: beta1  = 1._r8     ! coefficient of convective velocity (in computing W_*) [-]
    real(r8), parameter :: emg    = 0.97_r8   ! ground emissivity (0.97 for snow)
    real(r8), parameter :: zii    = 1000._r8  ! convective boundary height [m]
    real(r8) :: forc_th(1)              ! atmospheric potential temperature (Kelvin)
    real(r8) :: forc_vp(1)              ! atmospheric vapor pressure (Pa)
    real(r8) :: forc_rho(1)             ! density (kg/m**3)
    integer  :: i,fc,fp,g,c,p           ! do loop or array index
    integer  :: fncopy                  ! number of values in pft filter copy
    integer  :: fnold                   ! previous number of pft filter values
    integer  :: fpcopy(num_shlakep)     ! pft filter copy for iteration loop
    integer  :: iter                    ! iteration index
    integer  :: nmozsgn(lbp:ubp)        ! number of times moz changes sign
    integer  :: jtop(lbc:ubc)           ! top level for each column (no longer all 1)
   !real(r8) :: dtime                   ! land model time step (sec)
    real(r8) :: ax                      ! used in iteration loop for calculating t_grnd (numerator of NR solution)
    real(r8) :: bx                      ! used in iteration loop for calculating t_grnd (denomin. of NR solution)
    real(r8) :: degdT                   ! d(eg)/dT
    real(r8) :: dqh(lbp:ubp)            ! diff of humidity between ref. height and surface
    real(r8) :: dth(lbp:ubp)            ! diff of virtual temp. between ref. height and surface
    real(r8) :: dthv                    ! diff of vir. poten. temp. between ref. height and surface
    real(r8) :: dzsur(lbc:ubc)          ! 1/2 the top layer thickness (m)
    real(r8) :: eg                      ! water vapor pressure at temperature T [pa]
    real(r8) :: htvp(lbc:ubc)           ! latent heat of vapor of water (or sublimation) [j/kg]
    real(r8) :: obu(lbp:ubp)            ! monin-obukhov length (m)
    real(r8) :: obuold(lbp:ubp)         ! monin-obukhov length of previous iteration
    real(r8) :: qsatg(lbc:ubc)          ! saturated humidity [kg/kg]
    real(r8) :: qsatgdT(lbc:ubc)        ! d(qsatg)/dT
    real(r8) :: qstar                   ! moisture scaling parameter
    real(r8) :: ram(lbp:ubp)            ! aerodynamical resistance [s/m]
    real(r8) :: rah(lbp:ubp)            ! thermal resistance [s/m]
    real(r8) :: raw(lbp:ubp)            ! moisture resistance [s/m]
    real(r8) :: stftg3(lbp:ubp)         ! derivative of fluxes w.r.t ground temperature
    real(r8) :: temp1(lbp:ubp)          ! relation for potential temperature profile
    real(r8) :: temp12m(lbp:ubp)        ! relation for potential temperature profile applied at 2-m
    real(r8) :: temp2(lbp:ubp)          ! relation for specific humidity profile
    real(r8) :: temp22m(lbp:ubp)        ! relation for specific humidity profile applied at 2-m
    real(r8) :: tgbef(lbc:ubc)          ! initial ground temperature
    real(r8) :: thm(lbc:ubc)            ! intermediate variable (forc_t+0.0098*forc_hgt_t)
    real(r8) :: thv(lbc:ubc)            ! virtual potential temperature (kelvin)
    real(r8) :: thvstar                 ! virtual potential temperature scaling parameter
    real(r8) :: tksur                   ! thermal conductivity of snow/soil (w/m/kelvin)
    real(r8) :: tsur                    ! top layer temperature
    real(r8) :: tstar                   ! temperature scaling parameter
    real(r8) :: um(lbp:ubp)             ! wind speed including the stablity effect [m/s]
    real(r8) :: ur(lbp:ubp)             ! wind speed at reference height [m/s]
    real(r8) :: ustar(lbp:ubp)          ! friction velocity [m/s]
    real(r8) :: wc                      ! convective velocity [m/s]
    real(r8) :: zeta                    ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: zldis(lbp:ubp)          ! reference height "minus" zero displacement height [m]
    real(r8) :: displa(lbp:ubp)         ! displacement (always zero) [m]
   !real(r8) :: z0mg(lbp:ubp)           ! roughness length over ground, momentum [m]
   !real(r8) :: z0hg(lbp:ubp)           ! roughness length over ground, sensible heat [m]
   !real(r8) :: z0qg(lbp:ubp)           ! roughness length over ground, latent heat [m]
    real(r8) :: u10(1)                  ! 10-m wind (m/s) (for dust model)
    real(r8) :: fv(1)                   ! friction velocity (m/s) (for dust model)

    real(r8) :: fm(lbp:ubp)             ! needed for BGC only to diagnose 10m wind speed
    real(r8) :: bw                      ! partial density of water (ice + liquid)
    real(r8) :: t_grnd_temp             ! Used in surface flux correction over frozen ground
    real(r8) :: betaprime(lbc:ubc)      ! Effective beta: 1 for snow layers, beta(islak) otherwise
    character*256 :: message 
   ! This assumes all radiation is absorbed in the top snow layer and will need
   ! to be changed for CLM 4.

 !=== Variables For Updated Physical Part ===!
  !--- YangW: Updated Roughness Calculation ---!
    real(r8)             :: u2m                      ! 2 m wind speed (m/s)
    real(r8),intent(out) :: ws(1)                    ! surface friction velocity (m/s)
    real(r8),intent(out) :: ks(1)                    ! coefficient passed to ShalLakeTemperature
                                                     ! for calculation of decay of eddy diffusivity with depth
    real(r8),intent(out) :: z0mg(1)                  ! roughness length over ground, momentum      [m]
    real(r8),intent(out) :: z0hg(1)                  ! roughness length over ground, sensible heat [m]
    real(r8),intent(out) :: z0qg(1)                  ! roughness length over ground, latent heat   [m]
    real(r8),intent(in)  :: lakedepth(1)             ! column lake depth        [m]       
    real(r8),intent(in)  :: z_lake(1,nlevlake)       ! layer thickness for lake [m]         
    real(r8),intent(in)  :: lake_icefrac(1,nlevlake) ! lake ice fraction             
    !--- Local Temporary Pre-Calculated Variables For Subroutine roughness_lake ---
      integer              :: ihhh                             
      real(r8)             :: visa                  ! kinematic viscosity of dry air[m2/s]
      real(r8)             :: fetch(lbp:ubp)        ! lake fetch = 25*lakedepth [m] 
      real(r8)             :: cur(lbp:ubp)          ! Charnock parameter (-)  
      real(r8) , parameter :: cur0  = 0.01          ! min. Charnock parameter    0.01
      real(r8) , parameter :: curm  = 0.11          ! maximum Charnock parameter 0.11
      real(r8) , parameter :: fcrit = 22.           ! critical dimensionless fetch for Charnock parameter (Vickers & Mahrt 1997, orignial fcrit=22)
    !------------------------------------------------------------------------------
 !===========================================!

  !=== Some Variables Prepared For Coupled WRF-Lake namelist -> Subroutine Lake ===!
    logical, intent(in) :: lakedebug_yw
    real(r8)            :: beta(2)                              ! fraction solar rad absorbed at surface: depends on lake type
    real   , intent(in) :: beta_yw   

    !--- For Subroutine Lake -> LakeMain -> ShalLakeFluxes ---!
    logical, intent(in) :: roughness_flag    

    !--- For Subroutine Lake -> LakeMain -> ShalLakeTemperatures ---!
    logical, intent(in) :: diffusivity_flag       ! True For using the New diffusivity scheme
    integer, intent(in) :: diffusivity_index_yw   ! 1 for kw=mixing_factor*(km+ke+ked)
                                                  ! 2 for kw=km+mixing_factor*ke+ked
    real, intent(in)    :: tdmax
    real, intent(in)    :: mixing_factor    
    real, intent(in)    :: mixing_factor_ked    
    logical, intent(in) :: eta_flag
    real, intent(in)    :: eta_scale_yw
    real, intent(in)    :: eta_yw

    beta(:) = beta_yw
  !================================================================================!
  
! Constants for lake temperature model
!
    !---beta is the partial energy absorbed at the lake surface if no snow.
!    data beta/0.4_r8, 0.4_r8/  ! (deep lake, shallow lake)
!    data za  /0.6_r8, 0.5_r8/
!    data eta /0.1_r8, 0.5_r8/
!-----------------------------------------------------------------------
! dtime = get_step_size()
! Begin calculations

! Here Original Concurent
    forc_th(1)  = forc_t(1) * (forc_psrf(1)/ forc_pbot(1))**(rair/cpair)
    forc_vp(1)  = forc_q(1) * forc_pbot(1)/ (0.622 + 0.378 * forc_q(1))
    forc_rho(1) = (forc_pbot(1) - 0.378 * forc_vp(1)) / (rair * forc_t(1))

    do fc = 1, num_shlakec
       c = filter_shlakec(fc)
       g = cgridcell(c)

       ! Surface temperature and fluxes
       ! Find top layer
       if (snl(c) > 0 .or. snl(c) < -5) then
           if (lakedebug_yw) then            ! YangW: This part is designed in the original model
               WRITE(message,*)  'snl is not defined in ShalLakeFluxesMod'
               print *, message
           end if
          !CALL wrf_message(message)
          !CALL wrf_error_fatal("snl: out of range value")
       end if
!       if (snl(c) /= 0) then
!           write(6,*)'snl is not equal to zero in ShalLakeFluxesMod'
!           call endrun()
!       end if
       jtop(c) = snl(c) + 1

       if (snl(c) < 0) then
           betaprime(c) = 1._r8  !Assume all solar rad. absorbed at the surface of the top snow layer. 
           dzsur(c) = dz(c,jtop(c))/2._r8
       else
           betaprime(c) = beta(islak)
           dzsur(c) = dz_lake(c,1)/2._r8
       end if
       ! Originally this was 1*dz, but shouldn't it be 1/2?

       ! Saturated vapor pressure, specific humidity and their derivatives
       ! at lake surface
       CALL QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg(c), qsatgdT(c))

       ! Potential, virtual potential temperature, and wind speed at the
       ! reference height
       thm(c) = forc_t(g) + 0.0098_r8*forc_hgt_t(g)    ! intermediate variable
       thv(c) = forc_th(g)*(1._r8+0.61_r8*forc_q(g))   ! virtual potential T
    end do

! Here Original Concurent
    do fp = 1, num_shlakep
       p = filter_shlakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       nmozsgn(p) = 0
       obuold(p)  = 0._r8
       displa(p)  = 0._r8

      !Latent heat
       if (t_grnd(c) > tfrz) then
          htvp(c) = hvap
       else
          htvp(c) = hsub
       end if

       ! Zack Subin, 3/26/09: Shouldn't this be the ground temperature rather than the air temperature above?
       ! I'll change it for now.
       ! Initialize stability variables For Subroutine MoninObukIni to calcualte MO length and wind speed
       ur(p)    = max(1.0_r8,sqrt(forc_u(g)*forc_u(g)+forc_v(g)*forc_v(g)))  ! wind speed at reference height [m/s] 
       dth(p)   = thm(c)-t_grnd(c)                                           ! diff of virtual temp. between ref. height and surface
       dqh(p)   = forc_q(g)-qsatg(c)                                         ! diff of humidity between ref. height and surface
       dthv     = dth(p)*(1._r8+0.61_r8*forc_q(g))+0.61_r8*forc_th(g)*dqh(p) ! diff of vir. poten. temp. between ref. height and surface
       zldis(p) = forc_hgt_u(g) - 0._r8                                      ! reference height "minus" zero displacement height [m]

      !=== YangW: Updated Roughness Calculation ===!
        if (roughness_flag) then
           !=== Subin Roughness Lengths ===!
            ustar(p) = 0.06  
            wc       = 0.5   
            fetch(p) = 25*lakedepth(c) 
            !--- Charnock number---!
             cur(p) = cur0 + curm * exp(max(-(fetch(p)*grav/ur(p)/ur(p))**(1./3.)/fcrit, &  ! Fetch-limited	  
                                            -(z_lake(c,nlevlake)*grav)**0.5/ur(p)))         ! depth-limited
            !--- kinematic viscosity ---!
             visa = 1.326e-5*(1.+6.542e-3*(forc_t(g)-tfrz) &
                    + 8.301e-6*(forc_t(g)-tfrz)**2 - 4.84e-9*(forc_t(g)-tfrz)**3)
            !--- Initialized for roughness_lake ---!
            !--- Actually this part will be redone in MoninObykIni ---!
            !--- Why needs to initialize the z0mg ??? YangW---!
             if(dthv.ge.0.) then
                um(p)=max(ur(p),0.1)
             else
                um(p)=sqrt(ur(p)*ur(p)+wc*wc)
             endif
             do ihhh=1,5
                z0mg(p)  = 0.013*ustar(p)*ustar(p)/grav+0.11*visa/ustar(p)
                ustar(p) = vkc*um(p)/log(zldis(p)/z0mg(p))
             enddo
             CALL roughness_lake (snl(c),t_grnd(c),t_lake(c,1),lake_icefrac(c,1),forc_psrf(1),&
                                  cur(p),ustar(p),z0mg(p),z0hg(p),z0qg(p))
           !===============================!
        else
           !=== Original Roughness lengths ===!
            if (t_grnd(c) >= tfrz) then  ! for unfrozen lake
                z0mg(p) = 0.001_r8        !original 0.01
            else if(snl(c) == 0 ) then   ! for frozen lake
            ! Is this okay even if it is snow covered?  What is the roughness over
            ! non-veg. snow?
            ! Should ice roughness 0.04 more reasonalbe ?? YangW
                z0mg(p) = 0.005_r8        ! original 0.04, now for frozen lake without snow
            else                         ! for frozen lake with snow   
                z0mg(p) = 0.0024_r8
            end if
            z0hg(p) = z0mg(p)
            z0qg(p) = z0mg(p)
           !===================================! 
        end if
      !============================================!

      !--- Initialize Monin-Obukhov length and wind speed ---!
       CALL MoninObukIni(ur(p), thv(c), dthv, zldis(p), z0mg(p), um(p), obu(p))

    end do

    iter                  = 1
    fncopy                = num_shlakep
    fpcopy(1:num_shlakep) = filter_shlakep(1:num_shlakep)
    ! Begin stability iteration
    ITERATION : do while (iter <= niters .and. fncopy > 0)

       ! Determine friction velocity, and potential temperature and humidity
       ! profiles of the surface boundary layer
         CALL FrictionVelocity(pgridcell,forc_hgt,forc_hgt_u,          & ! i
                               forc_hgt_t,forc_hgt_q,                  & ! i
                               lbp, ubp, fncopy, fpcopy,               & ! i
                               displa, z0mg, z0hg, z0qg,               & ! i
                               obu,  um,                               & ! i
                               ustar,temp1, temp2, temp12m, temp22m)     ! o
                            !--- YangW: Following Variables For dustmodel -> Removed ---!
                            !u10, fv, iter, ur, fm)                    ! i&o

! Here Original Concurent
       do fp = 1, fncopy
          p = fpcopy(fp)
          c = pcolumn(p)
          g = pgridcell(p)

          tgbef(c) = t_grnd(c)
          if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0) then
             tksur = savedtke1(c)
             ! Set this to the eddy conductivity from the last
             ! timestep, as the molecular conductivity will be orders of magnitude too small.
             ! Will have to deal with first timestep.
             tsur = t_lake(c,1)
          else if (snl(c) == 0) then  !frozen but no snow layers
             tksur = tkice
             tsur  = t_lake(c,1)
          else
          !Need to calculate thermal conductivity of the top snow layer
             bw    = (h2osoi_ice(c,jtop(c))+h2osoi_liq(c,jtop(c)))/dz(c,jtop(c))
             tksur = tkairc + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkairc)
             tsur  = t_soisno(c,jtop(c))
          end if

          ! Determine aerodynamic resistances
          ram(p)  = 1._r8/(ustar(p)*ustar(p)/um(p))
          rah(p)  = 1._r8/(temp1(p)*ustar(p))
          raw(p)  = 1._r8/(temp2(p)*ustar(p))
          ram1(p) = ram(p)   !pass value to global variable

          ! Get derivative of fluxes with respect to ground temperature
          stftg3(p) = emg*sb*tgbef(c)*tgbef(c)*tgbef(c)

          ! Changed surface temperature from t_lake(c,1) to tsur.
          ! Also adjusted so that if there are snow layers present, all radiation is absorbed in the top layer.
          ax  = betaprime(c)*sabg(p) + emg*forc_lwrad(g) + 3._r8*stftg3(p)*tgbef(c)     &
                + forc_rho(g)*cpair/rah(p)*thm(c)                                       &
                - htvp(c)*forc_rho(g)/raw(p)*(qsatg(c)-qsatgdT(c)*tgbef(c) - forc_q(g)) &
                + tksur*tsur/dzsur(c)
          !Changed sabg(p) and to betaprime(c)*sabg(p).
          bx  = 4._r8*stftg3(p) + forc_rho(g)*cpair/rah(p)                              &
                + htvp(c)*forc_rho(g)/raw(p)*qsatgdT(c) + tksur/dzsur(c)

          t_grnd(c) = ax/bx

          ! Update htvp
          if (t_grnd(c) > tfrz) then
              htvp(c) = hvap
          else
              htvp(c) = hsub
          end if

          ! Surface fluxes of momentum, sensible and latent heat
          ! using ground temperatures from previous time step
          eflx_sh_grnd(p)  = forc_rho(g)*cpair*(t_grnd(c)-thm(c))/rah(p)
          qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-tgbef(c))-forc_q(g))/raw(p)

          ! Re-calculate saturated vapor pressure, specific humidity and their
          ! derivatives at lake surface
          CALL QSat(t_grnd(c), forc_pbot(g), eg, degdT, qsatg(c), qsatgdT(c))

          dth(p)  = thm(c)-t_grnd(c)
          dqh(p)  = forc_q(g)-qsatg(c)
          tstar   = temp1(p)*dth(p)
          qstar   = temp2(p)*dqh(p)
          thvstar = tstar*(1._r8+0.61_r8*forc_q(g)) + 0.61_r8*forc_th(g)*qstar
          zeta    = zldis(p)*vkc * grav*thvstar/(ustar(p)**2*thv(c))
          if (zeta >= 0._r8) then     !stable
             zeta  = min(2._r8,max(zeta,0.01_r8))
             um(p) = max(ur(p),0.1_r8)
          else                        !unstable
             zeta  = max(-100._r8,min(zeta,-0.01_r8))
             wc    = beta1*(-grav*ustar(p)*thvstar*zii/thv(c))**0.333_r8
             um(p) = sqrt(ur(p)*ur(p)+wc*wc)
          end if
          obu(p) = zldis(p)/zeta

          if (obuold(p)*obu(p) < 0._r8) nmozsgn(p) = nmozsgn(p)+1

          obuold(p) = obu(p)

       end do   ! end of filtered pft loop

       !=== YangW: Updated Roughness Calculation ===!
         if (roughness_flag) then
             CALL roughness_lake (snl(c),t_grnd(c),t_lake(c,1),lake_icefrac(c,1),forc_psrf(1),&
                                  cur(p),ustar(p),z0mg(p),z0hg(p),z0qg(p))
         end if
       !============================================!

       iter = iter + 1
       if (iter <= niters ) then
          ! Rebuild copy of pft filter for next pass through the ITERATION loop
          fnold  = fncopy
          fncopy = 0
          do fp = 1, fnold
             p  = fpcopy(fp)
             if (nmozsgn(p) < 3) then
                fncopy         = fncopy + 1
                fpcopy(fncopy) = p
             end if
          end do       ! end of filtered pft loop
       end if

    end do ITERATION   ! end of stability iteration

! Here Original Concurent
    do fp = 1, num_shlakep
       p = filter_shlakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       ! If there is snow on the ground and t_grnd > tfrz: reset t_grnd = tfrz.
       ! Re-evaluate ground fluxes.
       ! h2osno > 0.5 prevents spurious fluxes.
       ! note that qsatg and qsatgdT should be f(tgbef) (PET: not sure what this
       ! comment means)
       ! Zack Subin, 3/27/09: Since they are now a function of whatever t_grnd was before cooling
       !    to freezing temperature, then this value should be used in the derivative correction term.
       ! Should this happen if the lake temperature is below freezing, too? I'll assume that for now.
       ! Also, allow convection if ground temp is colder than lake but warmer than 4C, or warmer than 
       !    lake which is warmer than freezing but less than 4C.
       if ( (h2osno(c) > 0.5_r8 .or. t_lake(c,1) <= tfrz) .and. t_grnd(c) > tfrz) then
      !if ( t_lake(c,1) <= tfrz .and. t_grnd(c) > tfrz) then       ! A Optional Judgement ? Yang W
          t_grnd_temp      = t_grnd(c)
          t_grnd(c)        = tfrz
          eflx_sh_grnd(p)  = forc_rho(g)*cpair*(t_grnd(c)-thm(c))/rah(p)
          qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-t_grnd_temp) - forc_q(g))/raw(p)
               !--- YangW Updated : Using the passed-in parameter tdmax to replace 277 ---!
       else if ( (t_lake(c,1) > t_grnd(c) .and. t_grnd(c) > tdmax) .or. &
                 (t_lake(c,1) < t_grnd(c) .and. t_lake(c,1) > tfrz .and. t_grnd(c) < tdmax) ) then
                 ! Convective mixing will occur at surface
          t_grnd_temp      = t_grnd(c)
          t_grnd(c)        = t_lake(c,1)
          eflx_sh_grnd(p)  = forc_rho(g)*cpair*(t_grnd(c)-thm(c))/rah(p)
          qflx_evap_soi(p) = forc_rho(g)*(qsatg(c)+qsatgdT(c)*(t_grnd(c)-t_grnd_temp) - forc_q(g))/raw(p)
       end if

       ! Update htvp
       if (t_grnd(c) > tfrz) then
           htvp(c) = hvap
       else
           htvp(c) = hsub
       end if

       ! Net longwave from ground to atmosphere
       !eflx_lwrad_out(p) = (1._r8-emg)*forc_lwrad(g) + stftg3(p)*(-3._r8*tgbef(c)+4._r8*t_grnd(c))
       ! What is tgbef doing in this equation? Can't it be exact now? --Zack Subin, 4/14/09
       eflx_lwrad_out(p) = (1._r8-emg)*forc_lwrad(g) + emg*sb*t_grnd(c)**4

       ! Ground heat flux
       eflx_soil_grnd(p) = sabg(p) + forc_lwrad(g) - eflx_lwrad_out(p) - &
                           eflx_sh_grnd(p) - htvp(c)*qflx_evap_soi(p)
       !Why is this sabg(p) and not beta*sabg(p)??
       !I've kept this as the incorrect sabg so that the energy balance check will be correct.
       !This is the effective energy flux into the ground including the lake solar absorption
       !below the surface.  The variable eflx_gnet will be used to pass the actual heat flux
       !from the ground interface into the lake.
       taux(p) = -forc_rho(g)*forc_u(g)/ram(p)
       tauy(p) = -forc_rho(g)*forc_v(g)/ram(p)

       eflx_sh_tot(p)   = eflx_sh_grnd(p)
       qflx_evap_tot(p) = qflx_evap_soi(p)
       eflx_lh_tot(p)   = htvp(c)*qflx_evap_soi(p)
       eflx_lh_grnd(p)  = htvp(c)*qflx_evap_soi(p)

       ! 2 m height air temperature
       t_ref2m(p) = thm(c) + temp1(p)*dth(p)*(1._r8/temp12m(p) - 1._r8/temp1(p))

       ! 2 m height specific humidity
       q_ref2m(p) = forc_q(g) + temp2(p)*dqh(p)*(1._r8/temp22m(p) - 1._r8/temp2(p))

       ! Energy residual used for melting snow
       ! Effectively moved to ShalLakeTemp

       ! Prepare for lake layer temperature calculations below
       ! fin(c) = betaprime * sabg(p) + forc_lwrad(g) - (eflx_lwrad_out(p) + &
       !          eflx_sh_tot(p) + eflx_lh_tot(p))
       ! NOW this is just the net ground heat flux calculated below.

       eflx_gnet(p) = betaprime(c) * sabg(p) + forc_lwrad(g) - (eflx_lwrad_out(p) + &
                      eflx_sh_tot(p) + eflx_lh_tot(p))
       ! This is the actual heat flux from the ground interface into the lake, not including
       ! the light that penetrates the surface.

!       u2m = max(1.0_r8,ustar(p)/vkc*log(2._r8/z0mg(p)))
       ! u2 often goes below 1 m/s; it seems like the only reason for this minimum is to
       ! keep it from being zero in the ks equation below; 0.1 m/s is a better limit for
       ! stable conditions --ZS
       u2m   = max(0.1_r8,ustar(p)/vkc*log(2._r8/z0mg(p)))
       ws(c) = 1.2e-03_r8 * u2m
       ks(c) = 6.6_r8*sqrt(abs(sin(lat(g))))*(u2m**(-1.84_r8))

    end do
    !End of surface flux relevant code in original BiogeophysicsLakeMod until history loop.

    ! The following are needed for global average on history tape.

    ! Here Original Concurent
    do fp = 1, num_shlakep
       p = filter_shlakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)
      !t_veg(p) = forc_t(g)
      !This is an odd choice, since elsewhere t_veg = t_grnd for bare ground.
      !Zack Subin, 4/09
       t_veg(p)          = t_grnd(c)
       eflx_lwrad_net(p) = eflx_lwrad_out(p) - forc_lwrad(g)
       qflx_prec_grnd(p) = forc_rain(g) + forc_snow(g)
    end do

END SUBROUTINE ShalLakeFluxes
 
SUBROUTINE ShalLakeTemperature(t_grnd,h2osno,sabg,dz,dz_lake,z,zi,         & !i
                               z_lake,ws,ks,snl,eflx_gnet,lakedepth,       &
                               lake_icefrac,snowdp,                        & !i&o
                               eflx_sh_grnd,eflx_sh_tot,eflx_soil_grnd,    & !o
                               t_lake,t_soisno,h2osoi_liq,                 &
                               h2osoi_ice,savedtke1,                       &
                               frac_iceold,qflx_snomelt,imelt,             & 
                          !--- YangW: Updated Eddy Diffusivity Calculation---!
                          !--- Following are Local Variabls ---!
                               drhodz,n2,ri,km,ke,ked,kw,                      &
                          !--- YangW : Updated External Variables From the namelist ---!
                            lakedebug_yw, beta_yw      ,eta_flag         , eta_scale_yw, eta_yw,        &
                            roughness_flag , diffusivity_flag, diffusivity_index_yw, tdmax, mixing_factor, mixing_factor_ked)
!=======================================================================================================
! !DESCRIPTION:
! Calculates temperatures in the 20-25 layer column of (possible) snow,
! lake water, and soil beneath lake.
! Snow and soil temperatures are determined as in SoilTemperature, except
! for appropriate boundary conditions at the top of the snow (the flux is fixed
! to be the ground heat flux calculated in ShalLakeFluxes), the bottom of the snow
! (adjacent to top lake layer), and the top of the soil (adjacent to the bottom
! lake layer). Also, the soil is assumed to be always fully saturated (ShalLakeHydrology
! will have to insure this). The whole column is solved simultaneously as one tridiagonal matrix.
! Lake temperatures are determined from the Hostetler model as before, except now:
!    i) Lake water layers can freeze by any fraction and release latent heat; thermal
!       and mechanical properties are adjusted for ice fraction.
!   ii) Convective mixing (though not eddy diffusion) still occurs for frozen lakes.
!  iii) No sunlight is absorbed in the lake if there are snow layers.
!   iv) Light is allowed to reach the top soil layer (where it is assumed to be completely absorbed).
!    v) Lakes have variable depth, set ultimately in surface data set but now in initShalLakeMod.
!
! Eddy + molecular diffusion:
! d ts    d            d ts     1 ds
! ---- = -- [(km + ke) ----] + -- --
!  dt    dz             dz     cw dz
!
! where: ts = temperature (kelvin)
!         t = time (s)
!         z = depth (m)
!        km = molecular diffusion coefficient (m**2/s)
!        ke = eddy diffusion coefficient (m**2/s)
!        cw = heat capacity (j/m**3/kelvin)
!         s = heat source term (w/m**2)
!
!   Shallow lakes are allowed to have variable depth, set in _____.
!
!   For shallow lakes:    ke > 0 if unfrozen,
!       and convective mixing occurs WHETHER OR NOT frozen. (See e.g. Martynov...)
!
! Use the Crank-Nicholson method to set up tridiagonal system of equations to
! solve for ts at time n+1, where the temperature equation for layer i is
! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
!
! The solution conserves energy as:
!
! [For lake layers]
! cw*([ts(      1)] n+1 - [ts(      1)] n)*dz(      1)/dt + ... +
! cw*([ts(nlevlake)] n+1 - [ts(nlevlake)] n)*dz(nlevlake)/dt = fin
! But now there is phase change, so cv is not constant and there is
! latent heat.
!
! where:
! [ts] n   = old temperature (kelvin)
! [ts] n+1 = new temperature (kelvin)
! fin      = heat flux into lake (w/m**2)
!          = betaprime*sabg + forc_lwrad - eflx_lwrad_out - eflx_sh_tot - eflx_lh_tot
!          (This is now the same as the ground heat flux.)
!            + phi(1) + ... + phi(nlevlake) + phi(top soil level)
! betaprime = beta(islak) for no snow layers, and 1 for snow layers.
! This assumes all radiation is absorbed in the top snow layer and will need
! to be changed for CLM 4.
!
! WARNING: This subroutine assumes lake columns have one and only one pft.
!
! Outline:
! 1!) Initialization
! 2!) Lake density
! 3!) Diffusivity
! 4!) Heat source term from solar radiation penetrating lake
! 5!) Set thermal props and find initial energy content
! 6!) Set up vectors for tridiagonal matrix solution
! 7!) Solve tridiagonal and back-substitute
! 8!) (Optional) Do first energy check using temperature change at constant heat capacity.
! 9!) Phase change
! 9.5!) (Optional) Do second energy check using temperature change and latent heat, considering changed heat capacity.
!                  Also do soil water balance check.
!10!) Convective mixing 
!11!) Do final energy check to detect small numerical errors (especially from convection)
!     and dump small imbalance into sensible heat, or pass large errors to BalanceCheckMod for abort.
!
! REVISION HISTORY:
! Created by Zack Subin, 2009.
! Reedited by Hongping Gu, 2010.
!=========================================================================================================
!use TridiagonalMod     , only : Tridiagonal
    implicit none
!in:
    real(r8), intent(in)    :: t_grnd(1)                           ! ground temperature (Kelvin)
    real(r8), intent(inout) :: h2osno(1)                           ! snow water (mm H2O)
    real(r8), intent(in)    :: sabg(1)                             ! solar radiation absorbed by ground (W/m**2)
    real(r8), intent(in)    :: dz(1,-nlevsnow + 1:nlevsoil)        ! layer thickness for snow & soil (m)
    real(r8), intent(in)    :: dz_lake(1,nlevlake)                 ! layer thickness for lake (m)
    real(r8), intent(in)    :: z(1,-nlevsnow+1:nlevsoil)           ! layer depth for snow & soil (m)
    real(r8), intent(in)    :: zi(1,-nlevsnow+0:nlevsoil)          ! interface level below a "z" level (m)
                                                                   ! the other z and dz variables
    real(r8), intent(in)    :: z_lake(1,nlevlake)                  ! layer depth for lake (m)
    real(r8), intent(in)    :: ws(1)                               ! surface friction velocity (m/s)
    real(r8), intent(in)    :: ks(1)                               ! coefficient passed to ShalLakeTemperature
                                                                   ! for calculation of decay of eddy diffusivity with depth
    integer , intent(in)    :: snl(1)                              ! negative of number of snow layers
    real(r8), intent(inout) :: eflx_gnet(1)                        ! net heat flux into ground (W/m**2) at the surface interface
    real(r8), intent(in)    :: lakedepth(1)                        ! column lake depth (m)
   !real(r8), intent(in)    :: watsat(1,nlevsoil)                  ! volumetric soil water at saturation (porosity)
    real(r8), intent(inout) :: snowdp(1)                           ! snow height (m)
!out: 
    real(r8), intent(out)   :: eflx_sh_grnd(1)                     ! sensible heat flux from ground (W/m**2) [+ to atm]
    real(r8), intent(out)   :: eflx_sh_tot(1)                      ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), intent(out)   :: eflx_soil_grnd(1)                   ! heat flux into snow / lake (W/m**2) [+ = into soil]
                                                                   ! Here this includes the whole lake radiation absorbed.
    real(r8), intent(inout) :: t_lake(1,nlevlake)                  ! lake temperature (Kelvin)
    real(r8), intent(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)    ! soil (or snow) temperature (Kelvin)
    real(r8), intent(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)  ! liquid water (kg/m2) [for snow & soil layers]
    real(r8), intent(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)  ! ice lens (kg/m2) [for snow & soil layers]
    real(r8), intent(inout) :: lake_icefrac(1,nlevlake)            ! mass fraction of lake layer that is frozen
    real(r8), intent(out)   :: savedtke1(1)                        ! top level thermal conductivity (W/mK)
    real(r8), intent(out)   :: frac_iceold(1,-nlevsnow+1:nlevsoil) ! fraction of ice relative to the tot water
    real(r8), intent(out)   :: qflx_snomelt(1)                     ! snow melt (mm H2O /s)
    integer, intent(out)    :: imelt(1,-nlevsnow+1:nlevsoil)       ! flag for melting (=1), freezing (=2), Not=0 (new)

! OTHER LOCAL VARIABLES:
    integer , parameter  :: islak = 2   ! index of lake, 1 = deep lake, 2 = shallow lake
    real(r8), parameter  :: p0 = 1._r8  ! neutral value of turbulent prandtl number
    integer  :: i,j,fc,fp,g,c,p         ! do loop or array index
   !real(r8) :: dtime                   ! land model time step (sec)
   !--- YangW: beta will be Defined as a passin parameter ---!
   !real(r8) :: beta(2)                 ! fraction solar rad absorbed at surface: depends on lake type
    real(r8) :: za(2)                   ! base of surface absorption layer (m): depends on lake type
    real(r8) :: eta(2)                  ! light extinction coefficient (/m): depends on lake type
    real(r8) :: cwat                    ! specific heat capacity of water (j/m**3/kelvin)
    real(r8) :: cice_eff                ! effective heat capacity of ice (using density of
                                        ! water because layer depth is not adjusted when freezing
    real(r8) :: cfus                    ! effective heat of fusion per unit volume
                                        ! using water density as above
    real(r8) :: tkice_eff               ! effective conductivity since layer depth is constant
    real(r8) :: a(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)      ! "a" vector for tridiagonal matrix
    real(r8) :: b(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)      ! "b" vector for tridiagonal matrix
    real(r8) :: c1(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)     ! "c" vector for tridiagonal matrix
    real(r8) :: r(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)      ! "r" vector for tridiagonal solution
    real(r8) :: rhow(lbc:ubc,nlevlake)  ! density of water (kg/m**3)
    real(r8) :: phi(lbc:ubc,nlevlake)   ! solar radiation absorbed by layer (w/m**2)
    real(r8) :: rsfin                   ! relative flux of solar radiation into layer
    real(r8) :: rsfout                  ! relative flux of solar radiation out of layer
    real(r8) :: phi_soil(lbc:ubc)       ! solar radiation into top soil layer (W/m**2)
    real(r8) :: fin(lbc:ubc)            ! net heat flux into lake at ground interface (w/m**2)
    real(r8) :: ocvts(lbc:ubc)          ! (cwat*(t_lake[n  ])*dz
    real(r8) :: ncvts(lbc:ubc)          ! (cwat*(t_lake[n+1])*dz
    real(r8) :: zin                     ! depth at top of layer (m)
    real(r8) :: zout                    ! depth at bottom of layer (m)
    real(r8) :: num                     ! used in calculating ri
    real(r8) :: den                     ! used in calculating ri
    real(r8) :: tav_froz(lbc:ubc)       ! used in aver temp for convectively mixed layers (C)
    real(r8) :: tav_unfr(lbc:ubc)       ! "
    real(r8) :: nav(lbc:ubc)            ! used in aver temp for convectively mixed layers
    real(r8) :: phidum                  ! temporary value of phi
    real(r8) :: iceav(lbc:ubc)          ! used in calc aver ice for convectively mixed layers
    real(r8) :: qav(lbc:ubc)            ! used in calc aver heat content for conv. mixed layers
    integer  :: jtop(lbc:ubc)           ! top level for each column (no longer all 1)
    real(r8) :: cv (lbc:ubc,-nlevsnow+1:nlevsoil)            ! heat capacity of soil/snow [J/(m2 K)]
    real(r8) :: tk (lbc:ubc,-nlevsnow+1:nlevsoil)            ! thermal conductivity of soil/snow [W/(m K)]
                                                             ! (at interface below, except for j=0)
    real(r8) :: cv_lake (lbc:ubc,1:nlevlake)                 ! heat capacity [J/(m2 K)]
    real(r8) :: tk_lake (lbc:ubc,1:nlevlake)                 ! thermal conductivity at layer node [W/(m K)]
    real(r8) :: cvx (lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)  ! heat capacity for whole column [J/(m2 K)]
    real(r8) :: tkix(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)  ! thermal conductivity at layer interfaces
                                                             ! for whole column [W/(m K)]
    real(r8) :: tx(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)    ! temperature of whole column [K]
    real(r8) :: tktopsoillay(lbc:ubc)                        ! thermal conductivity [W/(m K)]
    real(r8) :: fnx(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)   ! heat diffusion through the layer interface below [W/m2]
    real(r8) :: phix(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)  ! solar source term for whole column [W/m**2]
    real(r8) :: zx(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil)    ! interface depth (+ below surface) for whole column [m]
    real(r8) :: dzm                                          ! used in computing tridiagonal matrix [m]
    real(r8) :: dzp                                          ! used in computing tridiagonal matrix [m]
    integer  :: jprime                                       ! j - nlevlake
    real(r8) :: factx(lbc:ubc,-nlevsnow+1:nlevlake+nlevsoil) ! coefficient used in computing tridiagonal matrix
    real(r8) :: t_lake_bef(lbc:ubc,1:nlevlake)               ! beginning lake temp for energy conservation check [K]
    real(r8) :: t_soisno_bef(lbc:ubc,-nlevsnow+1:nlevsoil)   ! beginning soil temp for E cons. check [K]
    real(r8) :: lhabs(lbc:ubc)        ! total per-column latent heat abs. from phase change  (J/m^2)
    real(r8) :: esum1(lbc:ubc)        ! temp for checking energy (J/m^2)
    real(r8) :: esum2(lbc:ubc)        ! ""
    real(r8) :: zsum(lbc:ubc)         ! temp for putting ice at the top during convection (m)
    real(r8) :: errsoi(1)             ! soil/lake energy conservation error (W/m**2)
    real(r8) :: eflx_snomelt(1)       ! snow melt heat flux (W/m**2)
    CHARACTER*256 :: message
!
! Constants for lake temperature model
   !--- YangW Updated: beta will be Defined as a passin parameter ---!
   !data beta/0.4_r8, 0.4_r8/  ! (deep lake, shallow lake)
    data za  /0.6_r8, 0.6_r8/
!   For now, keep beta and za for shallow lake the same as deep lake, until better data is found.
!   It looks like eta is key and that larger values give better results for shallow lakes.  Use
!   empirical expression from Hakanson (below). This is still a very unconstrained parameter
!   that deserves more attention.
!   Some radiation will be allowed to reach the soil.
!-----------------------------------------------------------------------

 !=== Variables For Updated Physical Part ===!
  !--- YangW: Why the sum soil water part Removed ??? ---
   !real(r8) :: wsum(lbc:ubc)              ! temp for checking water (kg/m^2)
   !real(r8) :: wsum_end(lbc:ubc)          ! temp for checking water (kg/m^2)
  !--- YangW: Updated Eddy Diffusivity Calculation ---!
    !--- For solar penetration scheme by ZJ---! 
      real(r8) :: rp
      real(r8) :: x1
      real(r8) :: x2
    !-----------------------------------------!

    real(r8) :: drhodz(nlevlake-1)          ! d[rhow]/dz(kg/m**4)
    real(r8) :: n2(nlevlake-1)              ! brunt-vaisala frequency (/s**2)
    real(r8) :: ri(nlevlake-1)              ! Richardson number
    real(r8),intent(out) :: km                   ! molecular diffusion coefficient (m**2/s)
    real(r8),intent(out) :: ke(nlevlake-1)       ! eddy diffusion coefficient (m**2/s)
    real(r8),intent(out) :: ked(nlevlake-1)      ! enhanced diffusion coefficient by unresolved mixing processes(m**2/s)
    real(r8),intent(out) :: kw(lbc:ubc,nlevlake) ! molecular + eddy + enhanced diffusivity coefficient (m**2/s)
 !===========================================!

  !=== Some Variables Prepared For Coupled WRF-Lake namelist -> Subroutine Lake ===!
    logical, intent(in) :: lakedebug_yw
    real(r8)            :: beta(2)                              ! fraction solar rad absorbed at surface: depends on lake type
    real   , intent(in) :: beta_yw   

    !--- For Subroutine Lake -> LakeMain -> ShalLakeFluxes ---!
    logical, intent(in) :: roughness_flag    

    !--- For Subroutine Lake -> LakeMain -> ShalLakeTemperatures ---!
    logical, intent(in) :: diffusivity_flag       ! True For using the New diffusivity scheme
    integer, intent(in) :: diffusivity_index_yw   ! 1 for kw=mixing_factor*(km+ke+ked)
                                                  ! 2 for kw=km+mixing_factor*ke+ked
    real, intent(in)    :: tdmax
    real, intent(in)    :: mixing_factor    
    real, intent(in)    :: mixing_factor_ked    
    logical, intent(in) :: eta_flag
    real, intent(in)    :: eta_scale_yw
    real, intent(in)    :: eta_yw

    beta(:) = beta_yw
  !================================================================================!

    ! 1!) Initialization
    ! Determine step size
    ! dtime = get_step_size()

    ! Initialize constants
    cwat     = cpliq*denh2o           ! water heat capacity per unit volume
    cice_eff = cpice*denh2o           ! use water density because layer depth is not adjusted
                                      ! for freezing
    cfus      = hfus*denh2o           ! latent heat per unit volume
    tkice_eff = tkice * denice/denh2o ! effective conductivity since layer depth is constant
    km        = tkwat/cwat            ! a constant (molecular diffusivity)

    ! Begin calculations
    do fc = 1, num_shlakec
       c = filter_shlakec(fc)
       ! Initialize Ebal quantities computed below
       ocvts(c) = 0._r8
       ncvts(c) = 0._r8
       esum1(c) = 0._r8
       esum2(c) = 0._r8
    end do

    ! Initialize set of previous time-step variables as in DriverInit,
    ! which is currently not called over lakes. This has to be done
    ! here because phase change will occur in this routine.
    ! Ice fraction of snow at previous time step
    do j = -nlevsnow+1,0
      do fc = 1, num_shlakec
         c = filter_shlakec(fc)
         if (j >= snl(c) + 1) then
            frac_iceold(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))
         end if
      end do
    end do

    !--- YangW: Why sum soil water part removed ??? ---!
    ! Sum soil water.
    !do j = 1, nlevsoil
    !   do fc = 1, num_shlakec
    !      c = filter_shlakec(fc)
    !      if (j == 1) wsum(c) = 0._r8
    !          wsum(c) = wsum(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
    !   end do
    !end do

    do fp = 1, num_shlakep
       p = filter_shlakep(fp)
       c = pcolumn(p)
       ! Prepare for lake layer temperature calculations below
       ! fin(c) = betaprime * sabg(p) + forc_lwrad(g) - (eflx_lwrad_out(p) + &
       !          eflx_sh_tot(p) + eflx_lh_tot(p)) 
       ! fin(c) now passed from ShalLakeFluxes as eflx_gnet
       fin(c) = eflx_gnet(p)
    end do

    ! 2!) Lake density
    do j = 1, nlevlake
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)
          rhow(c,j) = (1._r8 - lake_icefrac(c,j)) * & 
                      !1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,j)-277._r8))**1.68_r8 ) &
                      !--- YangW Updated : Using the passed-in parameter tdmax to replace 277 ---!
                       1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,j)-tdmax))**1.68_r8 ) &
                       + lake_icefrac(c,j)*denice
                    ! Allow for ice fraction; assume constant ice density.
                    ! Is this the right weighted average?
                    ! Using this average will make sure that surface ice is treated properly during
                    ! convective mixing.
       end do
    end do

    ! 3!) Diffusivity and implied thermal "conductivity" = diffusivity * cwat
    do j = 1, nlevlake-1
! Here Original Prevector
! Here Original Concurent
       do fc = 1, num_shlakec
          c         = filter_shlakec(fc)
          drhodz(j) = (rhow(c,j+1)-rhow(c,j)) / (z_lake(c,j+1)-z_lake(c,j))
          n2(j)     = grav / rhow(c,j) * drhodz(j)
          ! YangW: must set n2 >= 7.5e-5, or the ked will be infinity
          n2(j)     = max(n2(j),7.5e-5)                           
          ! Fixed sign error here: our z goes up going down into the lake, so no negative
          ! sign is needed to make this positive unlike in Hostetler. --ZS
          num    = 40._r8 * n2(j) * (vkc*z_lake(c,j))**2
          den    = max( (ws(c)**2) * exp(-2._r8*ks(c)*z_lake(c,j)), 1.e-10_r8 )
          ri(j)  = ( -1._r8 + sqrt( max(1._r8+num/den, 0._r8) ) ) / 20._r8

      !=== YangW: Updated Eddy Diffusivity & t_lake Calculation ===!
        if (diffusivity_flag) then  ! True For Adopting Eddy Diffusivity Calculation by Subin
           ked(j) = 1.04e-8_r8*n2(j)**(-0.43)                      ! For diffusivity_flag=true
           if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0) then
               ke(j)    = vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri(j)*ri(j)) 
               if (lakedepth(1) >= 25) then
                   if (diffusivity_index_yw == 1) then
                      kw(c,j) = km + ke(j) +ked(j) 
                      kw(c,j) = kw(c,j)*mixing_factor              ! Enhance mixing intensity lake deep than 25m
                   end if
                   if (diffusivity_index_yw == 2) then
                      kw(c,j) = km + ke(j)*mixing_factor +ked(j) 
                   end if
                   if (diffusivity_index_yw == 3) then
                      kw(c,j) = km + ke(j)*mixing_factor 
                   end if
                   if (diffusivity_index_yw == 4) then
                      if( t_lake(c,1) > tdmax ) then       
                         if (lakedepth(c) > 25.0 ) then 
                             kw(c,j) = km + ke(j)*mixing_factor 
                         else
                             kw(c,j) = km + ke(j)
                         end if
                      else 
                         if (lakedepth(c) > 25.0 ) then 
                             if (lakedepth(c) > 150.0 ) then 
                                 kw(c,j) = km + ke(j)*1.e+5_r8
                             else 
                                 kw(c,j) = km + ke(j)*1.e+4_r8
                             end if
                          else
                             kw(c,j) = km + ke(j)
                          end if
                      end if
                   end if
                   if (diffusivity_index_yw == 5) then
                      kw(c,j) = km + min(ke(j)*mixing_factor,0.01)+ked(j)*mixing_factor_ked 
                   end if
               else
                   kw(c,j)  = km + ke(j) +ked(j) 
               end if 
               tk_lake(c,j) = kw(c,j)*cwat
               ! If there is some ice in this layer (this should rarely happen because the surface
               ! is unfrozen and it will be unstable), still use the cwat to get out the tk b/c the eddy
               ! diffusivity equation assumes water.
           else                              ! frozen lake	 
               kw(c,j) = (km+ked(j))             
               if (lakedepth(1) >= 25) then  ! by Anning Huang	 
               ! kw(c,j) = kw(c,j)* 35.        
               end if
               tk_lake(c,j) = tkwat*tkice_eff /((1._r8-lake_icefrac(c,j))*tkice_eff &
                              + tkwat*lake_icefrac(c,j) )
               ! Assume the resistances add as for the calculation of conductivities at layer interfaces.
           end if
        else
           if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0) then
             ! ke = vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri*ri)
             ! if( t_lake(c,1) > 277.15_r8 ) then  ! Original Gu 
             !--- YangW Updated : Using the passed-in parameter tdmax to replace 277 ---!
              if( t_lake(c,1) > tdmax ) then       
                 if (lakedepth(c) > 15.0 ) then 
                    ke(j) = 1.e+2_r8*vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri(j)*ri(j))
                 else 
                    ke(j) = vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri(j)*ri(j))
                 endif
              else 
                 if (lakedepth(c) > 15.0 ) then 
                   if (lakedepth(c) > 150.0 ) then 
                       ke(j) = 1.e+5_r8*vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri(j)*ri(j))
                   else 
                       ke(j) =1.e+4_r8*vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri(j)*ri(j))
                   end if
                 else 
                   ke(j) = vkc*ws(c)*z_lake(c,j)/p0 * exp(-ks(c)*z_lake(c,j)) / (1._r8+37._r8*ri(j)*ri(j))
                 endif 
              end if
              kw(c,j)      = km + ke(j)
              tk_lake(c,j) = kw(c,j)*cwat
              ! If there is some ice in this layer (this should rarely happen because the surface
              ! is unfrozen and it will be unstable), still use the cwat to get out the tk b/c the eddy
              ! diffusivity equation assumes water.
           else
              kw(c,j)      = km
              tk_lake(c,j) = tkwat*tkice_eff / ( (1._r8-lake_icefrac(c,j))*tkice_eff &
                             + tkwat*lake_icefrac(c,j) )
              ! Assume the resistances add as for the calculation of conductivities at layer interfaces.
           end if
        end if
      !============================================================!
       end do         ! fc loop
    end do            ! j  loop ()

! Here Original Concurent
    do fc = 1, num_shlakec
       c               = filter_shlakec(fc)
       j               = nlevlake
       kw(c,nlevlake)  = kw(c,nlevlake-1)
       if (t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0) then
           tk_lake(c,j) = tk_lake(c,j-1)
       else
           tk_lake(c,j) = tkwat*tkice_eff / ( (1._r8-lake_icefrac(c,j))*tkice_eff &
                          + tkwat*lake_icefrac(c,j) )
       end if
       ! Use in surface flux calculation for next timestep.
       savedtke1(c) = kw(c,1)*cwat ! Will only be used if unfrozen
       ! set number of column levels for use by Tridiagonal below
       jtop(c) = snl(c) + 1
    end do

    ! 4!) Heat source term: unfrozen lakes only
    do j = 1, nlevlake
! Here Original Concurent
       do fp = 1, num_shlakep
          p = filter_shlakep(fp)
          c = pcolumn(p)

          !--- YangW Updated: eta will be Defined as a passin parameter ---!
            if (eta_flag) then
                if(eta_scale_yw > 0.) then
                   eta(:) = eta_scale_yw*1.1925_r8*lakedepth(c)**(-0.424)   ! Scale the eta
                 else
                   eta(:) = eta_yw
                 end if
            else
                ! Set eta(:), the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
                ! (regression of Secchi Depth with lake depth for small glacial basin lakes), and the
                ! Poole & Atkins expression for extinction coeffient of 1.7 / Secchi Depth (m).
                eta(:) = 1.1925_r8*lakedepth(c)**(-0.424)
            end if
            if(lakedebug_yw) then
               print *,"1.1925_r8*lakedepth(c)**(-0.424)=",1.1925_r8*lakedepth(c)**(-0.424)
               print *,"eta Now",eta(:)
            end if
          !----------------------------------------------------------------!

          !--- YangW Updated: change the solar penetration scheme Coded by ZJ ---!
           zin    = z_lake(c,j) - 0.5_r8*dz_lake(c,j)
           zout   = z_lake(c,j) + 0.5_r8*dz_lake(c,j)
           if (.false.) then
               rp     = 0.58               ! Clayton A. Paulson, Irradiance Measurements in the Upper Ocean, journal of physical oceanorgraphy
               x1     = 0.3
               x2     = 23
               rsfin  = rp*exp(-max(zin-za(islak),0._r8)/x1)+(1._r8-rp)*exp(-max(zin -za(islak),0._r8)/x2)
               rsfout = rp*exp(-max(zout-za(islak),0._r8)/x1)+(1._r8-rp)*exp(-max(zout-za(islak),0._r8)/x2)
           else
               rsfin  = exp( -eta(islak)*max(  zin-za(islak),0._r8 ) )
               rsfout = exp( -eta(islak)*max( zout-za(islak),0._r8 ) )
           end if
          !-----------------------------------------------------------------------!

          ! Let rsfout for bottom layer go into soil.
          ! This looks like it should be robust even for pathological cases,
          ! like lakes thinner than za.
          if(t_grnd(c) > tfrz .and. t_lake(c,1) > tfrz .and. snl(c) == 0) then ! For unfrozen lake
             phidum = (rsfin-rsfout) * sabg(p) * (1._r8-beta(islak))
             if(j == nlevlake) then
                phi_soil(c) = rsfout * sabg(p) * (1._r8-beta(islak))
             end if
          else if(j == 1 .and. snl(c) == 0) then                               ! if frozen but no snow layers
                  phidum = sabg(p) * (1._r8-beta(islak))
               else                                                            ! For Snow-Covered iced lake: radiation absorbed at surface
                  phidum = 0._r8
               if (j == nlevlake) phi_soil(c) = 0._r8
          end if
          phi(c,j) = phidum
       end do
    end do

    ! 5!) Set thermal properties and check initial energy content.
    ! For lake
    do j = 1, nlevlake
! Here Original Concurent
       do fc = 1, num_shlakec
          c            = filter_shlakec(fc)
          cv_lake(c,j) = dz_lake(c,j) * (cwat*(1._r8-lake_icefrac(c,j)) + cice_eff*lake_icefrac(c,j))
       end do
    end do

    ! For snow / soil
    CALL SoilThermProp_Lake (snl,dz,zi,z,t_soisno,h2osoi_liq,h2osoi_ice,    &
                             tk, cv, tktopsoillay)

    ! Sum cv*t_lake for energy check
    ! Include latent heat term, and correction for changing heat capacity with phase change.

    ! This will need to be over all soil / lake / snow layers. Lake is below.
    do j = 1, nlevlake
! Here Original Concurent
       do fc = 1, num_shlakec
          c        = filter_shlakec(fc)
         !ocvts(c) = ocvts(c) + cv_lake(c,j)*t_lake(c,j) &
          ocvts(c) = ocvts(c) + cv_lake(c,j)*(t_lake(c,j)-tfrz) &
                     + cfus*dz_lake(c,j)*(1._r8-lake_icefrac(c,j)) !&
                    !+ (cwat-cice_eff)*lake_icefrac(c)*tfrz*dz_lake(c,j) !enthalpy reconciliation term
          t_lake_bef(c,j) = t_lake(c,j)
       end do
    end do

    ! Now do for soil / snow layers
    do j = -nlevsnow + 1, nlevsoil
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)
          if (j >= jtop(c)) then
             !ocvts(c) = ocvts(c) + cv(c,j)*t_soisno(c,j) &
              ocvts(c) = ocvts(c) + cv(c,j)*(t_soisno(c,j)-tfrz) &
                         + hfus*h2osoi_liq(c,j) !&
                        !+ (cpliq-cpice)*h2osoi_ice(c,j)*tfrz !enthalpy reconciliation term
             if (j == 1 .and. h2osno(c) > 0._r8 .and. j == jtop(c)) then
                 ocvts(c) = ocvts(c) - h2osno(c)*hfus
             end if
             t_soisno_bef(c,j) = t_soisno(c,j)
             if(abs(t_soisno(c,j)-288) > 150)   then 
                if (lakedebug_yw) then            ! YangW: This part is designed in the original model
                    WRITE( message,* ) 'WARNING: Extreme t_soisno at c, level',c, j
                    print *, message
                end if
               !CALL wrf_error_fatal ( message )
             endif
          end if
       end do
    end do

    ! 6!) Set up vector r and vectors a, b, c1 that define tridiagonal matrix

    ! Heat capacity and resistance of snow without snow layers (<1cm) is ignored during diffusion,
    ! but its capacity to absorb latent heat may be used during phase change.

    ! Set up interface depths, zx, heat capacities, cvx, solar source terms, phix, and temperatures, tx.
    do j = -nlevsnow+1, nlevlake+nlevsoil
! Here Original Prevector
! Here Original Concurent
       do fc = 1,num_shlakec
          c      = filter_shlakec(fc)
          jprime = j - nlevlake
          if (j >= jtop(c)) then
             if (j < 1) then                ! snow layer
                 zx(c,j)   = z(c,j)
                 cvx(c,j)  = cv(c,j)
                 phix(c,j) = 0._r8
                 tx(c,j)   = t_soisno(c,j)
             else if (j <= nlevlake) then   ! lake layer
                 zx(c,j)   = z_lake(c,j)
                 cvx(c,j)  = cv_lake(c,j)
                 phix(c,j) = phi(c,j)
                 tx(c,j)   = t_lake(c,j)
             else                           ! soil layer
                zx(c,j)    = zx(c,nlevlake) + dz_lake(c,nlevlake)/2._r8 + z(c,jprime)
                cvx(c,j)   = cv(c,jprime)
                if (j == nlevlake + 1) then ! top soil layer
                    phix(c,j) = phi_soil(c)
                else                        ! middle or bottom soil layer
                    phix(c,j) = 0._r8
                end if
                tx(c,j) = t_soisno(c,jprime)
             end if
          end if
       end do
    end do

    ! Determine interface thermal conductivities, tkix
    do j = -nlevsnow+1, nlevlake+nlevsoil
! Here Original Prevector
! Here Original Concurent
       do fc = 1,num_shlakec
          c      = filter_shlakec(fc)
          jprime = j - nlevlake
          if (j >= jtop(c)) then
             if (j < 0) then              ! non-bottom snow layer
                 tkix(c,j) = tk(c,j)
             else if (j == 0) then        ! bottom snow layer
                 dzp       = zx(c,j+1) - zx(c,j)
                 tkix(c,j) = tk_lake(c,1)*tk(c,j)*dzp / &
                             (tk(c,j)*z_lake(c,1) + tk_lake(c,1)*(-z(c,j)) )
                             ! tk(c,0) is the conductivity at the middle of that layer, as defined in SoilThermProp_Lake
             else if (j < nlevlake) then  ! non-bottom lake layer
                 tkix(c,j) = ( tk_lake(c,j)*tk_lake(c,j+1) * (dz_lake(c,j+1)+dz_lake(c,j)) ) &
                             / ( tk_lake(c,j)*dz_lake(c,j+1) + tk_lake(c,j+1)*dz_lake(c,j) )
             else if (j == nlevlake) then ! bottom lake layer
                 dzp       = zx(c,j+1) - zx(c,j)
                 tkix(c,j) = (tktopsoillay(c)*tk_lake(c,j)*dzp / &
                             (tktopsoillay(c)*dz_lake(c,j)/2._r8 + tk_lake(c,j)*z(c,1) ) )
                             ! tktopsoillay is the conductivity at the middle of that layer, as defined in SoilThermProp_Lake
             else                         ! soil layer
                 tkix(c,j) = tk(c,jprime)
             end if
         end if
      end do 
   end do


    ! Determine heat diffusion through the layer interface and factor used in computing
    ! tridiagonal matrix and set up vector r and vectors a, b, c1 that define tridiagonal
    ! matrix and solve system
    do j = -nlevsnow+1, nlevlake+nlevsoil
! Here Original Prevector
! Here Original Concurent
       do fc = 1,num_shlakec
          c = filter_shlakec(fc)
          if (j >= jtop(c)) then
             if (j < nlevlake+nlevsoil) then ! top or interior layer
                 factx(c,j) = dtime/cvx(c,j)
                 fnx(c,j)   = tkix(c,j)*(tx(c,j+1)-tx(c,j))/(zx(c,j+1)-zx(c,j))
             else                            ! bottom soil layer
                 factx(c,j) = dtime/cvx(c,j)
                 fnx(c,j)   = 0._r8           ! not used
             end if
          end if
       enddo
    end do

    do j = -nlevsnow+1,nlevlake+nlevsoil
! Here Original Prevector
! Here Original Concurent
       do fc = 1,num_shlakec
          c = filter_shlakec(fc)
          if (j >= jtop(c)) then
             if (j == jtop(c)) then               ! top layer
                dzp     = zx(c,j+1)-zx(c,j)
                a(c,j)  = 0._r8
                b(c,j)  = 1+(1._r8-cnfac)*factx(c,j)*tkix(c,j)/dzp
                c1(c,j) =  -(1._r8-cnfac)*factx(c,j)*tkix(c,j)/dzp
                r(c,j)  = tx(c,j) + factx(c,j)*( fin(c) + phix(c,j) + cnfac*fnx(c,j) )
             else if (j < nlevlake+nlevsoil) then ! middle layer
                dzm     = (zx(c,j)-zx(c,j-1))
                dzp     = (zx(c,j+1)-zx(c,j))
                a(c,j)  =   - (1._r8-cnfac)*factx(c,j)* tkix(c,j-1)/dzm
                b(c,j)  = 1._r8+ (1._r8-cnfac)*factx(c,j)*(tkix(c,j)/dzp + tkix(c,j-1)/dzm)
                c1(c,j) =   - (1._r8-cnfac)*factx(c,j)* tkix(c,j)/dzp
                r(c,j)  = tx(c,j) + cnfac*factx(c,j)*( fnx(c,j) - fnx(c,j-1) ) + factx(c,j)*phix(c,j)
             else                                 ! bottom soil layer
                dzm     = (zx(c,j)-zx(c,j-1))
                a(c,j)  =   - (1._r8-cnfac)*factx(c,j)*tkix(c,j-1)/dzm
                b(c,j)  = 1._r8+ (1._r8-cnfac)*factx(c,j)*tkix(c,j-1)/dzm
                c1(c,j) = 0._r8
                r(c,j)  = tx(c,j) - cnfac*factx(c,j)*fnx(c,j-1)
             end if
          end if
       enddo
    end do

   ! 7!) Solve for tdsolution
    CALL Tridiagonal(lbc, ubc, -nlevsnow + 1, nlevlake + nlevsoil, jtop, num_shlakec, filter_shlakec, &
                     a, b, c1, r, tx)
 
    ! Set t_soisno and t_lake
    do j = -nlevsnow+1, nlevlake + nlevsoil
! Here Original Concurent
       do fc = 1, num_shlakec
          c      = filter_shlakec(fc)
          jprime = j - nlevlake
          ! Don't do anything with invalid snow layers.
          if (j >= jtop(c)) then
             if (j < 1) then              ! snow layer
                 t_soisno(c,j) = tx(c,j)
             else if (j <= nlevlake) then ! lake layer
                 t_lake(c,j)   = tx(c,j)
             else                         ! soil layer
                 t_soisno(c,jprime) = tx(c,j)
             end if
          end if
       end do
    end do

    ! 8!) Sum energy content and total energy into lake for energy check. Any errors will be from the
    !     Tridiagonal solution.


    ! 9!) Phase change
    CALL PhaseChange_Lake (snl,h2osno,dz,dz_lake,                        & ! i
                           t_soisno,h2osoi_liq,h2osoi_ice,               & ! i&o
                           lake_icefrac,t_lake, snowdp,                  & ! i&o
                           qflx_snomelt,eflx_snomelt,imelt,              & ! o  
                           cv, cv_lake,                                  & ! i&o
                           lhabs)                                          ! o


    ! 9.5!) Second energy check and water check.  Now check energy balance before and after phase
    !       change, considering the possibility of changed heat capacity during phase change, by
    !       using initial heat capacity in the first step, final heat capacity in the second step,
    !       and differences from tfrz only to avoid enthalpy correction for (cpliq-cpice)*melt*tfrz.
    !       Also check soil water sum.

    ! 10!) Convective mixing: make sure fracice*dz is conserved, heat content c*dz*T is conserved, and
    ! all ice ends up at the top. Done over all lakes even if frozen.
    ! Either an unstable density profile or ice in a layer below an incompletely frozen layer will trigger.

    !Recalculate density
    do j = 1, nlevlake
! Here Original Concurent
       do fc = 1, num_shlakec
          c         = filter_shlakec(fc)
          rhow(c,j) = (1._r8 - lake_icefrac(c,j)) * &
                     !1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,j)-277._r8))**1.68_r8 ) &
                     !--- YangW Updated : Using the passed-in parameter tdmax to replace 277 ---!
                      1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,j)-tdmax))**1.68_r8 ) &
                      + lake_icefrac(c,j)*denice
       end do
    end do

    do j = 1, nlevlake-1
! Here Original Concurent
       do fc = 1, num_shlakec
          c        = filter_shlakec(fc)
          qav(c)   = 0._r8
          nav(c)   = 0._r8
          iceav(c) = 0._r8
       end do

       do i = 1, j+1
! Here Original Concurent
          do fc = 1, num_shlakec
             c = filter_shlakec(fc)
             if (rhow(c,j) > rhow(c,j+1) .or. &
                (lake_icefrac(c,j) < 1._r8 .and. lake_icefrac(c,j+1) > 0._r8) ) then
                qav(c)   = qav(c) + dz_lake(c,i)*(t_lake(c,i)-tfrz) * & 
                           ((1._r8 - lake_icefrac(c,i))*cwat + lake_icefrac(c,i)*cice_eff)
               !tav(c)   = tav(c) + t_lake(c,i)*dz_lake(c,i)
                iceav(c) = iceav(c) + lake_icefrac(c,i)*dz_lake(c,i)
                nav(c)   = nav(c) + dz_lake(c,i)
             end if
          end do
       end do

! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)
          if (rhow(c,j) > rhow(c,j+1) .or. &
             (lake_icefrac(c,j) < 1._r8 .and. lake_icefrac(c,j+1) > 0._r8) ) then
             qav(c)   = qav(c)/nav(c)
             iceav(c) = iceav(c)/nav(c)
             !If the average temperature is above freezing, put the extra energy into the water.
             !If it is below freezing, take it away from the ice.
             if (qav(c) > 0._r8) then
                tav_froz(c) = 0._r8 ! Celsius
                tav_unfr(c) = qav(c) / ((1._r8 - iceav(c))*cwat)
             else if (qav(c) < 0._r8) then
                tav_froz(c) = qav(c) / (iceav(c)*cice_eff)
                tav_unfr(c) = 0._r8 ! Celsius
             else
                tav_froz(c) = 0._r8
                tav_unfr(c) = 0._r8
             end if
          end if
       end do

       do i = 1, j+1
! Here Original Concurent
          do fc = 1, num_shlakec
             c = filter_shlakec(fc)
             if (nav(c) > 0._r8) then
                !Put all the ice at the top.!
                !If the average temperature is above freezing, put the extra energy into the water.
                !If it is below freezing, take it away from the ice.
                !For the layer with both ice & water, be careful to use the average temperature
                !that preserves the correct total heat content given what the heat capacity of that
                !layer will actually be.
                if (i == 1) zsum(c) = 0._r8
                if ((zsum(c)+dz_lake(c,i))/nav(c) <= iceav(c)) then
                    lake_icefrac(c,i) = 1._r8
                    t_lake(c,i)       = tav_froz(c) + tfrz
                else if (zsum(c)/nav(c) < iceav(c)) then
                   lake_icefrac(c,i) = (iceav(c)*nav(c) - zsum(c)) / dz_lake(c,i)
                   ! Find average value that preserves correct heat content.
                   t_lake(c,i)       = ( lake_icefrac(c,i)*tav_froz(c)*cice_eff          &
                                        + (1._r8 - lake_icefrac(c,i))*tav_unfr(c)*cwat ) &
                                        / ( lake_icefrac(c,i)*cice_eff + (1-lake_icefrac(c,i))*cwat ) + tfrz
                else
                   lake_icefrac(c,i) = 0._r8
                   t_lake(c,i)       = tav_unfr(c) + tfrz
                end if
                zsum(c)              = zsum(c) + dz_lake(c,i)
                rhow(c,i)            = (1._r8 - lake_icefrac(c,i)) * & 
                                      !1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,i)-277._r8))**1.68_r8 ) &
                                      !--- YangW Updated : Using the passed-in parameter tdmax to replace 277 ---!
                                       1000._r8*( 1.0_r8 - 1.9549e-05_r8*(abs(t_lake(c,i)-tdmax))**1.68_r8 ) &
                                       + lake_icefrac(c,i)*denice
             end if
          end do
       end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!
    ! 11!) Re-evaluate thermal properties and sum energy content.
    ! For lake
    do j = 1, nlevlake
! Here Original Concurent
       do fc = 1, num_shlakec
          c            = filter_shlakec(fc)
          cv_lake(c,j) = dz_lake(c,j) * (cwat*(1._r8-lake_icefrac(c,j)) + cice_eff*lake_icefrac(c,j))
       end do
    end do

  !For snow / soil
  !call SoilThermProp_Lake(lbc, ubc, num_shlakec, filter_shlakec, tk, cv, tktopsoillay)
   CALL SoilThermProp_Lake (snl,dz,zi,z,t_soisno,h2osoi_liq,h2osoi_ice,    &
                            tk, cv, tktopsoillay)

    ! Do as above to sum energy content
    do j = 1, nlevlake
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)
         !ncvts(c) = ncvts(c) + cv_lake(c,j)*t_lake(c,j) &
          ncvts(c) = ncvts(c) + cv_lake(c,j)*(t_lake(c,j)-tfrz) &
                     + cfus*dz_lake(c,j)*(1._r8-lake_icefrac(c,j)) !&
                    !+ (cwat-cice_eff)*lake_icefrac(c)*tfrz*dz_lake(c,j) !enthalpy reconciliation term
          fin(c)   = fin(c) + phi(c,j)
       end do
    end do

    do j = -nlevsnow + 1, nlevsoil
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)
          if (j >= jtop(c)) then
             !ncvts(c) = ncvts(c) + cv(c,j)*t_soisno(c,j) &
              ncvts(c) = ncvts(c) + cv(c,j)*(t_soisno(c,j)-tfrz) &
                         + hfus*h2osoi_liq(c,j) !&
                        !+ (cpliq-cpice)*h2osoi_ice(c,j)*tfrz !enthalpy reconciliation term
             if (j == 1 .and. h2osno(c) > 0._r8 .and. j == jtop(c)) then
                 ncvts(c) = ncvts(c) - h2osno(c)*hfus
             end if
          end if
          if (j == 1) fin(c) = fin(c) + phi_soil(c)
       end do
    end do


    ! Check energy conservation.
    do fp = 1, num_shlakep
       p         = filter_shlakep(fp)
       c         = pcolumn(p)
       errsoi(c) = (ncvts(c)-ocvts(c)) / dtime - fin(c)
      !if (abs(errsoi(c)) < 0.10_r8) then ! else send to Balance Check and abort
       if (abs(errsoi(c)) < 10._r8) then  ! else send to Balance Check and abort
      !if (abs(errsoi(c)) < 1._r8) then   ! A Optional Judgement? YangW
          eflx_sh_tot(p)    = eflx_sh_tot(p) - errsoi(c)
          eflx_sh_grnd(p)   = eflx_sh_grnd(p) - errsoi(c)
          eflx_soil_grnd(p) = eflx_soil_grnd(p) + errsoi(c)
          eflx_gnet(p)      = eflx_gnet(p) + errsoi(c)
         !if (abs(errsoi(c)) > 1.e-3_r8) then
          if (abs(errsoi(c)) > 1.e-1_r8) then
              if (lakedebug_yw) then            ! YangW: This part is designed in the original model
                  write(message,*)'errsoi incorporated into sensible heat in ShalLakeTemperature: c, (W/m^2):', c, errsoi(c)
                  print *, message
              end if
            !CALL wrf_message(message)
          end if
          errsoi(c) = 0._r8
       end if
    end do
    ! This loop assumes only one point per column.

  end SUBROUTINE ShalLakeTemperature

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!BOP
!
! ROUTINE: SoilThermProp_Lake
!
! !INTERFACE:
  SUBROUTINE SoilThermProp_Lake (snl,dz,zi,z,t_soisno,h2osoi_liq,h2osoi_ice,    &
                                 tk, cv, tktopsoillay)

!
! !DESCRIPTION:
! Calculation of thermal conductivities and heat capacities of
! snow/soil layers
! (1) The volumetric heat capacity is calculated as a linear combination
!     in terms of the volumetric fraction of the constituent phases.
!
! (2) The thermal conductivity of soil is computed from the algorithm of
!     Johansen (as reported by Farouki 1981), and of snow is from the
!     formulation used in SNTHERM (Jordan 1991).
! The thermal conductivities at the interfaces between two neighboring
! layers (j, j+1) are derived from an assumption that the flux across
! the interface is equal to that from the node j to the interface and the
! flux from the interface to the node j+1.
!
! For lakes, the proper soil layers (not snow) should always be saturated.
!
! !USES:

    implicit none
!in
    integer , intent(in) :: snl(1)           ! number of snow layers
!    real(r8), intent(in) :: h2osno(1)        ! snow water (mm H2O)
   ! real(r8), intent(in) :: watsat(1,nlevsoil)      ! volumetric soil water at saturation (porosity)
   ! real(r8), intent(in) :: tksatu(1,nlevsoil)      ! thermal conductivity, saturated soil [W/m-K]
   ! real(r8), intent(in) :: tkmg(1,nlevsoil)        ! thermal conductivity, soil minerals  [W/m-K]
   ! real(r8), intent(in) :: tkdry(1,nlevsoil)       ! thermal conductivity, dry soil (W/m/Kelvin)
   ! real(r8), intent(in) :: csol(1,nlevsoil)        ! heat capacity, soil solids (J/m**3/Kelvin)
    real(r8), intent(in) :: dz(1,-nlevsnow+1:nlevsoil)          ! layer thickness (m)
    real(r8), intent(in) :: zi(1,-nlevsnow+0:nlevsoil)          ! interface level below a "z" level (m)
    real(r8), intent(in) :: z(1,-nlevsnow+1:nlevsoil)           ! layer depth (m)
    real(r8), intent(in) :: t_soisno(1,-nlevsnow+1:nlevsoil)    ! soil temperature (Kelvin)
    real(r8), intent(in) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)  ! liquid water (kg/m2)
    real(r8), intent(in) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)  ! ice lens (kg/m2)

!out
    real(r8), intent(out) :: cv(lbc:ubc,-nlevsnow+1:nlevsoil) ! heat capacity [J/(m2 K)]
    real(r8), intent(out) :: tk(lbc:ubc,-nlevsnow+1:nlevsoil) ! thermal conductivity [W/(m K)]
    real(r8), intent(out) :: tktopsoillay(lbc:ubc)          ! thermal conductivity [W/(m K)]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !CALLED FROM:
! subroutine ShalLakeTemperature in this module.
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/13/02, Peter Thornton: migrated to new data structures
! 7/01/03, Mariana Vertenstein: migrated to vector code
! 4/09, Zack Subin, adjustment for ShalLake code.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in scalars
!
!    integer , pointer :: clandunit(:)     ! column's landunit
!    integer , pointer :: ityplun(:)       ! landunit type
!
!EOP
! OTHER LOCAL VARIABLES:
    integer  :: l,c,j                     ! indices
    integer  :: fc                        ! lake filtered column indices
    real(r8) :: bw                        ! partial density of water (ice + liquid)
    real(r8) :: dksat                     ! thermal conductivity for saturated soil (j/(k s m))
    real(r8) :: dke                       ! kersten number
    real(r8) :: fl                        ! fraction of liquid or unfrozen water to total water
    real(r8) :: satw                      ! relative total water content of soil.
    real(r8) :: thk(lbc:ubc,-nlevsnow+1:nlevsoil) ! thermal conductivity of layer
    character*256 :: message 

! Thermal conductivity of soil from Farouki (1981)

    do j = -nlevsnow+1,nlevsoil
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)

          ! Only examine levels from 1->nlevsoil
          if (j >= 1) then
!             l = clandunit(c)
!             if (ityplun(l) /= istwet .AND. ityplun(l) /= istice) then
              ! This could be altered later for allowing this to be over glaciers.

          ! Soil should be saturated.
                satw = 1._r8
                   fl = h2osoi_liq(c,j)/(h2osoi_ice(c,j)+h2osoi_liq(c,j))
                   if (t_soisno(c,j) >= tfrz) then       ! Unfrozen soil
                      dke = max(0._r8, log10(satw) + 1.0_r8)
                      dksat = tksatu(c,j)
                   else                               ! Frozen soil
                      dke = satw
                      dksat = tkmg(c,j)*0.249_r8**(fl*watsat(c,j))*2.29_r8**watsat(c,j)
                   endif
                   thk(c,j) = dke*dksat + (1._r8-dke)*tkdry(c,j)
!             else
!                thk(c,j) = tkwat
!                if (t_soisno(c,j) < tfrz) thk(c,j) = tkice
!             endif
          endif

          ! Thermal conductivity of snow, which from Jordan (1991) pp. 18
          ! Only examine levels from snl(c)+1 -> 0 where snl(c) < 1
          if (snl(c)+1 < 1 .AND. (j >= snl(c)+1) .AND. (j <= 0)) then
             bw = (h2osoi_ice(c,j)+h2osoi_liq(c,j))/dz(c,j)
             thk(c,j) = tkairc + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkairc)
          end if

       end do
    end do

    ! Thermal conductivity at the layer interface

    ! Have to correct for the fact that bottom snow layer and top soil layer border lake.
    ! For the first case, the snow layer conductivity for the middle of the layer will be returned.
    ! Because the interfaces are below the soil layers, the conductivity for the top soil layer
    ! will have to be returned separately.
    do j = -nlevsnow+1,nlevsoil
! Here Original Concurent
       do fc = 1,num_shlakec
          c = filter_shlakec(fc)
          if (j >= snl(c)+1 .AND. j <= nlevsoil-1 .AND. j /= 0) then
             tk(c,j) = thk(c,j)*thk(c,j+1)*(z(c,j+1)-z(c,j)) &
                  /(thk(c,j)*(z(c,j+1)-zi(c,j))+thk(c,j+1)*(zi(c,j)-z(c,j)))
          else if (j == 0) then
             tk(c,j) = thk(c,j)
          else if (j == nlevsoil) then
             tk(c,j) = 0._r8
          end if
          ! For top soil layer.
          if (j == 1) tktopsoillay(c) = thk(c,j)
       end do
    end do

    ! Soil heat capacity, from de Vires (1963)

    do j = 1, nlevsoil
! Here Original Concurent
       do fc = 1,num_shlakec
          c = filter_shlakec(fc)
!          l = clandunit(c)
!          if (ityplun(l) /= istwet .AND. ityplun(l) /= istice) then
             cv(c,j) = csol(c,j)*(1-watsat(c,j))*dz(c,j) +   &
               (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
!          else
!             cv(c,j) = (h2osoi_ice(c,j)*cpice + h2osoi_liq(c,j)*cpliq)
!          endif
!          if (j == 1) then
!             if (snl(c)+1 == 1 .AND. h2osno(c) > 0._r8) then
!                cv(c,j) = cv(c,j) + cpice*h2osno(c)
!             end if
!          end if
       ! Won't worry about heat capacity for thin snow on lake with no snow layers.
       enddo
    end do

    ! Snow heat capacity

    do j = -nlevsnow+1,0
! Here Original Concurent
       do fc = 1,num_shlakec
          c = filter_shlakec(fc)
          if (snl(c)+1 < 1 .and. j >= snl(c)+1) then
             cv(c,j) = cpliq*h2osoi_liq(c,j) + cpice*h2osoi_ice(c,j)
          end if
       end do
    end do

  end SUBROUTINE SoilThermProp_Lake


!-----------------------------------------------------------------------
!BOP
!
! ROUTINE: PhaseChange_Lake
!
! !INTERFACE:
  SUBROUTINE PhaseChange_Lake (snl,h2osno,dz,dz_lake,                        & !i
                               t_soisno,h2osoi_liq,h2osoi_ice,               & !i&o
                               lake_icefrac,t_lake, snowdp,                  & !i&o
                               qflx_snomelt,eflx_snomelt,imelt,              & !o  
                               cv, cv_lake,                                  & !i&o
                               lhabs)                                          !o
!=============================================================================================
! !DESCRIPTION:
! Calculation of the phase change within snow, soil, & lake layers:
! (1) Check the conditions for which the phase change may take place,
!     i.e., the layer temperature is great than the freezing point
!     and the ice mass is not equal to zero (i.e. melting),
!     or the layer temperature is less than the freezing point
!     and the liquid water mass is greater than the allowable supercooled 
!    (i.e. freezing).
! (2) Assess the amount of phase change from the energy excess (or deficit)
!     after setting the layer temperature to freezing point, depending on
!     how much water or ice is available.
! (3) Re-adjust the ice and liquid mass, and the layer temperature: either to
!     the freezing point if enough water or ice is available to fully compensate,
!     or to a remaining temperature.
! The specific heats are assumed constant. Potential cycling errors resulting from
! this assumption will be trapped at the end of ShalLakeTemperature.
! !CALLED FROM:
! subroutine ShalLakeTemperature in this module
!
! !REVISION HISTORY:
! 04/2009 Zack Subin: Initial code
!==============================================================================================
! !USES:
!
    implicit none
!in: 
    integer , intent(in)    :: snl(1)                     ! number of snow layers
    real(r8), intent(inout) :: h2osno(1)                  ! snow water (mm H2O)
    real(r8), intent(in)    :: dz(1,-nlevsnow+1:nlevsoil) ! layer thickness (m)
    real(r8), intent(in)    :: dz_lake(1,nlevlake)        ! lake layer thickness (m)
    ! Needed in case snow height is less than critical value.
!inout: 
    real(r8), intent(inout) :: snowdp(1)                          ! snow height (m)
    real(r8), intent(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)   ! soil temperature (Kelvin)
    real(r8), intent(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil) ! liquid water (kg/m2)
    real(r8), intent(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil) ! ice lens (kg/m2)
    real(r8), intent(inout) :: lake_icefrac(1,nlevlake          ) ! mass fraction of lake layer that is frozen
    real(r8), intent(inout) :: t_lake(1,nlevlake)                 ! lake temperature (Kelvin)
!out: 
    real(r8), intent(out)   :: qflx_snomelt(1)                    ! snow melt (mm H2O /s)
    real(r8), intent(out)   :: eflx_snomelt(1)                    ! snow melt heat flux (W/m**2)
    integer, intent(out)    :: imelt(1,-nlevsnow+1:nlevsoil)      ! flag for melting (=1), freezing (=2), Not=0 (new)
                                                                  ! What's the sign of this? Is it just output?
    real(r8), intent(inout) :: cv(lbc:ubc,-nlevsnow+1:nlevsoil)   ! heat capacity [J/(m2 K)]
    real(r8), intent(inout) :: cv_lake (lbc:ubc,1:nlevlake)       ! heat capacity [J/(m2 K)]
    real(r8), intent(out)   :: lhabs(lbc:ubc)                     ! total per-column latent heat abs. (J/m^2)
! OTHER LOCAL VARIABLES:
    integer  :: j,c,g                              ! do loop index
    integer  :: fc                                 ! lake filtered column indices
   !real(r8) :: dtime                              ! land model time step (sec)
    real(r8) :: heatavail                          ! available energy for melting or freezing (J/m^2)
    real(r8) :: heatrem                            ! energy residual or loss after melting or freezing
    real(r8) :: melt                               ! actual melting (+) or freezing (-) [kg/m2]
    real(r8), parameter :: smallnumber = 1.e-7_r8  ! to prevent tiny residuals from rounding error
    logical  :: dophasechangeflag
!-----------------------------------------------------------------------
    !  dtime = get_step_size()
    ! Initialization

! Here Original Concurent
    do fc = 1,num_shlakec
       c = filter_shlakec(fc)
       qflx_snomelt(c) = 0._r8
       eflx_snomelt(c) = 0._r8
       lhabs(c)        = 0._r8
    end do

    do j = -nlevsnow+1,0
! Here Original Concurent
       do fc = 1,num_shlakec
          c = filter_shlakec(fc)
          if (j >= snl(c) + 1) imelt(c,j) = 0
       end do
    end do

    ! Check for case of snow without snow layers and top lake layer temp above freezing.

! Here Original Concurent
    do fc = 1,num_shlakec
       c = filter_shlakec(fc)
       if (snl(c) == 0 .and. h2osno(c) > 0._r8 .and. t_lake(c,1) > tfrz) then
          heatavail = (t_lake(c,1) - tfrz) * cv_lake(c,1)
          melt      = min(h2osno(c), heatavail/hfus)
          heatrem   = max(heatavail - melt*hfus, 0._r8)
                       !catch small negative value to keep t at tfrz
          t_lake(c,1)     = tfrz + heatrem/(cv_lake(c,1))
          snowdp(c)       = snowdp(c)*(1._r8 - melt/h2osno(c))
          h2osno(c)       = h2osno(c) - melt
          lhabs(c)        = lhabs(c) + melt*hfus
          qflx_snomelt(c) = qflx_snomelt(c) + melt
          ! Prevent tiny residuals
          if (h2osno(c) < smallnumber) h2osno(c) = 0._r8
          if (snowdp(c) < smallnumber) snowdp(c) = 0._r8
       end if
    end do

    ! Lake phase change

    do j = 1,nlevlake
! Here Original Concurent
       do fc = 1,num_shlakec
          c                 = filter_shlakec(fc)
          dophasechangeflag = .false.
          if (t_lake(c,j) > tfrz .and. lake_icefrac(c,j) > 0._r8) then ! melting
             dophasechangeflag = .true.
             heatavail         = (t_lake(c,j) - tfrz) * cv_lake(c,j)
             melt              = min(lake_icefrac(c,j)*denh2o*dz_lake(c,j), heatavail/hfus)
                                 !denh2o is used because layer thickness is not adjusted for freezing
             heatrem           = max(heatavail - melt*hfus, 0._r8)
                                 !catch small negative value to keep t at tfrz
          else if (t_lake(c,j) < tfrz .and. lake_icefrac(c,j) < 1._r8) then !freezing
             dophasechangeflag = .true.
             heatavail         = (t_lake(c,j) - tfrz) * cv_lake(c,j)
             melt              = max(-(1._r8-lake_icefrac(c,j))*denh2o*dz_lake(c,j), heatavail/hfus)
                                 !denh2o is used because layer thickness is not adjusted for freezing
             heatrem           = min(heatavail - melt*hfus, 0._r8)
                                 !catch small positive value to keep t at tfrz
          end if
          ! Update temperature and ice fraction.
          if (dophasechangeflag) then
             lake_icefrac(c,j) = lake_icefrac(c,j) - melt/(denh2o*dz_lake(c,j))
             lhabs(c)          = lhabs(c) + melt*hfus
          ! Update heat capacity
             cv_lake(c,j)      = cv_lake(c,j) + melt*(cpliq-cpice)
             t_lake(c,j)       = tfrz + heatrem/cv_lake(c,j)
             ! Prevent tiny residuals
             if (lake_icefrac(c,j) > 1._r8 - smallnumber) lake_icefrac(c,j) = 1._r8
             if (lake_icefrac(c,j) < smallnumber)         lake_icefrac(c,j) = 0._r8
          end if
       end do
    end do

    ! Snow & soil phase change
    do j = -nlevsnow+1,nlevsoil
! Here Original Concurent
       do fc = 1,num_shlakec
          c                 = filter_shlakec(fc)
          dophasechangeflag = .false.
          if (j >= snl(c) + 1) then
             if (t_soisno(c,j) > tfrz .and. h2osoi_ice(c,j) > 0._r8) then ! melting
                dophasechangeflag = .true.
                heatavail         = (t_soisno(c,j) - tfrz) * cv(c,j)
                melt              = min(h2osoi_ice(c,j), heatavail/hfus)
                heatrem           = max(heatavail - melt*hfus, 0._r8)
                                   !catch small negative value to keep t at tfrz
                if (j <= 0) then !snow
                   imelt(c,j)      = 1
                   qflx_snomelt(c) = qflx_snomelt(c) + melt
                end if
             else if (t_soisno(c,j) < tfrz .and. h2osoi_liq(c,j) > 0._r8) then !freezing
                dophasechangeflag = .true.
                heatavail         = (t_soisno(c,j) - tfrz) * cv(c,j)
                melt              = max(-h2osoi_liq(c,j), heatavail/hfus)
                heatrem           = min(heatavail - melt*hfus, 0._r8)
                                   !catch small positive value to keep t at tfrz
                if (j <= 0) then !snow
                   imelt(c,j)      = 2
                   qflx_snomelt(c) = qflx_snomelt(c) + melt
                   ! Does this works for both signs of melt in SnowHydrology? I think
                   ! qflx_snomelt(c) is just output.
                end if
             end if
             ! Update temperature and soil components.
             if (dophasechangeflag) then
                h2osoi_ice(c,j) = h2osoi_ice(c,j) - melt
                h2osoi_liq(c,j) = h2osoi_liq(c,j) + melt
                lhabs(c)        = lhabs(c) + melt*hfus
             ! Update heat capacity
                cv(c,j)         = cv(c,j) + melt*(cpliq-cpice)
                t_soisno(c,j)   = tfrz + heatrem/cv(c,j)
                ! Prevent tiny residuals
                if (h2osoi_ice(c,j) < smallnumber) h2osoi_ice(c,j) = 0._r8
                if (h2osoi_liq(c,j) < smallnumber) h2osoi_liq(c,j) = 0._r8
             end if
         end if
      end do
   end do

   ! Update eflx_snomelt(c)
! Here Original Concurent
    do fc = 1,num_shlakec
       c               = filter_shlakec(fc)
       eflx_snomelt(c) = qflx_snomelt(c)*hfus
    end do

   end SUBROUTINE PhaseChange_Lake


  SUBROUTINE ShalLakeHydrology(dz_lake,forc_rain,forc_snow,                      & !i
                               begwb,qflx_evap_tot,forc_t,do_capsnow,            &
                               t_grnd,qflx_evap_soi,                             &
                               qflx_snomelt,imelt,frac_iceold,                   & !i add by guhp
                               z,dz,zi,snl,h2osno,snowdp,lake_icefrac,t_lake,      & !i&o
                               endwb,snowage,snowice,snowliq,t_snow,             & !o
                               t_soisno,h2osoi_ice,h2osoi_liq,h2osoi_vol,        &
                               qflx_drain,qflx_surf,qflx_infl,qflx_qrgwl,        &
                               qcharge,qflx_prec_grnd,qflx_snowcap,              &
                               qflx_snowcap_col,qflx_snow_grnd_pft,              &
                               qflx_snow_grnd_col,qflx_rain_grnd,                &
                               qflx_evap_tot_col,soilalpha,zwt,fcov,             &
                               rootr_column,qflx_evap_grnd,qflx_sub_snow,        &
                               qflx_dew_snow,qflx_dew_grnd,qflx_rain_grnd_col)
                       
!==================================================================================
! !DESCRIPTION:
! Calculation of Shallow Lake Hydrology. Full hydrology of snow layers is
! done. However, there is no infiltration, and the water budget is balanced with 
! qflx_qrgwl. Lake water mass is kept constant. The soil is simply maintained at
! volumetric saturation if ice melting frees up pore space. Likewise, if the water
! portion alone at some point exceeds pore capacity, it is reduced. This is consistent
! with the possibility of initializing the soil layer with excess ice. The only
! real error with that is that the thermal conductivity will ignore the excess ice
! (and accompanying thickness change).
! 
! If snow layers are present over an unfrozen lake, and the top layer of the lake
! is capable of absorbing the latent heat without going below freezing, 
! the snow-water is runoff and the latent heat is subtracted from the lake.
!
! WARNING: This subroutine assumes lake columns have one and only one pft.
!
! Sequence is:
!  ShalLakeHydrology:
!    Do needed tasks from Hydrology1, Biogeophysics2, & top of Hydrology2.
!    -> SnowWater:             change of snow mass and snow water onto soil
!    -> SnowCompaction:        compaction of snow layers
!    -> CombineSnowLayers:     combine snow layers that are thinner than minimum
!    -> DivideSnowLayers:      subdivide snow layers that are thicker than maximum
!    Add water to soil if melting has left it with open pore space.
!    Cleanup and do water balance.
!    If snow layers are found above a lake with unfrozen top layer, whose top
!    layer has enough heat to melt all the snow ice without freezing, do so
!    and eliminate the snow layers.
!
! !REVISION HISTORY:
! Created by Zack Subin, 2009
!
!============================================================================================

! USES:
!
    implicit none

! in:

   ! integer , intent(in) :: clandunit(1)     ! column's landunit
   ! integer , intent(in) :: ityplun(1)       ! landunit type
   ! real(r8), intent(in) :: watsat(1,1:nlevsoil)      ! volumetric soil water at saturation (porosity)
    real(r8), intent(in) :: dz_lake(1,nlevlake)     ! layer thickness for lake (m)
    real(r8), intent(in) :: forc_rain(1)     ! rain rate [mm/s]
    real(r8), intent(in) :: forc_snow(1)     ! snow rate [mm/s]
    real(r8), intent(in) :: qflx_evap_tot(1) ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), intent(in) :: forc_t(1)        ! atmospheric temperature (Kelvin)
    logical , intent(in) :: do_capsnow(1)     ! true => do snow capping
    real(r8), intent(in) :: t_grnd(1)          ! ground temperature (Kelvin)
    real(r8), intent(in) :: qflx_evap_soi(1)   ! soil evaporation (mm H2O/s) (+ = to atm)
    real(r8), intent(in) :: qflx_snomelt(1)     !snow melt (mm H2O /s)
    integer,  intent(in) :: imelt(1,-nlevsnow+1:nlevsoil)        !flag for melting (=1), freezing (=2), Not=0

!inout:

    real(r8), intent(inout) :: begwb(1)         ! water mass begining of the time step

! inout: 

    
    real(r8), intent(inout) :: z(1,-nlevsnow+1:nlevsoil)           ! layer depth  (m)
    real(r8), intent(inout) :: dz(1,-nlevsnow+1:nlevsoil)          ! layer thickness depth (m)
    real(r8), intent(inout) :: zi(1,-nlevsnow+0:nlevsoil)          ! interface depth (m)
    integer , intent(inout) :: snl(1)           ! number of snow layers
    real(r8), intent(inout) :: h2osno(1)        ! snow water (mm H2O)
    real(r8), intent(inout) :: snowdp(1)        ! snow height (m)
    real(r8), intent(inout) :: lake_icefrac(1,nlevlake)  ! mass fraction of lake layer that is frozen
    real(r8), intent(inout) :: t_lake(1,nlevlake)        ! lake temperature (Kelvin)

    real(r8), intent(inout) :: frac_iceold(1,-nlevsnow+1:nlevsoil)      ! fraction of ice relative to the tot water
! out: 


    real(r8), intent(out) :: endwb(1)         ! water mass end of the time step
    real(r8), intent(out) :: snowage(1)       ! non dimensional snow age [-]
    real(r8), intent(out) :: snowice(1)       ! average snow ice lens
    real(r8), intent(out) :: snowliq(1)       ! average snow liquid water
    real(r8), intent(out) :: t_snow(1)        ! vertically averaged snow temperature
    real(r8), intent(out) :: t_soisno(1,-nlevsnow+1:nlevsoil)    ! snow temperature (Kelvin)
    real(r8), intent(out) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)  ! ice lens (kg/m2)
    real(r8), intent(out) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)  ! liquid water (kg/m2)
    real(r8), intent(out) :: h2osoi_vol(1,-nlevsnow+1:nlevsoil)  ! volumetric soil water (0<=h2osoi_vol<=watsat)[m3/m3]
    real(r8), intent(out) :: qflx_drain(1)    ! sub-surface runoff (mm H2O /s)
    real(r8), intent(out) :: qflx_surf(1)     ! surface runoff (mm H2O /s)
    real(r8), intent(out) :: qflx_infl(1)     ! infiltration (mm H2O /s)
    real(r8), intent(out) :: qflx_qrgwl(1)    ! qflx_surf at glaciers, wetlands, lakes
    real(r8), intent(out) :: qcharge(1)       ! aquifer recharge rate (mm/s)
    real(r8), intent(out) :: qflx_prec_grnd(1)     ! water onto ground including canopy runoff [kg/(m2 s)]
    real(r8), intent(out) :: qflx_snowcap(1)       ! excess precipitation due to snow capping (mm H2O /s) [+]
    real(r8), intent(out) :: qflx_snowcap_col(1)   ! excess precipitation due to snow capping (mm H2O /s) [+]
    real(r8), intent(out) :: qflx_snow_grnd_pft(1) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), intent(out) :: qflx_snow_grnd_col(1) ! snow on ground after interception (mm H2O/s) [+]
    real(r8), intent(out) :: qflx_rain_grnd(1)     ! rain on ground after interception (mm H2O/s) [+]
    real(r8), intent(out) :: qflx_evap_tot_col(1) !pft quantity averaged to the column (assuming one pft)
    real(r8) ,intent(out) :: soilalpha(1)     !factor that reduces ground saturated specific humidity (-)
    real(r8), intent(out) :: zwt(1)           !water table depth
    real(r8), intent(out) :: fcov(1)          !fractional area with water table at surface
    real(r8), intent(out) :: rootr_column(1,1:nlevsoil) !effective fraction of roots in each soil layer
    real(r8), intent(out) :: qflx_evap_grnd(1)  ! ground surface evaporation rate (mm H2O/s) [+]
    real(r8), intent(out) :: qflx_sub_snow(1)   ! sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), intent(out) :: qflx_dew_snow(1)   ! surface dew added to snow pack (mm H2O /s) [+]
    real(r8), intent(out) :: qflx_dew_grnd(1)   ! ground surface dew formation (mm H2O /s) [+]
    real(r8), intent(out) :: qflx_rain_grnd_col(1)   !rain on ground after interception (mm H2O/s) [+]

! Block of biogeochem currently not used.
    real(r8), pointer :: sucsat(:,:)      ! minimum soil suction (mm)
    real(r8), pointer :: bsw(:,:)         ! Clapp and Hornberger "b"
    real(r8), pointer :: bsw2(:,:)        ! Clapp and Hornberger "b" for CN code
    real(r8), pointer :: psisat(:,:)      ! soil water potential at saturation for CN code (MPa)
    real(r8), pointer :: vwcsat(:,:)      ! volumetric water content at saturation for CN code (m3/m3)
    real(r8), pointer :: wf(:)            ! soil water as frac. of whc for top 0.5 m
    real(r8), pointer :: soilpsi(:,:)     ! soil water potential in each soil layer (MPa)
    real(r8) :: psi,vwc,fsat               ! temporary variables for soilpsi calculation

! OTHER LOCAL VARIABLES:

    integer  :: p,fp,g,l,c,j,fc,jtop             ! indices
    integer  :: num_shlakesnowc                  ! number of column snow points
    integer  :: filter_shlakesnowc(ubc-lbc+1)    ! column filter for snow points
    integer  :: num_shlakenosnowc                ! number of column non-snow points
    integer  :: filter_shlakenosnowc(ubc-lbc+1)  ! column filter for non-snow points
!    real(r8) :: dtime                      ! land model time step (sec)
    integer  :: newnode                      ! flag when new snow node is set, (1=yes, 0=no)
    real(r8) :: dz_snowf                     ! layer thickness rate change due to precipitation [mm/s]
    real(r8) :: bifall                       ! bulk density of newly fallen dry snow [kg/m3]
    real(r8) :: fracsnow(lbp:ubp)            ! frac of precipitation that is snow
    real(r8) :: fracrain(lbp:ubp)            ! frac of precipitation that is rain
    real(r8) :: qflx_prec_grnd_snow(lbp:ubp) ! snow precipitation incident on ground [mm/s]
    real(r8) :: qflx_prec_grnd_rain(lbp:ubp) ! rain precipitation incident on ground [mm/s]
    real(r8) :: qflx_evap_soi_lim            ! temporary evap_soi limited by top snow layer content [mm/s]
    real(r8) :: h2osno_temp                  ! temporary h2osno [kg/m^2]
    real(r8), parameter :: snow_bd = 250._r8  !constant snow bulk density (only used in special case here) [kg/m^3]
    real(r8) :: sumsnowice(lbc:ubc)             ! sum of snow ice if snow layers found above unfrozen lake [kg/m&2]
    logical  :: unfrozen(lbc:ubc)            ! true if top lake layer is unfrozen with snow layers above
    real(r8) :: heatrem                      ! used in case above [J/m^2]
    real(r8) :: heatsum(lbc:ubc)             ! used in case above [J/m^2]
    real(r8) :: qflx_top_soil(1)     !net water input into soil from top (mm/s)
    character*256 :: message 

!-----------------------------------------------------------------------


    ! Determine step size

!    dtime = get_step_size()

    ! Add soil water to water balance.
    do j = 1, nlevsoil
! Here Original Concurent
      do fc = 1, num_shlakec
         c = filter_shlakec(fc)
         begwb(c) = begwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
      end do
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Do precipitation onto ground, etc., from Hydrology1.

! Here Original Concurent
    do fp = 1, num_shlakep
       p = filter_shlakep(fp)
       g = pgridcell(p)
!       l = plandunit(p)
       c = pcolumn(p)

       ! Precipitation onto ground (kg/(m2 s))
!       ! PET, 1/18/2005: Added new terms for mass balance correction
!       ! due to dynamic pft weight shifting (column-level h2ocan_loss)
!       ! Because the fractionation between rain and snow is indeterminate if
!       ! rain + snow = 0, I am adding this very small flux only to the rain
!       ! components.
       ! Not relevant unless PFTs are added to lake later.
!       if (frac_veg_nosno(p) == 0) then
          qflx_prec_grnd_snow(p) = forc_snow(g)
          qflx_prec_grnd_rain(p) = forc_rain(g) !+ h2ocan_loss(c)
!       else
!          qflx_prec_grnd_snow(p) = qflx_through_snow(p) + (qflx_candrip(p) * fracsnow(p))
!          qflx_prec_grnd_rain(p) = qflx_through_rain(p) + (qflx_candrip(p) * fracrain(p)) + h2ocan_loss(c)
!       end if
       qflx_prec_grnd(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)

       if (do_capsnow(c)) then
          qflx_snowcap(p) = qflx_prec_grnd_snow(p) + qflx_prec_grnd_rain(p)
          qflx_snow_grnd_pft(p) = 0._r8
          qflx_rain_grnd(p) = 0._r8
       else
          qflx_snowcap(p) = 0._r8
          qflx_snow_grnd_pft(p) = qflx_prec_grnd_snow(p)           ! ice onto ground (mm/s)
          qflx_rain_grnd(p)     = qflx_prec_grnd_rain(p)           ! liquid water onto ground (mm/s)
       end if
       ! Assuming one PFT; needed for below
       qflx_snow_grnd_col(c) = qflx_snow_grnd_pft(p)
       qflx_rain_grnd_col(c) = qflx_rain_grnd(p)

    end do ! (end pft loop)

    ! Determine snow height and snow water

! Here Original Concurent
    do fc = 1, num_shlakec
       c = filter_shlakec(fc)
!       l = clandunit(c)
       g = cgridcell(c)

       ! Use Alta relationship, Anderson(1976); LaChapelle(1961),
       ! U.S.Department of Agriculture Forest Service, Project F,
       ! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

       if (do_capsnow(c)) then
          dz_snowf = 0._r8
       else
          if (forc_t(g) > tfrz + 2._r8) then
             bifall=50._r8 + 1.7_r8*(17.0_r8)**1.5_r8
          else if (forc_t(g) > tfrz - 15._r8) then
             bifall=50._r8 + 1.7_r8*(forc_t(g) - tfrz + 15._r8)**1.5_r8
          else
             bifall=50._r8
          end if
          dz_snowf = qflx_snow_grnd_col(c)/bifall
          snowdp(c) = snowdp(c) + dz_snowf*dtime
          h2osno(c) = h2osno(c) + qflx_snow_grnd_col(c)*dtime  ! snow water equivalent (mm)
       end if

!       if (itype(l)==istwet .and. t_grnd(c)>tfrz) then
!          h2osno(c)=0._r8
!          snowdp(c)=0._r8
!          snowage(c)=0._r8
!       end if
       ! Take care of this later in function.

       ! When the snow accumulation exceeds 10 mm, initialize snow layer
       ! Currently, the water temperature for the precipitation is simply set
       ! as the surface air temperature

       newnode = 0    ! flag for when snow node will be initialized
       if (snl(c) == 0 .and. qflx_snow_grnd_col(c) > 0.0_r8 .and. snowdp(c) >= 0.01_r8) then
          newnode = 1
          snl(c) = -1
          dz(c,0) = snowdp(c)                       ! meter
          z(c,0) = -0.5_r8*dz(c,0)
          zi(c,-1) = -dz(c,0)
          snowage(c) = 0._r8                        ! snow age
          t_soisno(c,0) = min(tfrz, forc_t(g))      ! K
          h2osoi_ice(c,0) = h2osno(c)               ! kg/m2
          h2osoi_liq(c,0) = 0._r8                   ! kg/m2
          frac_iceold(c,0) = 1._r8
       end if

       ! The change of ice partial density of surface node due to precipitation.
       ! Only ice part of snowfall is added here, the liquid part will be added
       ! later.

       if (snl(c) < 0 .and. newnode == 0) then
          h2osoi_ice(c,snl(c)+1) = h2osoi_ice(c,snl(c)+1)+dtime*qflx_snow_grnd_col(c)
          dz(c,snl(c)+1) = dz(c,snl(c)+1)+dz_snowf*dtime
       end if

    end do

    ! Calculate sublimation and dew, adapted from HydrologyLake and Biogeophysics2.

! Here Original Concurent
    do fp = 1,num_shlakep
       p = filter_shlakep(fp)
       c = pcolumn(p)
       jtop = snl(c)+1

       ! Use column variables here
       qflx_evap_grnd(c) = 0._r8
       qflx_sub_snow(c) = 0._r8
       qflx_dew_snow(c) = 0._r8
       qflx_dew_grnd(c) = 0._r8

       if (jtop <= 0) then ! snow layers
          j = jtop
          ! Assign ground evaporation to sublimation from soil ice or to dew
          ! on snow or ground

          if (qflx_evap_soi(p) >= 0._r8) then
          ! for evaporation partitioning between liquid evap and ice sublimation, 
          ! use the ratio of liquid to (liquid+ice) in the top layer to determine split
          ! Since we're not limiting evap over lakes, but still can't remove more from top
          ! snow layer than there is there, create temp. limited evap_soi.
             qflx_evap_soi_lim = min(qflx_evap_soi(p), (h2osoi_liq(c,j)+h2osoi_ice(c,j))/dtime)
             if ((h2osoi_liq(c,j)+h2osoi_ice(c,j)) > 0._r8) then
                qflx_evap_grnd(c) = max(qflx_evap_soi_lim*(h2osoi_liq(c,j)/(h2osoi_liq(c,j)+h2osoi_ice(c,j))), 0._r8)
             else
                qflx_evap_grnd(c) = 0._r8
             end if
             qflx_sub_snow(c) = qflx_evap_soi_lim - qflx_evap_grnd(c)     
          else
             if (t_grnd(c) < tfrz) then
                qflx_dew_snow(c) = abs(qflx_evap_soi(p))
             else
                qflx_dew_grnd(c) = abs(qflx_evap_soi(p))
             end if
          end if
          ! Update the pft-level qflx_snowcap
          ! This was moved in from Hydrology2 to keep all pft-level
          ! calculations out of Hydrology2
          if (do_capsnow(c)) qflx_snowcap(p) = qflx_snowcap(p) + qflx_dew_snow(c) + qflx_dew_grnd(c)

       else ! No snow layers: do as in HydrologyLake but with actual clmtype variables
          if (qflx_evap_soi(p) >= 0._r8) then
             ! Sublimation: do not allow for more sublimation than there is snow
             ! after melt.  Remaining surface evaporation used for infiltration.
             qflx_sub_snow(c) = min(qflx_evap_soi(p), h2osno(c)/dtime)
             qflx_evap_grnd(c) = qflx_evap_soi(p) - qflx_sub_snow(c)
          else
             if (t_grnd(c) < tfrz-0.1_r8) then
                qflx_dew_snow(c) = abs(qflx_evap_soi(p))
             else
                qflx_dew_grnd(c) = abs(qflx_evap_soi(p))
             end if
          end if

          ! Update snow pack for dew & sub.
          h2osno_temp = h2osno(c)
          if (do_capsnow(c)) then
             h2osno(c) = h2osno(c) - qflx_sub_snow(c)*dtime
             qflx_snowcap(p) = qflx_snowcap(p) + qflx_dew_snow(c) + qflx_dew_grnd(c)
          else
             h2osno(c) = h2osno(c) + (-qflx_sub_snow(c)+qflx_dew_snow(c))*dtime
          end if
          if (h2osno_temp > 0._r8) then
             snowdp(c) = snowdp(c) * h2osno(c) / h2osno_temp
          else
             snowdp(c) = h2osno(c)/snow_bd !Assume a constant snow bulk density = 250.
          end if

          h2osno(c) = max(h2osno(c), 0._r8)

       end if

    qflx_snowcap_col(c) = qflx_snowcap(p)

    end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Determine initial snow/no-snow filters (will be modified possibly by
    ! routines CombineSnowLayers and DivideSnowLayers below

    CALL BuildSnowFilter(lbc, ubc, num_shlakec, filter_shlakec,snl,       &            !i
         num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc) !o

    ! Determine the change of snow mass and the snow water onto soil

    CALL SnowWater(lbc, ubc, num_shlakesnowc, filter_shlakesnowc,         & !i 
                   num_shlakenosnowc, filter_shlakenosnowc,               & !i 
                   snl,do_capsnow,qflx_snomelt,qflx_rain_grnd,            & !i 
                   qflx_sub_snow,qflx_evap_grnd,                          & !i   
                   qflx_dew_snow,qflx_dew_grnd,dz,                        & !i   
                   h2osoi_ice,h2osoi_liq,                                 & !i&o 
                   qflx_top_soil)                                           !o                        


    ! Determine soil hydrology
    ! Here this consists only of making sure that soil is saturated even as it melts and 10%
    ! of pore space opens up. Conversely, if excess ice is melting and the liquid water exceeds the
    ! saturation value, then remove water.

    do j = 1,nlevsoil
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)

          if (h2osoi_vol(c,j) < watsat(c,j)) then
             h2osoi_liq(c,j) = (watsat(c,j)*dz(c,j) - h2osoi_ice(c,j)/denice)*denh2o
          ! h2osoi_vol will be updated below, and this water addition will come from qflx_qrgwl
          else if (h2osoi_liq(c,j) > watsat(c,j)*denh2o*dz(c,j)) then
             h2osoi_liq(c,j) = watsat(c,j)*denh2o*dz(c,j)
          end if

       end do
    end do
!!!!!!!!!!

!    if (.not. is_perpetual()) then
    if (1==1) then

       ! Natural compaction and metamorphosis.

       CALL SnowCompaction(lbc, ubc, num_shlakesnowc, filter_shlakesnowc,   &!i
                           snl,imelt,frac_iceold,t_soisno,                  &!i
                           h2osoi_ice,h2osoi_liq,                           &!i
                           dz)                                               !&o

       ! Combine thin snow elements

       CALL CombineSnowLayers(lbc, ubc,                            & !i
                              num_shlakesnowc, filter_shlakesnowc, & !i&o
                              snl,h2osno,snowdp,dz,zi,             & !i&o
                              t_soisno,h2osoi_ice,h2osoi_liq,      & !i&o
                              z)  !o                              


       ! Divide thick snow elements

       CALL DivideSnowLayers(lbc, ubc,                             & !i
                             num_shlakesnowc, filter_shlakesnowc,  & !i&o
                             snl,dz,zi,t_soisno,                   & !i&o
                             h2osoi_ice,h2osoi_liq,                & !i&o
                             z)  !o


    else

       do fc = 1, num_shlakesnowc
          c = filter_shlakesnowc(fc)
          h2osno(c) = 0._r8
       end do
       do j = -nlevsnow+1,0
          do fc = 1, num_shlakesnowc
             c = filter_shlakesnowc(fc)
             if (j >= snl(c)+1) then
                h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
             end if
          end do
       end do

    end if

    ! Check for snow layers above lake with unfrozen top layer.  Mechanically,
    ! the snow will fall into the lake and melt or turn to ice.  If the top layer has
    ! sufficient heat to melt the snow without freezing, then that will be done.
    ! Otherwise, the top layer will undergo freezing, but only if the top layer will
    ! not freeze completely.  Otherwise, let the snow layers persist and melt by diffusion.
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)

          if (t_lake(c,1) > tfrz .and. lake_icefrac(c,1) == 0._r8 .and. snl(c) < 0) then
             unfrozen(c) = .true.
          else
             unfrozen(c) = .false.
          end if
       end do

    do j = -nlevsnow+1,0
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)

          if (unfrozen(c)) then
             if (j == -nlevsnow+1) then
                sumsnowice(c) = 0._r8
                heatsum(c) = 0._r8
             end if
             if (j >= snl(c)+1) then
                sumsnowice(c) = sumsnowice(c) + h2osoi_ice(c,j)
                heatsum(c) = heatsum(c) + h2osoi_ice(c,j)*cpice*(tfrz - t_soisno(c,j)) &
                           + h2osoi_liq(c,j)*cpliq*(tfrz - t_soisno(c,j))
             end if
          end if
       end do
    end do

! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)

          if (unfrozen(c)) then
             heatsum(c) = heatsum(c) + sumsnowice(c)*hfus
             heatrem = (t_lake(c,1) - tfrz)*cpliq*denh2o*dz_lake(c,1) - heatsum(c)

             if (heatrem + denh2o*dz_lake(c,1)*hfus > 0._r8) then            
                ! Remove snow and subtract the latent heat from the top layer.
                h2osno(c) = 0._r8
                snl(c) = 0
                ! The rest of the bookkeeping for the removed snow will be done below.
                if (heatrem > 0._r8) then ! simply subtract the heat from the layer
                   t_lake(c,1) = t_lake(c,1) - heatrem/(cpliq*denh2o*dz_lake(c,1))
                else !freeze part of the layer
                   t_lake(c,1) = tfrz
                   lake_icefrac(c,1) = -heatrem/(denh2o*dz_lake(c,1)*hfus)
                end if
             end if
          end if
       end do
!!!!!!!!!!!!

    ! Set snow age to zero if no snow

! Here Original Concurent
    do fc = 1, num_shlakesnowc
       c = filter_shlakesnowc(fc)
       if (snl(c) == 0) then
          snowage(c) = 0._r8
       end if
    end do

    ! Set empty snow layers to zero

    do j = -nlevsnow+1,0
! Here Original Concurent
       do fc = 1, num_shlakesnowc
          c = filter_shlakesnowc(fc)
          if (j <= snl(c) .and. snl(c) > -nlevsnow) then
             h2osoi_ice(c,j) = 0._r8
             h2osoi_liq(c,j) = 0._r8
             t_soisno(c,j) = 0._r8
             dz(c,j) = 0._r8
             z(c,j) = 0._r8
             zi(c,j-1) = 0._r8
          end if
       end do
    end do

    ! Build new snow filter

    CALL BuildSnowFilter(lbc, ubc, num_shlakec, filter_shlakec, snl,&   !i
         num_shlakesnowc, filter_shlakesnowc, num_shlakenosnowc, filter_shlakenosnowc) !o

    ! Vertically average t_soisno and sum of h2osoi_liq and h2osoi_ice
    ! over all snow layers for history output

! Here Original Concurent
    do fc = 1, num_shlakesnowc
       c = filter_shlakesnowc(fc)
       t_snow(c)  = 0._r8
       snowice(c) = 0._r8
       snowliq(c) = 0._r8
    end do
! Here Original Concurent
    do fc = 1, num_shlakenosnowc
       c = filter_shlakenosnowc(fc)
       t_snow(c)  = spval
       snowice(c) = spval
       snowliq(c) = spval
    end do

    do j = -nlevsnow+1, 0
! Here Original Concurent
       do fc = 1, num_shlakesnowc
          c = filter_shlakesnowc(fc)
          if (j >= snl(c)+1) then
             t_snow(c)  = t_snow(c) + t_soisno(c,j)
             snowice(c) = snowice(c) + h2osoi_ice(c,j)
             snowliq(c) = snowliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Determine ending water balance and volumetric soil water

! Here Original Concurent
    do fc = 1, num_shlakec
       
       c = filter_shlakec(fc)
       if (snl(c) < 0) t_snow(c) = t_snow(c)/abs(snl(c))
       endwb(c) = h2osno(c)
    end do

    do j = 1, nlevsoil
! Here Original Concurent
       do fc = 1, num_shlakec
          c = filter_shlakec(fc)
          endwb(c) = endwb(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) + h2osoi_ice(c,j)/(dz(c,j)*denice)
       end do
    end do


!!!!!!!!!!!!!
    ! Do history variables and set special landunit runoff (adapted from end of HydrologyLake)
! Here Original Concurent
    do fp = 1,num_shlakep
       p = filter_shlakep(fp)
       c = pcolumn(p)
       g = pgridcell(p)

       qflx_infl(c)      = 0._r8
       qflx_surf(c)      = 0._r8
       qflx_drain(c)     = 0._r8
       rootr_column(c,:) = spval
       soilalpha(c)      = spval
       zwt(c)            = spval
       fcov(c)           = spval
       qcharge(c)        = spval
!       h2osoi_vol(c,:)   = spval

       ! Insure water balance using qflx_qrgwl
       qflx_qrgwl(c)     = forc_rain(g) + forc_snow(g) - qflx_evap_tot(p) - (endwb(c)-begwb(c))/dtime

       ! The pft average must be done here for output to history tape
       qflx_evap_tot_col(c) = qflx_evap_tot(p)
    end do

!!!!!!!!!!!!!
!For now, bracket off the remaining biogeochem code.  May need to bring it back
!to do soil carbon and methane beneath lakes.

  end SUBROUTINE ShalLakeHydrology

  SUBROUTINE QSat (T, p, es, esdT, qs, qsdT)
!
! !DESCRIPTION:
! Computes saturation mixing ratio and the change in saturation
! mixing ratio with respect to temperature.
! Reference:  Polynomial approximations from:
!             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
!             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: T        ! temperature (K)
    real(r8), intent(in)  :: p        ! surface atmospheric pressure (pa)
    real(r8), intent(out) :: es       ! vapor pressure (pa)
    real(r8), intent(out) :: esdT     ! d(es)/d(T)
    real(r8), intent(out) :: qs       ! humidity (kg/kg)
    real(r8), intent(out) :: qsdT     ! d(qs)/d(T)
!
! !CALLED FROM:
! subroutine Biogeophysics1 in module Biogeophysics1Mod
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine CanopyFluxesMod CanopyFluxesMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8) :: T_limit
    real(r8) :: td,vp,vp1,vp2
!
! For water vapor (temperature range 0C-100C)
!
    real(r8), parameter :: a0 =  6.11213476
    real(r8), parameter :: a1 =  0.444007856
    real(r8), parameter :: a2 =  0.143064234e-01
    real(r8), parameter :: a3 =  0.264461437e-03
    real(r8), parameter :: a4 =  0.305903558e-05
    real(r8), parameter :: a5 =  0.196237241e-07
    real(r8), parameter :: a6 =  0.892344772e-10
    real(r8), parameter :: a7 = -0.373208410e-12
    real(r8), parameter :: a8 =  0.209339997e-15
!
! For derivative:water vapor
!
    real(r8), parameter :: b0 =  0.444017302
    real(r8), parameter :: b1 =  0.286064092e-01
    real(r8), parameter :: b2 =  0.794683137e-03
    real(r8), parameter :: b3 =  0.121211669e-04
    real(r8), parameter :: b4 =  0.103354611e-06
    real(r8), parameter :: b5 =  0.404125005e-09
    real(r8), parameter :: b6 = -0.788037859e-12
    real(r8), parameter :: b7 = -0.114596802e-13
    real(r8), parameter :: b8 =  0.381294516e-16
!
! For ice (temperature range -75C-0C)
!
    real(r8), parameter :: c0 =  6.11123516
    real(r8), parameter :: c1 =  0.503109514
    real(r8), parameter :: c2 =  0.188369801e-01
    real(r8), parameter :: c3 =  0.420547422e-03
    real(r8), parameter :: c4 =  0.614396778e-05
    real(r8), parameter :: c5 =  0.602780717e-07
    real(r8), parameter :: c6 =  0.387940929e-09
    real(r8), parameter :: c7 =  0.149436277e-11
    real(r8), parameter :: c8 =  0.262655803e-14
!
! For derivative:ice
!
    real(r8), parameter :: d0 =  0.503277922
    real(r8), parameter :: d1 =  0.377289173e-01
    real(r8), parameter :: d2 =  0.126801703e-02
    real(r8), parameter :: d3 =  0.249468427e-04
    real(r8), parameter :: d4 =  0.313703411e-06
    real(r8), parameter :: d5 =  0.257180651e-08
    real(r8), parameter :: d6 =  0.133268878e-10
    real(r8), parameter :: d7 =  0.394116744e-13
    real(r8), parameter :: d8 =  0.498070196e-16
!-----------------------------------------------------------------------

    T_limit = T - tfrz
    if (T_limit > 100.0) T_limit=100.0
    if (T_limit < -75.0) T_limit=-75.0

    td       = T_limit
    if (td >= 0.0) then
       es   = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
            + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
       esdT = b0 + td*(b1 + td*(b2 + td*(b3 + td*(b4 &
            + td*(b5 + td*(b6 + td*(b7 + td*b8)))))))
    else
       es   = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
            + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
       esdT = d0 + td*(d1 + td*(d2 + td*(d3 + td*(d4 &
            + td*(d5 + td*(d6 + td*(d7 + td*d8)))))))
    endif

    es    = es    * 100.            ! pa
    esdT  = esdT  * 100.            ! pa/K

    vp    = 1.0   / (p - 0.378*es)
    vp1   = 0.622 * vp
    vp2   = vp1   * vp

    qs    = es    * vp1             ! kg/kg
    qsdT  = esdT  * vp2 * p         ! 1 / K

  end SUBROUTINE QSat


  SUBROUTINE Tridiagonal (lbc, ubc, lbj, ubj, jtop, numf, filter, &
                          a, b, c, r, u)
!
! !DESCRIPTION:
! Tridiagonal matrix solution
!
! !USES:
  !  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)    :: lbc, ubc               ! lbinning and ubing column indices
    integer , intent(in)    :: lbj, ubj               ! lbinning and ubing level indices
    integer , intent(in)    :: jtop(lbc:ubc)          ! top level for each column
    integer , intent(in)    :: numf                   ! filter dimension
    integer , intent(in)    :: filter(1:numf)         ! filter
    real(r8), intent(in)    :: a(lbc:ubc, lbj:ubj)    ! "a" left off diagonal of tridiagonal matrix
    real(r8), intent(in)    :: b(lbc:ubc, lbj:ubj)    ! "b" diagonal column for tridiagonal matrix
    real(r8), intent(in)    :: c(lbc:ubc, lbj:ubj)    ! "c" right off diagonal tridiagonal matrix
    real(r8), intent(in)    :: r(lbc:ubc, lbj:ubj)    ! "r" forcing term of tridiagonal matrix
    real(r8), intent(inout) :: u(lbc:ubc, lbj:ubj)    ! solution
!
! !CALLED FROM:
! subroutine BiogeophysicsLake in module BiogeophysicsLakeMod
! subroutine SoilTemperature in module SoilTemperatureMod
! subroutine SoilWater in module HydrologyMod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!  1 July 2003: Mariana Vertenstein; modified for vectorization
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: j,ci,fc                   !indices
    real(r8) :: gam(lbc:ubc,lbj:ubj)      !temporary
    real(r8) :: bet(lbc:ubc)              !temporary
!-----------------------------------------------------------------------
! Solve the matrix

! Here Original Concurent
    do fc = 1,numf
       ci = filter(fc)
       bet(ci) = b(ci,jtop(ci))
    end do

    do j = lbj, ubj
! Here Original Prevector
! Here Original Concurent
       do fc = 1,numf
          ci = filter(fc)
          if (j >= jtop(ci)) then
             if (j == jtop(ci)) then
                u(ci,j) = r(ci,j) / bet(ci)
             else
                gam(ci,j) = c(ci,j-1) / bet(ci)
                bet(ci)   = b(ci,j) - a(ci,j) * gam(ci,j)
                u(ci,j)   = (r(ci,j) - a(ci,j)*u(ci,j-1)) / bet(ci)
             end if
          end if
       end do
    end do

!Cray X1 unroll directive used here as work-around for compiler issue 2003/10/20
! Here Original unroll
    do j = ubj-1,lbj,-1
! Here Original Prevector
! Here Original Concurent
       do fc = 1,numf
          ci = filter(fc)
          if (j >= jtop(ci)) then
             u(ci,j) = u(ci,j) - gam(ci,j+1) * u(ci,j+1)
          end if
       end do
    end do

  end SUBROUTINE Tridiagonal


  SUBROUTINE SnowWater(lbc, ubc, num_snowc, filter_snowc,           & ! i
                       num_nosnowc, filter_nosnowc,                 & ! i 
                       snl,do_capsnow,qflx_snomelt,qflx_rain_grnd,  & ! i
                       qflx_sub_snow,qflx_evap_grnd,                & ! i   
                       qflx_dew_snow,qflx_dew_grnd,dz,              & ! i   
                       h2osoi_ice,h2osoi_liq,                       & ! i&o 
                       qflx_top_soil)                                 ! o                        
!===============================================================================
! !DESCRIPTION:
! Evaluate the change of snow mass and the snow water onto soil.
! Water flow within snow is computed by an explicit and non-physical
! based scheme, which permits a part of liquid water over the holding
! capacity (a tentative value is used, i.e. equal to 0.033*porosity) to
! percolate into the underlying layer.  Except for cases where the
! porosity of one of the two neighboring layers is less than 0.05, zero
! flow is assumed. The water flow out of the bottom of the snow pack will
! participate as the input of the soil water and runoff.  This subroutine
! uses a filter for columns containing snow which must be constructed prior
! to being called.
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 15 November 2000: Mariana Vertenstein
! 2/26/02, Peter Thornton: Migrated to new data structures.
!=============================================================================
! !USES:
  !  use clmtype

    implicit none

!in:
    integer, intent(in) :: lbc, ubc                    ! column bounds
    integer, intent(in) :: num_snowc                   ! number of snow points in column filter
    integer, intent(in) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
    integer, intent(in) :: num_nosnowc                 ! number of non-snow points in column filter
    integer, intent(in) :: filter_nosnowc(ubc-lbc+1)   ! column filter for non-snow points

    integer , intent(in) :: snl(1)              !number of snow layers
    logical , intent(in) :: do_capsnow(1)       !true => do snow capping
    real(r8), intent(in) :: qflx_snomelt(1)     !snow melt (mm H2O /s)
    real(r8), intent(in) :: qflx_rain_grnd(1)   !rain on ground after interception (mm H2O/s) [+]
    real(r8), intent(in) :: qflx_sub_snow(1)    !sublimation rate from snow pack (mm H2O /s) [+]
    real(r8), intent(in) :: qflx_evap_grnd(1)   !ground surface evaporation rate (mm H2O/s) [+]
    real(r8), intent(in) :: qflx_dew_snow(1)    !surface dew added to snow pack (mm H2O /s) [+]
    real(r8), intent(in) :: qflx_dew_grnd(1)    !ground surface dew formation (mm H2O /s) [+]
    real(r8), intent(in) :: dz(1,-nlevsnow+1:nlevsoil)             !layer depth (m)


!inout: 

    real(r8), intent(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)     !ice lens (kg/m2)
    real(r8), intent(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)     !liquid water (kg/m2)

!out:

    real(r8), intent(out) :: qflx_top_soil(1)     !net water input into soil from top (mm/s)


! OTHER LOCAL VARIABLES:

    integer  :: c, j, fc                           !do loop/array indices
    real(r8) :: qin(lbc:ubc)                       !water flow into the elmement (mm/s)
    real(r8) :: qout(lbc:ubc)                      !water flow out of the elmement (mm/s)
    real(r8) :: wgdif                              !ice mass after minus sublimation
    real(r8) :: vol_liq(lbc:ubc,-nlevsnow+1:0)      !partial volume of liquid water in layer
    real(r8) :: vol_ice(lbc:ubc,-nlevsnow+1:0)      !partial volume of ice lens in layer
    real(r8) :: eff_porosity(lbc:ubc,-nlevsnow+1:0) !effective porosity = porosity - vol_ice
!-----------------------------------------------------------------------
    ! Renew the mass of ice lens (h2osoi_ice) and liquid (h2osoi_liq) in the
    ! surface snow layer resulting from sublimation (frost) / evaporation (condense)

! Here Original Concurent
    do fc = 1,num_snowc
       c = filter_snowc(fc)
       if (do_capsnow(c)) then
          wgdif = h2osoi_ice(c,snl(c)+1) - qflx_sub_snow(c)*dtime
          h2osoi_ice(c,snl(c)+1) = wgdif
          if (wgdif < 0.) then
             h2osoi_ice(c,snl(c)+1) = 0.
             h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
          end if
          h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) - qflx_evap_grnd(c) * dtime
       else
          wgdif = h2osoi_ice(c,snl(c)+1) + (qflx_dew_snow(c) - qflx_sub_snow(c)) * dtime
          h2osoi_ice(c,snl(c)+1) = wgdif
          if (wgdif < 0.) then
             h2osoi_ice(c,snl(c)+1) = 0.
             h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) + wgdif
          end if
          h2osoi_liq(c,snl(c)+1) = h2osoi_liq(c,snl(c)+1) +  &
               (qflx_rain_grnd(c) + qflx_dew_grnd(c) - qflx_evap_grnd(c)) * dtime
       end if
       h2osoi_liq(c,snl(c)+1) = max(0._r8, h2osoi_liq(c,snl(c)+1))
    end do

    ! Porosity and partial volume

    do j = -nlevsnow+1, 0
! Here Original Concurent
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             vol_ice(c,j) = min(1._r8, h2osoi_ice(c,j)/(dz(c,j)*denice))
             eff_porosity(c,j) = 1. - vol_ice(c,j)
             vol_liq(c,j) = min(eff_porosity(c,j),h2osoi_liq(c,j)/(dz(c,j)*denh2o))
          end if
       end do
    end do

    ! Capillary forces within snow are usually two or more orders of magnitude
    ! less than those of gravity. Only gravity terms are considered.
    ! the genernal expression for water flow is "K * ss**3", however,
    ! no effective parameterization for "K".  Thus, a very simple consideration
    ! (not physically based) is introduced:
    ! when the liquid water of layer exceeds the layer's holding
    ! capacity, the excess meltwater adds to the underlying neighbor layer.

    qin(:) = 0._r8

    do j = -nlevsnow+1, 0
! Here Original Concurent
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             h2osoi_liq(c,j) = h2osoi_liq(c,j) + qin(c)
             if (j <= -1) then
                ! No runoff over snow surface, just ponding on surface
                if (eff_porosity(c,j) < wimp .OR. eff_porosity(c,j+1) < wimp) then
                   qout(c) = 0._r8
                else
                   qout(c) = max(0._r8,(vol_liq(c,j)-ssi*eff_porosity(c,j))*dz(c,j))
                   qout(c) = min(qout(c),(1.-vol_ice(c,j+1)-vol_liq(c,j+1))*dz(c,j+1))
                end if
             else
                qout(c) = max(0._r8,(vol_liq(c,j) - ssi*eff_porosity(c,j))*dz(c,j))
             end if
             qout(c) = qout(c)*1000.
             h2osoi_liq(c,j) = h2osoi_liq(c,j) - qout(c)
             qin(c) = qout(c)
          end if
       end do
    end do

! Here Original Concurent
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       ! Qout from snow bottom
       qflx_top_soil(c) = qout(c) / dtime
    end do

! Here Original Concurent
    do fc = 1, num_nosnowc
       c = filter_nosnowc(fc)
       qflx_top_soil(c) = qflx_rain_grnd(c) + qflx_snomelt(c)
    end do

  end SUBROUTINE SnowWater

  SUBROUTINE SnowCompaction(lbc, ubc, num_snowc, filter_snowc,   &!i  
                           snl,imelt,frac_iceold,t_soisno,                  &!i  
                           h2osoi_ice,h2osoi_liq,                           &!i  
                           dz)                                               !i&o   


!================================================================================
! !DESCRIPTION:
! Determine the change in snow layer thickness due to compaction and
! settling.
! Three metamorphisms of changing snow characteristics are implemented,
! i.e., destructive, overburden, and melt. The treatments of the former
! two are from SNTHERM.89 and SNTHERM.99 (1991, 1999). The contribution
! due to melt metamorphism is simply taken as a ratio of snow ice
! fraction after the melting versus before the melting.
!
! CALLED FROM:
! subroutine Hydrology2 in module Hydrology2Mod
!
! REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures
!==============================================================================
! USES:
  !  use clmtype
!
! !ARGUMENTS:
    implicit none

!in:
    integer, intent(in) :: lbc, ubc                ! column bounds
    integer, intent(in) :: num_snowc               ! number of column snow points in column filter
    integer, intent(in) :: filter_snowc(ubc-lbc+1) ! column filter for snow points
    integer,  intent(in) :: snl(1)             !number of snow layers
    integer,  intent(in) :: imelt(1,-nlevsnow+1:nlevsoil)        !flag for melting (=1), freezing (=2), Not=0
    real(r8), intent(in) :: frac_iceold(1,-nlevsnow+1:nlevsoil)  !fraction of ice relative to the tot water
    real(r8), intent(in) :: t_soisno(1,-nlevsnow+1:nlevsoil)     !soil temperature (Kelvin)
    real(r8), intent(in) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)   !ice lens (kg/m2)
    real(r8), intent(in) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)   !liquid water (kg/m2)

!inout:

    real(r8), intent(inout) :: dz(1,-nlevsnow+1:nlevsoil)           !layer depth (m)

! OTHER LOCAL VARIABLES:

    integer :: j, c, fc                   ! indices
    real(r8), parameter :: c2 = 23.e-3    ! [m3/kg]
    real(r8), parameter :: c3 = 2.777e-6  ! [1/s]
    real(r8), parameter :: c4 = 0.04      ! [1/K]
    real(r8), parameter :: c5 = 2.0       !
    real(r8), parameter :: dm = 100.0     ! Upper Limit on Destructive Metamorphism Compaction [kg/m3]
    real(r8), parameter :: eta0 = 9.e+5   ! The Viscosity Coefficient Eta0 [kg-s/m2]
    real(r8) :: burden(lbc:ubc) ! pressure of overlying snow [kg/m2]
    real(r8) :: ddz1   ! Rate of settling of snowpack due to destructive metamorphism.
    real(r8) :: ddz2   ! Rate of compaction of snowpack due to overburden.
    real(r8) :: ddz3   ! Rate of compaction of snowpack due to melt [1/s]
    real(r8) :: dexpf  ! expf=exp(-c4*(273.15-t_soisno)).
    real(r8) :: fi     ! Fraction of ice relative to the total water content at current time step
    real(r8) :: td     ! t_soisno - tfrz [K]
    real(r8) :: pdzdtc ! Nodal rate of change in fractional-thickness due to compaction [fraction/s]
    real(r8) :: void   ! void (1 - vol_ice - vol_liq)
    real(r8) :: wx     ! water mass (ice+liquid) [kg/m2]
    real(r8) :: bi     ! partial density of ice [kg/m3]

!-----------------------------------------------------------------------


    ! Begin calculation - note that the following column loops are only invoked if snl(c) < 0

    burden(:) = 0._r8

    do j = -nlevsnow+1, 0
! Here Original Concurent
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then

             wx = h2osoi_ice(c,j) + h2osoi_liq(c,j)
             void = 1. - (h2osoi_ice(c,j)/denice + h2osoi_liq(c,j)/denh2o) / dz(c,j)

             ! Allow compaction only for non-saturated node and higher ice lens node.
             if (void > 0.001 .and. h2osoi_ice(c,j) > .1) then
                bi = h2osoi_ice(c,j) / dz(c,j)
                fi = h2osoi_ice(c,j) / wx
                td = tfrz-t_soisno(c,j)
                dexpf = exp(-c4*td)

                ! Settling as a result of destructive metamorphism

                ddz1 = -c3*dexpf
                if (bi > dm) ddz1 = ddz1*exp(-46.0e-3*(bi-dm))

                ! Liquid water term

                if (h2osoi_liq(c,j) > 0.01*dz(c,j)) ddz1=ddz1*c5

                ! Compaction due to overburden

                ddz2 = -burden(c)*exp(-0.08*td - c2*bi)/eta0

                ! Compaction occurring during melt

                if (imelt(c,j) == 1) then
                   ddz3 = - 1./dtime * max(0._r8,(frac_iceold(c,j) - fi)/frac_iceold(c,j))
                else
                   ddz3 = 0._r8
                end if

                ! Time rate of fractional change in dz (units of s-1)

                pdzdtc = ddz1 + ddz2 + ddz3

                ! The change in dz due to compaction

                dz(c,j) = dz(c,j) * (1.+pdzdtc*dtime)
             end if

             ! Pressure of overlying snow

             burden(c) = burden(c) + wx

          end if
       end do
    end do

  end SUBROUTINE SnowCompaction

  SUBROUTINE CombineSnowLayers(lbc, ubc,                            & !i
                              num_snowc, filter_snowc, & !i&o
                              snl,h2osno,snowdp,dz,zi,             & !i&o
                              t_soisno,h2osoi_ice,h2osoi_liq,      & !i&o
                              z)  !o
!==========================================================================
! !DESCRIPTION:
! Combine snow layers that are less than a minimum thickness or mass
! If the snow element thickness or mass is less than a prescribed minimum,
! then it is combined with a neighboring element.  The SUBROUTINE
! clm\_combo.f90 then executes the combination of mass and energy.
! !CALLED FROM:
! SUBROUTINE Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures.
!=========================================================================
! !USES:
  !  use clmtype
!
! !ARGUMENTS:
    implicit none
!in:
    integer, intent(in)    :: lbc, ubc                    ! column bounds
   ! integer, intent(in) :: clandunit(1)       !landunit index for each column
   ! integer, intent(in) :: ityplun(1)         !landunit type

!inout:
    integer, intent(inout) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(inout) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
    integer , intent(inout) :: snl(1)            !number of snow layers
    real(r8), intent(inout) :: h2osno(1)         !snow water (mm H2O)
    real(r8), intent(inout) :: snowdp(1)         !snow height (m)
    real(r8), intent(inout) :: dz(1,-nlevsnow+1:nlevsoil)           !layer depth (m)
    real(r8), intent(inout) :: zi(1,-nlevsnow+0:nlevsoil)           !interface level below a "z" level (m)
    real(r8), intent(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)     !soil temperature (Kelvin)
    real(r8), intent(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)   !ice lens (kg/m2)
    real(r8), intent(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)   !liquid water (kg/m2)

!out:

    real(r8), intent(out) :: z(1,-nlevsnow+1:nlevsoil)            !layer thickness (m)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer :: c, fc                 ! column indices
    integer :: i,k                   ! loop indices
    integer :: j,l                   ! node indices
    integer :: msn_old(lbc:ubc)      ! number of top snow layer
    integer :: mssi(lbc:ubc)         ! node index
    integer :: neibor                ! adjacent node selected for combination
    real(r8):: zwice(lbc:ubc)        ! total ice mass in snow
    real(r8):: zwliq (lbc:ubc)       ! total liquid water in snow
    real(r8):: dzmin(5)              ! minimum of top snow layer

    data dzmin /0.010, 0.015, 0.025, 0.055, 0.115/
!-----------------------------------------------------------------------

    ! Check the mass of ice lens of snow, when the total is less than a small value,
    ! combine it with the underlying neighbor.

! Here Original Concurent
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       msn_old(c) = snl(c)
    end do

    ! The following loop is NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)
   !    l = clandunit(c)                                                    
       do j = msn_old(c)+1,0
          if (h2osoi_ice(c,j) <= .1) then
           !  if (ityplun(l) == istsoil) then                                
           !     h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)        
           !     h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)       
           !  else if (ityplun(l) /= istsoil .and. j /= 0) then               
             h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) + h2osoi_liq(c,j)
             h2osoi_ice(c,j+1) = h2osoi_ice(c,j+1) + h2osoi_ice(c,j)
           !  end if 

             ! shift all elements above this down one.
             if (j > snl(c)+1 .and. snl(c) < -1) then
                do i = j, snl(c)+2, -1
                   t_soisno(c,i)   = t_soisno(c,i-1)
                   h2osoi_liq(c,i) = h2osoi_liq(c,i-1)
                   h2osoi_ice(c,i) = h2osoi_ice(c,i-1)
                   dz(c,i)         = dz(c,i-1)
                end do
             end if
             snl(c) = snl(c) + 1
          end if
       end do
    end do

! Here Original Concurent
    do fc = 1, num_snowc
       c = filter_snowc(fc)
       h2osno(c) = 0._r8
       snowdp(c) = 0._r8
       zwice(c)  = 0._r8
       zwliq(c)  = 0._r8
    end do

    do j = -nlevsnow+1,0
! Here Original Concurent
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             h2osno(c) = h2osno(c) + h2osoi_ice(c,j) + h2osoi_liq(c,j)
             snowdp(c) = snowdp(c) + dz(c,j)
             zwice(c)  = zwice(c) + h2osoi_ice(c,j)
             zwliq(c)  = zwliq(c) + h2osoi_liq(c,j)
          end if
       end do
    end do

    ! Check the snow depth - all snow gone
    ! The liquid water assumes ponding on soil surface.

! Here Original Concurent
    do fc = 1, num_snowc
       c = filter_snowc(fc)
      ! l = clandunit(c)                                         
       if (snowdp(c) < 0.01 .and. snowdp(c) > 0.) then
          snl(c) = 0
          h2osno(c) = zwice(c)
          if (h2osno(c) <= 0.) snowdp(c) = 0._r8
      !    if (ityplun(l) == istsoil) h2osoi_liq(c,1) = h2osoi_liq(c,1) + zwliq(c)    !change by guhp
       end if
    end do

    ! Check the snow depth - snow layers combined
    ! The following loop IS NOT VECTORIZED

    do fc = 1, num_snowc
       c = filter_snowc(fc)

       ! Two or more layers

       if (snl(c) < -1) then

          msn_old(c) = snl(c)
          mssi(c) = 1

          do i = msn_old(c)+1,0
             if (dz(c,i) < dzmin(mssi(c))) then

                if (i == snl(c)+1) then
                   ! If top node is removed, combine with bottom neighbor.
                   neibor = i + 1
                else if (i == 0) then
                   ! If the bottom neighbor is not snow, combine with the top neighbor.
                   neibor = i - 1
                else
                   ! If none of the above special cases apply, combine with the thinnest neighbor
                   neibor = i + 1
                   if ((dz(c,i-1)+dz(c,i)) < (dz(c,i+1)+dz(c,i))) neibor = i-1
                end if

                ! Node l and j are combined and stored as node j.
                if (neibor > i) then
                   j = neibor
                   l = i
                else
                   j = i
                   l = neibor
                end if

                CALL Combo (dz(c,j), h2osoi_liq(c,j), h2osoi_ice(c,j), &
                            t_soisno(c,j), dz(c,l), h2osoi_liq(c,l), h2osoi_ice(c,l), t_soisno(c,l) )

                ! Now shift all elements above this down one.
                if (j-1 > snl(c)+1) then
                   do k = j-1, snl(c)+2, -1
                      t_soisno(c,k) = t_soisno(c,k-1)
                      h2osoi_ice(c,k) = h2osoi_ice(c,k-1)
                      h2osoi_liq(c,k) = h2osoi_liq(c,k-1)
                      dz(c,k) = dz(c,k-1)
                   end do
                end if

                ! Decrease the number of snow layers
                snl(c) = snl(c) + 1
                if (snl(c) >= -1) EXIT

             else

                ! The layer thickness is greater than the prescribed minimum value
                mssi(c) = mssi(c) + 1

             end if
          end do

       end if

    end do

    ! Reset the node depth and the depth of layer interface

    do j = 0, -nlevsnow+1, -1
! Here Original Concurent
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c) + 1) then
             z(c,j) = zi(c,j) - 0.5*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

  end SUBROUTINE CombineSnowLayers

  SUBROUTINE DivideSnowLayers(lbc, ubc,                             & !i
                             num_snowc, filter_snowc,  & !i&o
                             snl,dz,zi,t_soisno,                   & !i&o
                             h2osoi_ice,h2osoi_liq,                & !i&o
                             z)  !o


!============================================================================
! !DESCRIPTION:
! Subdivides snow layers if they exceed their prescribed maximum thickness.
! !CALLED FROM:
! SUBROUTINE Hydrology2 in module Hydrology2Mod
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 2/28/02, Peter Thornton: Migrated to new data structures.
!============================================================================
! !USES:
 !   use clmtype
!
! !ARGUMENTS:
    implicit none

!in:
    integer, intent(in)    :: lbc, ubc                    ! column bounds

!inout:

    integer, intent(inout) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(inout) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
    integer , intent(inout) :: snl(1)            !number of snow layers
    real(r8), intent(inout) :: dz(1,-nlevsnow+1:nlevsoil)           !layer depth (m)
    real(r8), intent(inout) :: zi(1,-nlevsnow+0:nlevsoil)           !interface level below a "z" level (m)
    real(r8), intent(inout) :: t_soisno(1,-nlevsnow+1:nlevsoil)     !soil temperature (Kelvin)
    real(r8), intent(inout) :: h2osoi_ice(1,-nlevsnow+1:nlevsoil)   !ice lens (kg/m2)
    real(r8), intent(inout) :: h2osoi_liq(1,-nlevsnow+1:nlevsoil)   !liquid water (kg/m2)

!out: 

    real(r8), intent(out) :: z(1,-nlevsnow+1:nlevsoil)            !layer thickness (m)



! OTHER LOCAL VARIABLES:

    integer  :: j, c, fc               ! indices
    real(r8) :: drr                    ! thickness of the combined [m]
    integer  :: msno                   ! number of snow layer 1 (top) to msno (bottom)
    real(r8) :: dzsno(lbc:ubc,nlevsnow) ! Snow layer thickness [m]
    real(r8) :: swice(lbc:ubc,nlevsnow) ! Partial volume of ice [m3/m3]
    real(r8) :: swliq(lbc:ubc,nlevsnow) ! Partial volume of liquid water [m3/m3]
    real(r8) :: tsno(lbc:ubc ,nlevsnow) ! Nodel temperature [K]
    real(r8) :: zwice                  ! temporary
    real(r8) :: zwliq                  ! temporary
    real(r8) :: propor                 ! temporary
!-----------------------------------------------------------------------

    ! Begin calculation - note that the following column loops are only invoked
    ! for snow-covered columns

    do j = 1,nlevsnow
! Here Original Concurent
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j <= abs(snl(c))) then
             dzsno(c,j) = dz(c,j+snl(c))
             swice(c,j) = h2osoi_ice(c,j+snl(c))
             swliq(c,j) = h2osoi_liq(c,j+snl(c))
             tsno(c,j)  = t_soisno(c,j+snl(c))
          end if
       end do
    end do

! Here Original Concurent
    do fc = 1, num_snowc
       c = filter_snowc(fc)

       msno = abs(snl(c))

       if (msno == 1) then
          ! Specify a new snow layer
          if (dzsno(c,1) > 0.03) then
             msno = 2
             dzsno(c,1) = dzsno(c,1)/2.
             swice(c,1) = swice(c,1)/2.
             swliq(c,1) = swliq(c,1)/2.
             dzsno(c,2) = dzsno(c,1)
             swice(c,2) = swice(c,1)
             swliq(c,2) = swliq(c,1)
             tsno(c,2)  = tsno(c,1)
          end if
       end if

       if (msno > 1) then
          if (dzsno(c,1) > 0.02) then
             drr = dzsno(c,1) - 0.02
             propor = drr/dzsno(c,1)
             zwice = propor*swice(c,1)
             zwliq = propor*swliq(c,1)
             propor = 0.02/dzsno(c,1)
             swice(c,1) = propor*swice(c,1)
             swliq(c,1) = propor*swliq(c,1)
             dzsno(c,1) = 0.02

             CALL Combo (dzsno(c,2), swliq(c,2), swice(c,2), tsno(c,2), drr, &
                         zwliq, zwice, tsno(c,1))

             ! Subdivide a new layer
             if (msno <= 2 .and. dzsno(c,2) > 0.07) then
                msno = 3
                dzsno(c,2) = dzsno(c,2)/2.
                swice(c,2) = swice(c,2)/2.
                swliq(c,2) = swliq(c,2)/2.
                dzsno(c,3) = dzsno(c,2)
                swice(c,3) = swice(c,2)
                swliq(c,3) = swliq(c,2)
                tsno(c,3)  = tsno(c,2)
             end if
          end if
       end if

       if (msno > 2) then
          if (dzsno(c,2) > 0.05) then
             drr = dzsno(c,2) - 0.05
             propor = drr/dzsno(c,2)
             zwice = propor*swice(c,2)
             zwliq = propor*swliq(c,2)
             propor = 0.05/dzsno(c,2)
             swice(c,2) = propor*swice(c,2)
             swliq(c,2) = propor*swliq(c,2)
             dzsno(c,2) = 0.05

             CALL Combo (dzsno(c,3), swliq(c,3), swice(c,3), tsno(c,3), drr, &
                  zwliq, zwice, tsno(c,2))

             ! Subdivided a new layer
             if (msno <= 3 .and. dzsno(c,3) > 0.18) then
                msno =  4
                dzsno(c,3) = dzsno(c,3)/2.
                swice(c,3) = swice(c,3)/2.
                swliq(c,3) = swliq(c,3)/2.
                dzsno(c,4) = dzsno(c,3)
                swice(c,4) = swice(c,3)
                swliq(c,4) = swliq(c,3)
                tsno(c,4)  = tsno(c,3)
             end if
          end if
       end if

       if (msno > 3) then
          if (dzsno(c,3) > 0.11) then
             drr = dzsno(c,3) - 0.11
             propor = drr/dzsno(c,3)
             zwice = propor*swice(c,3)
             zwliq = propor*swliq(c,3)
             propor = 0.11/dzsno(c,3)
             swice(c,3) = propor*swice(c,3)
             swliq(c,3) = propor*swliq(c,3)
             dzsno(c,3) = 0.11

             CALL Combo (dzsno(c,4), swliq(c,4), swice(c,4), tsno(c,4), drr, &
                  zwliq, zwice, tsno(c,3))

             ! Subdivided a new layer
             if (msno <= 4 .and. dzsno(c,4) > 0.41) then
                msno = 5
                dzsno(c,4) = dzsno(c,4)/2.
                swice(c,4) = swice(c,4)/2.
                swliq(c,4) = swliq(c,4)/2.
                dzsno(c,5) = dzsno(c,4)
                swice(c,5) = swice(c,4)
                swliq(c,5) = swliq(c,4)
                tsno(c,5)  = tsno(c,4)
             end if
          end if
       end if

       if (msno > 4) then
          if (dzsno(c,4) > 0.23) then
             drr = dzsno(c,4) - 0.23
             propor = drr/dzsno(c,4)
             zwice = propor*swice(c,4)
             zwliq = propor*swliq(c,4)
             propor = 0.23/dzsno(c,4)
             swice(c,4) = propor*swice(c,4)
             swliq(c,4) = propor*swliq(c,4)
             dzsno(c,4) = 0.23

             CALL Combo (dzsno(c,5), swliq(c,5), swice(c,5), tsno(c,5), drr, &
                  zwliq, zwice, tsno(c,4))
          end if
       end if

       snl(c) = -msno

    end do

    do j = -nlevsnow+1,0
! Here Original Concurent
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             dz(c,j)         = dzsno(c,j-snl(c))
             h2osoi_ice(c,j) = swice(c,j-snl(c))
             h2osoi_liq(c,j) = swliq(c,j-snl(c))
             t_soisno(c,j)   = tsno(c,j-snl(c))
          end if
       end do
    end do

    do j = 0, -nlevsnow+1, -1
! Here Original Concurent
       do fc = 1, num_snowc
          c = filter_snowc(fc)
          if (j >= snl(c)+1) then
             z(c,j)    = zi(c,j) - 0.5*dz(c,j)
             zi(c,j-1) = zi(c,j) - dz(c,j)
          end if
       end do
    end do

  end SUBROUTINE DivideSnowLayers

  SUBROUTINE Combo(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
!
! !DESCRIPTION:
! Combines two elements and returns the following combined
! variables: dz, t, wliq, wice.
! The combined temperature is based on the equation:
! the sum of the enthalpies of the two elements =
! that of the combined element.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)    :: dz2   ! nodal thickness of 2 elements being combined [m]
    real(r8), intent(in)    :: wliq2 ! liquid water of element 2 [kg/m2]
    real(r8), intent(in)    :: wice2 ! ice of element 2 [kg/m2]
    real(r8), intent(in)    :: t2    ! nodal temperature of element 2 [K]
    real(r8), intent(inout) :: dz    ! nodal thickness of 1 elements being combined [m]
    real(r8), intent(inout) :: wliq  ! liquid water of element 1
    real(r8), intent(inout) :: wice  ! ice of element 1 [kg/m2]
    real(r8), intent(inout) :: t     ! nodel temperature of elment 1 [K]
!
! !CALLED FROM:
! SUBROUTINE CombineSnowLayers in this module
! SUBROUTINE DivideSnowLayers in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8) :: dzc   ! Total thickness of nodes 1 and 2 (dzc=dz+dz2).
    real(r8) :: wliqc ! Combined liquid water [kg/m2]
    real(r8) :: wicec ! Combined ice [kg/m2]
    real(r8) :: tc    ! Combined node temperature [K]
    real(r8) :: h     ! enthalpy of element 1 [J/m2]
    real(r8) :: h2    ! enthalpy of element 2 [J/m2]
    real(r8) :: hc    ! temporary
!-----------------------------------------------------------------------

    dzc = dz+dz2
    wicec = (wice+wice2)
    wliqc = (wliq+wliq2)
    h = (cpice*wice+cpliq*wliq) * (t-tfrz)+hfus*wliq
    h2= (cpice*wice2+cpliq*wliq2) * (t2-tfrz)+hfus*wliq2

    hc = h + h2
    if(hc < 0.)then
       tc = tfrz + hc/(cpice*wicec + cpliq*wliqc)
    else if (hc.le.hfus*wliqc) then
       tc = tfrz
    else
       tc = tfrz + (hc - hfus*wliqc) / (cpice*wicec + cpliq*wliqc)
    end if

    dz = dzc
    wice = wicec
    wliq = wliqc
    t = tc

  end SUBROUTINE Combo

  SUBROUTINE BuildSnowFilter(lbc, ubc, num_nolakec, filter_nolakec,snl, & !i
                             num_snowc, filter_snowc, &                   !o
                             num_nosnowc, filter_nosnowc)                 !o
!
! !DESCRIPTION:
! Constructs snow filter for use in vectorized loops for snow hydrology.
!
! !USES:
!    use clmtype
!
! !ARGUMENTS:
    implicit none
    integer, intent(in)  :: lbc, ubc                    ! column bounds
    integer, intent(in)  :: num_nolakec                 ! number of column non-lake points in column filter
    integer, intent(in)  :: filter_nolakec(ubc-lbc+1)   ! column filter for non-lake points
    integer, intent(in)  :: snl(1)                        ! number of snow layers
    integer, intent(out) :: num_snowc                   ! number of column snow points in column filter
    integer, intent(out) :: filter_snowc(ubc-lbc+1)     ! column filter for snow points
    integer, intent(out) :: num_nosnowc                 ! number of column non-snow points in column filter
    integer, intent(out) :: filter_nosnowc(ubc-lbc+1)   ! column filter for non-snow points
!
! !CALLED FROM:
! SUBROUTINE Hydrology2 in Hydrology2Mod
! SUBROUTINE CombineSnowLayers in this module
!
! !REVISION HISTORY:
! 2003 July 31: Forrest Hoffman
!
! !LOCAL VARIABLES:
! local pointers to implicit in arguments
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer  :: fc, c
!-----------------------------------------------------------------------


    ! Build snow/no-snow filters for other SUBROUTINEs

    num_snowc = 0
    num_nosnowc = 0
    do fc = 1, num_nolakec
       c = filter_nolakec(fc)
       if (snl(c) < 0) then
          num_snowc = num_snowc + 1
          filter_snowc(num_snowc) = c
       else
          num_nosnowc = num_nosnowc + 1
          filter_nosnowc(num_nosnowc) = c
       end if
    end do

  end SUBROUTINE BuildSnowFilter

SUBROUTINE FrictionVelocity(pgridcell,forc_hgt,forc_hgt_u,          & ! i
                            forc_hgt_t,forc_hgt_q,                  & ! i
                           !lbp, ubp, fncopy, fpcopy,               & ! i Special dummy argument declaration For this Subroutine YangW
                            lbp, ubp, fn, filterp,                  & ! i
                           !displa, z0mg, z0hg, z0qg,               & ! i Special dummy argument declaration For this Subroutine YangW
                            displa, z0m, z0h, z0q,                  & ! i
                            obu,  um,                               & ! i
                            ustar,temp1, temp2, temp12m, temp22m)     ! o
                            !--- YangW: Following Variables For dustmodel -> Removed ---!
                            !u10, fv, iter, ur, fm)                    ! i&o
!=============================================================================
! !DESCRIPTION:
! Calculation of the friction velocity, relation for potential
! temperature and humidity profiles of surface boundary layer.
! The scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
! Vol. 11, 2628-2644.
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 12/19/01, Peter Thornton
! Added arguments to eliminate passing clm derived type into this function.
! Created by Mariana Vertenstein
!============================================================================
! !USES:
! use clmtype
! use clm_atmlnd, only : clm_a2l
!
   implicit none
!in:
   integer , intent(in)  :: pgridcell(1)     ! pft's gridcell index
   real(r8), intent(in)  :: forc_hgt(1)      ! atmospheric reference height (m)
   real(r8), intent(in)  :: forc_hgt_u(1)    ! observational height of wind [m]
   real(r8), intent(in)  :: forc_hgt_t(1)    ! observational height of temperature [m]
   real(r8), intent(in)  :: forc_hgt_q(1)    ! observational height of humidity [m]
   integer , intent(in)  :: lbp, ubp         ! pft array bounds
   integer , intent(in)  :: fn               ! number of filtered pft elements
   integer , intent(in)  :: filterp(fn)      ! pft filter
   real(r8), intent(in)  :: displa(lbp:ubp)  ! displacement height (m)
   real(r8), intent(in)  :: z0m(lbp:ubp)     ! roughness length over vegetation, momentum [m]
   real(r8), intent(in)  :: z0h(lbp:ubp)     ! roughness length over vegetation, sensible heat [m]
   real(r8), intent(in)  :: z0q(lbp:ubp)     ! roughness length over vegetation, latent heat [m]
   real(r8), intent(in)  :: obu(lbp:ubp)     ! monin-obukhov length (m)
   real(r8), intent(in)  :: um(lbp:ubp)      ! wind speed including the stablity effect [m/s]
!out:
   real(r8), intent(out) :: ustar(lbp:ubp)   ! friction velocity [m/s]
   real(r8), intent(out) :: temp1(lbp:ubp)   ! relation for potential temperature profile
   real(r8), intent(out) :: temp12m(lbp:ubp) ! relation for potential temperature profile applied at 2-m
   real(r8), intent(out) :: temp2(lbp:ubp)   ! relation for specific humidity profile
   real(r8), intent(out) :: temp22m(lbp:ubp) ! relation for specific humidity profile applied at 2-m
!inout:

!OTHER LOCAL VARIABLES:
   real(r8), parameter :: zetam = 1.574_r8   ! transition point of flux-gradient relation (wind profile)
   real(r8), parameter :: zetat = 0.465_r8   ! transition point of flux-gradient relation (temp. profile)
   integer :: f                              ! pft-filter index
   integer :: p                              ! pft index
   integer :: g                              ! gridcell index
   real(r8):: zldis(lbp:ubp)                 ! reference height "minus" zero displacement heght [m]
   real(r8):: zeta(lbp:ubp)                  ! dimensionless height used in Monin-Obukhov theory

!--- YangW: Following Variables For dustmodel -> Removed ---!
  !integer,  intent(in)  :: iter             ! iteration number
  !real(r8), intent(in)  :: ur(lbp:ubp)      ! wind speed at reference height [m/s]
  !real(r8), intent(out) :: u10(1)           ! 10-m wind (m/s) (for dust model)
  !real(r8), intent(out) :: fv(1)            ! friction velocity (m/s) (for dust model)
  !real(r8), intent(inout) :: fm(lbp:ubp)    ! needed for DGVM only to diagnose 10m wind
!------------------------------------------------------------------------------
! Adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

   do f = 1, fn
      p = filterp(f)
      g = pgridcell(p)

    !--- Wind profile---!
      zldis(p) = forc_hgt_u(g)-displa(p)
      zeta(p)  = zldis(p)/obu(p)
      if (zeta(p) < -zetam) then
          ustar(p) = vkc*um(p)/(log(-zetam*obu(p)/z0m(p))&
                     - StabilityFunc1(-zetam)            &
                     + StabilityFunc1(z0m(p)/obu(p))     &
                     + 1.14_r8*((-zeta(p))**0.333_r8-(zetam)**0.333_r8))
      else if (zeta(p) < 0._r8) then
          ustar(p) = vkc*um(p)/(log(zldis(p)/z0m(p))     &
                     - StabilityFunc1(zeta(p))           &
                     + StabilityFunc1(z0m(p)/obu(p)))
      else if (zeta(p) <=  1._r8) then
          ustar(p) = vkc*um(p)/(log(zldis(p)/z0m(p)) + 5._r8*zeta(p) -5._r8*z0m(p)/obu(p))
      else
          ustar(p) = vkc*um(p)/(log(obu(p)/z0m(p))+5._r8-5._r8*z0m(p)/obu(p) &
                     +(5._r8*log(zeta(p))+zeta(p)-1._r8))
      end if

    !--- Temperature profile ---!
      zldis(p) = forc_hgt_t(g)-displa(p)
      zeta(p)  = zldis(p)/obu(p)
      if (zeta(p) < -zetat) then
          temp1(p) = vkc/(log(-zetat*obu(p)/z0h(p))  &
                     - StabilityFunc2(-zetat)        &
                     + StabilityFunc2(z0h(p)/obu(p)) &
                     + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
      else if (zeta(p) < 0._r8) then
          temp1(p) = vkc/(log(zldis(p)/z0h(p))       &
                     - StabilityFunc2(zeta(p))       &
                     + StabilityFunc2(z0h(p)/obu(p)))
      else if (zeta(p) <=  1._r8) then
          temp1(p) = vkc/(log(zldis(p)/z0h(p)) + 5._r8*zeta(p) - 5._r8*z0h(p)/obu(p))
      else
          temp1(p) = vkc/(log(obu(p)/z0h(p)) + 5._r8 - 5._r8*z0h(p)/obu(p) &
                     + (5._r8*log(zeta(p))+zeta(p)-1._r8))
      end if

    !--- Humidity profile ---!
      if (forc_hgt_q(g) == forc_hgt_t(g) .and. z0q(p) == z0h(p)) then
          temp2(p) = temp1(p)
      else
          zldis(p) = forc_hgt_q(g)-displa(p)
          zeta(p)  = zldis(p)/obu(p)
          if (zeta(p) < -zetat) then
              temp2(p) = vkc/(log(-zetat*obu(p)/z0q(p))  &
                         - StabilityFunc2(-zetat)        &
                         + StabilityFunc2(z0q(p)/obu(p)) &
                         + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
          else if (zeta(p) < 0._r8) then
              temp2(p) = vkc/(log(zldis(p)/z0q(p))       &
                         - StabilityFunc2(zeta(p))       &
                         + StabilityFunc2(z0q(p)/obu(p)))
          else if (zeta(p) <=  1._r8) then
              temp2(p) = vkc/(log(zldis(p)/z0q(p)) + 5._r8*zeta(p)-5._r8*z0q(p)/obu(p))
          else
              temp2(p) = vkc/(log(obu(p)/z0q(p)) + 5._r8 - 5._r8*z0q(p)/obu(p) &
                         + (5._r8*log(zeta(p))+zeta(p)-1._r8))
         end if
      endif

    !--- Temperature profile applied at 2-m ---!
      zldis(p) = 2.0_r8 + z0h(p)
      zeta(p)  = zldis(p)/obu(p)
      if (zeta(p) < -zetat) then
          temp12m(p) = vkc/(log(-zetat*obu(p)/z0h(p))  &
                       - StabilityFunc2(-zetat)        &
                       + StabilityFunc2(z0h(p)/obu(p)) &
                       + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
      else if (zeta(p) < 0._r8) then
          temp12m(p) = vkc/(log(zldis(p)/z0h(p))       &
                       - StabilityFunc2(zeta(p))       &
                       + StabilityFunc2(z0h(p)/obu(p)))
      else if (zeta(p) <=  1._r8) then
          temp12m(p) = vkc/(log(zldis(p)/z0h(p)) + 5._r8*zeta(p) - 5._r8*z0h(p)/obu(p))
      else
          temp12m(p) = vkc/(log(obu(p)/z0h(p)) + 5._r8 - 5._r8*z0h(p)/obu(p) &
                       + (5._r8*log(zeta(p))+zeta(p)-1._r8))
      end if

    !--- Humidity profile applied at 2-m ---!
      if (z0q(p) == z0h(p)) then
          temp22m(p) = temp12m(p)
      else
          zldis(p) = 2.0_r8 + z0q(p)
          zeta(p)  = zldis(p)/obu(p)
          if (zeta(p) < -zetat) then
              temp22m(p) = vkc/(log(-zetat*obu(p)/z0q(p))                           &
                           - StabilityFunc2(-zetat) + StabilityFunc2(z0q(p)/obu(p)) &
                           + 0.8_r8*((zetat)**(-0.333_r8)-(-zeta(p))**(-0.333_r8)))
          else if (zeta(p) < 0._r8) then
              temp22m(p) = vkc/(log(zldis(p)/z0q(p))                                &
                           - StabilityFunc2(zeta(p))+StabilityFunc2(z0q(p)/obu(p)))
         else if (zeta(p) <=  1._r8) then
              temp22m(p) = vkc/(log(zldis(p)/z0q(p)) + 5._r8*zeta(p)-5._r8*z0q(p)/obu(p))
         else
              temp22m(p) = vkc/(log(obu(p)/z0q(p)) + 5._r8 - 5._r8*z0q(p)/obu(p)    &
                           + (5._r8*log(zeta(p))+zeta(p)-1._r8))
         end if
      end if
   end do
   end SUBROUTINE FrictionVelocity

! !IROUTINE: StabilityFunc
!
! !INTERFACE:
   real(r8) function StabilityFunc1(zeta)
!
! !DESCRIPTION:
! Stability function for rib < 0.
!
! !USES:
!      use shr_const_mod, only: SHR_CONST_PI
!Zack Subin, 7/8/08
!
! !ARGUMENTS:
      implicit none
      real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
!
! !CALLED FROM:
! SUBROUTINE FrictionVelocity in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
      real(r8) :: chik, chik2
!------------------------------------------------------------------------------

      chik2 = sqrt(1._r8-16._r8*zeta)
      chik = sqrt(chik2)
      StabilityFunc1 = 2._r8*log((1._r8+chik)*0.5_r8) &
!Changed to pie, Zack Subin, 7/9/08
           + log((1._r8+chik2)*0.5_r8)-2._r8*atan(chik)+pie*0.5_r8

    end function StabilityFunc1

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: StabilityFunc2
!
! !INTERFACE:
   real(r8) function StabilityFunc2(zeta)
!
! !DESCRIPTION:
! Stability function for rib < 0.
!
! !USES:
!Removed by Zack Subin, 7/9/08
!     use shr_const_mod, only: SHR_CONST_PI
!
! !ARGUMENTS:
     implicit none
     real(r8), intent(in) :: zeta  ! dimensionless height used in Monin-Obukhov theory
!
! !CALLED FROM:
! SUBROUTINE FrictionVelocity in this module
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
     real(r8) :: chik2
!------------------------------------------------------------------------------

     chik2 = sqrt(1._r8-16._r8*zeta)
     StabilityFunc2 = 2._r8*log((1._r8+chik2)*0.5_r8)

   end function StabilityFunc2

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: MoninObukIni
!
! !INTERFACE:
  SUBROUTINE MoninObukIni (ur, thv, dthv, zldis, z0m, um, obu)
!
! !DESCRIPTION:
! Initialization of the Monin-Obukhov length.
! The scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate,
! Vol. 11, 2628-2644.
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in)  :: ur    ! wind speed at reference height [m/s]
    real(r8), intent(in)  :: thv   ! virtual potential temperature (kelvin)
    real(r8), intent(in)  :: dthv  ! diff of vir. poten. temp. between ref. height and surface
    real(r8), intent(in)  :: zldis ! reference height "minus" zero displacement heght [m]
    real(r8), intent(in)  :: z0m   ! roughness length, momentum [m]
    real(r8), intent(out) :: um    ! wind speed including the stability effect [m/s]
    real(r8), intent(out) :: obu   ! monin-obukhov length (m)
!
! !CALLED FROM:
! SUBROUTINE BareGroundFluxes in module BareGroundFluxesMod.F90
! SUBROUTINE BiogeophysicsLake in module BiogeophysicsLakeMod.F90
! SUBROUTINE CanopyFluxes in module CanopyFluxesMod.F90
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
!
!EOP
!
! !LOCAL VARIABLES:
!
    real(r8) :: wc    ! convective velocity [m/s]
    real(r8) :: rib   ! bulk Richardson number
    real(r8) :: zeta  ! dimensionless height used in Monin-Obukhov theory
    real(r8) :: ustar ! friction velocity [m/s]
!-----------------------------------------------------------------------

    ! Initial values of u* and convective velocity

    ustar=0.06_r8
    wc=0.5_r8
    if (dthv >= 0._r8) then
       um=max(ur,0.1_r8)
    else
       um=sqrt(ur*ur+wc*wc)
    endif

    rib=grav*zldis*dthv/(thv*um*um)

    if (rib >= 0._r8) then      ! neutral or stable
       zeta = rib*log(zldis/z0m)/(1._r8-5._r8*min(rib,0.19_r8))
       zeta = min(2._r8,max(zeta,0.01_r8 ))
    else                     ! unstable
       zeta=rib*log(zldis/z0m)
       zeta = max(-100._r8,min(zeta,-0.01_r8 ))
    endif

    obu=zldis/zeta

 end SUBROUTINE MoninObukIni


SUBROUTINE lakeini(IVGTYP,         ISLTYP,          HT,              SNOW,           & !i
                   lake_min_elev,     restart,        lakedepth_default, lake_depth, &
                   lakedepth2d,    savedtke12d,     snowdp2d,        h2osno2d,       & !o
                   snl2d,          t_grnd2d,        t_lake3d,        lake_icefrac3d, &
                   z_lake3d,       dz_lake3d,       t_soisno3d,      h2osoi_ice3d,   &
                   h2osoi_liq3d,   h2osoi_vol3d,    z3d,             dz3d,           &
                   zi3d,           watsat3d,        csol3d,          tkmg3d,         &
                   iswater,        xice,            xice_threshold,  xland,   tsk,   &
                   lakemask,       lakeflag,                                         &
                   lake_depth_flag, use_lakedepth,                                   &
                   tkdry3d,        tksatu3d,        lake,        its, ite, jts, jte, &
                   ims,ime, jms,jme,                                                 &
                  !--- YangW: Updated tlake_init ---!
                   tlake_init_flag, tlake_init_value,                                &
                  !--- YangW: Updated InOut_put Variables For albedo process ---!
                   sabg_lakeout, cosz_lakeout, albedo_sols_yw, albedo_soll_yw, albedo_solsd_yw, albedo_solld_yw, &  ! YangW: Output 4-kind albedo  
                   z0mg_yw,z0hg_yw,z0qg_yw,km_yw,ke_yw,ked_yw,kw_yw)    
!==============================================================================
! This subroutine was first edited by Hongping Gu for coupling
! 07/20/2010
!==============================================================================
  USE module_wrf_error
  implicit none

  INTEGER, INTENT(IN) :: iswater
  REAL,    INTENT(IN) :: xice_threshold
  REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: XICE
  REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN)    :: TSK
  REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT) :: XLAND
  REAL, DIMENSION( ims:ime , jms:jme ) :: LAKEMASK
  INTEGER , INTENT(IN)                 :: lakeflag
  INTEGER , INTENT(INOUT) :: lake_depth_flag
  INTEGER , INTENT(IN)    :: use_lakedepth
  LOGICAL , INTENT(IN)    :: restart
  INTEGER,  INTENT(IN)    :: ims,ime, jms,jme
  INTEGER,  INTENT(IN)    :: its,ite, jts,jte
  INTEGER, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: IVGTYP,       &
                                                           ISLTYP
  REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)    :: HT
  REAL,    DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: SNOW
  real,    intent(in)                                   :: lakedepth_default,lake_min_elev
  real,    dimension(ims:ime,jms:jme ),intent(out)      :: lakedepth2d,   &
                                                           savedtke12d
  real,    dimension(ims:ime,jms:jme ),intent(out)      :: snowdp2d,      &
                                                           h2osno2d,      &
                                                           snl2d,         &
                                                           t_grnd2d
                                                                              
  real,    dimension( ims:ime,1:nlevlake, jms:jme ),INTENT(out)           :: t_lake3d,       &
                                                                             lake_icefrac3d, &
                                                                             z_lake3d,       &
                                                                             dz_lake3d
  real,    dimension( ims:ime,-nlevsnow+1:nlevsoil, jms:jme ),INTENT(out) :: t_soisno3d,     &
                                                                             h2osoi_ice3d,   &
                                                                             h2osoi_liq3d,   &
                                                                             h2osoi_vol3d,   &
                                                                             z3d,            &
                                                                             dz3d
  real,    dimension( ims:ime,1:nlevsoil, jms:jme ),INTENT(out)           :: watsat3d,       &
                                                                             csol3d,         &
                                                                             tkmg3d,         &
                                                                             tkdry3d,        &
                                                                             tksatu3d
  real,    dimension( ims:ime,-nlevsnow+0:nlevsoil, jms:jme ),INTENT(out) :: zi3d            
  LOGICAL, DIMENSION( ims:ime, jms:jme ),intent(out)                      :: lake
  REAL, OPTIONAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN)            :: lake_depth

  real,    dimension( ims:ime,1:nlevsoil, jms:jme )  :: bsw3d,    &
                                                        bsw23d,   &
                                                        psisat3d, &
                                                        vwcsat3d, &
                                                        watdry3d, &
                                                        watopt3d, &
                                                        hksat3d,  &
                                                        sucsat3d, &
                                                        clay3d,   &
                                                        sand3d   
  integer  :: n,i,j,k,ib,lev,bottom                             ! indices
  real(r8),dimension(ims:ime,jms:jme )    :: bd2d               ! bulk density of dry soil material [kg/m^3]
  real(r8),dimension(ims:ime,jms:jme )    :: tkm2d              ! mineral conductivity
  real(r8),dimension(ims:ime,jms:jme )    :: xksat2d            ! maximum hydraulic conductivity of soil [mm/s]
  real(r8),dimension(ims:ime,jms:jme )    :: depthratio2d       ! ratio of lake depth to standard deep lake depth 
  real(r8),dimension(ims:ime,jms:jme )    :: clay2d             ! temporary
  real(r8),dimension(ims:ime,jms:jme )    :: sand2d             ! temporary
  real(r8)                                :: scalez  = 0.025_r8 ! Soil layer thickness discretization (m)
  logical,parameter                       :: arbinit = .true.
  real,parameter                          :: defval  = -999.0
  integer                                 :: isl
  integer                                 :: numb_lak           ! for debug
  character*256 :: message
  !--- YangW: Updated Output Variables For albedo process ---!
    logical, intent(in) :: tlake_init_flag                      ! true for adopting the default tsk to initialize tlake
    real,    intent(in) :: tlake_init_value                     ! The initizlized value for tlake
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: sabg_lakeout       ! solar radiation absorbed by ground
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: cosz_lakeout       ! the final coszen output: 
                                                                          ! cosz_lakeout=coszen_yw  (when coszen_flag=.true.)
                                                                          ! cosz_lakeout=coszen_lsm (when coszen_flag=.true.)
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: albedo_sols_yw     ! albedo for vis direct  solar radiation
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: albedo_soll_yw     ! albedo for nir direct  solar radiation
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: albedo_solsd_yw    ! albedo for vis diffuse solar radiation
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: albedo_solld_yw    ! albedo for nir diffuse solar radiation
  !----------------------------------------------------------!
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: z0mg_yw            ! roughness length over ground, momentum      (m)
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: z0hg_yw            ! roughness length over ground, sensible heat (m)
    real, dimension( ims:ime,jms:jme ) ,intent(out) :: z0qg_yw            ! roughness length over ground, latent heat   (m)

    real, dimension( ims:ime,jms:jme )              , intent(out)  :: km_yw   ! molecular diffusion coefficient (m**2/s)
    real, dimension( ims:ime,1:nlevlake-1, jms:jme ), intent(out)  :: ke_yw   ! eddy diffusion coefficient (m**2/s) 
    real, dimension( ims:ime,1:nlevlake-1, jms:jme ), intent(out)  :: ked_yw  ! enhanced diffusion coefficient (m**2/s) 
    real, dimension( ims:ime,1:nlevlake, jms:jme ),   intent(out)  :: kw_yw   ! km+ke+ked or km+ke coefficient (m**2/s) 

  IF ( RESTART ) RETURN 

  !--- YangW: Updated Output Variables For albedo process & roughness length ---!
    DO j = jts,jte
       DO i = its,ite
          sabg_lakeout(i,j)    = defval
          cosz_lakeout(i,j)    = defval
          albedo_sols_yw(i,j)  = defval
          albedo_soll_yw(i,j)  = defval
          albedo_solsd_yw(i,j) = defval
          albedo_solld_yw(i,j) = defval

          z0mg_yw(i,j)         = defval
          z0hg_yw(i,j)         = defval
          z0qg_yw(i,j)         = defval
       ENDDO
    ENDDO
  !------------------------------------------------------------------------------!

  DO j = jts,jte
     DO i = its,ite
        snowdp2d(i,j) = snow(i,j)*0.005 ! SNOW in kg/m^2 and snowdp in m
        h2osno2d(i,j) = snow(i,j)       ! mm 
     ENDDO
  ENDDO

! initialize all the grid with default value 
  DO j = jts,jte
     DO i = its,ite
        lakedepth2d(i,j)            = defval
        snl2d(i,j)                  = defval
        do k = -nlevsnow+1,nlevsoil
           h2osoi_liq3d(i,k,j)      = defval
           h2osoi_ice3d(i,k,j)      = defval
           t_soisno3d(i,k,j)        = defval
           z3d(i,k,j)               = defval 
           dz3d(i,k,j)              = defval                           
        enddo
        do k = 1,nlevlake 
           t_lake3d(i,k,j)          = defval
           lake_icefrac3d(i,k,j)    = defval
           z_lake3d(i,k,j)          = defval
           dz_lake3d(i,k,j)         = defval

           kw_yw(i,k,j)             = defval ! YangW: Total diffusivity
        enddo
        !--- YangW : Initialize extra diffusivity ---!
          km_yw(i,j)               = defval
          do k = 1,nlevlake-1
             ke_yw(i,k,j)          = defval
             ked_yw(i,k,j)         = defval
          enddo
        !--------------------------------------------!
     ENDDO
  ENDDO

! judge whether the grid is lake grid
   numb_lak = 0
       do i=its,ite
         do j=jts,jte
         IF (lakeflag.eq.0) THEN    
            if(ht(i,j)>=lake_min_elev) then 
              if ( xice(i,j).gt.xice_threshold) then   !mchen
                   ivgtyp(i,j) = iswater
                   xland(i,j) = 2. 
                   lake_icefrac3d(i,1,j) = xice(i,j)
                   xice(i,j)=0.0
               endif
            endif

            if(ivgtyp(i,j)==iswater.and.ht(i,j)>=lake_min_elev) then 
                lake(i,j)  = .true.
                lakemask(i,j) = 1
                numb_lak   = numb_lak + 1
            else 
                lake(i,j)  = .false.
                lakemask(i,j) = 0
            end if
         ELSE
            if(lakemask(i,j).eq.1) then 
                lake(i,j)  = .true.
                numb_lak   = numb_lak + 1
                if ( xice(i,j).gt.xice_threshold) then   !mchen
                   ivgtyp(i,j) = iswater
                   xland(i,j) = 2. 
                   lake_icefrac3d(i,1,j) = xice(i,j)
                   xice(i,j)=0.0
                endif
             else  
                lake(i,j)  = .false.
             endif
         ENDIF   ! end if lakeflag=0
        end do
       end do
    write(message,*) "the total number of lake grid is :", numb_lak
   !CALL wrf_message(message)
! CALL LakeDebug(msg)
! initialize lake grid 

  DO j = jts,jte
     DO i = its,ite
        if( lake(i,j) ) then
          ! t_soisno3d(i,:,j) = tsk(i,j)  
          ! t_lake3d(i,:,j)   = tsk(i,j)
          ! t_grnd2d(i,j)     = tsk(i,j)
           z3d(i,:,j)            = 0.0
           dz3d(i,:,j)           = 0.0
           zi3d(i,:,j)           = 0.0
           h2osoi_liq3d(i,:,j)   = 0.0
           h2osoi_ice3d(i,:,j)   = 0.0
           lake_icefrac3d(i,:,j) = 0.0
           h2osoi_vol3d(i,:,j)   = 0.0
           snl2d(i,j)            = 0.0
          if ( use_lakedepth.eq.1 .and.lake_depth_flag.eq.0 ) then !mchen
             !call wrf_error_fatal ( 'STOP: You need lake-depth information. Rerun WPS or set use_lakedepth = 0')
          end if
          if ( use_lakedepth.eq.0 .and.lake_depth_flag.eq.1 ) then !mchen
              lake_depth_flag = 0 
          end if
          !--- YangW Updated: Make sure we use the actual lake depth for the sake of ndown run ---!
          !--- We come into the problem that the lakedepth is equal to 50 after the ndown run ---!
          lake_depth_flag = 1
          if ( lake_depth_flag.eq.1 ) then
             if (lake_depth(i,j) > 0.0) then 
                 lakedepth2d(i,j)   = lake_depth(i,j)
             else
                if ( lakedepth_default  > 0.0 ) then
                    lakedepth2d(i,j)   = lakedepth_default
                else 
                    lakedepth2d(i,j)   = spval
                endif
             endif
          else
             if ( lakedepth_default  > 0.0 ) then
                 lakedepth2d(i,j)   = lakedepth_default
             else 
                 lakedepth2d(i,j)   = spval
             endif
          endif
       endif
     ENDDO
  ENDDO 
  
  if( nlevlake.eq.10 ) then
     dzlak(1) = 0.1_r8
     dzlak(2) = 1._r8
     dzlak(3) = 2._r8
     dzlak(4) = 3._r8
     dzlak(5) = 4._r8
     dzlak(6) = 5._r8
     dzlak(7) = 7._r8
     dzlak(8) = 7._r8
     dzlak(9) = 10.45_r8
     dzlak(10)= 10.45_r8

     zlak(1) =  0.05_r8
     zlak(2) =  0.6_r8
     zlak(3) =  2.1_r8
     zlak(4) =  4.6_r8
     zlak(5) =  8.1_r8
     zlak(6) = 12.6_r8
     zlak(7) = 18.6_r8
     zlak(8) = 25.6_r8
     zlak(9) = 34.325_r8
     zlak(10)= 44.775_r8
  endif
  if( nlevlake.eq.25 ) then
     dzlak(1)  = 0.1_r8
     dzlak(2)  = 0.25_r8
     dzlak(3)  = 0.25_r8
     dzlak(4)  = 0.25_r8
     dzlak(5)  = 0.25_r8
     dzlak(6)  = 0.5_r8
     dzlak(7)  = 0.5_r8
     dzlak(8)  = 0.5_r8
     dzlak(9)  = 0.5_r8
     dzlak(10) = 0.75_r8
     dzlak(11) = 0.75_r8
     dzlak(12) = 0.75_r8
     dzlak(13) = 0.75_r8
     dzlak(14) = 2_r8
     dzlak(15) = 2_r8
     dzlak(16) = 2.5_r8
     dzlak(17) = 2.5_r8
     dzlak(18) = 3.5_r8
     dzlak(19) = 3.5_r8
     dzlak(20) = 3.5_r8
     dzlak(21) = 3.5_r8
     dzlak(22) = 5.225_r8
     dzlak(23) = 5.225_r8
     dzlak(24) = 5.225_r8
     dzlak(25) = 5.225_r8

     zlak(1) = dzlak(1)/2._r8
     do k = 2,nlevlake
        zlak(k) = zlak(k-1) + (dzlak(k-1)+dzlak(k))/2._r8
     end do
  endif

   ! "0" refers to soil surface and "nlevsoil" refers to the bottom of model soil

   do j = 1, nlevsoil
      zsoi(j) = scalez*(exp(0.5_r8*(j-0.5_r8))-1._r8)    !node depths
   enddo

   dzsoi(1) = 0.5_r8*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
   do j = 2,nlevsoil-1
      dzsoi(j)= 0.5_r8*(zsoi(j+1)-zsoi(j-1))
   enddo
   dzsoi(nlevsoil) = zsoi(nlevsoil)-zsoi(nlevsoil-1)

   zisoi(0) = 0._r8
   do j = 1, nlevsoil-1
      zisoi(j) = 0.5_r8*(zsoi(j)+zsoi(j+1))         !interface depths
   enddo
   zisoi(nlevsoil) = zsoi(nlevsoil) + 0.5_r8*dzsoi(nlevsoil)
!!!!!!!!!!!!!!!!!!begin to initialize lake variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO j = jts,jte
  DO i = its,ite
      
     if ( lake(i,j) ) then

                             ! Soil hydraulic and thermal properties
         isl = ISLTYP(i,j)   
         if (isl == 14 ) isl = isl + 1 
         do k = 1,nlevsoil
            sand3d(i,k,j)  = sand(isl)
            clay3d(i,k,j)  = clay(isl)
         enddo

         do k = 1,nlevsoil
            clay2d(i,j)     = clay3d(i,k,j)
            sand2d(i,j)     = sand3d(i,k,j)
            watsat3d(i,k,j) = 0.489_r8 - 0.00126_r8*sand2d(i,j)
            bd2d(i,j)       = (1._r8-watsat3d(i,k,j))*2.7e3_r8
            xksat2d(i,j)    = 0.0070556_r8 *( 10._r8**(-0.884_r8+0.0153_r8*sand2d(i,j)) )          ! mm/s
            tkm2d(i,j)      = (8.80_r8*sand2d(i,j)+2.92_r8*clay2d(i,j))/(sand2d(i,j)+clay2d(i,j))  ! W/(m K)
            bsw3d(i,k,j)    = 2.91_r8 + 0.159_r8*clay2d(i,j)
            bsw23d(i,k,j)   = -(3.10_r8 + 0.157_r8*clay2d(i,j) - 0.003_r8*sand2d(i,j))
            psisat3d(i,k,j) = -(exp((1.54_r8 - 0.0095_r8*sand2d(i,j) + 0.0063_r8*(100.0_r8-sand2d(i,j)  &
                              -clay2d(i,j)))*log(10.0_r8))*9.8e-5_r8)
            vwcsat3d(i,k,j) = (50.5_r8 - 0.142_r8*sand2d(i,j) - 0.037_r8*clay2d(i,j))/100.0_r8
            hksat3d(i,k,j)  = xksat2d(i,j)
            sucsat3d(i,k,j) = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand2d(i,j)) )
            tkmg3d(i,k,j)   = tkm2d(i,j) ** (1._r8- watsat3d(i,k,j))
            tksatu3d(i,k,j) = tkmg3d(i,k,j)*0.57_r8**watsat3d(i,k,j)
            tkdry3d(i,k,j)  = (0.135_r8*bd2d(i,j) + 64.7_r8) / (2.7e3_r8 - 0.947_r8*bd2d(i,j))
            csol3d(i,k,j)   = (2.128_r8*sand2d(i,j)+2.385_r8*clay2d(i,j)) / (sand2d(i,j)+clay2d(i,j))*1.e6_r8  ! J/(m3 K)
            watdry3d(i,k,j) = watsat3d(i,k,j) * (316230._r8/sucsat3d(i,k,j)) ** (-1._r8/bsw3d(i,k,j))
            watopt3d(i,k,j) = watsat3d(i,k,j) * (158490._r8/sucsat3d(i,k,j)) ** (-1._r8/bsw3d(i,k,j))
         end do
         if (lakedepth2d(i,j) == spval) then
            lakedepth2d(i,j)          = zlak(nlevlake) + 0.5_r8*dzlak(nlevlake)
            z_lake3d(i,1:nlevlake,j)  = zlak(1:nlevlake)
            dz_lake3d(i,1:nlevlake,j) = dzlak(1:nlevlake)
         else
            depthratio2d(i,j)         = lakedepth2d(i,j) / (zlak(nlevlake) + 0.5_r8*dzlak(nlevlake)) 
            z_lake3d(i,1,j)           = zlak(1)
            dz_lake3d(i,1,j)          = dzlak(1)
            dz_lake3d(i,2:nlevlake,j) = dzlak(2:nlevlake)*depthratio2d(i,j)
            z_lake3d(i,2:nlevlake,j)  = zlak(2:nlevlake)*depthratio2d(i,j) + dz_lake3d(i,1,j)*(1._r8 - depthratio2d(i,j))
         end if

        !--- YangW: Initilize the Temperature Profile ---!
        if (tlake_init_flag) then
            t_soisno3d(i,1,j)      = tsk(i,j)
            t_lake3d(i,1,j)        = tsk(i,j)
            t_grnd2d(i,j)          = 277.0
            do k = 2, nlevlake
               if(z_lake3d(i,k,j).le.depth_c) then 
                  t_soisno3d(i,k,j)=tsk(i,j)+(277.0-tsk(i,j))/depth_c*z_lake3d(i,k,j)
                  t_lake3d(i,k,j)  =tsk(i,j)+(277.0-tsk(i,j))/depth_c*z_lake3d(i,k,j)
               else
                  t_soisno3d(i,k,j)      = 277.0
                  t_lake3d(i,k,j)        = 277.0
               end if 
            enddo
        else
         !--- #ZZZ: 修改湖温初始化方案 --- !
            if (tlake_init_value < 273.15+3) then
               t_soisno3d(i,:,j)      = tlake_init_value
               t_lake3d(i,:,j)        = tlake_init_value
               t_grnd2d(i,j)          = t_lake3d(i,1,j)
            else
               t_soisno3d(i,1,j)      = tlake_init_value
               t_lake3d(i,1,j)        = tlake_init_value
               t_grnd2d(i,j)          = 277.0
               do k = 2, nlevlake
                  if(z_lake3d(i,k,j).le.depth_c) then 
                     t_soisno3d(i,k,j)=tlake_init_value+(277.0-tlake_init_value)/depth_c*z_lake3d(i,k,j)
                     t_lake3d(i,k,j)  =tlake_init_value+(277.0-tlake_init_value)/depth_c*z_lake3d(i,k,j)
                  else
                     t_soisno3d(i,k,j)      = 277.0
                     t_lake3d(i,k,j)        = 277.0
                  end if 
            end if
        end if

        z3d(i,1:nlevsoil,j)  = zsoi(1:nlevsoil)
        zi3d(i,0:nlevsoil,j) = zisoi(0:nlevsoil)
        dz3d(i,1:nlevsoil,j) = dzsoi(1:nlevsoil)
        savedtke12d(i,j)     = tkwat ! Initialize for first timestep.

        if (snowdp2d(i,j) < 0.01_r8) then
           snl2d(i,j) = 0
           dz3d(i,-nlevsnow+1:0,j) = 0._r8
           z3d (i,-nlevsnow+1:0,j) = 0._r8
           zi3d(i,-nlevsnow+0:0,j) = 0._r8
        else
           if ((snowdp2d(i,j) >= 0.01_r8) .and. (snowdp2d(i,j) <= 0.03_r8)) then
              snl2d(i,j) = -1
              dz3d(i,0,j)  = snowdp2d(i,j)
           else if ((snowdp2d(i,j) > 0.03_r8) .and. (snowdp2d(i,j) <= 0.04_r8)) then
              snl2d(i,j) = -2
              dz3d(i,-1,j) = snowdp2d(i,j)/2._r8
              dz3d(i, 0,j) = dz3d(i,-1,j)
           else if ((snowdp2d(i,j) > 0.04_r8) .and. (snowdp2d(i,j) <= 0.07_r8)) then
              snl2d(i,j) = -2
              dz3d(i,-1,j) = 0.02_r8
              dz3d(i, 0,j) = snowdp2d(i,j) - dz3d(i,-1,j)
           else if ((snowdp2d(i,j) > 0.07_r8) .and. (snowdp2d(i,j) <= 0.12_r8)) then
              snl2d(i,j) = -3
              dz3d(i,-2,j) = 0.02_r8
              dz3d(i,-1,j) = (snowdp2d(i,j) - 0.02_r8)/2._r8
              dz3d(i, 0,j) = dz3d(i,-1,j)
           else if ((snowdp2d(i,j) > 0.12_r8) .and. (snowdp2d(i,j) <= 0.18_r8)) then
              snl2d(i,j) = -3
              dz3d(i,-2,j) = 0.02_r8
              dz3d(i,-1,j) = 0.05_r8
              dz3d(i, 0,j) = snowdp2d(i,j) - dz3d(i,-2,j) - dz3d(i,-1,j)
           else if ((snowdp2d(i,j) > 0.18_r8) .and. (snowdp2d(i,j) <= 0.29_r8)) then
              snl2d(i,j) = -4
              dz3d(i,-3,j) = 0.02_r8
              dz3d(i,-2,j) = 0.05_r8
              dz3d(i,-1,j) = (snowdp2d(i,j) - dz3d(i,-3,j) - dz3d(i,-2,j))/2._r8
              dz3d(i, 0,j) = dz3d(i,-1,j)
           else if ((snowdp2d(i,j) > 0.29_r8) .and. (snowdp2d(i,j) <= 0.41_r8)) then
              snl2d(i,j) = -4
              dz3d(i,-3,j) = 0.02_r8
              dz3d(i,-2,j) = 0.05_r8
              dz3d(i,-1,j) = 0.11_r8
              dz3d(i, 0,j) = snowdp2d(i,j) - dz3d(i,-3,j) - dz3d(i,-2,j) - dz3d(i,-1,j)
           else if ((snowdp2d(i,j) > 0.41_r8) .and. (snowdp2d(i,j) <= 0.64_r8)) then
              snl2d(i,j) = -5
              dz3d(i,-4,j) = 0.02_r8
              dz3d(i,-3,j) = 0.05_r8
              dz3d(i,-2,j) = 0.11_r8
              dz3d(i,-1,j) = (snowdp2d(i,j) - dz3d(i,-4,j) - dz3d(i,-3,j) - dz3d(i,-2,j))/2._r8
              dz3d(i, 0,j) = dz3d(i,-1,j)
           else if (snowdp2d(i,j) > 0.64_r8) then
              snl2d(i,j) = -5
              dz3d(i,-4,j) = 0.02_r8
              dz3d(i,-3,j) = 0.05_r8
              dz3d(i,-2,j) = 0.11_r8
              dz3d(i,-1,j) = 0.23_r8
              dz3d(i, 0,j)=snowdp2d(i,j)-dz3d(i,-4,j)-dz3d(i,-3,j)-dz3d(i,-2,j)-dz3d(i,-1,j)
           endif
        end if
 
        do k = 0, snl2d(i,j)+1, -1
           z3d(i,k,j)    = zi3d(i,k,j) - 0.5_r8*dz3d(i,k,j)
           zi3d(i,k-1,j) = zi3d(i,k,j) - dz3d(i,k,j)
        end do

! 3:subroutine makearbinit
        if (snl2d(i,j) < 0) then
           do k = snl2d(i,j)+1, 0
                ! Be careful because there may be new snow layers with bad temperatures like 0 even if
                ! coming from init. con. file.
              if(arbinit .or. t_soisno3d(i,k,j) > 300 .or. t_soisno3d(i,k,j) < 200) t_soisno3d(i,k,j) = 250._r8
           enddo
        end if

        do k = 1, nlevsoil
           if(arbinit .or. t_soisno3d(i,k,j) > 1000 .or. t_soisno3d(i,k,j) < 0) t_soisno3d(i,k,j) = t_lake3d(i,nlevlake,j)
        end do

        do k = 1, nlevlake
           if(arbinit .or. lake_icefrac3d(i,k,j) > 1._r8 .or. lake_icefrac3d(i,k,j) < 0._r8) then
              if(t_lake3d(i,k,j) >= tfrz) then
                 lake_icefrac3d(i,k,j) = 0._r8
              else
                 lake_icefrac3d(i,k,j) = 1._r8
              end if
           end if
        end do
        
        do k = 1,nlevsoil
           if (arbinit .or. h2osoi_vol3d(i,k,j) > 10._r8 .or. h2osoi_vol3d(i,k,j) < 0._r8) h2osoi_vol3d(i,k,j) = 1.0_r8
           h2osoi_vol3d(i,k,j) = min(h2osoi_vol3d(i,k,j),watsat3d(i,k,j))

             ! soil layers
           if (t_soisno3d(i,k,j) <= tfrz) then
              h2osoi_ice3d(i,k,j)  = dz3d(i,k,j)*denice*h2osoi_vol3d(i,k,j)
              h2osoi_liq3d(i,k,j) = 0._r8
           else
              h2osoi_ice3d(i,k,j) = 0._r8
              h2osoi_liq3d(i,k,j) = dz3d(i,k,j)*denh2o*h2osoi_vol3d(i,k,j)
           endif
        enddo

        do k = -nlevsnow+1, 0
           if (k > snl2d(i,j)) then
              h2osoi_ice3d(i,k,j) = dz3d(i,k,j)*bdsno
              h2osoi_liq3d(i,k,j) = 0._r8
           end if
        end do

    end if   !lake(i,j)
  ENDDO
  ENDDO

  END SUBROUTINE lakeini

 SUBROUTINE roughness_lake (snl,t_grnd,t_lake,lake_icefrac,forc_psrf,&
                            cur,ustar,z0mg,z0hg,z0qg)
    IMPLICIT NONE
    integer,  INTENT(in) :: snl             ! number of snow layers
    real(r8), INTENT(in) :: t_grnd          ! ground temperature
    real(r8), INTENT(in) :: t_lake(1)       ! surface lake layer temperature       [K]
    real(r8), INTENT(in) :: lake_icefrac(1) ! surface lake layer ice mass fraction [0-1]
    real(r8), INTENT(in) :: forc_psrf       ! atmosphere pressure at the surface   [pa]
    real(r8), INTENT(in) :: cur             ! Charnock parameter (-)
    real(r8), INTENT(in) :: ustar           ! u* in similarity theory              [m/s]
    real(r8), parameter :: cus  = 0.1       ! empirical constant for roughness under smooth flow
    real(r8), parameter :: kva0 = 1.51e-5   ! kinematic viscosity of air (m^2/s) at 20C and 1.013e5 Pa
    real(r8), parameter :: prn  = 0.713     ! Prandtl # for air at neutral stability
    real(r8), parameter :: sch  = 0.66      ! Schmidt # for water in air at neutral stability

   !--- Output Vaariables ---!
    real(r8), INTENT(out) :: z0mg           ! roughness length over ground, momentum      [m]
    real(r8), INTENT(out) :: z0hg           ! roughness length over ground, sensible heat [m]
    real(r8), INTENT(out) :: z0qg           ! roughness length over ground, latent heat   [m]

    real(r8) kva    ! kinematic viscosity of air at ground temperature and forcing pressure
    real(r8) sqre0  ! root of roughness Reynolds number

      if ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0).or. lake_icefrac(1) < 0.001) then ! for unfrozen lake
          kva  = kva0 * (t_grnd/293.15)**1.5 * 1.013e5/forc_psrf    ! kinematic viscosity of air
          z0mg = max(cus*kva/max(ustar,1.e-4),cur*ustar*ustar/grav) ! momentum roughness length
          z0mg = max(z0mg, 1.0e-5)                                  ! This limit is redundant with current values.

          sqre0 = (max(z0mg*ustar/kva,0.1))**0.5                    ! square root of roughness Reynolds number	
          z0hg  = z0mg * exp( -vkc/prn*( 4.*sqre0 - 3.2) )          ! SH roughness length
          z0qg  = z0mg * exp( -vkc/sch*( 4.*sqre0 - 4.2) )          ! LH roughness length
          z0qg  = max(z0qg, 1.0e-5)                                 ! Minimum allowed roughness length for unfrozen lakes
          z0hg  = max(z0hg, 1.0e-5)                                 ! set low so it is only to avoid floating point exceptions
     
      else if (snl == 0) then ! frozen lake with ice without snow snow cover			           
          z0mg = 0.001        ! original set up	0.001	    		 
          z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.51e-5)**0.45)
          z0qg = z0hg
      else                    ! use roughness over snow
          z0mg = 0.0024             
          z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.51e-5)**0.45)
          z0qg = z0hg
      end if  
  end SUBROUTINE roughness_lake

SUBROUTINE orb_coszen(coszen,calday,dlon,dlat)
!-------------------------------------------------------------------------------
! FUNCTION to return the cosine of the solar zenith angle. Assumes 365.0 days/year.
! Compute earth/orbit parameters using formula suggested by
! Duane Thresher. Use formulas from Berger, Andre 1978: Long-Term Variations of Daily
! Insolation and Quaternary Climatic Changes. J. of the Atmo. Sci. 35:2362-2367.
!
! Original version:  Erik Kluzek, Oct/1997, Brian Kauffman, Jan/98
! CCSM2.0 standard
! yongjiu dai (07/23/2002)
!-------------------------------------------------------------------------------
 implicit none
 
 real(r8), intent(in)  :: calday        ! Julian cal day (1.xx to 365.xx)
 real(r8), intent(in)  :: dlat          ! Centered latitude (radians)
 real(r8), intent(in)  :: dlon          ! Centered longitude (radians)
 real(r8)              :: coszen        ! cosine of the solar zenith angle 
!--- Local variables ---!
 real(r8) declin                        ! Solar declination (radians)
 real(r8) eccf                          ! Earth-sun distance factor (ie. (1/r)**2)
 real(r8) lambm                         ! Lambda m, mean long of perihelion (rad)
 real(r8) lmm                           ! Intermediate argument involving lambm
 real(r8) lamb                          ! Lambda, the earths long of perihelion
 real(r8) invrho                        ! Inverse normalized sun/earth distance
 real(r8) sinl                          ! Sine of lmm
 real(r8) pi                            ! 3.14159265358979323846...
 real(r8), parameter :: &
           dayspy=365.0,              & ! days per year
           ve=80.5,                   & ! Calday of vernal equinox assumes Jan 1 = calday 1
           eccen=1.672393084E-2,      & ! Eccentricity
           obliqr=0.409214646,        & ! Earths obliquity in radians
           lambm0=-3.2625366E-2,      & ! Mean long of perihelion at the vernal equinox (radians)
           mvelpp=4.92251015            ! moving vernal equinox longitude of 
                                        ! perihelion plus pi (radians)
!-------------------------------------------------------------------------------
   pi     = 4.*atan(1.)
   lambm  = lambm0 + (calday - ve)*2.*pi/dayspy
   lmm    = lambm  - mvelpp
   sinl   = sin(lmm)
   lamb   = lambm + eccen*(2.*sinl + eccen*(1.25*sin(2.*lmm) &
            + eccen*((13.0/12.0)*sin(3.*lmm) - 0.25*sinl)))
   invrho = (1. + eccen*cos(lamb - mvelpp)) / (1. - eccen*eccen)
   declin = asin(sin(obliqr)*sin(lamb))
   eccf   = invrho*invrho

   coszen = sin(dlat)*sin(declin) &
            - cos(dlat)*cos(declin)*cos(calday*2.0*pi+dlon)

 end SUBROUTINE orb_coszen

SUBROUTINE ticktime (deltim, idate)
!=======================================================================
! Desciptions: Get the jday & msec at the next Integration Step         !!
! Original author: Yongjiu Dai, September 15, 1999
! Notes: Greenwich time
!=======================================================================
  implicit none

  real(r8), INTENT(in)   :: deltim   ! model time step (second)
  integer, INTENT(inout) :: idate(3) ! year and julian day, sec, respectively
  integer maxday                     ! days of one year (365 or 366)

  idate(3) = idate(3) + nint(deltim)
  if(idate(3)>86400)then
     idate(3) = idate(3) - 86400
     if((mod(idate(1),4)==0 .AND. mod(idate(1),100)/=0) .OR. &
                                  mod(idate(1),400)==0)then
         maxday = 366
     else
         maxday = 365
     endif
     idate(2) = idate(2) + 1
     if(idate(2)>maxday) then
        idate(1) = idate(1) + 1
        idate(2) = 1
     endif
  endif

 END SUBROUTINE ticktime

SUBROUTINE rd_forcing(srad,coszen,forc_sols,forc_soll,forc_solsd,forc_solld)
!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
! Notes: Greenwich time
! Desciptions: Get the Four solar radiation 
!            : forc_sols   ! atm vis direct beam solar rad onto srf
!              forc_soll   ! atm nir direct beam solar rad onto srf [W/m2]
!              forc_solsd  ! atm vis diffuse solar rad onto srf     [W/m2]
!              forc_solld  ! atm nir diffuse solar rad onto srf     [W/m2]
!=======================================================================
 IMPLICIT NONE  
 real(r8), intent(in) :: srad                  ! solar is the total incoming radiation srad=swdown
 real(r8)             :: coszen                ! cosine of the solar zenith angle
 real(r8)             :: hsolar, ratio_rvrf    ! Used when IRAD=true
 real(r8)             :: cloud, difrat, vnrat  ! Used when IRAD=false
 logical              :: IRAD=.false.
 real(r8)             :: forc_sols             ! atm visible direct beam solar rad onto srf             [W/m2]	 
 real(r8)             :: forc_soll             ! atm Near-infrared (NIR) direct beam solar rad onto srf [W/m2]
 real(r8)             :: forc_solsd            ! atm vis diffuse solar rad onto srf [W/m2]	   
 real(r8)             :: forc_solld            ! atm nir diffuse solar rad onto srf [W/m2] 
 
!--- YangW: Present IRAD=false ---!
 if(IRAD)then
    !---------------------------------------------------------------
    ! NOTE: codes from CLM4.5-CESM1.2.0
    ! relationship between incoming NIR or VIS radiation and ratio of
    ! direct to diffuse radiation calculated based on one year's worth of
    ! hourly CAM output from CAM version cam3_5_55
    !---------------------------------------------------------------
       hsolar = srad*0.5_R8
      ! NIR (dir, diff)
       ratio_rvrf = min(0.99_R8,max(0.29548_R8 + 0.00504_R8*hsolar  &
                    -1.4957e-05_R8*hsolar**2 + 1.4881e-08_R8*hsolar**3,0.01_R8))
       forc_soll  = ratio_rvrf*hsolar
       forc_solld = (1._R8 - ratio_rvrf)*hsolar
     ! VIS (dir, diff)
       ratio_rvrf = min(0.99_R8,max(0.17639_R8 + 0.00380_R8*hsolar  &
                    -9.0039e-06_R8*hsolar**2 + 8.1351e-09_R8*hsolar**3,0.01_R8))
       forc_sols  = ratio_rvrf*hsolar
       forc_solsd = (1._R8 - ratio_rvrf)*hsolar
 else
    !---------------------------------------------------------------
    ! as the downward solar is in full band, an empirical expression
    ! will be used to divide fractions of band and incident
    ! (visible, near-infrad, dirct, diffuse)
    ! Julian calday (1.xx to 365.xx)
    !---------------------------------------------------------------
       cloud = (1160.*coszen-srad)/(963.*max(coszen,0.0001))
       cloud = max(cloud,0.)
       cloud = min(cloud,1.)
       cloud = max(0.58,cloud)  

       difrat = 0.0604/(coszen-0.0223)+0.0683
       if(difrat.lt.0.) difrat = 0.
       if(difrat.gt.1.) difrat = 1.

       difrat = difrat+(1.0-difrat)*cloud
       vnrat  = (580.-cloud*464.)/((580.-cloud*499.)+(580.-cloud*464.))

       forc_sols   = srad*(1.0-difrat)*vnrat         ! atm visible direct beam solar rad onto srf [W/m2] 
       forc_soll   = srad*(1.0-difrat)*(1.0-vnrat)   ! atm Near-infrared (NIR) direct beam solar rad onto srf [W/m2]
       forc_solsd  = srad*difrat*vnrat               ! atm vis diffuse solar rad onto srf [W/m2]   
       forc_solld  = srad*difrat*(1.0-vnrat)         ! atm nir diffuse solar rad onto srf [W/m2]
endif 
END SUBROUTINE rd_forcing

SUBROUTINE snowage_bats ( deltim,tg,scv,scvold,sag )
!=======================================================================
! Original version: Robert Dickinson
! Update snow cover and snow age, based on BATS code
!=======================================================================
  implicit none
!-------------------------- Dummy Argument -----------------------------
  real(r8), INTENT(in)    :: deltim ! seconds in a time step            [second] deltim=dtime=dt
  real(r8), INTENT(in)    :: tg     ! temperature surface               [K]      tg=t_grnd
  real(r8), INTENT(in)    :: scv    ! snow cover, water equivalent      [mm]     scv=h2osno
  real(r8), INTENT(in)    :: scvold ! snow cover for previous time step [mm]     scvold=h2osno_old
  real(r8), INTENT(inout) :: sag    ! non dimensional snow age [-]               sag=sag_yw    
!-------------------------- Local variables ----------------------------
  real(r8) :: age1   ! snow aging factor due to crystal growth [-]
  real(r8) :: age2   ! snow aging factor due to surface growth [-]
  real(r8) :: age3   ! snow aging factor due to accum of other particles [-]
  real(r8) :: arg    ! temporary variable used in snow age calculation [-]
  real(r8) :: arg2   ! temporary variable used in snow age calculation [-]
  real(r8) :: dela   ! temporary variable used in snow age calculation [-]
  real(r8) :: dels   ! temporary variable used in snow age calculation [-]
  real(r8) :: sge    ! temporary variable used in snow age calculation [-]
!-----------------------------------------------------------------------
      if(scv <= 0.) then
         sag = 0.
! Over antarctica
      else if (scv > 800.) then
         sag = 0.
! Away from antarctica
      else
         age3  = 0.3
         arg   = 5.e3*(1./tfrz-1./tg)
         arg2  = min(0.,10.*arg)
         age2  = exp(arg2)
         age1  = exp(arg)
         dela  = 1.e-6*deltim*(age1+age2+age3)
         dels  = 0.1*max(0.0,scv-scvold)
         sge   = (sag+dela)*(1.0-dels)
         sag   = max(0.0,sge)
      end if
end SUBROUTINE snowage_bats

SUBROUTINE alblake_yw(snowfraction_yw,snowage_yw,coszen,h2osno_yw,t_grnd_yw,snl_yw,alb_yw,lakedebug_yw)
!=======================================================================
! Original Concepts: Subin
! Brief Decription : Calculate albedo on different landtypes and solar radiation (visible vs infrared) 
!                    Snow albedo calculation is based on snowage when h2osno > 0 
!                    Return the averged albedo  based on snowfraction
!=======================================================================
IMPLICIT NONE
!--- Input Variables ---!
real(r8), intent(in)  :: snowfraction_yw ! non-dimensional snow age
real(r8), intent(in)  :: snowage_yw      ! snow fraction defined by YW
real(r8), intent(in)  :: coszen          ! cosine of the solar zenith angle
real(r8), intent(in)  :: h2osno_yw       ! snow water [mm]
real(r8), intent(in)  :: t_grnd_yw       ! ground temperature
integer,  intent(in)  :: snl_yw          ! number of snow layers
real(r8)              :: alb_yw(2,2)     ! averaged albedo according to different landtype and wave lendths of solar radiation [-]
logical,  intent(in)  :: lakedebug_yw

!--- Local Variables For Calculating snow albedo---!
real(r8) age,       & ! factor to reduce visible snow alb due to snow age [-]
         cff,       & ! snow alb correction factor for zenith angle > 60 [-]
         conn,      & ! constant (=0.5) for visible snow alb calculation [-]
         cons,      & ! constant (=0.2) for nir snow albedo calculation [-]
         czf,       & ! solar zenith correction for new snow albedo [-]
         dfalbl,    & ! snow albedo for diffuse nir radiation [-]
         dfalbs,    & ! snow albedo for diffuse visible solar radiation [-]
         dralbl,    & ! snow albedo for visible radiation [-]
         dralbs,    & ! snow albedo for near infrared radiation [-]
         sl,        & ! factor that helps control alb zenith dependence [-]
         snal0,     & ! alb for visible,incident on new snow (zen ang<60) [-]
         snal1        ! alb for NIR,    incident on new snow (zen angle<60) [-] 
real(r8) :: albsno_yw(2,2)   ! snow albedo according to different landtype and wave lendths of solar radiation [-]

real(r8) :: axx
! ----------------------------------------------------------------------
! 1. Initial set
! ----------------------------------------------------------------------
  ! short and long wave albedo for new snow
      snal0 = 0.85     ! shortwave
      snal1 = 0.65     ! long wave 
! ----------------------------------------------------------------------
! set default snow and final averaged albedo
      alb_yw(:,:)    = 0.  ! averaged albedo
      albsno_yw(:,:) = 0.  ! snow     albedo
     
! ----------------------------------------------------------------------
! 2. albedo for snow cover.
!    snow albedo depends on snow-age, zenith angle, and thickness
!    of snow age gives reduction of visible radiation
! ----------------------------------------------------------------------
    if (h2osno_yw>0.) then
        cons = 0.2
        conn = 0.5
        sl   = 2.0                          ! sl helps control albedo zenith dependence
       !correction for snow age
        age    = 1.-1./(1.+snowage_yw)      ! correction for snow age
        dfalbs = snal0*(1.-cons*age)        ! snow albedo for diffuse visible radiation
       !czf corrects albedo of new snow for solar zenith
        cff    = ((1.+1./sl)/(1.+coszen*2.*sl )- 1./sl)
        cff    = max(cff,0.)
        czf    = 0.4*cff*(1.-dfalbs)
        dralbs = dfalbs+czf          ! snow albedo for direct visible radiation
        dfalbl = snal1*(1.-conn*age) ! snow albedo for diffuse infrared radiation
        czf    = 0.4*cff*(1.-dfalbl)
        dralbl = dfalbl+czf          ! snow albedo for direct infrared radiation

        albsno_yw(1,1) = dralbs
        albsno_yw(2,1) = dralbl
        albsno_yw(1,2) = dfalbs
        albsno_yw(2,2) = dfalbl
    endif

!-----------------------------------------------------------------------
!3. albedo for unfrozen, frozen lake without snowcover
!-----------------------------------------------------------------------
  !--- albedo For unfrozen lake (also wetland,swamps,rice paddies) ---!
    alb_yw(:,1) = 0.05/(coszen+0.15)      ! for vis & nir direct radiation
    alb_yw(:,2) = 0.1                     ! for vis & nir diffuse radiation
 
  !--- albedo For frozen lake without snow layers ---!
     if (t_grnd_yw<=tfrz.and.snl_yw.eq.0) then  
        axx         = EXP(95*(t_grnd_yw-tfrz)/tfrz)
        alb_yw(1,1) = max(alb_yw(1,1),0.6*(1-axx)+0.1*axx)   ! for vis direct 
        alb_yw(1,2) = max(alb_yw(1,2),0.6*(1-axx)+0.1*axx)   ! for vis diffuse
        alb_yw(2,1) = max(alb_yw(2,1),0.4*(1-axx)+0.1*axx)   ! for nir direct 
        alb_yw(2,2) = max(alb_yw(2,2),0.4*(1-axx)+0.1*axx)   ! for nir diffuse
     endif 

    if (lakedebug_yw) then
        print *,"******************"
        print *,"t_grnd_yw=",t_grnd_yw
        print *,"snl_yw=",snl_yw
    endif
  !--- albedo For frozen lake with snow layers ---!
  !--- YangW: Shouldn't the following 4 lines be annotated ?? ---!
     if (t_grnd_yw<=tfrz.and.snl_yw<0) then   
         alb_yw(1,:) = 0.6                  ! for vis direct and diffuse radiation
         alb_yw(2,:) = 0.4                  ! for nir direct and diffuse radiation
     endif
 
!-----------------------------------------------------------------------
! 4 correction due to snow cover
!-----------------------------------------------------------------------
   alb_yw(:,:) = (1.-snowfraction_yw)*alb_yw(:,:) + snowfraction_yw*albsno_yw(:,:)          
   if (lakedebug_yw) then
       print *,"snowfraction_yw=",snowfraction_yw
       print *, "albedo_sols=",alb_yw(1,1),&
                "; albedo_soll=",alb_yw(2,1),&
                "; albedo_solsd=",alb_yw(1,2),&
                "; albedo_solld=",alb_yw(2,2)
   endif

end SUBROUTINE  alblake_yw

END MODULE module_sf_lake

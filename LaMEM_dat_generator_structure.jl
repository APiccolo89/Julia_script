"""
LaMEM generation setup 
Structure -> Grid 
Structure -> Timestepping 
Structure -> Free Surface
Structure -> Controls
Structure -> Visualisation
Structure -> Phase 
Structure -> Marker
Structure -> Solver Options (basic)
File structure 
File function (print .dat file)
flag 2 not printing none, -max(Float/Int)



"""

using Parameters


@with_kw struct Scaling
    units::String              = "geo"
	unit_temperature::Float64 = 1.0
	unit_length::Float64      = 1e3
	unit_viscosity ::Float64  = 1e18
	unit_stress::Float64      = 1e6
	unit_density::Float64     = 1e3
end

@with_kw struct Output_Options
    """
    -INT MAX/-FLOAT MAX are not active by default
    """
	time_end:: Float64  = 1.0     # simulation end time
	dt::Float64        = 0.05     # time step
	dt_min ::Float64   = 0.01     # minimum time step (declare divergence if lower value is attempted)
	dt_max::Float64    = 0.2      # maximum time step
	dt_out ::Float64   = 0.2      # output step (output at least at fixed time intervals)
	inc_dt::Float64    = 0.1      # time step increment per time step (fraction of unit)
	CFL ::Float64      = 0.5      # CFL (Courant-Friedrichs-Lewy) criterion
	CFLMAX ::Float64   = 0.8      # CFL criterion for elasticity
	nstep_max::Int32 = 50         # maximum allowed number of steps (lower bound: time_end/dt_max)
	nstep_out::Int32 = 1           # save output every n steps
	nstep_ini::Int32 = -max_value(Int32)     # save output for n initial steps
	nstep_rdb::Int32 = 10         # save restart database every n steps
	time_tol::Float64  = 1e-8      # relative tolerance for time comparisons
    
end 

@with_kw struct Grid_Options
    nel_x::Array{Int32} = 100
    seg_x:: Int32 = 1
    bias_x :: Array{Float64} = 1.0 
    coord_x:: Array{Float64} = [0.0 1.0] 
    nel_y::Array{Int32} = 100
    seg_y:: Int64 = 1
    bias_y :: Array{Float64} = 1.0 
    coord_y:: Array{Float64} = [0.0 1.0]
    nel_z:: Array{Int32} = 100
    seg_z:: Int64 = 1
    bias_z :: Array{Float64} = 1.0 
    coord_z:: Array{Float64} = [0.0 1.0]
end

@with_kw struct Free_Surface
    surf_use::Int32           = 0              # free surface activation flag
	surf_corr_phase::Int32    = 1              # air phase ratio correction flag (due to surface position)
	surf_level::Float64         = 0.5          # initial level
	surf_air_phase::Int32    = 0               # phase ID of sticky air layer
	surf_max_angle::Float64   = 45.0           # maximum angle with horizon (smoothed if larger)
	surf_topo_file::String     = "none"        # initial topography file (redundant)
end 


@with_kw struct Eros   
    erosion_model::Int32     = 2                # erosion model [0-none (default), 1-infinitely fast, 2-prescribed rate with given level]
    er_num_phases ::Int32     = 3                # number of erosion phases
    er_time_delims::Array{Float64}    = [0.5   2.5]      # erosion time delimiters (one less than number)
    er_rates ::Array{Float64}          = [0.2 0.1 0.2 ]     # constant erosion rates in different time periods
    er_levels ::Array{Float64}         = [1   2   1 ]       # levels above which we apply constant erosion rates in different time periods
end 

@with_kw struct Sedimentation 
    sediment_model::Int32      = 1                # sedimentation model [0-none (dafault), 1-prescribed rate with given level, 2-cont. margin]
    sed_num_layers::Int32     = 3                # number of sediment layers
    sed_time_delims::Array{Float64}   = [0.5   2.5]        # sediment layers time delimiters (one less than number)
    sed_rates::Array{Float64}         = [0.3 0.2 0.3]     # sediment rates in different time periods
    sed_levels::Array{Float64}         = [-5  -2  -2 ]      # levels below which we apply constant sediment rates in different time periods 
    sed_phases         = [1   2   3]       # sediment layers phase numbers in different time periods
end


@with_kw struct Cont_Margin
    marginO::Array{Float64}           = [0.0 0.0]          # lateral coordinates of continental margin - origin
    marginE::Array{Float64}           = [10.0 10.0        # lateral coordinates of continental margin - 2nd point
    hUp::Float64              = 1.5              # up dip thickness of sediment cover (onshore)
    hDown::Float64             = 0.1              # down dip thickness of sediment cover (off shore)
    dTrans::Float64            = 1.0              # half of transition zone
end

@with_kw struct Control_Values
    gravity:: Array{Float64}         = [0.0 0.0 -10.0]  # gravity vector
	FSSA:: Float64                   = 1.0            # free surface stabilization parameter [0 - 1]
	shear_heat_eff:: Float64         = 0.0            # shear heating efficiency parameter   [0 - 1]
	Adiabatic_Heat:: Float64         = 0.0            # Adiabatic Heating activaction flag and efficiency. [0.0 - 1.0] (e.g. 0.5 means that only 50% of the potential adiabatic heating affects the energy equation)   
	act_temp_diff:: Int32            = 0              # temperature diffusion activation flag
	act_therm_exp:: Int32            = 0              # thermal expansion activation flag
	act_steady_temp:: Int32          = 0              # steady-state temperature initial guess activation flag
	steady_temp_t::Float64           = 0.0            # time for (quasi-)steady-state temperature initial guess
	nstep_steady:: Int32             = 0              # number of steps for (quasi-)steady-state temperature initial guess (default = 1)
	init_lith_pres:: Int32           = 0              # initial pressure with lithostatic pressure (stabilizes compressible setups in the first steps)
	init_guess:: Int32               = 0             # initial guess flag
	p_litho_visc:: Int32             = 0              # use lithostatic pressure for creep laws
	p_litho_plast:: Int32            = 0              # use lithostatic pressure for plasticity
	p_lim_plast:: Int32              = 1              # limit pressure at first iteration for plasticity
	p_shift:: Int32 		         = 0              # constant [MPa] added to the total pressure field, before evaluating plasticity (e.g., when the domain is located @ some depth within the crust)  	
	act_p_shift :: Int32             = 0              # pressure shift activation flag (enforce zero pressure on average in the top cell layer); note: this overwrites p_shift above!
	eta_min::Float64                 = 1e18           # viscosity lower bound [Pas]
	eta_max::Float64                 = 1e25           # viscosity upper limit [Pas]
	eta_ref ::Float64                = 1e20           # reference viscosity (initial guess) [Pas]
	T_ref::Float64                   = 20.0           # reference temperature [C]
	RUGC ::Float64                   = 8.31           # universal gas constant (required only for non-dimensional setups)
	min_cohes ::Float64              = 2e7            # cohesion lower bound  [Pa]
	min_fric  ::Float64              = 5.0            # friction lower bound  [degree]
	tau_ult  ::Float64               = 1e9            # ultimate yield stress [Pa]
	rho_fluid ::Float64              = -max(Float64)            # fluid density for depth-dependent density model
	gw_level_type::String            = "none"            # ground water level type for pore pressure computation (see below)
	gw_level ::Float64               = -max(Float64)           # ground water level at the free surface (if defined)
	biot      ::Float64              = -max(Float64)            # Biot pressure parameter
	get_permea::Int32                = 0              # effective permeability computation activation flag
	rescal::Int32                    = 0             # stencil rescaling flag (for internal constraints, for example while computing permeability)
	mfmax::Float64                   = -max(Float64)            # maximum melt fraction affecting viscosity reduction
	lmaxit::Int32                    = 25             # maximum number of local rheology iterations 
	lrtol::Float64                   = 1e-6           # local rheology iterations relative tolerance
	act_dike::Int32                  = 0              # dike activation flag (additonal term in divergence)
	useTk::Int32                     = 0              # switch to use T-dependent conductivity, 0: not active
	Compute_velocity_gradient::Int32 = 0    # compute the velocity gradient tensor 1: active, 0: not active. If active, it automatically activates the output in the .pvd file
end

#####################################################
# Boundary Conditions                               #
#####################################################
@ks_args struct Basic_BoundaryC

    open_top_bound::Int32 = 1
    open_bot_bound::Int32 = 0
    permeable_phase_inflow::Int32 = 1 # Phase of the inflow material from the bottom (The temperature of the inflow phase it is the same of the bottom boundary)
    noslip::Array{Int32} = [0 0 0 0 0 0]
    fix_phase::Int32 = -max(Int32)
    fix_cell::Int32 = -max(Int32)
    fix_cell_file::String = "none" 
    temp_top::Float64               =   0.0
    temp_bot::Float64                =   1300.0                              # if only one value is given, it is assumed to be constant with time
    temp_bot_num_periods::Float64   =   2                                   # How many periods with different temp_bot do we have? 
    temp_bot_time_delim::Float64    =   1.2                                 # At which time do we switch from one to the next period?
    init_temp::Int32 = 1;
    pres_top::Float64 = 0.0
	pres_bot::Float64 = 10.0
	init_pres::Int32  = 1
end

@kw_args struct Markers 
    msetup::String  = "geom"            # setup type
	nmark_x::Int32         = 2                 # markers per cell in x-direction
	nmark_y::Int32        = 2                 # ...                 y-direction
	nmark_z::Int32        = 2                 # ...                 z-direction
	rand_noise::Int32     = 1                 # random noise flag
	rand_noiseGP::Int32    = 1                 # random noise flag, subsequently applied to geometric primitives
	bg_phase::Int32        = 1                 # background phase ID
	save_mark::Int32      = 0                # save marker to disk flag
	mark_load_file::String  = "./markers/mdb"     # marker input file (extension is .xxxxxxxx.dat)
	mark_save_file::String  = "none"     # marker output file (extension is .xxxxxxxx.dat)
	poly_file::String      = "none"  # polygon geometry file    (redundant)
	temp_file::String       = "none"  # initial temperature file (redundant)
	advect::String        = "basic"            # advection scheme
	interp::String         = "stag"              # velocity interpolation scheme
	stagp_a::Float64         = 0.7               # STAG_P velocity interpolation parameter
	mark_ctrl::String       = "none"              # marker control type
	nmark_lim ::Array{Int32}     = [10 100]            # min/max number per cell (marker control)
	nmark_avd ::Array{Int32}      = [3 3 3]             # x-y-z AVD refinement factors (avd marker control)
	nmark_sub ::Array{Int32}      = 1                 # max number of same phase markers per subcell (subgrid marker control)
end

####################################################################
# Geometric primitives
###############################################################
@kw_args struct Sphere 
    phase::Int32     = 1
	radius::Float64      = 1.5
	center::Float64      = [1.0 2.0 3.0]
	Temperature::String = "constant" # optional: Temperature of the sphere. possibilities: [constant]
	cstTemp::Float64     = 1000     # required in case of [constant]: temperature value [in Celcius in case of GEO units]
end


@kw_args struct Elispes 
    phase::Int32     = 1
	radius::Float64      = 1.5
	center::Array{Float64}     = [1.0 2.0 3.0]
    base::Array{Float64}        = [1.0 2.0 3.0]
	cap::Array{Float64}         = [3.0 5.0 7.0]
	Temperature::String = "constant" # optional: Temperature of the sphere. possibilities: [constant]
	cstTemp::Float64     = 1000     # required in case of [constant]: temperature value [in Celcius in case of GEO units]
end

@kw_args struct Cylinder 
    phase::Int32     = 1
	radius::Float64      = 1.5
	center::Array{Float64}      = [1.0 2.0 3.0]
	Temperature::String = "constant" # optional: Temperature of the sphere. possibilities: [constant]
	cstTemp::Float64     = 1000     # required in case of [constant]: temperature value [in Celcius in case of GEO units]
end

@kw_args struct Box
    phase::Int32      = 1
    bounds:: Array{Float64}      = [1.0 2.0 1.0 2.0 1.0 2.0] # box bound coordinates: left, right, front, back, bottom, top
    Temperature::String = "linear"  # optional: Temperature structure. possibilities: [constant, linear, halfspace]
    cstTemp::Float64    = 1000    # required in case of [constant]: temperature value [in Celcius in case of GEO units]
    topTemp::Float64    = 0       # required in case of [linear,halfspace]: temperature @ top [in Celcius in case of GEO units]
    botTemp::Float64    = 1300    # required in case of [linear,halfspace]: temperature @ bottom [in Celcius in case of GEO units]
    thermalAge::Float64  = 70      # required in case of [halfspace]: thermal age of lithosphere [in Myrs if GEO units are used]
end

@kw_args struct Ridge 

    phase::Int32       = 1
    bounds::Array{Float64}       = [1.0 2.0 1.0 2.0 1.0 2.0 ]  # box bound coordinates: left, right, front, back, bottom, top [top is seafloor]
    ridgeseg_x::Array{Float64}   = [1.5 1.5]                   # coordinate order: left, right [can be different for oblique ridge]
    ridgeseg_y::Array{Float64}   = [1.0 2.0 ]                  # coordinate order: front, back [needs to be the same as the front and back of bounds]
    topTemp::Float64    = 0                         # required: temperature @ top [in Celcius in case of GEO units]
    botTemp::Float64    = 1300                      # required: temperature @ bottom [in Celcius in case of GEO units]
    Temperature::String = "halfspace_age"             # initial temperature structure [ridge must be set to halfspace_age --> setTemp=4]
    age0::Float64        = 0.001                     # minimum age of seafloor at ridge [in Myr in case of GEO units]
    maxAge::Float64      = 60                        # [optional] parameter that indicates the maximum thermal age of a plate 
    v_spread::Float64   = 1                         # [optional] parameter that indicates the spreading velocity of the plate; if not defined it uses bvel_velin specified above

end


@kw_args struct Layer 
	phase::Int32       = 1
	top::Float64        = 5.0
	bottom::Float64     = 3.0
    cosine::Int32     = 0         # optional: add a cosine perturbation on top of the interface (if 1)
	wavelength::Float64  = 1         # required if cosine: wavelength in x-direction
	amplitude::Float64   = 0.1       # required if cosine: amplitude of perturbation         
	Temperature::String = "halfspace" # optional: Temperature structure. possibilities: [constant, linear, halfspace]
	cstTemp::Float64    = 1000      # required in case of [constant]: temperature value [in Celcius in case of GEO units]
	topTemp::Float64    = 0         # required in case of [linear,halfspace]: temperature @ top [in Celcius in case of GEO units]
	botTemp::Float64     = 1300      # required in case of [linear,halfspace]: temperature @ bottom [in Celcius in case of GEO units]
	thermalAge::Float64  = 70        # required in case of [halfspace]: thermal age of lithosphere [in Myrs if GEO units are used]
end













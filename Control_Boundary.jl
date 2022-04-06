
using Parameters; 
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
	rho_fluid ::Float64              = -typemax(Float64)            # fluid density for depth-dependent density model
	gw_level_type::String            = "none"            # ground water level type for pore pressure computation (see below)
	gw_level ::Float64               = -typemax(Float64)           # ground water level at the free surface (if defined)
	biot      ::Float64              = -typemax(Float64)            # Biot pressure parameter
	get_permea::Int32                = 0              # effective permeability computation activation flag
	rescal::Int32                    = 0             # stencil rescaling flag (for internal constraints, for example while computing permeability)
	mfmax::Float64                   = -typemax(Float64)            # maximum melt fraction affecting viscosity reduction
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
    fix_phase::Int32 = -typemax(Int32)
    fix_cell::Int32 = -typemax(Int32)
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

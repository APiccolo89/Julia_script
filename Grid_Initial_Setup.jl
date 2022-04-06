
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
	nstep_ini::Int32 = -typemax(Int32)     # save output for n initial steps
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

#######################################################
# Markers 
######################################################

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

@with_kw struct Passive_Tracers
    Passive_Tracer::Int32             = 1                               # Activate passive tracers?
    PassiveTracer_Box::Array{Float64}           =  [-600 600 -1 1 -300 -50]         # Dimensions of Box in which we disttribute passive tracers  
    PassiveTracer_Resolution::Array{Int32}  =  [100 1 100]                      # The number of passive tracers in every direction
    PassiveTracer_ActiveType::String    =  "Always"                  		  # Under which condition are they activated? []  
    PassiveTracer_ActiveValue::Float64   =  0.1                            # When this value is exceeded the tracers are being activated 
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


############################################################################################
#  Output Structures
############################################################################################

@kw_args struct Grid_Output
    out_file_name::String       = "output" # output file name
    out_pvd::Int32             = 0      # activate writing .pvd file
    out_phase::Int32           = 0
    out_density::Int32        = 0
    out_visc_total::Int32     = 0
    out_visc_creep::Int32      = 0
    out_velocity::Int32        = 0
    out_pressure::Int32        = 0
    out_tot_press::Int32       = 0
    out_eff_press::Int32       = 0
    out_over_press::Int32      = 0
    out_litho_press::Int32     = 0
    out_pore_pres::Int32s      = 0
    out_temperature::Int32     = 0
    out_dev_stress::Int32      = 0
    out_j2_dev_stress::Int32   = 0
    out_strain_rate::Int32     = 0
    out_j2_strain_rate::Int32  = 0
    out_shmax ::Int32          = 0
    out_ehmax ::Int32          = 0
    out_yield ::Int32          = 0
    out_rel_dif_rate::Int32    = 0     # relative proportion of diffusion creep strainrate
    out_rel_dis_rate::Int32    = 0     # relative proportion of dislocation creep strainrate
    out_rel_prl_rate::Int32    = 0     # relative proportion of peierls creep strainrate
    out_rel_pl_rate::Int32     = 0     # relative proportion of plastic strainrate
    out_plast_strain::Int32   = 0     # accumulated plastic strain
    out_plast_dissip::Int32    = 0     # plastic dissipation
    out_tot_displ::Int32       = 0
    out_moment_res::Int32      = 0
    out_cont_res::Int32        = 0
    out_energ_res::Int32       = 0
    out_melt_fraction::Int32   = 0
    out_fluid_density::Int32   = 0
    out_conductivity::Int32    = 0
    out_vel_gr_tensor::Int32   = 0
    out_avd::Int32     = 1 # activate AVD phase output
	out_avd_pvd::Int32 = 1 # activate writing .pvd file
	out_avd_ref::Int32 = 3 # AVD grid refinement factor
end


@kw_args struct Surf_Output
    out_surf::Int32         = 0 # activate surface output
    out_surf_pvd::Int32        = 0 # activate writing .pvd file
    out_surf_velocity::Int32    = 0
    out_surf_topography::Int32  = 0
    out_surf_amplitude::Int32   = 0
end

@kw_args struct Passive_Tracers_Output
    out_ptr::Int32              = 0    # activate  
    out_ptr_ID::Int32           = 0    # ID of the passive tracers
    out_ptr_phase::Int32        = 0    # phase of the passive tracers
    out_ptr_Pressure::Int32     = 0    # interpolated pressure
    out_ptr_Temperature::Int32  = 0    # temperature
    out_ptr_MeltFraction::Int32 = 0    # melt fraction computed using P-T of the marker 
    out_ptr_Active::Int32       = 0    # option that highlight the marker that are currently active
    out_ptr_Grid_Mf::Int32      = 0    # option that allow to store the melt fraction seen within the cell 
end

@kw_args struct Phase_Aggregator
    name::String     = "crust"  # phase aggregate output vector name
    numPhase::Int32 = 3      # number of phases to aggregate
    phaseID::Array{Int32}  = [1 5 15] # list of phase IDs to aggregate
end



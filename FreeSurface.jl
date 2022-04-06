using Parameters; 

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

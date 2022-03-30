#=
Main script
-> In the future, the path, path_S and Test_Name and name_pvtr are arguments
->
=#
#using GeoParams
include("Vtk_read_lamem.jl")
include("Slab_D_Functions.jl")


main()
C,D,F_,FSpec,TStep = initialize_test(); 





function _Test_Loop_Output(C::Coord_Model,D::DYN,F_::Field_Dyn,FSpec::Files_specification,TStep::Output_list)
"""
Input: 
C -> Coordinate of the test, Coord_model
D -> DYN Fields
F-> Fields 
FSpec ->Specification file
Tstep -> Output list
Output: 
-> strcture with relevant data (place holder)

Function: 
a) -> loop over the time step
b) -> Update the structure Surf, Dyn, and X 
c) -> Collect relevant data 
d) -> Update relevant data 
e) -> plot 
"""
for it=1:TStep.nTs





end




end


function initialize_test()
    # Function Initialization 

    Field = ["density [kg/m^3]",
             "visc_creep [Pa*s]",
             "velocity [cm/yr]",
            "dev_stress [MPa]",
             "strain_rate [1/s]",
            "j2_dev_stress [MPa]",
            "j2_strain_rate [1/s]",
            "tot_displ [km]",
            "_Continental_Lit_1 [ ]",
            "_Continental_Lit_2 [ ]",
            "_Oceanic_Lit_ [ ]",
            "_Continental_Crust_ [ ]"]
    F_ = Field_Dyn(Field)        
    path       = "/mnt/c/Users/Andrea Piccolo/Desktop"
    path_S     = "/mnt/c/Users/Andrea Piccolo/Dropbox/Bayreuth_Marcel/Tests"
    Test_Name  = "T_I0_VA18_V21_LM1_NE_FS"
    name_pvtr  = "SDet"
    FSpec      = Set_up_Fspec(path,path_S,name_pvtr,Test_Name);
    Tstep      = read_pvd_file(FSpec.file_pvd);
    #Find elegant and beauty way to do it
    C= _get_coordinate_(FSpec,Tstep,true,(-400.0,400.0),(0.0,0.0),(-400.0,0.0),1);
    D = DYN(Array{Float64}(undef,C.nz,C.nx));

    return C,D,F_,FSpec,TStep;

end

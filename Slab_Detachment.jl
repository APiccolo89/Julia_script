#=
Naming convection
-----------------------------------------------------------------------------------------------------------
OL (Output List)
Fsp (Test Specification)
D   (Dynamic properties: - the properties that are saved in the main computational grid)
F_d (Field Dynamic, field of the structure and relative dictionary)
S   (Surf properties - the properties that are saved in the free surface )
S_d (surf field, and relative dictionary)
Ph  (Phase properties - the properties that are saved in the free surface )
C    (Coordinate Model reference system, coordinates)
Sl    (Slab data type)
F_sl  (Slab properties and relative dictionary)
------------------------------------------------------------------------------------------------------------


=#
#using GeoParams
include("Vtk_read_lamem.jl")
include("Slab_D_Functions.jl")


main()
C,D,F_d,Fsp,OL = initialize_test(); 
_Test_Loop_Output(C,D,F_d,Fsp,OL)




function _Test_Loop_Output(C::Coord_Model,D::DYN,F_d::Field_Dyn,Fsp::Files_specification,OL::Output_list)
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
for it=1:OL.nTs
    _update_DYN!(D,Fsp,OL,C,it,F_d)




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

    dictionary = Dict("density [kg/m^3]"=>("ρ",1), 
                       "visc_creep [Pa*s]"=>("η_creep",1),
                       "velocity [cm/yr]"=>(["vx","vz"],3),
                       "dev_stress [MPa]"=>(["τ_xx","τ_zz","τ_xz"],9),
                       "strain_rate [1/s]"=>(["ϵ_xx","ϵ_zz","ϵ_xz"],9),
                       "j2_dev_stress [MPa]"=>("τ",1),
                       "j2_strain_rate [1/s]"=>("ϵ",1),
                       "_Continental_Lit_1 [ ]"=>("C₁",1),
                       "_Continental_Lit_2 [ ]"=>("C₂",1),
                       "_Oceanic_Lit_ [ ]"=>("Oₚ",1),
                       "_Continental_Crust_ [ ]"=>("Cont",1),
                       "tot_displ [km]" =>(["dx","dz"],3)
                       )

    F_d = Field_Dyn(Field,dictionary)        
    path       = "/mnt/c/Users/Andrea Piccolo/Desktop"
    path_S     = "/mnt/c/Users/Andrea Piccolo/Dropbox/Bayreuth_Marcel/Tests"
    Test_Name  = "T_I0_VA18_V21_LM1_NE_FS"
    name_pvtr  = "SDet"
    Fsp      = Set_up_Fspec(path,path_S,name_pvtr,Test_Name);
    OL      = read_pvd_file(Fsp.file_pvd);
    #Find elegant and beauty way to do it
    C= _get_coordinate_(Fsp,OL,true,(-400.0,400.0),(0.0,0.0),(-400.0,0.0),1);
    D = DYN(Array{Float64}(undef,length(C.z),length(C.x)));

    return C,D,F_d,Fsp,OL;

end

#=
Main script
-> In the future, the path, path_S and Test_Name and name_pvtr are arguments
->
=#
#using GeoParams

include("Vtk_read_lamem.jl")
include("Slab_D_Functions.jl")
path       = "C:\\Users\\Andrea Piccolo\\Desktop"
path_S     = "C:\\Users\\Andrea Piccolo\\Dropbox\\Bayreuth_Marcel\\Tests"
Test_Name  = "T_I0_VA18_V21_LM1_NE_FS"
name_pvtr  = "SDet"
FSpec      = Set_up_Fspec(path,path_S,name_pvtr,Test_Name);
Tstep      = read_pvd_file(FSpec.file_pvd);
#Find elegant and beauty way to do it 
C= _get_coordinate_(FSpec,true,(-400.0,400.0),(0.0,0.0),(-700.0,20.0),0)
#function main()

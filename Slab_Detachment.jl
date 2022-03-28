#=
Main script
=#
#using GeoParams



include("Vtk_read_lamem.jl")
include("Slab_D_Functions.jl")

path       = "C:\\Users\\Andrea Piccolo\\Desktop"
path_S     = "C:\\Users\\Andrea Piccolo\\Dropbox\\Bayreuth_Marcel\\Tests"
Test_Name  = "T_I0_VA18_V21_LM1_NE_FS"
name_pvtr  = "SDet"
FSpec      = Set_up_Fspec(path,path_S,name_pvtr,Test_Name)

Tstep = read_pvd_file(Fspec)

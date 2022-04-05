#=
Vtk utilities and general plotting tools
module read_vtk_LaMEM
=#
using Conda,PyCall,Plots,Printf,LinearAlgebra, PyPlot

vtk     =   pyimport("vtk")
dsa     =   pyimport("vtk.numpy_interface.dataset_adapter");



struct Files_specification
    path::String       #the path to the test
    pathS::String     # Path to save the output and database
    pathST::String
    file_pvd::String  # Path to the pvd file
    Test_Name::String #Test Name
    name_dyn::String   #name of the specific output
    name_surf:: String  #name of the specific output
    name_phase::String # name of the phase file
    name_PasTr::String # name of the passive tracers
end

struct Output_list    # Structure to handle the folder
    fLS:: Vector{String}      # Folder list of the test
    time:: Vector{Float64}    # time vector
    nTs:: Int64       # number of timestep
end

struct Coord_Model
    x:: Vector{Float64}   #x coordinate
    y:: Vector{Float64}   #y coordinate
    z:: Vector{Float64}   #z coordinate
    nx:: Int64            #number nodes along x
    ny:: Int64            #        ""         y
    nz:: Int64            #        ""         z
    ix:: Tuple{Int64,Int64}   # Touple with the ix_beg ix_end
    iy:: Tuple{Int64,Int64}
    iz:: Tuple{Int64,Int64}
    D2:: Bool                 # if true, we are dealing with 2D 
end

struct Field_Dyn
    Field:: Vector{String}
    DIC :: Dict 
    #Cmap::Vector{String}
end

mutable struct DYN
    vx:: Array{Float64}     #  Vel vector field (if 2D vel[1(x)2(z)])
    vz:: Array{Float64}
    dx:: Array{Float64}   #
    dz:: Array{Float64}   #
    τ ::Array{Float64}    # second invariant
    ϵ:: Array{Float64}     # second invariant
    τ_xx:: Array{Float64}  # Full tensor (if 2D -> T[1(xx)2(zz)3(xz),:,:])
    τ_zz:: Array{Float64}  # Full tensor (if 2D -> T[1(xx)2(zz)3(xz),:,:])
    τ_xz:: Array{Float64}  # Full tensor (if 2D -> T[1(xx)2(zz)3(xz),:,:])  
    ϵ_xx :: Array{Float64}
    ϵ_zz :: Array{Float64}
    ϵ_xz :: Array{Float64}
    Ψ::Array{Float64}
    η::Array{Float64}
    η_creep:: Array{Float64}
    ρ::Array{Float64}
    Oₚ::Array{Float64}
    C₁::Array{Float64}
    C₂::Array{Float64}
    Cont::Array{Float64}
    function DYN(b)
        """
        Initialize the structure with a random buffer  vector
        """
        return new(b, b, b, b, b, b ,b ,b ,b ,b ,b ,b ,b ,b ,b, b, b, b, b, b)
    end
end


struct SURF
    v  :: Array{Float64}   # velocity vector
    A :: Array{Float64}   # Amplitude
    T :: Array{Float64}   #Topography 
    function SURF(b)
        return new(b,b,b)
    end
end


function Set_up_Fspec(path::String,path_S::String,name_pvtr::String,Test_Name::String)
"""
    f(::String,::String,::String,::String)
    function that accept as argument the path containing the test, the path where to save the output
    the name of the .pvd file and the test name 
    Input: 
    path::String
    path_S::String
    name_pvtr::String
    Test_Name::String
    Output: 
    buf -> Type:: file specification where the path are saved and used 
    Function: 
    1) Generate the string containing the file specification (i.e. SDet (pvtr)+_surf.pvts)
    2) Generate the path where to look for the LaMEM output of a specific Test
    3) Generate the path where to save the output of the post_processing 
        3a) If the main save folder is not existing, create
        3b) If the folder of the test is not existing create
    Additional note: This script assumes that tests are grouped in projects and tests are performed for a general task
    i.e. if a user wants to explore the slab detachment perform a set of test that aims to solve this particolar task and each 
    tests is grouped accordingly. Path_S is the master folder of the group test: in this folder are contained - eventually - the 
    database of the tests, some picture that compare them, and the folder of each tests where the selected properties are saved in separated 
    folder
"""
    Filepvd    =   string(name_pvtr, ".pvd")
    Dyn        =   string(name_pvtr, ".pvtr")
    Surf       =   string(name_pvtr,"_surf" ,".pvtr")
    Phase      =   string(name_pvtr,"_phase",".pvts")
    PasTr      =   string(name_pvtr,"_passive_tracers",".pvtu")
    file_pvd   =   joinpath(path,Test_Name,Filepvd)
    path_ST    =   joinpath(path_S,Test_Name)
    #=
    Create the folder to save the output if it does not exist
    =#
    if isdir(path_S) == false
       mkdir(path_S)
    end

    if isdir(path_ST) == false
       mkdir(path_ST)
    end
    #=
    Fill the structure with the file specification
    =#
    buf = Files_specification(path,
                              path_S,
                              path_ST,
                              file_pvd,
                              Test_Name,
                              Dyn,
                              Surf,
                              Phase,
                              PasTr)
    return buf
end

function read_pvd_file(Fname::String)
   
    """
    Input: FName ::String
    Output: structure containing Time step, Time vector and total number of time
    steps
    buf = is a temporary buffer
    Function: 
    A) Create a local dictionary, and initialize two empty lists
    B) Read pvd file and collect: 1) Folder 2) Time 3) Total amount of folder
    C) Save in the data type Output_list
    """
    timestep_rgx = r"file=\".+?\/"
    time_rgx           = r"timestep=\".+?\""
    folder_list    = String[];
    time           = Float64[];
    f = open(Fname);
    lines = readlines(f);
    counter = 0;
    for l in lines
        Folder   = match(timestep_rgx,l);
        Timestep =  match(time_rgx,l);
        if Folder != nothing

            Time_s   = Timestep.match[11:end-1];
            f_s      = string(Folder.match[7:end-1]);
            println(counter);
            append!(time,parse(Float64,Time_s));
            push!(folder_list,f_s);
            counter += 1;
        end

    end

    close(f)
    buf = Output_list(folder_list,time,counter)
    return buf;
end


function _get_coordinate_(Fsp::Files_specification,OL::Output_list,D2::Bool,x_r::Tuple{Float64, Float64},y_r::Tuple{Float64, Float64},z_r::Tuple{Float64, Float64},istep::Int64)
   
    """
    f(fn,OL,D2,x_r,y_r,z_r,istep)
    Input
    fn ::File_specification  Information of the test
    OL :: Output_list        List of output 
    D2 ::Bool                Flag 2D or 3D 
    x_r :: Touple            zoom: -> delimiting area along x (y,z) to plot or do the post processing 
    .                        zoom:
    .                        zoom:
    istep ::Int64
    Output
    buf -> Coord_Model Data type (structure) Containing all the information related to the grid 
    Function: 
    A) Create the path to look for the file 
    B) Invoke the vtk utilities
    C) read the coordinates
    -> if 2D 
       select x-z direction
       if exist touple containing zoom
        ->create the new array 
        ->save the node where they start 
        ->update number of element
    D) Fill up the structure
    """

    fname =  joinpath(Fsp.path,Fsp.Test_Name,OL.fLS[istep],Fsp.name_dyn)

    reader  = vtk.vtkXMLPRectilinearGridReader()
    reader.SetFileName(fname)
    reader.Update()
    data    = reader.GetOutput()
    #=
    Extract the coordinates from the first timestep
    =#
    x               =   data.GetXCoordinates();
    y               =   data.GetYCoordinates();
    z               =   data.GetZCoordinates();
    x = convert(Array{Float64},x);
    y = convert(Array{Float64},y);
    z = convert(Array{Float64},z);
    nx              =   length(x);
    ny              =   length(y);
    nz              =   length(z);

    if D2 == true
        if(x_r[2]-x_r[1]>0.0)
            # Find all the nodes within the segment of investigation
            id = findall(b->(b>=x_r[1] && b<=x_r[2]),x);
            #find the first and last node of the segment
            ix = (id[1],id[end]);
            # update x, nx
            x  = x[id];
        else
            ix = (1,length(x))
        end
        if(z_r[2]-z_r[1]>0.0)

            # Find all the nodes within the segment of investigation
            id = findall(b->(b>=z_r[1] && b<=z_r[2]),z);
            #find the first and last node of the segment
            iz = (id[1],id[end]);
            # update x, nx
            z  = z[id];
        else
            iz = (1,length(z))
        end

        buf = Coord_Model(x,         # x coordinate
                          y,         # y coordinate
                          z,         # z coordinate
                          nx,        # nx along x
                          ny,         # ny along y
                          nz,        # nz along z
                          ix,        # izoom x
                          (0,0),     # izoom y
                          iz,        # izoom z
                          D2);       # Dimension

    else
    #=
    Place holder for the future
    =#

   end
 return buf;
end


function _update_DYN!(D::DYN,Fsp::Files_specification,OL::Output_list,C::Coord_Model,istep::Int64,F_d::Field_Dyn)

    fname =  joinpath(Fsp.path,Fsp.Test_Name,OL.fLS[istep],Fsp.name_dyn)

    reader  = vtk.vtkXMLPRectilinearGridReader()
    reader.SetFileName(fname)
    reader.Update()
    data    = reader.GetOutput()
    
    is = length(F_d.Field)
    for s=1:is 
        information = F_d.DIC[F_d.Field[s]]        
        buf      =   data.GetPointData().GetArray(F_d.Field[s]);
        b=_shape_array(buf,C,information);
        b      =  convert(Array{Float64},b);

        if information[2] == 3
            setfield!(D,Symbol(information[1][1]),b[1,:,:]);
            setfield!(D,Symbol(information[1][2]),b[2,:,:]);
        elseif information[2] == 9
            setfield!(D,Symbol(information[1][1]),b[1,:,:]);
            setfield!(D,Symbol(information[1][2]),b[3,:,:]);
            setfield!(D,Symbol(information[1][3]),b[2,:,:]);
        else
            setfield!(D,Symbol(information[1]),b);
        end
    
    end

end

function _shape_array(buf,C::Coord_Model,information)
    if information[2]>1
        nco = (information[2])
        if C.D2==true && nco==9
            b = Array{Float32}(undef,3,length(C.x),length(C.z));

        elseif C.D2==true && nco==3
            b = Array{Float32}(undef,2,length(C.x),length(C.z));

        end

        for i = 1:nco
            buf1 = buf[:,i];
            buf1 = reshape(buf1,(C.nx,C.ny,C.nz));
            if C.D2 == true 
                buf1 = buf1[:,1,:];
            end
            if i == 1 
                b[1,:,:] = buf1[C.ix[1]:C.ix[2],C.iz[1]:C.iz[2]];
            elseif i == 3
                b[2,:,:] = buf1[C.ix[1]:C.ix[2],C.iz[1]:C.iz[2]];
            elseif i == 9
                b[3,:,:] = buf1[C.ix[1]:C.ix[2],C.iz[1]:C.iz[2]];

            end
        end
        return b ;
    else 
        b    = reshape(buf,(C.nx,C.ny,C.nz));
        if C.D2 == true
            b = b[:,1,:];
            b = b[C.ix[1]:C.ix[2],C.iz[1]:C.iz[2]];
            return b;
            #place holder
        end
    end

end



function plot_grid_properties(D::DYN,Fsp::Files_specification,C::Coord_Model,field::Symbol,title::String,F::String,istep::Int64,OL::Output_list)

    fs = joinpath(Fsp.pathST,F)
    time = OL.time[istep]
    val = transpose(getfield(D,field))
    if log==true
        val = log10(val)
    end
    if isdir(fs) == false
        mkdir(fs)
    end
    # Contour the main unit 

    #

end
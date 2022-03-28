#=
Vtk utilities and general plotting tools
module read_vtk_LaMEM
=#

using PyCall
vtk     =   pyimport_conda("vtk","vtk")
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
    ix:: Vector{Int32}    # Touple with the ix_beg ix_end
    iy:: Vector{Int32}
    iz:: Vector{Int32}
end


function Set_up_Fspec(path,path_S,name_pvtr,Test_Name)
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

function read_pvd_file(Fname)
    #=
    Function that read the pvd files, and extract the folder that are needed for the post proccesing of the P-T-t path
    Introducing two regular expression that are not greedy ("?" allows to stop the regular expression at the first occurence of any characters)
    then opening a the file, then looping over the line of the pvd file to find the occurence.
    =#
    "
    Input: Filename
    Output: structure containing Time step, Time vector and total number of time
    steps
    buf = is a temporary buffer
    "
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
            println(counter);
        end

    end

    println("The total amount of timestep is $counter /n")
    close(f)
    buf = Output_list(folder_list,time,counter)
    return buf;
end

"
function _get_coordinate_(fn::Files_specification,2D::Bool==true,x_r::Float64=(0.0,0.0),y_r::Float64=(0.0,0.0),z_r::Float64=(0.0,0.0),istep::In64=0)
    #=
    Read the coordinate of the file. Adjust as a function of the
    fn   : Filespecification type structure
    2D   : bool (2D or 3D numerical experiment)
    x_r  : touple (xmin,xmax) if -> (0.0,0.0) No bounds
    y_r  : touple (ymin,ymax) ""
    z_r  : touple (zmin,zmax) ""
    =#

    fname =  joinpath(fn.path,fn.Test_Name,fn.fLS[1],fn.dyn)

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

    nx              =   length(x);
    ny              =   length(y);
    nz              =   length(z);

    if 2D == true
        if(x_r[2]-x_r[1]>0.0)
            # Find all the nodes within the segment of investigation
            id = findall(x->x>=x_r[1] & x<=x_r[2],x);
            #find the first and last node of the segment
            ix = (id[1],id[end]);
            # update x, nx
            x  = x[id];
            nx = length(x);
        end
        if(z_r[2]-z_r[1]>0.0)
            # Find all the nodes within the segment of investigation
            id = findall(x->x>=z_r[1] & x<=z_r[2],z);
            #find the first and last node of the segment
            iz = (id[1],id[end]);
            # update x, nx
            z  = z[id];
            nz = length(x);
        end

        buf = Coord_Model(x,         # x coordinate
                          (0.0,0.0), # y coordinate
                          z,         # z coordinate
                          nx,        # nx along x
                          0,         # ny along y
                          nz,        # nz along z
                          ix,        # izoom x
                          (0,0),     # izoom y
                          iz         # izoom z
                          );

    else
    #=
    Place holder for the future
    =#

 end
 return buf;
end
"

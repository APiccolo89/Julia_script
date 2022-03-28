#=
Vtk utilities and general plotting tools
module read_vtk_LaMEM
=#

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
            println(counter)
            append!(time,parse(Float64,Time_s))
            push!(folder_list,f_s)
            counter += 1;
            println(counter)
        end

    end

    println("The total amount of timestep is $counter /n")
    close(f)
    buf = Output_list(folder_list,time,counter)
    return buf;
end

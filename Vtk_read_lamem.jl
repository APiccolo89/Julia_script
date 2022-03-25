"""
Vtk utilities and general plotting tools
module read_vtk_LaMEM
"""
function read_pvd_file(Filename)
    """
    Function that read the pvd files, and extract the folder that are needed for the post proccesing of the P-T-t path
    Introducing two regular expression that are not greedy ("?" allows to stop the regular expression at the first occurence of any characters)
    then opening a the file, then looping over the line of the pvd file to find the occurence.
    """
    timestep_rgx = r"file=\".+?\/"
    time_rgx           = r"timestep=\".+?\""
    folder_list    = String[];
    time           = Float64[];
    f = open(Filename);
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
    return folder_list, time, counter;
end

"""
General outlines of what i'm going to do 
Project Julia Hackaton: convert my python script into julia
1) Produce fast plotting routine
2) Produce Data Base slab detachment
3) Parallel writing of database
-> Aims:
I want to write a script capable to post processes
several test at the same time, and update a Data
Base that contains the relevant data for my  project
to later use.
?) Object oriented
    -> How to produce classes of data that perform
    basic function?
    ->How to produce plot?
    ->How to optimize?
Structure:
-> Dyn (pvd file)     [General]
-> Surf (pvts file)   [General]
-> Ph_ag (-> Slab     [Specific]
          -> Lithos
          -> Crust)
"""
"""
OBJECT 1  (General)
-> path to the folder of thest (ptF)
-> path to the folder in which the output are saved (ptS)
-> Test Name
-> path where to save the specific output test (ptST)
OBJECT 2 (General)
-> lTs  : list of folder where find the input
-> nTs  : list of
-> time :
-> name_pvtr
-> name_pvts
OBJECT 3 (General)
->x    : real size of the model
->z
->nx   : number of nodes full sim
->nz
->xz   : area of interest
->zz
->nxz  : node of area of interest
->nzz
->Update : flag that warns if the size of the model is changing with time

OBJECT 4 Dyn (object size) (General)
-> List of properties
-> function that read input

OBJECT 5 FS (object size) (General/Specific)
-> List of properties
-> function that read input

OBJECT 5_b (specific)
-> list of properties and time-x maps data
-> function that saves .txt file

OBJECT 6 FS Slab (specific)
-> List of properties
-> Function that identify slab and save the averages values within the slab
-> Function that detect slab detachment and saves relevant properties
-> Function that plot the averages with time

Object 7 Non_dimensional (specific)
-> characteristic values
-> function that reads the input data and non dimensionalize the output

-> Function that creates the data set to use for additional analysis  (specific)



Specific (All function/object that belongs to the project of slab detachment)
General  (All the object that can be used outside the problem that I'm interested)
(Potential loop over test with potential to parallelize)
Main Block
Initialisation
-> paths
-> test name
-> creation folder
Call object 1
Call object 2
Call object 3
Call object 5_b
Call object 6
Call object 7

Loop over the time step

-> retrieve FS data
-> retrieve DYN data
-> Update Slab & plot 1D plot
-> Update Lithosphere Object
-> Update free surface & plot 1D plot
---------------------------------
-> Plot time-x/z maps
-> Update data base
-> Finish test X
"""

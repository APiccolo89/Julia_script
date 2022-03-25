"""
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
1st Block
OBJECT 1
-> path to the folder of thest (ptF)
-> path to the folder in which the output are saved (ptS)
-> Test Name
-> path where to save the specific output test (ptST)
OBJECT 2
-> lTs  : list of folder where find the input
-> nTs  : list of
-> name_pvtr
-> name_pvts
OBJECT 3
->x    : real size of the model
->z
->nx   : number of nodes full sim
->nz
->xz   : area of interest
->zz
->nxz  : node of area of interest
->nzz

OBJECT 4 Dyn (object size)
-> List of properties
-> function that read input

OBJECT 5 FS (object size)
-> List of properties
-> function that read input
-> function that saves .txt file

OBJECT 6 FS Slab
-> List of properties
-> Function that identify slab and save the




"""

#include("julia_first_setup.jl");                                     # packages 

include("DDAfunctions.jl");                                           # set of Julia functions

include("make_data_7_systems.jl");                                    # make data
include("run_DDA_Roessler.jl");                                       # run DDA
include("Roessler_ShowResults.jl");                                   # show results


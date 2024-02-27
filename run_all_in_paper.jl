include("DDAfunctions.jl");                                            # set of Julia functions

WL=2000;WS=500;WN=500;                                                 # assign window parameters
#WL=4000;WS=1000; WN=2000;                                             # other window settings to try

FN=@sprintf("CD_DDA/CD_DDA_data_NoNoise__WL%d_WS%d_WN%d.ascii",
             WL,WS,WN);                                                # noise free data file
if !isfile(FN)
   include("make_data_7_systems.jl");
end

include("run_DDA_NoNoise.jl");

include("run_DDA_15dB.jl");


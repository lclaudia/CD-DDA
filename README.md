This software reproduces the Roessler part of ``Network-motif delay differential analysis of brain activity during seizures", Chaos 33(12):123136; 2023. https://doi.org/10.1063/5.0165904

It was tested on linux and julia 1.10.1 (https://julialang.org/).

Use it at your own risk.


Do not forget to chmod for the executables:

chmod +x i_ODE_general_BIG
chmod +x run_DDA_ASCII


Single_Roessler_Example.jl   is an example on how to integrate a single Roessler system

run_all_in_paper.jl          does all computations in the paper
make_data_7_systems.jl       makes the data used in the paper 
run_DDA_NoNoise.jl           reproduces the analysis of noise free data 
run_DDA_15dB.jl              reproduces the analysis data with added white noise 


In case you want to try the software on other data:

I am happy to answer questions if you use the software on human data or data away from biology. 
For data from animals I do not want to participate in any way. Please respect life and be ethical! Animals have feelings and fears!

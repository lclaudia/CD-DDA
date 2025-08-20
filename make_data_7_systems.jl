include("DDAfunctions.jl");                                            # set of Julia functions
#
WL=2000;WS=500;WN=500;                                                 # assign window parameters
##WL=4000;WS=1000; WN=2000; 

DDA_DIR  = "DDA";  dir_exist(DDA_DIR);                                 # folders
DATA_DIR = "DATA"; dir_exist(DATA_DIR); 
FIG_DIR  = "FIG";  dir_exist(FIG_DIR); 

NrSyst=7;                                                              # 7 coupled systems
ROS=[[0  0 2];                                                         # single Roessler system
     [0  0 3];
     [1  0 1];
     [1  0 2];
     [2  0 0];
     [2  0 3];
     [2  1 3]
    ];
(MOD_nr,DIM,ODEorder,P) = make_MOD_nr(ROS,NrSyst);                     # encoding of the 7 Roessler systems
                                                                       # function defined in DDAfunctions.jl

a123= 0.21;                                                            # model parameters
a456= 0.20;
a7  = 0.18;
b1  = 0.2150;
b2  = 0.2020; 
b3  = 0.2041; 
b4  = 0.4050; 
b5  = 0.3991; 
b6  = 0.4100; 
b7  = 0.5000;
c   = 5.7;
c7  = 6.8;
MOD_par=[
         -1 -1 1 a123  b1 -c  1 
         -1 -1 1 a123  b2 -c  1
         -1 -1 1 a123  b3 -c  1
         -1 -1 1 a456  b4 -c  1
         -1 -1 1 a456  b5 -c  1
         -1 -1 1 a456  b6 -c  1
         -1 -1 1 a7    b7 -c7 1
        ];
MOD_par=reshape(MOD_par',size(ROS,1)*NrSyst)';

FromTo1=[];

FromTo2=[[4 0  0 1   7 0  0 1];                                       # from 4th system 1st Eq. variable 1 
                                                                      #   to 7th system 1st Eq. variable 1 
         [5 0  0 1   7 0  0 1];  
         [6 0  0 1   7 0  0 1]];

FromTo3=[[7 0  0 1   4 0  0 1];  
         [7 0  0 1   5 0  0 1];  
         [7 0  0 1   6 0  0 1]];
         
I2=make_MOD_nr_Coupling(FromTo2,DIM,P);                               # MOD_nr part for coupling; case (ii)
I3=make_MOD_nr_Coupling(FromTo3,DIM,P);                               # MOD_nr part for coupling; case (iii)
                                                                      # function defined in DDAfunctions.jl
II=[Int[],I2,I3];
        
epsilon=0.15;                                                         # coupling strength

MOD_par_add2=repeat([epsilon -epsilon],size(FromTo2,1),1)'[:]';       # MOD_par for coupling part
MOD_par_add3=repeat([epsilon -epsilon],size(FromTo3,1),1)'[:]';       # MOD_par for coupling part

MOD_par_add=[Float64[], MOD_par_add2, MOD_par_add3];

TAU=[32 9]; TM=maximum(TAU); dm=4;                                    # DDA parameters

LL=[WS*(WN-1)+WL+TM+dm; 
    WS*WN; 
    WS*WN+dm-1]; 

DELTA=2;                                                              # every second data point
CH_list=1:DIM:DIM*NrSyst;                                             # only x
TRANS=20000;
dt=0.05;

CASE=["i";"ii";"iii"];

for n_CASE=1:length(CASE)
    FN=@sprintf("%s%sCD_DDA_data__WL%d_WS%d_WN%d__case_%s.ascii", 
                 DATA_DIR,SL,WL,WS,WN,CASE[n_CASE]);                  # data file
    if !isfile(FN)
       X0 = rand(DIM*NrSyst,1);                                       # initial conditions
       if length(II[n_CASE])>0
          M1=[MOD_nr II[n_CASE]]; M2=[MOD_par MOD_par_add[n_CASE]];
       else
          M1=MOD_nr; M2=MOD_par;
       end
       integrate_ODE_general_BIG(M1,M2,                               # encoding of the coupled systems
                                 dt,                                  # step size of num. integration
                                 LL[n_CASE],                          # length 
                                 DIM*NrSyst,ODEorder,X0,              # parameters
                                 FN,                                  # output file
                                 CH_list,DELTA,                       # only x; every second point
                                 TRANS);                              # transient
    end
end

global X=Matrix{Any}(undef,1,NrSyst);
for n_CASE=1:length(CASE)
    FN=@sprintf("%s%sCD_DDA_data__WL%d_WS%d_WN%d__case_%s.ascii", 
                 DATA_DIR,SL,WL,WS,WN,CASE[n_CASE]);                  # data file
    if n_CASE == 1
       global X=readdlm(FN);
    else
       global X=vcat(X,readdlm(FN));
    end
end

SG = plot(layout = (length(CASE),NrSyst),size=(2100,800));            # make plot of delay embeddings
for n_CASE=1:length(CASE)
    for n_SYST=1:NrSyst
        plot!(SG,subplot=(n_CASE-1)*NrSyst+n_SYST,
              X[((20000:24000) .+ (n_CASE-1)*LL[n_CASE]),      n_SYST],
              X[((20000:24000) .+ (n_CASE-1)*LL[n_CASE]) .- 10,n_SYST],
              legend=false
             )
        end
    end
display(SG)
savefig(SG,@sprintf("%s%sRoessler_7syst_NoNoise.png",DATA_DIR,SL));

###  add noise  ###

SNR=15;                                                                # signal-to-noise ratio in dB

Y=X .* 1;                                                                 
for n_CASE=1:length(CASE)
    for n_SYST=1:NrSyst                                                # add noise
        range = (1:LL[n_CASE]) .+ (n_CASE-1)*LL[n_CASE];
        Y[range,n_SYST] = add_noise(Y[range,n_SYST],SNR); 
    end
end

SG = plot(layout = (length(CASE),NrSyst),size=(2100,800));            # make plot of delay embeddings
for n_CASE=1:length(CASE)
    for n_SYST=1:NrSyst
        plot!(SG,subplot=(n_CASE-1)*NrSyst+n_SYST,
              Y[((20000:24000) .+ (n_CASE-1)*LL[n_CASE]),      n_SYST],
              Y[((20000:24000) .+ (n_CASE-1)*LL[n_CASE]) .- 10,n_SYST],
              legend=false
             )
        end
    end
display(SG)

display(SG)
savefig(SG,@sprintf("%s%sRoessler_7syst_15dB.png",DATA_DIR,SL));



FN=@sprintf("%s%sCD_DDA_data_NoNoise__WL%d_WS%d_WN%d.ascii",DATA_DIR,SL,WL,WS,WN);
writedlm(FN, map(number_to_string, X),' ');                            # save data 

FN=@sprintf("%s%sCD_DDA_data_15dB__WL%d_WS%d_WN%d.ascii",DATA_DIR,SL,WL,WS,WN);
writedlm(FN, map(number_to_string, Y),' ');                            # save data 



Y = nothing; GC.gc();
X = nothing; GC.gc();

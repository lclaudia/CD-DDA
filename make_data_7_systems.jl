include("DDAfunctions.jl");                                           # set of Julia functions

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

a123=0.21;                                                             # model parameters
a456=0.20;
a7  =0.18;
b1 = 0.2150;
b2 = 0.2020; 
b3 = 0.2041; 
b4 = 0.4050; 
b5 = 0.3991; 
b6 = 0.4100; 
b7 = 0.5000;
c =5.7;
c7=6.8;
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

FromTo2=[[4 0  0 1   7 0  0 1];                                        # from 4th system 1st Eq. variable 1 
                                                                       #   to 7th system 1st Eq. variable 1 
         [5 0  0 1   7 0  0 1];  
         [6 0  0 1   7 0  0 1]];

FromTo3=[[7 0  0 1   4 0  0 1];  
         [7 0  0 1   5 0  0 1];  
         [7 0  0 1   6 0  0 1]];
         
I2=make_MOD_nr_Coupling(FromTo2,DIM,P);                                # MOD_nr part for coupling; case (ii)
I3=make_MOD_nr_Coupling(FromTo3,DIM,P);                                # MOD_nr part for coupling; case (iii)
                                                                       # function defined in DDAfunctions.jl
        
epsilon=0.15;                                                          # coupling strength

MOD_par_add_23=[epsilon -epsilon epsilon -epsilon epsilon -epsilon];   # MOD_par for coupling part

TAU=[32 9]; TM=maximum(TAU); dm=4;                                     # DDA parameters
WL=2000;WS=500; 

WN=100;                                                                # assign window number for each case
L1=WS*(WN-1)+WL+TM+dm;                                                 # ajust integration length to have 
L2=WS*WN;                                                              # equal number of windows for each case
L3=WS*WN+dm;

TRANS=20000;
dt=0.05;

X0 = rand(DIM*NrSyst,1);                                               # initial conditions

X=integrate_ODE_general_BIG(MOD_nr,MOD_par,                            # encoding of the 7 systems, case (i)
                            dt,                                        # step size of num. integration
                            L1*2,                                      # length * 2
                            DIM*NrSyst,ODEorder,X0,                    # parameters
                            TRANS);                                    # transient

X0=X[end,:];                                                           # initial conditions for (ii)

Y=integrate_ODE_general_BIG([MOD_nr I2],[MOD_par MOD_par_add_23],      # encoding of the coupled systems (ii)
                            dt,                                        # step size of num. integration
                            L2*2,                                      # length * 2
                            DIM*NrSyst,ODEorder,X0,                    # parameters
                            0);                                        # no transient

X=[X; Y];                                                              # concatenate (i) and (ii)
Y = nothing; GC.gc();                                                  # clear variable Y
X0=X[end,:];                                                           # initial conditions for (iii)

Y=integrate_ODE_general_BIG([MOD_nr I3],[MOD_par MOD_par_add_23],      # encoding of the coupled systems (iii)
                            dt,                                        # step size of num. integration
                            L3*2,                                      # length * 2
                            DIM*NrSyst,ODEorder,X0,                    # parameters
                            0);                                        # no transient

X=[X;Y];                                                               # concatenate (i), (ii), and (iii)
Y = nothing; GC.gc();                                                  # clear variable Y

X = X[1:2:end,1:3:end];                                                # take every second data point 
                                                                       #   and only u_{1,n}


SG = plot(layout = (3,7),size=(2100,800));                             # make plot of delay embeddings
for k1=1:3
    for k2=1:7
        plot!(SG,subplot=(k1-1)*7+k2,
              X[((20000:24000) .+ (k1-1)*L2),k2],
              X[((20000:24000) .+ (k1-1)*L2) .- 10,k2],
              legend=false
             )
        end
    end
display(SG)
savefig("Roessler_7syst_NoNoise.pdf")

FN="CD_DDA_data_NoNoise.ascii";                                        # noise free data file
writedlm(FN, map(number_to_string, X),' '); 

###  add noise  ###

SNR=15;                                                                # signal-to-noise ratio in dB
Y=X;                                                                   # Note, that
                                                                       # Julia arrays are not copied when 
                                                                       # assigned to another variable. 
                                                                       # After A = B, changing elements of
                                                                       # B will modify A as well.
for k=1:size(X,2)                                                      # add noise
    Y[:,k]=add_noise(Y[:,k],SNR); 
end

SG = plot(layout = (3,7),size=(2100,800));                             # make plot of delay embeddings
for k1=1:3
    for k2=1:7
        plot!(SG,subplot=(k1-1)*7+k2,
              Y[((20000:24000) .+ (k1-1)*L2),k2],
              Y[((20000:24000) .+ (k1-1)*L2) .- 10,k2],
              legend=false
             )
        end
    end
display(SG)
savefig("Roessler_7syst_15dB.pdf")

FN="CD_DDA_data_15dB.ascii";  
writedlm(FN, map(number_to_string, Y),' ');                            # save data 

Y = nothing; GC.gc();
X = nothing; GC.gc();

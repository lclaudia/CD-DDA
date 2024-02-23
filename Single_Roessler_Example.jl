include("DDAfunctions.jl");                                           # set of Julia functions

NrSyst=1;                                                             # 1 single system
ROS=[[0  0 2];                                                        # single Roessler system
     [0  0 3];
     [1  0 1];
     [1  0 2];
     [2  0 0];
     [2  0 3];
     [2  1 3]
    ];
 (MOD_nr,DIM,ODEorder,P) = make_MOD_nr(ROS,NrSyst);                    # encoding of the Roessler system
                                                                       # function defined in DDAfunctions.jl

a=.2; c=5.7;                
dt=.05; X0=rand(DIM,1);                                                # choice of parameters
L=10000; TRANS=5000;                                                   # integration length and transient

b=0.45;
MOD_par=[-1 -1 1 a b -c 1];                                            # parameters
# DO NOT FORGET: "chmod +x i_ODE_general_BIG" in linux!
X = integrate_ODE_general_BIG(MOD_nr,MOD_par,dt,
                              L,DIM,ODEorder,X0,TRANS);                # integrate system
                                                                       # function defined in DDAfunctions.jl
plot(X[:,1],X[:,2],X[:,3],                                             # plot the attractor
     color=:blue,legend=false,
     xlabel=L"x",ylabel=L"y",zlabel=L"z")

plot!(size=(500,500))  
savefig("Roessler_0.45.pdf")

b=1; 
MOD_par=[-1 -1 1 a b -c 1];                                            # parameters
X = integrate_ODE_general_BIG(MOD_nr,MOD_par,dt,
                              L,DIM,ODEorder,X0,TRANS);                # integrate system
                                                                       # function defined in DDAfunctions.jl
plot(X[:,1],X[:,2],X[:,3],                                             # plot the attractor
     color=:blue,legend=false,
     xlabel=L"x",ylabel=L"y",zlabel=L"z")
plot!(size=(500,500))
savefig("Roessler_1.pdf")

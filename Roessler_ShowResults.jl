#include("run_DDA_Roessler.jl");
include("DDAfunctions.jl");                                            # set of Julia functions

WL=2000;WS=500;                                                        # window length and shift
global WN=500;

DDA_DIR  = "DDA";  dir_exist(DDA_DIR);                                 # folders
DATA_DIR = "DATA";
FIG_DIR  = "FIG";  dir_exist(FIG_DIR); 

NOISE=["NoNoise";"15dB"];

NrSyst=7; DIM=3;                                                       # 7 3D systems 
NrCH=NrSyst; CH=collect(1:NrCH);                                       # x-components of 7 systems are channels

DDAmodel=[[0 0 1];                                                     #    a_1 v_1 +
          [0 0 2];                                                     #    a_2 v_2 +
          [1 1 1]];                                                    #    a_3 v_1^3
(MODEL, L_AF, DDAorder)=make_MODEL(DDAmodel);                          # DDA model encoding for DDA code

LIST=collect(combinations(CH,2));                                      # pairwise combinations of channels
LL1=vcat(LIST...)'; LIST=reduce(hcat,LIST)';

for n_NOISE = 1:length(NOISE)
    noise=NOISE[n_NOISE];

    FN_DDA=@sprintf("%s%s%s__WL%d_WS%d_WN%d.DDA",
                     DDA_DIR,SL,noise,WL,WS,WN);                       # DDA file

    E=fill(NaN,WN,NrSyst,NrSyst,3);                                    # dynamical ergodicity matrix
    C=fill(NaN,WN,NrSyst,NrSyst,3);                                    # causality matrix

    ST=readdlm(join([FN_DDA,"_ST"]));                                  # read ST DDA output file
    T=ST[:,1:2]; ST=ST[:,3:end];                                       # first 2 numbers in each line are 
                                                                       #    start and end of window
    ST=ST[:,L_AF:L_AF:end];                                            # need only error
    ST=reshape(ST,WN,3,NrSyst);                                        # reshape matrix: 3 cases, 7 systems

    CT=readdlm(join([FN_DDA,"_CT"]));                                  # read CT DDA output file
    CT=CT[:,3:end];                                                    # first 2 numbers in each line are 
                                                                       #    start and end of window
    CT=CT[:,L_AF:L_AF:end];                                            # need only error
    CT=reshape(CT,WN,3,size(LIST,1));                                  # reshape matrix: 3 cases,  
                                                                       #    length(LIST) combinations
    for l=1:size(LIST,1)
        ch1=LIST[l,1];ch2=LIST[l,2];
        E[:,ch1,ch2,:] = abs.( dropdims(mean(ST[:,:,[ch1,ch2]],dims=3),dims=3) ./ CT[:,:,l] .- 1 );
        E[:,ch2,ch1,:] = E[:,ch1,ch2,:];
    end    

    CD=readdlm(join([FN_DDA,"_CD_DDA_ST"]));                           # read CD-DDA output file
    CD=CD[:,3:end];                                                    # first 2 numbers in each line are 
                                                                       #    start and end of window
    CD=reshape(CD,WN,3,2,size(LIST,1));                                # reshape matrix: 3 cases,  
                                                                       #    length(LIST) combinations
    for l=1:size(LIST,1)
        ch1=LIST[l,1];ch2=LIST[l,2];
        C[:,ch1,ch2,:] = CD[:,:,2,l];
        C[:,ch2,ch1,:] = CD[:,:,1,l];
    end    

###   plot results

    l=@layout[a{0.7h} ; b c d];                                  
    SG = plot(layout = l,size=(1500,1500));
    CG= cgrad([:white, RGB(1,0.97,0.86), RGB(0.55,0.27,0.07)],
              [0,0.1],scale=:linear);

    e=reshape(E,size(E,1),NrSyst^2,3);
    e=permutedims(e,[1,3,2]);
    e=reshape(e,WN*3,NrSyst^2)';

    N=tril(reshape(1:NrSyst^2,NrSyst,NrSyst),-1)[:]; 
    N=filter(x -> x != 0, N);
    S=[join(string.(x), " ") for x in eachrow(LIST)];

    heatmap!(SG,subplot=1,
             e[N,:],
             c=CG,
             xtickfont=font(12), ytickfont=font(12),
             colorbar = true,
             yticks=(1:21,S),
             xticks=(100," ")
             ) 

    e=dropdims(mean(E[20:end,:,:,:],dims=1),dims=1);           
    for k=1:3
            heatmap!(SG,
                subplot = k+1,
                e[:,:,k],
                c=:jet,
                colorbar = true,
                xtickfont=font(12), ytickfont=font(12),
                xlims=(0.5, 7.5), ylims=(0.5, 7.5),
                title=@sprintf("(%d)",k),
                aspect_ratio = :equal
                )
    end
    display(SG)

    print("Make png file and continue? ");
    readline()
    savefig(SG,@sprintf("%s%sE__WL%d_WS%d_WN%d_%s.png",FIG_DIR,SL,WL,WS,WN,noise));   

###

    l=@layout[a{0.5h} ; b c d; e f g];                                 
    SG = plot(layout = l,size=(1500,1500));
    CG= cgrad([:white, RGB(1,0.97,0.86), RGB(0.55,0.27,0.07)],
              [0,0.1],scale=:linear);

    c=reshape(C,size(C,1),NrSyst^2,3);
    c=permutedims(c,[1,3,2]);
    c=reshape(c,WN*3,NrSyst^2)';
    c .= c .- minimum(filter(!isnan,c[:]));                                       
    c .= c ./ maximum(filter(!isnan,c[:]));
    
    N=setdiff(1:NrSyst^2, diagind(C[1,:,:,1]));
    S=collect(permutations(CH,2));
    S=reduce(hcat,S)';
    S=[join(string.(x), " ") for x in eachrow(S)];

    heatmap!(SG,subplot=1,
             c[N,:],
             c=CG,
             xtickfont=font(12), ytickfont=font(12),
             colorbar = true,
             yticks=(1:42,S),
             xticks=(100," ")
             ) 

    c=dropdims(mean(C[50:end,:,:,:],dims=1),dims=1);                    
    c .= c .- minimum(filter(!isnan,c[:]));                             
    c .= c ./ maximum(filter(!isnan,c[:]));

    CG= cgrad([RGB(0.9,0.9,0.9), RGB(0.3,.3,0.3), :magenta, :cyan],    
              [0.0, 0.25, 0.2501, 0.635, 1],scale=:linear);             

    h = [heatmap!(SG,subplot=k+1,
                  c[:,:,k],                                           
                  c = CG, clim=(0,1),
                  colorbar = true,
                  title=@sprintf("(%d)",k),
                  xtickfont=font(12), ytickfont=font(12),
                  xlims=(0.5, 7.5), ylims=(0.5, 7.5),
                  aspect_ratio = :equal
                )
         for k in 1:3];

    q=reshape(c,NrSyst*NrSyst,3);                                       
    for k=1:3
        q[:,k] .= q[:,k] .- minimum(filter(!isnan,q[:,k]));
        q[:,k] .= q[:,k] ./ maximum(filter(!isnan,q[:,k]));
    end
    q=reshape(c,NrSyst,NrSyst,3);

    GR.setarrowsize(0.5);

    MS = [1,1,1,2,2,2,3];
    colors = [colorant"plum2", colorant"mistyrose1", colorant"lavender"];

    for k=1:3
        A=q[:,:,k];
        A[A .< 0.25] .= 0;
        A[isnan.(A)] .= 0;
    
        graphplot!(SG,subplot=k+4,A,
                  method=:circular,nodeshape=:circle,
                  names=1:7,
                  markersize=0.15,
                  fontsize=20,
                  linewidth=3,
                  linealpha=1,
                  markercolor = colors[MS],
                  nodestrokecolor=colors[MS],
                  arrow=arrow(:closed,10),
                  )
    end
    display(SG)

    print("Make png file and continue? ");
    readline()
    savefig(SG,@sprintf("%s%sC__WL%d_WS%d_WN%d_%s.png",FIG_DIR,SL,WL,WS,WN,noise));   

###

    CE=C .* E;                                                             # causality * ergodicity

    l=@layout[a{0.5h} ; b c d; e f g];                                  
    SG = plot(layout = l,size=(1500,1500));
    CG= cgrad([:white, RGB(1,0.97,0.86), RGB(0.55,0.27,0.07)],
              [0,0.1],scale=:linear);
    
    c=reshape(CE,size(CE,1),NrSyst^2,3);
    c=permutedims(c,[1,3,2]);
    c=reshape(c,WN*3,NrSyst^2)';
    c .= c .- minimum(filter(!isnan,c[:]));                            
    c .= c ./ maximum(filter(!isnan,c[:]));

    heatmap!(SG,subplot=1,
             c[N,:],
             c=CG,
             xtickfont=font(12), ytickfont=font(12),
             colorbar = true,
             yticks=(1:42,S),
             xticks=(100," ")
             ) 

    c=dropdims(mean(CE[50:end,:,:,:],dims=1),dims=1);                  
    c .= c .- minimum(filter(!isnan,c[:]));                            
    c .= c ./ maximum(filter(!isnan,c[:]));

    CG= cgrad([RGB(0.9,0.9,0.9), RGB(0.3,.3,0.3), :magenta, :cyan],     
              [0.0, 0.25, 0.2501, 0.635, 1],scale=:linear);             

    h = [heatmap!(SG,subplot=k+1,
                  c[:,:,k],                                             
                  c = CG, clim=(0,1),
                  colorbar = true,
                  title=@sprintf("(%d)",k),
                  xtickfont=font(12), ytickfont=font(12),
                  xlims=(0.5, 7.5), ylims=(0.5, 7.5),
                  aspect_ratio = :equal
                )
         for k in 1:3];

    q=reshape(c,NrSyst*NrSyst,3);                                      
    for k=1:3
        q[:,k] .= q[:,k] .- minimum(filter(!isnan,q[:,k]));
        q[:,k] .= q[:,k] ./ maximum(filter(!isnan,q[:,k]));
    end
    q=reshape(c,NrSyst,NrSyst,3);
    
    for k=1:3
        A=q[:,:,k];
        A[A .< 0.25] .= 0;
        A[isnan.(A)] .= 0;
    
        graphplot!(SG,subplot=k+4,A,
                  method=:circular,nodeshape=:circle,
                  names=1:7,
                  markersize=0.15,
                  fontsize=20,
                  linewidth=3,
                  linealpha=1,
                  markercolor = colors[MS],
                  nodestrokecolor=colors[MS],
                  arrow=arrow(:closed,10),
                  init=0
                  )
    end
    display(SG)

    print("Make png file and continue? ");
    readline()
    savefig(SG,@sprintf("%s%sCE__WL%d_WS%d_WN%d_%s.png",FIG_DIR,SL,WL,WS,WN,noise));  
end


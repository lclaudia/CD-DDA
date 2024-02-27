#include("DDAfunctions.jl");                                            # set of Julia functions
#
#WL=2000;WS=500;WN=500;                                                 # window length and shift
##WL=4000;WS=1000;WN=2000;

NrSyst=7; DIM=3;                                                       # 7 3D systems 

NrCH=NrSyst; CH=collect(1:NrCH);                                       # x-components of 7 systems are channels

LIST=collect(combinations(CH,2));                                      # pairwise combinations of channels
LL1=vcat(LIST...)';
LIST=reduce(hcat,LIST)';

nr_delays=2; dm=4;                                                     # DDA parameters
                                                                       # encoding of DDA model
                                                                       # \dot{v} =  
DDAmodel=[[0 0 1];                                                     #    a_1 v_1 +
          [0 0 2];                                                     #    a_2 v_2 +
          [1 1 1]];                                                    #    a_3 v_1^3
(MODEL, L_AF, DDAorder)=make_MODEL(DDAmodel);                          # DDA model encoding for DDA code



TAU=[32 9]; TM=maximum(TAU);                                           # delays

dir_exist("CD_DDA");                                                   # make folder

FN_data=@sprintf("CD_DDA/CD_DDA_data_15dB__WL%d_WS%d_WN%d.ascii",
                  WL,WS,WN);                                           # noise free data file
FN_DDA=@sprintf("CD_DDA/CD_DDA_data_15dB__WL%d_WS%d_WN%d.DDA",
                 WL,WS,WN);                                            # DDA file

if !isfile(join([FN_DDA,"_ST"]))
   CMD = "./run_DDA_ASCII -ASCII";                                     # DDA executable
   CMD = "$CMD -MODEL $(join(MODEL," "))";                             # model
   CMD = "$CMD -TAU $(join(TAU," "))";                                 # delays       
   CMD = "$CMD -dm $dm -order $DDAorder -nr_tau $nr_delays";           # DDA parameters
   CMD = "$CMD -DATA_FN $FN_data -OUT_FN $FN_DDA";                     # input and output files
   CMD = "$CMD -WL $WL -WS $WS";                                       # window length and shift
   CMD = "$CMD -SELECT 1 1 0 0";                                       # ST and CT DDA
   CMD = "$CMD -WL_CT 2 -WS_CT 2";                                     # take 2 channels for CT DDA
   CMD = "$CMD -CT_CH_list $(join(LL1," "))";                          # list of channel pairs for CT DDA
   
   run(Cmd(string.(split(CMD, " "))));                                 # run ST and CT DDA
end

ST=readdlm(join([FN_DDA,"_ST"]));                                      # read ST DDA output file
T=ST[:,1:2]; ST=ST[:,3:end];                                           # first 2 numbers in each line are 
                                                                       #    start and end of window
WN=Int(size(T,1)/3);                                                   # window number = number of lines
ST=ST[:,L_AF:L_AF:end];                                                # need only error
ST=reshape(ST,WN,3,7);                                                 # reshape matrix: 3 cases, 7 systems

CT=readdlm(join([FN_DDA,"_CT"]));                                      # read CT DDA output file
CT=CT[:,3:end];                                                        # first 2 numbers in each line are 
                                                                       #    start and end of window
CT=CT[:,L_AF:L_AF:end];                                                # need only error
CT=reshape(CT,WN,3,size(LIST,1));                                      # reshape matrix: 3 cases,  
                                                                       #    length(LIST) combinations

E=fill(NaN,WN,7,7,3);                                                  # dynamical ergodicity matrix
for l=1:size(LIST,1)
    ch1=LIST[l,1];ch2=LIST[l,2];
    E[:,ch1,ch2,:] = abs.( dropdims(mean(ST[:,:,[ch1,ch2]],dims=3),dims=3) ./ CT[:,:,l] .- 1 );
    E[:,ch2,ch1,:] = E[:,ch1,ch2,:];
end    


l=@layout[a{0.7h} ; b c d];                                            # plot results
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


e=dropdims(mean(E[20:end,:,:,:],dims=1),dims=1);                       # mean over windows

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
                                                                       
print("Make pdf file and continue? ");
readline()
                                                                       
savefig(@sprintf("PDFs/CD_DDA_E_15dB__WL%d_WS%d_WN%d.pdf",WL,WS,WN));      # save figure as pdf

###

if !isfile(join([FN_DDA,"_CD_DDA_ST"]))
   CMD = "./run_DDA_ASCII -ASCII";                                     # DDA executable
   CMD = "$CMD -MODEL $(join(MODEL," "))";                             # model
   CMD = "$CMD -TAU $(join(TAU," "))";                                 # delays       
   CMD = "$CMD -dm $dm -order $DDAorder -nr_tau $nr_delays";           # DDA parameters
   CMD = "$CMD -DATA_FN $FN_data -OUT_FN $FN_DDA";                     # input and output files
   CMD = "$CMD -WL $WL -WS $WS";                                       # window length and shift
   CMD = "$CMD -SELECT 0 0 1 0";                                       # run CD-DDA
   CMD = "$CMD -PAIRS $(join(LL1," "))";                               # pairs for CD-DDA
   
   run(Cmd(string.(split(CMD, " "))));                                 # run CD-DDA
end
                                                                        
CD=readdlm(join([FN_DDA,"_CD_DDA_ST"]));                               # read CD-DDA output file
CD=CD[:,3:end];                                                        # first 2 numbers in each line are 
                                                                       #    start and end of window
CD=reshape(CD,WN,3,2,size(LIST,1));                                    # reshape matrix: 3 cases,  
                                                                       #    length(LIST) combinations

C=fill(NaN,WN,7,7,3);                                                  # causality matrix
for l=1:size(LIST,1)
    ch1=LIST[l,1];ch2=LIST[l,2];
    C[:,ch1,ch2,:] = CD[:,:,2,l];
    C[:,ch2,ch1,:] = CD[:,:,1,l];
end    
C[isnan.(C)].=0;


l=@layout[a{0.5h} ; b c d; e f g];                                     # plot results
SG = plot(layout = l,size=(1500,1500));
CG= cgrad([:white, RGB(1,0.97,0.86), RGB(0.55,0.27,0.07)],
          [0,0.1],scale=:linear);

c=reshape(C,size(C,1),NrSyst^2,3);
c=permutedims(c,[1,3,2]);
c=reshape(c,WN*3,NrSyst^2)';
c[isnan.(c)].=0;
c .= c .- minimum(c[:]);                                               # normalize to 0 and 1
c .= c ./ maximum(c[:]);

N=reshape(1:NrSyst^2,NrSyst,NrSyst);
k=[N[i, i] for i in 1:NrSyst];
N=(1:NrSyst^2)[setdiff(1:NrSyst^2,k)];
S=collect(permutations(CH,2));
S=[join(string.(x), " ") for x in eachrow(S)];

heatmap!(SG,subplot=1,
         c[N,:],
         c=CG,
         xtickfont=font(12), ytickfont=font(12),
         colorbar = true,
         yticks=(1:42,S),
         xticks=(100," ")
         ) 

c=dropdims(mean(C[50:end,:,:,:],dims=1),dims=1);                       # mean over windows
c[isnan.(c)].=0;
c .= c .- minimum(c[:]);                                               # normalize to 0 and 1
c .= c ./ maximum(c[:]);

CG= cgrad([RGB(0.9,0.9,0.9), RGB(0.3,.3,0.3), :magenta, :cyan],        # define the color map
          [0.0, 0.25, 0.2501, 0.635, 1],scale=:linear);                #   lower quarter in grays

h = [heatmap!(SG,subplot=k+1,
              c[:,:,k],                                                # heatmaps
              c = CG, clim=(0,1),
              colorbar = true,
              title=@sprintf("(%d)",k),
              xtickfont=font(12), ytickfont=font(12),
              xlims=(0.5, 7.5), ylims=(0.5, 7.5),
              aspect_ratio = :equal
            )
     for k in 1:3];

q=reshape(c,NrSyst*NrSyst,3);                                          # normalize each case to 0 and 1
for k=1:3
    q[:,k] .= q[:,k] .- minimum(q[:,k]);
    q[:,k] .= q[:,k] ./ maximum(q[:,k]);
end
q=reshape(c,NrSyst,NrSyst,3);
q[q .< 0.25] .= 0;                                                     # disregard lower quarter in graphs

[graphplot!(SG,subplot=k+4,
            DiGraph(q[:,:,k]),                                         # make digraphs
            markersize = .2,
            names = 1:NrSyst,
            fontsize = 12,
            linecolor = :darkgrey,
            method=:chorddiagram,
            fillcolor=:white
           ) 
     for k in 1:3];

display(SG)

print("Make pdf file and continue? ");
readline()
                                                                       
savefig(@sprintf("PDFs/CD_DDA_C_15dB__WL%d_WS%d_WN%d.pdf",WL,WS,WN))

###

CE=C .* E;                                                               # causality * ergodicity

l=@layout[a{0.5h} ; b c d; e f g];                                     # plot results
SG = plot(layout = l,size=(1500,1500));
CG= cgrad([:white, RGB(1,0.97,0.86), RGB(0.55,0.27,0.07)],
          [0,0.1],scale=:linear);

c=reshape(CE,size(CE,1),NrSyst^2,3);
c=permutedims(c,[1,3,2]);
c=reshape(c,WN*3,NrSyst^2)';
c[isnan.(c)].=0;
c .= c .- minimum(c[:]);                                               # normalize to 0 and 1
c .= c ./ maximum(c[:]);

heatmap!(SG,subplot=1,
         c[N,:],
         c=CG,
         xtickfont=font(12), ytickfont=font(12),
         colorbar = true,
         yticks=(1:42,S),
         xticks=(100," ")
         ) 

c=dropdims(mean(CE[50:end,:,:,:],dims=1),dims=1);                       # mean over windows
c[isnan.(c)].=0;
c .= c .- minimum(c[:]);                                               # normalize to 0 and 1
c .= c ./ maximum(c[:]);

CG= cgrad([RGB(0.9,0.9,0.9), RGB(0.3,.3,0.3), :magenta, :cyan],        # define the color map
          [0.0, 0.25, 0.2501, 0.635, 1],scale=:linear);                #   lower quarter in grays

h = [heatmap!(SG,subplot=k+1,
              c[:,:,k],                                                # heatmaps
              c = CG, clim=(0,1),
              colorbar = true,
              title=@sprintf("(%d)",k),
              xtickfont=font(12), ytickfont=font(12),
              xlims=(0.5, 7.5), ylims=(0.5, 7.5),
              aspect_ratio = :equal
            )
     for k in 1:3];

q=reshape(c,NrSyst*NrSyst,3);                                          # normalize each case to 0 and 1
for k=1:3
    q[:,k] .= q[:,k] .- minimum(q[:,k]);
    q[:,k] .= q[:,k] ./ maximum(q[:,k]);
end
q=reshape(c,NrSyst,NrSyst,3);
q[q .< 0.25] .= 0;                                                     # disregard lower quarter in graphs

[graphplot!(SG,subplot=k+4,
            DiGraph(q[:,:,k]),                                         # make digraphs
            markersize = .2,
            names = 1:NrSyst,
            fontsize = 12,
            linecolor = :darkgrey,
            method=:chorddiagram,
            fillcolor=:white
           ) 
     for k in 1:3];

display(SG)

print("Make pdf file and continue? ");
readline()
                                                                       
savefig(@sprintf("PDFs/CD_DDA_CE_15dB__WL%d_WS%d_WN%d.pdf",WL,WS,WN))

###

C[isnan.(C)] .= 0;
C=C./maximum(C[:]);
C=permutedims(C,[1,4,2,3]);
C=reshape(C,WN*3,NrCH,NrCH);

WLsvd=100; WSsvd=1;
WNsvd=Int(1+floor((size(C,1)-WLsvd)/(WSsvd)));

UU=fill(NaN,WNsvd,NrCH^2); SS=fill(NaN,WNsvd);
for wn=1:WNsvd
    (u,s,v)=svd(reshape(C[(1:WLsvd) .+ (wn-1)*WSsvd,:,:],WLsvd,NrCH^2)');
    UU[wn,:]=u[:,1];
    SS[wn]=s[1,1]; 
end

t=collect(1:WNsvd) .+ WLsvd;

plot(t,SS,label=L"{\cal C}",xticks=(WN:WN:2*WN),legendfontsize=18)

###

CE[isnan.(CE)] .= 0;
CE=CE./maximum(CE[:]);
CE=permutedims(CE,[1,4,2,3]);
CE=reshape(CE,WN*3,NrCH,NrCH);

WLsvd=100; WSsvd=1;
WNsvd=Int(1+floor((size(CE,1)-WLsvd)/(WSsvd)));

UU=fill(NaN,WNsvd,NrCH^2); SS=fill(NaN,WNsvd);
for wn=1:WNsvd
    (u,s,v)=svd(reshape(CE[(1:WLsvd) .+ (wn-1)*WSsvd,:,:],WLsvd,NrCH^2)');
    UU[wn,:]=u[:,1];
    SS[wn]=s[1,1]; 
end

t=collect(1:WNsvd) .+ WLsvd;

plot!(t,SS,label=L"{\cal C} \star {\cal E}",xticks=(WN:WN:2*WN),legendfontsize=18)

print("Make pdf file and continue? ");
readline()
                                                                       
savefig(@sprintf("PDFs/CD_DDA_SVD_15dB__WL%d_WS%d_WN%d.pdf",WL,WS,WN));


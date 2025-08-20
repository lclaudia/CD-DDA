using DataFrames
using Combinatorics
using LinearAlgebra
using Printf
using Random
using JLD2
using Statistics
using DelimitedFiles
using StatsBase
using Plots
using LaTeXStrings
using Graphs
using GraphRecipes
using Colors



if Sys.iswindows()
   SL="\\";
else
   SL="/";
end

function dir_exist(DIR)
    if !isdir(DIR)
        mkdir(DIR)
    end
end

function number_to_string(n::Number)
   return @sprintf("%.15f", n);
end

function integrate_ODE_general_BIG(MOD_nr,MOD_par,dt,L,DIM,ODEorder,X0,FNout,CH_list,DELTA,TRANS=nothing)
  if TRANS===nothing
     TRANS=0;
  end

  if Sys.iswindows()
     if !isfile("i_ODE_general_BIG.exe")
        cp("i_ODE_general_BIG","i_ODE_general_BIG.exe");
     end

     CMD=".\\i_ODE_general_BIG.exe";
  else
     CMD="./i_ODE_general_BIG";
  end

  MOD_NR = join(MOD_nr, " ");
  CMD = "$CMD -MODEL $MOD_NR";
  MOD_PAR = join(MOD_par, " ");
  CMD = "$CMD -PAR $MOD_PAR";
  ANF=join(X0," ");
  CMD = "$CMD -ANF $ANF";
  CMD = "$CMD -dt $dt";
  CMD = "$CMD -L $L";
  CMD = "$CMD -DIM $DIM";
  CMD = "$CMD -order $ODEorder";
  if TRANS>0
     CMD = "$CMD -TRANS $TRANS";
  end
  if length(FNout)>0
     CMD = "$CMD -FILE $FNout";
  end
  CMD = "$CMD -DELTA $DELTA";
  CMD = "$CMD -CH_list $(join(CH_list," "))";

  if length(FNout)>0
     if Sys.iswindows()
        run(Cmd(string.(split(CMD, " "))));
     else
        run(`sh -c $CMD`);
     end
  else
     if Sys.iswindows()
       X = read(Cmd(string.(split(CMD, " "))),String);
     else
       X = read(`sh -c $CMD`,String);
     end
     X = split(strip(X), '\n');
     X = hcat([parse.(Float64, split(row)) for row in X]...)';

     return X
  end
end

function index(DIM, ORDER)
    B = ones(DIM^ORDER, ORDER)
    
    if DIM > 1
        for i = 2:(DIM^ORDER)
            if B[i-1, ORDER] < DIM
                B[i, ORDER] = B[i-1, ORDER] + 1
            end
            
            for i_DIM = 1:ORDER-1
                if round((i/DIM^i_DIM - floor(i/DIM^i_DIM))*DIM^i_DIM) == 1
                    if B[i-DIM^i_DIM, ORDER-i_DIM] < DIM
                        for j = 0:DIM^i_DIM-1
                            B[i+j, ORDER-i_DIM] = B[i+j-DIM^i_DIM, ORDER-i_DIM] + 1
                        end
                    end
                end
            end
        end
        
        i_BB = 1
        BB = Vector{Int}[]
        for i = 1:size(B,1)
            jn = 1
            for j = 2:ORDER
                if B[i, j] >= B[i, j-1]
                    jn += 1
                end
            end
            if jn == ORDER
                push!(BB, B[i, :])
                i_BB += 1
            end
        end
    else
        println("DIM=1!!!")
    end
    
    return hcat(BB...)
end

function monomial_list(nr_delays, order)
    # monomials
    P = index(nr_delays+1, order)'
    P = P .- ones(Int64,size(P))
    
    P = P[2:size(P,1),:];
    
    return P
end

function make_MODEL(SYST)
   order=size(SYST,2);
   nr_delays=2;

   P=monomial_list(nr_delays,order);
   
   MODEL=fill(0,size(SYST,1))';
   for i=1:size(SYST,1)
       II=SYST[i,:]';
   
       MODEL[i] = findall( sum(abs.(repeat(II,size(P,1),1)-P),dims=2)' .== 0 )[1][2];
   end
   L_AF=length(MODEL)+1;

   return MODEL, L_AF, order
end


function make_MOD_nr(SYST,NrSyst)
   DIM=length(unique(SYST[:,1]));
   order=size(SYST,2)-1;

   P=[[0 0]; monomial_list(DIM*NrSyst,order)];
   
   MOD_nr=fill(0,size(SYST,1)*NrSyst,2);
   for n=1:NrSyst
       for i=1:size(SYST,1)
           II=SYST[i,2:end]';
           II[II .> 0] .+= DIM*(n-1);

           Nr=i+size(SYST,1)*(n-1);
           MOD_nr[Nr,2]=findall( sum(abs.(repeat(II,size(P,1),1)-P),dims=2)' .== 0 )[1][2] - 1;
           MOD_nr[Nr,1]=SYST[i,1]+DIM*(n-1);
       end
       #P[MOD_nr[1:size(SYST,1),2].+1,1:2]
   end
   MOD_nr=reshape(MOD_nr',size(SYST,1)*NrSyst*2)';
   
   return MOD_nr,DIM,order,P
end

function make_MOD_nr_Coupling(FromTo,DIM,P)
   order=size(P,2);
   II=fill(0,size(FromTo,1),4);
   for j=1:size(II,1)
       n1=FromTo[j,1]; k1=FromTo[j,2]+1; range1=3:3+order-1;
       n2=FromTo[j,1+range1[end]]; k2=FromTo[j,2+range1[end]]+1; range2=range1 .+ range1[end];

       JJ=FromTo[j,range1]'; JJ[JJ .> 0] .+= DIM * (n1-1); 
       II[j,4] = findall( sum(abs.(repeat(JJ,size(P,1),1)-P),dims=2)' .== 0 )[1][2] - 1;

       JJ=FromTo[j,range2]'; JJ[JJ .> 0] .+= DIM * (n2-1); 
       II[j,2] = findall( sum(abs.(repeat(JJ,size(P,1),1)-P),dims=2)' .== 0 )[1][2] - 1;

       II[j,1] = DIM*n2-(DIM-k2)-1;
       II[j,3] = DIM*n2-(DIM-k1)-1;
   end
   II=reshape(II',length(II))';

   return II
end

function add_noise(s,SNR)
   # check the length of the noise free signal
   N = length(s);

   # n is the noise realization, make it zero mean and unit variance
   n = randn(N);
   n .= (n.-mean(n))./std(n);
   # c is given  from SNR = 10*log10( var(s)/var(c*n) )
   c= sqrt( var(s)*10^-(SNR/10) );

   s_out = (s+c.*n);

   return s_out
end

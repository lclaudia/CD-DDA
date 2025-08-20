include("DDAfunctions.jl");                                            # set of Julia functions
#
WL=2000;WS=500;                                                        # window length and shift
global WN=500;
##WL=4000;WS=1000;WN=2000;

DDA_DIR  = "DDA";  dir_exist(DDA_DIR);                                 # folders
DATA_DIR = "DATA";
FIG_DIR  = "FIG";  dir_exist(FIG_DIR); 

NOISE=["NoNoise";"15dB"];

NrSyst=7; DIM=3;                                                       # 7 3D systems 
NrCH=NrSyst; CH=collect(1:NrCH);                                       # x-components of 7 systems are channels

LIST=collect(combinations(CH,2));                                      # pairwise combinations of channels
LL1=vcat(LIST...)'; LIST=reduce(hcat,LIST)';

nr_delays=2; dm=4;                                                     # DDA parameters
                                                                       # encoding of DDA model
                                                                       # \dot{v} =  
DDAmodel=[[0 0 1];                                                     #    a_1 v_1 +
          [0 0 2];                                                     #    a_2 v_2 +
          [1 1 1]];                                                    #    a_3 v_1^3
(MODEL, L_AF, DDAorder)=make_MODEL(DDAmodel);                          # DDA model encoding for DDA code

TAU=[32 9]; TM=maximum(TAU);                                           # delays

for n_NOISE = 1:length(NOISE)
    noise=NOISE[n_NOISE];

    FN_data=@sprintf("%s%sCD_DDA_data_%s__WL%d_WS%d_WN%d.ascii",
                      DATA_DIR,SL,noise,WL,WS,WN);                     # data file
    FN_DDA=@sprintf("%s%s%s__WL%d_WS%d_WN%d.DDA",
                     DDA_DIR,SL,noise,WL,WS,WN);                       # DDA file

    if !isfile(join([FN_DDA,"_ST"]))
       if Sys.iswindows()
          if !isfile("run_DDA_AsciiEdf.exe")
             cp("run_DDA_AsciiEdf","run_DDA_AsciiEdf.exe");
          end
    
          CMD=".\\run_DDA_AsciiEdf.exe";
       else
          CMD="./run_DDA_AsciiEdf";
       end
       
       CMD = "$CMD -MODEL $(join(MODEL," "))";                         # model
       CMD = "$CMD -TAU $(join(TAU," "))";                             # delays       
       CMD = "$CMD -dm $dm -order $DDAorder -nr_tau $nr_delays";       # DDA parameters
       CMD = "$CMD -DATA_FN $FN_data -OUT_FN $FN_DDA";                 # input and output files
       CMD = "$CMD -WL $WL -WS $WS";                                   # window length and shift
       CMD = "$CMD -SELECT 1 1 1 0";                                   # ST, CT, and CD DDA
       CMD = "$CMD -WL_CT 2 -WS_CT 2";                                 # take pairwise channels for CT and CD DDA
       CMD = "$CMD -CH_list $(join(LL1," "))";                         # list of channel pairs 
            
       if Sys.iswindows()                                              # run ST, CT, and CD DDA
          run(Cmd(string.(split(CMD, " "))));
       else
          run(`sh -c $CMD`);
       end
    
       rm(@sprintf("%s.info",FN_DDA));     
    end
end



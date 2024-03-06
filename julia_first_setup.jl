
import Pkg; Pkg.add("Combinatorics")
import Pkg; Pkg.add("DataFrames")
import Pkg; Pkg.add("LinearAlgebra")
import Pkg; Pkg.add("Printf")
import Pkg; Pkg.add("Random")
import Pkg; Pkg.add("JLD2")
import Pkg; Pkg.add("Statistics")
import Pkg; Pkg.add("DelimitedFiles")
import Pkg; Pkg.add("Plots")
import Pkg; Pkg.add("StatsBase")

import Pkg; Pkg.add("LaTeXStrings")
import Pkg; Pkg.add("Graphs")
import Pkg; Pkg.add("GraphRecipes")
import Pkg; Pkg.add("Colors")

import Pkg; Pkg.add("MAT")



if Sys.islinux()

   run(`cp i_ODE_general_BIG.linux64 i_ODE_general_BIG`);
   run(`chmod +x i_ODE_general_BIG`);

   run(`cp run_DDA_ASCII.linux64 run_DDA_ASCII`);
   run(`chmod +x run_DDA_ASCII`);

end


if Sys.isapple()
   unm=readchomp(`uname -m`);

   if unm == "arm64"
      run(`cp i_ODE_general_BIG.arm64 i_ODE_general_BIG`);
      run(`cp run_DDA_ASCII.arm64 run_DDA_ASCII`);
   end

   if unm == "x86_64"
      run(`cp i_ODE_general_BIG.x86_64 i_ODE_general_BIG`);
      run(`cp run_DDA_ASCII.x86_64 run_DDA_ASCII`);
   end

   run(`chmod +x i_ODE_general_BIG`);
   run(`chmod +x run_DDA_ASCII`);
end



# MIT License

# Copyright (c) 2021 Fraunhofer-Chalmers Centre

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

########################################
# Reads the results from the "result" directory (assuming that the CPLEX results from http://soa.iti.es/files/RCmax.7z has been formatted and placed there also)
# output results the data profiles and the result_table.md
# The data tables are generated from the results in results.7z
########################################
using DelimitedFiles, Plots, Statistics, Latexify, Printf
cd(@__DIR__ )
folders = [name  for  name in readdir("results") if isnothing(findfirst('.',name))]
##
function data_profile(n_task=nothing)
    T_lr = []
    T_gur = []
    for folder in folders
        res_path = "results/"*folder*"/res.txt"
        res,header = readdlm(res_path,',',header=true)
        t = res[:,(header[:].=="time")]
        if !isnothing(n_task)
            t = t[res[:,(header[:].=="n")].<=n_task]
        end
        T_lr = vcat(T_lr,t)

        res_path = "results/"*folder*"/res_gur_cut.txt"
        res,header = readdlm(res_path,',',header=true)
        t = res[:,(header[:].=="time")]
        if !isnothing(n_task)
            t = t[res[:,(header[:].=="n")].<=n_task]
        end
        T_gur = vcat(T_gur,t)
    end
    
    T = [T_lr T_gur]

    τ = 10 .^ range(-2.1, log10(120), length = 100) 
    n, n_solv = size(T)     
    ρ =  [[count(T[:,i] .< t) for t in τ] ./ n for i in 1:n_solv ]
    plotly()
    title_label = "Data profiles"
    if !isnothing(n_task)
        title_label *= ", n ≤ "*string(n_task)
    end
    plot(τ, 100*ρ, ylims=[0,100], title=title_label, label = ["LR" "Gur"], xlabel="CPU time [s]", ylabel ="Solved [%]", xaxis=:log, xticks = 10.0 .^(-2:2),legend=:bottomright)
end
plt = data_profile()
savefig(plt, "data_profile.png")
display(plt)

plt = data_profile(200)
savefig(plt, "data_profile200.png")
display(plt)
##
function getvalue(res,header,instance,label)
    row = findfirst(res[:,1].==instance)
    if isnothing(row)
        return NaN
    end
    return res[row, header[:].==label]    
end
function markbest(z,zmin)
    if z == zmin 
        return "**"*string(z)*"**"
    else 
        return z
    end
end
function calcmean(res, header, m, n, label)
    rows = (res[:,(header[:].=="m")] .== m) .& (res[:,(header[:].=="n")] .==n)
    if !any(rows)
        return 0
    end
    return  mean(res[rows[:],(header[:].==label)])    
end
function n_unsolved(res, header, m, n)
    rows = (res[:,(header[:].=="m")] .== m) .& (res[:,(header[:].=="n")] .==n)
    if !any(rows)
        return 0
    end
    return Int(sum(1 .-res[rows[:],(header[:].=="solved")])  ) 
end

function data_tables()
    # for each size compute a mean time, iter, nodes, solved, obj 
    machines = [10, 20,30,40,50]
    tasks = [100, 200,500,1000]
    outfile = "result_tables.md"
    io = open(outfile,"w")
    write(io, "#Result summary tables\n")
    write(io, "f <- solution not proved to be optimal\n")
    for folder in folders
        res_path = "results/"*folder*"/res.txt"

        out_table =["m" "n" "LR [s]" "LR #f" "LR iter" "LR nodes" "LR z" "LR root gap [%]" "Gur [s]" "Gur #f" "Gur iter" "Gur nodes" "Gur z" "Gur root gap [%]" "Gur gap [%]" "Cplex z"]

        LRres, LRheader = readdlm(res_path,',',header=true)      
        res_path = "results/"*folder*"/res_gur_cut.txt"  
        GURres, GURheader = readdlm(res_path,',',header=true)
        CplexRes, CplexHeader = readdlm("results/"*folder*"/Cplex2horas"*folder*".csv",';',header=true)
        CplexHeader = [CplexHeader "n" "m"]
        CplexRes = Any[CplexRes Int.(floor.(CplexRes[:,1] .-1,sigdigits=1)) Int.(floor.(CplexRes[:,1].-1;digits=-1)) .% 100]

        for m in machines
            for n in tasks
                LRtime = @sprintf("%0.2f", calcmean(LRres, LRheader, m, n, "time"))
                LRunsolved = n_unsolved(LRres, LRheader, m, n)
                LRiter = round(Int,calcmean(LRres, LRheader, m, n, "iter"))
                LRnodes = round(Int,calcmean(LRres, LRheader, m, n, "nodes"))
                LRz = calcmean(LRres, LRheader, m, n, "z")
                LRrootgap = @sprintf("%0.2f", calcmean(LRres, LRheader, m, n, "lbgaproot")*100)
                Gurtime =  @sprintf("%0.2f",calcmean(GURres, GURheader, m, n, "time"))
                Gurunsolved = n_unsolved(GURres, GURheader, m, n)
                Guriter =round(Int, calcmean(GURres, GURheader, m, n, "iter"))
                Gurnodes =round(Int, calcmean(GURres, GURheader, m, n, "nodes"))
                Gurz = round(calcmean(GURres, GURheader, m, n, "z"), digits=1)
                Gurrootgap = @sprintf("%0.2f", calcmean(GURres, GURheader, m, n, "lbrootgap")*100)
                Gurrgap =@sprintf("%0.2f",  calcmean(GURres, GURheader, m, n, "gap")*100  )
                Cplexz = calcmean(CplexRes, CplexHeader, m, n, "Objective")
                out_table = vcat(out_table, [m n LRtime LRunsolved LRiter LRnodes LRz LRrootgap Gurtime Gurunsolved Guriter Gurnodes Gurz Gurrootgap Gurrgap Cplexz])       
            end
        end
        write(io, "## "*folder*"\n\n")
        write(io,  string(md(out_table, latex=false)))
        write(io, "\n")
    end

    out_table =["type" "instance" "LR value" "LR opt" "Gurobi value" "Gurobi opt" "Cplex value" "Cplex opt"]
    for folder in folders
        res_path = "results/"*folder*"/res.txt"
        LRres, LRheader = readdlm(res_path,',',header=true)     
        res_path = "results/"*folder*"/res_gur_cut.txt"  
        GURres, GURheader = readdlm(res_path,',',header=true)
        GURres = round.(Int,GURres)
        CplexRes, CplexHeader = readdlm("results/"*folder*"/Cplex2horas"*folder*".csv",';',header=true)
        CplexRes = round.(Int,CplexRes)
        for instance in LRres[:,1]
            LRz = getvalue(LRres,LRheader,instance,"z")
            LRopt = getvalue(LRres,LRheader,instance,"solved")
            GURz = getvalue(GURres,GURheader,instance,"z")
            Guropt = getvalue(GURres,GURheader,instance,"solved")
            Cplexz = getvalue(CplexRes,CplexHeader,instance,"Objective")
            Cplexopt = getvalue(CplexRes,CplexHeader,instance,"Proven opt.")
            out_table = vcat(out_table,[folder instance LRz LRopt GURz Guropt Cplexz Cplexopt] )
        end
    end
    
    write(io, "## All instances summary \n\n")
    A = copy(out_table[2:end,3:end])
    A[isnan.(A)].=Inf
    z = A[:,1:2:end]
    best_z = Int.(minimum(z, dims=2))
    has_best =  z .<= best_z
    has_strict_best = reduce(hcat,[z[:,j] .< minimum(z[:,setdiff(1:3,j)],dims=2) for j in 1:3])
    has_opt = A[:,2:2:end] .== 1
    unique_opt = [has_opt[i,j] && sum(has_opt[i,:])==1  for i in 1:size(has_opt,1),j in 1:3]

    write(io,  string(md( [["" "LR" "Gurobi" "Cplex"];
     hcat("Best sol [%]", round.(100*sum(has_best,dims=1)/size(has_best,1),digits=2));
     hcat("Strict best sol [%]", round.(100*sum(has_strict_best,dims=1)/size(has_best,1),digits=2));
     hcat("Proved opt [%]",round.(100*sum(has_opt,dims=1)/size(has_opt,1),digits=2));
     hcat("Uniquely proved opt [%]",round.(100*sum(unique_opt,dims=1)/size(has_opt,1),digits=2))], latex=false)))
    
    writedlm("res_all_instances.csv", out_table,",")

    write(io, "## All instances \n\n")
    for (i, j) in Tuple.(findall(has_best))
        out_table[i+1,1+2*j] = "**"*string( out_table[i+1,1+2*j] )*"**"

    end

    write(io,  string(md(out_table, latex=false)))

    close(io)

end
data_tables()
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
# Runs the LR algorithm and Gurobi on all categories of instances in the "instance" folder and writes the result data to the "result" folder
# Results from our run is already in the results.7z folder
########################################
using DelimitedFiles
cd(@__DIR__ )
folders = [name  for  name in readdir("instances") if isnothing(findfirst('.',name))]
instance_names = Dict()
n_tot = 0
T_limit = 300
for folder in folders
    global  n_tot
    instance_name = []
    for file in readdir("instances/"*folder)
        if file[end-3:end] != ".csv"
            continue
        end
        inst_path = folder*"/"*file
        push!(instance_name,file[1:end-4])
        n_tot += 1
    end
    instance_names[folder]=instance_name
end
cp("../binary/RCmax.exe","RCmax.exe")
##
n_solved = 0
start = time()
for folder in folders
    global  n_solved
    out_path = "results/"*folder*"/res.txt"
    if isfile(out_path)
        rm(out_path)
    end

    for name in instance_names[folder]
        inst_path = "instances/"*folder*"/"*name*".csv"
        A = readdlm(inst_path,',')
        (m,n) = size(A)
        println(inst_path)
        if n <= 200 
            run(`RCmax.exe $inst_path $T_limit -o $out_path`)
        else
            run(`RCmax.exe $inst_path $T_limit -o $out_path -nL -nS`)
        end
        n_solved += 1
        println(time()-start,"s ", n_solved,"/",n_tot)
    end
    name_v = vcat("name", instance_names[folder])
    A = readdlm(out_path, ',')
    writedlm(out_path, [name_v A],',')
end

## with Gurobi 
using JuMP, Gurobi
function solveRCmax(inst_path)
    p = readdlm(inst_path,',')
    m,n = size(p)
    model = Model(Gurobi.Optimizer)
    @variable(model,z)
    @objective(model,Min,z)
    @variable(model, x[1:m,1:n], Bin)
    @constraint(model, [i in 1:m], p[i,:]'*x[i,:] <= z)
    @constraint(model, [j in 1:n], sum(x[:,j]) == 1)
    set_time_limit_sec(model,T_limit)
    set_silent(model)
    set_optimizer_attribute(model,"Threads", 1)
    set_optimizer_attribute(model, "Cuts", 2)
    set_optimizer_attribute(model,"MIPGap", 1e-8)
    set_optimizer_attribute(model,"MIPGapAbs", 0.999) # integer costs, ensure that it is not a round off error
    optimize!(model)
    
    zv = objective_value(model)
    zb = objective_bound(model)
    st = solve_time(model)
    si = simplex_iterations(model)
    nc = node_count(model)
    solved = st < T_limit

    unset_binary.(x)
    set_lower_bound.(x,0.0)
    optimize!(model)
    z_lp = objective_value(model)
    
    gap = (zv-zb)/zv 
    lb_root_gap = (zv-z_lp)/zv

    return [m n nc si st solved gap lb_root_gap zv]
end

start = time()
n_solved = 0
for folder in folders
    global  n_solved
    out_path = "results/"*folder*"/res_gur_cut.txt"
    A = ["name" "m" "n" "nodes" "iter" "time" "solved" "gap" "lbrootgap" "z"]
    for name in instance_names[folder]
        inst_path = "instances/"*folder*"/"*name*".csv"
        println(inst_path)
        res = solveRCmax(inst_path)
        n_solved += 1
        println(time()-start,"s ", n_solved,"/",n_tot)
        A = vcat(A, hcat(name, res))
        writedlm(out_path, A, ",")
    end
end
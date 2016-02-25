using FloorLayout, JuMP, CPLEX, Gurobi

filename = split(splitdir(@__FILE__)[2], ".")[1] # pick up current file name, strip off directory and file postfix
const fp = open(joinpath(pwd(),"$(filename).txt"), "w+")
@show dir = splitdir(@__FILE__)[1]
const fp = open(joinpath(dir,"results","$(filename).txt"), "w+")

maxlevels = 5

include("instances-1DFLP.jl")

for bench in [LW11, S8, S8H, S9, S9H, S10, S11, P15, P17, P18, AMI33]

    p = bench[:p]
    if istriu(p)
        p = p + p'
    else
        @assert issym(p)
    end
    l = bench[:l]
    N = length(l)
    L = 2sum(l)

    @assert size(p,1) == size(p,2) == length(l)

    # model = Model(solver=GurobiSolver())
    model = Model(solver=CplexSolver(CPX_PARAM_TILIM=2*60*60))

    @defVar(model, 0 ≤ pos[1:N] ≤ L)
    @defVar(model, obj[i=1:N,j=(i+1):N] ≥ 0)
    @defVar(model, bin[i=1:N,j=1:N;i!=j], Bin)

    @setObjective(model, Min, sum{p[r,s]*obj[r,s], (r,s)=combinations(1:N,2)})

    for i in 1:N, j in (i+1):N
        @addConstraints(model, begin
            pos[i] - pos[j] + l[i] + l[j] ≤ L*(1-bin[i,j])
            pos[j] - pos[i] + l[i] + l[j] ≤ L*(1-bin[j,i])

            bin[i,j] + bin[j,i] == 1

            obj[i,j] ≥  pos[i] - pos[j]
            obj[i,j] ≥ -pos[i] + pos[j]
        end)
    end

    solvetime = @elapsed (stat = solve(model))
    @assert stat == :Optimal

    println(fp, "$(bench[:name]): $(getObjectiveValue(model)) ($solvetime sec)")
    flush(fp)
end

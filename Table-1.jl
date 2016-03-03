using FloorLayout, JuMP, LiftedHierarchies, Gurobi, Mosek

filename = split(splitdir(@__FILE__)[2], ".")[1] # pick up current file name, strip off directory and file postfix
dir = splitdir(@__FILE__)[1]
const fp = open(joinpath(dir,"results","$(filename).txt"), "w+")

println(fp, "# Table 1 Results")
println(fp, "## Tests lifted LP and SDP representations for 1D-FLP. Runtimes can be compared against the level-two combinatorial bound, presented in Table 2.")
println(fp)

include("instances-1DFLP.jl")

const UB = Dict(:LW11  => 13867.0,
                :S8    => 1602.0,
                :S8H   => 4649.0,
                :S9    => 4939.0,
                :S9H   => 9391.0,
                :S10   => 5563.0,
                :S11   => 13867.0,
                :P15   => 12610.0,
                :P17   => 18508.0,
                :P18   => 21355.0,
                :AMI33 => 173628.47645759638)

gap(U,L) = 100*(U-L) / U

for bench in [S8, S8H, S9, S9H, S10, S11, LW11, P15, P17, P18, AMI33]

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

    println(fp, "Bench = $(bench[:name])")
    let LSmodel = Model(solver=GurobiSolver(TimeLimit=10*60))
        @defVar(LSmodel, 0 ≤ pos[1:N] ≤ L)
        @defVar(LSmodel, obj[i=1:N,j=(i+1):N] ≥ 0)
        @defVar(LSmodel, bin[i=1:N,j=1:N;i!=j], Bin)

        @setObjective(LSmodel, Min, sum{p[r,s]*obj[r,s], (r,s)=combinations(1:N,2)})

        for i in 1:N, j in (i+1):N
            @addConstraints(LSmodel, begin
                pos[i] - pos[j] + l[i] + l[j] ≤ L*(1-bin[i,j])
                pos[j] - pos[i] + l[i] + l[j] ≤ L*(1-bin[j,i])

                bin[i,j] + bin[j,i] == 1

                obj[i,j] ≥  pos[i] - pos[j]
                obj[i,j] ≥ -pos[i] + pos[j]
            end)
        end

        # @show LiftedHierarchies.sherali_adams!(LSmodel, 2)
        @show LiftedHierarchies.lovasz_schrijver!(LSmodel)
        solvetime = @elapsed solve(LSmodel, relaxation=true)

        upper = get(UB, bench, NaN)
        lower = getObjectiveValue(LSmodel)
        println(fp, "    LP: $lower lowerbound ($(gap(upper,lower)) gap)")
        println(fp, "        $(solvetime) sec")
        flush(fp)
    end

    let Lasmodel = Model(solver=MosekSolver())
        @defVar(Lasmodel, 0 ≤ pos[1:N] ≤ L)
        @defVar(Lasmodel, obj[i=1:N,j=(i+1):N] ≥ 0)
        @defVar(Lasmodel, bin[i=1:N,j=1:N;i!=j], Bin)

        @setObjective(Lasmodel, Min, sum{p[r,s]*obj[r,s], (r,s)=combinations(1:N,2)})

        for i in 1:N, j in (i+1):N
            @addConstraints(Lasmodel, begin
                pos[i] - pos[j] + l[i] + l[j] ≤ L*(1-bin[i,j])
                pos[j] - pos[i] + l[i] + l[j] ≤ L*(1-bin[j,i])

                bin[i,j] + bin[j,i] == 1

                obj[i,j] ≥  pos[i] - pos[j]
                obj[i,j] ≥ -pos[i] + pos[j]
            end)
        end

        # LiftedHierarchies.lasserre!(Lasmodel, 2)
        LiftedHierarchies.lovasz_schrijver_plus!(Lasmodel)
        buildInternalModel(Lasmodel)
            putdouparam(Lasmodel.internalModel.task, MSK_DPAR_OPTIMIZER_MAX_TIME, 10*60.0)
        solvetime = @elapsed solve(Lasmodel, relaxation=true)

        upper = get(UB, bench, NaN)
        lower = getObjectiveValue(Lasmodel)
        println(fp, "    SDP: $lower lowerbound ($(gap(upper,lower)) gap)")
        println(fp, "         $(solvetime) sec")
        flush(fp)
    end
end

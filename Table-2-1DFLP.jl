using FloorLayout, JuMP, Gurobi

filename = split(splitdir(@__FILE__)[2], ".")[1] # pick up current file name, strip off directory and file postfix
dir = splitdir(@__FILE__)[1]
const fp = open(joinpath(dir,"results","$(filename).txt"), "w+")

println(fp, "# Table 2 1D-FLP Results")
println(fp, "## Compares levels 2-5 of the combinatorial bounding scheme, applied to the 1D-FLP, against Gurobi with the same run-time budget.")
println(fp)

const maxlevels = 5

include("instances-1DFLP.jl")

const UB = Dict(:LW11 => 13867.0,
                :S8 => 1602.0,
                :S8H => 4649.0,
                :S9 => 4939.0,
                :S9H => 9391.0,
                :S10 => 5563.0,
                :S11 => 13867.0,
                :P15 => 12610.0,
                :P17 => 18508.0,
                :P18 => 21355.0,
                :AMI33 => 173628.47645759638)

gap(U,L) = 100*(U-L) / U

LEVELS = [2,2,3,4,5] # first run is just to ensure compilation

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

    println(fp, "Bench = $(bench[:name])")
    for K in LEVELS
        model = Model(solver=GurobiSolver())

        @defVar(model, t[1:N,1:N] >= 0)
        @setObjective(model, Min, sum{p[r,s]*t[r,s], (r,s)=combinations(1:N, 2)})

        println(fp, "    K = $K")
        function solve_multi_box(boxes::Vector{Int})
            if length(boxes) == 2 # short-circuit if we have the closed-form expression
                ii, jj = boxes[1], boxes[2]
                return p[ii,jj] * (l[ii] + l[jj])
            end

            innermod = Model(solver=GurobiSolver(OutputFlag=0))

            p′ = p[boxes,boxes]
            l′ = l[boxes]
            n′ = length(l′)
            L′ = 2sum(l′)
            @defVar(innermod, 0 ≤ x[1:n′] ≤ L′)
            @defVar(innermod, d[i=1:n′,j=(i+1):n′] ≥ 0)
            @defVar(innermod, w[i=1:n′,j=1:n′;i!=j], Bin)

            for i in 1:n′, j in (i+1):n′
                @addConstraints(innermod, begin
                    x[i] + l′[i] ≤ x[j] - l′[j] + L′*(1-w[i,j])
                    x[j] + l′[j] ≤ x[i] - l′[i] + L′*(1-w[j,i])

                    w[i,j] + w[j,i] == 1

                    d[i,j] ≥  x[i] - x[j]
                    d[i,j] ≥ -x[i] + x[j]
                end)
            end
            @setObjective(innermod, Min, sum{p′[r,s]*d[r,s], (r,s)=combinations(1:n′, 2)})

            stat = solve(innermod)
            @assert stat == :Optimal

            return getObjectiveValue(innermod)
        end

        elapse = 0.0
        totalsets = 0
        trimmed = 0
        for kk in 2:K, group in combinations(1:N, kk)
            # solve 3-box MIP
            cutout = false
            if kk > 2
                for x in group
                    cutout |= all(zz -> p[x,zz] <= 0, setdiff(group,x))
                    cutout && break
                end
            end
            totalsets += 1
            cutout && (trimmed += 1)
            cutout && continue
            elapse += @elapsed (rhs = solve_multi_box(group))
            @addConstraint(model, sum{p[r,s]*t[r,s], (r,s)=combinations(group,2)} >= rhs)
        end

        solvetime = @elapsed (stat = solve(model))
        @assert stat == :Optimal

        upper = get(UB, bench[:name], NaN)
        lower = getObjectiveValue(model)
        println(fp, "        MIP: $lower lowerbound ($(gap(upper,lower)) gap)")
        println(fp, "             $totalsets subproblems ($trimmed trimmed)")
        println(fp, "             $elapse elapsed, $solvetime solvetime")
        println(fp)

        # get MIP lowerbound
        bnbsolver = GurobiSolver(TimeLimit=(elapse+solvetime))
        bnbmod = Model(solver=bnbsolver)

        @defVar(bnbmod, 0 ≤ pos[1:N] ≤ L)
        @defVar(bnbmod, obj[i=1:N,j=(i+1):N] ≥ 0)
        @defVar(bnbmod, bin[i=1:N,j=1:N;i!=j], Bin)

        @setObjective(bnbmod, Min, sum{p[r,s]*obj[r,s], (r,s)=combinations(1:N,2)})

        for i in 1:N, j in (i+1):N
            @addConstraints(bnbmod, begin
                pos[i] - pos[j] + l[i] + l[j] ≤ L*(1-bin[i,j])
                pos[j] - pos[i] + l[i] + l[j] ≤ L*(1-bin[j,i])

                bin[i,j] + bin[j,i] == 1

                obj[i,j] ≥  pos[i] - pos[j]
                obj[i,j] ≥ -pos[i] + pos[j]
            end)
        end

        solve(bnbmod)

        lower = MathProgBase.getobjbound(bnbmod.internalModel)
        nodes = MathProgBase.getnodecount(bnbmod.internalModel)
        println(fp, "        B&B: $lower lowerbound ($(gap(upper,lower)) gap)")
        println(fp, "             $(elapse+solvetime) sec")
        println(fp, "             $nodes nodes")
        println(fp)
        flush(fp)
    end
end

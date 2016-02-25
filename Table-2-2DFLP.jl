using FloorLayout, JuMP, CPLEX, Gurobi

filename = split(splitdir(@__FILE__)[2], ".")[1] # pick up current file name, strip off directory and file postfix
dir = splitdir(@__FILE__)[1]
const fp = open(joinpath(dir,"results","$(filename).txt"), "w+")

println(fp, "# Table 2 2D-FLP Results")
println(fp, "## Tests lifted LP and SDP representations for 2D-FLP. Runtimes can be compared against the level-two combinatorial bound, presented in Table 2.")
println(fp)

β = 5
form = Unary(SOC())
maxlevel = 5

const UB = Dict(:hp => 62105.380137346525,
                :apte => 188631.01205865975,
                :xerox => 352437.03500702174,
                :Camp91 => 18522.78606519656,
                :Bozer97_1 => 221.72921422344973,
                :Bazaraa75_1 => 7883.477535194013,
                :Bazaraa75_2 => 6035.22990,
                :Bozer97_2 => 131.82764508093453,
                :Bozer91 => 23090.180383161554,
                :Armour62_1 => 3565.78203,
                :Armour62_2 => 249323.162)

gap(U,L) = 100*(U-L) / U

for bench in [:hp,:apte,:xerox,:Camp91,:Bozer97_1,:Bozer97_2,:Bazaraa75_1,:Bazaraa75_2,:Bozer91,:Armour62_1,:Armour62_2]

    prob = get_problem_data(bench, β)
    n = prob.N
    W = max(prob.H, sum(prob.wub), sum(prob.hub))
    H = max(prob.H, sum(prob.wub), sum(prob.hub))
    wlb, hlb = prob.wlb, prob.hlb
    wub, hub = prob.wub, prob.hub
    c = prob.c

    LEVELS = [2,2,3,4,5] # do two levels of 2 for compilation

    println(fp, "Bench = $bench")
    for K in LEVELS
        model = Model(solver=GurobiSolver())

        indx_set = combinations(1:n,2)

        @defVar(model, tx[1:n,1:n] >= 0)
        @defVar(model, ty[1:n,1:n] >= 0)

        @setObjective(model, :Min, sum{c[r,s]*(tx[r,s]+ty[r,s]), (r,s)=indx_set})

        println(fp, "    K = $K")
        function solve_multi_box(boxes::Vector{Int})
            if length(boxes) == 2 # short-circuit if we have the closed-form expression
                ii, jj = boxes[1], boxes[2]
                return 0.5c[ii,jj] * min(wlb[ii]+wlb[jj], hlb[ii]+hlb[jj])
            end

            innerprob = Problem(length(boxes),
                                W,
                                H,
                                wlb[boxes],
                                hlb[boxes],
                                wub[boxes],
                                hub[boxes],
                                prob.area[boxes],
                                prob.aspect[boxes],
                                c[boxes,boxes])

            innermod = Model(solver=GurobiSolver(OutputFlag=0))
            base_model(innermod, innerprob, form, symbreak=false)

            stat = solve(innermod)

            @assert stat == :Optimal

            return getObjectiveValue(innermod)
        end

        elapse = 0.0
        totalsets = 0
        trimmed = 0
        for kk in 2:K, group in combinations(1:n, kk)
            # solve 3-box MIP
            cutout = false
            if kk > 2
                for x in group
                    cutout |= all(zz -> c[x,zz] <= 0, setdiff(group,x))
                    cutout && break
                end
            end
            totalsets += 1
            cutout && (trimmed += 1)
            cutout && continue
            elapse += @elapsed (rhs = solve_multi_box(group))
            @addConstraint(model, sum{c[r,s]*(tx[r,s]+ty[r,s]), (r,s)=combinations(group,2)} >= rhs)
        end

        solvetime = @elapsed (stat = solve(model))
        @assert stat == :Optimal

        upper = get(UB, bench, NaN)
        lower = getObjectiveValue(model)
        println(fp, "        MIP: $lower lowerbound ($(gap(upper,lower)) gap)")
        println(fp, "             $totalsets subproblems ($trimmed trimmed)")
        println(fp, "             $elapse elapsed, $solvetime solvetime")
        flush(fp)

        # get MIP lowerbound
        bnbsolver = GurobiSolver(TimeLimit=(elapse+solvetime))
        bnbmod = Model(solver=bnbsolver)
        base_model(bnbmod, prob, form, symbreak=false)

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

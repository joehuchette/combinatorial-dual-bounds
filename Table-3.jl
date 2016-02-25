using FloorLayout, JuMP, CPLEX, Mosek, Gurobi

filename = split(splitdir(@__FILE__)[2], ".")[1] # pick up current file name, strip off directory and file postfix
dir = splitdir(@__FILE__)[1]
const fp = open(joinpath(dir,"results","$(filename).txt"), "w+")

println(fp, "# Table 3 Results")
println(fp, "## Tests Takouda SDP formulation for the 2D-FLP. Runtimes can be compared against the level-two combinatorial bound, presented in Table 2.")
println(fp)

const UB = Dict(:hp => 62105.380137346525,
                :apte => 188631.01205865975,
                :xerox => 352437.03500702174,
                :Camp91 => 18522.78606519656,
                :Bozer97_1 => 221.72921422344973,
                :Bazaraa75_1 => 7883.477535194013,
                :Bazaraa75_2 => 13213.552538505586,
                :Bozer97_2 => 131.82764508093453,
                :Bozer91 => 23090.180383161554,
                :Armour62_1 => 22679.140100826913,
                :Armour62_2 => 1.8652032684550043e6)

gap(U,L) = 100*(U-L) / U

const β = 5

for bench in [:hp,:apte,:xerox,:Camp91,:Bozer97_1,:Bozer97_2,:Bazaraa75_1,:Bazaraa75_2,:Bozer91,:Armour62_1,:Armour62_2]
    prob = get_problem_data(bench, β)
    prob.W = max(sum(prob.wub), sum(prob.hub))
    prob.H = max(sum(prob.wub), sum(prob.hub))

    elapse = @elapsed (obj = solveTakouda(prob))

    println(fp, "Bench = $bench")
    println(fp, "    $obj objective")
    println(fp, "    $elapse sec")
    flush(fp)
end

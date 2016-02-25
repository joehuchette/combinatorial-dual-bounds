## Computational materials from "Beating the SDP bound for the floor layout problem: A simple combinatorial idea"

The computational trials are conducted using [JuMP](), an algebraic modeling language embedded in Julia. To run the trials, download a copy of [Julia v0.4](). You can download the required packages with the following Julia script:
```jl
Pkg.add("JuMP")
Pkg.checkout("JuMP", "release-0.11")
Pkg.add("Gurobi")
Pkg.add("Mosek")
Pkg.clone("https://github.com/joehuchette/FloorLayout.jl.git")
Pkg.clone("https://github.com/joehuchette/LiftedHierarchies.jl.git")
```

Run ``compute-optimal-cost-1DFLP.jl`` to compute optimal costs for each of the 1D-FLP instances. The problem data for the 1D-FLP is available in ``instances-1DFLP.jl``; the data for the 2D-FLP instances are available in the ``FloorLayout.jl`` package.

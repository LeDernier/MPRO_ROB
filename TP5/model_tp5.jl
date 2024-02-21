"""
Problème de planification culturale durable

Par :
Leder Yhivert AGUIRRE RAMIREZ
Christian EL MAALOULY
"""


using CPLEX
using JuMP

# model and settings
model = Model(CPLEX.Optimizer)
set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0) # Remove the solver output

#parameters
l_max = 2
a_max = 2
t_max = 10
p_max = 40

E = [((2, 0, 0), (2, 0, 0)), ((2, 0, 0), (2, 1, 1)), ((2, 1, 1), (2, 2, 1)), ((2, 1, 1), (1, 0, 0)), ((2, 2, 1), (1, 0, 0)), ((1, 0, 0), (2, 0, 0)), ((1, 0, 0), (1, 1, 1)), ((1, 1, 1), (1, 0, 0)), ((1, 1, 1), (1, 2, 1)), ((1, 2, 1), (1, 0, 0))]
C = [1, 2]

L = collect(1:l_max)
A = collect(0:a_max)
T = collect(1:t_max)
P = collect(1:p_max)
R = Dict((1, 1, (0, 1)) => 72,
    (2, 1, (0, 1)) => 120,
    (1, 1, (0, 2)) => 54,
    (2, 1, (0, 2)) => 90,
    (1, 2, (2, 1)) => 54,
    (2, 2, (2, 1)) => 90,
    (1, 2, (1, 2)) => 39,
    (2, 2, (1, 2)) => 65)            # rendement
D = Dict((1, 1) => 1200,
    (1, 2) => 0,
    (1, 3) => 1200,
    (1, 4) => 0,
    (1, 5) => 1200,
    (1, 6) => 0,
    (1, 7) => 1200,
    (1, 8) => 0,
    (1, 9) => 1200,
    (1, 10) => 0,
    (2, 1) => 0,
    (2, 2) => 400,
    (2, 3) => 0,
    (2, 4) => 400,
    (2, 5) => 0,
    (2, 6) => 400,
    (2, 7) => 0,
    (2, 8) => 400,
    (2, 9) => 0,
    (2, 10) => 400)

# variables
@variable(model, x[p in P, e in E, t in T], Bin)
@variable(model, z[p in P], Bin)

# objective
@objective(model, Min, sum(z[p] for p in P))

#constraints
for t in T
    if t % 2 == 0
        @constraint(model, sum(x[p, e, t] * R[(e[2][1], e[2][2], (e[1][3] * 1, 2))] for p in P for e in E if e[2][2] != 0) >= D[(2, t)]) #demande
    else
        @constraint(model, sum(x[p, e, t] * R[e[2][1], e[2][2], (e[1][3] * 2, 1)] for p in P for e in E if e[2][2] != 0) >= D[(1, t)]) #demande
    end
end

@constraint(model, [p in P], x[p, ((l_max, 0, 0), (l_max, 1, 1)), 1] + x[p, ((l_max, 0, 0), (l_max, 0, 0)), 1] == z[p]) # flow conservation constraint in the node (j,1) / culture j, temps 0->1, assuming the initial state is (l,a,j)=(2,0,0)
@constraint(model, [p in P, t in collect(2:t_max)], sum(x[p, e, t] for e in E) == z[p])  # un seul arc/une seule arête par parcelle et temps 

for e in E
    @constraint(model, [p in P, t in collect(1:t_max-1)], x[p, e, t] <= sum(x[p, new_e, t+1] for new_e in E if new_e[1] == e[2])) # flow conservation constraint
end

optimize!(model)

if termination_status(model) == MOI.OPTIMAL
    x_val = value.(x)
    println("objective_value: ", objective_value(model))

    # Verifier les chemins
    for p in P
        println("p:", p)
        tmp_arc = []
        for t in T
            println("t:", t)
            for e in E
                if x_val[p, e, t] > 0
                    if tmp_arc != [] && e[1] != tmp_arc
                        println(" erreur:", e[1], " != ", tmp_arc)
                    end
                    println(" e:", e)
                    tmp_arc = e[2]
                end
            end
        end
        println()
    end

    # Verifier les demandes
    for t in T
        if t % 2 == 0
            val = sum(x_val[p, e, t] * R[(e[2][1], e[2][2], (e[1][3] * 1, 2))] for p in P for e in E if e[2][2] != 0)
            println("t:", t, " D[2, ", t, "]:", D[(2, t)], " val:", val)
        else
            val = sum(x_val[p, e, t] * R[e[2][1], e[2][2], (e[1][3] * 2, 1)] for p in P for e in E if e[2][2] != 0)
            println("t:", t, " D[1, ", t, "]:", D[(1, t)], " val:", val)
        end
    end
    #println("x_val:", x_val)

else
    println("Model infeasible.")
end
"""
RODD - TD7
Problème de planification culturale durable

Elève: Leder Yhivert AGUIRRE RAMIREZ
"""


using CPLEX
using JuMP

# model and settings
model = Model(CPLEX.Optimizer)
#set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)

#parameters
l_max = 2
a_max = 2
t_max = 2
p_max = 4
C = ("ri","ha")
V = ("ri","ha", 0)
E = [("ri","ha"), ("ha","ri"),("ri",0),("ha",0),(0,0),(0,"ha"), (0,"ri")]

L = collect(1:l_max)
A = collect(0:a_max)
T = collect(1:t_max)
P = collect(1:p_max)
R = Dict((1,1,(0,"ri")) => 72,
        (2,1,(0,"ri")) => 120,
        (1,1,(0,"ha")) => 54,
        (2,1,(0,"ha")) => 90,
        (1,2,("ha","ri")) => 54,
        (2,2,("ha","ri")) => 90,
        (1,2,("ri","ha")) => 39,
        (2,2,("ri","ha")) => 65)            # rendement
D = Dict(("ri",1) => 1200,
        ("ri",2) => 0,
        ("ri",3) => 1200,
        ("ri",4) => 0,
        ("ri",5) => 1200,
        ("ri",6) => 0,
        ("ri",7) => 1200,
        ("ri",8) => 0,
        ("ri",9) => 1200,
        ("ri",10) => 0,
        ("ha",1) => 0,
        ("ha",2) => 400,
        ("ha",3) => 0,
        ("ha",4) => 400,
        ("ha",5) => 0,
        ("ha",6) => 400,
        ("ha",7) => 0,
        ("ha",8) => 400,
        ("ha",9) => 0,
        ("ha",10) => 400)



# variables
@variable(model, x[p in P, l in L, a in A, e in E, t in T], Bin)
@variable(model, z[p in P], Bin)

# objective
@objective(model, Min, sum(z[p] for p in P))

#constraints
@constraint(model, [j in C, t in T], sum(x[p,l,a,(i,j),t]*R[(l,a,(i,j))] for p in P for l in L for a in A for i in V if (l,a,(i,j)) in keys(R)) >= D[(j,t)]/20) #demande
@constraint(model, [p in P], sum(x[p,l_max,1,(0,k),1] for k in C) == z[p]) # flow conservation constraint in the node (j,1) / culture j, temps 0->1, assuming the initial state is (l,a,j)=(2,0,0)
@constraint(model, [p in P, t in collect(1:t_max)], sum(x[p,l,a,e,t] for l in L for a in A for e in E) == z[p])  # un seul arc/une seule arête par parcelle et temps 
#@constraint(model, [p in P], sum(x[p,l,a,e,t] for l in L for a in A for e in E for t in T) <= z[p]*l_max*a_max*t_max*length(V)^2)   # relation entre les variables x et z

# flow-conservation constraints
#= @constraint(model, [p in P, l in L, a in collect(0:a_max-1), j in C, t in collect(1:t_max-1)], sum(x[p,l,a,(i,j),t] for i in V if (i,j) in E) == sum(x[p,l,a+1,(j,k),t+1] for k in C if (j,k) in E && a < a_max) + x[p,1,0,(j,0),t+1]) # flow conservation constraint in the node (j,t) / culture j, temps t->t+1
##@constraint(model, [p in P, j in C, t in collect(1:t_max-1)], sum(x[p,l,a_max,(i,j),t] for l in L for i in V if (i,j) in E) == x[p,1,0,(j,0),t+1]) # flow conservation constraint in the node (j,t) / culture j, temps t->t+1; a=a_max
@constraint(model, [p in P, l in collect(1:l_max-1), t in collect(1:t_max-1)], sum(x[p,l,0,(i,0),t] for i in V) == sum(x[p,l+1,0,(0,0),t+1]) + sum(x[p,l,1,(0,k),t+1] for k in C)) # flow conservation constraint in the node (j,t) / jachère j, temps t->t+1
 =#

@constraint(model, [p in P, l in collect(1:l_max-1), t in collect(1:t_max-1)], sum(x[p,l,0,(k,0),t] for k in V) <= x[p,l+1,0,(0,0),t+1])                    # arête (i,j) (=0,=0)
@constraint(model, [p in P, l in L, j in C, t in collect(1:t_max-1)], sum(x[p,l,0,(k,0),t] for k in V) <= x[p,l,1,(0,j),t+1])                                 # arête (i,j) (=0,!=0)

@constraint(model, [p in P, i in C, t in collect(1:t_max-1)], sum(x[p,l,a,(k,i),t] for k in V if (k,i) in E for l in L for a in collect(1:a_max)) <= x[p,1,0,(i,0),t+1])    # arête (i,j) (!=0,=0)

for i in C
    for j in C
        if (i,j) in E
            @constraint(model, [p in P, l in L, a in collect(1:a_max), t in collect(1:t_max-1)], sum(x[p,l,a-1,(k,i),t] for k in V if (k,i) in E) <= x[p,l,a,(i,j),t+1])                      # arête (i,j) (!=0,!=0)
        end
    end
end

println("MODEL:\n",model)
optimize!(model)
println("\nSOLVED.")

isModelFeasible = primal_status(model) == MOI.FEASIBLE_POINT
obj_val = undef
if isModelFeasible
    obj_val = objective_value(model)
    x_val = value.(model[:x])
    println("objective_value: $(obj_val)")
    println("x_val:")

    x_val_corr = [] #(p,t,l,a,(i,j))
    for key in keys(x_val)
        if x_val[key] > 0
            push!(x_val_corr, ((key[1],key[5]),key[2],key[3],key[4]))
        end
    end
    sort!(x_val_corr, by=x->x[1])
    for l in x_val_corr
        println("((p,t),l,a,(i,j),): ", l)
    end

    for j in C
        for t in T
            r_jt = sum(x_val[p,l,a,(i,j),t]*R[(l,a,(i,j))] for p in P for l in L for a in A for i in V if (l,a,(i,j)) in keys(R))
            println("r_{$(j),$(t)}: ", r_jt)
        end
    end
else
    println("Model infeasible.")
end
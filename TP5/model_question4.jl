"""
RODD - TD7
Problème de planification culturale durable

Elève: Leder Yhivert AGUIRRE RAMIREZ
"""


using CPLEX
using JuMP

# model and settings
model = Model(CPLEX.Optimizer)

# Désactive les sorties de CPLEX (optionnel)
#set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)

# Désactive le presolve (simplification automatique du modèle)
#set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", 0)

# Désactive la génération de coupes automatiques
set_optimizer_attribute(model, "CPXPARAM_MIP_Limits_CutsFactor", 0)

# Désactive la recherche dynamique avant de passer au Branch and Bound
set_optimizer_attribute(model, "CPXPARAM_MIP_Strategy_Search", 1)

#parameters
l_max = 2
a_max = 2
t_max = 10
p_max = 40
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


V_bar2 = []         # nodes of the auxiliary graph, of the form: (l,a,i), for i in V
E_bar2 = []         # edges of the auxiliary graph, of the form: ((l1,a1,i),(l2,a2,j)), for (i,j) in E

for e in E
    if e[1] == 0
        if e[2] == 0
            for l in collect(1:l_max-1)
                u = (l,0,e[1])
                v = (l+1,0,e[2])

                # add the new nodes
                if !(u in V_bar2)
                    push!(V_bar2,u)
                end
                if !(v in V_bar2)
                    push!(V_bar2,v)
                end

                # add the new edges
                push!(E_bar2,(u,v))
            end
        else
            for l in collect(1:l_max)
                u = (l,0,e[1])
                v = (l,1,e[2])

                # add the new nodes
                if !(u in V_bar2)
                    push!(V_bar2,u)
                end
                if !(v in V_bar2)
                    push!(V_bar2,v)
                end

                # add the new edges
                push!(E_bar2,(u,v))
            end
        end
    else
        if e[2] == 0
            for l in collect(1:l_max)
                for a in collect(1:a_max)
                    u = (l,a,e[1])
                    v = (1,0,e[2])

                    # add the new nodes
                    if !(u in V_bar2)
                        push!(V_bar2,u)
                    end
                    if !(v in V_bar2)
                        push!(V_bar2,v)
                    end

                    # add the new edges
                    push!(E_bar2,(u,v))
                end
            end
        else
            for l in collect(1:l_max)
                for a in collect(1:a_max-1)
                    u = (l,a,e[1])
                    v = (l,a+1,e[2])

                    # add the new nodes
                    if !(u in V_bar2)
                        push!(V_bar2,u)
                    end
                    if !(v in V_bar2)
                        push!(V_bar2,v)
                    end

                    # add the new edges
                    push!(E_bar2,(u,v))
                end
            end
        end
    end
end
    

# variables
@variable(model, x[p in P, e_bar2 in E_bar2, t in T], Bin)
@variable(model, z[p in P], Bin)

# objective
@objective(model, Min, sum(z[p] for p in P))

#constraints
@constraint(model, [j in C, t in T], sum(x[p,((l1,a1,i),(l2,a2,j)),t]*R[(l2,a2,(i,j))] for p in P for (l1,a1,i) in V_bar2 for l2 in L for a2 in A if ((l1,a1,i),(l2,a2,j)) in E_bar2) >= D[(j,t)]) # demand satisfaction
@constraint(model, [p in P], sum(x[p,((l_max,0,0),(l_max,1,k)),1] for k in C) == z[p]) # flow conservation constraint in the node (j,1) / crop j, time 0->1, assuming the initial state is (l,a,j)=(2,0,0)
@constraint(model, [p in P, t in collect(1:t_max)], sum(x[p,(e1,e2),t] for (e1,e2) in E_bar2) == z[p])  # only one edge by parcel and time 
@constraint(model, [p in P, t in collect(1:t_max-1), e in E_bar2], x[p,e,t] <= sum(x[p,(e[2],v),t+1] for v in V_bar2 if (e[2],v) in E_bar2))     # respecting precedence: x[p,(u,v),t] = 1 implies sum(x[p,(v,w),t+1] for w in V_bar2)

#println("MODEL:\n",model)
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
            push!(x_val_corr, ((key[1],key[3]),key[2]))
        end
    end
    sort!(x_val_corr, by=x->x[1])
    for val in x_val_corr
        println("((p,t),((l1,a1,i),(l2,a2,j))): ", val)
    end

    for j in C
        for t in T
            r_jt = sum(x_val[p,((l1,a1,i),(l2,a2,j)),t]*R[(l2,a2,(i,j))] for p in P for (l1,a1,i) in V_bar2 for l2 in L for a2 in A if ((l1,a1,i),(l2,a2,j)) in E_bar2)
            println("r_{$(j),$(t)}: ", r_jt)
        end
    end

    println("solution_time: ", solve_time(model))
    println("node_count: ", node_count(model))
else
    println("Model infeasible.")
end
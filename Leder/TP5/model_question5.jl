"""
RODD - TD7
Problème de planification culturale durable

Elève: Leder Yhivert AGUIRRE RAMIREZ
"""


using CPLEX
using JuMP

# model and settings
model = Model(CPLEX.Optimizer)
set_optimizer_attribute(model, "CPX_PARAM_SCRIND", 0)

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

#println("E_bar2:\n",E_bar2)

# N[u] := set of outgoing nodes from u=(l,a,j) for j in V
N = Dict();

for (u,v) in E_bar2
    if !(u in keys(N))
        N[u] = Vector{Tuple{Int64, Int64, Any}}()
        push!(N[u],v)
    else
        push!(N[u],v)
    end
end

println("N:\n",N)

# Rots[i] := i-th crop rotation of the form [u_1,...,u_{t_max}] where u_t = (l_t,a_t,j_t), j_t in V

Rots = []   # rotations
t=0
u = (2,0,0)

"""
    r: rotation of crops from 0 to t
    t: time step
"""
function dfs_recursion(rot=Vector{Tuple{Int64, Int64, Any}}([(2,0,0)]),t=0)
    # add the complete rotation
    if t == t_max
        push!(Rots, rot)
        return
    end

    # recursion
    t += 1
    u = last(rot)
    for v in N[u]
        rot_copy = deepcopy(rot)
        push!(rot_copy,v)
        dfs_recursion(rot_copy, t)
    end
end

# by calling dfs_recursion(), Rots is initialized: a list of rotations that start at (2,0,0)
dfs_recursion()
Rots_idx = collect(1:length(Rots))  # index set to bind the rotation i to the variable x_i
#= println("ROTATIONS:")
for r in Rots_idx
    println("$(r): ", Rots[r])
end =#
println("Number of crop rotations: ", length(Rots))

### REFORMULATION ###

# variables
@variable(model, x[r in Rots_idx] >= 0, Int)

# objective
@objective(model, Min, sum(x[r] for r in Rots_idx))

#constraints
# for the r-th rotation we have:
# Rots[Rots_idx[r]][t+1][1] = l_t
# Rots[Rots_idx[r]][t+1][2] = a_t
# Rots[Rots_idx[r]][t+1][3] = j_t
# Rots[Rots_idx[r]][t][3] = a_t
# ...so the demand satisfaction constraint becomes:
@constraint(model, [j in C, t in T], sum(x[r]*R[(Rots[Rots_idx[r]][t+1][1],Rots[Rots_idx[r]][t+1][2],(Rots[Rots_idx[r]][t][3],Rots[Rots_idx[r]][t+1][3]))] for r in Rots_idx if Rots[Rots_idx[r]][t+1][3] == j && (Rots[Rots_idx[r]][t+1][1],Rots[Rots_idx[r]][t+1][2],(Rots[Rots_idx[r]][t][3],Rots[Rots_idx[r]][t+1][3])) in keys(R)) >= D[(j,t)]) # demand satisfaction


#= @constraint(model, [p in P], sum(x[p,((l_max,0,0),(l_max,1,k)),1] for k in C) == z[p]) # flow conservation constraint in the node (j,1) / crop j, time 0->1, assuming the initial state is (l,a,j)=(2,0,0)
@constraint(model, [p in P, t in collect(1:t_max)], sum(x[p,(e1,e2),t] for (e1,e2) in E_bar2) == z[p])  # only one edge by parcel and time 
@constraint(model, [p in P, t in collect(1:t_max-1), e in E_bar2], x[p,e,t] <= sum(x[p,(e[2],v),t+1] for v in V_bar2 if (e[2],v) in E_bar2))     # respecting precedence: x[p,(u,v),t] = 1 implies sum(x[p,(v,w),t+1] for w in V_bar2)
 =#

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

    for r in keys(x_val)
        if x_val[r] > 0
            println("x[$(Rots[r[1]])]=", x_val[r])
        end
    end

    for j in C
        for t in T
            r_jt = sum(x_val[r]*R[(Rots[Rots_idx[r]][t+1][1],Rots[Rots_idx[r]][t+1][2],(Rots[Rots_idx[r]][t][3],Rots[Rots_idx[r]][t+1][3]))] for r in Rots_idx if Rots[Rots_idx[r]][t+1][3] == j && (Rots[Rots_idx[r]][t+1][1],Rots[Rots_idx[r]][t+1][2],(Rots[Rots_idx[r]][t][3],Rots[Rots_idx[r]][t+1][3])) in keys(R))
            println("r_{$(j),$(t)}: ", r_jt)
        end
    end
else
    println("Model infeasible.")
end
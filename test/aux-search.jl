
using Revise

using Tiger
import Tiger.build_base_model

using JuMP, Gurobi
using SparseArrays

import Base.convert


function get_exchange_rxns(cobra::Tiger.CobraModel; objective=true)
    ex_idx = (sum(abs.(cobra.S[:,1:length(cobra.rxns)]),dims=1)  .== 1)[:]
    if !objective
        ex_idx .&= abs.(cobra.c) .== 0.0
    end
    cobra.rxns[ex_idx]
end

indicator_name(varname) = "I__" * varname
function add_indicator(model::Model, var::VariableRef)
    ind = @variable(model, base_name=indicator_name(name(var)), binary=true)
    @constraint(model, ind => {var == 0.0})
    return ind
end

function add_indicator(model::Model, var::String)
    add_indicator(model, variable_by_name(model, var))
end

function create_primal_dual(cobra::Tiger.CobraModel; bound_vars::Vector{String}=[])
    c = cobra.c
    lb = cobra.lb
    ub = cobra.ub
    n = size(cobra.S,2)

    ieq = cobra.csense .== '='
    Aeq = cobra.S[ieq,:]
    beq = cobra.b[ieq]
    
    leq = cobra.csense .== '<'
    Ale = cobra.S[leq,:]
    ble = cobra.b[leq]

    geq = cobra.csense .== '>'
    Age = cobra.S[geq,:]
    bge = cobra.b[geq]

    A = [Ale; -Age]
    b = [ble; -bge]

    model = build_base_model(cobra)
    x = variable_by_name.(model, cobra.vars)

    meq = size(Aeq,1)
    mle = size(A,1)
    @variable(model, lambda[1:meq])
    @variable(model, 0 <= mu[1:mle])
    @variable(model, 0 <= muU[1:n])
    @variable(model, 0 <= muL[1:n])

    # z = muU .* ub
    @variable(model, 0 <= z[1:n])

    @constraint(model, Aeq'*lambda + A'*mu - muL + muU .== c)

    # strict duality
    @constraint(model, c[1:n]'*x[1:n] == lambda'*beq + mu'*b - muL'*lb + z'*ub)

    # bounding constraints
    nbound = length(bound_vars)
    @variable(model, Inds[1:nbound], binary=true)
    set_name.(Inds, "I__" .* bound_vars)

    for i = 1:n
        if name(x[i]) in bound_vars
            Ind = variable_by_name(model, "I__" * name(x[i]))
            @constraint(model, z[i] == muU[i])
            @constraint(model, !Ind => {muU[i] == 0.0})
            @constraint(model, !Ind => {z[i] == 0.0})
            @constraint(model, !Ind => {x[i] <= 0.0})
        else
            @constraint(model, z[i] == muU[i])
        end
    end

    @objective(model, Max, 0)

    return model, Inds
end

struct StandardLP
    # min sense * c'x
    # s.t.  Ax <= b
    #       Aeq x == beq
    #       lb <= x <= ub
    sense::Integer

    c::Vector{Real}

    A::SparseMatrixCSC{Float64,Int}
    b::Vector{Real}

    Aeq::SparseMatrixCSC{Float64,Int}
    beq::Vector{Real}

    lb::Vector{Real}
    ub::Vector{Real}

    vars::Vector{String}
end

function build_base_model(lp::StandardLP; optimizer=Gurobi.Optimizer, model=nothing)
    if isnothing(model)
        model = Model(optimizer)
    end

    n = length(lp.vars)
    x = @variable(model, [1:n])
    set_lower_bound.(x, lp.lb)
    set_upper_bound.(x, lp.ub)

    @constraint(model, lp.A*x .<= lp.b)
    @constraint(model, lp.Aeq*x .== lp.beq)

    if lp.sense == 1
        @objective(model, Min, lp.c'*x)
    else
        @objective(model, Max, lp.c'*x)
    end

    return model
end


function StandardLP(cobra::Tiger.CobraModel)
    cobra = deepcopy(cobra)
    sense = -1  # Cobra models maximize by default
    c = cobra.c
    lb = cobra.lb
    ub = cobra.ub
    n = size(cobra.S,2)

    ieq = cobra.csense .== '='
    Aeq = cobra.S[ieq,:]
    beq = cobra.b[ieq]
    
    leq = cobra.csense .== '<'
    Ale = cobra.S[leq,:]
    ble = cobra.b[leq]

    geq = cobra.csense .== '>'
    Age = cobra.S[geq,:]
    bge = cobra.b[geq]

    A = [Ale; -Age]
    b = [ble; -bge]

    StandardLP(
        sense,
        c,
        A,
        b,
        Aeq,
        beq,
        lb,
        ub,
        cobra.vars
    )
end

function find_conditional_media(cobra::Tiger.CobraModel, ko; min_wt_growth=1.0, max_ko_growth=0.0)
    gene_names = cobra.genes

    lp = StandardLP(cobra)
    n = length(lp.vars)

    model = Model(Gurobi.Optimizer)

    # add WT condition
    v_wt = @variable(model, [1:n])
    set_name.(v_wt, lp.vars .* "_WT")
    set_lower_bound.(v_wt, lp.lb)
    set_upper_bound.(v_wt, lp.ub)
    @constraint(model, lp.A * v_wt .<= lp.b)
    @constraint(model, lp.Aeq * v_wt .== lp.beq)

    # add KO condition with primal-dual expansion
    v_ko = @variable(model, [1:n])
    set_name.(v_ko, lp.vars .* "_KO")
    set_lower_bound.(v_ko, lp.lb)
    set_upper_bound.(v_ko, lp.ub)
    @constraint(model, lp.A * v_ko .<= lp.b)
    @constraint(model, lp.Aeq * v_ko .== lp.beq)

    @variable(model, lambda[1:size(lp.Aeq,1)])
    @variable(model, mu[1:size(lp.A,1)] >= 0)
    @variable(model, muU[1:n] >= 0)
    @variable(model, muL[1:n] >= 0)

    @constraint(model, lp.Aeq'*lambda + lp.A'*mu - muL + muU .== lp.sense*lp.c)

    ex_names = get_exchange_rxns(cobra; objective=false)
    nex = length(ex_names)
    @variable(model, I[1:nex], binary=true)
    for (i, ex) in enumerate(ex_names)
        # bind all variable except the KO variable
        if ex != ko
            j = findfirst(ex .== lp.vars)
            @constraint(model, !I[i] => {v_wt[j] == 0.0})
            @constraint(model, !I[i] => {v_ko[j] == 0.0})

            # transfer switching constraint to the dual except for
            # any reaction in the primal's objective function
            if lp.c[j] == 0
                @constraint(model, I[i] => {muU[j] == 0.0})
                @constraint(model, I[i] => {muL[j] == 0.0})
            end
        end
    end

    # strong duality
    i = [x in ex_names for x in lp.vars]
    @constraint(model, lambda'*lp.beq + mu'*lp.b - muL[.!i]'*lp.lb[.!i] + muU[.!i]'*lp.ub[.!i] - muL[i]'*(lp.lb[i] .* I) + muU[i]'*(lp.ub[i] .* I) == lp.sense*lp.c'*v_ko)

    # knockout KO variable
    set_upper_bound(variable_by_name(model, ko * "_KO"), 0.0)
    set_lower_bound(variable_by_name(model, ko * "_KO"), 0.0)

    # constraint WT to grow
    @constraint(model, lp.sense*lp.c' * v_wt >= min_wt_growth)

    # constraint KO to not grow
    @constraint(model, lp.sense*lp.c' * v_ko <= max_ko_growth)

    @objective(model, Min, sum(I))

    return model
end

function find_conditional_media_ex(cobra::Tiger.CobraModel, ko; min_wt_growth=1.0, max_ko_growth=0.0)
    gene_names = cobra.genes

    lp = StandardLP(cobra)
    n = length(lp.vars)

    model = Model(Gurobi.Optimizer)

    # add WT condition
    v_wt = @variable(model, [1:n])
    set_name.(v_wt, lp.vars .* "_WT")
    set_lower_bound.(v_wt, lp.lb)
    set_upper_bound.(v_wt, lp.ub)
    #@constraint(model, lp.A * v_wt .<= lp.b)
    @constraint(model, lp.Aeq * v_wt .== lp.beq)

    # add KO condition with primal-dual expansion
    ko_idx = findfirst(x->x==ko, lp.vars)
    ko_lb = copy(lp.lb)
    ko_ub = copy(lp.ub)
    ko_lb[ko_idx] = 0.0
    ko_ub[ko_idx] = 0.0
    v_ko = @variable(model, [1:n])
    set_name.(v_ko, lp.vars .* "_KO")
    set_lower_bound.(v_ko, ko_lb)
    set_upper_bound.(v_ko, ko_ub)
    #@constraint(model, lp.A * v_ko .<= lp.b)
    @constraint(model, lp.Aeq * v_ko .== lp.beq)

    @variable(model, lambda[1:size(lp.Aeq,1)])
    #@variable(model, mu[1:size(lp.A,1)] >= 0)
    @variable(model, muU[1:n] >= 0)
    @variable(model, muL[1:n] >= 0)

    @variable(model, I[1:n], binary=true)
    ex_names = get_exchange_rxns(cobra, objective=false)
    Ie = [x in ex_names for x in lp.vars]
    fix.(I[.!Ie], 1)

    #@constraint(model, lp.Aeq'*lambda + lp.A'*mu - muL + muU .== lp.sense*lp.c)
    @constraint(model, lp.Aeq'*lambda - muL + muU .== -lp.sense*lp.c)

    @constraint(model, v_wt[Ie] .<= lp.ub[Ie].*I[Ie])
    @constraint(model, v_ko[Ie] .<= ko_ub[Ie].*I[Ie])
    #@constraint(model, muI[Ie] .<= M*(1 .- I[Ie]))
    # for (i, ex) in enumerate(ex_names)
    #     # bind all variable except the KO variable
    #     if ex != ko
    #         j = findfirst(ex .== lp.vars)
    #         @constraint(model, !I[j] => {v_wt[j] == 0.0})
    #         @constraint(model, !I[j] => {v_ko[j] == 0.0})
    #     end

    #     # transfer switching constraint to the dual except for
    #     # any reaction in the primal's objective function
    #     if lp.c[j] == 0
    #         if lp.ub[i] > 0
    #             @constraint(model, I[j] => {muU[j] == 0.0})
    #         end
    #         if lp.lb[i] < 0
    #             @constraint(model, I[j] => {muL[j] == 0.0})
    #         end
    #     end
    # end

    # strong duality
    
    #@constraint(model, lambda'*lp.beq - muL[.!i]'*lp.lb[.!i] + muU[.!i]'*lp.ub[.!i] - muL[i]'*(lp.lb[i] .* I[i]) + muU[i]'*(lp.ub[i] .* I[i]) == -lp.sense*lp.c'*v_ko)
    @constraint(model, lambda'*lp.beq - muL'*ko_lb + muU'*(ko_ub.*I) == -lp.sense*lp.c'*v_ko)

    # constraint WT to grow
    @constraint(model, lp.c' * v_wt >= min_wt_growth)

    # constraint KO to not grow
    @constraint(model, lp.c' * v_ko <= max_ko_growth)

    @objective(model, Min, sum(I[Ie]))

    return model
end

function find_conditional_media_rxns(cobra::Tiger.CobraModel, ko; min_wt_growth=1.0, max_ko_growth=0.0)
    # we're assume a non-extended Cobra model, so only Sv=0 constraints

    m, n = size(cobra.S)

    # add WT condition
    model = Model(Gurobi.Optimizer)
    v_wt = @variable(model, [1:n])
    set_name.(v_wt, cobra.rxns .* "_WT")
    set_lower_bound.(v_wt, cobra.lb)
    set_upper_bound.(v_wt, cobra.ub)
    @constraint(model, cobra.S * v_wt .== cobra.b)

    # add KO condition with strong duality
    v_ko = @variable(model, [1:n])
    set_name.(v_ko, cobra.rxns .* "_KO")
    set_lower_bound.(v_ko, cobra.lb)
    set_upper_bound.(v_ko, cobra.ub)
    @constraint(model, cobra.S * v_ko .== cobra.b)

    @variable(model, lambda[1:m])
    @variable(model, muU[1:n] >= 0)
    @variable(model, muL[1:n] >= 0)

    @constraint(model, cobra.S'*lambda - muL + muU .== cobra.c)

    @variable(model, I[1:n], binary=true)
    for i = 1:n
        # bind all variable except the KO variable
        if name(v_wt[i]) != ko * "_WT"
            @constraint(model, !I[i] => {v_wt[i] == 0.0})
            @constraint(model, !I[i] => {v_ko[i] == 0.0})
        end
    end

    for i = 1:n
        # transfer switching constraint to the dual except for
        # any reaction in the primal's objective function
        if cobra.c[i] == 0
            @constraint(model, I[i] => {muU[i] == 0.0})
            @constraint(model, I[i] => {muL[i] == 0.0})
        end
    end

    # strong duality
    @constraint(model, lambda'*cobra.b - muL'*(cobra.lb.*I) + muU'*(cobra.ub.*I) == cobra.c'*v_ko)

    # knockout KO variable
    set_upper_bound(variable_by_name(model, ko * "_KO"), 0.0)
    set_lower_bound(variable_by_name(model, ko * "_KO"), 0.0)

    # constraint WT to grow
    @constraint(model, cobra.c' * v_wt >= min_wt_growth)

    # constraint KO to not grow
    @constraint(model, cobra.c' * v_ko <= max_ko_growth)

    @objective(model, Max, sum(I))

    return model

end


# =================== TESTING ===================


cobra = read_cobra("/Users/jensen/Dropbox/research/tiger.jl/Tiger/test/cobra_extended.mat", "cobra")
#cobra = read_cobra("/Users/jensen/Dropbox/research/tiger.jl/Tiger/test/cobra_small.mat", "cobra")

model = find_conditional_media_ex(cobra, "r5", min_wt_growth=0.1)
optimize!(model)

# cobra.ub[9] = 0.1
# cobra.ub[[8 10]] .= 0.0

# #cobra = extend_cobra_cnf(cobra, ub=1.0)
# #model = build_base_model(cobra)
# #optimize!(model)
# #fitness = single_deletions(model, variable_by_name.(model, cobra.genes))

# # model1 = find_conditional_media_ex(cobra, "r6", min_wt_growth=0.1)
# # optimize!(model1)
# # model2 = find_conditional_media_rxns(cobra, "r6", min_wt_growth=0.1)
# # optimize!(model2)

# ko = "r5"
# ko_idx = findfirst(x->x==ko, cobra.vars)
# cobra.lb[ko_idx] = 0.0
# cobra.ub[ko_idx] = 0.0

# #cobra.lb[1] = 0.0
# #cobra.ub[1] = 0.0

# lp = StandardLP(cobra)
# n = length(lp.vars)

# model = Model(Gurobi.Optimizer)

# # add KO condition with primal-dual expansion
# v_ko = @variable(model, [1:n])
# set_name.(v_ko, lp.vars .* "_KO")
# set_lower_bound.(v_ko, lp.lb)
# set_upper_bound.(v_ko, lp.ub)
# #@constraint(model, lp.A * v_ko .<= lp.b)
# @constraint(model, lp.Aeq * v_ko .== lp.beq)

# @variable(model, lambda[1:size(lp.Aeq,1)])
# #@variable(model, mu[1:size(lp.A,1)] >= 0)
# @variable(model, muU[1:n] >= 0)
# @variable(model, muL[1:n] >= 0)

# #@variable(model, I[1:n], binary=true)
# #@variable(model, muI[1:n] >= 0)
# #ex_names = get_exchange_rxns(cobra, objective=false)
# # Ie = [x in ex_names for x in lp.vars]
# # fix.(I[.!Ie], 1)

# #@constraint(model, lp.Aeq'*lambda + lp.A'*mu - muL + muU .== lp.sense*lp.c)
# @constraint(model, lp.Aeq'*lambda - muL + muU .== -lp.sense*lp.c)

# @variable(model, I[1:n], binary=true)
# for i in [1 2 3]
#     #@constraint(model, v_ko[i] <= I[i])
#     @constraint(model, !I[i] => {v_ko[i] == 0.0})
#     #@constraint(model, muU[i] <= 1 - I[i])
# end
# @constraint(model, I[4:10] .== 1)


# # M = 2
# # @constraint(model, v_wt[Ie] .<= M*I[Ie])
# # @constraint(model, v_ko[Ie] .<= M*I[Ie])
# # @constraint(model, muI[Ie] .<= M*(1 .- I[Ie]))

# # strong duality

# #@constraint(model, lambda'*lp.beq - muL[.!i]'*lp.lb[.!i] + muU[.!i]'*lp.ub[.!i] - muL[i]'*(lp.lb[i] .* I[i]) + muU[i]'*(lp.ub[i] .* I[i]) == -lp.sense*lp.c'*v_ko)
# #@constraint(model, lambda'*lp.beq - muL'*lp.lb + muU'*lp.ub + muI'*(M*I) == -lp.sense*lp.c'*v_ko)
# @constraint(model, lambda'*lp.beq - muL'*lp.lb + muU'*(lp.ub.*I) == -lp.sense*lp.c'*v_ko)

# #fix(variable_by_name(model, "r1_KO"), 0, force=true)

# @constraint(model, cobra.c' * v_ko <= 0)

# @objective(model, Max, sum(I))

# optimize!(model)
# #solution_summary(model, verbose=true)
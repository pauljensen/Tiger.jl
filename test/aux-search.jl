
using Revise

using Tiger
using JuMP, Gurobi

# =================== TESTING ===================

#cobra = read_cobra("/Users/jensen/Dropbox/research/tiger.jl/Tiger/test/cobra_extended.mat", "cobra")

cobra = read_cobra("/Users/jensen/Dropbox/research/tiger.jl/Tiger/test/cobra_small.mat", "cobra")
#cobra = extend_cobra_cnf(cobra, ub=1.0)
#model = build_base_model(cobra)
#optimize!(model)
#fitness = single_deletions(model, variable_by_name.(model, cobra.genes))

function get_exchange_rxns(cobra::Tiger.CobraModel)
    cobra.rxns[(sum(abs.(cobra.S[:,1:length(cobra.rxns)]),dims=1)  .== 1)[:]]
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

function find_conditional_media(cobra::Tiger.CobraModel, ko; min_wt_growth=1.0, max_ko_growth=0.0)
    ex_names = get_exchange_rxns(cobra)
    model, Iex = create_primal_dual(cobra, bound_vars=ex_names)
    v_ko = variable_by_name.(model, cobra.rxns)
    set_upper_bound(variable_by_name(model, ko), 0.0)
    vars_ko = variable_by_name.(model, cobra.vars)
    set_name.(vars_ko, name.(vars_ko) .* "__KO")

    # add the WT condition
    build_base_model(cobra; model)
    v_wt = variable_by_name.(model, cobra.rxns)
    vars_wt = variable_by_name.(model, cobra.vars)
    set_name.(vars_wt, name.(vars_wt) .* "__WT")
    ex_wt = variable_by_name.(model, ex_names .* "__WT")
    for i = 1:length(ex_wt)
        @constraint(model, !Iex[i] => {ex_wt[i] <= 0.0})
    end

    @constraint(model, cobra.c'*vars_wt >= min_wt_growth)
    @constraint(model, cobra.c'*vars_ko <= max_ko_growth)

    @objective(model, Min, sum(Iex))

    # add WT condition
    # add KO condition with strong duality
    # constraint WT to grow
    # constraint KO to not grow
    # set activity to zero for KO gene
    # add indicators to exchanges
    # set objective to minimize media

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
        if name(v_wt[i]) != ko * "_WT"
            @constraint(model, !I[i] => {v_wt[i] == 0.0})
            @constraint(model, !I[i] => {v_ko[i] == 0.0})
        end
    end
    # #@constraint(model, v_wt .<= cobra.ub .* I)
    # #@constraint(model, v_ko .<= cobra.ub .* I)

    # #@variable(model, z[1:n])
    for i = 1:n
        #@constraint(model, z[i] == muU[i])
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

    @objective(model, Min, sum(I))

    return model

end

#model = find_conditional_media(cobra, "g6")
model = find_conditional_media_rxns(cobra, "r4")
optimize!(model)

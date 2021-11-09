
using Revise

using Tiger
using JuMP

# =================== TESTING ===================

cobra = read_cobra("/Users/jensen/Dropbox/research/tiger.jl/Tiger/test/cobra_extended.mat", "cobra")
cobra = extend_cobra_cnf(cobra, ub=1.0)
model = build_base_model(cobra)
optimize!(model)
fitness = single_deletions(model, variable_by_name.(model, cobra.genes))

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

function create_primal_dual(cobra::Tiger.CobraModel)
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

    @constraint(model, Aeq'*lambda + A'*mu - muL + muU .== c)

    # strict duality
    @constraint(model, c[1:n]'*x[1:n] == lambda'*beq + mu'*b - muL'*lb + muU'*ub)

    @objective(model, Max, 0)

    return model
end

function find_conditional_media(cobra::Tiger.CobraModel, ko; min_wt_growth=0.9, max_ko_growth=0.0)
    model = create_primal_dual(cobra)
    v_ko = variable_by_name.(model, cobra.rxns)
    set_upper_bound(variable_by_name(model, ko), 0.0)
    vars_ko = variable_by_name.(model, cobra.vars)
    set_name.(vars_ko, name.(vars_ko) .* "__KO")

    # add the WT condition
    build_base_model(cobra; model)
    v_wt = variable_by_name.(model, cobra.rxns)
    vars_wt = variable_by_name.(model, cobra.vars)
    set_name.(vars_wt, name.(vars_wt) .* "__WT")

    ex_names = get_exchange_rxns(cobra)
    ex_wt = variable_by_name.(model, ex_names .* "__WT")
    ex_ko = variable_by_name.(model, ex_names .* "__KO")
    n = length(ex_names)
    @variable(model, Iex[1:n], binary=true)
    set_name.(Iex, "I__" .* ex_names)
    for i = 1:n
        # remove this loop using index notation?
        @constraint(model, Iex[i] => {ex_wt[i] == 0.0})
        @constraint(model, Iex[i] => {ex_ko[i] == 0.0})
    end

    @constraint(model, cobra.c'*vars_wt >= min_wt_growth)
    @constraint(model, cobra.c'*vars_ko <= max_ko_growth)

    @objective(model, Max, sum(Iex))

    # add WT condition
    # add KO condition with strong duality
    # constraint WT to grow
    # constraint KO to not grow
    # set activity to zero for KO gene
    # add indicators to exchanges
    # set objective to minimize media

    return model

end

model = find_conditional_media(cobra, "g4")

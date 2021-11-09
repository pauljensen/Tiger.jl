
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
    cobra.rxns[(sum(abs.(cobra.S),dims=1)  .== 1)[:]]
end

indicator_name(varname) = "I__" * varname
function add_indicator(model::Model, var::VariableRef)
    ind = @variable(model, base_name=indicator_name(name(var)), binary=true)
    @constraint(model, ind => {var == 0.0})
    ind
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

model = create_primal_dual(cobra)

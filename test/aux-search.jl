
using Revise

using Tiger
using JuMP

# =================== TESTING ===================

cobra = read_cobra("/Users/jensen/Dropbox/research/tiger.jl/Tiger/test/cobra_extended.mat", "cobra")
model = build_base_model(cobra)
optimize!(model)

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



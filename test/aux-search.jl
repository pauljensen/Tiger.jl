
using Revise

using Tiger
import Tiger.build_base_model

using JuMP, Gurobi
using SparseArrays

import Base.convert


function get_exchange_rxns(cobra::Tiger.CobraModel; objective=true)
    # Extending a CobraModel, e.g. by adding genes via CNF, adds columns to the S matrix.
    # We only want to search in the original variables at indices 1:length(cobra.rxns).
    ex_idx = (sum(abs.(cobra.S[:,1:length(cobra.rxns)]),dims=1)  .== 1)[:]
    if !objective
        ex_idx .&= abs.(cobra.c[1:length(cobra.rxns)]) .== 0.0
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

function find_conditional_media(cobra::Tiger.CobraModel, ko; min_wt_growth=1.0, max_ko_growth=0.0, exchanges=nothing, max_kos=0, optimize=true)
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

    # find the WT growth
    @objective(model, Max, cobra.c'*v_wt)
    optimize!(model)
    wt_growth = objective_value(model)

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
    @constraint(model, lp.A * v_ko .<= lp.b)
    @constraint(model, lp.Aeq * v_ko .== lp.beq)

    @variable(model, lambda[1:size(lp.Aeq,1)])
    @variable(model, mu[1:size(lp.A,1)] >= 0)
    @variable(model, muU[1:n] >= 0)
    @variable(model, muL[1:n] >= 0)

    @variable(model, I[1:n], binary=true)

    if isnothing(exchanges)
        ex_names = get_exchange_rxns(cobra, objective=false)
    else
        ex_names = exchanges
    end
    Ie = [x in ex_names for x in lp.vars]

    gene_names = cobra.genes
    Ig = [x in gene_names for x in lp.vars]

    fix.(I[.!(Ie .| Ig)], 1)

    @constraint(model, lp.Aeq'*lambda + lp.A'*mu - muL + muU .== -lp.sense*lp.c)

    @constraint(model, v_wt[Ie] .<= lp.ub[Ie].*I[Ie])
    @constraint(model, v_ko[Ie] .<= ko_ub[Ie].*I[Ie])
    if max_kos > 0
        @constraint(model, v_wt[Ig] .<= lp.ub[Ig].*I[Ig])
        @constraint(model, v_ko[Ig] .<= ko_ub[Ig].*I[Ig])
        @constraint(model, sum(Ig) - sum(I[Ig]) <= max_kos)
    end
    #@constraint(model, muI[Ie] .<= M*(1 .- I[Ie]))

    # strong duality
    @constraint(model, lambda'*lp.beq + mu'*lp.b - muL'*ko_lb + muU'*(ko_ub.*I) == -lp.sense*lp.c'*v_ko)

    # constraint WT to grow
    @constraint(model, lp.c' * v_wt >= min_wt_growth * wt_growth)

    # constraint KO to not grow
    @constraint(model, lp.c' * v_ko <= max_ko_growth * wt_growth)

    if max_kos > 0
        @objective(model, Max, sum(I[Ig]) - sum(I[Ie]))
    else
        @objective(model, Min, sum(I[Ie]))
    end

    if !optimize
        media = nothing
        kos = nothing
    else
        optimize!(model)
        media = exchanges[value.(I[Ie]) .≈ 1.0]
        kos = gene_names[value.(I[Ig]) .≈ 0.0]
    end

    return model, media, kos
end

# =================== TESTING ===================


# cobra = read_cobra("/Users/jensen/Dropbox/research/tiger.jl/Tiger/test/cobra_or.mat", "cobra")
# #cobra = read_cobra("/Users/jensen/Dropbox/research/tiger.jl/Tiger/test/cobra_small.mat", "cobra")

# cobra = extend_cobra_cnf(cobra, ub=1.0)
# model = find_conditional_media(cobra, "g6a", min_wt_growth=0.1)
# optimize!(model)


cobra = read_cobra("/Users/jensen/Dropbox/repos/COBRA_models/iSMU.mat", "iSMU")
cobra.lb[cobra.lb .> 0.0] .= 0.0  # remove NGAM
cobra = extend_cobra_cnf(cobra, ub=1000)

# we only want the exchanges with positive fluxes; these are
# the media exchanges in iSMU
ex_names = get_exchange_rxns(cobra, objective=false)
is_ex = [x in ex_names for x in cobra.vars]
ingredients = ex_names[nonzeros(cobra.S[:,is_ex]) .> 0]
is_ingredient = [x in ingredients for x in cobra.vars]

cobra.ub[is_ingredient] .= 10
cobra.ub[cobra.ub .>  1000] .=  1000
cobra.lb[cobra.lb .< -1000] .= -1000

optimize!(build_base_model(StandardLP(cobra)))

model, media, kos = find_conditional_media(cobra, "SMU_308", min_wt_growth=0.3, exchanges=ingredients)
#optimize!(model)
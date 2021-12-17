
module Tiger

import Base: print, isempty, convert, Iterators.flatten, Iterators.enumerate
import MAT.matread
using SparseArrays
using JuMP, Gurobi
using TOML

export parse_boolean,
       parse_gpr_string,
       read_cobra,
       build_base_model,
       read_media_file,
       set_media_bounds!,
       add_gprs_cnf!,
       single_deletions,
       extend_cobra_cnf,
       get_exchange_rxns,
       StandardLP,
       build_base_model

abstract type Boolean end
abstract type Junction <: Boolean end

struct Empty <: Boolean end

struct Atom <: Boolean
    name::String
end

struct And <: Junction
    left::Boolean
    right::Boolean
end

struct Or <: Junction
    left::Boolean
    right::Boolean
end

isempty(x::Empty) = true
isempty(x::Boolean) = false

convert(::Type{String}, a::Atom) = a.name

function print(x::Boolean)
    indent_step = "   "
    aux(a::Atom, indent) = println(indent * "\$" * a.name)
    aux(b::Or, indent) = begin
        println(indent * "Or")
        aux(b.left, indent * indent_step)
        aux(b.right, indent * indent_step)
    end
    aux(b::And, indent) = begin
        println(indent * "And")
        aux(b.left, indent * indent_step)
        aux(b.right, indent * indent_step)
    end
    aux(x, "")
end

function substitute_op(e)
    if isa(e, Expr) && e.head == :call
        if e.args[1] == :&
            e.args[1] = :And
        elseif e.args[1] == :|
            e.args[1] = :Or
        end

        for i = 2:length(e.args)
            e.args[i] = substitute_op(e.args[i])
        end
    elseif isa(e, Symbol) || isa(e, String)
        e = Atom(String(e))
    end

    return e
end

function parse_boolean(s::AbstractString)
    if isempty(s)
        return Empty()
    end
    ex = Meta.parse(s)
    return eval(substitute_op(ex))
end

insert(e::Union{Atom,And}, i::Or) = Or(i, e)
insert(e::Or, i::Or) = begin
    if isa(e.left, Or)
        insert(e.left, i)
    else
        Or(Or(i, e.left), e.right)
    end
end

dnf(e::Atom) = e
dnf(e::Or) = begin
    l = dnf(e.left)
    r = dnf(e.right)
    if isa(l, Atom) || isa(l, And)
        Or(r, l)
    elseif isa(r, Atom) || isa(r, And)
        Or(l, r)
    else
        insert(l, r)
    end
end
dnf(e::And) = begin
    l = dnf(e.left)
    r = dnf(e.right)
    if !isa(l, Or) && !isa(r, Or)
        e
    elseif isa(l, Or) && !isa(r, Or)
        Or(
            dnf(And(l.left, r)),
            And(l.right, r)
        )
    elseif isa(r, Or) && !isa(l, Or)
        Or(
            dnf(And(r.left, l)),
            And(r.right, l)
        )
    else
        # Or -- Or
        sub = Or(
            Or(
                And(l.left, r.right),
                And(l.left, r.left)
            ),
            And(l.right, r.left)
        )
        Or(
            dnf(sub),
            And(l.right, r.right)
        )
    end
end

swap_andor(b::Boolean) = b
swap_andor(b::And) = Or(swap_andor(b.left), swap_andor(b.right))
swap_andor(b::Or) = And(swap_andor(b.left), swap_andor(b.right))

cnf(b::Boolean) = b |> swap_andor |> dnf |> swap_andor

atoms(e::Empty) = []
atoms(a::Atom) = [a]
atoms(b::Boolean) = unique([atoms(b.left); atoms(b.right)])

name(a::Atom) = a.name
names(a::Atom) = [a.name]
names(b::Boolean) = name.(atoms(b))

right_names(b::Boolean) = [names(b)]
right_names(j::Junction) = push!(right_names(j.left), names(j.right))

dnf_groups(b::Boolean) = begin
    groups = []
    current = b
    while isa(current, Or)
        push!(groups, names(current.right))
        current = current.left
    end
    if !isempty(current)
        push!(groups, names(current))
    end
    return groups
end

cnf_groups(b::Boolean) = b |> swap_andor |> dnf_groups

# =================== Cobra Stuff ===================

function parse_gpr_string(s; and_str="and", or_str="or")
    if and_str != "&"
        s = replace(s, and_str => "&")
    end
    if or_str != "|"
        s = replace(s, or_str => "|")
    end
    parse_boolean(s)
end

struct CobraModel 
    description::String

    c::Vector{Real}
    S::SparseMatrixCSC{Float64,Int}
    b::Vector{Real}
    lb::Vector{Real}
    ub::Vector{Real}
    csense::Vector{Char}

    rev::Vector{Bool}

    mets::Vector{String}
    metNames::Vector{String}
    rxns::Vector{String}
    rxnNames::Vector{String}

    vars::Vector{String}

    grRules::Vector{String}
    gprs::Vector{Boolean}
    genes::Vector{String}
end

"""
    read_cobra(matfile::AbstractString, name::AbstractString)

Read a Cobra model from a Matlab .mat file.

`name` is the variable name in the MAT file containing the model structure.
"""
function read_cobra(matfile::AbstractString, name::AbstractString)
    vars = matread(matfile)[name]
    gprs = parse_gpr_string.(vars["grRules"][:])
    genes = unique(flatten(atoms.(gprs)))
    CobraModel(
        vars["description"],

        vars["c"][:],
        sparse(vars["S"]),
        vars["b"][:],
        vars["lb"][:],
        vars["ub"][:],
        repeat(['='], length(vars["b"])),

        vars["rev"][:],

        vars["mets"][:],
        vars["metNames"][:],
        vars["rxns"][:],
        vars["rxnNames"][:],
        vars["rxns"][:],

        vars["grRules"][:],
        gprs,
        genes
    )
end

"""
    expand_cobra(cobra, ncons=0, nvars=0)

Add `ncons` constraints and `nvars` variables to a Cobra model.

Only the fields c, S, b, lb, ub, and csense are extended. If a variable is added,
all entries are zero. If a constraint is added, it take the form `0=0`.
"""
function expand_cobra(cobra; ncons=0, nvars=0, vars=nothing)
    m, n = size(cobra.S)

    c = zeros(n+nvars)
    c[1:n] = cobra.c

    S = spzeros(m+ncons, n+nvars)
    S[1:m,1:n] = cobra.S

    b = zeros(m+ncons)
    b[1:m] = cobra.b

    lb = zeros(n+nvars)
    lb[1:n] = cobra.lb

    ub = zeros(n+nvars)
    ub[1:n] = cobra.ub

    csense = repeat(['='], m+ncons)
    csense[1:m] = cobra.csense

    vars_ = repeat([""], n+nvars)
    vars_[1:n] = cobra.vars
    if !isnothing(vars)
        vars_[n+1:end] = vars
    end

    CobraModel(
        cobra.description,
        c,
        S,
        b,
        lb,
        ub,
        csense,
        cobra.rev,
        cobra.mets,
        cobra.metNames,
        cobra.rxns,
        cobra.rxnNames,
        vars_,
        cobra.grRules,
        cobra.gprs,
        cobra.genes
    )
end

"""
    build_base_model(cobra::CobraModel, optimizer=Gurobi.Optimizer, model=nothing)

Convert a CobraModel object into a JuMP model.

Adds variables for each reaction and defines the problem
    max c'*v
    s.t.
        S*v <|=|> b
        lb <= v <= ub

If model=nothing, a new model is created; otherwise the constraints are added to
the existing model.
"""
function build_base_model(cobra::CobraModel; optimizer=Gurobi.Optimizer, model=nothing)
    if isnothing(model)
        model = Model(optimizer)
    end

    nv = length(cobra.lb)

    #v = @variable(model, cobra.lb[i] <= [i = 1:nv] <= cobra.ub[i])
    v = @variable(model, [1:nv])
    set_name.(v, cobra.vars)
    set_lower_bound.(v, cobra.lb)
    set_upper_bound.(v, cobra.ub)

    sense = cobra.csense .== '<'
    if any(sense)
        @constraint(model, cobra.S[sense,:]*v .<= cobra.b[sense])
    end

    sense = cobra.csense .== '='
    if any(sense)
        @constraint(model, cobra.S[sense,:]*v .== cobra.b[sense])
    end

    sense = cobra.csense .== '>'
    if any(sense)
        @constraint(model, cobra.S[sense,:]*v .>= cobra.b[sense])
    end

    @objective(model, Max, cobra.c' * v)

    return model
end

"""
    add_gprs_cnf!(model::Model, cobra::CobraModel)

Add GPRs to a JuMP model using a CNF encoding.

Adds a continuous, nonnegative variable for each gene and constraints
to enforce the GPR rules. Variables are named based on the strings in 
`cobra.genes`.

Setting a gene variable to zero is analogous to a gene knockout; otherwise
the value of each gene corresponds to the magnitude of the flux through the
the associated reactions. Note that the GPR logic is encoded. If two enzymes
can catalyze a reaction, the sum of the corresponding gene variables will 
equal the reaction flux.
"""
function add_gprs_cnf!(model::Model, cobra::CobraModel)
    ng = length(cobra.genes)
    for gene in cobra.genes
        @variable(model, base_name=gene, lower_bound=0.0)
    end

    for (i, rxn) in enumerate(cobra.rxns)
        flux = variable_by_name(model, rxn)
        rule = cobra.gprs[i]
        if isempty(rule) continue end
        groups = cnf_groups(rule)
        for group in groups
            activities = variable_by_name.(model, group)
            @constraint(model, flux <= sum(activities))
            @constraint(model, flux >= -sum(activities))
        end
    end
end

function extend_cobra_cnf(cobra::CobraModel; ub=Inf)
    ng = length(cobra.genes)
    m, n = size(cobra.S)
    groups = cnf_groups.(cobra.gprs)
    ncons = 2(sum(length.(groups)))

    cobra = expand_cobra(cobra; ncons, nvars=ng, vars=cobra.genes)
    cobra.ub[n+1:end] .= ub

    gene_index(g) = findfirst(x -> x==g, cobra.vars)
    current = m
    for i = 1:n
        if isempty(groups[i]) continue end
        for group in groups[i]
            gidxs = gene_index.(group)
            current += 1
            cobra.S[current,i] = 1
            cobra.S[current,gidxs] .= -1
            cobra.csense[current] = '<'
            cobra.b[current] = 0

            current += 1
            cobra.S[current,i] = 1
            cobra.S[current,gidxs] .= 1
            cobra.csense[current] = '>'
            cobra.b[current] = 0
        end
    end

    return cobra
end


function read_media_file(mediafile::String)
    return TOML.parsefile(mediafile)
end

function set_media_bounds!(model::Model, media::Dict{String, Any}; set_defaults=true)
    if set_defaults
        varnames = all_variables(model) .|> JuMP.name
        for nm in varnames
            if occursin(media["exchange_pattern"], nm)
                var = variable_by_name(model, nm)
                set_lower_bound(var, media["default_exchange_bounds"][1])
                set_upper_bound(var, media["default_exchange_bounds"][2])
            end
        end
    end

    bounds = media["reactions"]
    for nm in keys(bounds)
        var = variable_by_name(model, nm)
        set_lower_bound(var, bounds[nm][1])
        set_upper_bound(var, bounds[nm][2])
    end
end

function set_media_bounds!(model::CobraModel, media::Dict{String, Any}; set_defaults=true)
    exchanges = get_exchange_rxns(model)
    for (i,rxn) in enumerate(model.rxns)
        if !(rxn in exchanges || haskey(media["reactions"], rxn))
            continue
        end

        if haskey(media["reactions"], rxn)
            model.lb[i] = media["reactions"][rxn][1]
            model.ub[i] = media["reactions"][rxn][2]
        elseif set_defaults
            model.lb[i] = media["default_exchange_bounds"][1]
            model.ub[i] = media["default_exchange_bounds"][2]
        end
    end
end

"""
    set_media_bounds!(model, filename::AbstractString; set_defaults=true)

Set exchange bounds using values in a TOML file.

See `./test/SampleMedia.toml` for a description of the media file format.
The media file can also be read using `read_media_file`. The resulting Dict
can be passed instead of the filename.
"""
function set_media_bounds!(model, filename::AbstractString; set_defaults=true)
    set_media_bounds!(model, read_media_file(filename); set_defaults)
end

struct Bounds
    lb::Union{Nothing,Float64}
    ub::Union{Nothing,Float64}
end

"""
    save_bounds(var::VariableRef)

Create a Bounds object with `var`'s current bounds. The a bound is `nothing` if no bound is set.
"""
function save_bounds(var::VariableRef)
    Bounds(
        has_lower_bound(var) ? lower_bound(var) : nothing,
        has_upper_bound(var) ? upper_bound(var) : nothing
    )
end


"""
    restore_bounds!(var::VariableRef, bounds::Bounds)

Set the bounds of a JuMP variable. If either bound is `nothing`, the current bound is deleted.
"""
function restore_bounds!(var::VariableRef, bounds::Bounds)
    if isnothing(bounds.lb) && has_lower_bound(var)
        delete_lower_bound(var)
    else
        set_lower_bound(var, bounds.lb)
    end

    if isnothing(bounds.ub) && has_upper_bound(var)
        delete_upper_bound(var)
    else
        set_upper_bound(var, bounds.ub)
    end
end

"""
    single_deletions(model::Model, vars=all_variables(model))

Solve the model when each variable is constrainted to zero.

Returns vectors with the objective values for each deletion and the solver status.
"""
function single_deletions(model::Model, vars=all_variables(model))
    n = length(vars)
    objval = zeros(n)
    status = Array{Any}(undef, n)

    for (i, var) in enumerate(vars)
        prev_bounds = save_bounds(var)
        set_lower_bound(var, 0.0)
        set_upper_bound(var, 0.0)
        optimize!(model)
        objval[i] = objective_value(model)
        status[i] = termination_status(model)
        restore_bounds!(var, prev_bounds)
    end

    return objval, status
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

function get_exchange_rxns(cobra::Tiger.CobraModel; objective=true)
    # Extending a CobraModel, e.g. by adding genes via CNF, adds columns to the S matrix.
    # We only want to search in the original variables at indices 1:length(cobra.rxns).
    ex_idx = (sum(abs.(cobra.S[:,1:length(cobra.rxns)]),dims=1)  .== 1)[:]
    if !objective
        ex_idx .&= abs.(cobra.c[1:length(cobra.rxns)]) .== 0.0
    end
    cobra.rxns[ex_idx]
end

end  # module

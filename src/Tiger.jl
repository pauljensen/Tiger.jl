
module Tiger

import Base.print
import MAT.matread
using SparseArrays
using JuMP, Gurobi

abstract type Boolean end
abstract type Junction <: Boolean end

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
    elseif isa(e, Symbol)
        e = Atom(String(e))
    end

    return e
end

function parse_boolean(s)
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
    push!(groups, names(current))
    return groups
end

cnf_groups(b::Boolean) = b |> swap_andor |> dnf_groups

# =================== Cobra Stuff ===================

struct CobraModel 
    description::String

    c::Vector{Real}
    S::SparseMatrixCSC{Real,Integer}
    b::Vector{Real}
    lb::Vector{Real}
    ub::Vector{Real}

    rev::Vector{Bool}

    mets::Vector{String}
    metNames::Vector{String}
    rxns::Vector{String}
    rxnNames::Vector{String}

    grRules::Vector{String}
    genes::Vector{String}
end

function read_cobra(matfile, name)
    vars = matread(matfile)[name]
    CobraModel(
        vars["description"],

        vars["c"][:],
        sparse(vars["S"]),
        vars["b"][:],
        vars["lb"][:],
        vars["ub"][:],

        vars["rev"][:],

        vars["mets"][:],
        vars["metNames"][:],
        vars["rxns"][:],
        vars["rxnNames"][:],

        vars["grRules"][:],
        vars["genes"][:]
    )
end

function build_base_model(cobra::CobraModel, optimizer=Gurobi.Optimizer)
    model = Model(optimizer)

    nr = length(cobra.lb)
    @variable(model, cobra.lb[i] <= v[i = 1:nr] <= cobra.ub[i])
    @constraint(model, cobra.S * v .== cobra.b)
    @objective(model, Max, cobra.c' * v)

    return model
end

function add_gprs_cnf!(model, cobra::CobraModel)
    gene_ub = 1e10

    # add gene variables to model

    for i = 1:length(cobra.rxns)
        if isempty(cobra.grRules[i])
            continue
        end
        1
    end
end

end  # module

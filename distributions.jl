using StatsBase
using Distributions


abstract type PM end


struct Discr <: PM
    # given n_atoms constructs type discrete measure which has atoms and weights
    # sample space is [-1,1]
    atoms::Vector{Float64}
    weights::Vector{Float64}

    function Discr(n_atoms::Int)
        atoms = 2*rand(n_atoms).-1
        weights = rand(n_atoms)
        weights = weights./sum(weights)
        new(atoms, weights)
    end

    function Discr(atoms::Vector{Float64}, weights::Vector{Float64})
        new(atoms, weights)
    end
end

struct Gammarv<: PM
    # given α shape and θ scale parameters and b  constructs measure on [0,b]
    # firstly simulate gamma r.v. and then if it's larger than b set it uniformly on [0,b]
    α::Float64
    θ::Float64
    b::Float64
end

struct Uniformrv<: PM
    # uniform pm on [a,b]
    a::Float64
    b::Float64
end


function Base.copy(d::Discr)
    return Discr(d.atoms, d.weights)
end

function convert_measure(atoms::Vector{Float64})
    # this function converts discrete measure into an another same discrete 
    # measure which has unique atoms
    uniq_atoms = unique(atoms)
    n = length(atoms)
    counts = countmap(atoms)
    weights = [counts[val]/n for val in uniq_atoms]
    return uniq_atoms, weights
end


function emp_distr(pm::Discr, n::Int)
    # returns the empirical distribution of discrete p.m., i.e. unique atoms and associated weights from
    # simulation
    dist = Categorical(pm.weights)
    x = pm.atoms[rand(dist,n)]
    return x
end

function emp_distr(pm::Gammarv, n::Int)
    # returns the empirical distribution of gamma p.m., i.e. unique atoms and associated weights from
    # simulation
    dist = Gamma(pm.α,pm.θ)
    x = rand(dist,n)
    for i in 1:length(x)
        if x[i] > pm.b
            x[i] = rand(Uniform(0,pm.b))
        end
    end
    return x
end

function emp_distr(pm::Uniformrv, n::Int)
    # returns the empirical distribution of uniform(a,b) p.m., i.e. unique atoms and associated weights from
    # simulation
    dist = Uniform(pm.a,pm.b)
    x = rand(dist,n)
    return x
end










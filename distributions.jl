using StatsBase
using Distributions


struct Discr
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

function Base.copy(d::Discr)
    return Discr(d.atoms, d.weights)
end

convert_measure = function(atoms::Vector{Float64})
    # this function converts discrete measure into an another same discrete 
    # measure which has unique atoms
    uniq_atoms = unique(atoms)
    n = length(atoms)
    counts = countmap(atoms)
    weights = [counts[val]/n for val in uniq_atoms]
    return uniq_atoms, weights
end


emp_distr = function(pm::Discr, n::Int)
    # returns the empirical distribution of discrete p.m., i.e. unique atoms and associated weights from
    # simulation
    dist = Categorical(pm.weights)
    x = pm.atoms[rand(dist,n)]
    return convert_measure(x)
end









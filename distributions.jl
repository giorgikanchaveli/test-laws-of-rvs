using StatsBase
using Distributions


abstract type PM end


struct Discr <: PM
    # given n_atoms constructs type discrete measure which has atoms and weights
    # sample space is [a,b]
    atoms::Vector{Float64}
    weights::Vector{Float64}
    a::Float64
    b::Float64

    function Discr(n_atoms::Int,a::Float64,b::Float64)
        atoms = rand(Uniform(a,b),n_atoms)
        weights = rand(n_atoms)
        weights = weights ./ sum(weights)
        new(atoms, weights,a,b)
    end

    function Discr(atoms::Vector{Float64}, weights::Vector{Float64},a::Float64,b::Float64)
        new(atoms, weights,a,b)
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

struct Betarv<:PM
    # Beta distribution on [0,1]
    α::Float64
    β::Float64
    a::Float64
    b::Float64
    function Betarv(α::Float64,β::Float64)
        new(α,β,0.0,1.0)
    end
end




function Base.copy(d::Discr)
    return Discr(d.atoms, d.weights,d.a,d.b)
end

function unique_measure(atoms::Vector{Float64})
    # this function converts discrete measure into an another same discrete 
    # measure which has unique atoms
    uniq_atoms = unique(atoms)
    n = length(atoms)
    counts = countmap(atoms)
    weights = [counts[val]/n for val in uniq_atoms]
    return uniq_atoms, weights
end

function project_grid(x::Vector{Float64},a::Float64,b::Float64)
    # We project atoms on grid and updage weights accordingly
    M = 200
    Δ = (b-a)/M
    atoms = collect(a:Δ:b)
    weights = zeros(M+1)
    for i in 1:length(x)
        l = floor(Int,(x[i]-a)/Δ+1)
        if x[i] - (l-1)*Δ - a > Δ/2
            j = l+1
        else
            j = l
        end
        weights[j]+=1
    end
    weights = weights ./ sum(weights)
    return atoms, weights
end

function emp_distr(pm::Discr, n::Int)
    # returns the empirical distribution of discrete p.m., i.e. unique atoms and associated weights from
    # simulation
    dist = Categorical(pm.weights)
    x = pm.atoms[rand(dist,n)]
    return unique_measure(x)
end

function emp_distr(pm::Betarv, n)
    dist = Beta(pm.α,pm.β)
    x = rand(dist, n)
    if n<=200
        return x, fill(1/n,n)
    else
        return project_grid(x, 0.0,1.0)
    end
end



# unda davamato projection on grid
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
    return x,fill(1/n,n)
end

function emp_distr(pm::Uniformrv, n::Int)
    # returns the empirical distribution of uniform(a,b) p.m., i.e. unique atoms and associated weights from
    # simulation
    dist = Uniform(pm.a,pm.b)
    x = rand(dist,n)
    if n<=200
        return x,fill(1/n,n)
    else
        return project_grid(x, pm.a,pm.b)
    end
end


function CDF(x::Float64, pm::Discr)
    # CDF of discrete distribution
    y = 0.0
    if x>=pm.b
        return 1.0
    elseif x<= pm.a
        return 0.0
    else
        for i in 1:length(pm.atoms)
            if pm.atoms[i]<=x
                y+=pm.weights[i]
            end
        end
    end
    return y
end



function CDF(x::Float64, pm::Uniformrv)
    # CDF of discrete distribution
    if x <= pm.a
        return 0.0
    elseif x>= pm.b
        return 1.0
    else
        return (x-pm.a)/(pm.b-pm.a)
    end
end


function CDF(x::Float64, pm::Betarv)
    # CDF of discrete distribution
    d = Beta(pm.α,pm.β)
    if x <= 0.0
        return 0.0
    elseif x>1.0
        return 1.0
    else
        return quadgk(x->pdf(d,x),0.0,x)[1]
    end
end








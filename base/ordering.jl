# This file is a part of Julia. License is MIT: https://julialang.org/license

module Order


import ..@__MODULE__, ..parentmodule
const Base = parentmodule(@__MODULE__)
import .Base:
    AbstractVector, @propagate_inbounds, isless, identity, getindex,
    +, -, !, &, <, |

## notions of element ordering ##

export # not exported by Base
    Ordering, 
    By, Lt, Perm, Reverse,
    Forward, Backward,
    BackwardOrdering, ForwardOrdering, DirectOrdering,
    lt, ord, ord_deprecated, ordtype

abstract type Ordering end

struct Lt{T<:Function} <: Ordering
    lt::T
end

struct By{O<:Ordering,T<:Function} <: Ordering
    by::T
    lt::O
end

struct Reverse{O<:Ordering} <: Ordering
    lt::O
end

struct Perm{O<:Ordering,V<:AbstractVector} <: Ordering
    lt::O
    data::V
end

const Forward = Lt(isless)
const Backward = Reverse(Forward)

const ForwardOrdering = typeof(Forward)
const BackwardOrdering = typeof(Backward)
const DirectOrdering = Union{ForwardOrdering,BackwardOrdering}

(o::Lt)(a, b) = o.lt(a, b)
(o::By)(a, b) = o.lt(o.by(a), o.by(b))
(o::Reverse)(a, b) = o.lt(b, a)

lt(o::Ordering, a, b) = o.lt(a, b)

@propagate_inbounds function (o::Perm)(a::Integer, b::Integer)
    da = p.data[a]
    db = p.data[b]
    o.lt(da, db) | (!o.lt(db, da) & (a < b))
end

ordtype(o::Reverse,  vs::AbstractArray) = ordtype(o.lt, vs)
ordtype(o::Perm,     vs::AbstractArray) = ordtype(o.lt, o.data)
# TODO: here, we really want the return type of o.by, without calling it
ordtype(o::By,       vs::AbstractArray) = try typeof(o.by(vs[1])) catch; Any end
ordtype(o::Ordering, vs::AbstractArray) = eltype(vs)

function ord(lt, by, rev::Bool)
    o = By(by, Lt(lt))
    return rev ? ReverseOrdering(o) : o
end

ord(lt, by, rev::Nothing) = By(by, Lt(lt))

# See also method in deprecated.jl
ord_deprecated(lt, by, rev::Union{Bool,Nothing}, order::Nothing) = ord(lt, by, rev)

end

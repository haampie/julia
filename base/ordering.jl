# This file is a part of Julia. License is MIT: https://julialang.org/license

module Order


import ..@__MODULE__, ..parentmodule
const Base = parentmodule(@__MODULE__)
import .Base: AbstractVector, @propagate_inbounds, @inline, 
    isless, identity, !, &, <, |

## notions of element ordering ##

export # not exported by Base
    Ordering, Forward, Backward,
    By, Less, Reverse, Perm,
    ForwardOrdering, BackwardOrdering, DirectOrdering

abstract type Ordering end

struct Less{T<:Function} <: Ordering
    isless::T
end

const Forward = Less(isless)

struct By{T<:Function,O<:Ordering} <: Ordering
    by::T
    isless::O
end

struct Reverse{O<:Ordering} <: Ordering
    isless::O
end

struct Perm{O<:Ordering,V<:AbstractVector} <: Ordering
    isless::O
    data::V
end

@inline Reverse(o::Reverse) = o.isless
@inline By(by::typeof(identity), o::Ordering) = o
@inline By(by::Function) = By(by, Forward)
@inline Perm(data::AbstractVector) = Perm(Forward, data)

const Backward = Reverse(Forward)

const ForwardOrdering = typeof(Forward)
const BackwardOrdering = typeof(Backward)
const DirectOrdering = Union{ForwardOrdering,BackwardOrdering}

@inline (o::Less)(a, b) = o.isless(a, b)
@inline (o::By)(a, b) = o.isless(o.by(a), o.by(b))
@inline (o::Reverse)(a, b) = o.isless(b, a)

@propagate_inbounds function (o::Perm)(a, b)
    da, db = o.data[a], o.data[b]
    o.isless(da, db) | (!o.isless(db, da) & (a < b))
end

end

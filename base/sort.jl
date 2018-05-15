# This file is a part of Julia. License is MIT: https://julialang.org/license

module Sort

import ..@__MODULE__, ..parentmodule
const Base = parentmodule(@__MODULE__)
using .Base.Order
using .Base: copymutable, LinearIndices, IndexStyle, viewindexing, IndexLinear, _length, (:),
    eachindex, axes, first, last, similar, start, next, done, zip, @views, OrdinalRange,
    AbstractVector, @inbounds, AbstractRange, @eval, @inline, Vector, @noinline,
    AbstractMatrix, AbstractUnitRange, isless, identity, eltype, >, <, <=, >=, |, +, -, *, !,
    extrema, sub_with_overflow, add_with_overflow, oneunit, div, getindex, setindex!,
    length, resize!, fill

using .Base: >>>, !==

import .Base:
    sort,
    sort!,
    issorted,
    sortperm,
    Slice,
    to_indices

export # also exported by Base
    # order-only:
    issorted,
    searchsorted,
    searchsortedfirst,
    searchsortedlast,
    # order & algorithm:
    sort,
    sort!,
    sortperm,
    sortperm!,
    partialsort,
    partialsort!,
    partialsortperm,
    partialsortperm!,
    sortrows,
    sortcols,
    # algorithms:
    InsertionSort,
    QuickSort,
    MergeSort,
    PartialQuickSort

export # not exported by Base
    Algorithm,
    DEFAULT_UNSTABLE,
    DEFAULT_STABLE,
    SMALL_ALGORITHM,
    SMALL_THRESHOLD


## functions requiring only ordering ##

"""
    issorted(v, order::Ordering=Forward)

Test whether a vector is in sorted order.

# Examples
```jldoctest
julia> issorted([1, 2, 3])
true

julia> issorted([(1, "b"), (2, "a")], By(x -> x[1]))
true

julia> issorted([(1, "b"), (2, "a")], By(x -> x[2]))
false

julia> issorted([(1, "b"), (2, "a")], Reverse(By(x -> x[2])))
true
```
"""
function issorted(itr, lt::Ordering = Forward)
    state = start(itr)
    done(itr,state) && return true
    prev, state = next(itr, state)
    while !done(itr, state)
        this, state = next(itr, state)
        lt(this, prev) && return false
        prev = this
    end
    return true
end

"""
    partialsort!(v, k, o::Ordering=Forward)

Partially sort the vector `v` in place, according to the order specified by `o` so that the
value at index `k` (or range of adjacent values if `k` is a range) occurs
at the position where it would appear if the array were fully sorted via a non-stable
algorithm. If `k` is a single index, that value is returned; if `k` is a range, an array of
values at those indices is returned. Note that `partialsort!` does not fully sort the input
array.

# Examples
```jldoctest
julia> a = [1, 2, 4, 3, 4]
5-element Array{Int64,1}:
 1
 2
 4
 3
 4

julia> partialsort!(a, 4)
4

julia> a
5-element Array{Int64,1}:
 1
 2
 3
 4
 4

julia> a = [1, 2, 4, 3, 4]
5-element Array{Int64,1}:
 1
 2
 4
 3
 4

julia> partialsort!(a, 4, Backward)
2

julia> a
5-element Array{Int64,1}:
 4
 4
 3
 2
 1
```
"""
function partialsort!(v::AbstractVector, k::Union{Int,OrdinalRange}, o::Ordering = Forward)
    inds = axes(v, 1)
    sort!(v, first(inds), last(inds), PartialQuickSort(k), o)
    maybeview(v, k)
end

maybeview(v, k) = view(v, k)
maybeview(v, k::Integer) = v[k]

"""
    partialsort(v, k, o::Ordering = Forward)

Variant of [`partialsort!`](@ref) which copies `v` before partially sorting it, thereby returning the
same thing as `partialsort!` but leaving `v` unmodified.
"""
partialsort(v::AbstractVector, k::Union{Int,OrdinalRange}, o::Ordering = Forward) =
    partialsort!(copymutable(v), k, o)


# reference on sorted binary search:
#   http://www.tbray.org/ongoing/When/200x/2003/03/22/Binary

# index of the first value of vector a that is greater than or equal to x;
# returns length(v)+1 if x is greater than all values in v.
function searchsortedfirst(v::AbstractVector, x, lo::Int, hi::Int, lt::Ordering = Forward)
    lo = lo-1
    hi = hi+1
    @inbounds while lo < hi-1
        m = (lo+hi)>>>1
        if lt(v[m], x)
            lo = m
        else
            hi = m
        end
    end
    return hi
end

# index of the last value of vector a that is less than or equal to x;
# returns 0 if x is less than all values of v.
function searchsortedlast(v::AbstractVector, x, lo::Int, hi::Int, lt::Ordering = Forward)
    lo = lo-1
    hi = hi+1
    @inbounds while lo < hi-1
        m = (lo+hi)>>>1
        if lt(x, v[m])
            hi = m
        else
            lo = m
        end
    end
    return lo
end

# returns the range of indices of v equal to x
# if v does not contain x, returns a 0-length range
# indicating the insertion point of x
function searchsorted(v::AbstractVector, x, ilo::Int, ihi::Int, lt::Ordering = Forward)
    lo = ilo-1
    hi = ihi+1
    @inbounds while lo < hi-1
        m = (lo+hi)>>>1
        if lt(v[m], x)
            lo = m
        elseif lt(x, v[m])
            hi = m
        else
            a = searchsortedfirst(v, x, max(lo,ilo), m, lt)
            b = searchsortedlast(v, x, m, min(hi,ihi), lt)
            return a : b
        end
    end
    return (lo + 1) : (hi - 1)
end

function searchsortedlast(a::AbstractRange{<:Real}, x::Real, lt::DirectOrdering = Forward)
    if step(a) == 0
        lt(x, first(a)) ? 0 : length(a)
    else
        n = round(Integer, clamp((x - first(a)) / step(a) + 1, 1, length(a)))
        lt(x, a[n]) ? n - 1 : n
    end
end

function searchsortedfirst(a::AbstractRange{<:Real}, x::Real, lt::DirectOrdering = Forward)
    if step(a) == 0
        lt(first(a), x) ? length(a) + 1 : 1
    else
        n = round(Integer, clamp((x - first(a)) / step(a) + 1, 1, length(a)))
        lt(a[n] ,x) ? n + 1 : n
    end
end

function searchsortedlast(a::AbstractRange{<:Integer}, x::Real, lt::DirectOrdering = Forward)
    if step(a) == 0
        lt(x, first(a)) ? 0 : length(a)
    else
        clamp( fld(floor(Integer, x) - first(a), step(a)) + 1, 0, length(a))
    end
end

function searchsortedfirst(a::AbstractRange{<:Integer}, x::Real, lt::DirectOrdering = Forward)
    if step(a) == 0
        lt(first(a), x) ? length(a)+1 : 1
    else
        clamp(-fld(floor(Integer, -x) + first(a), step(a)) + 1, 1, length(a) + 1)
    end
end

function searchsortedfirst(a::AbstractRange{<:Integer}, x::Unsigned, lt::DirectOrdering = Forward)
    if lt(first(a), x)
        if step(a) == 0
            length(a) + 1
        else
            min(cld(x - first(a), step(a)), length(a)) + 1
        end
    else
        1
    end
end

function searchsortedlast(a::AbstractRange{<:Integer}, x::Unsigned, lt::DirectOrdering = Forward)
    if lt(x, first(a))
        0
    elseif step(a) == 0
        length(a)
    else
        min(fld(x - first(a), step(a)) + 1, length(a))
    end
end

searchsorted(a::AbstractRange{<:Real}, x::Real, lt::DirectOrdering = Forward) =
    searchsortedfirst(a, x, lt) : searchsortedlast(a, x, lt)

for s in (:searchsortedfirst, :searchsortedlast, :searchsorted)
    @eval $s(v::AbstractVector, x, o::Ordering = Forward) = (inds = axes(v, 1); $s(v,x,first(inds),last(inds),o))
end

"""
    searchsorted(a, x, order::Ordering = Forward)

Return the range of indices of `a` which compare as equal to `x` (using binary search)
according to the order specified by the `by`, `lt` and `rev` keywords, assuming that `a`
is already sorted in that order. Return an empty range located at the insertion point
if `a` does not contain values equal to `x`.

# Examples
```jldoctest
julia> a = [4, 3, 2, 1]
4-element Array{Int64,1}:
 4
 3
 2
 1

julia> searchsorted(a, 4)
5:4

julia> searchsorted(a, 4, Backward)
1:1
```
""" searchsorted

"""
    searchsortedfirst(a, x, order::Ordering = Forward)

Return the index of the first value in `a` greater than or equal to `x`, according to the
specified order. Return `length(a) + 1` if `x` is greater than all values in `a`.
`a` is assumed to be sorted.

# Examples
```jldoctest
julia> searchsortedfirst([1, 2, 4, 5, 14], 4)
3

julia> searchsortedfirst([1, 2, 4, 5, 14], 4, Backward)
1

julia> searchsortedfirst([1, 2, 4, 5, 14], 15)
6
```
""" searchsortedfirst

"""
    searchsortedlast(a, x, order::Ordering = Forward)

Return the index of the last value in `a` less than or equal to `x`, according to the
specified order. Return `0` if `x` is less than all values in `a`. `a` is assumed to
be sorted.

# Examples
```jldoctest
julia> searchsortedlast([1, 2, 4, 5, 14], 4)
3

julia> searchsortedlast([1, 2, 4, 5, 14], 4, Backward)
5

julia> searchsortedlast([1, 2, 4, 5, 14], -1)
0
```
""" searchsortedlast


## sorting algorithms ##

abstract type Algorithm end

struct InsertionSortAlg <: Algorithm end
struct QuickSortAlg     <: Algorithm end
struct MergeSortAlg     <: Algorithm end

struct PartialQuickSort{T <: Union{Int,OrdinalRange}} <: Algorithm
    k::T
end

Base.first(a::PartialQuickSort{Int}) = 1
Base.last(a::PartialQuickSort{Int}) = a.k
Base.first(a::PartialQuickSort) = first(a.k)
Base.last(a::PartialQuickSort) = last(a.k)

const InsertionSort = InsertionSortAlg()
const QuickSort     = QuickSortAlg()
const MergeSort     = MergeSortAlg()

const DEFAULT_UNSTABLE = QuickSort
const DEFAULT_STABLE   = MergeSort
const SMALL_ALGORITHM  = InsertionSort
const SMALL_THRESHOLD  = 20

function sort!(v::AbstractVector, lo::Int, hi::Int, ::InsertionSortAlg, lt::Ordering = Forward)
    @inbounds for i = lo+1:hi
        j = i
        x = v[i]
        while j > lo
            if lt(x, v[j-1])
                v[j] = v[j-1]
                j -= 1
                continue
            end
            break
        end
        v[j] = x
    end
    return v
end

# selectpivot!
#
# Given 3 locations in an array (lo, mi, and hi), sort v[lo], v[mi], v[hi]) and
# choose the middle value as a pivot
#
# Upon return, the pivot is in v[lo], and v[hi] is guaranteed to be
# greater than the pivot

@inline function selectpivot!(v::AbstractVector, lo::Int, hi::Int, lt::Ordering = Forward)
    @inbounds begin
        mi = (lo+hi)>>>1

        # sort the values in v[lo], v[mi], v[hi]

        if lt(v[mi], v[lo])
            v[mi], v[lo] = v[lo], v[mi]
        end
        if lt(v[hi], v[mi])
            if lt(v[hi], v[lo])
                v[lo], v[mi], v[hi] = v[hi], v[lo], v[mi]
            else
                v[hi], v[mi] = v[mi], v[hi]
            end
        end

        # move v[mi] to v[lo] and use it as the pivot
        v[lo], v[mi] = v[mi], v[lo]
        pivot = v[lo]
    end

    # return the pivot
    return pivot
end

# partition!
#
# select a pivot, and partition v according to the pivot

function partition!(v::AbstractVector, lo::Int, hi::Int, lt::Ordering = Forward)
    pivot = selectpivot!(v, lo, hi, lt)
    # pivot == v[lo], v[hi] > pivot
    i, j = lo, hi
    @inbounds while true
        i += 1; j -= 1
        while lt(v[i], pivot); i += 1; end;
        while lt(pivot, v[j]); j -= 1; end;
        i >= j && break
        v[i], v[j] = v[j], v[i]
    end
    v[j], v[lo] = pivot, v[j]

    # v[j] == pivot
    # v[k] >= pivot for k > j
    # v[i] <= pivot for i < j
    return j
end

function sort!(v::AbstractVector, lo::Int, hi::Int, a::QuickSortAlg, lt::Ordering = Forward)
    @inbounds while lo < hi
        hi-lo <= SMALL_THRESHOLD && return sort!(v, lo, hi, SMALL_ALGORITHM, lt)
        j = partition!(v, lo, hi, lt)
        if j-lo < hi-j
            # recurse on the smaller chunk
            # this is necessary to preserve O(log(n))
            # stack space in the worst case (rather than O(n))
            lo < (j-1) && sort!(v, lo, j-1, a, lt)
            lo = j+1
        else
            j+1 < hi && sort!(v, j+1, hi, a, lt)
            hi = j-1
        end
    end
    return v
end

function sort!(v::AbstractVector, lo::Int, hi::Int, a::MergeSortAlg, lt::Ordering, t=similar(v,0))
    @inbounds if lo < hi
        hi-lo <= SMALL_THRESHOLD && return sort!(v, lo, hi, SMALL_ALGORITHM, lt)

        m = (lo+hi)>>>1
        (length(t) < m-lo+1) && resize!(t, m-lo+1)

        sort!(v, lo,  m,  a, lt, t)
        sort!(v, m+1, hi, a, lt, t)

        i, j = 1, lo
        while j <= m
            t[i] = v[j]
            i += 1
            j += 1
        end

        i, k = 1, lo
        while k < j <= hi
            if lt(v[j], t[i])
                v[k] = v[j]
                j += 1
            else
                v[k] = t[i]
                i += 1
            end
            k += 1
        end
        while k < j
            v[k] = t[i]
            k += 1
            i += 1
        end
    end

    return v
end

## TODO: When PartialQuickSort is parameterized by an Int, this version of sort
##       has one less comparison per loop than the version below, but enabling
##       it causes return type inference to fail for sort/sort! (#12833)
##
# function sort!(v::AbstractVector, lo::Int, hi::Int, a::PartialQuickSort{Int},
#                o::Ordering)
#     @inbounds while lo < hi
#         hi-lo <= SMALL_THRESHOLD && return sort!(v, lo, hi, SMALL_ALGORITHM, o)
#         j = partition!(v, lo, hi, o)
#         if j >= a.k
#             # we don't need to sort anything bigger than j
#             hi = j-1
#         elseif j-lo < hi-j
#             # recurse on the smaller chunk
#             # this is necessary to preserve O(log(n))
#             # stack space in the worst case (rather than O(n))
#             lo < (j-1) && sort!(v, lo, j-1, a, o)
#             lo = j+1
#         else
#             (j+1) < hi && sort!(v, j+1, hi, a, o)
#             hi = j-1
#         end
#     end
#     return v
# end


function sort!(v::AbstractVector, lo::Int, hi::Int, a::PartialQuickSort,
               lt::Ordering = Forward)
    @inbounds while lo < hi
        hi-lo <= SMALL_THRESHOLD && return sort!(v, lo, hi, SMALL_ALGORITHM, lt)
        j = partition!(v, lo, hi, lt)

        if j <= first(a)
            lo = j+1
        elseif j >= last(a)
            hi = j-1
        else
            if j-lo < hi-j
                lo < (j-1) && sort!(v, lo, j-1, a, lt)
                lo = j+1
            else
                hi > (j+1) && sort!(v, j+1, hi, a, lt)
                hi = j-1
            end
        end
    end
    return v
end


## generic sorting methods ##

defalg(v::AbstractArray) = DEFAULT_STABLE
defalg(v::AbstractArray{<:Number}) = DEFAULT_UNSTABLE

function sort!(v::AbstractVector, alg::Algorithm, order::Ordering = Forward)
    inds = axes(v,1)
    sort!(v,first(inds),last(inds),alg,order)
end

"""
    sort!(v, order::Ordering=Forward; alg::Algorithm=defalg(v))

Sort the vector `v` in place. `QuickSort` is used by default for numeric arrays while
`MergeSort` is used for other arrays. You can specify an algorithm to use via the `alg`
keyword (see Sorting Algorithms for available algorithms).

# Examples
```jldoctest
julia> v = [3, 1, 2]; sort!(v); v
3-element Array{Int64,1}:
 1
 2
 3

julia> v = [3, 1, 2]; sort!(v, Backward); v
3-element Array{Int64,1}:
 3
 2
 1

julia> v = [(1, "c"), (3, "a"), (2, "b")]; sort!(v, By(x -> x[1])); v
3-element Array{Tuple{Int64,String},1}:
 (1, "c")
 (2, "b")
 (3, "a")

julia> v = [(1, "c"), (3, "a"), (2, "b")]; sort!(v, By(x -> x[2])); v
3-element Array{Tuple{Int64,String},1}:
 (3, "a")
 (2, "b")
 (1, "c")
```
"""
function sort!(v::AbstractVector, order::Ordering = Forward; alg::Algorithm=defalg(v))
    if order === Forward && isa(v,Vector) && eltype(v)<:Integer
        n = _length(v)
        if n > 1
            min, max = extrema(v)
            (diff, o1) = sub_with_overflow(max, min)
            (rangelen, o2) = add_with_overflow(diff, oneunit(diff))
            if !o1 && !o2 && rangelen < div(n,2)
                return sort_int_range!(v, rangelen, min)
            end
        end
    end
    sort!(v, alg, order)
end

# sort! for vectors of few unique integers
function sort_int_range!(x::Vector{<:Integer}, rangelen, minval)
    offs = 1 - minval
    n = length(x)

    where = fill(0, rangelen)
    @inbounds for i = 1:n
        where[x[i] + offs] += 1
    end

    idx = 1
    @inbounds for i = 1:rangelen
        lastidx = idx + where[i] - 1
        val = i-offs
        for j = idx:lastidx
            x[j] = val
        end
        idx = lastidx + 1
    end

    return x
end

"""
    sort(v, order::Ordering=Forward; alg::Algorithm=defalg(v))

Variant of [`sort!`](@ref) that returns a sorted copy of `v` leaving `v` itself unmodified.

# Examples
```jldoctest
julia> v = [3, 1, 2];

julia> sort(v)
3-element Array{Int64,1}:
 1
 2
 3

julia> v
3-element Array{Int64,1}:
 3
 1
 2
```
"""
sort(v::AbstractVector, o::Ordering = Forward; kws...) = sort!(copymutable(v), o; kws...)

## partialsortperm: the permutation to sort the first k elements of an array ##

"""
    partialsortperm(v, k; alg=<algorithm>)

Return a partial permutation of the vector `v`, according to the order specified by
`by`, `lt` and `rev`, so that `v[output]` returns the first `k` (or range of adjacent values
if `k` is a range) values of a fully sorted version of `v`. If `k` is a single index,
the index in `v` of the value which would be sorted at position `k` is returned;
if `k` is a range, an array with the indices in `v` of the values which would be sorted in
these positions is returned.

Note that this is equivalent to, but more efficient than, calling `sortperm(...)[k]`.
"""
partialsortperm(v::AbstractVector, k::Union{Integer,OrdinalRange}, order::Ordering=Forward; kwargs...) =
    partialsortperm!(similar(Vector{eltype(k)}, axes(v,1)), v, k, order; kwargs..., initialized=false)

"""
    partialsortperm!(ix, v, k, order::Ordering; initialized=false)

Like [`partialsortperm`](@ref), but accepts a preallocated index vector `ix`. If `initialized` is `false`
(the default), `ix` is initialized to contain the values `1:length(ix)`.
"""
function partialsortperm!(ix::AbstractVector{<:Integer}, v::AbstractVector,
                          k::Union{Int, OrdinalRange}, order::Ordering=Forward;
                          initialized::Bool=false)
    if !initialized
        @inbounds for i = axes(ix,1)
            ix[i] = i
        end
    end

    # do partial quicksort
    sort!(ix, PartialQuickSort(k), Perm(order, v))

    maybeview(ix, k)
end

## sortperm: the permutation to sort an array ##

"""
    sortperm(v, order::Ordering=Forward; alg::Algorithm=DEFAULT_UNSTABLE)

Return a permutation vector of indices of `v` that puts it in sorted order. Specify `alg` to
choose a particular sorting algorithm (see Sorting Algorithms). `MergeSort` is used by
default, and since it is stable, the resulting permutation will be the lexicographically
first one that puts the input array into sorted order – i.e. indices of equal elements
appear in ascending order. If you choose a non-stable sorting algorithm such as `QuickSort`,
a different permutation that puts the array into order may be returned. The order is
specified using the same keywords as `sort!`.

See also [`sortperm!`](@ref).

# Examples
```jldoctest
julia> v = [3, 1, 2];

julia> p = sortperm(v)
3-element Array{Int64,1}:
 2
 3
 1

julia> v[p]
3-element Array{Int64,1}:
 1
 2
 3
```
"""
function sortperm(v::AbstractVector, order::Ordering=Forward;
                  alg::Algorithm=DEFAULT_UNSTABLE)
    if order === Forward && isa(v,Vector) && eltype(v)<:Integer
        n = _length(v)
        if n > 1
            min, max = extrema(v)
            (diff, o1) = sub_with_overflow(max, min)
            (rangelen, o2) = add_with_overflow(diff, oneunit(diff))
            if !o1 && !o2 && rangelen < div(n,2)
                return sortperm_int_range(v, rangelen, min)
            end
        end
    end
    p = similar(Vector{Int}, axes(v, 1))
    for (i,ind) in zip(eachindex(p), axes(v, 1))
        p[i] = ind
    end
    sort!(p, alg, Perm(order,v))
end


"""
    sortperm!(ix, v, order::Ordering=Forward; alg::Algorithm=DEFAULT_UNSTABLE, initialized::Bool=false)

Like [`sortperm`](@ref), but accepts a preallocated index vector `ix`.  If `initialized` is `false`
(the default), `ix` is initialized to contain the values `1:length(v)`.

# Examples
```jldoctest
julia> v = [3, 1, 2]; p = zeros(Int, 3);

julia> sortperm!(p, v); p
3-element Array{Int64,1}:
 2
 3
 1

julia> v[p]
3-element Array{Int64,1}:
 1
 2
 3
```
"""
function sortperm!(x::AbstractVector{<:Integer}, v::AbstractVector, order::Ordering=Forward;
                   alg::Algorithm=DEFAULT_UNSTABLE,
                   initialized::Bool=false)
    if axes(x,1) != axes(v,1)
        throw(ArgumentError("index vector must have the same indices as the source vector, $(axes(x,1)) != $(axes(v,1))"))
    end
    if !initialized
        @inbounds for i = axes(v,1)
            x[i] = i
        end
    end
    sort!(x, alg, Perm(order,v))
end

# sortperm for vectors of few unique integers
function sortperm_int_range(x::Vector{<:Integer}, rangelen, minval)
    offs = 1 - minval
    n = length(x)

    where = fill(0, rangelen+1)
    where[1] = 1
    @inbounds for i = 1:n
        where[x[i] + offs + 1] += 1
    end

    #cumsum!(where, where)
    @inbounds for i = 2:length(where)
        where[i] += where[i-1]
    end

    P = Vector{Int}(undef, n)
    @inbounds for i = 1:n
        label = x[i] + offs
        P[where[label]] = i
        where[label] += 1
    end

    return P
end

## sorting multi-dimensional arrays ##

"""
    sort(A, o::Ordering = Forward; dims::Integer, alg::Algorithm=DEFAULT_UNSTABLE)

Sort a multidimensional array `A` along the given dimension.
See [`sort!`](@ref) for a description of possible
keyword arguments.

# Examples
```jldoctest
julia> A = [4 3; 1 2]
2×2 Array{Int64,2}:
 4  3
 1  2

julia> sort(A, dims = 1)
2×2 Array{Int64,2}:
 1  2
 4  3

julia> sort(A, dims = 2)
2×2 Array{Int64,2}:
 3  4
 1  2
```
"""
function sort(A::AbstractArray, order::Ordering = Forward;
              dims::Integer,
              alg::Algorithm=DEFAULT_UNSTABLE,
              initialized::Union{Bool,Nothing}=nothing)
    dim = dims
    if initialized !== nothing
        Base.depwarn("`initialized` keyword argument is deprecated", :sort)
    end
    n = length(axes(A, dim))
    if dim != 1
        pdims = (dim, setdiff(1:ndims(A), dim)...)  # put the selected dimension first
        Ap = permutedims(A, pdims)
        Av = vec(Ap)
        sort_chunks!(Av, n, alg, order)
        permutedims(Ap, invperm(pdims))
    else
        Av = A[:]
        sort_chunks!(Av, n, alg, order)
        reshape(Av, axes(A))
    end
end

@noinline function sort_chunks!(Av, n, alg, order)
    inds = LinearIndices(Av)
    for s = first(inds):n:last(inds)
        sort!(Av, s, s+n-1, alg, order)
    end
    Av
end


"""
    sortrows(A, order::Ordering = Forward; alg::Algorithm=DEFAULT_UNSTABLE)

Sort the rows of matrix `A` lexicographically.
See [`sort!`](@ref) for a description of possible
keyword arguments.

# Examples
```jldoctest
julia> sortrows([7 3 5; -1 6 4; 9 -2 8])
3×3 Array{Int64,2}:
 -1   6  4
  7   3  5
  9  -2  8

julia> sortrows([7 3 5; -1 6 4; 9 -2 8], By((x,y) -> isless(x[2],y[2])))
3×3 Array{Int64,2}:
  9  -2  8
  7   3  5
 -1   6  4

julia> sortrows([7 3 5; -1 6 4; 9 -2 8], Backward)
3×3 Array{Int64,2}:
  9  -2  8
  7   3  5
 -1   6  4
```
"""
function sortrows(A::AbstractMatrix, o::Ordering = Forward; kws...)
    inds = axes(A,1)
    T = slicetypeof(A, inds, :)
    rows = similar(A, T, axes(A, 1))
    for i in inds
        rows[i] = view(A, i, :)
    end
    p = sortperm(rows, o; kws...)
    A[p,:]
end

"""
    sortcols(A, order::Ordering = Forward; alg::Algorithm=DEFAULT_UNSTABLE)

Sort the columns of matrix `A` lexicographically.
See [`sort!`](@ref) for a description of possible
keyword arguments.

# Examples
```jldoctest
julia> sortcols([7 3 5; 6 -1 -4; 9 -2 8])
3×3 Array{Int64,2}:
  3   5  7
 -1  -4  6
 -2   8  9

julia> sortcols([7 3 5; 6 -1 -4; 9 -2 8], By((x,y)->isless(x[2],y[2])), alg=InsertionSort)
3×3 Array{Int64,2}:
  5   3  7
 -4  -1  6
  8  -2  9

julia> sortcols([7 3 5; 6 -1 -4; 9 -2 8], Backward)
3×3 Array{Int64,2}:
 7   5   3
 6  -4  -1
 9   8  -2
```
"""
function sortcols(A::AbstractMatrix, order::Ordering = Forward; kws...)
    inds = axes(A,2)
    T = slicetypeof(A, :, inds)
    cols = similar(A, T, axes(A, 2))
    for i in inds
        cols[i] = view(A, :, i)
    end
    p = sortperm(cols, order; kws...)
    A[:,p]
end

function slicetypeof(A::AbstractArray{T}, i1, i2) where T
    I = map(slice_dummy, to_indices(A, (i1, i2)))
    fast = isa(IndexStyle(viewindexing(I), IndexStyle(A)), IndexLinear)
    SubArray{T,1,typeof(A),typeof(I),fast}
end
slice_dummy(S::Slice) = S
slice_dummy(::AbstractUnitRange{T}) where {T} = oneunit(T)

## fast clever sorting for floats ##

module Float
using ..Sort
using ...Order
using ..Base: @inbounds, AbstractVector, Vector, last, axes

import Core.Intrinsics: slt_int
import ..Sort: sort!
import ...Order: DirectOrdering

const Floats = Union{Float32,Float64}

struct Left <: Ordering end
struct Right <: Ordering end

left(::DirectOrdering) = Left()
right(::DirectOrdering) = Right()

left(o::Perm) = Perm(left(o.isless), o.data)
right(o::Perm) = Perm(right(o.isless), o.data)

(o::Left)(x::T, y::T) where {T<:Floats} = slt_int(y, x)
(o::Right)(x::T, y::T) where {T<:Floats} = slt_int(x, y)

isnan(o::DirectOrdering, x::Floats) = (x!=x)
isnan(o::Perm, i::Int) = isnan(o.isless,o.data[i])

function nans2left!(v::AbstractVector, o::Ordering, lo::Int=first(axes(v,1)), hi::Int=last(axes(v,1)))
    i = lo
    @inbounds while i <= hi && isnan(o,v[i])
        i += 1
    end
    j = i + 1
    @inbounds while j <= hi
        if isnan(o,v[j])
            v[i], v[j] = v[j], v[i]
            i += 1
        end
        j += 1
    end
    return i, hi
end
function nans2right!(v::AbstractVector, o::Ordering, lo::Int=first(axes(v,1)), hi::Int=last(axes(v,1)))
    i = hi
    @inbounds while lo <= i && isnan(o,v[i])
        i -= 1
    end
    j = i - 1
    @inbounds while lo <= j
        if isnan(o,v[j])
            v[i], v[j] = v[j], v[i]
            i -= 1
        end
        j -= 1
    end
    return lo, i
end

nans2end!(v::AbstractVector, o::ForwardOrdering) = nans2right!(v,o)
nans2end!(v::AbstractVector, o::BackwardOrdering) = nans2left!(v,o)
nans2end!(v::AbstractVector{Int}, o::Perm{<:ForwardOrdering}) = nans2right!(v,o)
nans2end!(v::AbstractVector{Int}, o::Perm{<:BackwardOrdering}) = nans2left!(v,o)

issignleft(lt::ForwardOrdering, x::Floats) = lt(x, zero(x))
issignleft(lt::BackwardOrdering, x::Floats) = lt(x, -zero(x))
issignleft(o::Perm, i::Int) = issignleft(o.isless, o.data[i])

function fpsort!(v::AbstractVector, a::Algorithm, o::Ordering)
    i, j = lo, hi = nans2end!(v,o)
    @inbounds while true
        while i <= j &&  issignleft(o,v[i]); i += 1; end
        while i <= j && !issignleft(o,v[j]); j -= 1; end
        i <= j || break
        v[i], v[j] = v[j], v[i]
        i += 1; j -= 1
    end
    sort!(v, lo, j,  a, left(o))
    sort!(v, i,  hi, a, right(o))
    return v
end


fpsort!(v::AbstractVector, a::Sort.PartialQuickSort, o::Ordering) =
    sort!(v, first(axes(v,1)), last(axes(v,1)), a, o)

sort!(v::AbstractVector{<:Floats}, a::Algorithm, o::DirectOrdering) = fpsort!(v,a,o)
sort!(v::Vector{Int}, a::Algorithm, o::Perm{<:DirectOrdering,<:Vector{<:Floats}}) = fpsort!(v,a,o)

end # module Sort.Float

end # module Sort

# This file is a part of Julia. License is MIT: https://julialang.org/license

# The tuples and types that do not include 128 bit sizes are necessary to handle
# certain issues on 32-bit machines, and also to simplify promotion rules, as
# they are also used elsewhere where Int128/UInt128 support is separated out,
# such as in hashing2.jl

const BitSigned32_types      = (Int8, Int16, Int32)
const BitUnsigned32_types    = (UInt8, UInt16, UInt32)
const BitInteger32_types     = (BitSigned32_types..., BitUnsigned32_types...)

const BitSigned64_types      = (BitSigned32_types..., Int64)
const BitUnsigned64_types    = (BitUnsigned32_types..., UInt64)
const BitInteger64_types     = (BitSigned64_types..., BitUnsigned64_types...)

const BitSigned_types        = (BitSigned64_types..., Int128)
const BitUnsigned_types      = (BitUnsigned64_types..., UInt128)
const BitInteger_types       = (BitSigned_types..., BitUnsigned_types...)

const BitSignedSmall_types   = Int === Int64 ? ( Int8,  Int16,  Int32) : ( Int8,  Int16)
const BitUnsignedSmall_types = Int === Int64 ? (UInt8, UInt16, UInt32) : (UInt8, UInt16)
const BitIntegerSmall_types  = (BitSignedSmall_types..., BitUnsignedSmall_types...)

const BitSigned32      = Union{BitSigned32_types...}
const BitUnsigned32    = Union{BitUnsigned32_types...}
const BitInteger32     = Union{BitInteger32_types...}

const BitSigned64      = Union{BitSigned64_types...}
const BitUnsigned64    = Union{BitUnsigned64_types...}
const BitInteger64     = Union{BitInteger64_types...}

const BitSigned        = Union{BitSigned_types...}
const BitUnsigned      = Union{BitUnsigned_types...}
const BitInteger       = Union{BitInteger_types...}

const BitSignedSmall   = Union{BitSignedSmall_types...}
const BitUnsignedSmall = Union{BitUnsignedSmall_types...}
const BitIntegerSmall  = Union{BitIntegerSmall_types...}

const BitSigned64T     = Union{Type{Int8}, Type{Int16}, Type{Int32}, Type{Int64}}
const BitUnsigned64T   = Union{Type{UInt8}, Type{UInt16}, Type{UInt32}, Type{UInt64}}

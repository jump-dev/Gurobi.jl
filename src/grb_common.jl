# Common stuff

## convenient types and type conversion functions

typealias GChars Union{Cchar, Char}
typealias IVec Vector{Cint}
typealias FVec Vector{Float64}
typealias CVec Vector{Cchar}

typealias GCharOrVec Union{Cchar, Char, Vector{Cchar}, Vector{Char}}

typealias Bounds{T<:Real} Union{T, Vector{T}}
typealias CoeffMat Union{Matrix{Float64}, SparseMatrixCSC{Float64}}

cchar(c::Cchar) = c
cchar(c::Char) = convert(Cchar, c)

ivec(v::IVec) = v
fvec(v::FVec) = v
cvec(v::CVec) = v

ivec(v::Vector) = convert(IVec, v)
fvec(v::Vector) = convert(FVec, v)
cvec(v::Vector) = convert(CVec, v)

# cvecx(v, n) and fvecx(v, n)
# converts v into a vector of Cchar or Float64 of length n,
# where v can be either a scalar or a vector of length n.

_chklen(v, n::Integer) = (length(v) == n || error("Inconsistent argument dimensions."))

cvecx(c::GChars, n::Integer) = fill(cchar(c), n)
cvecx(c::Vector{Cchar}, n::Integer) = (_chklen(c, n); c)
cvecx(c::Vector{Char}, n::Integer) = (_chklen(c, n); convert(Vector{Cchar}, c))

fvecx(v::Real, n::Integer) = fill(Float64(v), n)
fvecx(v::Vector{Float64}, n::Integer) = (_chklen(v, n); v)
fvecx{T<:Real}(v::Vector{T}, n::Integer) = (_chklen(v, n); convert(Vector{Float64}, v))

# empty vector & matrix (for the purpose of supplying default arguments)

const emptyfvec = Array(Float64, 0)
const emptyfmat = Array(Float64, 0, 0)

# macro to call a Gurobi C function
macro grb_ccall(func, args...)
    f = "GRB$(func)"
    args = map(esc,args)
    is_unix() && return quote
        ccall(($f,libgurobi), $(args...))
    end
    is_windows() && VERSION < v"0.6-" && return quote
        ccall(($f,libgurobi), stdcall, $(args...))
    end
    is_windows() && VERSION >= v"0.6-" && return quote
        ccall(($f,libgurobi), $(esc(:stdcall)), $(args...))
    end
end


# Gurobi library version
function getlibversion()
    _major = Cint[0]
    _minor = Cint[0]
    _tech = Cint[0]
    @grb_ccall(version, Void, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), _major, _minor, _tech)
    return VersionNumber(_major[1], _minor[1], _tech[1])        
end

# version need not be export
# one can write Gurobi.version to get the version numbers
const version = getlibversion()


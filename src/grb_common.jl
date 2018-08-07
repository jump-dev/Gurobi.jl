# Common stuff

## convenient types and type conversion functions

const GChars          = Union{Cchar, Char}
const IVec            = Vector{Cint}
const FVec            = Vector{Float64}
const CVec            = Vector{Cchar}

const GCharOrVec      = Union{Cchar, Char, Vector{Cchar}, Vector{Char}}

const Bounds{T<:Real} = Union{T, Vector{T}}
const CoeffMat        = Union{Matrix{Float64}, SparseMatrixCSC{Float64}}

cchar(c::Cchar) = c
cchar(c::Char) = convert(Cchar, c)

ivec(v::IVec) = v
fvec(v::FVec) = v
cvec(v::CVec) = v

ivec(v::AbstractVector) = convert(IVec, v)
fvec(v::AbstractVector) = convert(FVec, v)
cvec(v::AbstractVector) = convert(CVec, v)

ivec(s::Integer) = Cint[s]

# cvecx(v, n) and fvecx(v, n)
# converts v into a vector of Cchar or Float64 of length n,
# where v can be either a scalar or a vector of length n.

_chklen(v, n::Integer) = (length(v) == n || error("Inconsistent argument dimensions."))

cvecx(c::GChars, n::Integer) = fill(cchar(c), n)
cvecx(c::Vector{Cchar}, n::Integer) = (_chklen(c, n); c)
cvecx(c::Vector{Char}, n::Integer) = (_chklen(c, n); convert(Vector{Cchar}, c))

fvecx(v::Real, n::Integer) = fill(Float64(v), n)
fvecx(v::Vector{Float64}, n::Integer) = (_chklen(v, n); v)
fvecx(v::Vector{T}, n::Integer) where {T<:Real} = (_chklen(v, n); convert(Vector{Float64}, v))

# empty vector & matrix (for the purpose of supplying default arguments)

const emptyfvec = Array{Float64}(undef, 0)
const emptyfmat = Array{Float64}(undef, 0, 0)

# macro to call a Gurobi C function
macro grb_ccall(func, args...)
    f = "GRB$(func)"
    args = map(esc,args)
    if Compat.Sys.isunix()
        return quote
            ccall(($f,libgurobi), $(args...))
        end
    elseif Compat.Sys.iswindows()
        return quote
            ccall(($f,libgurobi), $(esc(:stdcall)), $(args...))
        end
    end
    error("System not recognised.s")
end


# Gurobi library version
function getlibversion()
    _major = Cint[0]
    _minor = Cint[0]
    _tech = Cint[0]
    @grb_ccall(version, Nothing, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}), _major, _minor, _tech)
    return VersionNumber(_major[1], _minor[1], _tech[1])
end

# version need not be export
# one can write Gurobi.version to get the version numbers
const version = getlibversion()

const GRB_INFINITY = 1e100
const GRB_BOUNDMAX = 1e30
function checkvalue(x::Real, bound)
    abs(x) > bound && abs(x) != Inf
end
checkvalue(x::Vector, bound) = any(checkvalue(c, bound) for c in x)

_objwarning(c) = Compat.@warn("Gurobi will silently silently truncate " *
    "objective coefficients >$(GRB_INFINITY) or <-$(GRB_INFINITY). Current " *
    "objective coefficient extrema: $(extrema(c))")

_boundwarning(lb, ub) = Compat.@warn("Gurobi has implicit variable bounds of " *
    "[-1e30, 1e30]. Settings variable bounds outside this can cause " *
    "infeasibility or unboundedness. Current lower bound extrema: " *
    "$(extrema(lb)). Current upper bound extrema: $(extrema(ub))")

# Common stuff

## convenient types and type conversion functions

typealias GChars Union(Cchar, Char)
typealias IVec Vector{Cint}
typealias FVec Vector{Float64}

cchar(c::Cchar) = c
cchar(c::Char) = convert(Cchar, c)

ivec{I<:Integer}(v::Vector{I}) = convert(IVec, v)
fvec{T<:Real}(v::Vector{T}) = convert(FVec, v)

# macro to call a Gurobi C function
macro grb_ccall(func, args...)
    f = "GRB$(func)"
    quote
        ccall(($f,libgurobi), $(args...))
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





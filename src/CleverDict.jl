"""
    CleverDict{K, V}

A smart storage type for managing sequential objects with non-decreasing integer
indices.

Provided no keys are deleted, the backing storage is a `Vector{V}`. Once a key
has been deleted, the backing storage switches to a standard Julia `Dict{K,
V}`.

Use `new_key(c::CleverDict)` to obtain the next key in the sequential order.

When `empty!` is called, the largest index is stored in `offset` so that new
keys are generated from `offset + 1`.

Overload the functions `index_to_key` and `key_to_index` to enable mappings
between the integer index of the vector (plus the offset) and the dictionary key.
"""
mutable struct CleverDict{K, V}
    last_index::Int
    vector::Union{Nothing, Vector{V}}
    offset::Int
    dict::Union{Nothing, Dict{K, V}}
    CleverDict{K, V}() where {K, V} = new{K, V}(0, V[], 0, nothing)
end

"""
    index_to_key(::Type{K}, index::Int)

Create a new key associated with the integer value `index`.
"""
function index_to_key(::Type{K}, index::Int) where {K}
    error("You need to define `index_to_key`.")
end

"""
    key_to_index(key::K)

Map `key` to an integer valued index.
"""
function key_to_index(key::K) where {K}
    error("You need to define `key_to_index`.")
end

"""
    new_key(c::CleverDict{K, V}) where {K, V}

Return the next available key in the ordered sequence.
"""
function new_key(c::CleverDict{K, V}) where {K, V}
    c.last_index += 1
    return index_to_key(K, c.last_index)
end

function Base.empty!(c::CleverDict{K, V}) where {K, V}
    c.vector = V[]
    c.offset = c.last_index
    c.dict = nothing
    return
end

function Base.getindex(c::CleverDict{K, V}, key::K) where {K, V}
    if c.vector !== nothing
        # Case I) no call to `MOI.delete!`, so return the offsetted element:
        return c.vector[key_to_index(key) - c.offset]
    end
    # Case II) `MOI.delete!` must have been called, so return the element
    #          from the dictionary.
    return c.dict[key]
end

function Base.setindex!(c::CleverDict{K, V}, val::V, key::K) where {K, V}
    if c.vector !== nothing
        push!(c.vector, val)
    else
        c.dict[key] = val
    end
    return val
end

function Base.delete!(c::CleverDict{K, V}, key::K) where {K, V}
    if c.vector !== nothing
        c.dict = Dict{K, V}()
        for (i, info) in enumerate(c.vector)
            if key.value == i + c.offset
                continue
            end
            c.dict[MOI.VariableIndex(i + c.offset)] = info
        end
        c.vector = nothing
    else
        delete!(c.dict, key)
    end
    return
end

function Base.length(c::CleverDict)
    if c.vector !== nothing
        return length(c.vector)
    else
        return length(c.dict)
    end
end

function Base.iterate(c::CleverDict)
    if c.vector !== nothing
        next = iterate(c.vector)
        if next === nothing
            return nothing
        else
            (val, new_state) = next
            return (MOI.VariableIndex(1 + c.offset), val), new_state
        end
    else
        return iterate(c.dict)
    end
end

function Base.iterate(c::CleverDict, state)
    if c.vector !== nothing
        next = iterate(c.vector, state)
        if next === nothing
            return nothing
        else
            (val, new_state) = next
            return MOI.VariableIndex(state + c.offset) => val, new_state
        end
    else
        return iterate(c.dict, state)
    end
end

Base.haskey(::CleverDict, key) = false
function Base.haskey(c::CleverDict{K, V}, key::K) where {K, V}
    if c.vector !== nothing
        return 1 <= key_to_index(key) - c.offset <= length(c.vector)
    else
        return haskey(c.dict, key)
    end
end

function Base.values(c::CleverDict)
    return c.vector !== nothing ? c.vector : values(c.dict)
end

function Base.keys(c::CleverDict{MOI.VariableIndex, V}) where {V}
    if c.vector !== nothing
        return MOI.VariableIndex.((1:length(c.vector)) .+ c.offset)
    end
    return keys(c.dict)
end

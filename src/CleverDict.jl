"""
    CleverDict{K, V}

A smart storage type for managing sequential objects with non-decreasing integer
indices.

Provided no keys are deleted, the backing storage is a `Vector{V}`. Once a key
has been deleted, the backing storage switches to a standard Julia `Dict{K,
V}`.

The i'th ordered element can be obtained with `c[LinearIndex(i)]`.

Use `new_key(c::CleverDict)` to obtain the next key in the sequential order.

Overload the functions `index_to_key` and `key_to_index` to enable mappings
between the integer index of the vector and the dictionary key.
"""
mutable struct CleverDict{K, V}
    last_index::Int
    vector::Vector{V}
    dict::Union{Nothing, Dict{K, V}}
    CleverDict{K, V}() where {K, V} = new{K, V}(0, V[], nothing)
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
    empty!(c.vector)
    c.last_index = 0
    c.dict = nothing
    return
end

function Base.getindex(c::CleverDict{K, V}, key::K) where {K, V}
    if c.dict === nothing
        # Case I) no call to `Base.delete!`, so return the element:
        return c.vector[key_to_index(key)]
    end
    # Case II) `Base.delete!` must have been called, so return the element
    #          from the dictionary.
    return c.dict[key]
end

function Base.setindex!(c::CleverDict{K, V}, val::V, key::K) where {K, V}
    push!(c.vector, val)
    if c.dict !== nothing
        c.dict[key] = val
    end
    return val
end

struct LinearIndex
    i::Int
end
Base.getindex(c::CleverDict, index::LinearIndex) = c.vector[index.i]

function Base.delete!(c::CleverDict{K, V}, key::K) where {K, V}
    if c.dict === nothing
        c.dict = Dict{K, V}()
        for (i, info) in enumerate(c.vector)
            if key_to_index(key) == i
                continue
            end
            c.dict[index_to_key(K, i)] = info
        end
        splice!(c.vector, key_to_index(key))
    else
        for (i, v) in enumerate(c.vector)
            if v == key
                splice!(c.vector, i)
                break
            end
        end
        delete!(c.dict, key)
    end
    return
end

Base.length(c::CleverDict) = length(c.vector)

function Base.iterate(c::CleverDict{K, V}) where {K, V}
    if c.dict === nothing
        next = iterate(c.vector)
        if next === nothing
            return nothing
        else
            (val, new_state) = next
            return index_to_key(K, 1) => val, new_state
        end
    else
        return iterate(c.dict)
    end
end

function Base.iterate(c::CleverDict{K, V}, state) where {K, V}
    if c.dict === nothing
        next = iterate(c.vector, state)
        if next === nothing
            return nothing
        else
            (val, new_state) = next
            return index_to_key(K, state) => val, new_state
        end
    else
        return iterate(c.dict, state)
    end
end

Base.haskey(::CleverDict, key) = false
function Base.haskey(c::CleverDict{K, V}, key::K) where {K, V}
    if c.dict === nothing
        return 1 <= key_to_index(key) <= length(c.vector)
    else
        return haskey(c.dict, key)
    end
end

function Base.values(c::CleverDict)
    return c.dict === nothing ? c.vector : values(c.dict)
end

function Base.keys(c::CleverDict{K, V}) where {K, V}
    if c.dict === nothing
        return index_to_key.(K, 1:length(c.vector))
    end
    return keys(c.dict)
end

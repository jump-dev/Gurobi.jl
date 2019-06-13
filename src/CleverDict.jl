import OrderedCollections

"""
    CleverDict{K, V}

A smart storage type for managing sequential objects with non-decreasing integer
indices.

Provided no keys are deleted, the backing storage is a `Vector{V}`. Once a key
has been deleted, the backing storage switches to an `OrderedDict{K, V}`.

The i'th ordered element can be obtained with `c[LinearIndex(i)]`.

Note that querying a `LinearIndex` immediately after deleting a key via
`delete!` is very slow. (It requires a rebuild of an ordered list of variables.)

Use `new_key(c::CleverDict)` to obtain the next key in the sequential order.

Overload the functions `index_to_key` and `key_to_index` to enable mappings
between the integer index of the vector and the dictionary key.
"""
mutable struct CleverDict{K, V}
    last_index::Int
    vector::Union{Nothing, Vector{V}}
    dict::Union{Nothing, OrderedCollections.OrderedDict{K, V}}
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

Map `key` to an integer valued index, assuming that there have been no
deletions.
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
    else
        # Case II) `Base.delete!` must have been called, so return the element
        #          from the dictionary.
        return c.dict[key]
    end
end

function Base.setindex!(c::CleverDict{K, V}, val::V, key::K) where {K, V}
    if c.vector !== nothing
        push!(c.vector, val)
    end
    if c.dict !== nothing
        c.dict[key] = val
    end
    return val
end

struct LinearIndex
    i::Int
end

function Base.getindex(c::CleverDict{K, V}, index::LinearIndex) where {K, V}
    # Get the `index` linear element. If `c.vector` is currently `nothing`
    # (i.e., there has been a deletion), rebuild `c.vector`. This is a
    # trade-off: We could ensure `c.vector` is always updated, but this requires
    # a `splice!` in `delete!`, making deletions costly. However, it makes this
    # `getindex` operation trival because we would never have to rebuild the
    # vector.
    # The current implemented approach offers quick deletions, but an expensive
    # rebuild the first time you query a `LinearIndex` after a deletion. Once
    # the rebuild is done, there are quick queries until the next deletion.
    # Thus, the worst-case is a user repeatedly deleting a key and then
    # querying a LinearIndex (e.g., getting the MOI objective function).
    if c.vector === nothing
        c.vector = Vector{V}(undef, length(c))
        for (i, val) in values(c.dict)
            c.vector[i] = val
        end
    end
    return c.vector[index.i]
end

function Base.delete!(c::CleverDict{K, V}, key::K) where {K, V}
    if c.dict === nothing
        c.dict = OrderedCollections.OrderedDict{K, V}()
        for (i, info) in enumerate(c.vector)
            c.dict[index_to_key(K, i)] = info
        end
    end
    delete!(c.dict, key)
    c.vector = nothing
    return
end

function Base.length(c::CleverDict)
    if c.vector !== nothing
        return length(c.vector)
    else
        return length(c.dict)
    end
end

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

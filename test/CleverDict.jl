using Gurobi, Test

struct MyKey
    x::Int
end
Gurobi.key_to_index(key::MyKey) = key.x
Gurobi.index_to_key(::Type{MyKey}, index::Int) = MyKey(index)

@testset "CleverDict" begin
    @testset "get/set" begin
        d = Gurobi.CleverDict{MyKey, String}()
        key = Gurobi.add_item(d, "first")
        @test key == MyKey(1)
        @test d[key] == "first"
        @test haskey(d, key) == true
        @test_throws KeyError d[MyKey(2)]
        delete!(d, key)
        @test_throws KeyError d[key]
        @test_throws KeyError d[key] = "key"
        @test haskey(d, key) == false
        key2 = Gurobi.add_item(d, "second")
        @test key2 == MyKey(2)
        @test d[key2] == "second"
        @test d.vector === nothing
        @test d.dict !== nothing
        d[key2] = "third"
        @test d[key2] == "third"

        empty!(d)

        key = Gurobi.add_item(d, "first")
        @test key == MyKey(1)
        @test d[key] == "first"
        d[key] = "zeroth"
        @test d[key] == "zeroth"
        @test haskey(d, key) == true
        @test_throws KeyError d[MyKey(2)]
        delete!(d, key)
        @test_throws KeyError d[key]
        @test_throws KeyError d[key] = "key"
        @test haskey(d, key) == false
        key2 = Gurobi.add_item(d, "second")
        @test key2 == MyKey(2)
        @test d[key2] == "second"
        @test d.vector === nothing
        @test d.dict !== nothing
    end

    @testset "LinearIndex" begin
        d = Gurobi.CleverDict{MyKey, String}()
        key = Gurobi.add_item(d, "first")
        @test d[Gurobi.LinearIndex(1)] == "first"
        key2 = Gurobi.add_item(d, "second")
        @test d[Gurobi.LinearIndex(2)] == "second"
        @test length(d) == 2
        delete!(d, key)
        @test d.vector === nothing
        @test d[Gurobi.LinearIndex(1)] == "second"
        @test_throws KeyError d[Gurobi.LinearIndex(2)]
        @test length(d) == 1
        @test d.vector !== nothing
    end

    @testset "keys/values" begin
        d = Gurobi.CleverDict{MyKey, String}()
        key = Gurobi.add_item(d, "first")
        key2 = Gurobi.add_item(d, "second")
        @test collect(keys(d)) == [MyKey(1), MyKey(2)]
        @test collect(values(d)) == ["first", "second"]
        delete!(d, key)
        key3 = Gurobi.add_item(d, "third")
        @test collect(keys(d)) == [MyKey(2), MyKey(3)]
        @test collect(values(d)) == ["second", "third"]
    end

    @testset "iterate" begin
        d = Gurobi.CleverDict{MyKey, String}()
        key = Gurobi.add_item(d, "first")
        key2 = Gurobi.add_item(d, "second")
        my_keys = MyKey[]
        my_values = String[]
        for (k, v) in d
           push!(my_keys, k)
           push!(my_values, v)
        end
        @test my_keys == [MyKey(1), MyKey(2)]
        @test my_values == ["first", "second"]
        delete!(d, key)
        key3 = Gurobi.add_item(d, "third")
        my_keys = MyKey[]
        my_values = String[]
        for (k, v) in d
           push!(my_keys, k)
           push!(my_values, v)
        end
        @test my_keys == [MyKey(2), MyKey(3)]
        @test my_values == ["second", "third"]
    end

    @testset "iterate ii" begin
        d = Gurobi.CleverDict{MyKey, String}()
        key = Gurobi.add_item(d, "first")
        key2 = Gurobi.add_item(d, "second")
        my_keys = MyKey[]
        my_values = String[]
        for (k, v) in d
            push!(my_keys, k)
            push!(my_values, v)
        end
        @test my_keys == [MyKey(1), MyKey(2)]
        @test my_values == ["first", "second"]
        delete!(d, key)
        @test d[Gurobi.LinearIndex(1)] == "second"
        key3 = Gurobi.add_item(d, "third")
        my_keys = MyKey[]
        my_values = String[]
        for (k, v) in d
            push!(my_keys, k)
            push!(my_values, v)
        end
        @test my_keys == [MyKey(2), MyKey(3)]
        @test my_values == ["second", "third"]
    end

    @testset "delete!" begin
        d = Gurobi.CleverDict{MyKey, String}()
        @test length(d) == 0
        @test delete!(d, MyKey(0)) == nothing
        k1 = Gurobi.add_item(d, "a")
        k2 = Gurobi.add_item(d, "b")
        d[Gurobi.LinearIndex(2)] == "b"
        delete!(d, k1)
        d[Gurobi.LinearIndex(1)] == "b"
        k3 = Gurobi.add_item(d, "c")
        @test d[k3] == "c"
        @test d[Gurobi.LinearIndex(1)] == "b"
        @test d[Gurobi.LinearIndex(2)] == "c"
    end
end

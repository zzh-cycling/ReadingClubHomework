include("BP.jl")
using Test
using Graphs

function is_tree(g::AbstractGraph)
    n = nv(g)  # 顶点数
    e = ne(g)  # 边数

    # 边数不等于 n-1 直接返回 false
    e != n - 1 && return false

    # 检查连通性（空图或单顶点图直接视为连通）
    n ≤ 1 && return true

    # 使用 BFS/DFS 判断连通性（避免依赖特定包的实现）
    visited = falses(n)
    queue = [1]  # 从第一个顶点开始遍历
    visited[1] = true
    count = 1

    while !isempty(queue)
        v = popfirst!(queue)
        for u in neighbors(g, v)
            if !visited[u]
                visited[u] = true
                push!(queue, u)
                count += 1
            end
        end
    end

    # 所有顶点都被访问过则为连通
    return count == n
end

clauses1 = SATformula([[-2, -3, -4, 5], [-1, -5, 6], [-5, 7], [-1, -6, -7], [-1, -2, 5], [-1, -3, 5], [-1, -4, 5], [-1, 2, 3, 4, 5, -6]])
FG = FactorGraph(clauses1)




@testset "FactorGraph" begin
    clauses2 = SATformula([[1, -2], [-1, 2, 3], [3, 4, -5], [1]])
    FG=FactorGraph(clauses2)
    # 1 -> 5 var index, 6->9 clause idx, or a, b, c, d
    g = FG.g
    @test g == SimpleGraph{Int64}(9, [[6, 7, 9], [6, 7], [7, 8], [8], [8], [1, 2], [1, 2, 3], [3, 4, 5], [1]])
    @test FG.num_vars == 5
end

function Set_messages(FG::FactorGraph)
    # Set messages in the factor graph
    messages = Dict{Tuple{Int, Int}, Vector{Float64}}()
    g = FG.g
    for e in edges(g)
        messages[(src(e), dst(e))] = [0.5, 0.5]  # Initialize messages to uniform distribution
        messages[(dst(e), src(e))] = [0.5, 0.5]  # Initialize messages in the opposite direction
    end
    
    return messages
end

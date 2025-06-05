using Random
using Statistics
using Graphs, SimpleWeightedGraphs
using LinearAlgebra
using SparseArrays

struct SATformula
    clauses::Vector{Vector{Int}}
    num_vars::Int
    num_clauses::Int
end

struct SATsolution
    solutions::Vector{Symbol}
    num_vars::Int
    clauses::Vector{Vector{Int}}  # Store the clauses for reference
end

# function SATformula(clauses::Vector{Vector{Int}})::SATformula
#     num_vars = maximum(abs.(reduce(vcat, clauses)))
#     num_clauses = length(clauses)
#     return SATformula(clauses, num_vars, num_clauses)
# end

Base.:(==)(s1::SATformula, s2::SATformula) = s1.clauses == s2.clauses && s1.num_vars == s2.num_vars && s1.num_clauses== s2.num_clauses
Base.:(==)(s1::SATsolution, s2::SATsolution) = s1.solutions == s2.solutions && s1.num_vars == s2.num_vars && s1.clauses == s2.clauses

function SATformula(clauses::Vector{T})::SATformula where T 
    if isempty(clauses)
        return SATformula(clauses, 0, 0)
    end
    non_empty_clauses = filter(!isempty, clauses)
    if isempty(non_empty_clauses)
        num_vars = 0
    else
        num_vars = maximum(abs.(reduce(vcat, non_empty_clauses)))
    end
    num_clauses = length(clauses)
    return SATformula(clauses, num_vars, num_clauses)
end

function SATsolution(solutions::Vector{Symbol},formula::SATformula)::SATsolution
    if length(solutions) != formula.num_vars
        throw(ArgumentError("Solutions length must match number of variables."))
    end
    clauses = formula.clauses
    num_vars = formula.num_vars
    return SATsolution(solutions, num_vars, clauses)
end

# factor graph, the first num_vars are the variables (the indices), the rest are the factors (the tensors)
struct FactorGraph{T}
    g::SimpleWeightedGraph{T}
    num_vars::Int
    function FactorGraph(g::SimpleWeightedGraph{T}, num_vars::Int) where T
        for e in edges(g)
            s, d = src(e), dst(e)
            # neighbor of a variable is a factor, and vice versa
            @assert ((s ≤ num_vars) && (d > num_vars)) || ((s > num_vars) && (d ≤ num_vars))
        end
        new{T}(g, num_vars)
    end
end

Base.show(io::IO, fg::FactorGraph) = print(io, "FactorGraph{variables: $(fg.num_vars), factors: $(nv(fg.g) - fg.num_vars)}")
Base.copy(fg::FactorGraph) = FactorGraph(copy(fg.g), fg.num_vars)

Graphs.edges(fg::FactorGraph) = edges(fg.g)
Graphs.vertices(fg::FactorGraph) = vertices(fg.g)
Graphs.nv(fg::FactorGraph) = nv(fg.g)
Graphs.ne(fg::FactorGraph) = ne(fg.g)
Graphs.has_edge(fg::FactorGraph, s, d) = has_edge(fg.g, s, d)
Graphs.has_vertex(fg::FactorGraph, v) = has_vertex(fg.g, v)
Graphs.rem_edge!(fg::FactorGraph, s, d) = rem_edge!(fg.g, s, d)
Graphs.rem_vertex!(fg::FactorGraph, v) = rem_vertex!(fg.g, v)
Graphs.rem_vertices!(fg::FactorGraph, vs) = rem_vertices!(fg.g, vs)
Graphs.add_edge!(fg::FactorGraph, s, d) = add_edge!(fg.g, s, d)
Graphs.add_vertex!(fg::FactorGraph) = add_vertex!(fg.g)
Graphs.add_vertices!(fg::FactorGraph, n) = add_vertices!(fg.g, n)
Graphs.neighbors(fg::FactorGraph, v) = neighbors(fg.g, v)
is_factor(fg::FactorGraph, v) = v > fg.num_vars
is_variable(fg::FactorGraph, v) = v ≤ fg.num_vars

function FactorGraph(code::SATformula)::FactorGraph
    num_vars = code.num_vars
    num_clauses = code.num_clauses
    g = SimpleWeightedGraph(num_vars + num_clauses)
    
    # Add edges between variables and clauses
    for (i, clause) in enumerate(code.clauses)
        clause_idx = num_vars + i  # Clause index in the graph
        for var in clause
            # Add edge from variable to clause
            add_edge!(g, abs(var), clause_idx, sign(var))
        end
    end
    
    return FactorGraph(g, num_vars)
end

function is_variable(v::Int, num_vars::Int)::Bool
    """
    Check if the vertex v is a variable (i.e., its index is less than or equal to num_vars).

    """
    return v<= num_vars
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

function seperate_messages(messages::Dict{Tuple{Int, Int}, TA}, num_vars::Int64) where TA
    """
    Seperate messages into two dictionaries: one for variable to factor messages and one for factor to variable messages.
    """
    messages_v2f = Dict{Tuple{Int, Int}, TA}()
    messages_f2v = Dict{Tuple{Int, Int}, TA}()
    
    for (key, value) in messages
        if is_variable(key[1], num_vars)  # Variable to factor message
            messages_v2f[key] = value
        else  # Factor to variable message
            messages_f2v[key] = value
        end
    end
    
    return messages_v2f, messages_f2v
end

function BP_update!(update_messages::Dict{Tuple{Int, Int}, TA}, num_vars::Int64, edge_types::Dict{Tuple{Int, Int}, Int}, damping_factor::Float64=0.5, tolerance::Float64=1e-6 ,randomvalue::Bool=true) where TA
    # Input: Set of all update_messages arriving onto each variable node j ∈ V(a)\i
    messages_v2f, messages_f2v = seperate_messages(messages, num_vars)
    bond = collect(keys(messages))
    fuse(x,y)=[x...,y...]
    vertices = Set(foldl(fuse, bond))
    len = length(messages)

    for vertex in vertices

    end

    new_messages = Dict{Tuple{Int, Int}, Vector{Float64}}()
    
    V_plus = Dict{Int, Vector{Tuple{Int, Int}}}()  # 变量j的原变量边 (a,j) where J=-1
    V_minus = Dict{Int, Vector{Tuple{Int, Int}}}() # 变量j的否定边 (a,j) where J=1
    
    I, J, V = findnz(edge_types)  # I=行索引, J=列索引, V=元素值

    for idx in eachindex(V)
        i = I[idx]
        j = J[idx]
        val = V[idx]

        # 注意：Julia 稀疏矩阵默认列优先存储
        if val == -1.0
            # 添加到 V_plus (j对应键, (i,j)是坐标)
            arr = get!(() -> Tuple{Int, Int}[], V_plus, j)
            push!(arr, (i, j))
        elseif val == 1.0
            # 添加到 V_minus
            arr = get!(() -> Tuple{Int, Int}[], V_minus, j)
            push!(arr, (i, j))
        end
    end
    
    # 遍历每条消息边 (a -> i)
    for (a, i) in keys(messages)
        # 获取 V(a)\i：所有与a相连的变量节点，排除i
        neighbors = [j for (a_prime, j) in keys(messages) if a_prime == a && j != i]
        isempty(neighbors) && continue  # 跳过无邻居的情况

        γ_product = 1.0  # 初始化γ乘积

        # 对每个邻居j ∈ V(a)\i 计算γ_j→a
        for j in neighbors
            # 获取边(a,j)的类型
            J_aj = get(edge_types, (a, j), 0)
            @assert J_aj in [1, -1] "Edge type for ($a, $j) must be 1 or -1"

            # 确定Vu(j→a)和Vs(j→a)的边集合
            if J_aj == -1
                # J_aj=-1: 原变量边，按公式Vu(j)=V_plus[j]\a，Vs=V_minus[j]\a
                Vu = [b for (b, j_node) in V_plus[j] if b != a]
                Vs = [b for (b, j_node) in V_minus[j] if b != a]
            else
                # J_aj=1: 否定边，按公式Vu(j)=V_minus[j]\a，Vs=V_plus[j]\a
                Vu = [b for (b, j_node) in V_minus[j] if b != a]
                Vs = [b for (b, j_node) in V_plus[j] if b != a]
            end

            # 计算Pu_j→a: ∏_{b ∈ Vu} (1 - δ_b→j)
            Pu = 1.0
            for (b, j_node) in Vu
                # 根据边类型选择消息分量：原变量边取[1]，否定边取[2]
                J_bj = edge_types[(b, j_node)]
                δ = (J_bj == -1) ? messages[(b, j_node)][1] : messages[(b, j_node)][2]
                Pu *= (1 - δ)
            end
            Pu = isempty(Vu) ? 1.0 : Pu  # 处理空集

            # 计算Ps_j→a: ∏_{b ∈ Vs} (1 - δ_b→j)
            Ps = 1.0
            for (b, j_node) in Vs
                J_bj = edge_types[(b, j_node)]
                δ = (J_bj == -1) ? messages[(b, j_node)][1] : messages[(b, j_node)][2]
                Ps *= (1 - δ)
            end
            Ps = isempty(Vs) ? 1.0 : Ps

            # 计算γ_j→a = Ps / (Pu + Ps)
            γ = Ps / (Pu + Ps + eps())  # 加eps()防止除零
            γ_product *= γ
        end

        # 根据边类型(a,i)确定更新哪个消息分量
        J_ai = edge_types[(a, i)]
        new_δ = (J_ai == -1) ? [γ_product, messages[(a, i)][2]] : [messages[(a, i)][1], γ_product]
        
        # 应用阻尼（可选）
        if damping_factor < 1.0
            old_δ = messages[(a, i)]
            new_δ = damping_factor * new_δ + (1 - damping_factor) * old_δ
        end

        new_messages[(a, i)] = new_δ
    end

    # 更新原消息字典
    for (k, v) in new_messages
        messages[k] = v
    end

    return messages
end

function BP(FG::FactorGraph, max_iter::Int=1000, randomvalue::Bool=true,  tol::Float64=1e-6)
    """
    Belief Propagation algorithm for solving SAT problems represented as a FactorGraph.
    Parameters:
        g : FactorGraph representing the SAT problem
        max_iter : maximum number of iterations
        tol : tolerance for convergence
    Returns:
        "UN-CONVERGED" or all messages.
    """
    # Initialize messages
    messages = Set_messages(FG)
    g = FG.g
    # Initialize variable assignments
    
    edge_types = g.weights  # Edge types: -1 for original variable edges, 1 for negated variable edges
    
    
    
    messages_v2f, messages_f2v = seperate_messages(messages, FG.num_vars)
    t = collect(keys(messages_f2v))
    order = randomvalue ? t[sortperm(rand(length(t)))] : t

    
    

    for iter in 1:max_iter
        # Update messages from variables to clauses
        for (a, i) in order
            # 获取 V(a)\i：所有与a相连的变量节点，排除i
            neighbors = [j for (a_prime, j) in keys(messages) if a_prime == a && j != i]
            isempty(neighbors) && continue  # 跳过无邻居的情况
            input_message = Dict((a, k) => messages[(a, k)] for k in neighbors if haskey(messages, (a, k)))
            output_message = BP_update!(input_message, FG.num_vars, edge_types, 0.5, tol, randomvalue)

        end     


        for v in 1:g.num_vars
            neighbors = g.neighbors[v]
            for c in neighbors
                msg = 1.0  # Start with a neutral message
                for other_v in neighbors
                    if other_v != v
                        msg *= messages[(other_v, c)]
                    end
                end
                messages[(v, c)] = msg
            end
        end
        
        # Update messages from clauses to variables
        for c in 1:g.num_clauses
            neighbors = g.clause_neighbors[c]
            for v in neighbors
                msg = 1.0  # Start with a neutral message
                for other_c in g.variable_neighbors[v]
                    if other_c != c
                        msg *= messages[(v, other_c)]
                    end
                end
                messages[(c, v)] = msg
            end
        end
        
        # Check convergence (optional)
        if iter > 1 && maximum(abs.(values(messages) .- prev_messages)) < tol
            break
        end
        
        prev_messages = copy(messages)
    end
    
    # Extract variable assignments from messages (simplified logic)
    for v in 1:g.num_vars
        assignment[v] = if rand() < 0.5 :t else :f end  # Random assignment as placeholder
    end
    return SATsolution(assignment, g.num_vars, g.clauses)
    
end


# function compute_entropy(alpha; k=3, N=100_000, num_iter=100, num_samples=10_000)
#     """
#     BP algorithm for calculating entropy density in random 3-SAT
#     Parameters:
#         alpha : clause density
#         k     : SAT problem order (default 3-SAT)
#         N     : population size
#         num_iter : iterations for population dynamics
#         num_samples : Monte Carlo samples for entropy calculation
#     """
#     # Initialize populations with uniform distributions
#     rng = MersenneTwister()  # Separate RNG for thread safety
#     eta = rand(rng, N)
#     eta_s = rand(rng, N)
    
#     # Population dynamics iterations
#     for _ in 1:num_iter
#         # Update ηs population (clause to variable messages)
#         for i in 1:N
#             selected = rand(rng, 1:N, k-1)
#             new_eta_s = 1.0 - prod(1.0 .- eta[selected])
#             replace_idx = rand(rng, 1:N)
#             eta_s[replace_idx] = new_eta_s
#         end
        
#         # Update η population (variable to clause messages)
#         for i in 1:N
#             p = rand(rng, Poisson(k * alpha / 2))
#             q = rand(rng, Poisson(k * alpha / 2))
#             total = p + q
            
#             if total == 0
#                 new_eta = 0.5
#             else
#                 selected = rand(rng, 1:N, total)
#                 pos_neg = eta_s[selected]
#                 pos = @view pos_neg[1:p]
#                 neg = @view pos_neg[p+1:end]
                
#                 # Stabilize computation using logarithms
#                 log_prod1 = sum(log, pos) + sum(x -> log(1 - x), neg)
#                 log_prod2 = sum(x -> log(1 - x), pos) + sum(log, neg)
#                 max_log = max(log_prod1, log_prod2)
#                 log_total = max_log + log(exp(log_prod1 - max_log) + exp(log_prod2 - max_log))
#                 new_eta = exp(log_prod1 - log_total)
#             end
            
#             replace_idx = rand(rng, 1:N)
#             eta[replace_idx] = new_eta
#         end
#     end
    
#     # Calculate entropy components
#     term1 = zeros(num_samples)
#     term2 = zeros(num_samples)
#     term3 = zeros(num_samples)
    
#     # Parallel computation using threads
#     Threads.@threads for i in 1:num_samples
#         rng = Random.MersenneTwister()  # Thread-local RNG
        
#         # Term 1 calculation
#         indices = rand(rng, 1:N, k)
#         product = prod(eta[indices])
#         term1[i] = log(1 - product + 1e-10)
        
#         # Term 2 calculation
#         p = rand(rng, Poisson(k * alpha / 2))
#         q = rand(rng, Poisson(k * alpha / 2))
#         total = p + q
        
#         if total == 0
#             term2[i] = log(2.0)
#         else
#             selected = rand(rng, 1:N, total)
#             pos_neg = eta_s[selected]
#             pos = @view pos_neg[1:p]
#             neg = @view pos_neg[p+1:end]
            
#             log_prod1 = sum(log, pos) + sum(x -> log(1 - x), neg)
#             log_prod2 = sum(x -> log(1 - x), pos) + sum(log, neg)
#             max_log = max(log_prod1, log_prod2)
#             term2[i] = max_log + log(exp(log_prod1 - max_log) + exp(log_prod2 - max_log))
#         end
        
#         # Term 3 calculation
#         idx_eta = rand(rng, 1:N)
#         idx_eta_s = rand(rng, 1:N)
#         term3[i] = log(eta_s[idx_eta_s] * eta[idx_eta] + 
#                       (1 - eta_s[idx_eta_s]) * (1 - eta[idx_eta]) + 1e-10)
#     end
    
#     # Final entropy calculation
#     s = alpha * mean(term1) + mean(term2) - (k * alpha / 2) * mean(term3)
#     return s
# end

# # Example usage
# let alpha = 4.0  # Clause density
#     entropy = compute_entropy(alpha)
#     println("Entropy density: ", entropy)
# end
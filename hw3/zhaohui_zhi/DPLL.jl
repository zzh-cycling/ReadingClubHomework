# Only used for Conjunctive Normal Form (CNF)
# Do not consider Empty clauses input, but empty clause is valid.
struct SATformula
    clauses::Vector{Vector{Int}}
    num_vars::Int
    num_clauses::Int
end

struct SATsolution
    solutions::Vector{Bool}
    num_vars::Int
    clauses::Vector{Vector{Int}}  # Store the clauses for reference
end

function SATformula(clauses::Vector{Vector{Int}})::SATformula
    num_vars = maximum(abs.(reduce(vcat, clauses)))
    num_clauses = length(clauses)
    return SATformula(clauses, num_vars, num_clauses)
end

function SATformula(clauses::Vector{T})::SATformula where T   
    if isempty(clauses)
        throw(ArgumentError("Clauses cannot be empty."))
    end
    if any(length(c) == 0 for c in clauses)
        throw(ArgumentError("Clauses cannot contain empty vectors."))
    end
    if any(!all(isinteger, c) for c in clauses)
        throw(ArgumentError("All literals in clauses must be integers."))
    end
    num_vars = maximum(abs.(reduce(vcat, clauses)))
    num_clauses = length(clauses)
    return SATformula(clauses, num_vars, num_clauses)
end

function SATsolution(solutions::Vector{Bool},formula::SATformula)::SATsolution
    if length(solutions) != formula.num_vars
        throw(ArgumentError("Solutions length must match number of variables."))
    end
    clauses = formula.clauses
    num_vars = formula.num_vars
    return SATsolution(solutions, num_vars, clauses)
end

function is_unit_clause(clause::Vector{Int}, assignment::Vector{Union{Nothing, Bool}})
    # a unit clause has exactly one literal unassigned, will others not satisfied.
    index = abs.(clause)
    try solution=assignment[index]
        for literal in clause
            if literal > 0 && solution[literal] == true
                return true
            elseif literal < 0 && assignment[-literal] == false
                return true
            end
        end
    catch e
        throw(ArgumentError("Clause has additional literal beyond assignment."))
    end
    return length(clause) == 1
end

function is_empty_clause(clause::Vector{Int})
    return length(clause) == 0
end

function is_satisfied(clause::Vector{Int}, assignment::Vector{Bool})
    # once one literal is true, thus clause is satisfied
    if isempty(clause)
        return true
    end
    for literal in clause
        if literal > 0 && assignment[literal] == true
            return true
        elseif literal < 0 && assignment[-literal] == false
            return true
        end
    end
    return false
end

function is_satisfied_all(clauses::SATformula, assignment::SATsolution)
    clauses_set= clauses.clauses
    solution = assignment.solutions
    for clause in clauses_set
        if !is_satisfied(clause, solution)
            return false
        end
    end
    return true
end

function is_unsatisfied(clauses::SATformula)
    clauses_set = clauses.clauses
    for clause in clauses_set
        if is_empty_clause(clause)
            return true
        end
    end
    return false
end

function unit_propagate(clauses::SATformula, assignment::SATsolution)
    clauses_set = clauses.clauses
    solution = assignment.solutions
    while true
        unit_clauses = filter(is_unit_clause, clauses_set)
        if isempty(unit_clauses)
            break
        end
        for unit_clause in unit_clauses
            literal = unit_clause[1]
            if literal > 0
                solution[literal] = true
            else
                solution[-literal] = false
            end
            clauses = filter(c -> !is_satisfied(c, solution), clauses_set)
        end
    end
    return SATformula(clauses_set), SATsolution(solution, clauses_set)
end

function choose_literal(clauses::SATformula, assignment::SATsolution)
    clauses_set = clauses.clauses
    solution = assignment.solutions
    for clause in clauses
        for literal in clause
            if literal > 0 && assignment[literal] === nothing
                return literal
            elseif literal < 0 && assignment[-literal] === nothing
                return -literal
            end
        end
    end
    return nothing
end

function dpll(clauses::SATformula, assignment::SATsolution)
    clauses, assignment = unit_propagate(clauses, assignment)
    if is_unsatisfied(clauses, assignment)
        return false
    end
    if is_satisfied_all(clauses, assignment)
        return true
    end
    literal = choose_literal(clauses, assignment)
    if literal == nothing
        return true
    end
    assignment[literal] = true
    result = dpll(clauses, assignment)
    if result == true
        return true
    end
    assignment[literal] = false
    return dpll(clauses, assignment)
end

# function solve_sat(formula::SATformula)
#     assignment = fill(nothing, formula.num_vars)
#     result = dpll(formula.clauses, assignment)
#     if result == true
#         return SATsolution(formula.clauses, formula.num_vars)
#     else
#         return nothing
#     end
# end

# function print_solution(solution::SATsolution)
#     if solution == nothing
#         println("No solution found.")
#     else
#         println("Solution found:")
#         for i in 1:solution.num_vars
#             if solution.clauses[i] == true
#                 println("Variable $i: true")
#             elseif solution.clauses[i] == false
#                 println("Variable $i: false")
#             else
#                 println("Variable $i: undefined")
#             end
#         end
#     end
# end
# function main()
#     clauses = [[1, -2], [-1, 2], [2, 3], [-3]]
#     formula = SATformula(clauses)
#     solution = solve_sat(formula)
#     print_solution(solution)
# end
# main()
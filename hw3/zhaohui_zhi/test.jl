using Test
include("DPLL.jl")

@testset "is_satisfied" begin
    # Test case 1: Clause satisfied by positive literal
    clause1 = [1, -2]
    assignment1 = [true, false]
    @test is_satisfied(clause1, assignment1) == true

    # Test case 2: Clause satisfied by negative literal
    clause2 = [1, -2]
    assignment2 = [true, true]
    @test is_satisfied(clause2, assignment2) == true

    # Test case 3: Clause not satisfied
    clause3 = [1, -2]
    assignment3 = [false, true]
    @test is_satisfied(clause3, assignment3) == false

    # Test case 4: Empty clause (should be satisfied)
    clause4 = []
    assignment4 = [true, false]
    @test is_satisfied(clause4, assignment4) == false
end

@testset "unit_propagate" begin
    # Test case 1: Simple unit propagation
    clauses1 = [[1, -2], [-1, 2], [2]]
    assignment1 = [nothing, nothing, nothing]
    new_clauses1, new_assignment1 = unit_propagate(clauses1, assignment1)
    @test new_assignment1 == [true, false, true]
    @test new_clauses1 == []

    # Test case 2: No unit clauses
    clauses2 = [[1, -2], [-3]]
    assignment2 = [nothing, nothing, nothing]
    new_clauses2, new_assignment2 = unit_propagate(clauses2, assignment2)
    @test new_assignment2 == [nothing, nothing, nothing]
    @test new_clauses2 == [[1, -2], [-3]]

    # Test case 3: Empty clauses
    clauses3 = []
    assignment3 = [nothing]
    new_clauses3, new_assignment3 = unit_propagate(clauses3, assignment3)
    @test new_assignment3 == [nothing]
    @test new_clauses3 == []
end


# Test cases for DPLL algorithm
@testset "DPLL Tests" begin
    # Test case 1: Simple satisfiable formula
    formula1 = [[1, -2], [-1, 2], [2, 3]]
    result1 = dpll(formula1)
    @test result1 == true

    # Test case 2: Simple unsatisfiable formula
    formula2 = [[1, -2], [-1, 2], [2, -3]]
    result2 = dpll(formula2)
    @test result2 == false

    # Test case 3: More complex satisfiable formula
    formula3 = [[1, 2], [-1, -3], [3]]
    result3 = dpll(formula3)
    @test result3 == true

    # Test case 4: Empty clause (unsatisfiable)
    formula4 = [[1, -2], [-1, 2], []]
    result4 = dpll(formula4)
    @test result4 == false

    # Test case 5: All variables are true
    formula5 = [[1], [2], [3]]
    result5 = dpll(formula5)
    @test result5 == true

end
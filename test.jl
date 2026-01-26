using Serialization
using CSV
using NamedArrays
using SparseArrays
using DataFrames
using LinearAlgebra


path_input_julia = "./src/test_micro_it2/inputs_julia"

max_flows_N = deserialize("$(path_input_julia)/max_flows_N.jls")



function jl_return_constraints_for_micro_iteration(subproblem_id::Ptr{UInt8})

        println("entered in return_constraints_for_micro_iteration")
        F_N_values = csv_to_Dict("./src/micro_iteration_1_constraints_sub_0_benders_1.csv")
        sub_id = unsafe_strin(subproblem_id) 
        println("sub_id $(sub_id)")

        # dict_results_overflow_N = get_overflows_N(unsafe_string(subproblem_id),F_N_values,vars_dict["inc_by_sub"][unsafe_string(subproblem_id)])

        # dict_results_overflow_N_K = get_overflows_N_K(unsafe_string(subproblem_id), F_N_values)

        # constraints_to_add, N_constraints_micro_it = sort_results_and_return_constraints(dict_results_overflow_N,dict_results_overflow_N_K)
        # append!(vars_dict["inc_by_sub"][unsafe_string(subproblem_id)],N_constraints_micro_it)
        println("got to the end !!!!!")
end 


jl_return_constraints_for_micro_iteration("sub/sub_0.mps")



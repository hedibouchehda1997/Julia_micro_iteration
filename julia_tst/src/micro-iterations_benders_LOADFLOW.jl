using Serialization
using CSV
using DataFrames

################################################################
###################### HELPER FUNCTIONS ########################
################################################################
function csv_to_Dict(path::AbstractString)
    df = CSV.read(path, DataFrame, header = false,stringtype=String)
    @assert ncol(df) ≥ 2 "Le CSV doit contenir au moins deux colonnes"
    return Dict(df[!, 1] .=> df[!, 2])
end

############################################
######### INPUTS FROM STUDY FILE ###########
############################################
""" These inputs are loaded from the study subfile inputs_julia.
They can be loaded in C++ code and converted into Julia objects, so they are not re-read at each Benders iteration."""

path_input_julia = "mps_files_debug/test_micro_it0/inputs_julia" 

dict_incident_outage_AC_branches = deserialize("$(path_input_julia)/dict_incident_outage_AC_branches.jls")
dict_incident_HVDC_branches = deserialize("$(path_input_julia)/dict_incident_HVDC_branches.jls")
max_flows_N = deserialize("$(path_input_julia)/max_flows_N.jls")
max_flows_N_K = deserialize("$(path_input_julia)/max_flows_N_K.jls")


############################################
######### INPUTS FROM PTDF SCRIPT ###########
############################################
""" Inputs obtained at the beginning of each benders iteration with the script micro-iterations_benders_PTDF.jl.
Here they have been stored as julia objects but they should be stored as C++ objects."""

path_output_julia = "mps_files_debug/test_micro_it0/output_julia_PTDF"

all_monitored_branches = deserialize("$(path_output_julia)/all_monitored_branches.jls")
dict_incident_factors = deserialize("$(path_output_julia)/dict_incident_factors.jls")
HVDC_new_dict = deserialize("$(path_output_julia)/HVDC_new_dict.jls")

################################################################
######################## FUNCTIONS #############################
################################################################

function get_overflows_N(max_flows_N, all_monitored_branches,v,F_N_values,N_constraints_added)
    """ This function estimates the overflows in N situation a regarding the monitored line + the branches invested in.

    Arguments : 
    - max_flows_N -> dictionnary of maximum flow on monitored branch (branch,v) => value
    - all_monitored_branches -> vector of monitored branches, including monitored branch from the initial network + invested branch at this Benders iteration
    - v -> variant = hour of the subproblem
    - F_N_values -> flow values obtained by solving the subproblem
    - N_constraints_added -> list of N constraints previously added to the subproblem

    Return :
    - dict_results_overflow_N -> dictionnary with the overflows in N
    """


    #Initialize overflows dictionnaries
    dict_results_overflow_N = Dict{String, Float64}() #Branche monitorée

    # Compute N overflows
    for monitored in setdiff(all_monitored_branches, N_constraints_added) # on ne regarde pas les branches déjà rajouter en N précédemment
        overflow = max((abs(F_N_values[monitored])) - abs(max_flows_N[monitored,v]),0)
        if overflow > 0 
            dict_results_overflow_N[monitored] = overflow
        end
    end

    return dict_results_overflow_N
end

function get_overflows_N_K(max_flows_N_K, all_monitored_branches, dict_incident_outage_AC_branches, dict_incident_HVDC_branches, dict_incident_factors, HVDC_new_dict, v, F_N_values)
    """ This function estimates the overflows in N-k situation, regarding the monitored line + the branches invested in.
    To compute the flow through a monitored branch during an incident involving AC branches and HVDC branches, we use incident_factors and HVDC_sensi factors :
    F_N_K(monitored) = F_N(monitored) + incident_factor(monitored,incident,outage_AC_branch) * F_N(outage_AC_branch) + HVDC_sensi(outage_HVDC_branch,monitored)

    Arguments : 
    - max_flows_N_K -> dictionnary of maximum flow on monitored branch (branch,v) => value
    - all_monitored_branches -> vector of monitored branches, including monitored branch from the initial network + invested branch at this Benders iteration
    - dict_incident_outage_AC_branches -> dictionnary of AC branches involved in each incident
    - dict_incident_HVDC_branches -> dictionnary of HVDC branches involved in each incident
    - dict_incident_factors-> dictionary of incident factors computed for the Benders iteration topology: (monitored_branch,incident,AC_branch)
    - HVDC_new_dict -> dictionnary with HVDC sensitivity factors computed for the Benders iteration topology : (HVDC_branch,monitored_branch)
    - F_N_values -> flow values obtained by solving the subproblem

    Return :
    - dict_results_overflow_N_K -> dictionnary with the overflows in N-k
    """

    #Initialize overflows dictionnaries
    dict_results_overflow_N_K= Dict{Tuple{String, String}, Float64}() # (branche monitorée,incident)

    incidents_vector = keys(dict_incident_outage_AC_branches)

    # Compute N-K overflows : iterative process on each incident
    for incident in incidents_vector

        ## Retrieve incident information # 
        outage_AC_branches = dict_incident_outage_AC_branches[incident]
        outage_hvdcs       = dict_incident_HVDC_branches[incident]

        # Compute N-K overflows for monitored lines + candidates invested in
        for monitored in all_monitored_branches
            if !(monitored in outage_AC_branches)
                F_N_k_loadflow = F_N_values[monitored] + sum((F_N_values[outage_AC_branch] * dict_incident_factors[monitored,incident,outage_AC_branch] for outage_AC_branch in outage_AC_branches);init=0.0) - sum((F_N_values[outage_hvdc] * HVDC_new_dict[outage_hvdc, monitored] for outage_hvdc in outage_hvdcs); init= 0.0)
                overflow = max((abs(F_N_k_loadflow) - abs(max_flows_N_K[monitored, v])),0) 
                if overflow > 0 
                    dict_results_overflow_N_K[monitored,incident] = overflow
                end 
            end
        end
    end

    dict_results_overflow_N_K = sort(dict_results_overflow_N_K, by=x->x[2], rev=true)
    return dict_results_overflow_N_K
end

function sort_results_and_return_constraints(dict_results_overflow_N,dict_results_overflow_N_K)
    """ This function sorts the results given by the two dictionnaries of overflows and keeps the most violated constraints. 

    Arguments : 
    - dict_results_overflow_N -> dictionnary with the overflows in N
    - dict_results_overflow_N_K -> dictionnary with the overflows in N-k
    

    Return :
    - max_threat_N_K -> Dictionnary which gives the incident that causes the most important overflow for each monitored branch
    - N_K_overflow_in_N -> Dictionnary which indicates if a branch with overflow in a N-k situation has also overflow in N situation (so that N situation wil be addressed first)
    - lignes_N -> List of branches with overflow in N situation
    """
    
    


    constraints_to_add = [] # Vector which will be sent to C++
    N_constraints_micro_it = [] # Vector used to track N_constraints added during this micro-iteration
    N_K_constraints_micro_it = [] # Vector used to track N_K_constraints added during this micro-iteration

    
    for monitored in keys(dict_results_overflow_N)
        append!(constraints_to_add,["branch_$(monitored)"])
        append!(N_constraints_micro_it,[monitored])
    end

    for overflow in keys(dict_results_overflow_N_K)
        monitored = overflow[1]  #overflow de la forme (monitored,incident) => value
        if !(monitored in N_constraints_micro_it) & !(monitored in N_K_constraints_micro_it) # on ajoute d'abord la contrainte en N et seulement un incident par lignes surveillées
            incident = overflow[2] #on ajoute le plus gros incident (max threat)
            append!(constraints_to_add,["branch_$(monitored)_inc_$(incident)"])
            append!(constraints_to_add,["inc_$(incident)"])
            append!(N_K_constraints_micro_it,[monitored])
        end
    end
    
    return constraints_to_add, N_constraints_micro_it 
end

################################################################
######################## INPUTS FROM BENDERS ###################
################################################################
""" Here the dictionnary is loaded from a csv but it will be sent by Benders in a CPP object."""
F_N_values = csv_to_Dict("inputs_benders/micro_iteration_1_constraints_sub_0_benders_1.csv")

#1ere itération : initialize N_constraints_added
N_constraints_added = []

# other iterations :
# N_constraints_added charged

# name of subproblem : 
v = "sub_0.mps" # ou sub/sub_0.mps ?

#####################################################################################
##################################### COMPUTATION ###################################
#####################################################################################

function return_constraints_for_micro_iteration(max_flows_N,max_flows_N_K,dict_incident_outage_AC_branches,dict_incident_HVDC_branches, # from study inputs
                                               all_monitored_branches,dict_incident_factors,HVDC_new_dict, # from benders iteration - PTDF script
                                                v, F_N_values,N_constraints_added) # from benders subproblem
    dict_results_overflow_N = get_overflows_N(max_flows_N,all_monitored_branches,v,F_N_values,N_constraints_added)

    dict_results_overflow_N_K = get_overflows_N_K(max_flows_N_K, all_monitored_branches, dict_incident_outage_AC_branches, dict_incident_HVDC_branches, dict_incident_factors, HVDC_new_dict, v, F_N_values)

    constraints_to_add, N_constraints_micro_it = sort_results_and_return_constraints(dict_results_overflow_N,dict_results_overflow_N_K)

    append!(N_constraints_added,N_constraints_micro_it)

    return constraints_to_add, N_constraints_added

end

constraints_to_add, N_constraints_added = return_constraints_for_micro_iteration(max_flows_N,max_flows_N_K,dict_incident_outage_AC_branches,dict_incident_HVDC_branches, all_monitored_branches,dict_incident_factors,HVDC_new_dict, v, F_N_values,N_constraints_added) 
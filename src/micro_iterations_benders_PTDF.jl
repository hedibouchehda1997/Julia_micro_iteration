using Serialization
using CSV
using NamedArrays
using SparseArrays
using DataFrames
using LinearAlgebra

################################################################
###################### HELPER FUNCTIONS ########################
################################################################
"""Used to open csv file and change type"""
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

path_input_julia = "./test_micro_it2/inputs_julia" 

B_inv = deserialize("$(path_input_julia)/B_inv.jls")
Ab = deserialize("$(path_input_julia)/Ab.jls")
Yl = deserialize("$(path_input_julia)/Yl.jls")
A_hvdc = deserialize("$(path_input_julia)/A_hvdc.jls")
branches_to_candidates_dict = deserialize("$(path_input_julia)/branches_to_candidates_dict.jls")
dict_incident_outage_AC_branches = deserialize("$(path_input_julia)/dict_incident_outage_AC_branches.jls")
max_flows_N = deserialize("$(path_input_julia)/max_flows_N.jls")
max_flows_N_K = deserialize("$(path_input_julia)/max_flows_N_K.jls")
n_side1_dict = deserialize("$(path_input_julia)/n_side1_dict.jls")
n_side2_dict = deserialize("$(path_input_julia)/n_side2_dict.jls")

################################################################
######################## FUNCTIONS #############################
################################################################

function get_invested_branches(z_dict::Dict, branches_to_candidates_dict::Dict)
    """ This fonction gives the list of invested and not invested branches at the previous master iteration, based on the binary value of z of the associated candidates. 
    !!! One candidate can represent two branches, in the case of LAdouble for example. That's why we don't use directly the id of candidates. 
    Arguments : z_dict : type to determine (charged from CPP)
    Returns : branches_invested (Vector) and branches_not_invested (Vector)"""
    branches_invested = collect(k for (k,v) in branches_to_candidates_dict if z_dict[v]== 1)   
    branches_not_invested = collect(k for (k,v) in branches_to_candidates_dict if z_dict[v]== 0)
    
    return branches_invested, branches_not_invested
end

function update_PTDF_after_lines_removal(B_inv, Ab, Yl, branches_not_invested)
    """This function computes the new PTDF from the PTDF of the complete network, thanks to the Woodbury formula. 
    Arguments :
    - B_inv : inverse of admittance matrix of the complete network (NamedMatrix)
    - Ab : Incidence matrix of the complete network (Named Sparse Matrix)
    - Yl : Diagonal matrix of admittance of the complete network (Named Sparse Matrix)
    - branches_not_invested = vector of lines not invested

    Returns :
    - PTDF_new : PTDF of the network with branches not invested removed (NamedMatrix)
    """
    k = length(branches_not_invested) # number of deleted lines on the network 
    n = size(B_inv, 1) # number of nodes without slack nodes


    # Construire A_k = [√b_l * a_l_sans_slack]
    A_k = zeros(n, k)
    for (j, br) in enumerate(branches_not_invested)
        # vecteur d’incidence complet
        a = Ab[br, :] |> collect
        # suppression du slack
        y_l = Yl[br, br]
        A_k[:, j] = sqrt(y_l) * a
    end

    # Calcul de la matrice intermédiaire M qu'il faut ensuite inverser
    M = I - A_k' * (B_inv * A_k)
    M_inv = inv(Matrix(M))  # si îlot possible : utiliser pinv(M)

    # Mise à jour de B_inv
    BinvA = B_inv * A_k
    Bnew_inv = B_inv + BinvA * M_inv * (BinvA')

    # Mettre à jour Ab et Yl (suppression des lignes retirées)
    Ab_new = Ab[setdiff(names(Ab, 1), branches_not_invested), :]
    Yl_new = Yl[setdiff(names(Yl, 1), branches_not_invested), setdiff(names(Yl, 2), branches_not_invested)]

    # Nouvelle PTDF
    PTDF_new = Yl_new * Ab_new * Bnew_inv


    return PTDF_new
end

function update_HVDC_sensi_after_lines_removal(A_hvdc, PTDF_new)
    """This function computes the new HVDC sensitivity factors. 
    Arguments :
    - A_hvdc (Named Sparse Matrix)
    - PTDF_new : PTDF Matrix (NamedMatrix)

    Returns :
    - HVDC_sensitivity_matrix_new : HVDC sensitivity factors of the network with branches not invested removed (NamedMatrix)
    """
    HVDC_sensitivity_matrix_new = PTDF_new * A_hvdc

    return HVDC_sensitivity_matrix_new
end

function compute_sensis_after_lines_removal(B_inv, Ab, Yl, A_hvdc, branches_not_invested)
    """This function computes sensitivity factors for the new network. 
    Arguments :
    - B_inv : inverse of admittance matrix of the complete network
    - Ab : Incidence matrix of the complete network
    - Yl : Diagonal matrix of admittance of the complete network
    - branches_not_invested = vector of lines not invested
    - all_monitored_branches = 

    Returns :
    - PTDF_new : PTDF of the network with branches not invested removed (NamedMatrix)
    - HVDC_sensi_new : HVDC sensitivity factors of the network with branches not invested removed (NamedMatrix)
    - PTDF_new_dict : dictionnary with PTDF_new values
    - HVDC_sensi_new_dict : dictionnary with HVDC_sensi_new values

    """
    # Compute matrix (PTDF and HVDC sensis)
    PTDF_new = update_PTDF_after_lines_removal(B_inv, Ab, Yl, branches_not_invested)
    HVDC_sensi_new = update_HVDC_sensi_after_lines_removal(A_hvdc,PTDF_new)

    # Pruning 
    threshold_ptdf = 1e-5
    PTDF_new[abs.(PTDF_new) .< threshold_ptdf] .= 0.

    threshold_hvdc_sensi = 1e-5
    HVDC_sensi_new[abs.(HVDC_sensi_new) .< threshold_hvdc_sensi] .= 0.0
   
    return PTDF_new, HVDC_sensi_new
end

function helper_convert_sensitivity_array_to_dict(named_array::NamedArray{Float64, 2})
    """This helper function takes as argument a named array, and converts it to a dictionnary of type Dict{Tuple{eltype(col_labels), eltype(row_labels)}, eltype(named_array)} , 
    where the first element of the tuple key is the row, and the second element is the column.
    Returns : the dictionnary containing the NamedArray data"""

    row_labels = names(named_array, 1)
    col_labels = names(named_array, 2)
    sensi_dict = Dict{Tuple{eltype(col_labels), eltype(row_labels)}, eltype(named_array)}(
        (col_labels[j], row_labels[i]) => named_array[i, j]
        for i in eachindex(row_labels), j in eachindex(col_labels)
    )

    return sensi_dict
end

function compute_incident_factors_after_lines_removal(all_monitored_branches, dict_incident_outage_AC_branches, branches_invested, PTDF_new, n_side1_dict, n_side2_dict)
    """This function computes the incident factors, used to get flow in N-K situation for the new network :
    if the incident i corresponds to X AC outages o, the flow in line monitored m is F_N_k[m,i] = F_N[m] + sum(F_N[o] * incident_factor[m,i,o] for o in i)

    Arguments :
    - all_monitored_branches : Vector of monitored branches
    - dict_incident_outage_AC_branches : 
    - PTDF_new : PTDF Matrix
    - n_side1_dict : dictionnary which gives the origin node of a branch
    - n_side2_dict : dictionnary which gives the extremity node of a branch

    Returns :
    - dict_incident_factors : dictionnary with incident values : [monitored,incident,outage] => incident_factor 
    """
    ## PRE TREATMENTS ##

    #Initialise a dictionnary which will contain all the FPM matrices
    #It will contain NamedArrays, or nothing if the incident contains a branch breaking synchronicity
    dict_incident_factors = Dict{Tuple{String,String,String}, Float64}()

    PTDF_matrix = Array(PTDF_new)

    PTDF_new_dict = helper_convert_sensitivity_array_to_dict(PTDF_new[intersect(names(PTDF_new, 1), all_monitored_branches), :])

    ## MAIN LOOP ##

    #Iterate over all incidents to compute the incident factors

    for incident_id in keys(dict_incident_outage_AC_branches)

        #Get the involved AC branches
        outage_AC_branches = dict_incident_outage_AC_branches[incident_id]

        #Get the origin (side 1) and extremity (side 2) nodes of the AC outage branches
        side1_nodes = [n_side1_dict[branch] for branch in outage_AC_branches]
        side2_nodes = [n_side2_dict[branch] for branch in outage_AC_branches]
        
        #Extract the part of the PTDF corresponding to the side 1 nodes and outaged branches
        side1_nodes_indexes = map(node -> findfirst(isequal(node), names(PTDF_new, 2)), side1_nodes)
        outage_AC_branches_indexes = map(branch -> findfirst(isequal(branch), names(PTDF_new, 1)), outage_AC_branches)
        PTDF_side1_nodes_outage_lines = PTDF_matrix[outage_AC_branches_indexes, side1_nodes_indexes]
        
        #Extract the part of the PTDF corresponding to the side 2 nodes and outaged branches
        side2_nodes_indexes = map(node -> findfirst(isequal(node), names(PTDF_new, 2)), side2_nodes)
        PTDF_side2_nodes_outage_lines = PTDF_matrix[outage_AC_branches_indexes, side2_nodes_indexes]

        
        #Create corresponding identity matrix (used below to compute the FPM)
        nb_outage_AC_branches = length(outage_AC_branches)
        I_matrix = Matrix{Float64}(I, nb_outage_AC_branches, nb_outage_AC_branches)
        
        #Compute the FPM (fictitious powers matrix), such that fictitious_powers_vector = FPM * flow_outage_lines

        #We first compute FPM^-1 , which is the matrix such that flow_outage_lines = FPM^-1 * fictitious_powers_vector
        FPM_inv_matrix = I_matrix - (Matrix(PTDF_side1_nodes_outage_lines) - Matrix(PTDF_side2_nodes_outage_lines))
        FPM_matrix = inv(FPM_inv_matrix)
        FPM = NamedArray(FPM_matrix, (outage_AC_branches, outage_AC_branches), ("branch_fictitious_power","outage_branch_id"))

        for monitored_line in all_monitored_branches, outage_line in outage_AC_branches
            incident_factor = 0

            for outage_line_bis in outage_AC_branches
                node1_outage = n_side1_dict[outage_line_bis]
                node2_outage = n_side2_dict[outage_line_bis]

                PTDF_node1 = PTDF_new_dict[node1_outage,monitored_line]
                PTDF_node2 = PTDF_new_dict[node2_outage,monitored_line]

                incident_factor += FPM[outage_line_bis,outage_line]*(PTDF_node1 - PTDF_node2)
            end

            dict_incident_factors[monitored_line,incident_id,outage_line] = incident_factor
        end
    end

    # Data filtering 
    threshold_factors = 1e-5

    dict_incident_factors = Dict(
        k => (abs(v) >= threshold_factors ? v : 0)
        for (k, v) in dict_incident_factors
    )


    #Do not return anything (in place modification)
    return dict_incident_factors

end



################################################################
######################## INPUTS FROM BENDERS #############################
################################################################
""" Here the dictionnary is loaded from a csv but it will be sent by Benders in a CPP object."""
z_dict = csv_to_Dict("./master_iter_1.csv") 

#####################################################################################
##################################### COMPUTATION ###################################
#####################################################################################
function compute_factors_for_microiterations(B_inv,Ab,Yl,A_hvdc,n_side1_dict,n_side2_dict,
                                             max_flows_N, max_flows_N_K, dict_incident_outage_AC_branches, branches_to_candidates_dict,
                                            z_dict)

    # Get vectors of invested and not invested branches
    branches_invested,branches_not_invested = get_invested_branches(z_dict,branches_to_candidates_dict)

    # Get vector of monitored lines and incidents lines
    monitored_N = unique(first.(keys(max_flows_N)))
    monitored_N_K = unique(first.(keys(max_flows_N_K)))
    all_monitored_branches = setdiff(monitored_N ∪ monitored_N_K,branches_not_invested)

    branches_with_incidents = unique(vcat(values(dict_incident_outage_AC_branches)...))

    # Compute PTDF and HVDC sensi matrix
    PTDF_new , HVDC_new = compute_sensis_after_lines_removal(B_inv, Ab, Yl, A_hvdc, branches_not_invested)

    # Write HVDC matrix into a dictionnary
    branches_to_keep_in_sensi_dicts = collect(Set([all_monitored_branches; branches_with_incidents])) #to get rid of duplicates
    HVDC_new_dict = helper_convert_sensitivity_array_to_dict(HVDC_new[intersect(names(HVDC_new, 1), branches_to_keep_in_sensi_dicts), :])

    # Compute incident_factors 
    dict_incident_factors = compute_incident_factors_after_lines_removal(all_monitored_branches, dict_incident_outage_AC_branches, branches_invested, PTDF_new, n_side1_dict, n_side2_dict)

    return HVDC_new_dict, dict_incident_factors, all_monitored_branches
end

HVDC_new_dict, dict_incident_factors, all_monitored_branches = compute_factors_for_microiterations(B_inv,Ab,Yl,A_hvdc,n_side1_dict,n_side2_dict,
                                            max_flows_N, max_flows_N_K, dict_incident_outage_AC_branches, branches_to_candidates_dict,
                                            z_dict)


df_1 = DataFrame(
    key1 = [k[1] for k in keys(HVDC_new_dict)],
    key2 = [k[2] for k in keys(HVDC_new_dict)],
    value = collect(values(HVDC_new_dict))
)
CSV.write("./HVDC_new_dict_1.csv",df_1)
df_2 = DataFrame(
        key1 = [k[1] for k in keys(dict_incident_factors)],
        key2 = [k[2] for k in keys(dict_incident_factors)],
        key3 = [k[3] for k in keys(dict_incident_factors)],
        value = collect(values(dict_incident_factors))
    )  
CSV.write("./dict_incident_factors_1.csv",df_2)
df_3 = DataFrame(value=all_monitored_branches)
CSV.write("./all_monitored_branches_1.csv",df_3) 


#####################################################################################
##################################### LOAD AND SEND RESULTS TO JULIA LOADFLOW #######
#####################################################################################

# The objets HVDC_new_dict, dict_incident_factors and all_monitored_branches should be sent to the second julia script : micro_iterations_benders_LOADFLOW.jl

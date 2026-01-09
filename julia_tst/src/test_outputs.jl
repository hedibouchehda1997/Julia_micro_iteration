

using CSV 
    using DataFrames

all_monitored_branches_csv = CSV.read("./all_monitored_branches.csv",DataFrame)
all_monitored_branches_1_csv = CSV.read("./all_monitored_branches_1.csv",DataFrame)

if (isequal(sort(all_monitored_branches_csv), sort(all_monitored_branches_1_csv)))  
    println("all_monitored_branches_csv is good !! ")

end 

dict_inicident_factors = CSV.read("./dict_incident_factors.csv",DataFrame)
dict_inicident_factors_1 = CSV.read("./dict_incident_factors_1.csv",DataFrame)
if (isequal(sort(dict_inicident_factors), sort(dict_inicident_factors_1)))  
    println("dict_inicident_factors is good !! ")

end 

HVDC_new_dict = CSV.read("./HVDC_new_dict.csv",DataFrame)
HVDC_new_dict_1 = CSV.read("./HVDC_new_dict_1.csv",DataFrame)
if (isequal(sort(HVDC_new_dict), sort(HVDC_new_dict_1)))  
    println("HVDC_new_dict_1 is good !! ")

end 

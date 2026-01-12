
include("MyLib.jl")


subs = ["sub_0.mps","sub_1.mps","sub_2.mps"]
subs_ptr = pointer.(subs)
ptr = pointer(subs_ptr)
sub_ids = MyLib.SubProblemsIds(ptr,Cint(3)) 

MyLib.jl_load_variables(sub_ids)




res_1 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-107-108_0"),Cint(0))
res_2 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-108-203"),Cint(0))
res_3 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-115-116_0"),Cint(0))
res_4 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-116-117_0"),Cint(0))
res_5 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-116-118"),Cint(0))
res_6 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-116-121"),Cint(0))
res_7 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-128-203_0"),Cint(0))
res_8 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-216-218"),Cint(0))
res_9 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-216-219_0"),Cint(0))
res_10 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-217-219"),Cint(0))
res_11 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-316-318"),Cint(0))
res_12 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-330-213"),Cint(0))
res_13 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-330-315"),Cint(0))
res_14 = MyLib.CandidateLineMasterIterationResult(pointer("LINE-330-319"),Cint(0))

println("test zebbbi ")

candidates = [res_1,res_2,res_3,res_4,res_5,res_6,res_7,res_8,res_9,res_10,res_11,res_12,res_13,res_14] 
ptr = pointer(candidates) 

master_benders_input = MyLib.MasterBendersInput(ptr,Cint(14))

dict = MyLib.build_dict_from_master_iter_result(master_benders_input)

println("printing the built dictionary") 
println(dict)

z_dict = MyLib.csv_to_Dict("./src/master_iter_1.csv")


if (z_dict == dict) 
    println("build dict from struct is all good !! ")
end 


MyLib.jl_compute_factors_for_microiterations(master_benders_input)
str = pointer("sub_0.mps")

flow_1 = MyLib.FlowN(pointer("LINE-103-109"),Cdouble(49.3576)) 
flow_2 = MyLib.FlowN(pointer("LINE-107-108,"),Cdouble(-29.0993))
flow_3 = MyLib.FlowN(pointer("LINE-107-108_0"),Cdouble(-29.0993))

flows =[flow_1, flow_2, flow_3] 
flows_ptr = pointer(flows)

flows_struct = MyLib.FlowNList(flows_ptr,Cint(3))

MyLib.jl_return_constraints_for_micro_iteration(str,flows_struct)

### Testing compute factors function 

#include <iostream>
#include <vector>


extern "C" {
    #include "libmylib/include/julia_init.h"
    #include "MyLib.h"
}



int main() {

    init_julia(0, NULL);

    char* sub_0 = "sub_0.mps" ; 
    char* sub_1 = "sub_1.mps" ; 
    char* sub_2 = "sub_2.mps" ; 

    char** subs = (char**) malloc(3*sizeof(char*)) ;
    subs[0] = sub_0 ; 
    subs[1] = sub_1 ; 
    subs[2] = sub_2 ; 
    
    auto sub_probs_ids = SubProblemIds{subs,3} ; 

    jl_load_variables(sub_probs_ids) ; 

    char* line_1 = "LINE-107-108_0" ; 
    auto res_1 = CandidateLineMasterIterationResult{line_1,0} ; 

    char* line_2 = "LINE-108-203" ; 
    auto res_2 = CandidateLineMasterIterationResult{line_2,0} ; 
    
    char* line_3 = "LINE-115-116_0" ; 
    auto res_3 = CandidateLineMasterIterationResult{line_3,0} ; 
    
    char* line_4 = "LINE-116-117_0" ; 
    auto res_4 = CandidateLineMasterIterationResult({line_4,0}) ; 
    
    char* line_5 = "LINE-116-118" ; 
    auto res_5 = CandidateLineMasterIterationResult({line_5,0}) ; 
    
    char* line_6 = "LINE-116-121" ; 
    auto res_6 = CandidateLineMasterIterationResult({line_6,0}) ; 

    char* line_7 = "LINE-128-203_0" ; 
    auto res_7 = CandidateLineMasterIterationResult({line_7,0}) ;

    char* line_8 = "LINE-216-218" ; 
    auto res_8 = CandidateLineMasterIterationResult({line_8,0}) ;
    
    char* line_9 = "LINE-216-219_0" ; 
    auto res_9 = CandidateLineMasterIterationResult({line_9,0}) ;

    char* line_10 = "LINE-217-219" ; 
    auto res_10 = CandidateLineMasterIterationResult({line_10,0}) ;
    
    char* line_11 = "LINE-316-318" ; 
    auto res_11 = CandidateLineMasterIterationResult({line_11,0}) ;

    char* line_12 = "LINE-330-213" ; 
    auto res_12 = CandidateLineMasterIterationResult({line_12,0}) ;

    char* line_13 = "LINE-330-315" ; 
    auto res_13 = CandidateLineMasterIterationResult({line_13,0}) ;

    char* line_14 = "LINE-330-319" ; 
    auto res_14 = CandidateLineMasterIterationResult({line_14,0}) ;


    CandidateLineMasterIterationResult* candidates_res = (CandidateLineMasterIterationResult*) malloc(sizeof(CandidateLineMasterIterationResult) *15) ; 
    candidates_res[0] = res_1 ; 
    candidates_res[1] = res_2 ; 
    candidates_res[2] = res_3 ; 
    candidates_res[3] = res_4 ; 
    candidates_res[4] = res_5 ; 
    candidates_res[5] = res_6 ; 
    candidates_res[6] = res_7 ; 
    candidates_res[7] = res_8 ; 
    candidates_res[8] = res_9 ; 
    candidates_res[9] = res_10 ; 
    candidates_res[10] = res_11 ; 
    candidates_res[11] = res_12 ; 
    candidates_res[12] = res_13 ; 
    candidates_res[13] = res_14 ; 

    auto master_benders_input = MasterBendersInput{candidates_res,14} ; 
    
    
   

    jl_compute_factors_for_microiterations(master_benders_input) ;

    std::cout<<"ended jl_compute_factors_for_microiterations"<<std::endl; 

    char* sub = "sub_0.mps" ; 
    
    char* flow_line_1 = "LINE-103-109" ; 
    auto flow_1 = FlowN{flow_line_1,49.3576} ; 
    
    char* flow_line_2 = "LINE-107-108" ; 
    auto flow_2 = FlowN{flow_line_2,-29.0993} ;

    char* flow_line_3 = "LINE-107-108_0" ; 
    auto flow_3 = FlowN{flow_line_3,-29.0993} ; 

    FlowN* flows = (FlowN*) malloc(sizeof(FlowN)*3) ; 
    flows[0] = flow_1 ; 
    flows[1] = flow_2 ;
    flows[2] = flow_3 ; 

    auto flow_list = FlowNList{flows,3} ; 

    jl_return_constraints_for_micro_iteration(sub,flow_list) ; 

    jl_test() ; 

    shutdown_julia(0);
    return 0;
}
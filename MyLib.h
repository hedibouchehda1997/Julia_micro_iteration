#ifndef MYLIB_H
#define MYLIB_H

#ifdef __cplusplus
extern "C" {
#endif

struct CandidateLineMasterIterationResult 
{
    char* candidate_line_id ; 
    int is_invested ; 
}; 

struct MasterBendersInput
{
    CandidateLineMasterIterationResult* candidates_res; 
    int size; 
} ; 

struct SubProblemIds 
{
    char** subProblems_ids ; 
    int n_subproblems ; 
} ; 

struct FlowN
{
    char* line_id ; 
    double values ; 
}; 

struct FlowNList
{
    FlowN* flows ; 
    int size;
};

void jl_load_variables(SubProblemIds) ; 
void jl_test() ; 
void jl_compute_factors_for_microiterations(MasterBendersInput) ; 
void jl_return_constraints_for_micro_iteration(char*,FlowNList) ; 

#ifdef __cplusplus
}
#endif

#endif // MYLIB_H
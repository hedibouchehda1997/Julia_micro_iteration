#include <dlfcn.h>
#include <iostream>


extern "C" {
    #include "libmylib/include/julia_init.h"
    #include "MyLib.h"
}

using init_julia_FUNC = void (*) (int, char*); 
using shut_down_julia_FUNC = void(*)(int) ; 
using jl_test_FUNC = void(*)() ; 
using jl_load_variables_FUNC = void(*) (SubProblemIds) ; 


int main() 
{
    void* handle = dlopen("./libmylib/lib/libmylib.so",RTLD_NOW) ; 

    if (handle) 
    {
        
        std::cout<<"handle to the so file is not null "<<std::endl ;
        
        init_julia_FUNC init_julia = (init_julia_FUNC) dlsym(handle,"init_julia") ; 
        init_julia(0,NULL) ;
        
        jl_test_FUNC jl_test = (jl_test_FUNC) dlsym(handle,"jl_test") ; 
        jl_test() ; 

        int n_subproblems = 4 ; 
        char** subproblems_ids = (char**) malloc(n_subproblems * sizeof(char*)) ; 
        
        subproblems_ids[0] = "sub_0.mps" ; 
        subproblems_ids[1] = "sub_1.mps" ; 
        subproblems_ids[2] = "sub_2.mps" ; 
        subproblems_ids[3] = "sub_3.mps" ; 

        SubProblemIds sub_ids = SubProblemIds{subproblems_ids, n_subproblems} ;  
        jl_load_variables_FUNC jl_load_variables = (jl_load_variables_FUNC) dlsym(handle,"jl_load_variables") ; 
        jl_load_variables(sub_ids) ; 


        shut_down_julia_FUNC shut_down_julia = (shut_down_julia_FUNC) dlsym(handle,"shutdown_julia") ; 
        shut_down_julia(0) ;


    }
}
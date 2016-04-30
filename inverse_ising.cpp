#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "lbfgs.h"
#include "dynamic_inverse_ising.hpp"

static lbfgsfloatval_t evaluate(void *instance,
				const double * x,
				double *g,
				const int n,
				const lbfgsfloatval_t step) {
    DynamicInverseIsing * dynamic_inverse_ising = 
	static_cast<DynamicInverseIsing *>(instance);
    double SL = dynamic_inverse_ising->evaluate(n, x, g);
    printf("%.10f, %.10f, %.10f \n", SL, g[n-1], x[n-1]);
    //getchar();
    return SL;
}

static int progress(void *instance,
		    const double *s,
		    const double *g,
		    const double SL,
		    const lbfgsfloatval_t xnorm,
		    const lbfgsfloatval_t gnorm,
		    const lbfgsfloatval_t step,
		    int n,
		    int k,
		    int ls) {
    //if (k%10 == 0) {
	printf("Iteration %d:  ",k);
	printf("Object function = %16.15f  ", SL);
	printf(" = %16.15f  step = %16.15f\n", gnorm, step);
	//}
    return 0;
}

int main()
{
    lbfgs_parameter_t param;
    double * x;
    double SL;

    DynamicInverseIsing dynamic_inverse_ising;
    dynamic_inverse_ising.init();
    dynamic_inverse_ising.getEvidence();

    int N = dynamic_inverse_ising.L + 2;
    x = new double[N];
    if(x==NULL)
    {
	std::cout<<"Allocating storage for WeightMatrix FAILED!"<< "\n";
	return -1;
    }
    //dynamic_inverse_ising.make_initial(x);
    for (int i=0; i<N; i++)
    {
        x[i] = 0.0;
    }

    lbfgs_parameter_init(&param);
    param.m = 10;
    //param.epsilon = 1e-5;
    param.max_iterations = 20000;
    param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_WOLFE;
  
    //int status = lbfgs(L+1,ss,&SL,evaluate,progress,NULL,&param);
    int status = lbfgs(N,x,&SL,evaluate,progress,&dynamic_inverse_ising,&param);
    printf("L-BFGS optimization terminated with status code = %d, lambda=%f\n",status, x[N-2]);

    dynamic_inverse_ising.update_using_x(x);
    dynamic_inverse_ising.output_result(x);
    for (int i=0; i<N; i++)
    {
        printf("%f \n", x[i]);
    }
    return 0;
}	
	


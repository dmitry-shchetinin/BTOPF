/*=================================================================
 * This function solves LP with one linear constraint by greedy algorithm:
 * max c'x
 * s.t. xl<= x <=xu
 *      bl<=a'x<=bu
 * output info shows: 0 - constraint not binding, 1 - lower bound of
 * constraint is binding, 2 - upper bound of constraint is binding, 3 - problem
 * is infeasible
 *=================================================================*/


#include "mex.h"
#include "math.h"
#include "matrix.h"

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{
    //initialize variables
    int *ind_change, *ind_new, info, i, index, num_change;
    mxLogical *ind_1, *ind_2;
    double *c, *a, *xl, *xu, *xu_out, *ratios, b_change, obj, 
            sum_a, eps, bl, bu, b_viol, a_temp, *new_order;
    mwSize n;
    mxArray *array1, *array2, *array3, *array4, *array5[1], *array6, *lhs1[2];
    
    //get pointers to input variables
    c=mxGetPr(prhs[0]);
    a=mxGetPr(prhs[1]);
    xl=mxGetPr(prhs[2]);
    xu=mxGetPr(prhs[3]);
    bl=mxGetScalar(prhs[4]);
    bu=mxGetScalar(prhs[5]);
    eps=1e-14;
    
    //get number of variables
    n=mxGetNumberOfElements(prhs[0]);  
    
    //initialize intermediate arrays
    array1=mxCreateNumericMatrix(n, 1, mxDOUBLE_CLASS, mxREAL);
    xu_out=mxGetPr(array1);
    array2=mxCreateLogicalMatrix(n, 1);
    ind_1=mxGetLogicals(array2);
    array3=mxCreateLogicalMatrix(n, 1);
    ind_2=mxGetLogicals(array3);
    
    //initialize scalar values
    obj=0;
    sum_a=0;
    b_change=0;
    num_change=0;
    
    //loop over all variables and record necessary values
    for (i=0; i<n; i++) {
        //change upper bound s.t. xl=0
        xu_out[i]=xu[i]-xl[i];
        //record corresponding change on constraint bounds
        b_change=b_change+a[i]*xl[i];
        //record contribution of xl to objective function
        obj=obj+c[i]*xl[i];
        
        //deal with entries with positive c
        if (c[i]>0) {
            if (a[i]>eps)
                ind_1[i]=1;
            else if (a[i]<-eps)
                ind_2[i]=1;
            //record contribution to objective and constraint
            obj=obj+c[i]*xu_out[i];
            sum_a=sum_a+a[i]*xu_out[i];
        }
    }  
    
    //update constraint bounds
    bl=bl-b_change;
    bu=bu-b_change;
    
    //check if a'x constraint is binding at the optimum solution
    if (sum_a>=bl && sum_a<=bu) {
        //record optimal solution
        plhs[0]=mxCreateDoubleScalar(obj);
        plhs[1]=mxCreateDoubleScalar(0);
        //free memory
        mxDestroyArray(array1);
        mxDestroyArray(array2);
        mxDestroyArray(array3);
        return;
    }
    else if (sum_a>bu) { //the value is too big
        //we need to either remove some x with a>0,c>0 or add some x with a<0,c<=0 until a'x=bu
        for (i=0; i<n; i++) {
            if (ind_1[i] || (c[i]<=0 && a[i]<-eps)) {
                ind_1[i]=1;
                num_change++;
            }
            else
                ind_1[i]=0;
        }
        //compute the violation of the constraint
        b_viol=sum_a-bu;
        info=2;
    }
    else { //the value is too small
        //we need to either remove some x with a<0,c>0 or add some x with a>0,c<=0 until a'x=bl
        for (i=0; i<n; i++) {
            if (ind_2[i] || (c[i]<=0 && a[i]>eps)) {
                ind_1[i]=1;
                num_change++;
            }
            else
                ind_1[i]=0;
        }
        //compute the violation of the constraint
        b_viol=bl-sum_a;
        info=1;
    }
    
    //check if the problem is infeasible
    if (num_change==0)
    {
        plhs[0]=mxCreateDoubleScalar(-1e20);
        plhs[1]=mxCreateDoubleScalar(3);
        //free memory
        mxDestroyArray(array1);
        mxDestroyArray(array2);
        mxDestroyArray(array3);
        return;
    } 
    
    //extract elements which can help make problem feasible and their |c/a| ratio
    array4=mxCreateNumericMatrix(num_change, 1, mxINT32_CLASS, mxREAL);
    ind_new=(int *)mxGetData(array4);
    array5[0]=mxCreateNumericMatrix(num_change, 1, mxDOUBLE_CLASS, mxREAL);
    ratios=mxGetPr(array5[0]);
    index=0;
    for (i=0; i<n; i++) {
        if (ind_1[i]) {
            ind_new[index]=i;
            ratios[index]=fabs(c[i]/a[i]);
            index++;
        }
    } 
    
    //sort ratios in increasing order
    mexCallMATLAB(2, lhs1, 1, array5, "sort");
    new_order=mxGetPr(lhs1[1]);
    
    //reorder indices correspondingly
    array6=mxCreateNumericMatrix(num_change, 1, mxINT32_CLASS, mxREAL);
    ind_change=(int *)mxGetData(array6);
    for (i=0; i<num_change; i++)
        ind_change[i]=ind_new[(int)new_order[i]-1];  
 
    //add and/or remove elements by a greedy strategy until the constraint is satisfied
    for (i=0; i<num_change; i++)
    {
        index=ind_change[i];
        //check if we need to remove the entire x entry
        a_temp=fabs(a[index]*xu_out[index]);
        if (b_viol>=a_temp)
        {
            b_viol=b_viol-a_temp;
            obj=obj-fabs(c[index]*xu_out[index]);
        }
        else
        {
            obj=obj-fabs(c[index]*b_viol/a[index]);
            b_viol=0;
            break;
        }
    } 
            
    //record output
    if (fabs(b_viol)>1e-6)
        plhs[0]=mxCreateDoubleScalar(-1e20);
    else
        plhs[0]=mxCreateDoubleScalar(obj);
    plhs[1]=mxCreateDoubleScalar(info);

    //free memory
    mxDestroyArray(array1);
    mxDestroyArray(array2);
    mxDestroyArray(array3);
    mxDestroyArray(array4);
    mxDestroyArray(array5[0]);
    mxDestroyArray(array6);
    mxDestroyArray(lhs1[0]);
    mxDestroyArray(lhs1[1]);
}
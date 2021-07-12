/*   Copyright 2005-2015 The MathWorks, Inc. */
/*   Written by Peter Lindskog. */

/* Include libraries. */
#include "mex.h"
#include <math.h>

/* Specify the number of outputs here. */
#define NY 1

/* State equations. */
void compute_dx(double *dx, double t, double *x, double *u, double **p,
                const mxArray *auxvar)
{
    // Declare model parameters
    double *a1, *a2, *a3, *Cm1, *Rm1, *R0dm1, *R0cm1;
    // Declare intermediate variables
    double u_l = u;
    double u_oc;
    double Q_e;
    double i_b;
    
    // Retrieve model parameters
    a1 = p[0]; a2 = p[1]; a3 = p[2]; // Parameters for the U_oc(SoC)
    Cm1 = p[3]; Rm1 = p[4]; R0dm1 = p[5]; R0cm1 = p[6];
    
    
    /* Declare model parameters. */
    double *Rs, *Cp, *g_ch,*g_pk, *lam1, *lam2, *lam3, *c1, *c2, *c3, *c4, *c5, *c6, *v_r, *tau_x, *g_i;
    /* Declare intermediate variables. */
    double i;   /* injected current.    */     
    double q;   /* injected charge.     */
    double alpha_g;
    double beta_g;
    double z;
    double g_x;
    double g_p_inf;
    double y;   /* output (voltage over skin) */

    
    i = u[0];
    q = u[1];
        
    /* Retrieve model parameters. */
    Rs = p[0];   
    Cp  = p[1];   
    g_ch = p[2];   
    g_pk  = p[3];   
    lam1 = p[4];   
    lam2 = p[5];   
    lam3 = p[6];
    c1 = p[7];
    c2 = p[8];
    c3 = p[9];
    c4 = p[10];
    c5 = p[11];
    c6 = p[12];
    v_r = p[13];
    tau_x = p[14];
    g_i = p[15];
    
    /* Calculate voltage over skin */
    y = x[0] + Rs[0]*i;
    
    if (fabs(i)<0.1e-3 && fabs(y)<0.05){
        // if system at rest (both current and voltage 0: Rp does not change)
        dx[1] = 0;
    }else{
        /* Intermediate equations */
        alpha_g = c2[0]*(1-exp(c3[0]*(y+v_r[0])));
        beta_g = c4[0]/(exp(c5[0]*(y+v_r[0]+c6[0]))+1);
        z = beta_g + alpha_g*(c1[0]-beta_g);
        g_x = (g_ch[0]+i*g_i[0]-q)/tau_x[0];
        g_p_inf = g_pk[0]*(lam1[0]+lam2[0]*(1-exp(lam3[0]*i)))/(g_pk[0]-(lam1[0]+lam2[0]*(1-exp(lam3[0]*i))));
        dx[1] = (g_p_inf-x[1])*z+g_x;
    }
    
    
    /* x[0]: Capacitor Voltage u_cp. */
    /* x[1]: Parallel Conductance g_p. */
    dx[0] = (-x[1]*x[0]+i)/Cp[0];
    
//       printf("Completed Calculation of dx @ t = %d\n",t[0]);
}

/* Output equation. */
void compute_y(double *y, double t, double *x, double *u, double **p,
               const mxArray *auxvar)
{
    double *Rs;
    Rs = p[0];   

    /* y[0]: Skin Voltage. */
    y[0] = x[0]+Rs[0]*u[0];
    
//     printf("Completed Calculation of y @ t = %d\n",t);
}



/*----------------------------------------------------------------------- *
   DO NOT MODIFY THE CODE BELOW UNLESS YOU NEED TO PASS ADDITIONAL
   INFORMATION TO COMPUTE_DX AND COMPUTE_Y
 
   To add extra arguments to compute_dx and compute_y (e.g., size
   information), modify the definitions above and calls below.
 *-----------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Declaration of input and output arguments. */
    double *x, *u, **p, *dx, *y, *t;
    int     i, np;
    size_t  nu, nx;
    const mxArray *auxvar = NULL; /* Cell array of additional data. */
    
    if (nrhs < 3) {
        mexErrMsgIdAndTxt("IDNLGREY:ODE_FILE:InvalidSyntax",
        "At least 3 inputs expected (t, u, x).");
    }
    
    /* Determine if auxiliary variables were passed as last input.  */
    if ((nrhs > 3) && (mxIsCell(prhs[nrhs-1]))) {
        /* Auxiliary variables were passed as input. */
        auxvar = prhs[nrhs-1];
        np = nrhs - 4; /* Number of parameters (could be 0). */
    } else {
        /* Auxiliary variables were not passed. */
        np = nrhs - 3; /* Number of parameters. */
    }
    
    /* Determine number of inputs and states. */
    nx = mxGetNumberOfElements(prhs[1]); /* Number of states. */
    nu = mxGetNumberOfElements(prhs[2]); /* Number of inputs. */
    
    /* Obtain double data pointers from mxArrays. */
    t = mxGetPr(prhs[0]);  /* Current time value (scalar). */
    x = mxGetPr(prhs[1]);  /* States at time t. */
    u = mxGetPr(prhs[2]);  /* Inputs at time t. */
    
    p = mxCalloc(np, sizeof(double*));
    for (i = 0; i < np; i++) {
        p[i] = mxGetPr(prhs[3+i]); /* Parameter arrays. */
    }
    
    /* Create matrix for the return arguments. */
    plhs[0] = mxCreateDoubleMatrix(nx, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(NY, 1, mxREAL);
    dx      = mxGetPr(plhs[0]); /* State derivative values. */
    y       = mxGetPr(plhs[1]); /* Output values. */
    
    /*
      Call the state and output update functions.
      
      Note: You may also pass other inputs that you might need,
      such as number of states (nx) and number of parameters (np).
      You may also omit unused inputs (such as auxvar).
      
      For example, you may want to use orders nx and nu, but not time (t)
      or auxiliary data (auxvar). You may write these functions as:
          compute_dx(dx, nx, nu, x, u, p);
          compute_y(y, nx, nu, x, u, p);
    */
    
    /* Call function for state derivative update. */
    compute_dx(dx, t[0], x, u, p, auxvar);
    
    /* Call function for output update. */
    compute_y(y, t[0], x, u, p, auxvar);
    
    /* Clean up. */
    mxFree(p);
}
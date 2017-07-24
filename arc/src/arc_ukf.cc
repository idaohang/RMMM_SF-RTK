/*********************************************************************************
 *  ARC-SRTK - Single Frequency RTK Pisitioning Library
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  Created on: July 19, 2017
 *      Author: SuJingLan
 *********************************************************************************/

/**
 * @brief ARC-SRTK UKF Functions
 * @author sujinglan
 */
#include <rtklib.h>
#include "arc.h"

#define MAXSTATES 100
/* Macro for defining an exception----------------------------------------*/
ARC_DEFINE_EXCEPTION(Exception, std::runtime_error);

/*-------------------------------------------------------------------------
 * Computes cholesky decomposition of A
 * @param A: a n x n matrix. On exit, the lower triangle
 * of A is overwritten with the cholesky decomposition
 * @param n: the dimension of the square matrix A
 * @param sigma: the vector of diagonal elements of the cholesky
 * decomposition
 *------------------------------------------------------------------------*/
static void arc_ukf_cholesky_decomposition(double *A, unsigned n, double *sigma)
{
    arc_log(ARC_INFO,"arc_ukf_cholesky_decomposition : \n");

    unsigned i,j,k;
    double t;

    for (i=0;i<n;i++) {
        if (i>0) {
            for(j=i;j<n;j++) {
                t=0.0;
                for(k=0;k<i;k++)
                    t+=A[j*n+k]*A[i*n+k];
                A[j*n+i]-=t;
            }
        }
        if(A[i*n+i]<=0.0) {
            return;
        }
        sigma[i]=sqrt(A[i*n+i]);
        t=1.0/sigma[i];
        for(j=i;j<n;j++) A[j*n+i]*=t;
    }
}
/* --------------------------------------------------------------------------
 * Solve the linear system AX^t=B^t given the cholesky decomposition of A
 * @param A: the cholesky decomposition of A, as given by the
 * routine ukf_cholesky_decomposition
 * @param n: the size of the problem
 * @param sigma: the vector of diagonal elements of A
 * @param B: the right hand side of the linear system. Size n x m
 * @param m: the number of simultaneou systems to solve
 * @param X: a matrix holding the result
 * --------------------------------------------------------------------------*/
static void arc_ukf_cholesky_solve(double *A, unsigned n, double *sigma, double *B,
                                   unsigned m, double *X)
{
    arc_log(ARC_INFO,"arc_ukf_cholesky_solve : \n");

    int i,j,k;
    double t;

    for (i =0;i<m;i++) { /* iterate over the lines of B */
        for (j=0;j<n;j++) { /* solve Ly=B */
            t=B[i*n+j];
            for (k=j-1;k>=0;k--)
                t-=A[j*n+k]*X[i*n+k];
            X[i*n+j]=t/sigma[j];
        }

        for (j=n-1;j>=0;j--) { /* solve Ltx=y */
            t=X[i*n+j];
            for (k=j+1;k<n;k++)
                t-=A[k*n+j]*X[i*n+k];
            X[i*n+j]=t/sigma[j];
        }
    }
}
/* --------------------------------------------------------------------
 * Lower triangle of the inverse of the cholesky decomposition.
 * On exit, the lower triangle of A is overwritten with the lower
 * triangle of the inverse.
 * @param A: the cholesky decomposition of A, as given by the
 * routine ukf_cholesky_decomposition
 * @param n: the size of the problem
 * @param sigma: the vector of diagonal elements of A
 *---------------------------------------------------------------------*/
static void arc_ukf_cholesky_invert(double *A, unsigned n, double *sigma)
{
    arc_log(ARC_INFO,"arc_ukf_cholesky_invert : \n");

    double t;
    int i,j,k;

    for (i=0;i<n;i++) {
        A[i*n+i]=1.0/sigma[i];
        for (j=i+1;j<n;j++) {
            t=0.0;
            for (k=i;k<j;k++)
                t-=A[j*n+k]*A[k*n+i];
            A[j*n+i]=t/sigma[j];
        }
    }
}
/* -----------------------------------------------------------------
 * creates a new unscented kalman filter structure
 * @param state_dim: size of state variable
 * @param mesaure_dim: size of measurement
 * @param Q: additive model noise covariance matrix
 * @param R: additive measurement noise covariance matrix
 * @param ffun: the funcion describing the filter evolution equation
 * @param mfun: the measurement function
 * ----------------------------------------------------------------*/
extern ukf_t* arc_ukf_filter_new(unsigned int state_dim,
                                 unsigned int measure_dim,
                                 double *Q,
                                 double *R,
                                 filter_function ffun,
                                 measure_function mfun)
{
    arc_log(ARC_INFO,"arc_ukf_filter_new : \n");

    ukf_t *filter;
    int Size;
    unsigned err=0;
    /* nothing to do if no state or measurement !*/
    if(state_dim==0||measure_dim==0) return NULL;
    /* alloc new structure */
    if (!(filter=(ukf_t*)malloc(sizeof(ukf_t)))) return NULL;

    filter->state_dim=state_dim;
    filter->measure_dim=measure_dim;
    filter->ffun=ffun;
    filter->mfun=mfun;

    filter->x=(double*)malloc(state_dim*sizeof(double));
    err|=(filter->x==NULL);

    filter->y=(double*)malloc(measure_dim*sizeof(double));
    err|=(filter->y==NULL);

    Size=state_dim*state_dim;

    filter->P=(double*)malloc(Size*sizeof(double));
    err|=(filter->P==NULL);

    Size=2*state_dim+1; /* 2n+1 */

    filter->wm=(double*)malloc(Size*sizeof(double));
    err|=(filter->wm==NULL);

    filter->wc=(double*)malloc(Size*sizeof(double));
    err|=(filter->wc==NULL);

    filter->sigma_point=(double*)malloc(Size*state_dim*sizeof(double));
    err|=(filter->sigma_point==NULL);

    Size=filter->state_dim;

    filter->sigma=(double*)malloc(Size*sizeof(double));
    err|=(filter->sigma==NULL);

    filter->PM=(double*)malloc(Size*Size*sizeof(double));
    err|=(filter->PM==NULL);

    filter->PM_save=(double*)malloc(Size*Size*sizeof(double));
    err|=(filter->PM==NULL);

    filter->xm=(double*)malloc(Size*sizeof(double));
    err|=(filter->xm==NULL);

    filter->ym=(double*)malloc(filter->measure_dim *sizeof(double));
    err|=(filter->ym==NULL);

    Size=2*filter->state_dim+1;
    filter->khi=(double*)malloc(2*Size*filter->state_dim*sizeof(double));
    err|=(filter->khi==NULL);

    filter->khi_y=(double*)malloc(Size*filter->measure_dim*sizeof(double));
    err|=(filter->khi_y==NULL);

    Size=filter->measure_dim;
    filter->Pyy=(double*)malloc(Size*Size*sizeof(double));
    err|=(filter->Pyy==NULL);

    filter->Pxy=(double*)malloc(Size*filter->state_dim*sizeof(double));
    err|=(filter->Pxy==NULL);

    filter->dx=(double*)malloc(filter->state_dim*sizeof(double));
    err|=(filter->dx==NULL);

    filter->dy=(double*)malloc(filter->measure_dim*sizeof(double));
    err|=(filter->dy==NULL);

    filter->gain=(double*)malloc(filter->state_dim*filter->measure_dim*sizeof(double));
    err|=(filter->gain==NULL);

    filter->sigma_y=(double*)malloc(filter->measure_dim*sizeof(double));
    err|=(filter->sigma_y==NULL);

    filter->KL=(double*)malloc(filter->state_dim*filter->measure_dim*sizeof(double));
    err|=(filter->KL==NULL);

    if(err!=0) return 0;
    
    Size=filter->state_dim;
    filter->Q=(double*)malloc(Size*Size*sizeof(double));
    if(filter->Q==NULL) {
        return 0;
    }
    memcpy(filter->Q,Q,Size*Size*sizeof(double));

    Size=filter->measure_dim;
    filter->R=(double*)malloc(Size*Size*sizeof(double));
    if(filter->R==NULL) {
        return 0;
    }
    memcpy(filter->R,R,Size*Size*sizeof(double));
    return filter;
}
/* -------------------------------------------------------------------
 * set filter weight using default procedure
 * @param alpha: spread parameter
 * @param k: scaling parameter
 * @param beta: distribution fitting parameter (2 for Gaussian)
 * -----------------------------------------------------------------*/
extern void arc_ukf_filter_compute_weights(ukf_t  *filter,double alpha,
                                           double ZCount,double beta)
{
    arc_log(ARC_INFO,"arc_ukf_filter_compute_weights : \n");

    double l;
    double lam;
    unsigned YCount;

    if(!filter) return;
    l=(double)filter->state_dim;  /* states numbers */

    /* lambda parameter */
    lam=alpha*alpha*(l+ZCount)-l;

    filter->wm[0]=lam/(lam+l);
    filter->wc[0]=filter->wm[0]+(1.0-alpha*alpha+beta);
    for(YCount=1;YCount<=2*filter->state_dim;YCount++) {
        filter->wm[YCount]=0.5/(lam+l);
        filter->wc[YCount]=0.5/(lam+l);
    }
    filter->gamma=alpha*sqrt(l+ZCount);
}
/* --------------------------------------------------------------------
 * Reset filter to new state
 * @param x0: initial state
 * @param P0: initial state covariance matrix
 * -------------------------------------------------------------------*/
extern void arc_ukf_filter_reset(ukf_t* filter,double *x0,double *P0)
{
    arc_log(ARC_INFO,"arc_ukf_filter_reset : \n");

    if(filter) {
        /* state of the filter */
        if (x0) arc_matcpy(filter->x,x0,filter->state_dim,1);
        if (P0) arc_matcpy(filter->P,P0,filter->state_dim,filter->state_dim);
    }
}
/* --------------------------------------------------------------------
 * Get filter state
 * @param x: A vector that will hold the state
 * @param P: A vector that will hold the state covariance
 *--------------------------------------------------------------------*/
extern void arc_ukf_filter_get_state(ukf_t *filter, double *x, double* P)
{
    arc_log(ARC_INFO,"arc_ukf_filter_get_state : \n");

    if(filter) {
        if (x) arc_matcpy(x,filter->x,filter->state_dim,1);
        if (P) arc_matcpy(P,filter->P,filter->state_dim,filter->state_dim);
    }
}
/* --------------------------------------------------------------------
 * Update filter using a measure
 * @param y: The measure vector
 * @param u: the command
 * -------------------------------------------------------------------*/
extern int arc_ukf_filter_update(ukf_t *filter, double *y, double *u,
                                 double*F,double *G)
{
    int l=filter->state_dim; /* numbers of states */
    int m=filter->measure_dim; /* numbers of measures */
    int i,j,k;
    double t;
    /* propagate measurements and gotten measurements max difference for check */
    static const double MAXDY=5.0;

    arc_log(ARC_INFO,"arc_ukf_filter_update : update filter using a measure \n");

    /* cholesky decomposition of the state covariance matrix */
    arc_ukf_cholesky_decomposition(filter->P,l,filter->sigma);
    
    arc_log(ARC_INFO,"arc_ukf_filter_update,cholesky decomposition P=\n ");
    arc_tracemat(ARC_MATPRINTF,filter->P,filter->state_dim,filter->state_dim,10,4);

    /* ================================= */
    /* compute sigma points */
    for (j=0;j<l;j++) {
        filter->sigma_point[j]=filter->x[j];
    }
    for (i=0;i<l;i++) { /* states numbers */
        for (j=0;j<i;j++) {
            filter->sigma_point[(i+1)*l+j]  =filter->x[j] ;
            filter->sigma_point[(i+1+l)*l+j]=filter->x[j] ;
        }
        filter->sigma_point[(i+1)*l+i]  =filter->x[i]+filter->gamma*filter->sigma[i];
        filter->sigma_point[(i+1+l)*l+i]=filter->x[i]-filter->gamma*filter->sigma[i];
        for (j=i+1;j<l;j++) { /* state numbers */
            filter->sigma_point[(i+1)*l+j]  =filter->x[j]+filter->gamma*filter->P[j*l+i];
            filter->sigma_point[(i+1+l)*l+j]=filter->x[j]-filter->gamma*filter->P[j*l+i];
        }
    }
    /* ================================= */
    /* propagate sigma points:y=f(xi) ,return khi[]*/
    for (i=0;i<2*l+1;i++) {  /* sigma point numbers */
        filter->ffun(filter->state_dim,&(filter->sigma_point[i*l]),&(filter->khi[i*l]));
    }
    /* check propagate sigma point whether is good */
    if (arc_norm(filter->khi,filter->state_dim)<=0.0) {
        arc_log(ARC_WARNING,"propagate sigma points failed \n");
        return 0;
    }
    arc_log(ARC_INFO,"sigma points:\n");
    arc_tracemat(ARC_MATPRINTF,filter->sigma_point,l,2*l+1,16,4);

    /* compute state prediction xm */
    for (i=0;i<l;i++) {  /* states numbers */
        filter->xm[i]=filter->wm[0]*filter->khi[i];
        for (j=1;j<2*l+1;j++) {  /* sigma point numbers */
            filter->xm[i]+=filter->wm[j]*filter->khi[j*l+i];
        }
    }
    arc_log(ARC_INFO,"arc_ukf_filter_update : propagate sigma points,"
            "its mean sigma point is : \n");
    arc_tracemat(ARC_MATPRINTF,filter->xm,l,1,16,4);
    
    /* ================================ */
    /* time update */
    /* start with state covariance matrix */

    /* noise model for states */
    for (i=0;i<l*l;i++) {
        filter->PM[i]=filter->Q[i];
    }
    /* accumulate covariances */
    for (i=0;i<2*l+1;i++) {
        for (j=0;j<l;j++) {
            filter->dx[j]=filter->khi[i*l+j]-filter->xm[j];
        }
        for (j=0;j<l;j++) {
            for (k=0;k<l;k++) {
                filter->PM[j*l+k]+=filter->wc[i]*filter->dx[j]*filter->dx[k];
            }
        }
    }
    /* save PM matrix */
    for(i=0;i<l*l;i++) {
        filter->PM_save[i]=filter->PM[i];
    }
    /* redraw sigma points */
    arc_ukf_cholesky_decomposition(filter->PM,l,filter->sigma);

    for (j=0;j<l;j++) {
        filter->sigma_point[j]=filter->xm[j];
    }
    for (i=0;i<l;i++) {
        for (j = 0;j<i;j++) {
            filter->sigma_point[(i+1)*l+j]  =filter->xm[j];
            filter->sigma_point[(i+1+l)*l+j]=filter->xm[j];
        }
        filter->sigma_point[(i+1)*l+i]  =filter->xm[i]+filter->gamma*filter->sigma[i];
        filter->sigma_point[(i+1+l)*l+i]=filter->xm[i]-filter->gamma*filter->sigma[i];
        for (j=i+1;j<l;j++) {
            filter->sigma_point[(i+1)*l+j]  =filter->xm[j]+filter->gamma*filter->PM[j*l+i];
            filter->sigma_point[(i+1+l)*l+j]=filter->xm[j]-filter->gamma*filter->PM[j*l+i];
        }
    }
    /* ================================= */
    /* propagate measurement */
    
    arc_log(ARC_INFO,"arc_ukf_filter_update : propagate measurement \n");
    
    for (i=0;i<2*l+1;i++) {
        filter->mfun(&(filter->sigma_point[i*l]),&(filter->khi_y[i*m]));
    }
    /* check propagate measurement */
    if (arc_norm(filter->khi_y,filter->measure_dim)<=0.0) {
        arc_log(ARC_WARNING,"ukf propagate measurement failed \n");
        return 0;
    }
    arc_log(ARC_INFO,"ukf propagate measurement = \n");
    arc_tracemat(ARC_MATPRINTF,filter->khi_y,m,2*l+1,16,4);

    /* measurement prediction */
    for (i=0;i<m;i++) {
        filter->ym[i]=filter->wm[0]*filter->khi_y[i];
        for(j=1;j<2*l+1;j++) {
            filter->ym[i]+=filter->wm[j]*filter->khi_y[j*m+i];
        }
    }
    arc_log(ARC_INFO,"measurement prediction meanings = \n");
    arc_tracemat(ARC_MATPRINTF,filter->ym,filter->measure_dim,1,13,4);

    /* measurement update
    /* Pyy matrix
    /* start with measure covariance matrix */
    for(i=0;i<m*m;i++) filter->Pyy[i]=filter->R[i];

    /* accumulate covariances */
    for(i=0;i<2*l+1;i++) {
        for(j=0;j<m;j++)
            filter->dy[j]=filter->khi_y[i*m+j]-filter->ym[j];
        for(j=0;j<m;j++)
            for(k=0;k<m;k++) {
                filter->Pyy[j*m+k]+=filter->wc[i]*filter->dy[j]*filter->dy[k];
            }
    }
    arc_log(ARC_INFO,"arc_ukf_filter_update : Pyy matrix \n");
    arc_tracemat(ARC_MATPRINTF,filter->Pyy,m,m,10,4);

    /* Pxy matrix */
    for(i=0;i<m*l;i++) filter->Pxy[i]=0.0;

    /* accumulate covariances */
    for(i=0;i<2*l+1;i++) {
        for(j=0;j<m;j++) {
            filter->dy[j]=filter->khi_y[i*m+j]-filter->ym[j];
        }
        for(j=0;j<l;j++) {
            filter->dx[j]=filter->sigma_point[i*l+j]-filter->xm[j];
        }
        for (j=0;j<l;j++) {
            for (k=0;k<m;k++) {
                filter->Pxy[j*m+k]+=filter->wc[i]*filter->dx[j]*filter->dy[k];
            }
        }
    }
    arc_log(ARC_INFO,"arc_ukf_filter_update : Pxy matrix \n");
    arc_tracemat(ARC_MATPRINTF,filter->Pxy,m,l,10,4);

    /* gain de kalman */
    arc_ukf_cholesky_decomposition(filter->Pyy,m,filter->sigma_y);
    arc_ukf_cholesky_solve(filter->Pyy,m,filter->sigma_y,filter->Pxy,l,filter->gain);

    /* restore PM matrix */
    for(i=0;i<l*l;i++) {
        filter->P[i]=filter->PM_save[i];
    }
    arc_log(ARC_INFO,"arc_ukf_filter_update :ukf measurements =\n");
    arc_tracemat(ARC_MATPRINTF,y,filter->measure_dim,1,10,4);
    
    /* update state */
    for (j=0;j<m;j++) {
        filter->dy[j]=y[j]-filter->ym[j];
    }
    arc_log(ARC_INFO,"arc_ukf_filter_update : dy \n");
    arc_tracemat(ARC_MATPRINTF,filter->dy,1,m,10,4);

    for (i=0;i<l;i++) {
        filter->x[i]=filter->xm[i];
        t=0.0;
        for (j=0;j<m;j++) {
            t+=filter->gain[i*m+j]*filter->dy[j];
        }
        filter->x[i]+=t;
    }
    arc_log(ARC_INFO,"arc_ukf_filter_update : x \n");
    arc_tracemat(ARC_MATPRINTF,filter->x,filter->state_dim,1,10,4);

    for(i=0;i<l;i++) {
        for(j=0;j<m;j++) {
            t=0.0;
            for(k=j;k<m;k++) {
                t+=filter->gain[i*m+k]*filter->Pyy[k*m+j];
            }
            filter->KL[i*m+j]=t;
        }
    }
    for (i=0;i<l;i++) {
        for(j=0;j<l;j++) {
            t=0.0;
            for (k=0;k<m;k++) {
                t+=filter->KL[i*m+k]*filter->KL[j*m+k];
            }
            filter->P[i*l+j]-=t;
        }
    }
    arc_log(ARC_INFO,"arc_ukf_filter_update : P \n");
    arc_tracemat(ARC_MATPRINTF,filter->P,filter->state_dim,filter->state_dim,10,4);
    return 1;
    /* finished with kalman iteration ! */
}
/* ukf_filter_delete -----------------------------------------------*/
extern void arc_ukf_filter_delete(ukf_t *filter)
{
    arc_log(ARC_INFO,"arc_ukf_filter_delete : \n");

    /* free ukf */
    filter->measure_dim=filter->state_dim=0;  /* todo:this function have some unknown bugs,must to fix */
    filter->gamma=0.0;

    if (filter->x)           free(filter->x);
    if (filter->y)           free(filter->y);
    if (filter->P)           free(filter->P);
    if (filter->Q)           free(filter->Q);
    if (filter->R)           free(filter->R);
    if (filter->wm)          free(filter->wm);
    if (filter->wc)          free(filter->wc);
    if (filter->sigma_point) free(filter->sigma_point);
    if (filter->sigma)       free(filter->sigma);
    if (filter->sigma_y)     free(filter->sigma_y);
    if (filter->PM)          free(filter->PM);
    if (filter->PM_save)     free(filter->PM_save);
    if (filter->xm)          free(filter->xm);
    if (filter->khi)         free(filter->khi);  /* todo:memory leak or memory overflow,need to fix */
    if (filter->ym)          free(filter->ym);
    if (filter->khi_y)       free(filter->khi_y);
    if (filter->Pyy)         free(filter->Pyy);
    if (filter->Pxy)         free(filter->Pxy);
    if (filter->dx)          free(filter->dx);
    if (filter->dy)          free(filter->dy);
    if (filter->gain)        free(filter->gain);
    if (filter->KL)          free(filter->KL);
    if (filter)              free(filter);  /* todo:here may have more better to complete */
    filter=NULL;
}

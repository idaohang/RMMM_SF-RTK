
#include "arc.h"

/* constants/macros ----------------------------------------------------------*/
#define LOOPMAX     10000           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<=(y)?(x):(y))

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
static int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A= arc_mat(n, n);

    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
        if ((D[i]=A[i+i*n])<=0.0) {info=-1; break;}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;

    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;

    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;

    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------*/
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=arc_zeros(n,n),*dist=arc_mat(n,1),
            *zb=arc_mat(n,1),*z=arc_mat(n,1),*step=arc_mat(n,1);

    k=n-1; dist[k]=0.0;
    zb[k]=zs[k];
    z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
    for (c=0;c<LOOPMAX;c++) {
        newdist=dist[k]+y*y/D[k];
        if (newdist<maxdist) {
            if (k!=0) {
                dist[--k]=newdist;
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); y=zb[k]-z[k]; step[k]=SGN(y);
            }
            else {
                if (nn<m) {
                    if (nn==0||newdist>s[imax]) imax=nn;
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];
                    s[nn++]=newdist;
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];
                }
                z[0]+=step[0]; y=zb[0]-z[0]; step[0]=-step[0]-SGN(step[0]);
            }
        }
        else {
            if (k==n-1) break;
            else {
                k++;
                z[k]+=step[k]; y=zb[k]-z[k]; step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {
            if (s[i]<s[j]) continue;
            SWAP(s[i],s[j]);
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);

    if (c>=LOOPMAX) {
        fprintf(stderr,"%s : search loop count overflow\n",__FILE__);
        return -1;
    }
    return 0;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int arc_lambda(int n, int m, const double *a, const double *Q, double *F,
                      double *s,double *Do,double *Lo)
{
    int info;
    double *L,*D,*Z,*z,*E;

    arc_log(ARC_INFO,"arc_lambda :\n");

    if (n<=0||m<=0) return -1;
    L=arc_zeros(n,n); D=arc_mat(n,1); Z=arc_eye(n); z=arc_mat(n,1); E=arc_mat(n,m);

    /* LD factorization */
    if (!(info=LD(n,Q,L,D))) {

        /* lambda reduction */
        reduction(n,L,D,Z);
        arc_matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */

        /* mlambda search */
        if (!(info=search(n,m,L,D,z,E,s))) {

            info=arc_solve("T",Z,E,n,m,F); /* F=Z'\E */
        }
    }
    if (Do) arc_matcpy(Do,D,n,1);
    if (Lo) arc_matcpy(Lo,L,n,n);

    arc_log(ARC_INFO,"L=\n");
    arc_tracemat(ARC_MATPRINTF,Lo,n,n,10,4);

    arc_log(ARC_INFO,"D=\n");
    arc_tracemat(ARC_MATPRINTF,Do,1,n,10,4);

    free(L); free(D); free(Z); free(z); free(E);
    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_reduction(int n, const double *Q, double *Z)
{
    double *L,*D;
    int i,j,info;

    if (n<=0) return -1;

    L=arc_zeros(n,n); D=arc_mat(n,1);

    for (i=0;i<n;i++) for (j=0;j<n;j++) {
            Z[i+j*n]=i==j?1.0:0.0;
        }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);

    free(L); free(D);
    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s)
{
    double *L,*D;
    int info;

    if (n<=0||m<=0) return -1;

    L=arc_zeros(n,n); D=arc_mat(n,1);

    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* mlambda search */
    info=search(n,m,L,D,a,F,s);

    free(L); free(D);
    return info;
}
static int arc_exctract_L(const double *L,int k,double *LL,int m,int n)
{
    int i,j;

    for (i=0;i<n;i++) LL[i]=L[m-1+k*i]; return n;
}
extern int arc_bootstrap(int n,const double *a, const double *Q, double *F,double *Ps)
{
    double *L,*D,*Z,*an,*ZT;
    double *af,*afc,*S,*zhat,*LL;
    int i,j,k,info;

    if (n<=0) return -1;

    L=arc_zeros(n,n); D=arc_mat(n,1); Z=arc_eye(n);

    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D); free(Z);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);

    if (Ps) *Ps=arc_amb_bs_success(D,n);

    zhat=arc_mat(n,1);
    an=arc_mat(n,1);
    ZT=arc_mat(n,n);

    for (i=0;i<n;i++) an[i]=a[i]-int(a[i]);

    arc_matmul("TN",n,1,n,1.0,Z,an,0.0,zhat); /* z=Z'*a */

    af=arc_zeros(n,1); afc=arc_zeros(n,1);
    S=arc_zeros(n,1); LL=arc_mat(n,1);

    afc[n-1]=zhat[n-1]; af[n-1]=ROUND(afc[n-1]);

    for (i=n-1;i>0;i--) {
        j=arc_exctract_L(L,n,LL,i+1,i);

        for (k=0;k<j;k++) S[k]=S[k]+(af[j]-afc[j])*LL[k];
        afc[j-1]=zhat[j-1]+S[j-1];
        af[j-1]=ROUND(afc[j-1]);
    }

    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) ZT[i+j*n]=Z[j+i*n];
    }
    if (!(info=arc_matinv(ZT,n))) {
        arc_matmul("NN",n,1,n,1.0,ZT,af,0.0,F);
        for (i=0;i<n;i++) F[i]+=int(a[i]);
    }

    free(L); free(D); free(Z);
    free(af); free(afc);
    free(S); free(LL);

    return info;
}
static int arc_exctract_vec(const double *D,int n,int k,double *Do)
{
    int i,j; for (i=k-1,j=0;i<n;i++) Do[j++]=D[i]; return j;
}
static void arc_exctract_mat(const double *A,int r,int c,int a1,int b1,
                             int a2,int b2,double *Ao)
{
    int i,j,rs,re,cs,ce,nr,nc;

    rs=MIN(a1,a2); re=MAX(a1,a2);
    cs=MIN(b1,b2); ce=MAX(b1,b2);
    nr=re-rs+1; nc=ce-cs+1;

    for (i=0;i<nc;i++) {
        for (j=0;j<nr;j++) Ao[i*nr+j]=A[(cs-1+i)*r+rs-1+j];
    }
}
static int arc_parsearch(const double *zhat,const double *Qzhat,const double *Z,
                         const double *L,const double *D,double P0,int n,int m,
                         double *F,double*s)
{
    int i,j,k=1,info,l;
    double ps,*DD,*zz,*LL,*Qzpar,*QP,*zpar;

    ps=arc_amb_bs_success(D,n);

    DD=arc_mat(1,n);
    zz=arc_mat(1,n); LL=arc_mat(n,n);

    while(ps<P0&&k<n) {
        k=k+1;
        j=arc_exctract_vec(D,n,k,DD);
        ps=arc_amb_bs_success(DD,j);
    }
    Qzpar=arc_mat(n,n);
    QP=arc_mat(n,n); zpar=arc_mat(m,n);

    if (ps>P0) {

        if (k==1) {
            search(n,m,L,D,zhat,F,s);
        }
        else {
            arc_exctract_vec(zhat,n,k,zz);
            arc_exctract_vec(D,n,k,DD);
            arc_exctract_mat(L,n,n,k,k,n,n,LL);

            search(n-k+1,m,LL,DD,zz,zpar,s);

            arc_exctract_mat(Qzhat,n,n,k,k,n,n,Qzpar);
            arc_exctract_mat(Qzhat,n,n,k-1,k,1,n,LL);

            if (!(info=arc_matinv(Qzpar,n-k+1))) {

                arc_matmul("NN",k-1,n-k+1,n-k+1,1.0,LL,Qzpar,0.0,QP);

                for (i=0;i<m;i++) {
                    for (j=0;j<n-k+1;j++) {
                        DD[j]=zz[j]-zpar[i*(n-k+1)+j];
                    }
                    arc_matmul("NN",k-1,1,n-k+1,1.0,QP,DD,0.0,LL);
                    for (j=0;j<k-1;j++) {
                        F[i*n+j]=ROUND(zhat[j]-LL[j]);
                    }
                    for (j=k-1,l=0;j<n;j++) {
                        F[i*n+j]=zpar[i*(n-k+1)+l++];
                    }
                }
            }
        }
    }

    free(DD); free(LL); free(zz);
    free(Qzpar); free(QP);
    free(zpar);

    return info;
}
extern int arc_par_lambda(const double *a,const double *Qa,int n,int m,double *F,
                          double *s,double p0)
{
    double *L,*D,*Z,*zhat,*Qz,*Q,*E;
    int info;

    if (n<=0) return -1;

    L=arc_zeros(n,n); D=arc_mat(n,1); Z=arc_eye(n);
    E=arc_mat(m,n);

    /* LD factorization */
    if ((info=LD(n,Qa,L,D))) {
        free(L); free(D); free(Z);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);

    zhat=arc_mat(n,1);
    Qz=arc_mat(n,n); Q=arc_mat(n,n);

    arc_matmul("TN",n,1,n,1.0,Z,a,0.0,zhat); /* z=Z'*a */
    arc_matmul("TN",n,n,n,1.0,Z,Qa,0.0,Q);
    arc_matmul("NN",n,n,n,1.0,Q,Z,0.0,Qz); /* Qz=Z'*Qa*Z */

    arc_parsearch(zhat,Qz,Z,L,D,p0,n,m,E,s);

    info=arc_solve("T",Z,E,n,m,F); /* F=Z'\E */

    free(L); free(D); free(Z);
    free(zhat); free(Qz); free(Q); free(E);

    return info;
}
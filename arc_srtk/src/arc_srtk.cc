
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
 *  Created on: July 07, 2017
 *********************************************************************************/

#include "arc.h"
#include <iomanip>
#include <fstream>
#include <arc.h>

using namespace std;

/* constants/macros ----------------------------------------------------------*/
#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define VAR_POS     SQR(10.0)       /* initial variance of receiver pos (m^2) */
#define VAR_GRA     SQR(0.001)      /* initial variance of gradient (m^2) */
#define VAR_AMB     SQR(10.0)       /* initial variance of ambguity (cycle^2) */
#define INIT_ZWD    0.15            /* initial zwd (m) */
#define VAR_POS_AMB SQR(2.0)        /* after ambiguity fix,the variance of position */

#define GAP_RESION  120             /* gap to reset ionosphere parameters (epochs) */
#define VAR_HOLDAMB 0.001           /* constraint to hold ambiguity (cycle^2) */
#define VAR_AMB_GROUP_EL SQR(0.1)   /* constraint to high elevation double-difference ambiguity */

#define TTOL_MOVEB  (1.0+2*DTTOL)
                                    /* time sync tolerance for moving-baseline (s) */
#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define AMB_THRES   0.5             /* cycles */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:6)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)      /* phase bias (s:satno,f:freq) */

/* constants/global variables for ukf problem--------------------------------------------*/
static int _UKF_NX_      =0;                        /* ukf problem statec numbers */
static int _UKF_NV_      =0;                        /* ukf problem measurements numbers */
static rtk_t *_UKF_RTK_  =NULL;                     /* ukf problem rtk data struct pointor */
static const nav_t *_UKF_NAV_=NULL;                 /* ukf problem navigation data struct pointor */
static const obsd_t*_UKF_OBS_=NULL;                 /* ukf problem observation data pointor */
static int _UKF_NU_      =0;                        /* ukf problem rover station observation data numbers */
static int _UKF_NR_      =0;                        /* ukf problem base station observation data numbers */
static double *_UKF_RS_  =NULL;                     /* ukf problem satellite position/velecity pointor */
static double *_UKF_DTS_ =NULL;                     /* ukf problem satellite clock/clock drift pointor */
static double *_UKF_AZEL_=NULL;                     /* ukf problem base and rover station - satellite azel */
static double *_UKF_E_   =NULL;                     /* ukf problem base and rover station - satellite's line of light */
static int *_UKF_SVH_    =NULL;                     /* ukf problem satellite health flag */
static double *_UKF_R_   =NULL;                     /* ukf problem measurements variance */
static double *_UKF_Y_   =NULL;                     /* ukf problem sigma point's single-difference measurement vector */
static double *_UKF_DY_  =NULL;                     /* ukf problem sigma point's double-difference measurement vector */
static int    *_UKF_IX_  =NULL;                     /* ukf problem record the active states index */
static int     _UKF_ANX_ =0;                        /* ukf problem active state numbers */
static double *_UKF_Q_   =NULL;                     /* ukf problem system noise matrix */
static double *_UKF_XP_  =NULL;                     /* ukf problem prior states value */
static double *_UKF_PP_  =NULL;                     /* ukf problem prior states covariance matrix */
static double *_UKF_MEAS_=NULL;                     /* ukf problem measurements vector */
static int *_UKF_SAT_    =NULL;                     /* ukf problem rover station and base station common satellite list */
static int *_UKF_IU_     =NULL;
static int *_UKF_IR_     =NULL;
static int  _UKF_NS_     =0;                        /* ukf problem rover station and base station common satellite numbers */

/* Macro for defining an exception------------------------------------------------------*/
ARC_DEFINE_EXCEPTION(Exception,std::runtime_error);

/* single-differenced observable ---------------------------------------------*/
static double arc_sdobs(const obsd_t *obs, int i, int j, int f)
{
    double pi=f<NFREQ?obs[i].L[f]:obs[i].P[f-NFREQ];
    double pj=f<NFREQ?obs[j].L[f]:obs[j].P[f-NFREQ];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* single-differenced measurement error variance -----------------------------*/
static double arc_varerr(int sat, int sys, double el, double bl, double dt, int f,
                         const prcopt_t *opt)
{
    double a,b,c=opt->err[3]*bl/1E4,d=CLIGHT*opt->sclkstab*dt,fact=1.0;
    double sinel=sin(el),geo_factor=1.0;
    int i=sys==SYS_GLO?1:(sys==SYS_GAL?2:0),nf=NF(opt);

    /* bds geo */
    if (sys==SYS_CMP&&arc_is_bds_geo(sat)) geo_factor=EFACT_GEO; /* todo:this may be no effect,now just for test */

    /* extended error model */
    if (f>=nf&&opt->exterr.ena[0]) { /* code */
        a=opt->exterr.cerr[i][  (f-nf)*2];
        b=opt->exterr.cerr[i][1+(f-nf)*2];
        if (sys==SYS_SBS) {a*=EFACT_SBS; b*=EFACT_SBS;}
    }
    else if (f<nf&&opt->exterr.ena[1]) { /* phase */
        a=opt->exterr.perr[i][  f*2];
        b=opt->exterr.perr[i][1+f*2];
        if (sys==SYS_SBS) {a*=EFACT_SBS; b*=EFACT_SBS;}
    }
    else { /* normal error model */
        if (f>=nf) fact=opt->eratio[f-nf];
        if (fact<=0.0)  fact=opt->eratio[0];
        fact*=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
        a=fact*opt->err[1];
        b=fact*opt->err[2];
    }
    return 2.0*geo_factor*(opt->ionoopt==IONOOPT_IFLC?3.0:1.0)*(a*a+b*b/sinel/sinel+c*c)+d*d;
}
/* baseline length -----------------------------------------------------------*/
static double arc_baseline(const double *ru,const double *rb,double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return arc_norm(dr, 3);
}
/* initialize state and covariance -------------------------------------------*/
static void arc_initx(rtk_t *rtk,double xi,double var,int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}
static void arc_diff_pr_initx(double *x,double *P,double xi,double var,int i,
                              int nx)
{
    int j; x[i]=xi; for (j=0;j<nx;j++) P[i+j*nx]=P[j+i*nx]=i==j?var:0.0;
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void arc_testeclipse(const obsd_t *obs,int n,const nav_t *nav,double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;

    arc_log(3,"testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    arc_sunmoonpos(gpst2utc(obs[0].time), erpv, rsun, NULL, NULL);
    arc_normv3(rsun,esun);

    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;

        if ((r=arc_norm(rs+i*6,3))<=0.0) continue;

        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;

        /* sun-earth-satellite angle */
        cosa=arc_dot(rs+i*6,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);

        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;

        arc_log(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
                obs[i].sat);
        for (j=0;j<3;j++) rs[j+i*6]=0.0;
    }
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double arc_yaw_nominal(double beta,double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* shadow-crossing GPS IIA ---------------------------------------------------*/
static int arc_yaw_shadow_IIA(double beta,double mu,double eps0,double R,
                              double mudot,double *yaw)
{
    double mu_s,mu_e;

    mu_s=-sqrt(SQR(eps0)-SQR(beta));
    mu_e=-mu_s;

    if (mu_s<=mu&&mu<mu_e) {
        *yaw=atan2(-tan(beta),sin(mu_s))+R*(mu-mu_s)/mudot;
    }
    else if (mu_e<=mu&&mu<mu_e+T_POSTSHADOW*mudot) {
        return 0;
    }
    return 1;
}
/* shadow-crossing GLONASS-M -------------------------------------------------*/
static int arc_yaw_shadow_GLO(double beta,double mu,double eps0,double R,
                              double mudot,double *yaw)
{
    double mu_s,mu_e,mu_f,tan_beta,sin_mu_s;

    if (beta<0) R=-R;
    tan_beta=tan(beta);

    mu_s=-acos(cos(eps0)/cos(beta));
    mu_e=-mu_s;
    sin_mu_s=sin(mu_s);
    mu_f=mudot*(atan2(-tan_beta,-sin_mu_s)-atan2(-tan_beta,sin_mu_s))/R+mu_s;

    if (mu_s<=mu&&mu<mu_f) {
        *yaw=atan2(-tan_beta,sin_mu_s)+R*(mu-mu_s)/mudot;
    }
    else if (mu_f<=mu&&mu<mu_e) {
        *yaw=atan2(-tan_beta,-sin_mu_s);
    }
    return 1;
}
/* noon-turn maneuver --------------------------------------------------------*/
static int arc_yaw_noon(double beta,double mu,double beta0,double R,
                        double mudot,double *yaw)
{
    double mu_s,y;

    if (beta>=0) R=-R;
    mu_s=PI-sqrt(beta0*fabs(beta)-SQR(beta));

    if (mu_s<=mu) {
        y=atan2(-tan(beta),sin(mu_s))+R*(mu-mu_s)/mudot;
        if ((beta>=0&&y>*yaw)||(beta<0&&y<*yaw)) *yaw=y;
    }
    return 1;
}
/* midnight-turn maneuver ----------------------------------------------------*/
static int arc_yaw_midnight(double beta,double mu,double beta0,double R,
                            double mudot,double *yaw)
{
    double mu_s,y;

    if (beta<0) R=-R;

    mu_s=-sqrt(beta0*fabs(beta)-SQR(beta));
    if (mu_s<=mu) {
        y=atan2(-tan(beta),sin(mu_s))+R*(mu-mu_s)/mudot;
        if ((beta>=0&&y<*yaw)||(beta<0&&y>*yaw)) *yaw=y;
    }
    return 1;
}
/* yaw-angle of GPS IIA (ref [8]) --------------------------------------------*/
static int arc_yaw_IIA(int sat,int opt,double beta,double mu,double *yaw)
{
    const double R_GPSIIA[]={
            0.1046,0.1230,0.1255,0.1249,0.1003,0.1230,0.1136,0.1169,0.1253,0.0999,
            0.1230,0.1230,0.1230,0.1230,0.1092,0.1230,0.1230,0.1230,0.1230,0.1230,
            0.1230,0.1230,0.1230,0.0960,0.0838,0.1284,0.1183,0.1230,0.1024,0.1042,
            0.1230,0.1100,0.1230
    };
    double R=R_GPSIIA[sat-1]*D2R,beta0=atan(MUDOT_GPS/R);

    *yaw=atan2(-tan(beta),sin(mu));

    if (opt==2) { /* precise yaw */
        if (mu<PI/2.0&&fabs(beta)<EPS0_GPS) {
            if (!arc_yaw_shadow_IIA(beta,mu,EPS0_GPS,R,MUDOT_GPS,yaw)) return 0;
        }
        else if (mu>PI/2.0&&fabs(beta)<beta0) {
            if (!arc_yaw_noon(beta,mu,beta0,R,MUDOT_GPS,yaw)) return 0;
        }
    }
    return 1;
}
/* yaw-angle of GPS IIR (ref [8]) --------------------------------------------*/
static int arc_yaw_IIR(int sat,int opt,double beta,double mu,double *yaw)
{
    const double R=0.2*D2R;
    double beta0=atan(MUDOT_GPS/R);

    *yaw=atan2(-tan(beta),sin(mu));

    if (opt==2) { /* precise yaw */
        if (mu<PI/2.0&&fabs(beta)<beta0) {
            if (!arc_yaw_midnight(beta,mu,beta0,R,MUDOT_GPS,yaw)) return 0;
        }
        else if (mu>PI/2.0&&fabs(beta)<beta0) {
            if (!arc_yaw_noon(beta,mu,beta0,R,MUDOT_GPS,yaw)) return 0;
        }
    }
    *yaw+=PI;
    return 1;
}
/* yaw-angle of GPS IIF (ref [9]) --------------------------------------------*/
static int arc_yaw_IIF(int sat,int opt,double beta,double mu,double *yaw)
{
    const double R0=0.06*D2R,R1=0.11*D2R;
    double beta0=atan(MUDOT_GPS/R1);

    *yaw=atan2(-tan(beta),sin(mu));

    if (opt==2) { /* precise yaw */
        if (fabs(mu)<EPS0_GPS&&fabs(beta)<EPS0_GPS) {
            if (!arc_yaw_shadow_GLO(beta,mu,EPS0_GPS,R0,MUDOT_GPS,yaw)) return 0;
        }
        else if (mu>PI/2.0&&fabs(beta)<beta0) {
            if (!arc_yaw_noon(beta,mu,beta0,R1,MUDOT_GPS,yaw)) return 0;
        }
    }
    return 1;
}
/* yaw-angle of Galileo (ref [11]) -------------------------------------------*/
static int arc_yaw_GAL(int sat,int opt,double beta,double mu,double *yaw)
{
    *yaw=arc_yaw_nominal(beta,mu);
    return 1;
}
static int arc_yaw_CMP(int sat,int opt,double beta,double mu,double *yaw)
{
    *yaw=0.0;
    return 1;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int arc_yaw_angle(int sat,const char *type,int opt,double beta,double mu,
                         double *yaw)
{
    if      (strstr(type,"BLOCK IIA")) return arc_yaw_IIA(sat,opt,beta,mu,yaw);
    else if (strstr(type,"BLOCK IIR")) return arc_yaw_IIR(sat,opt,beta,mu,yaw);
    else if (strstr(type,"BLOCK IIF")) return arc_yaw_IIF(sat,opt,beta,mu,yaw);
    else if (strstr(type,"Galileo"  )) return arc_yaw_GAL(sat,opt,beta,mu,yaw);
    else if (strstr(type,"BEIDOU"   )) return arc_yaw_CMP(sat,opt,beta,mu,yaw);
    return 0;
}
/* satellite attitude model --------------------------------------------------*/
static int arc_sat_yaw(gtime_t time,int sat,const char *type,int opt,
                       const double *rs,double *exs,double *eys)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0};
    int i;

    arc_sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

    /* beta and orbit angle */
    arc_matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    arc_cross3(ri,ri+3,n);
    arc_cross3(rsun,n,p);
    if (!arc_normv3(rs,es)||!arc_normv3(rsun,esun)||!arc_normv3(n,en)||
        !arc_normv3(p,ep)) return 0;
    beta=PI/2.0-acos(arc_dot(esun,en,3));
    E=acos(arc_dot(es,ep,3));
    mu=PI/2.0+(arc_dot(es,esun,3)<=0?-E:E);
    if      (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;

    /* yaw-angle of satellite */
    if (!arc_yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;

    /* satellite fixed x,y-vector */
    arc_cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/
static int arc_model_phw(gtime_t time,int sat,const char *type,int opt,
                         const double *rs,const double *rr,double *phw)
{
    double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
    double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
    int i;

    if (opt<=0) return 1; /* no phase windup */

    /* satellite yaw attitude model */
    if (!arc_sat_yaw(time,sat,type,opt,rs,exs,eys)) return 0;

    /* unit vector satellite to receiver */
    for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
    if (!arc_normv3(r,ek)) return 0;

    /* unit vectors of receiver antenna */
    ecef2pos(rr,pos);
    xyz2enu(pos,E);
    exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
    eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */

    /* phase windup effect */
    arc_cross3(ek,eys,eks);
    arc_cross3(ek,eyr,ekr);
    for (i=0;i<3;i++) {
        ds[i]=exs[i]-ek[i]*arc_dot(ek,exs,3)-eks[i];
        dr[i]=exr[i]-ek[i]*arc_dot(ek,exr,3)+ekr[i];
    }
    cosp=arc_dot(ds,dr,3)/arc_norm(ds,3)/arc_norm(dr,3);
    if      (cosp<-1.0) cosp=-1.0;
    else if (cosp> 1.0) cosp= 1.0;
    ph=acos(cosp)/2.0/PI;
    arc_cross3(ds,dr,drs);
    if (arc_dot(ek,drs,3)<0.0) ph=-ph;

    *phw=ph+floor(*phw-ph+0.5); /* in cycle */
    return 1;
}
/* select common satellites between rover and reference station --------------*/
static int arc_selsat(const obsd_t *obs,double *azel,int nu,int nr,
                      const prcopt_t *opt,int *sat,int *iu,int *ir)
{
    int i,j,k=0;

    arc_log(ARC_INFO, "nu=%d nr=%d\n", nu, nr);

    for (i=0,j=nu;i<nu&&j<nu+nr;i++,j++) {
        if      (obs[i].sat<obs[j].sat) j--;
        else if (obs[i].sat>obs[j].sat) i--;
        else if (azel[1+j*2]>=opt->elmin) { /* elevation at base station */
            sat[k]=obs[i].sat; iu[k]=i; ir[k++]=j;
            arc_log(4, "(%2d) sat=%3d iu=%2d ir=%2d\n", k-1,obs[i].sat,i,j);
        }
    }
    return k;
}
/* temporal update of position/velocity/acceleration -------------------------*/
static void arc_udpos(rtk_t *rtk,double tt)
{
    /* todo:dynamic mode is not good,cause of this is not clearly */
    int i,j;
    double var;

    arc_log(ARC_INFO,"arc_udpos   : tt=%.3f\n", tt);

    /* fixed mode */
    /* todo:this mode may useless */
    if (rtk->opt.mode==PMODE_FIXED) {
        for (i=0;i<3;i++) arc_initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (arc_norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        /* initial rover station velecity */
        if (rtk->opt.dynamics) for (i=3;i<6;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
    }

    /* update ceres solver problem active states index list,included velecity */
    if (rtk->opt.dynamics) for (i=0;i<6;i++) rtk->ceres_active_x[i]=1;  /* todo:set active states index in here may be not best */
    else for (i=0;i<3;i++) rtk->ceres_active_x[i]=1;

    /* static mode */
    if (rtk->opt.mode==PMODE_STATIC) return;

    /* reset rover station pistion and its variance */
    /* todo:using standard positioning to initial ukf prior states and its covariacne matrix,
     * todo:and may have more better methids to do this */

    /* kinamic without dynamic mode */
    if (!rtk->opt.dynamics) {
        if (rtk->opt.init_pnt) { /* use pnt to initial rtk */
            for (i=0;i<3;i++) {
                arc_initx(rtk,rtk->sol.rr[i],rtk->sol.qr[i],i); return;
            }
        }
        else if (rtk->opt.init_dc) { /* use dd-pseudorange to initial rtk */
            for (i=0;i<3;i++) rtk->P[i+i*rtk->nx]=rtk->sol.prsol.qr[i];
            rtk->P[1+rtk->nx*0]=rtk->P[0+rtk->nx*1]=rtk->sol.prsol.qr[3];
            rtk->P[2+rtk->nx*0]=rtk->P[0+rtk->nx*2]=rtk->sol.prsol.qr[5];
            rtk->P[1+rtk->nx*2]=rtk->P[2+rtk->nx*1]=rtk->sol.prsol.qr[4]; /* todo:may be have more better way to do this */
            return;
        }
        else { /* use prior information to initial rtk */
            if (rtk->sol.ratio>=rtk->opt.thresar[0]) {
                for (i=0;i<3;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_POS_AMB,i);
                return; /* here may more better than just use pnt to initial */
            }
            /* pnt to initial */
            for (i=0;i<3;i++) rtk->P[i+i*rtk->nx]=rtk->sol.qr[i];
            rtk->P[1+rtk->nx*0]=rtk->P[0+rtk->nx*1]=rtk->sol.qr[3];
            rtk->P[2+rtk->nx*0]=rtk->P[0+rtk->nx*2]=rtk->sol.qr[5];
            rtk->P[1+rtk->nx*2]=rtk->P[2+rtk->nx*1]=rtk->sol.qr[4]; /* todo:may be have more better way to do this */
            return;
        }
    }
    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=rtk->P[i+i*rtk->nx]; var/=3.0;

    if (var>VAR_POS||var<0.0) {
        /* reset position with large variance */
        for (i=0;i<3;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        arc_log(ARC_INFO,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
    /* state transition of position/velocity/acceleration for dynamic mode */
    double *F=arc_eye(rtk->nx),*FP=arc_mat(rtk->nx,rtk->nx),
           *xp=arc_mat(rtk->nx,1),pos[3]={0},Q[9]={0},Qv[9]={0};

    /* compute the state transition matrix */
    for (i=0;i<3;i++) F[i+(i+3)*rtk->nx]=tt;

    arc_log(ARC_INFO,"arc_udpos : state transition matrix =\n");
    arc_tracemat(ARC_MATPRINTF,F,rtk->nx,rtk->nx,10,4);

    arc_log(ARC_INFO,"arc_udpos : before state transition P =\n");
    arc_tracemat(ARC_MATPRINTF,rtk->P,rtk->nx,rtk->nx,10,4);

    /* x=F*x, P=F*P*F+Q */
    arc_matmul("NN",rtk->nx,1,rtk->nx,1.0,F,rtk->x,0.0,xp);
    arc_matcpy(rtk->x,xp,rtk->nx,1);
    arc_matmul("NN",rtk->nx,rtk->nx,rtk->nx,1.0,F,rtk->P,0.0,FP);
    arc_matmul("NT",rtk->nx,rtk->nx,rtk->nx,1.0,FP,F,0.0,rtk->P);

    arc_log(ARC_INFO,"arc_udpos : after state transition P = \n");
    arc_tracemat(ARC_MATPRINTF,rtk->P,rtk->nx,rtk->nx,10,4);

    /* process noise added to only velecity */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3]); Q[8]=SQR(rtk->opt.prn[4]);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        rtk->P[i+3+(j+3)*rtk->nx]+=Qv[i+j*3];  /* add system process noise */
    }
    free(F); free(FP); free(xp);
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void arc_udion(rtk_t *rtk,double tt,double bl,const int *sat,int ns)
{
    double el,fact;
    int i,j;

    arc_log(ARC_INFO, "arc_udion   : tt=%.1f bl=%.0f ns=%d\n",tt,bl,ns);

    for (i=1;i<=MAXSAT;i++) {
        j=II(i,&rtk->opt);
        if (rtk->x[j]!=0.0&&
            rtk->ssat[i-1].outc[0]>GAP_RESION&&rtk->ssat[i-1].outc[1]>GAP_RESION)
            rtk->x[j]=0.0;
    }
    for (i=0;i<ns;i++) {
        j=II(sat[i],&rtk->opt);

        if (rtk->x[j]==0.0) {
            arc_initx(rtk,1E-6,SQR(rtk->opt.std[1]*bl/1E4),j);
        }
        else {
            /* elevation dependent factor of process noise */
            el=rtk->ssat[sat[i]-1].azel[1];
            fact=cos(el);
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[1]*bl/1E4*fact)*tt;
        }
    }
}
/* temporal update of tropospheric parameters --------------------------------*/
static void arc_udtrop(rtk_t *rtk,double tt,double bl)
{
    int i,j,k;

    arc_log(ARC_INFO, "arc_udtrop  : tt=%.1f\n", tt);

    for (i=0;i<2;i++) {
        j=IT(i,&rtk->opt);

        if (rtk->x[j]==0.0) {
            arc_initx(rtk,INIT_ZWD,SQR(rtk->opt.std[2]),j); /* initial zwd */

            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) arc_initx(rtk,1E-6,VAR_GRA,++j);
            }
        }
        else {
            rtk->P[j+j*rtk->nx]+=SQR(rtk->opt.prn[2])*tt;

            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) {
                    rtk->P[++j*(1+rtk->nx)]+=SQR(rtk->opt.prn[2]*0.3)*fabs(rtk->tt);
                }
            }
        }
    }
}
/* detect cycle slip by double-difference inter-stations and epoches ------------*/
static int arc_detslp_ddre(rtk_t *rtk, const obsd_t *obs, const int *iu, const int *ir,
                           int ns, const nav_t *nav, const double *rs)
{
    double *H,*v,*er,*eb;
    double dx[4],Qx[4*4];
    double dr,db;
    int i,j,k,nv,stat=0;
    int sat[MAXSAT];

    arc_log(ARC_INFO,"arc_detslp_ddre : \n");

    if (ns<5) {
        arc_log(ARC_WARNING,"arc_detslp_ddre: lack of satellate,ns=%d",ns);
        return 0;
    }
    H=arc_mat(ns,4); v=arc_mat(ns,1); er=arc_mat(ns,3); eb=arc_mat(ns,3);
    for (nv=i=0;i<ns;i++) {
        if (obs[iu[i]].sat!=obs[ir[i]].sat) continue;  /* rover and base station common satellite */
        sat[nv]=obs[iu[i]].sat;

        if (rtk->x[IB(sat[nv],0,&rtk->opt)]==0) continue;
        if (rtk->ssat[sat[nv]-1].slip[0]) continue;
        if (rtk->ssat[sat[nv]-1].ph[0][0]==0
            ||rtk->ssat[sat[nv]-1].ph[1][0]==0) continue;

        if (rtk->ssat[sat[nv]-1].r0[0]==0
            ||rtk->ssat[sat[nv]-1].r0[1]==0) continue;

        dr=arc_geodist(rs+iu[i]*6,rtk->sol.rr,er+i*3); /* rover station to satellite distance */
        db=arc_geodist(rs+ir[i]*6,rtk->rb,eb+i*3); /* base station to satellite distance */

        for (j=0;j<3;j++) {
            H[nv*4+j]=-er[i*3+j]; /* fill design matrix */
        }
        H[nv*4+3]=1;

        v[nv]=((obs[iu[i]].L[0]-obs[ir[i]].L[0])-(rtk->ssat[sat[nv]-1].ph[0][0]-rtk->ssat[sat[nv]-1].ph[1][0]))
              *nav->lam[sat[nv]-1][0]-((dr-db)-(rtk->ssat[sat[nv]-1].r0[0]-rtk->ssat[sat[nv]-1].r0[1]));
        nv++;
    }
    if (nv<5) {
        arc_log(ARC_WARNING,"arc_detslp_ddre: lack of satellate,nv=%d",nv);
        free(H); free(v); free(er); free(eb);
        return 0;
    }
    if (arc_lsq(H,v,4,nv,dx,Qx)>0) {
        arc_log(ARC_WARNING,"arc_detslp_ddre: lsq error");
        free(H); free(v); free(er); free(eb);
        return 0;
    }
    arc_tracemat(ARC_MATPRINTF,v,1,nv,10,4);
    arc_matmul("TN",nv,1,4,1.0,H,dx,-1,v); /* v=H'*dx-v */

    double *Q; Q=arc_mat(4,4);
    arc_matmul("NT",4,4,nv,1.0,H,H,0.0,Q); /* Q=H*H' */
    arc_matinv(Q,4);
    double *QH; QH=arc_mat(nv,4);
    arc_matmul("TN",nv,4,4,1.0,H,Q,0.0,QH); /* QH=H'*Q-1 */
    double *Qvv; Qvv=arc_eye(nv);
    arc_matmul("NN",nv,nv,4,-1.0,QH,H,1,Qvv); /* Qvv=Qvv^-1-H'*Q-1*H */

    double delta;
    delta=sqrt(arc_dot(v,v,nv)/(nv-4));
    double *vv;
    vv=arc_mat(nv,1);
    for (i=0;i<nv;i++) {
        vv[i]=v[i]/(sqrt(Qvv[i*nv+i])*delta);
    }

    double maxv=0;
    for (k =i=0;i<nv;i++) {
        if (fabs(vv[i])>fabs(maxv)) {
            maxv=vv[i]; k=i;
        }
    }
    /* todo:maybe here have bugs(for example,v=[-991.2134 -989.5139 -990.1663 -991.1516 -988.6747 -990.3013]) */
    if (fabs(maxv)>2&&fabs(v[k])>0.02) stat=sat[k];

    free(H); free(v); free(er); free(eb);
    free(Q); free(QH); free(Qvv); free(vv);

    return stat;
}
/* ------------------------------------------------------------------------- */
static void arc_detsl_new(rtk_t *rtk, const obsd_t *obs, const int *iu, const int *ir,
                          int ns, const nav_t *nav, const double *rs)
{
    int sat=0,i;

    arc_log(ARC_INFO,"arc_detsl_new : \n");

    for (i=0;i<ns;i++) {
        sat=arc_detslp_ddre(rtk,obs,iu,ir,ns,nav,rs);
        if (sat!=0) {
            arc_log(ARC_WARNING,"arc_detsl_new: slip detected,sat=%d,%s",
                    sat,time_str(obs[0].time,3));
            rtk->ssat[sat-1].slip[0]|=1;
        }
        if (sat==1) break;
    }  /* todo:this condition may have bugs */
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void arc_detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv)
{
    unsigned int slip,LLI;
    int f=0,sat=obs[i].sat;

    arc_log(ARC_INFO, "arc_detslp_ll: i=%d rcv=%d\n",i,rcv);

    if (obs[i].L[f]==0.0) return;

    /* restore previous LLI */
    if (rcv==1) LLI=getbitu(&rtk->ssat[sat-1].slip[f],0,2); /* rover */
    else        LLI=getbitu(&rtk->ssat[sat-1].slip[f],2,2); /* base  */

    /* detect slip by cycle slip flag in LLI */
    if (rtk->tt>=0.0) { /* forward */
        if (obs[i].LLI[f]&1) {
            arc_log(ARC_WARNING, "arc_detslp_ll : "
                            "slip detected forward (sat=%2d rcv=%d F=%d LLI=%x)\n",
                    sat,rcv,f+1,obs[i].LLI[f]);
        }
        slip=obs[i].LLI[f];
    }
    else { /* backward */
        if (LLI&1) {
            arc_log(ARC_WARNING, "arc_detslp_ll : "
                            "slip detected backward (sat=%2d rcv=%d F=%d LLI=%x)\n",
                    sat,rcv,f+1,LLI);
        }
        slip=LLI;
    }
    /* detect slip by parity unknown flag transition in LLI */
    if (((LLI&2)&&!(obs[i].LLI[f]&2))||(!(LLI&2)&&(obs[i].LLI[f]&2))) {
        arc_log(ARC_WARNING, "arc_detslp_ll : "
                        "slip detected half-cyc (sat=%2d rcv=%d F=%d LLI=%x->%x)\n",
                sat,rcv,f+1,LLI,obs[i].LLI[f]);
        slip|=1;
    }
    /* save current LLI */
    if (rcv==1) setbitu(&rtk->ssat[sat-1].slip[f],0,2,obs[i].LLI[f]);
    else        setbitu(&rtk->ssat[sat-1].slip[f],2,2,obs[i].LLI[f]);

    /* save slip and half-cycle valid flag */
    rtk->ssat[sat-1].slip[f]|=(unsigned char)slip;
    rtk->ssat[sat-1].half[f]=(obs[i].LLI[f]&2)?0:1;
}
/* all ambiguity reset--------------------------------------------------------*/
static void arc_ubbias_all(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat,
                           const int *iu, const int *ir, int ns, const nav_t *nav)
{
    double cp,pr,*bias,offset,lami;
    int i,j,f=0;

    arc_log(ARC_INFO,"arc_udbias  : tt=%.1f ns=%d\n",tt,ns);

    bias=arc_zeros(ns,1);

    /* estimate approximate phase-bias by phase - code */
    for (i=j=0,offset=0.0;i<ns;i++) {

        cp=arc_sdobs(obs,iu[i],ir[i],f); /* cycle */
        pr=arc_sdobs(obs,iu[i],ir[i],f+NFREQ);
        lami=nav->lam[sat[i]-1][f];
        if (cp==0.0||pr==0.0||lami<=0.0) continue;

        bias[i]=cp-pr/lami;

        if (rtk->x[IB(sat[i],f,&rtk->opt)]!=0.0) {
            offset+=bias[i]-rtk->x[IB(sat[i],f,&rtk->opt)];
            j++;
        }
    }
    /* set initial states of phase-bias */
    for (i=0;i<ns;i++) {
        arc_initx(rtk,bias[i],SQR(rtk->opt.std[0]),IB(sat[i],f,&rtk->opt));
    }
    free(bias);
}
/* temporal update of phase biases -------------------------------------------*/
static void arc_udbias(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav,
                       const double *rs)
{
    double cp,pr,*bias,offset,lami;
    int i,j,f=0,slip,reset;

    arc_log(ARC_INFO,"arc_udbias  : tt=%.1f ns=%d\n",tt,ns);

    if (rtk->opt.reset_amb_all) {
        return arc_ubbias_all(rtk,tt,obs,sat,iu,ir,ns,nav);  /* todo:this method of ambiguity initial is bad */
    }
    for (i=0;i<ns;i++) {

        /* detect cycle slip by LLI */
        rtk->ssat[sat[i]-1].slip[f]&=0xFC;
        arc_detslp_ll(rtk,obs,iu[i],1);
        arc_detslp_ll(rtk,obs,ir[i],2);

        /* update half-cycle valid flag */
        rtk->ssat[sat[i]-1].half[f]=
                !((obs[iu[i]].LLI[f]&2)||(obs[ir[i]].LLI[f]&2));
    }
    arc_detsl_new(rtk,obs,iu,ir,ns,nav,rs);  /* detect cycle slip by new ways */

    /* reset phase-bias if instantaneous AR or expire obs outage counter */
    for (i=1;i<=MAXSAT;i++) {

        reset=++rtk->ssat[i-1].outc[0]>(unsigned int)rtk->opt.maxout;
        if (rtk->opt.modear==ARMODE_INST&&rtk->x[IB(i,0,&rtk->opt)]!=0.0) {
            arc_initx(rtk,0.0,VAR_AMB,IB(i,0,&rtk->opt));  /* each epoch is re initialized */
        }
        else if (reset&&rtk->x[IB(i,0,&rtk->opt)]!=0.0) {  /* reset ambiguity */
            arc_initx(rtk,0.0,VAR_AMB,IB(i,0,&rtk->opt));
            arc_log(ARC_INFO, "arc_udbias : obs outage counter overflow (sat=%3d L%d n=%d)\n",
                    i,f+1,rtk->ssat[i-1].outc[0]);
        }
        if (rtk->opt.modear!=ARMODE_INST&&reset) {
            rtk->ssat[i-1].lock[0]=-rtk->opt.minlock;
        }
    }
    /* reset phase-bias if detecting cycle slip */
    for (i=0;i<ns;i++) {
        j=IB(sat[i],0,&rtk->opt);
        rtk->P[j+j*rtk->nx]+=rtk->opt.prn[0]*rtk->opt.prn[0]*tt;
        slip=rtk->ssat[sat[i]-1].slip[0];
        if (rtk->opt.modear==ARMODE_INST||!(slip&1)) continue;
        rtk->x[j]=0.0;
        rtk->ssat[sat[i]-1].lock[0]=-rtk->opt.minlock;
    }
    bias=arc_zeros(ns,1);

    /* estimate approximate phase-bias by phase - code */
    for (i=j=0,offset=0.0;i<ns;i++) {

        cp=arc_sdobs(obs,iu[i],ir[i],0); /* cycle */
        pr=arc_sdobs(obs,iu[i],ir[i],0+NFREQ);
        lami=nav->lam[sat[i]-1][0];
        if (cp==0.0||pr==0.0||lami<=0.0) continue;

        /* single-difference ambiguity */
        bias[i]=cp-pr/lami;

        /* if is ARMODE_INST,then here it no process */
        if (rtk->x[IB(sat[i],0,&rtk->opt)]!=0.0) {
            offset+=bias[i]-rtk->x[IB(sat[i],0,&rtk->opt)];
            j++;
        }
    }
    /* correct phase-bias offset to enssure phase-code coherency */
    /* todo:this step may have problems,need to think */
    if (j>0) {
        for (i=1;i<=MAXSAT;i++) {
            if (rtk->x[IB(i,0,&rtk->opt)]!=0.0) rtk->x[IB(i,0,&rtk->opt)]+=offset/j;
        }
    }
    /* set initial states of phase-bias (ARMODE_INST) */
    for (i=0;i<ns;i++) {
        if (bias[i]==0.0||rtk->x[IB(sat[i],0,&rtk->opt)]!=0.0) continue;
        arc_initx(rtk,bias[i],SQR(rtk->opt.std[0]),IB(sat[i],0,&rtk->opt));
    }
    free(bias);
}
/* temporal update of states --------------------------------------------------*/
static void arc_udstate(rtk_t *rtk, const obsd_t *obs, const int *sat,
                        const int *iu, const int *ir, int ns, const nav_t *nav,
                        const double *rs)
{
    double tt=fabs(rtk->tt),bl,dr[3];

    arc_log(ARC_INFO, "arc_udstate : ns=%d\n", ns);

    /* temporal update of position/velocity/acceleration */
    arc_udpos(rtk,tt);

    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt>=IONOOPT_EST) {
        bl=arc_baseline(rtk->x,rtk->rb,dr);
        arc_udion(rtk,tt,bl,sat,ns);
    }
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt>=TROPOPT_EST) {
        arc_udtrop(rtk,tt,bl);
    }
    /* temporal update of phase-bias */
    if (rtk->opt.mode>PMODE_DGPS) {
        arc_udbias(rtk,tt,obs,sat,iu,ir,ns,nav,rs);
    }
}
/* undifferenced phase/code residual for satellite ---------------------------*/
static void arc_zdres_sat(int base, double r, const obsd_t *obs, const nav_t *nav,
                          const double *azel, const double *dant,double dion,
                          double vion,const prcopt_t *opt, double *y,const rtk_t *rtk,
                          double *ukf_y,int *nzd)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=0,nf=1;

    if (lam[i]==0.0) return;

    /* check snr mask */
    if (testsnr(base,i,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) return;

    /* residuals = observable - pseudorange,think about phase windup correction */
    if (y) {
        if (obs->L[i]!=0.0) y[i   ]=obs->L[i]*lam[i]-r-dant[i]+dion-rtk->ssat[obs->sat-1].phw*lam[i];
        if (obs->P[i]!=0.0) y[i+nf]=obs->P[i]       -r-dant[i]-dion;
    }
    /* save sigma point measurements to ukf strcut */
    if (ukf_y) {
        ukf_y[i   ]=r+dant[i]+rtk->ssat[obs->sat-1].phw*lam[i];
        ukf_y[i+nf]=r+dant[i];
    }
    /* record undifferenced residual numbers */
    if (nzd) *nzd+=2;
}
/* undifferenced phase/code residuals ----------------------------------------*/
static int arc_zdres(int base, const obsd_t *obs, int n, const double *rs,
                     const double *dts,const int *svh, const nav_t *nav,
                     const double *rr,const prcopt_t *opt, int index, double *y,
                     double *e,double *azel,rtk_t* rtk,double *ukf_y)
{
    double r,rr_[3],pos[3],dant[NFREQ]={0},disp[3];
    double zhd,zazel[]={0.0,90.0*D2R},dion,vion,*py=y,*pukfy=ukf_y;
    int i=0,nf=1,nzd=0;

    arc_log(ARC_INFO, "arc_zdres   : n=%d\n",n);

    /* reset undifferenced phase/code observation and residuals */
    if (y) for (i=0;i<n*nf*2;i++) y[i]=0.0;
    if (ukf_y) for (i=0;i<n*nf*2;i++) ukf_y[i]=0.0;

    if (arc_norm(rr,3)<=0.0) return 0; /* no receiver position */

    for (i=0;i<3;i++) rr_[i]=rr[i];

    /* earth tide correction */
    if (opt->tidecorr) {
        arc_tidedisp(gpst2utc(obs[0].time),rr_,opt->tidecorr,&nav->erp,
                     opt->odisp[base],disp);
        for (i=0;i<3;i++) rr_[i]+=disp[i];
    }
    ecef2pos(rr_,pos);

    for (i=0;i<n;i++) {
        /* compute geometric-range and azimuth/elevation angle */
        if ((r=arc_geodist(rs+i*6,rr_,e+i*3))<=0.0) continue;
        if (arc_satazel(pos,e+i*3,azel+i*2)<opt->elmin) continue;

        /* excluded satellite? */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;

        /* satellite clock-bias */
        r+=-CLIGHT*dts[i*2];

        /* troposphere delay model (hydrostatic) */
        zhd=arc_tropmodel(obs[0].time,pos,zazel,0.0);
        r+=arc_tropmapf(obs[i].time,pos,azel+i*2,NULL)*zhd;

        /* ionospheric corrections */
        if (!arc_ionocorr(obs[i].time,nav,obs[i].sat,pos,azel+i*2,
                          IONOOPT_BRDC,&dion,&vion)) continue;
        /* receiver antenna phase center correction */
        arc_antmodel(opt->pcvr+index,opt->antdel[index],azel+i*2,opt->posopt[1],dant);

        /* phase windup model */
        if (!arc_model_phw(rtk->sol.time,obs[i].sat,nav->pcvs[obs[i].sat-1].type,
                           opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[obs[i].sat-1].phw)) {
            continue;
        }
        /* todo:here may have some more better methods to do this */
        /* adjust measurements vector pointor */
        if (py!=NULL) py=y+i*nf*2;           /* for akf and ekf */

        if (pukfy!=NULL) pukfy=ukf_y+i*nf*2; /* for ukf */

        /* undifferenced phase/code residual for satellite */
        arc_zdres_sat(base,r,obs+i,nav,azel+i*2,dant,dion,vion,opt,py,rtk,pukfy,&nzd);
    }
    arc_log(ARC_INFO,"arc_zdres : rr_=%.3f %.3f %.3f\n",rr_[0],rr_[1],rr_[2]);
    arc_log(ARC_INFO,"arc_zdres : pos=%.9f %.9f %.3f\n",pos[0]*R2D,pos[1]*R2D,pos[2]);
    for (i=0;i<n;i++) {
        arc_log(ARC_INFO,"arc_zdres : sat=%2d %13.3f %13.3f %13.3f %13.10f %6.1f %5.1f\n",
                obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],dts[i*2],azel[i*2]*R2D,
                azel[1+i*2]*R2D);
    }
    if (y) {
        arc_log(ARC_INFO, "arc_zdres : y=\n");
        arc_tracemat(ARC_MATPRINTF,y,nf*2,n,13,3);
    }
    if (ukf_y) {
        arc_log(ARC_INFO, "arc_zdres : ukf y=\n");
        arc_tracemat(ARC_MATPRINTF,ukf_y,nf*2,n,13,3);
    }
    return nzd;
}
/* test valid observation data -----------------------------------------------*/
static int arc_validobs(int i,int j,int f,int nf,double *y)
{
    /* if no phase observable, psudorange is also unusable */
    return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
           (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
/* double-differenced measurement error covariance ---------------------------*/
static void arc_ddcov(const int *nb,int n,const double *Ri,const double *Rj,
                      int nv,double *R)
{
    int i,j,k=0,b=0;

    arc_log(ARC_INFO,"arc_ddcov   : n=%d\n",n);

    if (nv<=0) return; /* no double-difference measurements */

    for (i=0;i<nv*nv;i++) R[i]=0.0;
    for (b=0;b<n;k+=nb[b++]) {
        for (i=0;i<nb[b];i++) for (j=0;j<nb[b];j++) {
            R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
        }
    }
    arc_log(ARC_INFO,"R=\n");
    arc_tracemat(5,R,nv,nv,8,6);
}
/* precise tropspheric model -------------------------------------------------*/
static double arc_prectrop(gtime_t time, const double *pos, int r,
                           const double *azel, const prcopt_t *opt, const double *x,
                           double *dtdx)
{
    double m_w=0.0,cotz,grad_n,grad_e;
    int i=IT(r,opt);

    /* wet mapping function */
    arc_tropmapf(time,pos,azel,&m_w);

    if (opt->tropopt>=TROPOPT_ESTG&&azel[1]>0.0) {

        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[i+1]+grad_e*x[i+2];
        dtdx[1]=grad_n*x[i];
        dtdx[2]=grad_e*x[i];
    }
    else dtdx[1]=dtdx[2]=0.0;
    dtdx[0]=m_w;
    return m_w*x[i];
}
/* test navi system (m=0:gps/qzs/sbs,1:glo,2:gal,3:bds) ----------------------*/
static int arc_test_sys(int sys, int m)
{
    switch (sys) {
        case SYS_GPS: return m==0;
        case SYS_SBS: return m==0;
        case SYS_GAL: return m==2;
        case SYS_CMP: return m==3;
    }
    return 0;
}
/* given sat no and return its observation index------------------------------*/
static int arc_ukf_get_dd_obsind(const rtk_t* rtk,int *ir, int*iu,const obsd_t*obs,
                                 int rsat,int usat)
{
    arc_log(ARC_INFO,"arc_ukf_get_dd_obsind : \n");

    int i,j;
    for (i=0,j=0;i<_UKF_NR_+_UKF_NU_;i++) if (obs[i].sat==rsat) ir[j++]=i;
    for (i=0,j=0;i<_UKF_NR_+_UKF_NU_;i++) if (obs[i].sat==usat) iu[j++]=i;
    if (j==0) return 1;  /* no double-difference measurements */
    return 0;            /* success,return double-difference measurements pair */
}
/* ukf:double-difference phase/code measurements------------------------------*/
static int arc_ukf_ddmeas(rtk_t* rtk,const nav_t* nav,double dt,double *x,const double* P,
                          int *sat, int*iu, int*ir,int ns,double *azel,double *y,
                          double *ddy)
{
    prcopt_t *opt=&rtk->opt;
    double posu[3],posr[3],didxi=0.0,didxj=0.0,*im;
    double *tropr,*tropu,*dtdxr,*dtdxu,lami,lamj,fi,fj;
    int i,j,m,f=0,nv=0,nf=1,n=0,sati,satj;

    arc_log(ARC_INFO,"arc_ukf_ddmeas   : dt=%.1f nx=%d ns=%d\n", dt,rtk->nx,ns);
    arc_log(ARC_INFO,"arc_ukf_ddmeas   : np=%4d  nc=%4d \n",rtk->np,rtk->nc);

    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);

    im=arc_mat(ns,1);
    tropu=arc_mat(ns,1); tropr=arc_mat(ns,1);
    dtdxu=arc_mat(ns,3); dtdxr=arc_mat(ns,3);

    /* reset double-difference measurements for ukf */
    if (ddy) for (i=0;i<ns;i++) ddy[i]=0.0;

    /* compute factors of ionospheric and tropospheric delay */
    for (i=0;i<ns;i++) {
        if (opt->ionoopt>=IONOOPT_EST) {
            im[i]=(arc_ionmapf(posu,azel+iu[i]*2)+arc_ionmapf(posr,azel+ir[i]*2))/2.0;
        }
        if (opt->tropopt>=TROPOPT_EST) {
            tropu[i]=arc_prectrop(rtk->sol.time,posu,0,azel+iu[i]*2,opt,x,dtdxu+i*3);
            tropr[i]=arc_prectrop(rtk->sol.time,posr,1,azel+ir[i]*2,opt,x,dtdxr+i*3);
        }
    }
    for (m=0;m<2;m++) {  /* loop for phase and code */

        /* numbers of double-difference phase observations */
        if      (m==0) n=rtk->nc;
        /* numbers of double-difference code observations */
        else if (m==1) n=rtk->np;

        for (j=0;j<n;j++) {  /* loop for numbers of phase and code observations */

            sati=rtk->sat[2*nv  ];
            satj=rtk->sat[2*nv+1];

            lami=nav->lam[sati-1][0];
            lamj=nav->lam[satj-1][0];

            /* ukf sigma point measurement */
            if (y&&ddy) {
                ddy[nv]=(y[f+iu[rtk->obs_ind[2*nv  ]]*nf*2]-y[f+ir[rtk->obs_ind[2*nv  ]]*nf*2])-
                        (y[f+iu[rtk->obs_ind[2*nv+1]]*nf*2]-y[f+ir[rtk->obs_ind[2*nv+1]]*nf*2]);
            }
            /* double-differenced ionospheric delay term */
            if (opt->ionoopt==IONOOPT_EST) {
                fi=lami/lam_carr[0]; fj=lamj/lam_carr[0];
                didxi=(m==0?-1.0:1.0)*fi*fi*im[i];
                didxj=(m==0?-1.0:1.0)*fj*fj*im[j];

                /* compute the ukf sigma point measurements */
                if (ddy) ddy[nv]+=didxi*x[II(sati,opt)]-didxj*x[II(satj,opt)];
            }
            /* double-differenced tropospheric delay term */
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                /* compute the ukf sigma point measurements */
                if (ddy&&y) {
                    ddy[nv]+=(tropu[i]-tropu[j])-(tropr[i]-tropr[j]);
                }
            }
            /* double-differenced phase-bias term */
            if (m==0) {
                /* compute the ukf sigma point measurements */
                if (y&&ddy) {
                    ddy[nv]+=lami*x[IB(sati,f,opt)]-lamj*x[IB(satj,f,opt)];
                }
            }
            arc_log(ARC_INFO,"arc_ukf_ddmeas : sat=%3d-%3d %s%d  ddy=%14.4lf \n",
                    sati,satj,m==0?"L":"P",1,ddy[nv]);
            nv++;  /* increase double-residual index */
        }
    }
    free(im);free(tropu); free(tropr);
    free(dtdxu); free(dtdxr);
    return nv;
}
/* exctract actcive H of kalman filter----------------------------------------*/
static int arc_kalman_exct_xHP(const rtk_t *rtk,const double *Hi,double *Ho,
                                int nv,double *x,double *P)
{
    int i,j,k,nx=rtk->nx,*ix;

    arc_log(ARC_INFO,"arc_kalman_exct_H:\n");

    ix=arc_imat(nx,1); for (i=0,j=0;i<nx;i++) if (rtk->ceres_active_x[i]) ix[j++]=i;

    arc_log(ARC_INFO,"active sates index =\n");
    arc_tracemati(ARC_MATPRINTF,ix,1,j,2,0);

    if (x) for (i=0;i<j;i++) x[i]=rtk->x[ix[i]];
    if (P) for (i=0;i<j;i++) for (k=0;k<j;k++) P[i*j+k]=rtk->P[ix[i]*nx+ix[k]];

    if (Ho) for (i=0;i<nv;i++) for (k=0;k<j;k++) Ho[i*j+k]=Hi[i*nx+ix[k]];

    arc_log(ARC_INFO,"x=\n"); arc_tracemat(ARC_MATPRINTF,x,1,j,10,4);
    arc_log(ARC_INFO,"P=\n"); arc_tracemat(ARC_MATPRINTF,P,j,j,10,4);
    arc_log(ARC_INFO,"H=\n"); arc_tracemat(ARC_MATPRINTF,Ho,j,nv,10,4);

    free(ix);
    return j; /* numbers of active states */
}
/* normalization of innovation series of kalman filter------------------------*/
static int arc_kalman_norm_inno(const rtk_t *rtk,double *v,int nv,const double *H,
                                const double *R)
{
    int i,nx=rtk->nx,n,info=1;
    double *He,*Pe,*HP,*Qe,*L=NULL,*ve;

    arc_log(ARC_INFO,"arc_kalman_norm_inno :\n");

    He=arc_mat(nx,nv); Pe=arc_zeros(nx,nx);

    n=arc_kalman_exct_xHP(rtk,H,He,nv,NULL,Pe);

    HP=arc_mat(nv,n); Qe=arc_mat(nv,nv); ve=arc_mat(nv,1);
    arc_matmul("TN",nv,n,n,1.0,He,Pe,0.0,HP);
    arc_matcpy(Qe,R,nv,nv);
    arc_matmul("NN",nv,nv,n,1.0,HP,H,1.0,Qe);

    arc_log(ARC_INFO,"Qe=\n"); arc_tracemat(ARC_MATPRINTF,Qe,nv,nv,10,4);

    if (!(info=arc_matinv(Qe,nv))) {

        /* cholesky matrix decomposition */
        L=arc_cholesky(Qe,nv);

        arc_log(ARC_INFO,"L=\n"); arc_tracemat(ARC_MATPRINTF,L,nv,nv,10,4);

        /* normalization */
        if (L) arc_matmul("NN",nv,1,nv,1.0,L,v,0.0,ve);

        arc_log(ARC_INFO,"ve=\n"); arc_tracemat(ARC_MATPRINTF,ve,1,nv,10,4);
    }
    if (L) free(L); free(He); free(Pe); free(ve); free(Qe);
    return info;
}
/* covariance matrix of normalization of innovation series of kalman filte----*/
static void arc_kalman_norm_Qino(const double *v,int nv,double *Qv)
{
    arc_log(ARC_INFO,"arc_kalman_norm_Qino :\n");

    int i;
    double ave_v=0.0,*ve=arc_mat(nv,1);
    for (i=0;i<nv;i++) ave_v+=v[i]; ave_v/=nv;

    arc_log(ARC_INFO,"average of innovation series=%10.6lf\n",ave_v);

    for (i=0;i<nv;i++) ve[i]=v[i]-ave_v;
    if (Qv) arc_matmul("NT",nv,nv,1,1.0,ve,ve,0.0,Qv);

    arc_log(ARC_INFO,"Qv=\n"); arc_tracemat(ARC_MATPRINTF,Qv,nv,nv,10,4);

    free(ve);
}
/* kalman robust check function-----------------------------------------------*/
static double arc_robust_chk(int nv,double alpha,double rk)
{
    if (alpha<=0.0) {
        if (rk>=chisqr[nv]) return SQRT(chisqr[nv]/rk); /* alpha=0.001 */
    }
    else {
        if (rk>=arc_re_chi2(nv,alpha)) return SQRT(chisqr[nv]/rk);
    }
    return 1.0; /* todo:need to test */
}
/* kalman filter Phi matrix---------------------------------------------------*/
static void arc_kalman_robust_phi(const rtk_t *rtk,const double *vn,int nv,
                                  double *phi) /* todo:need to test */
{
    int i;
    double alpha=rtk->opt.kalman_robust_alpha;

    arc_log(ARC_INFO,"arc_kalman_phi :\n");

    for (i=0;i<nv;i++) phi[i+i*nv]=arc_robust_chk(nv,alpha,fabs(vn[i]));

    arc_log(ARC_INFO,"Phi=\n");
    arc_tracemat(ARC_MATPRINTF,phi,nv,nv,10,4);
}
/* double-differenced phase/code residuals -----------------------------------*/
static int arc_ddres(rtk_t *rtk,const nav_t *nav,double dt,const double *x,
                     const double *P,const int *sat,double *y,double *e,
                     double *azel,const int *iu,const int *ir,int ns,double *v,
                     double *H,double *R,int *vflg)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3],didxi=0.0,didxj=0.0,*im;
    double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,lami,lamj,fi,fj,*Hi=NULL;
    int i,j,k,m,f,ff=0,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=1;

    arc_log(ARC_INFO, "arc_ddres   : dt=%.1f nx=%d ns=%d\n", dt,rtk->nx,ns);

    bl=arc_baseline(x,rtk->rb,dr);
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);

    Ri=arc_mat(ns*nf*2+2,1); Rj=arc_mat(ns*nf*2+2,1); im=arc_mat(ns,1);
    tropu=arc_mat(ns,1); tropr=arc_mat(ns,1); dtdxu=arc_mat(ns,3); dtdxr=arc_mat(ns,3);

    /* initial satellite status informations */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].resp[0]=rtk->ssat[i].resc[0]=0.0;
    }
    /* reset hold the phase and pseudorange observations numbers */
    rtk->nc=rtk->np=0;

    /* compute factors of ionospheric and tropospheric delay */
    for (i=0;i<ns;i++) {
        if (opt->ionoopt>=IONOOPT_EST) {
            im[i]=(arc_ionmapf(posu,azel+iu[i]*2)+arc_ionmapf(posr,azel+ir[i]*2))/2.0;
        }
        if (opt->tropopt>=TROPOPT_EST) {
            tropu[i]=arc_prectrop(rtk->sol.time,posu,0,azel+iu[i]*2,opt,x,dtdxu+i*3);
            tropr[i]=arc_prectrop(rtk->sol.time,posr,1,azel+ir[i]*2,opt,x,dtdxr+i*3);
        }
    }
    for (m=0;m<4;m++) /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

    for (f=opt->mode>PMODE_DGPS?0:nf;f<nf*2;f++) {

        /* search reference satellite with highest elevation */
        for (i=-1,j=0;j<ns;j++) {
            sysi=rtk->ssat[sat[j]-1].sys;
            if (!arc_test_sys(sysi,m)) continue;

            /* y[i] may be zero,so must check out,this step is importance */
            if (y) if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;

            /* set the reference satellite index */
            if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
        }
        if (i<0) continue;  /* i is the reference satellite */

        /* make double difference */
        for (j=0;j<ns;j++) {
            if (i==j) continue;  /* sat[i] is reference satellite */
            sysi=rtk->ssat[sat[i]-1].sys;
            sysj=rtk->ssat[sat[j]-1].sys;
            if (!arc_test_sys(sysj,m)) continue;

            /* todo:may be have bugs on running,need to fix in future */
            if (y) if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;

            ff=f%nf;
            lami=nav->lam[sat[i]-1][ff];
            lamj=nav->lam[sat[j]-1][ff];
            if (lami<=0.0||lamj<=0.0) continue;
            if (H) {
                Hi=H+nv*rtk->nx;
                for (k=0;k<rtk->nx;k++) Hi[k]=0.0;
            }
            /* double-differenced residual */
            if (v) v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])
                        -(y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);

            /* partial derivatives by rover position */
            if (H) {
                for (k=0;k<3;k++) {
                    Hi[k]=-e[k+iu[i]*3]+e[k+iu[j]*3];
                }
            }
            /* double-differenced ionospheric delay term */
            if (opt->ionoopt==IONOOPT_EST) {
                fi=lami/lam_carr[0]; fj=lamj/lam_carr[0];
                didxi=(f<nf?-1.0:1.0)*fi*fi*im[i];
                didxj=(f<nf?-1.0:1.0)*fj*fj*im[j];
                if (v) v[nv]-=didxi*x[II(sat[i],opt)]-didxj*x[II(sat[j],opt)];

                /* design matrix */
                if (H) {
                    Hi[II(sat[i],opt)]= didxi;
                    Hi[II(sat[j],opt)]=-didxj;
                }
            }
            /* double-differenced tropospheric delay term */
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                if (v) v[nv]-=(tropu[i]-tropu[j])-(tropr[i]-tropr[j]);
                for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {

                    if (!H) continue;
                    /* design matrix of double-differenced tropospheric delay term */
                    Hi[IT(0,opt)+k]= (dtdxu[k+i*3]-dtdxu[k+j*3]);
                    Hi[IT(1,opt)+k]=-(dtdxr[k+i*3]-dtdxr[k+j*3]);
                }
            }
            /* double-differenced phase-bias term */
            if (f<nf) {
                if (v) v[nv]-=lami*x[IB(sat[i],f,opt)]-lamj*x[IB(sat[j],f,opt)];
                /* design matrix */
                if (H) {
                    Hi[IB(sat[i],f,opt)]= lami;
                    Hi[IB(sat[j],f,opt)]=-lamj;
                }
            }
            if (v) {
                if (f<nf) rtk->ssat[sat[j]-1].resc[f   ]=v[nv];
                else      rtk->ssat[sat[j]-1].resp[f-nf]=v[nv];
            }
            /* test innovation */
            if (v) if (opt->maxinno>0.0  /* todo:test innovation need to modify,have bugs for ukf */
                       &&fabs(v[nv])>opt->maxinno) {
                if (f<nf) {
                    rtk->ssat[sat[i]-1].rejc[f]++;
                    rtk->ssat[sat[j]-1].rejc[f]++;
                }
                arc_log(ARC_WARNING, "arc_ddres : outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                        sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
                continue;
            }
            /* save double-difference satellite pair index and its observation index*/
            rtk->sat[nv*2]=sat[i]; rtk->sat[nv*2+1]=sat[j];
            rtk->obs_ind[2*nv]=i;  rtk->obs_ind[2*nv+1]=j;

            /* hold the active states index */
            if (opt->ionoopt==IONOOPT_EST) {
                /* updates ceres solver active states index list */
                rtk->ceres_active_x[II(sat[i],opt)]=1;  /* this is important */
                rtk->ceres_active_x[II(sat[j],opt)]=1;
            }
            /* double-differenced tropospheric delay term */
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {
                    /* updates ceres solver active states index list */
                    rtk->ceres_active_x[IT(0,opt)+k]=1;
                    rtk->ceres_active_x[IT(1,opt)+k]=1;
                }
            }
            /* double-differenced phase-bias term */
            if (f<nf) {
                /* updates ceres solver active states index list */
                rtk->ceres_active_x[IB(sat[i],f,opt)]=1;
                rtk->ceres_active_x[IB(sat[j],f,opt)]=1;
                /* hold phase observation numbers */
                rtk->nc++;
            }
            else rtk->np++;  /* hold code observation numbers */

            /* single-differenced measurement error variances */
            if (R) {
                Rj[nv]=arc_varerr(sat[j],sysj,azel[1+iu[j]*2],bl,dt,f,opt);
                Ri[nv]=arc_varerr(sat[i],sysi,azel[1+iu[i]*2],bl,dt,f,opt);
            }
            /* set valid data flags */
            if (opt->mode>PMODE_DGPS) {
                if (f<nf) rtk->ssat[sat[i]-1].vsat[f]=rtk->ssat[sat[j]-1].vsat[f]=1;
            }
            else {
                rtk->ssat[sat[i]-1].vsat[f-nf]=rtk->ssat[sat[j]-1].vsat[f-nf]=1;
            }
            arc_log(ARC_INFO,"arc_ddres : sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n",sat[i],
                    sat[j],f<nf?"L":"P",f%nf+1,v[nv],R==NULL?-999.0:Ri[nv],R==NULL?-999.0:Rj[nv]);
            if (vflg) vflg[nv]=(sat[i]<<16)|(sat[j]<<8)|((f<nf?0:1)<<4)|(f%nf);
            nb[b]++; nv++;  /* increase double-residual index */
        }
        b++;
    }
    /* end of system loop */

    if (H&&nv) {
        arc_log(ARC_INFO, "arc_ddres : H=\n"); arc_tracemat(ARC_MATPRINTF, H,rtk->nx,nv,7,4);
    }
    /* double-differenced measurement error covariance */
    if (R&&nv) {
        arc_ddcov(nb,b,Ri,Rj,nv,R);
    }
    arc_log(ARC_INFO,"arc_ddres : active states index = \n");
    arc_tracemati(ARC_MATPRINTF,rtk->ceres_active_x,1,rtk->nx,2,1);

    free(Ri); free(Rj); free(im);
    free(tropu); free(tropr); free(dtdxu); free(dtdxr);
    return nv;
}
/* difference pseudorange for initialing--------------------------------------*/
static int arc_diff_pr_ddres(rtk_t *rtk,const nav_t *nav,double dt,const double *x,
                             const double *P,const int *sat,double *y,double *e,
                             double *azel,const int *iu,const int *ir,int ns,double *v,
                             double *H,double *R,int *vflg)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3];
    double *Ri,*Rj,lami,lamj,*Hi=NULL;
    int i,j,k,m,f,ff=0,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=1,nx=NP(&rtk->opt);

    arc_log(ARC_INFO,"arc_init_dc_res   : dt=%.1f nx=%d ns=%d\n",dt,nx,ns);

    bl=arc_baseline(x,rtk->rb,dr);
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);

    Ri=arc_mat(ns*nf*2+2,1); Rj=arc_mat(ns*nf*2+2,1);

    for (m=0;m<NUMOFSYS;m++) /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        for (f=nf;f<nf*2;f++) {

            /* search reference satellite with highest elevation */
            for (i=-1,j=0;j<ns;j++) {
                sysi=rtk->ssat[sat[j]-1].sys;
                if (!arc_test_sys(sysi,m)) continue;

                /* y[i] may be zero,so must check out,this step is importance */
                if (y) if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;

                /* set the reference satellite index */
                if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
            }
            if (i<0) continue;  /* i is the reference satellite */

            /* make double difference */
            for (j=0;j<ns;j++) {
                if (i==j) continue;  /* sat[i] is reference satellite */
                sysi=rtk->ssat[sat[i]-1].sys;
                sysj=rtk->ssat[sat[j]-1].sys;
                if (!arc_test_sys(sysj,m)) continue;

                /* todo:may be have bugs on running,need to fix in future */
                if (y) if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;

                ff=f%nf;
                lami=nav->lam[sat[i]-1][ff];
                lamj=nav->lam[sat[j]-1][ff];
                if (lami<=0.0||lamj<=0.0) continue;
                if (H) {
                    Hi=H+nv*nx;
                    for (k=0;k<nx;k++) Hi[k]=0.0;
                }
                /* double-differenced residual */
                if (v) v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])
                             -(y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);

                /* partial derivatives by rover position */
                if (H) {
                    for (k=0;k<3;k++) {
                        Hi[k]=-e[k+iu[i]*3]+e[k+iu[j]*3];
                    }
                }
                /* test innovation */
                if (v) if (opt->maxinno>0.0  /* todo:test innovation need to modify,have bugs for ukf */
                           &&fabs(v[nv])>opt->maxinno) {
                    arc_log(ARC_WARNING,"arc_init_dc_res : outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                            sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
                    continue;
                }
                /* single-differenced measurement error variances */
                if (R) {
                    Rj[nv]=arc_varerr(sat[j],sysj,azel[1+iu[j]*2],bl,dt,f,opt);
                    Ri[nv]=arc_varerr(sat[i],sysi,azel[1+iu[i]*2],bl,dt,f,opt);
                }

                arc_log(ARC_INFO,"arc_init_dc_res : sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n",
                        sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv],
                        R==NULL?-999.0:Ri[nv],R==NULL?-999.0:Rj[nv]);

                if (vflg) vflg[nv]=(sat[i]<<16)|(sat[j]<<8)|((f<nf?0:1)<<4)|(f%nf);

                nb[b]++;
                nv++;  /* increase double-residual index */
            }
            b++;
        }
    /* end of system loop */

    if (H) {
        arc_log(ARC_INFO,"arc_init_dc_res : H=\n"); arc_tracemat(ARC_MATPRINTF,H,nx,nv,7,4);
    }
    /* double-differenced measurement error covariance */
    if (R) arc_ddcov(nb,b,Ri,Rj,nv,R);

    free(Ri); free(Rj);
    return nv;
}
/* time-interpolation of residuals (for post-mission) ------------------------*/
static double arc_intpres(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                          rtk_t *rtk, double *y)
{
    static obsd_t obsb[MAXOBS];
    static double yb[MAXOBS*NFREQ*2],rs[MAXOBS*6],dts[MAXOBS*2],var[MAXOBS];
    static double e[MAXOBS*3],azel[MAXOBS*2];
    static int nb=0,svh[MAXOBS*2];
    prcopt_t *opt=&rtk->opt;
    double tt=timediff(time,obs[0].time),ttb,*p,*q;
    int i,j,k,nf=1;

    arc_log(ARC_INFO, "arc_intpres : n=%d tt=%.1f\n", n, tt);

    if (nb==0||fabs(tt)<DTTOL) {
        nb=n; for (i=0;i<n;i++) obsb[i]=obs[i];
        return tt;
    }
    ttb=timediff(time,obsb[0].time);
    if (fabs(ttb)>opt->maxtdiff*2.0||ttb==tt) return tt;

    arc_satposs(time,obsb,nb,nav,opt->sateph,rs,dts,var,svh);

    if (!arc_zdres(1,obsb,nb,rs,dts,svh,nav,rtk->rb,opt,1,yb,e,azel,rtk,NULL)) {
        return tt;
    }
    for (i=0;i<n;i++) {
        for (j=0;j<nb;j++) if (obsb[j].sat==obs[i].sat) break;
        if (j>=nb) continue;
        for (k=0,p=y+i*nf*2,q=yb+j*nf*2;k<nf*2;k++,p++,q++) {
            if (*p==0.0||*q==0.0) *p=0.0;
            else *p=(ttb*(*p)-tt*(*q))/(ttb-tt);
        }
    }
    return fabs(ttb)>fabs(tt)?ttb:tt;
}
/* single to double-difference transformation matrix (D') --------------------*/
static int arc_ddmat(rtk_t *rtk, double *D)
{
    int i,j,k,m,f,nb=0,nx=rtk->nx,na=rtk->na,nf=1,nofix,ref=-1;
    double el=rtk->opt.elmaskar;

    arc_log(ARC_INFO, "arc_ddmat   :\n");

    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].fix[0]=0; /* initial fix flag,this step is very importance */
        rtk->ssat[i].group=0;  /* initial state,all satellite is un-groupping */
    }
    for (i=0;i<NUMOFSYS;i++) rtk->amb_refsat[i]=0; /* initial for reference satellite */

    /* initial single-difference to double-diffrence transformation matrix */
    if (D) for (i=0;i<na;i++) D[i+i*nx]=1.0;

    for (m=0;m<4;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        nofix=(m==3&&rtk->opt.bdsmodear==0); /* this flag is importance */
        for (f=0,k=na;f<nf;f++,k+=MAXSAT) {  /* loop for frequency */

            for (i=k;i<k+MAXSAT;i++) { /* loop for ambiguity list */
                /* select the reference single-difference ambiguity */
                if (rtk->x[i]==0.0||!rtk->ceres_active_x[i]||!arc_test_sys(rtk->ssat[i-k].sys,m)||
                    !rtk->ssat[i-k].vsat[f]||!rtk->ssat[i-k].half[f]) {
                    continue;
                }
                if (rtk->ssat[i-k].lock[f]>0&&!(rtk->ssat[i-k].slip[f]&2)&&
                    rtk->ssat[i-k].azel[1]>=el&&!nofix&&rtk->ceres_active_x[i]) {
                    el=rtk->ssat[i-k].azel[1]; ref=i; /* index of reference single-difference ambiguity */
                    continue;
                }
            }
            /* set the reference ambiguity index */
            rtk->ssat[ref-k].fix[f]=2;  /* fix */
            rtk->ssat[ref-k].group=1;   /* first group */
            rtk->amb_refsat[m]=ref-k+1; /* sat no.todo:this step need to modified */

            for (j=k;j<k+MAXSAT;j++) {
                if (ref==j||rtk->x[j]==0.0||!arc_test_sys(rtk->ssat[j-k].sys,m)||
                    !rtk->ssat[j-k].vsat[f]||!rtk->ceres_active_x[j]) {
                    continue;
                }
                if (rtk->ssat[j-k].lock[f]>0&&!(rtk->ssat[j-k].slip[f]&2)&&
                    rtk->ssat[ref-k].vsat[f]&&rtk->ceres_active_x[j]&&
                    rtk->ssat[j-k].azel[1]>=rtk->opt.elmaskar&&!nofix) {
                    if (D) D[ref+(na+nb)*nx]= 1.0;  /* reference single-difference ambiguity */
                    if (D) D[j  +(na+nb)*nx]=-1.0;  /* other single-difference ambiguity */
                    rtk->amb_index[nb]=j-k+1;       /* double-difference ambiguity index (sat no),no included reference satellite */
                    nb++;                           /* numbers of double-difference ambiguity */
                    rtk->ssat[j-k].fix[f]=2;        /* fix this ambiguity*/
                    rtk->ssat[j-k].group =1;        /* first group */
                }
            }
        }
    }
    rtk->amb_nb=nb; /* save numbers of double-difference ambiguity */
    if (D) {
        arc_log(ARC_INFO, "D=\n"); arc_tracemat(ARC_MATPRINTF,D,nx,na+nb,2,0);
    }
    return nb; /* numbers of double-difference ambiguity */
}
/* adop---------------------------------------------------------------------- */
static double arc_amb_adop(const double *Qb,const double *xb,int nb)
{
    double det=0.0; arc_matdet(Qb,nb,&det); return det; /* todo:need to test */
}
/* sucess probability value of ambiguity -------------------------------------*/
static double arc_amb_success(double adop,int nb)
{
    double pb=0.0; pb=arc_norm_distri(1.0/(2.0*adop));
    return pow(2.0*pb-1.0,nb); /* todo:need to test */
}
/* compute the bootstrapped success rate ------------------------------------*/
static double arc_amb_bs_success(const double *D,int n)
{
    double s=1.0;
    int i; for (i=0;i<n;i++) s*=(2.0*arc_norm_distri(0.5/SQRT(D[i]))-1.0);
    return s;
}
/* FF-ratio------------------------------------------------------------------*/
static int arc_FFRatio(int namb,double ffailure,double ratio,double pf) /* todo:need to test */
{
    int i,j,n;
    const double *tp=NULL;
    double ffratio=0.0;

    arc_log(ARC_INFO,"arc_FFRatio :\n");

    if (pf<=ffailure) return 1; if (namb<=0) return 0;

    if      (ffailure==0.01 ) n=41,tp=&ff_ratio_table1[0];
    else if (ffailure==0.001) n=64,tp=&ff_ratio_table2[0];

    namb=MIN(namb,n-1);

    /* search ffratio in tables */
    for (i=0;i<31;i++) if (tp[i*n+namb-1]>pf) break;

    /* interpolation of pf value */
    ffratio=(tp[i*n+namb]-tp[(i-1)*n+namb])
            /(tp[i*n]-tp[(i-1)*n])*(pf-tp[(i-1)*n])+tp[(i-1)*n+namb];

    arc_log(ARC_INFO,"FF-ratio=%8.4lf \n",ffratio);

    /* check ffratio */
    if (ffratio<=0.0) return 0;
    if (ratio>=(1.0/ffratio)) return 1; /* valid ok */

    return 0; /* valid failed */
}
/* restore single-differenced ambiguity --------------------------------------*/
static void arc_restamb(rtk_t *rtk, const double *bias, double *xa,int group)
{
    int i,j,n,m,index[MAXSAT],nv=0,ref=-1,refsat;

    arc_log(ARC_INFO,"arc_restamb :\n");

    for (i=0;i<rtk->nx;i++) xa[i]=rtk->x [i];  /* ambiguity */
    for (i=0;i<rtk->na;i++) xa[i]=rtk->xa[i];  /* station position/trp/iono/clock and so on */

    for (j=1;j<(group==0?2:3);j++)  /* loop for first and second group */

        for (m=0;m<NUMOFSYS;m++) {  /* m==3 is bds */

        /* reference single-difference ambiguity */
        refsat=(group==0?rtk->amb_refsat[m]:
                         rtk->amb_group_refsat[m][j-1]);

        if (refsat<=0) continue; /* no double-difference satellites */

        arc_log(ARC_INFO,"m=%d,reference single-difference satellites : %d",m,refsat);

        for (n=i=0;i<MAXSAT;i++) {
            if (!arc_test_sys(rtk->ssat[i].sys,m)
                ||rtk->ssat[i].fix[0]!=2
                ||rtk->ssat[i].group!=j /* j=1 is first group,j=2 is second group */
                ||(refsat==i+1)) {
                continue; /* no included reference satellite */
            }
            index[n++]=IB(i+1,0,&rtk->opt); /* include reference single-differnce ambiguity */
        }
        if (n<1) continue;
        ref=IB(refsat,0,&rtk->opt); /* reference single-differnce ambiguity */

        if (arc_conffunc(ROUND(xa[ref]),xa[ref],
                         SQRT(rtk->P[ref+ref*rtk->nx]))>=rtk->opt.amb_ref_thres) {
            arc_log(ARC_WARNING,"arc_restamb : reference single-difference "
                    "ambiguity is round(sat=%4d)",refsat);
            xa[ref]=ROUND(xa[ref]);
        }
        for (i=0;i<n;i++) xa[index[i]]=xa[ref]-bias[nv++];
    }
}
/* hold integer ambiguity ----------------------------------------------------*/
static void arc_holdamb(rtk_t *rtk, const double *xa,int group)
{
    double *v,*H,*R;
    int i,j,n,m,info,index[MAXSAT],nb=rtk->nx-rtk->na,nv=0,ref=-1,refsat=0;

    arc_log(ARC_INFO,"arc_holdamb :\n");

    v=arc_mat(nb,1); H=arc_zeros(nb,rtk->nx);

    for (j=1;j<(group==0?2:3);j++) /* two group double-difference ambiguity */

    for (m=0;m<NUMOFSYS;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        /* reference single-difference satellite */
        refsat=group==0?rtk->amb_refsat[m]:
                        rtk->amb_group_refsat[m][j];

        /* no double-difference satellites */
        if (refsat<=0) continue;

        if (group) {
            arc_log(ARC_INFO,"group(j=%4d,m=%4d,refsat=%4d)",j,m,refsat);
        }
        for (n=i=0;i<MAXSAT;i++) {
            if (!arc_test_sys(rtk->ssat[i].sys,m)
                ||rtk->ssat[i].fix[0]!=2
                ||rtk->ssat[i].azel[1]<rtk->opt.elmaskhold
                ||rtk->ssat[i].group!=j
                ||(i+1==refsat)) { /* exclued reference single-difference ambiguity */
                continue;
            }
            index[n++]=IB(i+1,0,&rtk->opt); /* index of fixing single-difference ambiguity */
            rtk->ssat[i].fix[0]=3; /* hold this single-difference ambiguity */
        }
        /* states index of reference single-differnce ambiguity */
        ref=IB(refsat,0,&rtk->opt);

        /* constraint to fixed ambiguity */
        for (i=1;i<n;i++) {
            v[nv]=(xa[ref]-xa[index[i]])-(rtk->x[ref]-rtk->x[index[i]]);
            if (fabs(v[nv])>=AMB_THRES) continue;
            H[ref     +nv*rtk->nx]= 1.0; /* reference ambiguity */
            H[index[i]+nv*rtk->nx]=-1.0; /* fix this single-difference ambiguity */
            nv++;
        }
    }
    if (nv>0) {
        R=arc_zeros(nv,nv);
        for (i=0;i<nv;i++) R[i+i*nv]=VAR_HOLDAMB;

        /* update states with constraints */
        if ((info=arc_filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv,NULL))) {
            arc_log(ARC_WARNING,"filter error (info=%d)\n",info);
        }
        free(R);
    }
    free(v); free(H);
}
static int arc_cmpel(const void *p1,const void *p2)
{
    return *((double*)p1)>=*((double*)p2);
}
/* adjust satellites grouping-----------------------------------------------*/
static double arc_amb_adjust_el(const rtk_t *rtk,const int* sat,int ns)
{
    int i;
    double el=rtk->opt.amb_el_group,*satel=arc_mat(ns,1);

    for (i=0;i<ns;i++) satel[i]=rtk->ssat[sat[i]-1].azel[1];

    /* sort satellites elevations */
    qsort(satel,ns,sizeof(double),arc_cmpel);

    /* adjust group-el'value */
    if (satel[0]>=el||satel[ns-1]<=el) el=satel[ns/2];

    arc_log(ARC_INFO,"arc_amb_adjust_el : el=%.3f\n",el);
    return el;  /* ensure two-el-group have elements */
}
/* resolve part integer ambiguity by LAMBDA-----------------------------------
 * here just reduce the max variance of ambiguity---------------------------*/
static int arc_resamb_reduceQ(int _na_,int _nb_,int na,int nb,double *Qy,
                              double *Qb,double *Qab,int *ix)
{
    int i,j,k=0,*ixx=arc_imat(na+nb,1);
    double val=-9999.0;

    arc_log(ARC_INFO,"arc_resamb_reduceQ : \n");

    /* find max variance position */
    for (i=0;i<_nb_;i++) {
        if ((Qy[_na_+i+(_na_+i)*(_na_+_nb_)]>=val)&&ix[_na_+i]) {
            val=Qy[_na_+i+(_na_+i)*(_na_+_nb_)]; k=i;
        }
    }
    for (i=0,j=0;i<_nb_;i++) {
        if (k==i||!ix[_na_+i]) continue; ixx[na+(j++)]=i+_na_;
    }
    for (i=0;i<_na_;i++) ixx[i]=i;

    arc_log(ARC_INFO,"arc_resamb_reduceQ : index=\n");
    arc_tracemati(ARC_MATPRINTF,ixx,1,nb+na-1,10,4);

    /* fill Qb matrix */
    if (Qb) for (i=0;i<nb-1;i++) {
        for (j=0;j<nb-1;j++) Qb[i+(nb-1)*j]=Qy[ixx[i+_na_]+ixx[j+_na_]*(_na_+_nb_)];
    }
    /* fill Qab matrix */
    if (Qab) for (i=0;i<na;i++) {
        for (j=0;j<nb-1;j++) Qab[i+j*na]=Qy[ixx[i]+ixx[j+_na_]*(_na_+_nb_)];
    }
    ix[_na_+k]=0;
    free(ixx); return k;
}
/* lambda/mlambda integer least-square estimation-----------------------------*/
static int arc_amb_lambda(int nb,double *y,double *Qb,double *b,double *s,
                          double*ratio)
{
    int info=1; if (!(info=arc_lambda(nb,2,y,Qb,b,s,NULL,NULL))) {

        arc_log(ARC_INFO, "N(1)=");
        arc_tracemat(ARC_MATPRINTF,b,1,nb,10,3);
        arc_log(ARC_INFO, "N(2)=");
        arc_tracemat(ARC_MATPRINTF,b+nb,1,nb,10,3);

        *ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (*ratio>999.9) *ratio=999.9f;
    }
    else {
        info=1; arc_log(ARC_WARNING,"arc_amb_lambda : lambda failed");
    }
    return info;
}
/* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) -----*/
static void arc_amb_s2d(const rtk_t *rtk,const double *D,double *y,double *Qy,
                        int ny,int nx)
{
    double *DP=arc_mat(ny,nx);
    if (y) arc_matmul("TN",ny,1,nx,1.0,D,rtk->x,0.0,y);
    if (Qy) {
        arc_matmul("TN",ny,nx,nx,1.0,D,rtk->P,0.0,DP);
        arc_matmul("NN",ny,ny,nx,1.0,DP,D,0.0,Qy);
    }
    free(DP);
}
/* extract ambiguity covariance matrix----------------------------------------*/
static void arc_amb_Qb(int nb,int na,int ny,double* Qb,double*Qab,
                       const double* Qy)
{
    int i,j;
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]=Qy[na+i+(na+j)*ny];
    for (i=0;i<na;i++) for (j=0;j<nb;j++) Qab[i+j*na]=Qy[   i+(na+j)*ny];
}
/* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b))--------------------*/
static int arc_amb_fix_sol(rtk_t *rtk,int nb,int na,double *Qb,const double *y,
                           const double *Qab)
{
    int info=1;
    double *db=arc_mat(nb,1),*QQ=arc_mat(na,nb);

    if (!(info=arc_matinv(Qb,nb))) {
        arc_matmul("NN",nb,1,nb,1.0,Qb,y+na,0.0,db);
        arc_matmul("NN",na,1,nb,-1.0,Qab,db,1.0,rtk->xa);

        /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
        arc_matmul("NN",na,nb,nb,1.0,Qab,Qb,0.0,QQ);
        arc_matmul("NT",na,na,nb,-1.0,QQ,Qab,1.0,rtk->Pa);

        arc_log(ARC_INFO,"arc_amb_fix_sol : xa=\n");
        arc_tracemat(ARC_MATPRINTF,rtk->xa,na,1,13,4);
        arc_log(ARC_INFO,"arc_amb_fix_sol : Pa=\n");
        arc_tracemat(ARC_MATPRINTF,rtk->Pa,na,na,13,4);
    }
    free(db); free(QQ);
    return info;
}
/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
static int arc_resamb_LAMBDA(rtk_t *rtk,double *bias,double *xa)
{
    prcopt_t *opt=&rtk->opt;
    int i,j,ny,nb,nx=rtk->nx,na=rtk->na,k,n,*ix,_na_,_nb_,ok=0;
    double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*QQ,s[2],*yb,*DB=NULL;

    arc_log(ARC_INFO, "arc_resamb_LAMBDA : nx=%d\n", nx);

    rtk->sol.ratio=0.0; rtk->sol.dop.dops[4]=-999.0; rtk->sol.p_ar=(float)-999.0;

    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
        rtk->opt.thresar[0]<1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D=arc_zeros(nx,nx);
    if ((nb=arc_ddmat(rtk,D))<=0) {
        arc_log(ARC_WARNING, "arc_resamb_LAMBDA : no valid double-difference\n");
        free(D);
        return 0;
    }
    ny=na+nb; _na_=na; _nb_=nb;
    y=arc_mat(ny,1); Qy=arc_mat(ny,ny);DP=arc_mat(ny,nx);
    b=arc_mat(nb,2); db=arc_mat(nb,1); Qb=arc_mat(nb,nb);
    Qab=arc_mat(na,nb); QQ=arc_mat(na,nb);
    yb=arc_mat(nb+na,1);  /* save the reduce y value */
    ix=arc_imat(na+nb,1);  /* hold active ambiguity index */
    DB=arc_mat(1,nb);

    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    arc_matmul("TN",ny,1,nx,1.0, D,rtk->x,0.0,y);
    arc_matmul("TN",ny,nx,nx,1.0,D,rtk->P,0.0,DP);
    arc_matmul("NN",ny,ny,nx,1.0,DP,D,0.0,Qy);

    arc_log(ARC_INFO,"arc_resamb_LAMBDA : Qy=\n");
    arc_tracemat(ARC_MATPRINTF,Qy,na+nb,na+nb,10,4);

    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]=Qy[na+i+(na+j)*ny];
    for (i=0;i<na;i++) for (j=0;j<nb;j++) Qab[i+j*na]=Qy[   i+(na+j)*ny];
    for (i=0;i<na+nb;i++) ix[i]=1;

    arc_log(ARC_INFO, "arc_resamb_LAMBDA : N(0)=");
    arc_tracemat(ARC_MATPRINTF,y+na,1,nb,10,3);

    arc_log(ARC_INFO,"arc_resamb_LAMBDA : Qb= \n");
    arc_tracemat(ARC_MATPRINTF,Qb,nb,nb,10,3);

    arc_log(ARC_INFO,"arc_resamb_LAMBDA : Qab=\n");
    arc_tracemat(ARC_MATPRINTF,Qab,na,nb,10,4);

    /* updates adop */
    if (rtk->opt.amb_adop) {
        rtk->sol.dop.dops[4]=arc_amb_adop(Qb,y+na,nb);
        arc_log(ARC_INFO,"ADOP=%8.4lf \n",rtk->sol.dop.dops[4]); /* todo:need to test */
    }
    arc_matcpy(yb,y,nb+na,1);

    /* lambda/mlambda integer least-square estimation */
    for (n=0;n<(rtk->opt.amb_part?rtk->opt.amb_iter:1);n++) {

        if (nb==0) break;  /* no ambiguitys */

        if (!(arc_lambda(nb,2,yb+na,Qb,b,s,DB,NULL))) {

            arc_log(ARC_INFO, "N(1)=");
            arc_tracemat(ARC_MATPRINTF,b,1,nb,10,3);
            arc_log(ARC_INFO, "N(2)=");
            arc_tracemat(ARC_MATPRINTF,b+nb,1,nb,10,3);

            /* updates ar sucess probability value */
            rtk->sol.p_ar=arc_amb_bs_success(DB,nb);

            arc_log(ARC_INFO,"AR sucess probability=%8.4f \n",rtk->sol.p_ar);

            rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
            if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;

            /* validation by popular ratio-test */
            ok=(opt->amb_ffratio?arc_FFRatio(nb,opt->amb_ffailure,
                                             rtk->sol.ratio,1.0-rtk->sol.p_ar):
                                 (s[0]<=0.0||s[1]/s[0]>=(opt->thresar[0]))); /* todo:need to test */
            if (ok) {

                /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
                for (i=0;i<na;i++) {
                    rtk->xa[i]=rtk->x[i];
                    for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
                }
                for (i=0;i<nb;i++) {
                    bias[i]=b[i];   /* fixed solutions */
                    yb[na+i]-=b[i]; /* b0-b */
                }
                if (!arc_matinv(Qb,nb)) {
                    arc_matmul("NN",nb,1,nb,1.0,Qb,yb+na,0.0,db);
                    arc_matmul("NN",na,1,nb,-1.0,Qab,db,1.0,rtk->xa);

                    /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                    arc_matmul("NN",na,nb,nb,1.0,Qab,Qb,0.0,QQ);
                    arc_matmul("NT",na,na,nb,-1.0,QQ,Qab,1.0,rtk->Pa);

                    arc_log(ARC_INFO,"arc_resamb_LAMBDA : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                            nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);

                    /* restore single-differenced ambiguity */
                    arc_restamb(rtk,bias,xa,0); /* bias is the fixed solutions */
                    break; /* lambda solve success */
                }
                else nb=0;
            }
            else {
                if (rtk->opt.amb_part) { /* todo:this option is not good for ARMODE_FIXHOLD */
                    
                    if (n==((rtk->opt.amb_part?rtk->opt.amb_iter:1)-1)) {
                        nb=0; break;
                    }
                    /* find max ambiguity variance position */
                    k=arc_resamb_reduceQ(_na_,_nb_,na,nb,Qy,Qb,Qab,ix);

                    arc_log(ARC_INFO,"arc_resamb_LAMBDA : reduce Qb=\n");
                    arc_tracemat(ARC_MATPRINTF,Qb,nb-1,nb-1,10,4);
                    arc_log(ARC_INFO,"arc_resamb_LAMBDA : reduce Qab=\n");
                    arc_tracemat(ARC_MATPRINTF,Qab,na,nb-1,10,4);

                    arc_tracemati(ARC_MATPRINTF,ix,1,ny,2,1);

                    /* reset the value of ambiguity list */
                    for (i=0,j=0;i<_nb_;i++) {
                        if (!ix[i+na]) continue; yb[na+(j++)]=y[na+i];
                    }
                    arc_log(ARC_INFO,"arc_resamb_LAMBDA : reduce y=\n");
                    arc_tracemat(ARC_MATPRINTF,yb+na,1,j,10,4);

                    /* reset the numbers of ambiguitys */
                    nb=j;
                    /* reset fix flag */
                    rtk->ssat[rtk->amb_index[k]-1].fix[0]=0; /* no fix this double-differnce ambiguity,
                                                              * don't warn about reference ambiguity */
                    /* todo:may be the reference single-difference ambiguity,need to fix this bug */
                }
                else nb=0;
            }
        }
    }
    if (rtk->sol.ratio<opt->thresar[0]) {
        arc_log(ARC_WARNING, "arc_resamb_LAMBDA : ambiguity validation "
                        "failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",nb,s[1]/s[0],s[0],s[1]);
    }
    free(D); free(y); free(Qy); free(DP); free(ix); free(DB);
    free(b); free(db); free(Qb); free(Qab); free(QQ); free(yb);

    return nb; /* number of ambiguities */
}
/* resolve integer ambiguity by group-LAMBDA ---------------------------------*/
static int arc_resamb_group_LAMBDA(rtk_t *rtk,double *bias,double *xa)
{
    return arc_resamb_LAMBDA(rtk,bias,xa);
}
/* validation of solution ----------------------------------------------------*/
static int arc_valpos(rtk_t *rtk,const double *v,const double *R,const int *vflg,
                      int nv,double thres)
{
    double fact=thres*thres;
    int i,stat=1,sat1,sat2,type,freq;
    char stype[8];

    arc_log(ARC_INFO,"arc_valpos  : nv=%d thres=%.1f\n",nv,thres);

    /* post-fit residual test */
    for (i=0;i<nv;i++) {
        if (v[i]*v[i]<=fact*R[i+i*nv]) continue;
        sat1=(vflg[i]>>16)&0xFF;
        sat2=(vflg[i]>> 8)&0xFF;
        type=(vflg[i]>> 4)&0xF;
        freq=vflg[i]&0xF;
        strcpy(stype,type==0?"L":(type==1?"L":"C"));
        arc_log(ARC_WARNING,"arc_valpos : "
                        "large residual (sat=%2d-%2d %s%d v=%6.3f sig=%.3f)\n",
                sat1,sat2,stype,freq+1,v[i],SQRT(R[i+i*nv]));
    }
    return stat;
}
/* Q parameters for adaptive Kaman filter-------------------------------------*/
static int arc_adap_Q(const rtk_t *rtk,double *Q,int n)
{
    arc_log(ARC_INFO,"arc_adap_Q : \n");

    int i,k;
    /* rover station position noise */
    for (i=0;i<3;i++) Q[i+n*i]=SQR(rtk->opt.prn[5]);
    /* ambguity noise */
    for (i=3,k=0;i<rtk->nx;i++) {
        if (rtk->ceres_active_x[i]) Q[(3+k)*n+(3+(k++))]=SQR(rtk->opt.prn[0]);
    }
    /* Q is positive definite matrix,and its col is equal to n */
    if ((k+3)!=n) return 0;
    return 1;
}
/* C0 matrix for adaptive Kaman filter----------------------------------------*/
static int arc_adap_C0(const rtk_t* rtk,const double *v,double *C0,int m,
                       double lam)
{
    arc_log(ARC_INFO,"arc_adap_C0 : \n");

    static int first=1;
    static double lamk;

    if (first) {
        arc_matmul("NT",m,m,1,0.5,v,v,0.0,C0); return 1; first=0;
    }
    if (lam<1.0) return 0;
    lamk=lam/(1.0+lam);
    arc_matmul("NT",m,m,1,lamk,v,v,0.0,C0);
    return 1;
}
/* M matrix for adaptive Kaman filter-----------------------------------------*/
static int arc_adap_M(const rtk_t* rtk,const double *H,const double *P,
                      const double *R,int m,int n,double *M)
{
    arc_log(ARC_INFO,"arc_adap_M : \n");

    double *F;
    F=arc_mat(n,m);
    arc_matmul("NN",n,m,n,1.0,P,H,0.0,F);
    arc_matmul("TN",m,m,n,1.0,H,F,0.0,M);

    free(F); return 1;
}
/* N matrix for adaptive Kaman filter-----------------------------------------*/
static int arc_adap_N(const rtk_t* rtk,const double *H,const double *Q,
                      const double *R,const double *C0,int m,int n,double *N)
{
    arc_log(ARC_INFO,"arc_adap_N : \n");

    int i,j;
    double *F;

    F=arc_mat(n,m);
    arc_matcpy(N,R,m,m);
    arc_matmul("NN",n,m,n,1.0,Q,H,0.0,F);
    arc_matmul("TN",m,m,n,1.0,H,F,1.0,N);
    for (i=0;i<m;i++) for (j=0;j<m;j++) N[i+j*m]=C0[i+j*m]-N[i+j*m];
    return 1;
}
/* adaptive Kaman filter------------------------------------------------------*/
extern int adap_kaman_filter(rtk_t* rtk,double *x,double *P,const double *H,
                             const double *v,const double *R,int n,int m)
{
    arc_log(ARC_INFO,"adap_kaman_filter : \n");

    int i,j,*ix,k;
    double *C0,*M,*N,*Q,*H_,*P_;

    ix=arc_imat(n,1);
    for (i=0,k=0;i<rtk->nx;i++) if (rtk->ceres_active_x[i]) ix[k++]=i;
    Q=arc_zeros(k,k);
    if (!arc_adap_Q(rtk,Q,k)) {
        free(ix); free(Q); return 0;
    }
    C0=arc_mat(m,m);
    if (!arc_adap_C0(rtk,v,C0,m,rtk->lam)) {
        free(ix); free(Q); free(C0); return 0;
    }
    H_=arc_mat(k,m); M=arc_mat(m,m); P_=arc_mat(k,k);
    for (i=0;i<k;i++) {
        for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
        for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
    }
    if (!arc_adap_M(rtk,H_,P_,R,m,k,M)) {
        free(ix); free(Q); free(P_); free(C0); free(H_); free(M);
        return 0;
    }
    N=arc_mat(m,m);
    if (!arc_adap_N(rtk,H_,Q,R,C0,m,k,N)) {
        free(ix); free(Q); free(P_);
        free(C0); free(H_); free(M); free(N);
        return 0;
    }
    rtk->lam=MAX(1.0,arc_mattrace(N,m)/arc_mattrace(M,m));

    double *F=arc_mat(k,m),*_Q_=arc_mat(m,m),*K=arc_mat(k,m),
            *I=arc_eye(k),*xp=arc_mat(k,1),*Pp=arc_zeros(k,k);
    arc_matcpy(_Q_,R,m,m);
    arc_matcpy(Pp,P_,k,k);
    for (i=0;i<k;i++) xp[i]=x[ix[i]];

    arc_matmul("NN",k,m,k,1.0,P_,H_,0.0,F);            /* Q=H'*P*H+R */
    arc_matmul("TN",m,m,k,1.0,H_,F,1.0,_Q_);
    if (!(arc_matinv(_Q_,m))) {
        arc_matmul("NN",k,m,m,1.0,F,_Q_,0.0,K);        /* K=P*H*Q^-1 */
        arc_matmul("NN",k,1,m,rtk->lam,K,v,1.0,xp);    /* xp=x+K*v */
        arc_matmul("NT",k,k,m,-1.0,K,H_,1.0,I);        /* Pp=(I-K*H')*P */
        arc_matmul("NN",k,k,k,1.0,I,P_,0.0,Pp);
    }
    for (i=0;i<k;i++) {
        if (x) x[ix[i]]=xp[i];
        if (P) for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp[i+j*k];
    }
    free(F);  free(_Q_); free(K); free(I); free(xp); free(Pp);
    free(ix); free(Q); free(P_);
    free(C0); free(H_); free(M); free(N);
    return 1;
}
/* ukf: Propagate the sigma points through the dynamic model-----------------*/
static void arc_ukf_filterfunc(int states_dim,double *input_X,double *output_X)
{
    arc_log(ARC_INFO,"arc_ukf_filterfunc : \n");

    int i;
    for (i=0;i<states_dim;i++) output_X[i]=input_X[i];  /* todo:here have not think about dynamic model */
}
/* ukf: Propagate sigma points through the measurement model -----------------*/
static void arc_ukf_measfunc(double *input_X,double *output_Y)  /* todo:this function have some problem and bugs,need to fix */
{
    arc_log(ARC_INFO,"arc_ukf_measfunc : \n");

    int i,ndd=0;
    double dt,*xp=arc_zeros(_UKF_RTK_->nx,1),*ddy=arc_mat(_UKF_NV_,1);

    /* adjust input states dim,make it meet rtk_t struct */
    /* todo:this method is not best,and need to be modify in future */
    for (i=0;i<_UKF_ANX_;i++) xp[_UKF_IX_[i]]=input_X[i];

    /* base and rover station observation time difference */
    dt=timediff(_UKF_OBS_[0].time,_UKF_OBS_[_UKF_NU_].time);

    /* undifferenced residuals for rover */
    if (!(arc_zdres(0,_UKF_OBS_,_UKF_NU_,_UKF_RS_,_UKF_DTS_,_UKF_SVH_,
                    _UKF_NAV_,xp,&_UKF_RTK_->opt,0,NULL,_UKF_E_,_UKF_AZEL_,
                    _UKF_RTK_,_UKF_Y_))) {
        arc_log(ARC_WARNING, "arc_ukf_measfunc : rover initial position error\n");
        free(xp); free(ddy); return;
    }
    /* double-differenced measurements for ukf*/
    if ((ndd=arc_ukf_ddmeas(_UKF_RTK_,_UKF_NAV_,dt,xp,NULL,_UKF_SAT_,_UKF_IU_,_UKF_IR_,
                            _UKF_NS_,_UKF_AZEL_,_UKF_Y_,ddy))<1) {
        arc_log(ARC_WARNING, "arc_ukf_measfunc : no double-differenced residual\n");
        free(xp); free(ddy); return;
    }
    ARC_ASSERT_TRUE_DBG(Exception,ndd==_UKF_NV_,"Updates measurements Failed");

    if (output_Y) arc_matcpy(output_Y,ddy,_UKF_NV_,1);

    free(xp); free(ddy); return;
}
/* ukf: set the all needed paranmeter pointor,just for initial ukf process----*/
static void arc_ukf_initial(rtk_t *rtk,const nav_t* nav,const obsd_t* obs,int nu,int nr,
                            double *rs,double *dts,double *e,int *svh,int *vflag,double *azel,
                            double *kf_y)
{
    arc_log(ARC_INFO,"arc_ukf_initial :\n");

    /* before using ukf,must use this function to initial all data and struct */
    _UKF_RTK_  =rtk;       /* rtk data struct */
    _UKF_NAV_  =nav;       /* navigation data */
    _UKF_OBS_  =obs;       /* observation data */
    _UKF_NU_   =nu;        /* numbers of rover observations */
    _UKF_NR_   =nr;        /* numbers of base observations */
    _UKF_NX_   =rtk->nx;   /* numbers of states,todo:this need to re-think,because this states numbers is large */
    _UKF_RS_   =rs;        /* satellite position/velecity */
    _UKF_SVH_  =svh;       /* satellite health flag */
    _UKF_AZEL_ =azel;      /* satellite-rover station azel pointor */
    _UKF_E_    =e;         /* observation satellite - rover line of light */
    _UKF_DTS_  =dts;       /* satellite clock/dirft pointor */
    if (_UKF_R_==NULL) {   /* initial double-difference measurements variance */
        _UKF_R_=arc_zeros(MAXSAT,MAXSAT);   /* todo:may be allocate too large memory */
    }
    if (_UKF_Y_==NULL) {   /* initial sigma point single-difference measurements vector */
        _UKF_Y_=arc_zeros(MAXSAT*2*2,1);    /* todo:may be allocate too large memory */
    }
    if (_UKF_DY_==NULL) {  /* initial sigma point double-difference measurements vector */
        _UKF_DY_=arc_zeros(MAXSAT*2,1);     /* todo:may be allocate too large memory,can be modified */
    }
    if (_UKF_IX_==NULL) {  /* initial active states index list */
        _UKF_IX_=arc_imat(rtk->nx,1);
    }
    if (_UKF_Q_==NULL) {   /* initial states system noise matrix */
        _UKF_Q_=arc_zeros(rtk->nx,rtk->nx);
    }
    if (_UKF_XP_==NULL) {  /* initial prior states */
        _UKF_XP_=arc_zeros(rtk->nx,1);
    }
    if (_UKF_PP_==NULL) {  /* initial prior states covariance matrix */
        _UKF_PP_=arc_zeros(rtk->nx,rtk->nx);
    }
    if (_UKF_MEAS_==NULL) { /* initial measurements vector */
        _UKF_MEAS_=arc_zeros(MAXSAT,1);;
    }
}
/* ukf: set satellite comom list----------------------------------------------*/
static void arc_ukf_set_sat(int *sat,int *ir,int *iu,int ns)
{
    arc_log(ARC_INFO,"arc_ukf_set_sat : \n");

    _UKF_SAT_=sat;
    _UKF_IU_ =iu;
    _UKF_IR_ =ir;
    _UKF_NS_ =ns;
}
/* ukf: free ukf position variance--------------------------------------------*/
extern void arc_ukf_free_problem()  /* todo:may be have more better way to free ukf problem */
{
    arc_log(ARC_INFO,"arc_ukf_free_problem : \n");

    if (_UKF_MEAS_) free(_UKF_MEAS_);
    if (_UKF_PP_)   free(_UKF_PP_);
    if (_UKF_XP_)   free(_UKF_XP_);
    if (_UKF_Q_)    free(_UKF_Q_);
    if (_UKF_IX_)   free(_UKF_IX_);
    if (_UKF_DY_)   free(_UKF_DY_);
    if (_UKF_Y_)    free(_UKF_Y_);
    if (_UKF_R_)    free(_UKF_R_);
}
/* ukf: count active states index list----------------------------------------*/
static void arc_ukf_activex(const rtk_t *rtk)
{
    arc_log(ARC_INFO,"arc_ukf_activex : \n");

    int i,k=0; for (i=0;i<rtk->nx;i++) if (rtk->ceres_active_x[i]) _UKF_IX_[k++]=i;
    _UKF_ANX_=k;
}
/* get the active states index given raw index in rtk->x----------------------*/
static int arc_get_activex_index(const rtk_t *rtk,int index)
{
    arc_log(ARC_INFO,"arc_get_activex_index : \n");

    int i,k;
    if (index<0) return -1;
    for (i=0,k=0;i<rtk->nx;i++) {
        if (i==index&&rtk->ceres_active_x[i]) return k;
        if (rtk->ceres_active_x[i]) k++;  /* todo:may be have more better way to find index */
    }
    return -1;
}
/* ukf: states system noise matrix--------------------------------------------*/
static void arc_ukf_Q(const rtk_t *rtk)
{
    arc_log(ARC_INFO,"arc_ukf_Q : \n");

    int i,j,k,indx;

    /* reset states systen noise matrix,this step is very important */
    for (i=0;i<_UKF_ANX_;i++) {
        for (j=0;j<_UKF_ANX_;j++) _UKF_Q_[i*_UKF_ANX_+j]=0.0;
    }
    /* build the states system noise matrix */
    if (_UKF_Q_) {
        for (i=0;i<3;i++) _UKF_Q_[i+_UKF_ANX_*i]=SQR(rtk->opt.prn[5]);
        if (rtk->opt.tropopt==TROPOPT_EST
            ||rtk->opt.tropopt==TROPOPT_ESTG) {  /* troposphere model states */
            for (k=0;k<(rtk->opt.tropopt<TROPOPT_ESTG?1:3);k++) {
                /* updates Q matrix */
                if ((indx=arc_get_activex_index(rtk,IT(0,&rtk->opt)+k))<0)
                    continue;
                _UKF_Q_[indx+indx*_UKF_ANX_]=SQR(rtk->opt.prn[2]);
                if ((indx=arc_get_activex_index(rtk,IT(1,&rtk->opt)+k))<0)
                    continue;
                _UKF_Q_[indx+indx*_UKF_ANX_]=SQR(rtk->opt.prn[2]);
            }
        }
        if (rtk->opt.ionoopt==IONOOPT_EST) {
            /* updates inon states system noise */
            for (i=0;i<MAXSAT;i++) {
                if (rtk->ceres_active_x[II(i+1,&rtk->opt)]) {
                    if ((indx=(arc_get_activex_index(rtk,II(i+1,&rtk->opt))))<0)
                        continue;
                    /* iono states process noise */
                    _UKF_Q_[indx+indx*_UKF_ANX_]=SQR(rtk->opt.prn[1]);
                }
            }
        }
        /* ambguity process noise */
        for (i=0;i<MAXSAT;i++) {
            if (!((satsys(i+1,NULL)&_UKF_RTK_->opt.navsys)&&
                    _UKF_RTK_->opt.exsats[i]!=1))
                continue;
            if (rtk->ceres_active_x[IB(i+1,0,&rtk->opt)]==0)
                continue;
            if ((indx=arc_get_activex_index(rtk,IB(i+1,0,&rtk->opt)))<0)
                continue;
            _UKF_Q_[indx+indx*_UKF_ANX_]=SQR(rtk->opt.prn[0]);
        }
    }
}
/* ukf:set the states prior value and states covariance matrix----------------*/
static void arc_ukf_get_prior_XP(const rtk_t *rtk,double *xp,double *Pp)
{
    arc_log(ARC_INFO,"arc_ukf_get_prior_XP : \n");

    int i,j;
    for (i=0;i<_UKF_ANX_;i++) {
        if (xp) xp[i]=rtk->x[_UKF_IX_[i]];
        if (Pp) {
            for (j=0;j<_UKF_ANX_;j++)
                Pp[i*_UKF_ANX_+j]=rtk->P[rtk->nx*_UKF_IX_[i]+_UKF_IX_[j]];
        }
    }
}
/* ukf:get ukf updates states-------------------------------------------------*/
static void arc_ukf_get_updatax(ukf_t *ukf,double *xp,double *Pp)
{
    arc_log(ARC_INFO,"arc_ukf_get_updatax : \n");

    int i,j;
    double *xpp=arc_zeros(ukf->state_dim,1),
           *Ppp=arc_zeros(ukf->state_dim,ukf->state_dim);
    arc_ukf_filter_get_state(ukf,xpp,Ppp);  /* get the ukf update states and covariance matrix */
    for (i=0;i<_UKF_ANX_;i++) {
        if (xp) xp[_UKF_IX_[i]]=xpp[i];     /* todo:here may be modified */
        if (Pp) for (j=0;j<_UKF_ANX_;j++)
                Pp[_UKF_IX_[i]*_UKF_RTK_->nx+_UKF_IX_[j]]=Ppp[i*_UKF_ANX_+j];
    }
    free(xpp); free(Ppp);
}
/* make double-differnce measurements-----------------------------------------*/
static int arc_ukf_dd_meas()
{
    arc_log(ARC_INFO,"arc_ukf_dd_meas : \n"); 

    int i,k,iu[2]={0},ir[2]={0};

    if (_UKF_MEAS_==NULL) return 0;  /* ukf struct must be initial */

    for (i=0,k=0;i<_UKF_NV_;i++) {  /* todo:loop times is _UKF_NV_,here may be have some problems */
        if (arc_ukf_get_dd_obsind(_UKF_RTK_,ir,iu,_UKF_OBS_,
                                  _UKF_RTK_->sat[2*i],_UKF_RTK_->sat[2*i+1]))
            continue;
        /* carry phase double-difference observation */
        if (i<_UKF_RTK_->nc)
            _UKF_MEAS_[k++]=-(_UKF_OBS_[iu[0]].L[0]*_UKF_NAV_->lam[_UKF_OBS_[iu[0]].sat-1][0]
                             -_UKF_OBS_[iu[1]].L[0]*_UKF_NAV_->lam[_UKF_OBS_[iu[1]].sat-1][0]-
                             (_UKF_OBS_[ir[0]].L[0]*_UKF_NAV_->lam[_UKF_OBS_[ir[0]].sat-1][0]
                             -_UKF_OBS_[ir[1]].L[0]*_UKF_NAV_->lam[_UKF_OBS_[ir[1]].sat-1][0]));
        /* pseudorange double-difference observation */
        else if (i>=_UKF_RTK_->nc)
            _UKF_MEAS_[k++]=-(_UKF_OBS_[iu[0]].P[0]-_UKF_OBS_[iu[1]].P[0]-
                             (_UKF_OBS_[ir[0]].P[0]-_UKF_OBS_[ir[1]].P[0]));
    }
    ARC_ASSERT_TRUE_DBG(Exception,k==_UKF_NV_,"make double-differnce measurements failed");

    if (k!=_UKF_NV_) return 0;  /* todo:this condition may be not best */
    return k;
}
/* states updates of difference pseudorange positioning-----------------------*/
static void arc_diff_pr_update(const rtk_t* rtk,double *xp,double *Pp,double tt)
{
    int i,j,nx=NP(&rtk->opt);
    double var;

    arc_log(ARC_INFO,"arc_diff_pr_update   : tt=%.3f\n",tt);

    /* fixed mode */
    if (rtk->opt.mode==PMODE_FIXED) {
        for (i=0;i<3;i++) {
            arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],1E-6,i,nx);
        }
        return;
    }
    /* initialize position for first epoch */
    if (arc_norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) {
            arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_POS,i,nx);
        }
        /* initial rover station velecity */
        if (rtk->opt.dynamics) for (i=3;i<nx;i++) {
                arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_VEL,i,nx);
            }
    }
    /* static mode */
    if (rtk->opt.mode==PMODE_STATIC) return;

    /* todo:using standard positioning to initial ukf prior states and its covariacne matrix,
     * todo:and may have more better methids to do this */

    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
        if (rtk->opt.init_pnt) {
            for (i=0;i<3;i++) Pp[i+i*nx]=rtk->sol.prsol.qr[i];
            Pp[1+nx*0]=Pp[0+nx*1]=rtk->sol.prsol.qr[3];
            Pp[2+nx*0]=Pp[0+nx*2]=rtk->sol.prsol.qr[5];
            Pp[1+nx*2]=Pp[2+nx*1]=rtk->sol.prsol.qr[4]; /* todo:may be have more better way to do this */
            return;
        }
        else {
            if (rtk->sol.ratio>=rtk->opt.thresar[0]) {
                for (i=0;i<3;i++) arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_POS_AMB,i,nx);
                return;
            }
            for (i=0;i<3;i++) arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_POS,i,nx);
            return;
        }
    }
    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=rtk->P[i+i*rtk->nx]; var/=3.0;

    if (var>VAR_POS||var<0.0) {
        /* reset position with large variance */
        for (i=0;i<3;i++) arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_POS,i,nx);
        for (i=3;i<6;i++) arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_VEL,i,nx);
        arc_log(ARC_INFO,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
    /* state transition of position/velocity/acceleration for dynamic mode */
    double *F=arc_eye(nx),*FP=arc_mat(nx,nx),
            *xpp=arc_mat(rtk->nx,1),pos[3]={0},Q[9]={0},Qv[9]={0};

    /* compute the state transition matrix */
    for (i=0;i<3;i++) F[i+(i+3)*nx]=tt;

    arc_log(ARC_INFO,"arc_diff_pr_update : state transition matrix=\n");
    arc_tracemat(ARC_MATPRINTF,F,nx,nx,10,4);

    arc_log(ARC_INFO,"arc_diff_pr_update : before state transition P=\n");
    arc_tracemat(ARC_MATPRINTF,Pp,nx,nx,10,4);

    /* x=F*x, P=F*P*F+Q */
    arc_matmul("NN",nx,1,nx,1.0,F,xp,0.0,xpp);
    arc_matcpy(xp,xpp,nx,1);
    arc_matmul("NN",nx,nx,nx,1.0,F,Pp,0.0,FP);
    arc_matmul("NT",nx,nx,nx,1.0,FP,F,0.0,Pp);

    arc_log(ARC_INFO,"arc_diff_pr_update : after state transition P= \n");
    arc_tracemat(ARC_MATPRINTF,Pp,nx,nx,10,4);

    /* process noise added to only velecity */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3]); Q[8]=SQR(rtk->opt.prn[4]);
    ecef2pos(xp,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
            Pp[i+3+(j+3)*nx]+=Qv[i+j*3];  /* add system process noise */
        }
    free(F); free(FP); free(xpp);
}
/* difference pseudorange positioning-----------------------------------------*/
static int arc_diff_pr_relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                              const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,*v,*H,*R,*xp,*Pp,*bias,dt,*D=NULL,*Qv=NULL;
    int i,j,n=nu+nr,ns,ny,nv,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],niter,nx=NP(&rtk->opt);
    int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
    int nf=1;

    arc_log(ARC_INFO,"arc_diff_pr_relpos  : nx=%d nu=%d nr=%d\n",rtk->nx,nu,nr);

    dt=timediff(time,obs[nu].time);

    rs=arc_mat(6,n); dts=arc_mat(2,n);
    var=arc_mat(1,n); y=arc_mat(nf*2,n); e=arc_mat(3,n);
    azel=arc_zeros(2,n);

    /* initial satellite valid flag and snr */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        rtk->ssat[i].vsat[0]=0;
        rtk->ssat[i].snr [0]=0;
    }
    /* satellite positions/clocks */
    arc_satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        arc_testeclipse(obs,n,nav,rs);
    }
    /* undifferenced residuals for base station */
    if (!arc_zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,svh+nu,nav,rtk->rb,opt,1,
                   y+nu*nf*2,e+nu*3,azel+nu*2,rtk,NULL)) {
        arc_log(ARC_WARNING,"arc_diff_pr_relpos : initial base station position error\n");
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    /* time-interpolation of residuals (for post-processing) */
    if (opt->intpref) {
        dt=arc_intpres(time,obs+nu,nr,nav,rtk,y+nu*nf*2);
    }
    /* select common satellites between rover and base-station */
    if ((ns=arc_selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
        arc_log(ARC_WARNING,"arc_diff_pr_relpos : no common satellite\n");

        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }

    /* initial states for kalman filter */
    xp=arc_mat(nx,1); Pp=arc_zeros(nx,nx);
    arc_matcpy(xp,rtk->x,nx,1);
    for (i=0;i<nx;i++) for (j=0;j<nx;j++) Pp[i+j*nx]=rtk->P[i+rtk->nx*j];

    /* updates states */
    arc_diff_pr_update(rtk,xp,Pp,fabs(rtk->tt));

    arc_log(ARC_INFO,"arc_diff_pr_relpos: xp=\n");
    arc_tracemat(ARC_MATPRINTF,xp,1,nx,13,4);

    arc_log(ARC_INFO,"arc_diff_pr_relpos : Pp=\n");
    arc_tracemat(ARC_MATPRINTF,Pp,nx,nx,10,4);

    ny=ns*nf*2+2;
    v=arc_mat(ny,1); H=arc_zeros(nx,ny);
    R=arc_mat(ny,ny); bias=arc_mat(nx,1);

    /* add 2 iterations for baseline-constraint moving-base */
    niter=opt->niter+(opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0?2:0);

    for (i=0;i<niter;i++) {  /* iterations compute */

        /* undifferenced residuals for rover */
        if (!arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk,NULL)) {
            arc_log(ARC_WARNING,"arc_diff_pr_relpos : rover initial position error\n");
            break;
        }
        /* double-differenced residuals and partial derivatives */
        if ((nv=arc_diff_pr_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg))<1) {
            arc_log(ARC_WARNING,"arc_diff_pr_relpos : no double-differenced residual\n");
            break;
        }
        arc_log(ARC_INFO,"arc_diff_pr_relpos ：double-differenced residual vector : \n");
        arc_tracemat(ARC_MATPRINTF,v,nv,1,10,4);

        /* robust kalman filter */
        if (opt->kalman_robust) {
            D=arc_eye(nv); Qv=arc_mat(nv,nv);
            arc_kalman_norm_inno(rtk,v,nv,H,R); arc_kalman_norm_Qino(v,nv,Qv);
            arc_kalman_robust_phi(rtk,v,nv,D); /* todo:need to test */
        }
        /* adaptive kaman filter */
        if (opt->adapt_filter) {
            if (!adap_kaman_filter(rtk,xp,Pp,H,v,R,nx,nv)) {
                arc_log(ARC_WARNING,"arc_diff_pr_relpos : adaptive filter error (info=%d)\n",info);
                break;
            }
        }
        /* kalman filter measurement update */
        else {
            if ((info=arc_filter(xp,Pp,H,v,R,nx,nv,opt->kalman_robust?D:NULL))) {
                arc_log(ARC_WARNING,"arc_diff_pr_relpos : filter error (info=%d)\n",info);
                break;
            }
            arc_log(ARC_INFO,"arc_diff_pr_relpos : x(%d)=",i+1);
            if (opt->mode==PMODE_STATIC)                arc_tracemat(ARC_MATPRINTF,xp,3,1,10,4);
            if (opt->mode==PMODE_KINEMA&&opt->dynamics) arc_tracemat(ARC_MATPRINTF,xp,6,1,10,4);

            arc_log(ARC_INFO,"arc_diff_pr_relpos : P(%d)=",i+1);
            arc_tracemat(ARC_MATPRINTF,Pp,nx,nx,10,4);
        }
    }
    /* save solution status */
    arc_matcpy(rtk->sol.prsol.rr,xp,nx,1); /* position/velecity */

    for (j=0;j<nx;j++) rtk->sol.prsol.qr[j]=(float)Pp[j+j*nx];

    rtk->sol.prsol.qr[3]=(float)Pp[1];    /* cov xy */
    rtk->sol.prsol.qr[4]=(float)Pp[2+nx]; /* cov yz */
    rtk->sol.prsol.qr[5]=(float)Pp[2];    /* cov zx */

    /* todo:difference pseudorange positioning just for initialing for rtk-postioning */

    if (D) free(D); if (Qv) free(Qv);
    free(rs); free(dts); free(var); free(y); free(e); free(azel);
    free(xp); free(Pp);  free(v); free(H); free(R); free(bias);
    return 1;
}
/* relative positioning ------------------------------------------------------*/
#ifdef ARC_TEST
ofstream fp_ukf_ceres("/home/sujinglan/arc_rtk/arc_test/data/gps_bds/static/arc_ukf_pos");
static int I=0;
#endif
static int arc_relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                      const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,*v,*H,*R,*xp,*Pp,*xa,*bias,dt,*D=NULL,*Qv=NULL;
    int i,j,f,n=nu+nr,ns,ny,nv,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],niter;
    int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
    int stat=rtk->opt.mode<=PMODE_DGPS?SOLQ_DGPS:SOLQ_FLOAT;
    int nf=1;

#ifdef ARC_TEST
    I++;  /* just for debug */
#endif

    arc_log(ARC_INFO, "arc_relpos  : nx=%d nu=%d nr=%d\n", rtk->nx, nu, nr);

    dt=timediff(time,obs[nu].time);

    rs=arc_mat(6,n); dts=arc_mat(2,n);
    var=arc_mat(1,n); y=arc_mat(nf*2,n); e=arc_mat(3,n);
    azel=arc_zeros(2,n);

    /* initial satellite valid flag and snr */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        rtk->ssat[i].vsat[0]=0;
        rtk->ssat[i].snr [0]=0;
    }
    /* reset ceres problem solver active states index list */
    for (i=0;i<rtk->nx;i++) rtk->ceres_active_x[i]=0;

    /* satellite positions/clocks */
    arc_satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    /* initial ukf solution needed data struct */
    if (opt->ukf) arc_ukf_initial(rtk,nav,obs,nu,nr,rs,dts,e,svh,vflg,azel,y);

    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        arc_testeclipse(obs,n,nav,rs);
    }
    /* undifferenced residuals for base station */
    if (!arc_zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,svh+nu,nav,rtk->rb,opt,1,
                   y+nu*nf*2,e+nu*3,azel+nu*2,rtk,_UKF_Y_==NULL?NULL:_UKF_Y_+nu*nf*2)) {
        arc_log(ARC_WARNING, "arc_relpos : initial base station position error\n");
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    if (_UKF_Y_&&y) {
        arc_log(ARC_INFO,"base station undifferenced measurements for ukf :\n");
        arc_tracemat(ARC_MATPRINTF,_UKF_Y_+nu*nf*2,2,nr,10,4);
        arc_log(ARC_INFO,"base station undifferenced measurements for kf :\n");
        arc_tracemat(ARC_MATPRINTF,y+nu*nf*2,2,nr,10,4);
    }
    /* time-interpolation of residuals (for post-processing) */
    if (opt->intpref) {
        dt=arc_intpres(time,obs+nu,nr,nav,rtk,y+nu*nf*2);
    }
    /* select common satellites between rover and base-station */
    if ((ns=arc_selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
        arc_log(ARC_WARNING, "arc_relpos : no common satellite\n");

        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    /* set ukf common satellites between rover and base-station */
    if (opt->ukf) arc_ukf_set_sat(sat,ir,iu,ns);

    /* temporal update of states */
    arc_udstate(rtk,obs,sat,iu,ir,ns,nav,rs);

    xp=arc_mat(rtk->nx,1); Pp=arc_zeros(rtk->nx,rtk->nx);
    xa=arc_mat(rtk->nx,1);
    arc_matcpy(xp,rtk->x,rtk->nx,1);

    ny=ns*nf*2+2;
    v=arc_mat(ny,1); H=arc_zeros(rtk->nx,ny);
    R=arc_mat(ny,ny); bias=arc_mat(rtk->nx,1);

    /* add 2 iterations for baseline-constraint moving-base */
    niter=opt->niter+(opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0?2:0);

    if (opt->ceres==0) {
        if (opt->ukf==0) for (i=0;i<niter;i++) {  /* iterations compute */
            /* undifferenced residuals for rover */
            if (!arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk,NULL)) {
                arc_log(ARC_WARNING, "arc_relpos : rover initial position error\n");
                stat=SOLQ_NONE;
                break;
            }
            /* double-differenced residuals and partial derivatives */
            if ((nv=arc_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg))<1) {
                arc_log(ARC_WARNING, "arc_relpos : no double-differenced residual\n");
                stat=SOLQ_NONE;
                break;
            }
            arc_log(ARC_INFO,"arc_relpos ： double-differenced residual vector : \n");
            arc_tracemat(ARC_MATPRINTF,v,nv,1,10,4);

            /* robust kalman filter */
            if (opt->kalman_robust) {
                D=arc_eye(nv); Qv=arc_mat(nv,nv);
                arc_kalman_norm_inno(rtk,v,nv,H,R); arc_kalman_norm_Qino(v,nv,Qv);
                arc_kalman_robust_phi(rtk,v,nv,D); /* todo:need to test */
            }
            arc_matcpy(Pp,rtk->P,rtk->nx,rtk->nx);

            /* adaptive kaman filter */
            if (opt->adapt_filter) {
                if (!adap_kaman_filter(rtk,xp,Pp,H,v,R,rtk->nx,nv)) {
                    arc_log(ARC_WARNING,"arc_relpos : adaptive filter error (info=%d)\n",info);
                    stat=SOLQ_NONE;
                    break;
                }
            }
            /* kalman filter measurement update */
            else {
                if ((info=arc_filter(xp,Pp,H,v,R,rtk->nx,nv,opt->kalman_robust?D:NULL))) {
                    arc_log(ARC_WARNING, "arc_relpos : filter error (info=%d)\n",info);
                    stat=SOLQ_NONE;
                    break;
                }
                arc_log(ARC_INFO,"arc_relpos : x(%d)=",i+1);
                if (opt->mode==PMODE_STATIC)                arc_tracemat(ARC_MATPRINTF,xp,3,1,10,4);
                if (opt->mode==PMODE_KINEMA&&opt->dynamics) arc_tracemat(ARC_MATPRINTF,xp,6,1,10,4);

                arc_log(ARC_INFO,"arc_relpos : P(%d)=",i+1);
                arc_tracemat(ARC_MATPRINTF,Pp,rtk->nx,rtk->nx,10,4);
            }
        }
        else if (opt->ukf) {
            /* inital states covariance matrix */
            arc_matcpy(Pp,rtk->P,rtk->nx,rtk->nx);  /* todo:this step may reduce code running'speed */

            /* undifferenced residuals for rover,here is initial for ukf */
            if (!arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk,NULL)) {
                arc_log(ARC_WARNING, "arc_relpos : rover initial position error\n");
                stat=SOLQ_NONE;
            }
            /* double-differenced residuals and partial derivatives,here is initial for ukf */
            if ((_UKF_NV_=arc_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,
                                    ns,v,NULL,_UKF_R_,vflg))<1) {
                arc_log(ARC_WARNING, "arc_relpos : no double-differenced residual\n");
                stat=SOLQ_NONE;
            }

            /* start ukf */
            if (_UKF_NV_>=1) {  /* todo:this condition may not best */
                /* count active states index list */
                arc_ukf_activex(rtk);
                /* compute the states system noise matrix */
                arc_ukf_Q(rtk);

                arc_log(ARC_INFO,"arc_relpos : Q(ukf)=\n");
                arc_tracemat(ARC_MATPRINTF,_UKF_Q_,_UKF_ANX_,_UKF_ANX_,10,4);

                /* creates a new unscented kalman filter structure */
                ukf_t *ukf=arc_ukf_filter_new(_UKF_ANX_,_UKF_NV_,_UKF_Q_,_UKF_R_,
                                              arc_ukf_filterfunc,arc_ukf_measfunc);

                arc_log(ARC_INFO,"arc_relpos : R(ukf)=\n");
                arc_tracemat(ARC_MATPRINTF,_UKF_R_,_UKF_NV_,_UKF_NV_,10,4);

                /* set the prior states and its covariance matrix */
                arc_ukf_get_prior_XP(rtk,_UKF_XP_,_UKF_PP_);

                arc_log(ARC_INFO,"arc_relpos : x-prior(ukf)=\n");
                arc_tracemat(ARC_MATPRINTF,_UKF_XP_,_UKF_ANX_,1,13,4);
                arc_log(ARC_INFO,"arc_relpos : P-prior(ukf)=\n");
                arc_tracemat(ARC_MATPRINTF,_UKF_PP_,_UKF_ANX_,_UKF_ANX_,13,4);

                arc_ukf_filter_reset(ukf,_UKF_XP_,_UKF_PP_);

                /* compute weight of ukf */
                arc_ukf_filter_compute_weights(ukf,rtk->opt.ukf_alpha,
                                               _UKF_ANX_-3,rtk->opt.ukf_beta);
                /* make double-difference measurements */
                if (arc_ukf_dd_meas()!=_UKF_NV_) {  /* todo:here may have bugs */
                    arc_log(ARC_ERROR, "arc_relpos : no double-differenced measurements\n");
                    stat=SOLQ_NONE;
                }
                arc_log(ARC_INFO,"arc_relpos : double-differenced measurements(ukf)\n");
                arc_tracemat(ARC_MATPRINTF,_UKF_MEAS_,_UKF_NV_,1,10,4);

                /* states filter */
                if (arc_ukf_filter_update(ukf,_UKF_MEAS_,NULL,NULL,NULL)) {
                    /* updates states */
                    arc_ukf_get_updatax(ukf,xp,Pp);
                    arc_log(ARC_INFO,"ukf updates x : \n");
                    arc_tracemat(ARC_MATPRINTF,xp,3,1,10,4);
                }
                else {
                    arc_log(ARC_WARNING, "arc_relpos : ukf updates failed \n");
                    stat=SOLQ_NONE;
                }
#ifdef ARC_TEST
                fp_ukf_ceres<<setiosflags(ios::fixed)<<setprecision(10)
                            <<xp[0]<<"   "<<xp[1]<<"   "<<xp[2]<<"   "<<std::endl;
#endif
                /* free ukf problem */
                arc_ukf_filter_delete(ukf);  /* todo:why always free?may be have some better way to do */
            }
            else {
                arc_log(ARC_WARNING,"arc_relpos : no double-difference measurements for ukf \n");
                stat=SOLQ_NONE;
            }
        }
    }
    if (stat!=SOLQ_NONE&&arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk,NULL)) {

        /* post-fit residuals for float solution */
        nv=arc_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg);

        /* validation of float solution */
        if (arc_valpos(rtk,v,R,vflg,nv,ARC_SOLVALTHRES)) {

            /* update state and covariance matrix */
            arc_matcpy(rtk->x,xp,rtk->nx,1);
            arc_matcpy(rtk->P,Pp,rtk->nx,rtk->nx);

            /* update ambiguity control struct */
            rtk->sol.ns=0;
            for (i=0;i<ns;i++) for (f=0;f<nf;f++) {
                if (!rtk->ssat[sat[i]-1].vsat[f]) continue;
                rtk->ssat[sat[i]-1].lock[f]++;  /* lock counter of phase */
                rtk->ssat[sat[i]-1].outc[f]=0;  /* outage counter of phase is re-initial */
                if (f==0) rtk->sol.ns++; /* valid satellite count by L1 */
            }
            /* lack of valid satellites */
            if (rtk->sol.ns<4) stat=SOLQ_NONE;
        }
        else stat=SOLQ_NONE;
    }
    /* resolve integer ambiguity by LAMBDA */
    if (stat!=SOLQ_NONE&&(opt->amb_group?arc_resamb_group_LAMBDA(rtk,bias,xa)
                                        :arc_resamb_LAMBDA(rtk,bias,xa)>1)) {

        /* resolve integer ambiguity by LAMBDA success */
        if (arc_zdres(0,obs,nu,rs,dts,svh,nav,xa,opt,0,y,e,azel,rtk,NULL)) {

            /* post-fit reisiduals for fixed solution */
            nv=arc_ddres(rtk,nav,dt,xa,NULL,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg);

            /* validation of fixed solution */
            if (arc_valpos(rtk,v,R,vflg,nv,ARC_SOLVALTHRES)) {

                /* hold integer ambiguity */
                if (++rtk->nfix>=rtk->opt.minfix&&
                        (rtk->opt.modear==ARMODE_FIXHOLD||
                        rtk->opt.modear==ARMODE_INST)) {
                    arc_holdamb(rtk,xa,opt->amb_group);
                }
                stat=SOLQ_FIX;
            }
        }
    }
    /* save solution status */
    if (stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];  /* fix solutions */
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
    }
    else {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[1+2*rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];
        rtk->nfix=0;
    }
    for (i=0;i<n;i++) for (j=0;j<nf;j++) {
        if (obs[i].L[j]==0.0) continue;
        rtk->ssat[obs[i].sat-1].pt[obs[i].rcv-1][j]=obs[i].time;
        rtk->ssat[obs[i].sat-1].ph[obs[i].rcv-1][j]=obs[i].L[j];
    }
    for (i=0;i<n;i++) {  /* save rover and base station to satellite's distance */
        if (obs[i].rcv-1==0)
            rtk->ssat[obs[i].sat-1].r0[obs[i].rcv-1]=arc_geodist(rs+i*6,rtk->x, e+i*3);
        else
            rtk->ssat[obs[i].sat-1].r0[obs[i].rcv-1]=arc_geodist(rs+i*6,rtk->rb,e+i*3);
    }
    for (i=0;i<ns;i++) for (j=0;j<nf;j++) {
        /* output snr of rover receiver */
        rtk->ssat[sat[i]-1].snr[j]=obs[iu[i]].SNR[j];
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<nf;j++) {
        if (rtk->ssat[i].fix[j]==2&&stat==SOLQ_FIX) rtk->ssat[i].fix[j]=1;
        if (rtk->ssat[i].slip[j]&1) rtk->ssat[i].slipc[j]++;
    }
    if (D) free(D); if (Qv) free(Qv);
    free(rs); free(dts); free(var); free(y); free(e); free(azel);
    free(xp); free(Pp);  free(xa);  free(v); free(H); free(R); free(bias);

    if (stat!=SOLQ_NONE) rtk->sol.stat=stat;
    return stat!=SOLQ_NONE;
}
/* number of estimated states ------------------------------------------------*/
extern int arc_pppnx(const prcopt_t *opt)
{
    return NX(opt);
}
/* initialize rtk control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void arc_rtkinit(rtk_t *rtk, const prcopt_t *opt)
{
    sol_t sol0={{0}};
    ambc_t ambc0={{{0}}};
    ssat_t ssat0={0};
    int i;

    arc_log(ARC_INFO, "rtkinit :\n");

    rtk->sol=sol0;
    for (i=0;i<6;i++) rtk->rb[i]=0.0;
    rtk->nx=opt->mode<=PMODE_FIXED?NX(opt): arc_pppnx(opt);
    rtk->na=opt->mode<=PMODE_FIXED?NR(opt): arc_pppnx(opt);
    rtk->tt=0.0;
    rtk->x=arc_zeros(rtk->nx,1);
    rtk->P=arc_zeros(rtk->nx,rtk->nx);
    rtk->xa=arc_zeros(rtk->na,1);
    rtk->Pa=arc_zeros(rtk->na,rtk->na);
    rtk->nfix=rtk->neb=0;
    for (i=0;i<MAXSAT;i++) {
        rtk->ambc[i]=ambc0;
        rtk->ssat[i]=ssat0;
    }
    for (i=0;i<MAXERRMSG;i++) rtk->errbuf[i]=0;
    rtk->opt=*opt;

    /* ceres solver options */
    rtk->ceres_active_x=arc_imat(rtk->nx,1);  /* todo: may be use static variance */

    /* ambiguity solver options */
    for (i=0;i<MAXSAT;i++) rtk->amb_index[i]=0;
}
/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void arc_rtkfree(rtk_t *rtk)
{
    arc_log(ARC_INFO, "rtkfree :\n");

    rtk->nx=rtk->na=0;
    if (rtk->x)  free(rtk->x ); rtk->x =NULL;
    if (rtk->P)  free(rtk->P ); rtk->P =NULL;
    if (rtk->xa) free(rtk->xa); rtk->xa=NULL;
    if (rtk->Pa) free(rtk->Pa); rtk->Pa=NULL;
    if (rtk->ceres_active_x) {
        free(rtk->ceres_active_x); rtk->ceres_active_x=NULL;
    }
}
/* adjust option for single rtk precise positioning -------------------------*/
static void arc_adjust_option(rtk_t *rtk)
{
    if (rtk->opt.modear==ARMODE_FIXHOLD) rtk->opt.amb_part=ARC_AMB_NOPART;
    if (rtk->opt.mode  ==PMODE_KINEMA)   rtk->opt.modear  =ARMODE_INST;  /* todo:slip cycle often happen if rover station is kinemaed */
    if (rtk->opt.navsys==SYS_CMP)        rtk->opt.thresar[0]=1.5;  /* todo:when the Ratio value is 1.5,the effect is better than 3.0 */
}
/* arc single rtk precise positioning ---------------------------------------*/
extern int arc_srtkpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    sol_t solb={{0}};
    gtime_t time;
    int i,nu,nr;
    char msg[128]="";

    arc_log(ARC_INFO, "arc_srtkpos  : time=%s n=%d\n", time_str(obs[0].time, 3), n);
    arc_log(ARC_WARNING, "arc_srtkpos : obs=\n"); arc_traceobs(4,obs,n);

    /* set base staion position */
    if (opt->refpos<=POSOPT_RINEX&&opt->mode!=PMODE_SINGLE&&
        opt->mode!=PMODE_MOVEB) {
        for (i=0;i<6;i++) rtk->rb[i]=i<3?opt->rb[i]:0.0;
    }
    /* count rover/base station observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ; /* rover */
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ; /* base */

    time=rtk->sol.time; /* previous epoch */

    /* rover position by single point positioning */
    if (!arc_pntpos(obs,nu,nav,&rtk->opt,&rtk->sol,NULL,rtk->ssat,msg)) {
        arc_log(ARC_WARNING, "arc_srtkpos : point pos error (%s)\n", msg);
        if (!rtk->opt.dynamics) {
            return 0;
        }
    }
    if (time.time!=0) rtk->tt=timediff(rtk->sol.time,time);

    /* single point positioning */
    if (opt->mode==PMODE_SINGLE) {
        return 1;
    }
    /* suppress output of single solution */
    if (!opt->outsingle) {
        rtk->sol.stat=SOLQ_NONE;
    }
    /* check number of data of base station and age of differential */
    if (nr==0) {
        arc_log(ARC_ERROR, "arc_srtkpos : no base station observation data for rtk\n");
        return 1;
    }
    if (opt->mode==PMODE_MOVEB) { /* moving baseline */

        /* estimate position/velocity of base station */
        if (!arc_pntpos(obs+nu,nr,nav,&rtk->opt,&solb,NULL,NULL,msg)) {
            arc_log(ARC_WARNING, "arc_srtkpos : base station position error (%s)\n", msg);
            return 0;
        }
        rtk->sol.age=(float)timediff(rtk->sol.time,solb.time);

        if (fabs(rtk->sol.age)>TTOL_MOVEB) {
            arc_log(ARC_WARNING, "arc_srtkpos : time sync error "
                    "for moving-base (age=%.1f)\n", rtk->sol.age);
            return 0;
        }
        for (i=0;i<6;i++) rtk->rb[i]=solb.rr[i];

        /* time-synchronized position of base station */
        for (i=0;i<3;i++) rtk->rb[i]+=rtk->rb[i+3]*rtk->sol.age;
    }
    else {
        rtk->sol.age=(float)timediff(obs[0].time,obs[nu].time);

        if (fabs(rtk->sol.age)>opt->maxtdiff) {
            arc_log(ARC_WARNING, "arc_srtkpos : age of differential "
                    "error (age=%.1f)\n", rtk->sol.age);
            return 1;
        }
    }
    /* adjust option,this step just for optimize positioning process */
    if (opt->auto_ajust_opt) arc_adjust_option(rtk); /* todo:just for debug */

    /* difference pseudorange positioning for initialing */
    if (opt->init_dc) arc_diff_pr_relpos(rtk,obs,nu,nr,nav); /* todo:may be have no effects,need to test */

    /* relative potitioning */
    arc_relpos(rtk,obs,nu,nr,nav);
    return 1;
}

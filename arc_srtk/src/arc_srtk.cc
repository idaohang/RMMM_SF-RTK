
#include "arc.h"
#include <iomanip>
#include <fstream>
#include <rtklib.h>

using namespace std;

/* constants/macros ----------------------------------------------------------*/
#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)
#define MAXSTATES   80         /* max numbers of states */

#define VAR_POS     SQR(10.0)  /* initial variance of receiver pos (m^2) */
#define VAR_GRA     SQR(0.001) /* initial variance of gradient (m^2) */
#define VAR_AMB     SQR(10.0)  /* initial variance of ambguity (cycle^2) */
#define INIT_ZWD    0.15       /* initial zwd (m) */

#define GAP_RESION  120        /* gap to reset ionosphere parameters (epochs) */
#define VAR_HOLDAMB 0.001      /* constraint to hold ambiguity (cycle^2) */

#define TTOL_MOVEB  (1.0+2*DTTOL)
                               /* time sync tolerance for moving-baseline (s) */
#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
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

/* constants/global variables for ceres problem--------------------------------------------*/
#define  _WINDOWSSIZE_   10                         /* ceres problem sovler windows size */
static double *_H_=NULL;                            /* ceres problem jacobi matrix pointor */
static int _NX_=0;                                  /* ceres problem parameters block numbers*/
static int _NV_=0;                                  /* ceres problem residual block numbers */
static rtk_t *_RTK_;                                /* ceres problem rtk data type */
static const obsd_t*_OBS_;                          /* ceres problem observation data */
static const nav_t *_NAV_;                          /* ceres problem navigation data */
static int _NU_;                                    /* ceres problem rover observation data */
static int _NR_;                                    /* ceres problem base observation data */
static double *_RS_;                                /* ceres problem satellite position/velecity */
static double *_DTS_;                               /* ceres problem satellite clock drift */
static double *_Y_;                                 /* ceres problem undifferenced residuals for rover and base */
static double *_AZEL_;                              /* ceres problem satellite azel */
static double *_E_;                                 /* ceres problem rover station - satellite's line of light */
static int *_SVH_;                                  /* ceres problem satellite health flag */
static int *_VFLAG_;                                /* ceres problem satellite valid flag */
static int *_ParaBlock_=NULL;                       /* ceres problem parameters block list */
static double **_Para_=NULL;                        /* ceres problem parameters block pointor */
static int _Para_Const_List_[MAXSTATES]={-1};       /* ceres problem const parameters list */
static int _NCP_=0;                                 /* ceres problem const parameters numbers */
static double *_XP_=NULL;                           /* ceres problem prior states */
static double *_PP_=NULL;                           /* ceres problem prior covariance matrix */
static double *_X_=NULL;                            /* ceres problem states */
static double *_R_=NULL;                            /* ceres problem measurement covariance matrix */

static double _WINDOWS_X_[_WINDOWSSIZE_][_WINDOWSSIZE_*3+MAXSAT];
                                                    /* ceres problem solver windows sates */
static int _WINDOWS_FREAME_COUNT=0;                 /* ceres problem how many frame in cunrrent time */

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
    double sinel=sin(el);
    int i=sys==SYS_GLO?1:(sys==SYS_GAL?2:0),nf=NF(opt);

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
    return 2.0*(opt->ionoopt==IONOOPT_IFLC?3.0:1.0)*(a*a+b*b/sinel/sinel+c*c)+d*d;
}
/* baseline length -----------------------------------------------------------*/
static double arc_baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return arc_norm(dr, 3);
}
/* initialize state and covariance -------------------------------------------*/
static void arc_initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void arc_testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;

    arc_log(3,"testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
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
static double arc_yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* shadow-crossing GPS IIA ---------------------------------------------------*/
static int arc_yaw_shadow_IIA(double beta, double mu, double eps0, double R,
                              double mudot, double *yaw)
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
static int arc_yaw_shadow_GLO(double beta, double mu, double eps0, double R,
                              double mudot, double *yaw)
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
static int arc_yaw_noon(double beta, double mu, double beta0, double R,
                    double mudot, double *yaw)
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
static int arc_yaw_midnight(double beta, double mu, double beta0, double R,
                            double mudot, double *yaw)
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
static int arc_yaw_IIA(int sat, int opt, double beta, double mu, double *yaw)
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
static int arc_yaw_IIR(int sat, int opt, double beta, double mu, double *yaw)
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
static int arc_yaw_IIF(int sat, int opt, double beta, double mu, double *yaw)
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
static int arc_yaw_GAL(int sat, int opt, double beta, double mu, double *yaw)
{
    *yaw=arc_yaw_nominal(beta,mu);
    return 1;
}

static int arc_yaw_CMP(int sat, int opt, double beta, double mu, double *yaw)
{
    *yaw=0.0;
    return 1;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int arc_yaw_angle(int sat, const char *type, int opt, double beta, double mu,
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
static int arc_sat_yaw(gtime_t time, int sat, const char *type, int opt,
                       const double *rs, double *exs, double *eys)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0};
    int i;

    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);

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
static int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw)
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
static int arc_selsat(const obsd_t *obs, double *azel, int nu, int nr,
                      const prcopt_t *opt, int *sat, int *iu, int *ir)
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
static void arc_udpos(rtk_t *rtk, double tt)
{
    int i;

    arc_log(ARC_INFO, "arc_udpos   : tt=%.3f\n", tt);

    /* fixed mode */
    if (rtk->opt.mode==PMODE_FIXED) {
        for (i=0;i<3;i++) arc_initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (arc_norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_POS,i);
    }
    /* reset rover station pistion and its variance */
    for (i=0;i<3;i++) arc_initx(rtk,rtk->sol.rr[i],rtk->sol.qr[i],i);
    
    /* update ceres solver problem active states index list */
    for (i=0;i<3;i++) rtk->ceres_active_x[i]=1;

    /* add velecity noise to position variance */
    for (i=0;i<3;i++) rtk->P[i+i*rtk->nx]+=SQR(rtk->opt.prn[5])*tt;
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void arc_udion(rtk_t *rtk, double tt, double bl, const int *sat, int ns)
{
    double el,fact;
    int i,j;

    arc_log(ARC_INFO, "arc_udion   : tt=%.1f bl=%.0f ns=%d\n", tt, bl, ns);

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
static void arc_udtrop(rtk_t *rtk, double tt, double bl)
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
/* detect cycle slip by LLI --------------------------------------------------*/
static void arc_detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv)
{
    unsigned int slip,LLI;
    int f=0,sat=obs[i].sat;

    arc_log(ARC_INFO, "arc_detslp_ll: i=%d rcv=%d\n", i, rcv);

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
/* temporal update of phase biases -------------------------------------------*/
static void arc_udbias(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav)
{
    double cp,pr,*bias,offset,lami;
    int i,j,f=0,slip,reset;

    arc_log(ARC_INFO,"arc_udbias  : tt=%.1f ns=%d\n",tt,ns);

    for (i=0;i<ns;i++) {

        /* detect cycle slip by LLI */
        rtk->ssat[sat[i]-1].slip[f]&=0xFC;
        arc_detslp_ll(rtk,obs,iu[i],1);
        arc_detslp_ll(rtk,obs,ir[i],2);

        /* update half-cycle valid flag */
        rtk->ssat[sat[i]-1].half[f]=
                !((obs[iu[i]].LLI[f]&2)||(obs[ir[i]].LLI[f]&2));
    }
    /* reset phase-bias if instantaneous AR or expire obs outage counter */
    for (i=1;i<=MAXSAT;i++) {

        reset=++rtk->ssat[i-1].outc[f]>(unsigned int)rtk->opt.maxout;
        if (rtk->opt.modear==ARMODE_INST&&rtk->x[IB(i,f,&rtk->opt)]!=0.0) {
            arc_initx(rtk,0.0,VAR_AMB,IB(i,f,&rtk->opt));
        }
        else if (reset&&rtk->x[IB(i,f,&rtk->opt)]!=0.0) {
            arc_initx(rtk,0.0,VAR_AMB,IB(i,f,&rtk->opt));
            arc_log(ARC_INFO, "arc_udbias : obs outage counter overflow (sat=%3d L%d n=%d)\n",
                    i,f+1,rtk->ssat[i-1].outc[f]);
        }
        if (rtk->opt.modear!=ARMODE_INST&&reset) {
            rtk->ssat[i-1].lock[f]=-rtk->opt.minlock;
        }
    }
    /* reset phase-bias if detecting cycle slip */
    for (i=0;i<ns;i++) {
        j=IB(sat[i],f,&rtk->opt);
        rtk->P[j+j*rtk->nx]+=rtk->opt.prn[0]*rtk->opt.prn[0]*tt;
        slip=rtk->ssat[sat[i]-1].slip[f];
        if (rtk->opt.modear==ARMODE_INST||!(slip&1)) continue;
        rtk->x[j]=0.0;
        rtk->ssat[sat[i]-1].lock[f]=-rtk->opt.minlock;
    }
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
    /* correct phase-bias offset to enssure phase-code coherency */
    if (j>0) {
        for (i=1;i<=MAXSAT;i++) {
            if (rtk->x[IB(i,f,&rtk->opt)]!=0.0) rtk->x[IB(i,f,&rtk->opt)]+=offset/j;
        }
    }
    /* set initial states of phase-bias */
    for (i=0;i<ns;i++) {
        if (bias[i]==0.0||rtk->x[IB(sat[i],f,&rtk->opt)]!=0.0) continue;
        arc_initx(rtk,bias[i],SQR(rtk->opt.std[0]),IB(sat[i],f,&rtk->opt));
    }
    free(bias);
}
/* temporal update of states --------------------------------------------------*/
static void arc_udstate(rtk_t *rtk, const obsd_t *obs, const int *sat,
                        const int *iu, const int *ir, int ns, const nav_t *nav)
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
        arc_udbias(rtk,tt,obs,sat,iu,ir,ns,nav);
    }
}
/* undifferenced phase/code residual for satellite ---------------------------*/
static void arc_zdres_sat(int base, double r, const obsd_t *obs, const nav_t *nav,
                          const double *azel, const double *dant,double dion,
                          double vion,const prcopt_t *opt, double *y,const rtk_t *rtk)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=0,nf=1;

    if (lam[i]==0.0) return;

    /* check snr mask */
    if (testsnr(base,i,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) return;

    /* residuals = observable - pseudorange,think about phase windup correction */
    if (obs->L[i]!=0.0) y[i   ]=obs->L[i]*lam[i]-r-dant[i]+dion-rtk->ssat[obs->sat-1].phw*lam[i];
    if (obs->P[i]!=0.0) y[i+nf]=obs->P[i]       -r-dant[i]-dion;
}
/* undifferenced phase/code residuals ----------------------------------------*/
static int arc_zdres(int base, const obsd_t *obs, int n, const double *rs,
                     const double *dts, const int *svh, const nav_t *nav,
                     const double *rr, const prcopt_t *opt, int index, double *y,
                     double *e, double *azel, rtk_t* rtk)
{
    double r,rr_[3],pos[3],dant[NFREQ]={0},disp[3];
    double zhd,zazel[]={0.0,90.0*D2R},dion,vion;
    int i=0,nf=1;

    arc_log(ARC_INFO, "arc_zdres   : n=%d\n", n);

    for (i=0;i<n*nf*2;i++) y[i]=0.0;

    if (arc_norm(rr,3)<=0.0) return 0; /* no receiver position */

    for (i=0;i<3;i++) rr_[i]=rr[i];

    /* earth tide correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),rr_,opt->tidecorr,&nav->erp,
                 opt->odisp[base],disp);
        for (i=0;i<3;i++) rr_[i]+=disp[i];
    }
    ecef2pos(rr_,pos);

    for (i=0;i<n;i++) {
        /* compute geometric-range and azimuth/elevation angle */
        if ((r=geodist(rs+i*6,rr_,e+i*3))<=0.0) continue;
        if (satazel(pos,e+i*3,azel+i*2)<opt->elmin) continue;

        /* excluded satellite? */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;

        /* satellite clock-bias */
        r+=-CLIGHT*dts[i*2];

        /* troposphere delay model (hydrostatic) */
        zhd=tropmodel(obs[0].time,pos,zazel,0.0);
        r+=tropmapf(obs[i].time,pos,azel+i*2,NULL)*zhd;

        /* ionospheric corrections */
        if (!ionocorr(obs[i].time,nav,obs[i].sat,pos,azel+i*2,
                      IONOOPT_BRDC,&dion,&vion)) continue;
        /* receiver antenna phase center correction */
        antmodel(opt->pcvr+index,opt->antdel[index],azel+i*2,opt->posopt[1],dant);
        
        /* phase windup model */
        if (!model_phw(rtk->sol.time,obs[i].sat,nav->pcvs[obs[i].sat-1].type,
                       opt->posopt[2]?2:0,rs+i*6,rr,&rtk->ssat[obs[i].sat-1].phw)) {
            continue;
        }
        /* undifferenced phase/code residual for satellite */
        arc_zdres_sat(base,r,obs+i,nav,azel+i*2,dant,dion,vion,opt,y+i*nf*2,rtk);
    }
    arc_log(ARC_INFO, "arc_zdres : rr_=%.3f %.3f %.3f\n", rr_[0], rr_[1], rr_[2]);
    arc_log(ARC_INFO, "arc_zdres : pos=%.9f %.9f %.3f\n", pos[0] * R2D, pos[1] * R2D, pos[2]);
    for (i=0;i<n;i++) {
        arc_log(ARC_INFO, "arc_zdres : sat=%2d %13.3f %13.3f %13.3f %13.10f %6.1f %5.1f\n",
                obs[i].sat,rs[i*6],rs[1+i*6],rs[2+i*6],dts[i*2],azel[i*2]*R2D,
                azel[1+i*2]*R2D);
    }
    arc_log(ARC_INFO, "arc_zdres : y=\n");
    arc_tracemat(4,y,nf*2,n,13,3);

    return 1;
}
/* test valid observation data -----------------------------------------------*/
static int arc_validobs(int i, int j, int f, int nf, double *y)
{
    /* if no phase observable, psudorange is also unusable */
    return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
           (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
/* double-differenced measurement error covariance ---------------------------*/
static void arc_ddcov(const int *nb, int n, const double *Ri, const double *Rj,
                      int nv, double *R)
{
    int i,j,k=0,b=0;

    arc_log(ARC_INFO, "arc_ddcov   : n=%d\n", n);

    for (i=0;i<nv*nv;i++) R[i]=0.0;
    for (b=0;b<n;k+=nb[b++]) {
        for (i=0;i<nb[b];i++) for (j=0;j<nb[b];j++) {
            R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
        }
    }
    arc_log(ARC_INFO, "R=\n");
    arc_tracemat(5, R, nv, nv, 8, 6);
}
/* baseline length constraint ------------------------------------------------*/
static int arc_constbl(rtk_t *rtk, const double *x, const double *P, double *v,
                       double *H, double *Ri, double *Rj, int index)
{
    const double thres=0.1; /* threshold for nonliearity (v.2.3.0) */
    double xb[3],b[3],bb,var=0.0;
    int i;

    arc_log(ARC_INFO, "arc_constbl : \n");

    /* no constraint */
    if (rtk->opt.baseline[0]<=0.0) return 0;

    /* time-adjusted baseline vector and length */
    for (i=0;i<3;i++) {
        xb[i]=rtk->rb[i]+rtk->rb[i+3]*rtk->sol.age;
        b[i]=x[i]-xb[i];
    }
    bb= arc_norm(b, 3);

    /* approximate variance of solution */
    if (P) {
        for (i=0;i<3;i++) var+=P[i+i*rtk->nx];
        var/=3.0;
    }
    /* check nonlinearity */
    if (var>thres*thres*bb*bb) {
        arc_log(ARC_INFO, "arc_constbl : "
                "equation nonlinear (bb=%.3f var=%.3f)\n", bb, var);
        return 0;
    }
    /* constraint to baseline length */
    v[index]=rtk->opt.baseline[0]-bb;
    if (H) {
        for (i=0;i<3;i++) H[i+index*rtk->nx]=b[i]/bb;
    }
    Ri[index]=0.0;
    Rj[index]=SQR(rtk->opt.baseline[1]);

    arc_log(ARC_INFO, "baseline len   "
            "v=%13.3f R=%8.6f %8.6f\n", v[index], Ri[index], Rj[index]);
    return 1;
}
/* precise tropspheric model -------------------------------------------------*/
static double arc_prectrop(gtime_t time, const double *pos, int r,
                           const double *azel, const prcopt_t *opt, const double *x,
                           double *dtdx)
{
    double m_w=0.0,cotz,grad_n,grad_e;
    int i=IT(r,opt);

    /* wet mapping function */
    tropmapf(time,pos,azel,&m_w);

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
/* double-differenced phase/code residuals -----------------------------------*/
static int arc_ddres(rtk_t *rtk, const nav_t *nav, double dt, const double *x,
                     const double *P, const int *sat, double *y, double *e,
                     double *azel, const int *iu, const int *ir, int ns, double *v,
                     double *H, double *R, int *vflg)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3],didxi=0.0,didxj=0.0,*im;
    double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,lami,lamj,fi,fj,*Hi=NULL;
    int i,j,k,m,f,ff=0,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=1;

    arc_log(ARC_INFO, "arc_ddres   : dt=%.1f nx=%d ns=%d\n", dt, rtk->nx, ns);

    bl=arc_baseline(x,rtk->rb,dr);
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);

    Ri=arc_mat(ns*nf*2+2,1); Rj=arc_mat(ns*nf*2+2,1); im=arc_mat(ns,1);
    tropu=arc_mat(ns,1); tropr=arc_mat(ns,1); dtdxu=arc_mat(ns,3); dtdxr=arc_mat(ns,3);

    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].resp[0]=rtk->ssat[i].resc[0]=0.0;
    }
    /* compute factors of ionospheric and tropospheric delay */
    for (i=0;i<ns;i++) {
        if (opt->ionoopt>=IONOOPT_EST) {
            im[i]=(ionmapf(posu,azel+iu[i]*2)+ionmapf(posr,azel+ir[i]*2))/2.0;
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
            if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;
            if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
        }
        if (i<0) continue;

        /* make double difference */
        for (j=0;j<ns;j++) {
            if (i==j) continue;
            sysi=rtk->ssat[sat[i]-1].sys;
            sysj=rtk->ssat[sat[j]-1].sys;
            if (!arc_test_sys(sysj,m)) continue;
            if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;

            ff=f%nf;
            lami=nav->lam[sat[i]-1][ff];
            lamj=nav->lam[sat[j]-1][ff];
            if (lami<=0.0||lamj<=0.0) continue;
            if (H) {
                Hi=H+nv*rtk->nx;
                for (k=0;k<rtk->nx;k++) Hi[k]=0.0;
            }
            /* double-differenced residual */
            v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])-
                  (y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);

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
                v[nv]-=didxi*x[II(sat[i],opt)]-didxj*x[II(sat[j],opt)];
                if (H) {
                    Hi[II(sat[i],opt)]= didxi;
                    Hi[II(sat[j],opt)]=-didxj;
                }
                /* updates ceres solver active states index list */
                rtk->ceres_active_x[II(sat[i],opt)]=1;
                rtk->ceres_active_x[II(sat[j],opt)]=1;
            }
            /* double-differenced tropospheric delay term */
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                v[nv]-=(tropu[i]-tropu[j])-(tropr[i]-tropr[j]);
                for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {
                    if (!H) continue;
                    Hi[IT(0,opt)+k]= (dtdxu[k+i*3]-dtdxu[k+j*3]);
                    Hi[IT(1,opt)+k]=-(dtdxr[k+i*3]-dtdxr[k+j*3]);
                    /* updates ceres solver active states index list */
                    rtk->ceres_active_x[IT(0,opt)+k]=1;
                    rtk->ceres_active_x[IT(1,opt)+k]=1;
                }
            }
            /* double-differenced phase-bias term */
            if (f<nf) {
                v[nv]-=lami*x[IB(sat[i],f,opt)]-lamj*x[IB(sat[j],f,opt)];
                if (H) {
                    Hi[IB(sat[i],f,opt)]= lami;
                    Hi[IB(sat[j],f,opt)]=-lamj;
                }
                /* updates ceres solver active states index list */
                rtk->ceres_active_x[IB(sat[i],f,opt)]=1;
                rtk->ceres_active_x[IB(sat[j],f,opt)]=1;
            }
            if (f<nf) rtk->ssat[sat[j]-1].resc[f   ]=v[nv];
            else      rtk->ssat[sat[j]-1].resp[f-nf]=v[nv];

            /* test innovation */
            if (opt->ceres==0) if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                if (f<nf) {
                    rtk->ssat[sat[i]-1].rejc[f]++;
                    rtk->ssat[sat[j]-1].rejc[f]++;
                }
                arc_log(ARC_WARNING, "arc_ddres : outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                        sat[i], sat[j], f < nf ? "L" : "P", f % nf + 1, v[nv]);
                continue;
            }
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
            arc_log(ARC_INFO,"arc_ddres : sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n", sat[i],
                    sat[j],f<nf?"L":"P",f%nf+1,v[nv],Ri[nv],Rj[nv]);

            vflg[nv++]=(sat[i]<<16)|(sat[j]<<8)|((f<nf?0:1)<<4)|(f%nf);
            nb[b]++;
        }
        b++;
    }
    /* end of system loop */

    /* baseline length constraint for moving baseline */
    if (opt->mode==PMODE_MOVEB&&arc_constbl(rtk,x,P,v,H,Ri,Rj,nv)) {
        vflg[nv++]=3<<4;
        nb[b++]++;
    }
    if (H) {
        arc_log(ARC_INFO, "arc_ddres : H=\n");
        arc_tracemat(ARC_MATPRINTF, H,rtk->nx,nv,7,4);}

    /* double-differenced measurement error covariance */
    if (R) arc_ddcov(nb,b,Ri,Rj,nv,R);

    free(Ri); free(Rj); free(im);
    free(tropu); free(tropr); free(dtdxu); free(dtdxr);

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

    satposs(time,obsb,nb,nav,opt->sateph,rs,dts,var,svh);

    if (!arc_zdres(1,obsb,nb,rs,dts,svh,nav,rtk->rb,opt,1,yb,e,azel,rtk)) {
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
    int i,j,k,m,f,nb=0,nx=rtk->nx,na=rtk->na,nf=1,nofix;

    arc_log(ARC_INFO, "arc_ddmat   :\n");

    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].fix[0]=0;
    }
    for (i=0;i<na;i++) D[i+i*nx]=1.0;

    for (m=0;m<4;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        nofix=(m==1&&rtk->opt.glomodear==0)
              ||(m==3&&rtk->opt.bdsmodear==0);
        for (f=0,k=na;f<nf;f++,k+=MAXSAT) {

            for (i=k;i<k+MAXSAT;i++) {
                if (rtk->x[i]==0.0||!arc_test_sys(rtk->ssat[i-k].sys,m)||
                    !rtk->ssat[i-k].vsat[f]||!rtk->ssat[i-k].half[f]) {
                    continue;
                }
                if (rtk->ssat[i-k].lock[f]>0&&!(rtk->ssat[i-k].slip[f]&2)&&
                    rtk->ssat[i-k].azel[1]>=rtk->opt.elmaskar&&!nofix) {
                    rtk->ssat[i-k].fix[f]=2; /* fix */
                    break;
                }
                else rtk->ssat[i-k].fix[f]=1;
            }
            for (j=k;j<k+MAXSAT;j++) {
                if (i==j||rtk->x[j]==0.0||!arc_test_sys(rtk->ssat[j-k].sys,m)||
                    !rtk->ssat[j-k].vsat[f]) {
                    continue;
                }
                if (rtk->ssat[j-k].lock[f]>0&&!(rtk->ssat[j-k].slip[f]&2)&&
                    rtk->ssat[i-k].vsat[f]&&
                    rtk->ssat[j-k].azel[1]>=rtk->opt.elmaskar&&!nofix) {
                    D[i+(na+nb)*nx]= 1.0;
                    D[j+(na+nb)*nx]=-1.0;
                    nb++;
                    rtk->ssat[j-k].fix[f]=2; /* fix */
                }
                else rtk->ssat[j-k].fix[f]=1;
            }
        }
    }
    arc_log(ARC_INFO, "D=\n");
    arc_tracemat(5, D, nx, na + nb, 2, 0);
    return nb;
}
/* restore single-differenced ambiguity --------------------------------------*/
static void arc_restamb(rtk_t *rtk, const double *bias, double *xa)
{
    int i,n,m,f,index[MAXSAT],nv=0,nf=NF(&rtk->opt);

    arc_log(ARC_INFO, "arc_restamb :\n");

    for (i=0;i<rtk->nx;i++) xa[i]=rtk->x [i];
    for (i=0;i<rtk->na;i++) xa[i]=rtk->xa[i];

    for (m=0;m<4;m++) for (f=0;f<nf;f++) {

        for (n=i=0;i<MAXSAT;i++) {
            if (!arc_test_sys(rtk->ssat[i].sys,m)
                ||rtk->ssat[i].fix[f]!=2) {
                continue;
            }
            index[n++]=IB(i+1,f,&rtk->opt);
        }
        if (n<2) continue;

        xa[index[0]]=rtk->x[index[0]];

        for (i=1;i<n;i++) {
            xa[index[i]]=xa[index[0]]-bias[nv++];
        }
    }
}
/* hold integer ambiguity ----------------------------------------------------*/
static void arc_holdamb(rtk_t *rtk, const double *xa)
{
    double *v,*H,*R;
    int i,n,m,f,info,index[MAXSAT],nb=rtk->nx-rtk->na,nv=0,nf=NF(&rtk->opt);

    arc_log(ARC_INFO, "arc_holdamb :\n");

    v= arc_mat(nb,1); H=arc_zeros(nb,rtk->nx);

    for (m=0;m<4;m++) for (f=0;f<nf;f++) {

        for (n=i=0;i<MAXSAT;i++) {
            if (!arc_test_sys(rtk->ssat[i].sys,m)
                ||rtk->ssat[i].fix[f]!=2||
                rtk->ssat[i].azel[1]<rtk->opt.elmaskhold) {
                continue;
            }
            index[n++]=IB(i+1,f,&rtk->opt);
            rtk->ssat[i].fix[f]=3; /* hold */
        }
        /* constraint to fixed ambiguity */
        for (i=1;i<n;i++) {
            v[nv]=(xa[index[0]]-xa[index[i]])
                  -(rtk->x[index[0]]-rtk->x[index[i]]);
            H[index[0]+nv*rtk->nx]= 1.0;
            H[index[i]+nv*rtk->nx]=-1.0;
            nv++;
        }
    }
    if (nv>0) {
        R= arc_zeros(nv,nv);
        for (i=0;i<nv;i++) R[i+i*nv]=VAR_HOLDAMB;

        /* update states with constraints */
        if ((info= arc_filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv))) {
            arc_log(ARC_WARNING, "filter error (info=%d)\n", info);
        }
        free(R);
    }
    free(v); free(H);
}
/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
static int arc_resamb_LAMBDA(rtk_t *rtk, double *bias, double *xa)
{
    prcopt_t *opt=&rtk->opt;
    int i,j,ny,nb,info,nx=rtk->nx,na=rtk->na;
    double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*QQ,s[2];

    arc_log(ARC_INFO, "arc_resamb_LAMBDA : nx=%d\n", nx);

    rtk->sol.ratio=0.0;

    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
        rtk->opt.thresar[0]<1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D= arc_zeros(nx, nx);
    if ((nb=arc_ddmat(rtk,D))<=0) {
        arc_log(ARC_WARNING, "arc_resamb_LAMBDA : no valid double-difference\n");
        free(D);
        return 0;
    }
    ny=na+nb;
    y=arc_mat(ny,1); Qy=arc_mat(ny,ny);DP=arc_mat(ny,nx);
    b=arc_mat(nb,2); db=arc_mat(nb,1); Qb=arc_mat(nb,nb);
    Qab=arc_mat(na,nb); QQ=arc_mat(na,nb);

    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    arc_matmul("TN",ny,1,nx,1.0, D,rtk->x,0.0,y);
    arc_matmul("TN",ny,nx,nx,1.0,D,rtk->P,0.0,DP);
    arc_matmul("NN",ny,ny,nx,1.0,DP,D,0.0,Qy);

    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]=Qy[na+i+(na+j)*ny];
    for (i=0;i<na;i++) for (j=0;j<nb;j++) Qab[i+j*na]=Qy[   i+(na+j)*ny];

    arc_log(ARC_INFO, "arc_resamb_LAMBDA : N(0)=");
    arc_tracemat(ARC_MATPRINTF,y+na,1,nb,10,3);

    /* lambda/mlambda integer least-square estimation */
    if (!(info=lambda(nb,2,y+na,Qb,b,s))) {

        arc_tracemat(ARC_MATPRINTF,Qb,nb,nb,10,3);
        arc_log(ARC_INFO, "N(1)=");
        arc_tracemat(ARC_MATPRINTF,b,1,nb,10,3);
        arc_log(ARC_INFO, "N(2)=");
        arc_tracemat(ARC_MATPRINTF,b+nb,1,nb,10,3);

        rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;

        /* validation by popular ratio-test */
        if (s[0]<=0.0||s[1]/s[0]>=opt->thresar[0]) {

            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i=0;i<na;i++) {
                rtk->xa[i]=rtk->x[i];
                for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
            }
            for (i=0;i<nb;i++) {
                bias[i]=b[i];
                y[na+i]-=b[i];
            }
            if (!arc_matinv(Qb,nb)) {
                arc_matmul("NN",nb,1,nb,1.0,Qb,y+na,0.0,db);
                arc_matmul("NN",na,1,nb,-1.0,Qab,db,1.0,rtk->xa);

                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                arc_matmul("NN",na,nb,nb,1.0,Qab,Qb,0.0,QQ);
                arc_matmul("NT",na,na,nb,-1.0,QQ,Qab,1.0,rtk->Pa);

                arc_log(ARC_INFO, "arc_resamb : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                        nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);
                /* restore single-differenced ambiguity */
                arc_restamb(rtk,bias,xa);
            }
            else nb=0;
        }
        else { /* validation failed */
            arc_log(ARC_WARNING, "arc_resamb_LAMBDA : ambiguity validation "
                            "failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                    nb,s[1]/s[0],s[0],s[1]);
            nb=0;
        }
    }
    else {
        arc_log(ARC_WARNING, "lambda error (info=%d)\n", info);
    }
    free(D); free(y); free(Qy); free(DP);
    free(b); free(db); free(Qb); free(Qab); free(QQ);

    return nb; /* number of ambiguities */
}
/* validation of solution ----------------------------------------------------*/
static int arc_valpos(rtk_t *rtk, const double *v, const double *R, const int *vflg,
                      int nv, double thres)
{
    double fact=thres*thres;
    int i,stat=1,sat1,sat2,type,freq;
    char stype[8];

    arc_log(ARC_INFO, "arc_valpos  : nv=%d thres=%.1f\n", nv, thres);

    /* post-fit residual test */
    for (i=0;i<nv;i++) {
        if (v[i]*v[i]<=fact*R[i+i*nv]) continue;
        sat1=(vflg[i]>>16)&0xFF;
        sat2=(vflg[i]>> 8)&0xFF;
        type=(vflg[i]>> 4)&0xF;
        freq=vflg[i]&0xF;
        strcpy(stype,type==0?"L":(type==1?"L":"C"));
        arc_log(ARC_WARNING, "arc_valpos : "
                        "large residual (sat=%2d-%2d %s%d v=%6.3f sig=%.3f)\n",
                sat1,sat2,stype,freq+1,v[i],SQRT(R[i+i*nv]));
    }
    return stat;
}
/* Q parameters for adaptive Kaman filter-------------------------------------*/
static int arc_adap_Q(const rtk_t *rtk,double *Q,int n)
{
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
static int arc_adap_C0(const rtk_t* rtk ,const double *v,double *C0,int m,
                       double lam)
{
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
extern int adap_kaman_filter(rtk_t* rtk,double *x, double *P, const double *H,
                             const double *v,const double *R,int n,int m)
{
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
static void arc_ukf_filterfunc(ukf_t *ukf,double *input_X,double *output_X)
{
    int i;
    for (i=0;i<ukf->state_dim;i++) output_X[i]=input_X[i];
}
/* ukf: Propagate sigma points through the measurement model -----------------*/
static void arc_ukf_measfunc(ukf_t *ukf,double *input_X,double *output_Y)
{

}
/* ceres problem pointor asignment -------------------------------------------*/
static void arc_ceres_init(double *H,const prcopt_t* opt,int nv,double* rs,
                           double *dts,double *y,double *azel,int nu,int nr,
                           double *e,int *svh,int *vflag,rtk_t *rtk,
                           const obsd_t *obs,const nav_t* nav,double *R)
{
    int i;
    /* set the jacobi matric pointor to ceres problem */
    _H_=H;
    /* set the ceres parameters block numbers */
    _NX_=NX(opt);
    /* set the ceres residual block numbers */
    _NV_=nv;
    /* set satellite position/clock pointor */
    _RS_=rs; _DTS_=dts;
    /* set the satellite azimuth/elevation angle */
    _AZEL_=azel;
    /* set the undifferenced residuals for rover and base */
    _Y_=y;
    /* get the rover and base station obsevation numbers */
    _NU_=nu; _NR_=nr;
    /* rover station - satellite's line of light */
    _E_=e;
    /* satellite health flag */
    _SVH_=svh;
    /* satellite valid flag */
    _VFLAG_=vflag;
    /* observation data/navigation data pointor */
    _RTK_=rtk; _OBS_=obs; _NAV_=nav;

    /* measurements covariance matrix */
    _R_=R;

    /* initial parameters block list */
    if (_ParaBlock_==NULL) {
        _ParaBlock_=(int*)malloc(sizeof(int)*rtk->nx);
        for (i=0;i<rtk->nx;i++) _ParaBlock_[i]=1;
    }
    /* initial ceres problem prior information */
    if (_XP_==NULL) {
        _XP_=arc_mat(_NX_,1); _PP_=arc_mat(_NX_,_NX_);
    }
    /* initial ceres problem states */
    if (_X_==NULL) {
        _X_=arc_mat(_NX_,1);   /* allocate memory to ceres problem states */
    }
    for (i=0;i<_NX_;i++) _X_[i]=_RTK_->x[i]; /* ceres problem states asignment */
}
/* asignment parameters block ------------------------------------------------*/
static double** arc_ceres_para(const rtk_t* rtk)
{
    int i;
    double **para;

    /* allocate memeory */
    for (i=0;i<MAXSTATES;i++) para[i]=(double*)malloc(sizeof(double));

    /* asigne parameter pointor */
    for (i=0;i<_NX_;i++) para[i]=_X_+i;
    return para;
}
/* check whether there is some parameter can be set to be const parameters ----*/
static int arc_para_chk(const rtk_t* rtk,const double *H,const double *x)
{
    int i,k=0;
    for (i=0,k=0;i<_NX_;i++) {
        if (_RTK_->ceres_active_x[i]==0) _Para_Const_List_[k++]=i;
    }
    return k;
}
/* ceres solve problem cost function -----------------------------------------*/
static int arc_ceres_residual(void* m,double** parameters,double* residuals,
                              double** jacobians)
{
    prcopt_t *opt=&_RTK_->opt;
    gtime_t time=_OBS_[0].time;
    double *v,*R,dt,*L,*H;
    int ns,ny,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],i,j;

    /* adjust parameters positions */
    for (i=0;i<_NX_;i++) _X_[i]=parameters[i][0];

    /* select common satellites between rover and base-station */
    if ((ns=arc_selsat(_OBS_,_AZEL_,_NU_,_NR_,opt,sat,iu,ir))<=0) {
        arc_log(ARC_WARNING, "ceres_residual : no common satellite\n");
        return 0;
    }
    ny=ns*2+2;
    v=arc_mat(1,ny); R=arc_mat(ny,ny);
    dt=timediff(time,_OBS_[_NU_].time);

    /* undifferenced residuals for rover */
    if (!arc_zdres(0,_OBS_,_NU_,_RS_,_DTS_,_SVH_,_NAV_,_X_,opt,0,_Y_,_E_,_AZEL_,_RTK_)) {
        arc_log(ARC_WARNING, "ceres_residual : rover initial position error\n");
        free(v); free(R);
        return 0;
    }
    /* double-differenced residuals and partial derivatives */
    if ((_NV_=arc_ddres(_RTK_,_NAV_,dt,_X_,NULL,sat,_Y_,_E_,
                        _AZEL_,iu,ir,ns,v,_H_,R,_VFLAG_))<1) {
        arc_log(ARC_WARNING, "ceres_residual : no double-differenced residual\n");
        free(v); free(R);
        return 0;
    }
    if (opt->ceres_cholesky) {
        /* cholesky decomposition */
        L=arc_cholesky(R,_NV_); H=arc_mat(_NV_,_NX_);
        /* asign the double-difference residual */
        if (residuals) arc_matmul("NN",_NV_,1,_NV_,1.0,L,v,0.0,residuals);
        /* get the double-difference jacobi matrix */
        arc_matmul("NT",_NV_,_NX_,_NV_,1.0,L,_H_,0.0,H);
    }
    else {
        H=_H_;
        if (residuals) arc_matcpy(residuals,v,_NV_,1);
    }
    if (jacobians) {
        if (opt->ceres_cholesky) {
            /* asignment jacobi matrix , this case is thinking about the covariance */
            for (i=0;i<_NX_;i++) for (j=0;j<_NV_;j++) {
                if (jacobians[i]) jacobians[i][j]=-H[i*_NV_+j];
            }
        }
        else for (i=0;i<_NX_;i++) {
            /* asignment jacobi matrix , don't think about the covariance */
            if (jacobians[i]) for (j=0;j<_NV_;j++) jacobians[i][j]=-H[j*_NX_+i];
        }
    }
    if (opt->ceres_cholesky) {
        free(L); free(H);
    }
    free(v); free(R);
    return 1;
}
/* ceres sovel problem covariance matrix -------------------------------------*/
static int arc_ceres_cov(ceres_summary_t* summay,double *P)
{
    int i,j,k,col,row,*ix=NULL;
    double *J=arc_ceres_get_jacobis(summay,&row,&col),*Pp=NULL,*JR=NULL,*R=NULL;

    ix=arc_imat(_NX_,1); for (i=0,k=0;i<_NX_;i++) if (_RTK_->ceres_active_x[i]) ix[k++]=i;

    /* compute the covariance matrix */
    Pp=arc_mat(row,row);
    if (_RTK_->opt.ceres_cholesky) {
        arc_matmul("TN",col,col,row,1.0,J,J,0.0,Pp);
    }
    else if (_RTK_->opt.ceres_cholesky==0) {
        JR=arc_mat(col,row); R=arc_mat(_NV_,_NV_);
        arc_matcpy(R,_R_,_NV_,_NV_);
        if (!arc_matinv(R,_NV_)) {
            arc_matmul("TN",col,row,row,1.0,J,_R_,0.0,JR);
            arc_matmul("NN",col,col,row,1.0,JR,J,0.0,Pp);
        }
    }
    for (i=0;i<k;i++) for (j=0;j<k;j++) P[ix[i]*_NX_+ix[j]]=Pp[i*k+j];
    if (ix) free(ix); if (Pp) free(Pp);
    if (JR) free(JR); if (R)  free(R);
    return 1;
}
/* relative positioning ------------------------------------------------------*/
ofstream fp_pf_ceres("/home/sujinglan/arc_rtk/arc_test/data/gps_bds/static/arc_ceres_pos");
static int arc_relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                      const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,*v,*H,*R,*xp,*Pp,*xa,*bias,dt;
    int i,j,f,n=nu+nr,ns,ny,nv,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],niter;
    int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
    int stat=rtk->opt.mode<=PMODE_DGPS?SOLQ_DGPS:SOLQ_FLOAT;
    int nf=1;

    arc_log(ARC_INFO, "arc_relpos  : nx=%d nu=%d nr=%d\n", rtk->nx, nu, nr);

    dt=timediff(time,obs[nu].time);

    rs=arc_mat(6,n); dts=arc_mat(2,n);
    var=arc_mat(1,n); y=arc_mat(nf*2,n); e=arc_mat(3,n);
    azel= arc_zeros(2,n);

    /* initial satellite valid flag and snr */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        rtk->ssat[i].vsat[0]=0;
        rtk->ssat[i].snr [0]=0;
    }
    /* reset ceres problem solver active states index list */
    for (i=0;i<rtk->nx;i++) rtk->ceres_active_x[i]=0;

    /* satellite positions/clocks */
    satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        arc_testeclipse(obs,n,nav,rs);
    }

    /* undifferenced residuals for base station */
    if (!arc_zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,svh+nu,nav,rtk->rb,opt,1,
               y+nu*nf*2,e+nu*3,azel+nu*2,rtk)) {
        arc_log(ARC_WARNING, "arc_relpos : initial base station position error\n");

        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
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
    /* temporal update of states */
    arc_udstate(rtk,obs,sat,iu,ir,ns,nav);

    xp=arc_mat(rtk->nx,1); Pp=arc_zeros(rtk->nx,rtk->nx);
    xa=arc_mat(rtk->nx,1);
    arc_matcpy(xp,rtk->x,rtk->nx,1);

    ny=ns*nf*2+2;
    v=arc_mat(ny,1); H=arc_zeros(rtk->nx,ny);
    R=arc_mat(ny,ny); bias=arc_mat(rtk->nx,1);

    /* add 2 iterations for baseline-constraint moving-base */
    niter=opt->niter+(opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0?2:0);

    if (opt->ceres==0) for (i=0;i<niter;i++) {
        /* undifferenced residuals for rover */
        if (!arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk)) {
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
        arc_matcpy(Pp,rtk->P,rtk->nx,rtk->nx);
        /* adaptive Kaman filter */
        if (opt->adapt_filter) {
            if (!adap_kaman_filter(rtk,xp,Pp,H,v,R,rtk->nx,nv)) {
                arc_log(ARC_WARNING, "arc_relpos : adaptive filter error (info=%d)\n", info);
                stat=SOLQ_NONE;
                break;
            }
        }
        /* kalman filter measurement update */
        else {
            if ((info=arc_filter(xp,Pp,H,v,R,rtk->nx,nv))) {
                arc_log(ARC_WARNING, "arc_relpos : filter error (info=%d)\n", info);
                stat=SOLQ_NONE;
                break;
            }
            arc_log(ARC_INFO, "arc_relpos : x(%d)=", i+1);
        }
    }
    else if (opt->ceres==ARC_CERES_SINGLE) {

        /* initial ceres problem covariance matrix */
        arc_matcpy(Pp,rtk->P,rtk->nx,rtk->nx);

        /* undifferenced residuals for rover */
        if (!arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk)) {
            arc_log(ARC_WARNING, "arc_relpos : rover initial position error\n");
            stat=SOLQ_NONE;
        }
        /* double-differenced residuals and partial derivatives */
        if ((nv=arc_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg))<1) {
            arc_log(ARC_WARNING, "arc_relpos : no double-differenced residual\n");
            stat=SOLQ_NONE;
        }
        /* build a ceres solve problem */
        ceres_problem_t *ceres_problem=arc_ceres_create_problem();
        ceres_option_t  *ceres_option=arc_ceres_create_option();
        ceres_summary_t *ceres_summary=arc_ceres_create_summary();

        /* ceres problem initial */
        arc_ceres_init(H,opt,nv,rs,dts,y,azel,nu,nr,e,svh,vflg,rtk,obs,nav,R);

        /* get the ceres problem parameter pointor */
        _Para_=arc_ceres_para(rtk);
        /* add parameters block to ceres problem */
        arc_ceres_add_para_block(ceres_problem,_NX_,_Para_);

        /* check parameter list ,and get some parameters can be const */
        if ((_NCP_=arc_para_chk(rtk,H,_X_))>0) {
            for (i=0;i<_NCP_;i++) arc_ceres_set_para_const(ceres_problem,_Para_[_Para_Const_List_[i]]);
        }
        /* create a ceres solve problem to add all the double-difference residuals. */
        arc_ceres_problem_add_residual_block(
                ceres_problem,
                arc_ceres_residual,    /* Cost function */
                NULL,                  /* Points to the measurement */
                NULL,                  /* No loss function */
                ceres_create_huber_loss_function_data(1.0),
                                       /* loss function user data */
                _NV_,                  /* Number of residuals */
                _NX_,                  /* Number of parameter blocks */
                _ParaBlock_,           /* NUmber of parameter size */
                _Para_);               /* Parameters pointor */
        /* solve the problem */
        arc_ceres_solvex(ceres_problem,ceres_summary,ceres_option);

        /* copy ceres problem states to xp and covariance matrix to Pp */
        if (xp) arc_matcpy(xp,_X_,_NX_,1);
        if (Pp) arc_ceres_cov(ceres_summary,Pp);

        arc_tracemat(ARC_MATPRINTF,Pp,rtk->nx,rtk->nx,10,4);

        /* free ceres problem */
        arc_ceres_free_problem(ceres_problem);
        arc_ceres_free_option(ceres_option);
        arc_ceres_free_summary(ceres_summary);

        fp_pf_ceres<<setiosflags(ios::fixed)<<setprecision(10);
        for (int i=0;i<3;i++) fp_pf_ceres<<_X_[i]<<"  ";
        fp_pf_ceres<<std::endl;
    }
    else if (opt->ceres_windows==ARC_CERES_WINDOWS) {
        

    }
    if (stat!=SOLQ_NONE&&arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk)) {

        /* post-fit residuals for float solution */
        nv=arc_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg);

        /* validation of float solution */
        if (arc_valpos(rtk,v,R,vflg,nv,4.0)) {

            /* update state and covariance matrix */
            arc_matcpy(rtk->x,xp,rtk->nx,1);
            arc_matcpy(rtk->P,Pp,rtk->nx,rtk->nx);

            /* update ambiguity control struct */
            rtk->sol.ns=0;
            for (i=0;i<ns;i++) for (f=0;f<nf;f++) {
                if (!rtk->ssat[sat[i]-1].vsat[f]) continue;
                rtk->ssat[sat[i]-1].lock[f]++;
                rtk->ssat[sat[i]-1].outc[f]=0;
                if (f==0) rtk->sol.ns++; /* valid satellite count by L1 */
            }
            /* lack of valid satellites */
            if (rtk->sol.ns<4) stat=SOLQ_NONE;
        }
        else stat=SOLQ_NONE;
    }
    /* resolve integer ambiguity by LAMBDA */
    if (stat!=SOLQ_NONE&&arc_resamb_LAMBDA(rtk,bias,xa)>1) {

        if (arc_zdres(0,obs,nu,rs,dts,svh,nav,xa,opt,0,y,e,azel,rtk)) {

            /* post-fit reisiduals for fixed solution */
            nv=arc_ddres(rtk,nav,dt,xa,NULL,sat,y,e,azel,iu,ir,ns,v,NULL,R,vflg);

            /* validation of fixed solution */
            if (arc_valpos(rtk,v,R,vflg,nv,ARC_SOLVALTHRES)) {

                /* hold integer ambiguity */
                if (++rtk->nfix>=rtk->opt.minfix&&
                    rtk->opt.modear==ARMODE_FIXHOLD) {
                    arc_holdamb(rtk,xa);  /* todo:ambguity fix function need to improved */
                }
                stat=SOLQ_FIX;
            }
        }
    }
    /* save solution status */
    if (stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];
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
    for (i=0;i<ns;i++) for (j=0;j<nf;j++) {

        /* output snr of rover receiver */
        rtk->ssat[sat[i]-1].snr[j]=obs[iu[i]].SNR[j];
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<nf;j++) {
        if (rtk->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) rtk->ssat[i].fix[j]=1;
        if (rtk->ssat[i].slip[j]&1) rtk->ssat[i].slipc[j]++;
    }
    free(rs); free(dts); free(var); free(y); free(e); free(azel);
    free(xp); free(Pp);  free(xa);  free(v); free(H); free(R); free(bias);

    if (stat!=SOLQ_NONE) rtk->sol.stat=stat;
    return stat!=SOLQ_NONE;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
    return NX(opt);
}
/* initialize rtk control ------------------------------------------------------
* initialize rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkinit(rtk_t *rtk, const prcopt_t *opt)
{
    sol_t sol0={{0}};
    ambc_t ambc0={{{0}}};
    ssat_t ssat0={0};
    int i;

    arc_log(ARC_INFO, "rtkinit :\n");

    rtk->sol=sol0;
    for (i=0;i<6;i++) rtk->rb[i]=0.0;
    rtk->nx=opt->mode<=PMODE_FIXED?NX(opt):pppnx(opt);
    rtk->na=opt->mode<=PMODE_FIXED?NR(opt):pppnx(opt);
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
    rtk->ceres_active_x=arc_imat(rtk->nx,1);
}
/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkfree(rtk_t *rtk)
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
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;   /* rover */
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;   /* base */

    time=rtk->sol.time; /* previous epoch */

    /* rover position by single point positioning */
    if (!pntpos(obs,nu,nav,&rtk->opt,&rtk->sol,NULL,rtk->ssat,msg)) {
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
    if (opt->mode==PMODE_MOVEB) { /*  moving baseline */

        /* estimate position/velocity of base station */
        if (!pntpos(obs+nu,nr,nav,&rtk->opt,&solb,NULL,NULL,msg)) {
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
    /* relative potitioning */
    arc_relpos(rtk,obs,nu,nr,nav);
    return 1;
}

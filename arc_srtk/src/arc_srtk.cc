
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
#include <rtklib.h>

using namespace std;

/* constants/macros ----------------------------------------------------------*/
#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define VAR_POS     SQR(30.0)       /* initial variance of receiver pos (m^2) */
#define VAR_VEL     SQR(10.0)
#define VAR_GRA     SQR(0.001)      /* initial variance of gradient (m^2) */
#define INIT_ZWD    0.15            /* initial zwd (m) */
#define VAR_POS_AMB SQR(5.0)        /* after ambiguity fix,the variance of position */
#define VAR_CLKDRI  SQR(100.0)      /* initial variance of rover station clock drift  */

#define GAP_RESION  120             /* gap to reset ionosphere parameters (epochs) */
#define VAR_HOLDAMB 0.01            /* constraint to hold ambiguity (cycle^2) */
#define VAR_DOPPLER SQR(30.0)

#define TTOL_MOVEB  (1.0+2*DTTOL)
/* time sync tolerance for moving-baseline (s) */
#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define AMB_THRES   0.5             /* cycles */
#define FACTOR_RIRJ 6.0             /* thres of ratio of reference ambiguity variance and others */
#define MAXLARGERESC 10             /* max of large double-differecen residuals */
#define REVARFACTOR  1.1            /* reset the large double-difference residuals variance */
#define MAXAMBSQ     200
#define GOOGRATIO    5.0            /* a good ratio from fix ambiguity solutions */
#define FIXCOUNTC    300            /* min fix count to no-ambiguity double-difference solutions */
#define MAXDIFFAMB   300.0          /* max difference of precious epoch and current epoch for amb-fix */
#define NEWSIZE      10
#define MINTROP      10

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:6)
#define NPDC(opt)   ((opt)->dynamics_dc==0?3:6)
#define NCLK(opt)   ((opt)->est_doppler==0?0:1)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt)+NCLK(opt))
#define NX(opt)     (NR(opt)+NB(opt))
#define NXDC(opt)   (NPDC(opt)+NCLK(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+NCLK(opt)+(s)-1)       /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NCLK(opt)+NI(opt)+NT(opt)/2*(r))
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
/* tropos (r:0=rov,1:ref) */
#define IC(opt)     (NP(opt))                       /* rover station clock drift */
#define ICDC(opt)   (NPDC(opt))                     /* rover station clock drift */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)      /* phase bias (s:satno,f:freq) */

/* Macro for defining an exception------------------------------------------------------*/
ARC_DEFINE_EXCEPTION(Exception,std::runtime_error);

/* global variables for debuger---------------------------------------------------------*/
static int EPOCH=0;                                 /* epoch no. */

/* global variables ----------------------------------------------------------*/
static int statlevel=0;                             /* rtk status output level (0:off) */
static FILE *fp_stat=NULL;                          /* rtk status file pointer */
static char file_stat[1024]="";                     /* rtk status file original path */
static gtime_t time_stat={0};                       /* rtk status file time */
/* open solution status file ---------------------------------------------------
* open solution status file and set output level
* args   : char     *file   I   rtk status file
*          int      level   I   rtk status level (0: off)
* return : status (1:ok,0:error)
* notes  : file can constain time keywords (%Y,%y,%m...) defined in reppath().
*          The time to replace keywords is based on UTC of CPU time.
* output : solution status file record format
*
*   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          posx/posy/posz    : position x/y/z ecef (m) float
*          posxf/posyf/poszf : position x/y/z ecef (m) fixed
*
*   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          vele/veln/velu    : velocity e/n/u (m/s) float
*          acce/accn/accu    : acceleration e/n/u (m/s^2) float
*          velef/velnf/veluf : velocity e/n/u (m/s) fixed
*          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
*
*   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          clk1     : receiver clock bias GPS (ns)
*          clk2     : receiver clock bias GLO-GPS (ns)
*          clk3     : receiver clock bias GAL-GPS (ns)
*          clk4     : receiver clock bias BDS-GPS (ns)
*
*   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          sat      : satellite id
*          az/el    : azimuth/elevation angle(deg)
*          ion      : vertical ionospheric delay L1 (m) float
*          ion-fixed: vertical ionospheric delay L1 (m) fixed
*
*   $TROP,week,tow,stat,rcv,ztd,ztdf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          rcv      : receiver (1:rover,2:base station)
*          ztd      : zenith total delay (m) float
*          ztdf     : zenith total delay (m) fixed
*
*   $HWBIAS,week,tow,stat,frq,bias,biasf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          frq      : frequency (1:L1,2:L2,...)
*          bias     : h/w bias coefficient (m/MHz) float
*          biasf    : h/w bias coefficient (m/MHz) fixed
*
*   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
*          week/tow : gps week no/time of week (s)
*          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
*          az/el    : azimuth/elevation angle (deg)
*          resp     : pseudorange residual (m)
*          resc     : carrier-phase residual (m)
*          vsat     : valid data flag (0:invalid,1:valid)
*          snr      : signal strength (dbHz)
*          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
*          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
*          lock     : carrier-lock count
*          outc     : data outage count
*          slipc    : cycle-slip count
*          rejc     : data reject (outlier) count
*
*-----------------------------------------------------------------------------*/
extern int rtkopenstat(const char *file, int level)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];

    arc_log(ARC_INFO,"rtkopenstat: file=%s level=%d\n",file,level);

    if (level<=0) return 0;

    reppath(file,path,time,"","");

    if (!(fp_stat=fopen(path,"w"))) {
        arc_log(ARC_WARNING,"rtkopenstat: file open error path=%s\n",path);
        return 0;
    }
    strcpy(file_stat,file);
    time_stat=time;
    statlevel=level;
    return 1;
}
/* close solution status file --------------------------------------------------
* close solution status file
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkclosestat(void)
{
    arc_log(ARC_INFO,"rtkclosestat:\n");

    if (fp_stat) fclose(fp_stat);
    fp_stat=NULL;
    file_stat[0]='\0';
    statlevel=0;
}
/* write solution status to buffer -------------------------------------------*/
extern int rtkoutstat(rtk_t *rtk, char *buff)
{
    ssat_t *ssat;
    double tow,pos[3],vel[3],acc[3],vela[3]={0},acca[3]={0},xa[3];
    int i,j,week,est,nfreq,nf=NF(&rtk->opt);
    char id[32],*p=buff;

    if (rtk->sol.stat<=SOLQ_NONE) {
        return 0;
    }
    est=rtk->opt.mode>=PMODE_DGPS;
    nfreq=est?nf:1;
    tow=time2gpst(rtk->sol.time,&week);

    /* receiver position */
    if (est) {
        for (i=0;i<3;i++) xa[i]=i<rtk->na?rtk->xa[i]:0.0;
        p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                   rtk->sol.stat,rtk->x[0],rtk->x[1],rtk->x[2],xa[0],xa[1],
                   xa[2]);
    }
    else {
        p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                   rtk->sol.stat,rtk->sol.rr[0],rtk->sol.rr[1],rtk->sol.rr[2],
                   0.0,0.0,0.0);
    }
    /* receiver velocity and acceleration */
    if (est&&rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->x+3,vel);
        ecef2enu(pos,rtk->x+6,acc);
        if (rtk->na>=6) ecef2enu(pos,rtk->xa+3,vela);
        if (rtk->na>=9) ecef2enu(pos,rtk->xa+6,acca);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                   week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],acc[0],acc[1],
                   acc[2],vela[0],vela[1],vela[2],acca[0],acca[1],acca[2]);
    }
    else {
        ecef2pos(rtk->sol.rr,pos);
        ecef2enu(pos,rtk->sol.rr+3,vel);
        p+=sprintf(p,"$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                   week,tow,rtk->sol.stat,vel[0],vel[1],vel[2],
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    }
    /* receiver clocks */
    p+=sprintf(p,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
               week,tow,rtk->sol.stat,1,rtk->sol.dtr[0]*1E9,rtk->sol.dtr[1]*1E9,
               rtk->sol.dtr[2]*1E9,rtk->sol.dtr[3]*1E9);

    /* ionospheric parameters */
    if (est&&rtk->opt.ionoopt==IONOOPT_EST) {
        for (i=0;i<MAXSAT;i++) {
            ssat=rtk->ssat+i;
            if (!ssat->vs) continue;
            satno2id(i+1,id);
            j=II(i+1,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,id,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
                       rtk->x[j],xa[0]);
        }
    }
    /* tropospheric parameters */
    if (est&&(rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG)) {
        for (i=0;i<2;i++) {
            j=IT(i,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,i+1,rtk->x[j],xa[0]);
        }
    }
    /* receiver h/w bias */
    if (est&&rtk->opt.glomodear==2) {
        for (i=0;i<nfreq;i++) {
            j=IL(i,&rtk->opt);
            xa[0]=j<rtk->na?rtk->xa[j]:0.0;
            p+=sprintf(p,"$HWBIAS,%d,%.3f,%d,%d,%.4f,%.4f\n",week,tow,
                       rtk->sol.stat,i+1,rtk->x[j],xa[0]);
        }
    }
    return (int)(p-buff);
}
/* swap solution status file -------------------------------------------------*/
static void swapsolstat(void)
{
    gtime_t time=utc2gpst(timeget());
    char path[1024];

    if ((int)(time2gpst(time     ,NULL)/INT_SWAP_STAT)==
        (int)(time2gpst(time_stat,NULL)/INT_SWAP_STAT)) {
        return;
    }
    time_stat=time;

    if (!reppath(file_stat,path,time,"","")) {
        return;
    }
    if (fp_stat) fclose(fp_stat);

    if (!(fp_stat=fopen(path,"w"))) {
        arc_log(ARC_WARNING,"swapsolstat: file open error path=%s\n",path);
        return;
    }
    arc_log(ARC_INFO,"swapsolstat: path=%s\n",path);
}
/* output solution status ----------------------------------------------------*/
static void outsolstat(rtk_t *rtk)
{
    ssat_t *ssat;
    double tow;
    char buff[MAXSOLMSG+1],id[32];
    int i,j,n,week,nfreq,nf=NF(&rtk->opt);

    if (statlevel<=0||!fp_stat||!rtk->sol.stat) return;

    arc_log(ARC_INFO,"outsolstat:\n");

    /* swap solution status file */
    swapsolstat();

    /* write solution status */
    n=rtkoutstat(rtk,buff); buff[n]='\0';

    fputs(buff,fp_stat);

    if (rtk->sol.stat==SOLQ_NONE||statlevel<=1) return;

    tow=time2gpst(rtk->sol.time,&week);
    nfreq=rtk->opt.mode>=PMODE_DGPS?nf:1;

    /* write residuals and status */
    for (i=0;i<MAXSAT;i++) {
        ssat=rtk->ssat+i;
        if (!ssat->vs) continue;
        satno2id(i+1,id);
        for (j=0;j<nfreq;j++) {
            fprintf(fp_stat,"$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
                    week,tow,id,j+1,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
                    ssat->resp[j],ssat->resc[j],ssat->vsat[j],ssat->snr[j]*0.25,
                    ssat->fix[j],ssat->slip[j]&3,ssat->lock[j],ssat->outc[j],
                    ssat->slipc[j],ssat->rejc[j]);
        }
    }
}
/* single-differenced observable ---------------------------------------------*/
static double arc_sdobs(const obsd_t *obs, int i, int j, int f)
{
    double pi=f<NFREQ?obs[i].L[f]:obs[i].P[f-NFREQ];
    double pj=f<NFREQ?obs[j].L[f]:obs[j].P[f-NFREQ];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* snr measurement error variance model---------------------------------------*/
extern double arc_snr_varerr(const double dbhz,int f,int nf,const prcopt_t* opt)
{
#if 0
    return 0.0;
#endif
    static const double C=1.61E-2;
    if (f< nf) return C*pow(10.0,-dbhz/10.0)/(opt->eratio[f]*2.0);
    if (f>=nf) return C*pow(10.0,-dbhz/10.0);
}
/* single-differenced measurement error variance -----------------------------*/
static double arc_varerr(int sat, int sys, double el, double bl,double dt, int f,
                         const prcopt_t *opt)
{
    double a,b,c=opt->err[3]*bl/1E4,d=CLIGHT*opt->sclkstab*dt,fact=1.0;
    double sinel=sin(el),geo_factor=1.0;
    int i=sys==SYS_GLO?1:(sys==SYS_GAL?2:0),nf=NF(opt);

    /* bds geo */
    if (sys==SYS_CMP&&arc_is_bds_geo(sat)) geo_factor=EFACT_GEO;

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
/*----------------------------------------------------------------------------*/
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

    arc_log(ARC_INFO,"testeclipse:\n");

    /* unit vector of sun direction (ecef) */
    arc_sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
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

        arc_log(ARC_INFO,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
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
            arc_log(ARC_INFO,"(%2d) sat=%3d iu=%2d ir=%2d\n", k-1,obs[i].sat,i,j);
        }
    }
    return k;
}
/* change reference satellites------------------------------------------------*/
static int arc_chg_refsat(rtk_t *rtk,const nav_t *nav,int nv)
{
    arc_log(ARC_INFO,"arc_chg_refsat :\n");

    int i,chg=0;
    static char sys[NUMOFSYS][8]={"GPS","GLONASS","GAL","BDS"};

    for (i=0;i<NUMOFSYS;i++) { /* i=0:gps/qzs/sbs,1:glo,2:gal,3:bds */
        if (rtk->prefsat[i]==0) {
            rtk->prefsat[i]=rtk->refsat[i]; continue;
        }
        if (rtk->refsat[i]==rtk->prefsat[i]) continue;
        chg|=1;
        arc_log(ARC_WARNING,"arc_chg_refsat :%s,change reference satellite :%d->%d,epoch=%d",
                sys[i],rtk->prefsat[i],rtk->refsat[i],EPOCH);
    }
    return chg;
}
/*----------------------------------------------------------------------------*/
static void arc_chk_prefsat(rtk_t *rtk,const nav_t *nav,int nv,int *vf)
{
    int i,j,sat,sys;

    arc_log(ARC_INFO,"arc_chk_refsat :\n");

    for (i=0;i<NUMOFSYS;i++) vf[i]=-1;

    for (i=0;i<NUMOFSYS;i++) { /* i=0:gps/qzs/sbs,1:glo,2:gal,3:bds */
        for (j=0;j<nv;j++) {
            sat=rtk->sat[2*j+1]; sys=rtk->ssat[sat-1].sys;
            if (!arc_test_sys(sys,i)) continue; vf[i]=0;
            if (sat==rtk->prefsat[i]) {vf[i]=1;break;}
        }
    }
}
/* temporal update of position/velocity/acceleration -------------------------*/
static void arc_udpos(rtk_t *rtk,double tt)
{
    int i,j,nx=rtk->nx,ic=IC(&rtk->opt);
    double var=0.0;

    arc_log(ARC_INFO,"arc_udpos   : tt=%.3f\n",tt);

    /* fixed mode */
    if (rtk->opt.mode==PMODE_FIXED) {
        for (i=0;i<3;i++) arc_initx(rtk,rtk->opt.ru[i],1E-8,i);
        return;
    }
    /* initialize position for first epoch */
    if (arc_norm(rtk->x,3)<=0.0||rtk->sol.clk_jmp) {
        for (i=0;i<3;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        /* initial rover station velecity */
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        }
        /* initial rover station clock drift */
        if (rtk->opt.est_doppler) arc_initx(rtk,rtk->sol.clk_dri,
                                            VAR_CLKDRI,IC(&rtk->opt));
        return; /* initial compeleted */
    }
    /* dynamic mode for rover station */
    if (rtk->opt.dynamics) {
        /* add some system noice to states,here only velecity */
        for (i=3;i<6;i++) {
            if (rtk->x[i]==0.0) rtk->x[i]=rtk->opt.vel_prn[i-3];
        }
    }
    /* add clock drift noice */
    if (rtk->opt.est_doppler) {
        if (rtk->x[ic]==0.0) rtk->x[ic]=rtk->opt.clk_dri_prn;
    }
    /* update ceres solver problem active states index list,included velecity */
    if (rtk->opt.dynamics) {
        for (i=0;i<6;i++) rtk->ceres_active_x[i]=1;
    }
    else for (i=0;i<3;i++) rtk->ceres_active_x[i]=1;

    if (rtk->opt.est_doppler) rtk->ceres_active_x[ic]=1; /* clock drift */

    /* static mode */
    if (rtk->opt.mode==PMODE_STATIC) return;

    /* reset rover station pistion and its variance */
    /* todo:using standard positioning to initial ukf prior states and its covariacne matrix,
     * todo:and may have more better methids to do this */

    /* kinamic without dynamic mode */
    if (!rtk->opt.dynamics) { /* don't think about rover station velecity */

        if (rtk->opt.est_doppler) { /* clock drift initial */
            rtk->x[IC(&rtk->opt)]=rtk->sol.clk_dri;
        }
        if (rtk->opt.init_pnt) { /* use pnt to initial rtk */
            for (i=0;i<3;i++) {
                arc_initx(rtk,rtk->sol.rr[i],rtk->sol.qr[i],i); return;
            }
        }
        else if (rtk->opt.init_dc) { /* use dd-pseudorange to initial rtk */

            for (i=0;i<3;i++) rtk->P[i+i*nx]=rtk->sol.prsol.qr[i];
            rtk->P[1+nx*0]=rtk->P[0+nx*1]=rtk->sol.prsol.qr[3];
            rtk->P[2+nx*0]=rtk->P[0+nx*2]=rtk->sol.prsol.qr[5];
            rtk->P[1+nx*2]=rtk->P[2+nx*1]=rtk->sol.prsol.qr[4];
            arc_matcpy(rtk->x,rtk->sol.prsol.rr,3,1); /* rover station position */

            /* add clock drift system noice */
            if (rtk->opt.est_doppler) rtk->P[ic+ic*nx]+=SQR(rtk->opt.clk_dri_prn)*tt;

            arc_log(ARC_INFO,"arc_udpos: dd-pseudorange to initial rtk P=\n");
            arc_tracemat(ARC_MATPRINTF,rtk->P,nx,nx,10,4);
            return;
        }
        else {
            for (i=0;i<3;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_POS,i);

            arc_log(ARC_INFO,"arc_udpos: dd-pseudorange to initial rtk P=\n");
            arc_tracemat(ARC_MATPRINTF,rtk->P,nx,nx,10,4);
            return;
        }
    }
    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=rtk->P[i+i*nx]; var/=3.0;

    if (var>=VAR_POS||var<0.0) {
        /* reset position with large variance */
        for (i=0;i<3;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) arc_initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        arc_log(ARC_INFO,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
    /* state transition of position/velocity/acceleration for dynamic mode */
    double *F=arc_eye(nx),*FP=arc_mat(nx,nx),
            *xp=arc_mat(nx,1),pos[3]={0},Q[9]={0},Qv[9]={0};

    /* compute the state transition matrix */
    for (i=0;i<3;i++) F[i+(i+3)*nx]=tt;

    arc_log(ARC_INFO,"arc_udpos : state transition matrix =\n");
    arc_tracemat(ARC_MATPRINTF,F,nx,nx,10,4);

    arc_log(ARC_INFO,"arc_udpos : before state transition P =\n");
    arc_tracemat(ARC_MATPRINTF,rtk->P,nx,nx,10,4);

    /* x=F*x, P=F*P*F+Q */
    arc_matmul("NN",nx,1,nx,1.0,F,rtk->x,0.0,xp);
    arc_matcpy(rtk->x,xp,nx,1);
    arc_matmul("NN",nx,nx,nx,1.0,F,rtk->P,0.0,FP);
    arc_matmul("NT",nx,nx,nx,1.0,FP,F,0.0,rtk->P);

    arc_log(ARC_INFO,"arc_udpos : after state transition P = \n");
    arc_tracemat(ARC_MATPRINTF,rtk->P,nx,nx,10,4);

    /* process noise added to only velecity */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3]); Q[8]=SQR(rtk->opt.prn[4]);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
            rtk->P[i+3+(j+3)*rtk->nx]+=Qv[i+j*3];  /* add system process noise */
        }
    /* process noice of clock drift */
    if (rtk->opt.est_doppler) rtk->P[ic+ic*nx]+=SQR(rtk->opt.clk_dri_prn)*tt;

    free(F); free(FP); free(xp);
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void arc_udion(rtk_t *rtk,double tt,double bl,const int *sat,int ns)
{
    double el,fact;
    int i,j;

    arc_log(ARC_INFO,"arc_udion   : tt=%.1f bl=%.0f ns=%d\n",tt,bl,ns);

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
    double pos[3],azel[]={0.0,PI/2.0},zwd=INIT_ZWD,var;
    int i,j,k;
    static int trpc=0;

    arc_log(ARC_INFO,"udtrop  : tt=%.1f\n",tt);

    for (i=0;i<2;i++) {
        j=IT(i,&rtk->opt);

        if (rtk->x[j]<=0.0||trpc++>MINTROP) {

            if (i==0) ecef2pos(rtk->sol.rr,pos);
            if (i==1) ecef2pos(rtk->opt.rb,pos);
            sbstropcorr(rtk->sol.time,pos,azel,&var,&zwd);
            arc_initx(rtk,zwd,var,j); /* initial zwd */

            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) arc_initx(rtk,1E-6,VAR_GRA,++j);
            }
            trpc=0; /* reset counts of trop-holding */
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

    arc_matmul("TN",nv,1,4,1.0,H,dx,-1,v); /* v=H'*dx-v */
    arc_tracemat(ARC_MATPRINTF,v,1,nv,10,4);

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
    }
}
/* snr detect-----------------------------------------------------------------*/
static void arc_detsnr(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns)
{
    arc_log(ARC_INFO,"arc_detsnr :\n");

    int i,j=0,k=0,*ixb,*ixr;
    double snr_b[MAXSAT],snr_r[MAXSAT],asb=0.0,asr=0.0,r0_b,r0_r;
    static const double r=arc_re_norm(1.0-rtk->opt.snr_alpha/2.0);

    ixb=arc_imat(ns,1); ixr=arc_imat(ns,1);

    /* for base station */
    for (i=0;i<ns;i++) if (obs[ir[i]].SNR[0]) snr_b[j]=obs[ir[i]].SNR[0]*0.25,ixb[j++]=sat[i];
    /* for rover station */
    for (i=0;i<ns;i++) if (obs[iu[i]].SNR[0]) snr_r[k]=obs[iu[i]].SNR[0]*0.25,ixr[k++]=sat[i];

    if (j==0||k==0) {
        free(ixb); free(ixr); return; /* no observations */
    }

    for (i=0;i<j;i++) asb+=snr_b[i]; asb/=j;
    for (i=0;i<k;i++) asr+=snr_r[i]; asr/=k;

    for (i=0;i<j;i++) snr_b[i]-=asb; for (i=0;i<k;i++) snr_r[i]-=asr;

    arc_log(ARC_INFO,"snr_b=\n"); arc_tracemat(ARC_MATPRINTF,snr_b,1,j,10,4);
    arc_log(ARC_INFO,"snr_r=\n"); arc_tracemat(ARC_MATPRINTF,snr_r,1,k,10,4);

    arc_matmul("NT",1,1,j,1.0/(j-1),snr_b,snr_b,0.0,&r0_b);
    arc_matmul("NT",1,1,k,1.0/(k-1),snr_r,snr_r,0.0,&r0_r);

    /* detect snr for base station*/
    for (i=0;i<j;i++) {
        if ((fabs(snr_b[i])/SQRT(r0_b)>=r)) {
            arc_log(ARC_WARNING,"base station detect outlier: %3d \n",ixb[i]);
            rtk->ssat[ixb[i]-1].snrf[0]|=1; /* detected */
        }
    }
    /* detect snr for rover station */
    for (i=0;i<k;i++) {
        if ((fabs(snr_r[i])/SQRT(r0_r)>=r)) {
            arc_log(ARC_WARNING,"rover station detect outlier: %3d \n",ixr[i]);
            rtk->ssat[ixr[i]-1].snrf[0]|=1; /* detected */
        }
    }
    free(ixb); free(ixr);
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void arc_detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv)
{
    unsigned int slip,LLI;
    int f=0,sat=obs[i].sat;

    arc_log(ARC_INFO,"arc_detslp_ll: i=%d rcv=%d\n",i,rcv);

    if (obs[i].L[f]==0.0) return;

    /* restore previous LLI */
    if (rcv==1) LLI=getbitu(&rtk->ssat[sat-1].slip[f],0,2); /* rover */
    else        LLI=getbitu(&rtk->ssat[sat-1].slip[f],2,2); /* base  */

    /* detect slip by cycle slip flag in LLI */
    if (rtk->tt>=0.0) { /* forward */
        if (obs[i].LLI[f]&1) {
            arc_log(ARC_WARNING,"arc_detslp_ll : "
                            "slip detected forward (sat=%2d rcv=%d F=%d LLI=%x)\n",
                    sat,rcv,f+1,obs[i].LLI[f]);
        }
        slip=obs[i].LLI[f];
    }
    else { /* backward */
        if (LLI&1) {
            arc_log(ARC_WARNING,"arc_detslp_ll : "
                            "slip detected backward (sat=%2d rcv=%d F=%d LLI=%x)\n",
                    sat,rcv,f+1,LLI);
        }
        slip=LLI;
    }
    /* detect slip by parity unknown flag transition in LLI */
    if (((LLI&2)&&!(obs[i].LLI[f]&2))||(!(LLI&2)&&(obs[i].LLI[f]&2))) {
        arc_log(ARC_WARNING,"arc_detslp_ll : "
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

    arc_log(ARC_WARNING,"arc_detslp_ll : half-cycle, %2d \n",rtk->ssat[sat-1].half[f]);
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
    double cp,pr,*bias,offset=0.0,lami,*ba;
    int i,j,f=0,slip,reset,clk_jump=0;

    arc_log(ARC_INFO,"arc_udbias  : tt=%.1f ns=%d\n",tt,ns);

    /* snr check */
    if (rtk->opt.snr_det) arc_detsnr(rtk,obs,sat,iu,ir,ns);

    if (rtk->opt.reset_amb_all) {
        return arc_ubbias_all(rtk,tt,obs,sat,iu,ir,ns,nav);
    }
    for (i=0;i<ns;i++) {

        /* detect cycle slip by LLI */
        rtk->ssat[sat[i]-1].slip[f]&=0xFC; /* reset */
        arc_detslp_ll(rtk,obs,iu[i],1);
        arc_detslp_ll(rtk,obs,ir[i],2);

        /* update half-cycle valid flag */
        rtk->ssat[sat[i]-1].half[f]=
                !((obs[iu[i]].LLI[f]&2)||(obs[ir[i]].LLI[f]&2));
    }
    arc_detsl_new(rtk,obs,iu,ir,ns,nav,rs); /* detect cycle slip by new ways */

    arc_tracemat(ARC_MATPRINTF,rtk->x,1,rtk->nx,10,4);

    /* handle day-boundary clock jump */
    if (rtk->opt.posopt[5]) {
        clk_jump=(ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0);
    }
    clk_jump|=rtk->sol.clk_jmp; /* reciver clock jump handle */

    /* reset phase-bias if instantaneous AR or expire obs outage counter */
    for (i=1;i<=MAXSAT;i++) {

        reset=++rtk->ssat[i-1].outc[0]>(unsigned int)rtk->opt.maxout;
        if ((rtk->opt.modear==ARMODE_INST&&rtk->x[IB(i,0,&rtk->opt)]!=0.0)||
            clk_jump) {
            arc_initx(rtk,0.0,0.0,IB(i,0,&rtk->opt));  /* each epoch is re initialized */
        }
        else if (reset&&rtk->x[IB(i,0,&rtk->opt)]!=0.0) {  /* reset ambiguity */
            arc_initx(rtk,0.0,0.0,IB(i,0,&rtk->opt));
            arc_log(ARC_INFO,"arc_udbias : obs outage counter overflow (sat=%3d L%d n=%d)\n",
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
        if (rtk->opt.modear==ARMODE_INST||!(slip&1)||!rtk->ssat[sat[i]-1].snrf[0]) continue;
        rtk->x[j]=0.0;
        rtk->ssat[sat[i]-1].lock[0]=-rtk->opt.minlock;
        if (rtk->ssat[sat[i]-1].snrf[0]) rtk->ssat[sat[i]-1].snrc[0]++;
    }
    bias=arc_zeros(ns,1); ba=arc_zeros(ns,1);

    /* estimate approximate phase-bias by phase - code */
    for (i=j=0,offset=0.0;i<ns;i++) {

        cp=arc_sdobs(obs,iu[i],ir[i],0); /* cycle */
        pr=arc_sdobs(obs,iu[i],ir[i],0+NFREQ);
        lami=nav->lam[sat[i]-1][0];
        if (cp==0.0||pr==0.0||lami<=0.0) continue;

        /* single-difference ambiguity */
        bias[i]=cp-pr/lami;
        ba  [i]=rtk->x[IB(sat[i],0,&rtk->opt)];

        /* if is ARMODE_INST,then here it no process */
        if (rtk->x[IB(sat[i],f,&rtk->opt)]!=0.0) {
            offset+=bias[i]-rtk->x[IB(sat[i],f,&rtk->opt)];
            j++;
        }
    }
    arc_log(ARC_INFO,"bias=\n");
    arc_tracemat(ARC_MATPRINTF,bias,1,ns,10,4);
    arc_log(ARC_INFO,"ba=\n");
    arc_tracemat(ARC_MATPRINTF,ba,1,ns,10,4);

    /* correct phase-bias offset to enssure phase-code coherency */
    if (j>0) {
        for (i=1;i<=MAXSAT;i++) {
            if (rtk->x[IB(i,f,&rtk->opt)]!=0.0) rtk->x[IB(i,f,&rtk->opt)]+=offset/j;
        }
    }
    /* re-check ambiguity */
    if (rtk->opt.modear==ARMODE_FIXHOLD) for (i=0;i<ns;i++) {
        if (fabs(bias[i]-ba[i])>=rtk->opt.reset_hold) {
            arc_initx(rtk,bias[i],SQR(rtk->opt.std[0]),IB(sat[i],0,&rtk->opt));
        }
        arc_log(ARC_WARNING,"arc_udbias : re-check ambiguity \n");
    }
    /* set initial states of phase-bias */
    for (i=0;i<ns;i++) {
        if (bias[i]==0.0||rtk->x[IB(sat[i],0,&rtk->opt)]!=0.0) continue;
        arc_initx(rtk,bias[i],SQR(rtk->opt.std[0]),IB(sat[i],0,&rtk->opt));
    }
    arc_log(ARC_INFO,"arc_udbias : after ambiguity updates P= \n");
    arc_tracemat(ARC_MATPRINTF,rtk->P,rtk->nx,rtk->nx,10,4);
    free(bias); free(ba);
}
/* test valid observation data -----------------------------------------------*/
static int arc_validobs(int i,int j,int f,int nf,double *y)
{
    /* if no phase observable, psudorange is also unusable */
    return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
           (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
/* select the reference satellite for double-difference-----------------------*/
static int arc_sel_refsat(const rtk_t *rtk,const int *sat,int ns,const int *ir,
                          const int *iu,const double *azel,double *y,int *refsat)
{
    arc_log(ARC_INFO,"arc_sel_refsat :\n");

    int i,m,sys,j;

    for (m=0;m<NUMOFSYS;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        for (j=0,i=-1;j<ns;j++) { /* numbers of commom satellite */

            sys=rtk->ssat[sat[j]-1].sys;

            if (!arc_test_sys(sys,m)) continue;

            /* slip detect */
            if (rtk->ssat[sat[j]-1].slip[0]) continue;

            /* snr outlier */
            if (rtk->ssat[sat[j]-1].snrf[0]) continue;

            if (y) if (!arc_validobs(iu[j],ir[j],0,1,y)) continue;

            /* set the reference satellite index */
            if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
        }
        if (refsat) refsat[m]=i; /* reference satellite */
    }
}
/* get double-difference ambiguity -------------------------------------------*/
static ddamb_t *arc_get_ddamb(amb_t *bias,int sat1,int sat2)
{
    arc_log(ARC_INFO,"arc_get_ddamb: \n");

    int i; for (i=0;i<bias->nb;i++)
        if (sat1==bias->amb[i].sat1
            &&sat2==bias->amb[i].sat2) return bias->amb+i;
    return NULL;
}
/* add double-difference ambiguity-------------------------------------------*/
static  int arc_add_ddamb(amb_t *amb)
{
    ddamb_t *data;

    if (amb->nmax<=amb->nb) {
        if (amb->nmax<=0) amb->nmax=MAXSAT; else amb->nmax*=2;
        if (!(data=(ddamb_t *)realloc(amb->amb,sizeof(ddamb_t)*amb->nmax))) {
            arc_log(ARC_FATAL, "add double-difference ambiguity: memalloc error n=%dx%d\n",
                    sizeof(ddamb_t),amb->nmax);
            free(amb->amb); amb->amb=NULL; amb->nb=amb->nmax=0;
            return -1;
        }
        amb->amb=data;
    }
    return 1;
}
/* check P matrix and x matrix of rtk struct-----------------------------------*/
static int arc_chk_xP(rtk_t *rtk,int id)
{
    arc_log(ARC_INFO,"arc_chk_xP: \n");

    int i,j,pnx=rtk->nx,*ix;
    double *p=NULL;

    if (id<rtk->nx||!rtk->P||!rtk->x) return 1;

    arc_log(ARC_INFO,"old x=\n");
    arc_tracemat(ARC_MATPRINTF,rtk->x,rtk->nx,1,10,4);

    arc_log(ARC_INFO,"old P=\n");
    arc_tracemat(ARC_MATPRINTF,rtk->P,rtk->nx,rtk->nx,10,4);

    rtk->nx+=NEWSIZE;
    if (!(p=arc_zeros(1,rtk->nx))) {
        arc_log(ARC_FATAL,"add double-difference ambiguity: memalloc error n=%dx%d\n",
                sizeof(double),rtk->nx);
        free(rtk->x); rtk->x=NULL; rtk->nx=0;
        return -1;
    }
    arc_matcpy(p,rtk->x,1,pnx);
    free(rtk->x); rtk->x=p; /* new size */
    if (!(p=arc_zeros(rtk->nx,rtk->nx))) {
        arc_log(ARC_FATAL,"add double-difference ambiguity: memalloc error n=%dx%d\n",
                sizeof(double),rtk->nx);
        free(rtk->P); rtk->P=NULL; rtk->nx=0;
        return -1;
    }
    for (i=0;i<pnx;i++) for (j=0;j<pnx;j++) p[i+rtk->nx*j]=rtk->P[i+pnx*j];
    free(rtk->P); rtk->P=p; /* new size */

    if (!(ix=arc_imat(1,rtk->nx))) {
        arc_log(ARC_FATAL,"add double-difference ambiguity: memalloc error n=%dx%d\n",
                sizeof(double),rtk->nx);
        free(rtk->ceres_active_x); rtk->ceres_active_x=NULL; rtk->nx=0;
        return -1;
    }
    for (i=0;i<pnx;i++) ix[i]=rtk->ceres_active_x[i];
    free(rtk->ceres_active_x); rtk->ceres_active_x=ix; /* new size */

    arc_log(ARC_INFO,"new x=\n");
    arc_tracemat(ARC_MATPRINTF,rtk->x,rtk->nx,1,10,4);

    arc_log(ARC_INFO,"new P=\n");
    arc_tracemat(ARC_MATPRINTF,rtk->P,rtk->nx,rtk->nx,10,4);

    return 0; /* new size */
}
/* update double-difference ambiguity -----------------------------------------*/
static void arc_update_ddamb(rtk_t *rtk, const obsd_t *obs, const int *sat,
                             const int *iu, const int *ir, int ns, const nav_t *nav,
                             const double *rs,double *y,const double *azel)
{
    arc_log(ARC_INFO,"arc_update_ddamb :\n");

    int m,i,j,refsat[NUMOFSYS]={-1},sysi,sysj,ind=-1;
    ddamb_t *amb=NULL,amb0={0};
    double cp1,pr1,cp2,pr2,bias,lami,lamj;
    prcopt_t *opt=&rtk->opt;
    amb_t *pamb=&rtk->sol.bias; /* double-difference ambiguity pointor */

    for (i=0;i<ns;i++) {

        /* detect cycle slip by LLI */
        rtk->ssat[sat[i]-1].slip[0]&=0xFC; /* reset */
        arc_detslp_ll(rtk,obs,iu[i],1);
        arc_detslp_ll(rtk,obs,ir[i],2);

        /* update half-cycle valid flag */
        rtk->ssat[sat[i]-1].half[0]=
                !((obs[iu[i]].LLI[0]&2)||(obs[ir[i]].LLI[0]&2));
    }
    arc_detsl_new(rtk,obs,iu,ir,ns,nav,rs); /* detect cycle slip by new ways */

    /* select the reference satellite for double-difference */
    arc_sel_refsat(rtk,sat,ns,ir,iu,azel,y,refsat);

    /* updates double-difference ambiguity */
    for (m=0;m<NUMOFSYS;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */
        i=refsat[m]; /* reference satellite */
        if (i<0) continue;

        for (j=0;j<ns;j++) {
            if (i==j) continue;

            sysi=rtk->ssat[sat[i]-1].sys;
            sysj=rtk->ssat[sat[j]-1].sys;

            if (!arc_test_sys(sysj,m)) continue;
            if (!arc_test_sys(sysi,m)) continue;

            if (y) if (!arc_validobs(iu[j],ir[j],0,1,y)) continue;

            if (rtk->ssat[sat[j]-1].slip[0]) continue; /* slip detect */

            lami=nav->lam[sat[i]-1][0];
            lamj=nav->lam[sat[j]-1][0];

            cp1=arc_sdobs(obs,iu[i],ir[i],0); /* cycle */
            pr1=arc_sdobs(obs,iu[i],ir[i],0+NFREQ);
            cp2=arc_sdobs(obs,iu[j],ir[j],0); /* cycle */
            pr2=arc_sdobs(obs,iu[j],ir[j],0+NFREQ);

            if (cp1==0.0||pr1==0.0||cp2==0.0||pr2==0.0
                ||lami<=0.0||lamj<=0.0||lami!=lamj) continue;

            bias=(cp1-pr1/lami)-(cp2-pr2/lamj);

            if ((amb=arc_get_ddamb(pamb,sat[i],sat[j]))==NULL) {
                arc_add_ddamb(pamb);
                amb=&pamb->amb[pamb->nb]; /* new double-difference ambiguity */
                amb->id=pamb->nb++; /* id of this double-difference ambiguity */
                arc_chk_xP(rtk,amb->id+rtk->na); /* check x and P matrix */
            }
            amb->update=1; /* updates this ambiguity */
            amb->sat1=sat[i]; /* reference satellite */
            amb->sat2=sat[j]; /* another satellite */
            ind=amb->id+rtk->na; /*index of this ambiguity  */

            if (fabs(amb->b-bias)<2.0&&amb->c>=opt->minfix) {
                bias=amb->b; /* use last time fixed-ambiguity */
                arc_initx(rtk,bias,SQR(opt->std[0])/SQR(50.0),ind); /* fix and hold */
            }
            else arc_initx(rtk,bias,SQR(opt->std[0]),ind); /* slip-cycle */
        }
    }
    for (i=0;i<pamb->nb;i++) {
        arc_log(ARC_INFO,"%s : %3d-%3d update=%2d id=%2d \n",time_str(pamb->amb[i].t,2),
                pamb->amb[i].sat1,pamb->amb[i].sat2,pamb->amb[i].update,pamb->amb[i].id);
    }
}
/* temporal update of states --------------------------------------------------*/
static void arc_udstate(rtk_t *rtk, const obsd_t *obs, const int *sat,
                        const int *iu, const int *ir, int ns, const nav_t *nav,
                        const double *rs,double*y,const double *azel)
{
    double tt=fabs(rtk->tt),bl,dr[3];

    arc_log(ARC_INFO,"arc_udstate : ns=%d\n",ns);

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
        if (!rtk->opt.use_dd_sol) { /* general mode for update ambiguity */
            arc_udbias(rtk,tt,obs,sat,iu,ir,ns,nav,rs);
        }
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

    arc_log(ARC_INFO,"arc_zdres   : n=%d\n",n);

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

    for (i=0;i<n;i++) { /* loop for numbers of observations */
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
        arc_log(ARC_INFO,"arc_zdres : y=\n");
        arc_tracemat(ARC_MATPRINTF,y,nf*2,n,13,3);
    }
    if (ukf_y) {
        arc_log(ARC_INFO,"arc_zdres : ukf y=\n");
        arc_tracemat(ARC_MATPRINTF,ukf_y,nf*2,n,13,3);
    }
    return nzd;
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
        for (i=0;i<nb[b];i++) {
            for (j=0;j<nb[b];j++) {
                R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
            }
        }
    }
    arc_log(ARC_INFO,"R=\n");
    arc_tracemat(ARC_MATPRINTF,R,nv,nv,8,6);
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
static int arc_kalman_norm_inno(const rtk_t *rtk,const double *v,int nv,const double *H,
                                const double *R,double*ve)
{
    int i,nx=rtk->nx,n,info=1;
    double *He,*Pe,*HP,*Qe,*L=NULL;

    arc_log(ARC_INFO,"arc_kalman_norm_inno :\n");

    arc_log(ARC_INFO,"v=\n"); arc_tracemat(ARC_MATPRINTF,v,nv,1,10,4);

    He=arc_zeros(nx,nv); Pe=arc_zeros(nx,nx);

    n=arc_kalman_exct_xHP(rtk,H,He,nv,NULL,Pe);

    HP=arc_mat(nv,n); Qe=arc_mat(nv,nv);
    arc_matmul("TN",nv,n,n,1.0,He,Pe,0.0,HP);
    arc_matcpy(Qe,R,nv,nv);
    arc_matmul("NN",nv,nv,n,1.0,HP,He,1.0,Qe);

    arc_log(ARC_INFO,"Qe=\n");
    arc_tracemat(ARC_MATPRINTF,Qe,nv,nv,10,4);

    if (!(info=arc_matinv(Qe,nv))) {

        arc_log(ARC_INFO,"Qe-inv=\n");
        arc_tracemat(ARC_MATPRINTF,Qe,nv,nv,10,4);

        /* cholesky matrix decomposition */
        L=arc_cholesky(Qe,nv);

        arc_log(ARC_INFO,"L=\n");
        arc_tracemat(ARC_MATPRINTF,L,nv,nv,10,4);

        /* normalization */
        if (L&&ve) arc_matmul("NN",nv,1,nv,1.0,L,v,0.0,ve);

        arc_log(ARC_INFO,"ve=\n");
        arc_tracemat(ARC_MATPRINTF,ve,1,nv,10,4);
    }
    if (L) free(L); free(He); free(Pe); free(Qe);
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
    arc_log(ARC_INFO,"arc_robust_chk : \n");

    static const double rr=arc_re_chi2(nv,alpha);

    if (alpha<=0.0) {
        if (rk>chisqr[nv]) return SQRT(chisqr[nv]/rk); /* alpha=0.001 */
    }
    else {
        if (rk>=rr) return 0.95;/* robust */
    }
    return 1.0;
}
/* kalman filter Phi matrix---------------------------------------------------*/
static void arc_kalman_robust_phi(const rtk_t *rtk,const double *vn,int nv,
                                  int nx,double *phi,const double *Qv)
{
    int i,n=nv;
    double alpha=rtk->opt.kalman_robust_alpha;

    arc_log(ARC_INFO,"arc_kalman_phi :\n");

    for (i=0;i<nv;i++) phi[i+i*nv]=arc_robust_chk(n-1,alpha,fabs(Qv[i+i*nv]));

    arc_log(ARC_INFO,"Phi=\n");
    arc_tracemat(ARC_MATPRINTF,phi,nv,nv,10,4);
}
/* doppler partial derivatives by rover position------------------------------*/
static void arc_doppler_dr(const double *rr,const double *rs,const double *e,
                           double *dr,int dynamic,int estclk,int ic)
{
    arc_log(ARC_INFO,"arc_doppler_dr\n");

    int i;
    double dv[3],dp[3],pr,t;

    for (i=0;i<3;i++) dp[i]=rs[i]-rr[i]; pr=arc_norm(dp,3);
    for (i=0;i<3;i++) dv[i]=rs[i+3]-rr[i+3];

    t=pow(SQR(pr),1.5);

    /* doppler partial derivatives of rover station velecity term */
    if (dynamic) for (i=0;i<3;i++) dr[3+i]=-e[i];

    /* doppler partial derivatives of clock drift */
    if (estclk) dr[ic]=1.0;

    /* doppler partial derivatives of rover station term */
    dr[0]=2.0*dv[2]*dp[2]*dp[0]/(2.0*t)-dv[0]/pr+2.0*dv[0]*dp[0]*dp[0]/(2.0*t)+
          2.0*dv[1]*dp[1]*dp[0]/(2.0*t);
    dr[1]=2.0*dv[1]*dp[1]*dp[0]/(2.0*t)-dv[1]/pr+2.0*dv[2]*dp[2]*dp[1]/(2.0*t)+
          2.0*dv[0]*dp[0]*dp[1]/(2.0*t);
    dr[2]=2.0*dv[0]*dp[0]*dp[2]/(2.0*t)-dv[1]/pr+2.0*dv[1]*dp[1]*dp[2]/(2.0*t)+
          2.0*dv[2]*dp[2]*dp[2]/(2.0*t);
}
/* doppler residuals----------------------------------------------------------*/
static int arc_doppler_res(rtk_t *rtk,const double *x,const int *sat,int ns,
                           double *v,double *H,double *R,int nv,int nx,int dynamic,
                           int estclk,int ic)
{
    int i,j;
    ssat_t *psat=NULL;

    arc_log(ARC_INFO,"arc_doppler_res: \n");

    for (i=0,j=0;i<ns;i++) {

        if (rtk->ssat[sat[i]-1].doppler[0]==0.0) continue;
        psat=rtk->ssat+(sat[i]-1);

        /* doppler partial derivatives */
        if (H) arc_doppler_dr(x,psat->rs,psat->e,H+nv*nx,dynamic,estclk,ic);

        if (v) v[nv]=psat->doppler_res; /* doppler residuals */
        if (R) R[j ]=VAR_DOPPLER; /* measurement variance */

        arc_log(ARC_INFO,"arc_doppler_res :sat=%2d,v=%8.4lf,R=%8.4lf \n",
                sat[i],v==NULL?-999.0:v[nv],R==NULL?-999.0:R[j]);
        nv++; j++;
    }
    psat=NULL; return j;
}
/* no ambiguity double-differenced phase/code residuals-----------------------*/
static int arc_ddres_noamb(rtk_t *rtk,const nav_t *nav,double dt,const double *x,
                           const int *sat,double *y,double *e,double *azel,
                           int ns,double *v,double *H,double *R,int *vflg,
                           int *ir,int *iu)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3];
    double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,lami,lamj,*Hi=NULL;
    int i,j,k,m,f,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=1;
    int sat1=0,sat2=0,l;

    arc_log(ARC_INFO,"arc_ddres_noamb   : dt=%.1f nx=%d ns=%d\n",dt,rtk->nx,ns);

    bl=arc_baseline(x,rtk->rb,dr);
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);

    Ri=arc_mat(ns*nf*2+2,1); Rj=arc_mat(ns*nf*2+2,1);
    tropu=arc_mat(ns,1); tropr=arc_mat(ns,1);
    dtdxu=arc_mat(ns,3); dtdxr=arc_mat(ns,3);

    /* compute factors of ionospheric and tropospheric delay */
    for (i=0;i<ns;i++) {
        if (opt->tropopt>=TROPOPT_EST) {
            tropu[i]=arc_prectrop(rtk->sol.time,posu,0,azel+iu[i]*2,opt,x,dtdxu+i*3);
            tropr[i]=arc_prectrop(rtk->sol.time,posr,1,azel+ir[i]*2,opt,x,dtdxr+i*3);
        }
    }
    /* double-differenced phase/code residuals */
    for (m=0;m<4;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        for (f=0;f<nf;f++) { /* loop for frequencies,here no use code */

            for (l=0;l<rtk->bias.nb;l++) { /* loop for satellite */

                sat1=rtk->bias.amb[l].sat1; /* reference satellite */
                sat2=rtk->bias.amb[l].sat2; /* another satellite */

                /* ratio test and fix count test*/
                if (rtk->bias.amb[l].ratio<GOOGRATIO||
                    rtk->bias.amb[l].c<FIXCOUNTC) continue;

                /* time diffrence check */
                if (fabs(timediff(rtk->sol.time,
                                  rtk->bias.amb[l].t))>=MAXDIFFAMB) continue;

                if (sat1==0||sat2==0||sat1==sat2) continue; /* no double-difference satellite */

                /* detect slip */
                if (rtk->ssat[sat1-1].slip[0]||rtk->ssat[sat2-1].slip[0])
                    continue;

                /* get navigation system */
                sysi=satsys(sat1,NULL); sysj=satsys(sat2,NULL);

                if (!arc_test_sys(sysi,m)) continue; /* test system */
                if (!arc_test_sys(sysj,m)) continue;

                /* observation index */
                for (i=-1,k=0;k<ns;k++) if (sat[k]==sat1) {i=k;break;}
                for (j=-1,k=0;k<ns;k++) if (sat[k]==sat2) {j=k;break;}
                if (i<0||j<0) continue;

                /* check carrier wave lengths */
                lami=nav->lam[sat1-1][0];
                lamj=nav->lam[sat2-1][0];
                if (lami<=0.0||lamj<=0.0||lami!=lamj) continue;
                if (H) {
                    Hi=H+nv*rtk->nx;
                    for (k=0;k<rtk->nx;k++) Hi[k]=0.0;
                }
                /* double-differenced residual */
                if (v) v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])
                             -(y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);

                /* double-difference ambiguitu residuals */
                if (v) v[nv]-=lami*rtk->bias.amb[l].b;

                /* partial derivatives by rover position */
                if (H) {
                    for (k=0;k<3;k++) {
                        Hi[k]=-e[k+iu[i]*3]+e[k+iu[j]*3];
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
                /* test innovation */
                if (v) {
                    if (opt->maxinno>0.0
                        &&fabs(v[nv])>opt->maxinno) {
                        arc_log(ARC_WARNING,"arc_ddres_noamb : outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                                sat1,sat2,"L",1,v[nv]);
                        continue;
                    }
                }
                /* single-differenced measurement error variances */
                if (R) {
                    Rj[nv]=arc_varerr(sat2,sysj,azel[1+iu[j]*2],bl,dt,f,opt)+
                           arc_snr_varerr(rtk->ssat[sat[j]-1].snr[0]*0.25,f,nf,opt);
                    Ri[nv]=arc_varerr(sat1,sysi,azel[1+iu[i]*2],bl,dt,f,opt)+
                           arc_snr_varerr(rtk->ssat[sat[i]-1].snr[0]*0.25,f,nf,opt);
                }
                arc_log(ARC_INFO,"arc_ddres_noamb : "
                                "sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n",
                        sat1,sat2,"L",f%nf+1,v[nv],R==NULL?-999.0:Ri[nv],R==NULL?-999.0:Rj[nv]);

                if (vflg) vflg[nv]=(sat1<<16)|(sat2<<8)|((f<nf?0:1)<<4)|(f%nf);
                nb[b]++; nv++; /* increase double-residual index */
            }
        }
        b++;
    } /* end of system loop */

    if (H&&nv) {
        arc_log(ARC_INFO,"arc_ddres_noamb : H=\n");
        arc_tracemat(ARC_MATPRINTF,H,rtk->nx,nv,7,4);
    }
    /* double-differenced measurement error covariance */
    if (R&&nv) {
        arc_ddcov(nb,b,Ri,Rj,nv,R);
    }
    free(Ri); free(Rj); free(tropu);
    free(tropr); free(dtdxu); free(dtdxr);
    return nv;
}
/* double-differenced phase/code residuals -----------------------------------*/
static int arc_ddres(rtk_t *rtk,const nav_t *nav,double dt,const double *x,
                     const double *P,const int *sat,double *y,double *e,
                     double *azel,const int *iu,const int *ir,int ns,double *v,
                     double *H,double *R,int *vflg)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3],didxi=0.0,didxj=0.0,*im;
    double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,lami,lamj,fi,fj,*Hi=NULL,*RR=NULL,R1,R2;
    int i,j,k,m,f,ff=0,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=1,nd=0;

    arc_log(ARC_INFO,"arc_ddres   : dt=%.1f nx=%d ns=%d\n",dt,rtk->nx,ns);

    bl=arc_baseline(x,rtk->rb,dr);
    if (bl<=0.0) {
        arc_log(ARC_WARNING,"base station and rover station position is error\n");
        return 0;
    }
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

                /* slip detect */
                if (rtk->ssat[sat[j]-1].slip[0]) continue;

                /* snr outlier */
                if (rtk->ssat[sat[j]-1].snrf[0]) continue;

                /* y[i] may be zero,so must check out,this step is importance */
                if (y) if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;

                /* set the reference satellite index */
                if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
            }
            if (i<0) continue; /* i is the reference satellite */
            rtk->refsat[m]=sat[i]; /* double-difference reference satellite */

            /* make double difference */
            for (j=0;j<ns;j++) {
                if (i==j) continue;  /* sat[i] is reference satellite */
                sysi=rtk->ssat[sat[i]-1].sys;
                sysj=rtk->ssat[sat[j]-1].sys;
                if (!arc_test_sys(sysj,m)) continue;

                if (y) if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;

                /* if detect snr outlier */
                if (rtk->ssat[sat[i]-1].snrf[0]||rtk->ssat[sat[j]-1].snrf[0]) continue;

                /* check factor of ratio of reference ambiguity variance and others */
                R1=arc_varerr(sat[j],sysj,azel[1+iu[j]*2],bl,dt,f,opt);
                R2=arc_varerr(sat[i],sysi,azel[1+iu[i]*2],bl,dt,f,opt);
                if ((R1/R2)>=FACTOR_RIRJ) {
                    arc_log(ARC_WARNING,"check factor of ratio of reference "
                                    "ambiguity variance and others is failed,sat=%3d-%3d,factor=%8.4lf\n",
                            sat[i],sat[j],R1/R2);
                    continue;
                }
                /* check carrier wave lengths */
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
                /* test innovation for carrier-phase and pseudorange */
                if (v) {
                    if (f<nf&&opt->maxinnoc>0.0
                        &&fabs(v[nv])>opt->maxinnoc) { /* for carrier-phase */

                        rtk->ssat[sat[i]-1].rejc[f]++;
                        rtk->ssat[sat[j]-1].rejc[f]++;

                        arc_log(ARC_WARNING,"arc_ddres : outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                                sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
                        continue;
                    }
                    else if (f==nf&&opt->maxinno>0.0&&
                             fabs(v[nv])>opt->maxinno) { /* for pseudorange */
                        arc_log(ARC_WARNING,"arc_ddres : outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                                sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
                        continue;
                    }
                }
                if (v) {
                    if (f<nf) rtk->ssat[sat[j]-1].resc[f   ]=v[nv];
                    else      rtk->ssat[sat[j]-1].resp[f-nf]=v[nv];
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
                    Rj[nv]=arc_varerr(sat[j],sysj,azel[1+iu[j]*2],bl,dt,f,opt)+
                           arc_snr_varerr(rtk->ssat[sat[j]-1].snr[0]*0.25,f,nf,opt);
                    Ri[nv]=arc_varerr(sat[i],sysi,azel[1+iu[i]*2],bl,dt,f,opt)+
                           arc_snr_varerr(rtk->ssat[sat[i]-1].snr[0]*0.25,f,nf,opt);
                }
                /* adjust for larger double-differecen resiudal's satellite */
                if (rtk->ssat[sat[j]-1].large_resc[0]) {
                    Rj[nv]*=REVARFACTOR;  Ri[nv]*=REVARFACTOR;
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
                nb[b]++; nv++; /* increase double-residual index */
            }
            b++;
        }
    /* end of system loop */

    if (H&&nv) {
        arc_log(ARC_INFO,"arc_ddres : H=\n");
        arc_tracemat(ARC_MATPRINTF,H,rtk->nx,nv,7,4);
    }
    /* double-differenced measurement error covariance */
    if (R&&nv) {
        arc_ddcov(nb,b,Ri,Rj,nv,R);
    }
    /* check whther change reference satellite */
    if (arc_chg_refsat(rtk,nav,nv)) {

        int vf[NUMOFSYS]; arc_chk_prefsat(rtk,nav,nv,vf);

        arc_log(ARC_INFO,"vf=\n");
        arc_tracemati(ARC_MATPRINTF,vf,1,NUMOFSYS,2,1);

        for (i=0;i<NUMOFSYS;i++) {
            if (vf[i]==-1) continue;
            if (vf[i]== 0) rtk->prefsat[i]=rtk->refsat[i];
            if (vf[i]== 1) {
                if ((rtk->ref_delay[i]++)>opt->amb_ref_delayc) {
                    rtk->prefsat  [i]=rtk->refsat[i];
                    rtk->ref_delay[i]=0; /* reset */
                    arc_log(ARC_WARNING,
                            "arc_ddres : change reference satellite when delay count get\n");
                }
            }
        }
        arc_log(ARC_INFO,"check whther change reference satellite is done\n");
    }
    arc_log(ARC_INFO,"arc_ddres : active states index=\n");
    arc_tracemati(ARC_MATPRINTF,rtk->ceres_active_x,1,rtk->nx,2,1);

    /* doppler partial derivatives */
    if (opt->est_doppler) {

        /* get doppler residuals */
        nd=arc_doppler_res(rtk,x,sat,ns,v,H,Ri,nv,rtk->nx,
                           rtk->opt.dynamics,rtk->opt.est_doppler,IC(opt));
        nv+=nd; /* numbers of residuals */

        /* rerset covariance of active states */
        if (R&&nv) {
            RR=arc_zeros(nv,nv);
            for (i=0;i<nv-nd;i++) {
                for (j=0;j<nv-nd;j++) RR[i+j*nv]=R[i+j*(nv-nd)];
            }
            for (i=0;i<nd;i++) RR[nv-nd+i+(nv-nd+i)*nv]=Ri[i];

            /* save measurements variance matrix */
            arc_matcpy(R,RR,nv,nv); free(RR);

            arc_log(ARC_INFO,"arc_ddres : add doppler,R=\n");
            arc_tracemat(ARC_MATPRINTF,R,nv,nv,10,4);
        }
        if (H&&nv) {
            arc_log(ARC_INFO,"arc_ddres : add doppler,H=\n");
            arc_tracemat(ARC_MATPRINTF,H,rtk->nx,nv,8,4);
        }
    }
    free(Ri); free(Rj); free(im);
    free(tropu); free(tropr); free(dtdxu); free(dtdxr);
    return nv;
}
/* double-difference residuals in estimating double-difference ambiguity-------*/
static int arc_ddres_ddamb(rtk_t *rtk,const nav_t *nav,double dt,const double *x,
                           const double *P,const int *sat,double *y,double *e,
                           double *azel,const int *iu,const int *ir,int ns,
                           double *v,double *H,double *R,int *vflg)
{
    arc_log(ARC_INFO,"arc_ddres_ddamb \n");

    prcopt_t *opt=&rtk->opt;
    amb_t *pamb=&rtk->sol.bias;
    double bl,dr[3],posu[3],posr[3],didxi=0.0,didxj=0.0,*im;
    double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,fi,fj,lami,lamj,*Hi=NULL,*RR=NULL;
    int i,j,k,f,m,ff=0,nv=0,nb[NFREQ*4*2+2]={0},b=0,nf=1;
    int sat1,sat2,pi,pj,ib=-1,sysi,sysj,nd=0;

    bl=arc_baseline(x,rtk->rb,dr);
    if (bl<=0.0) {
        arc_log(ARC_WARNING,"base station and rover station position is error\n");
        return 0;
    }
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);

    Ri=arc_mat(ns*nf*2+2,1); Rj=arc_mat(ns*nf*2+2,1); im=arc_mat(ns,1);
    tropu=arc_mat(ns,1); tropr=arc_mat(ns,1); dtdxu=arc_mat(ns,3); dtdxr=arc_mat(ns,3);

    /* initial satellite status informations */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].resp[0]=rtk->ssat[i].resc[0]=0.0;
    }
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
    /* construct double-difference residuals */
    for (m=0;m<4;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        for (f=opt->mode>PMODE_DGPS?0:nf;f<nf*2;f++) {

            for (i=0;i<pamb->nb;i++) {

                if (!pamb->amb[i].update) continue;
                sat1=pamb->amb[i].sat1; sat2=pamb->amb[i].sat2;

                if (sat1<=0||sat2<=0) continue;

                sysi=rtk->ssat[sat1-1].sys;
                sysj=rtk->ssat[sat2-1].sys;
                if (!arc_test_sys(sysi,m)) continue;
                if (!arc_test_sys(sysj,m)) continue;

                /* observation index */
                for (pi=-1,k=0;k<ns;k++) if (sat[k]==sat1) {pi=k;break;}
                for (pj=-1,k=0;k<ns;k++) if (sat[k]==sat2) {pj=k;break;}
                if (pi<0||pj<0) continue;

                /* check carrier wave lengths */
                ff=f%nf;
                lami=nav->lam[sat1-1][ff];
                lamj=nav->lam[sat2-1][ff];
                if (lami<=0.0||lamj<=0.0||lami!=lamj) continue;
                if (H) {
                    Hi=H+nv*rtk->nx;
                    for (k=0;k<rtk->nx;k++) Hi[k]=0.0;
                }
                /* double-differenced residual */
                if (v) v[nv]=(y[f+iu[pi]*nf*2]-y[f+ir[pi]*nf*2])
                            -(y[f+iu[pj]*nf*2]-y[f+ir[pj]*nf*2]);

                /* partial derivatives by rover position */
                if (H) {
                    for (k=0;k<3;k++) {
                        Hi[k]=-e[k+iu[pi]*3]+e[k+iu[pj]*3];
                    }
                }
                /* double-differenced ionospheric delay term */
                if (opt->ionoopt==IONOOPT_EST) {
                    fi=lami/lam_carr[0]; fj=lamj/lam_carr[0];
                    didxi=(f<nf?-1.0:1.0)*fi*fi*im[pi];
                    didxj=(f<nf?-1.0:1.0)*fj*fj*im[pj];
                    if (v) v[nv]-=didxi*x[II(sat1,opt)]-didxj*x[II(sat1,opt)];

                    /* design matrix */
                    if (H) {
                        Hi[II(sat1,opt)]= didxi;
                        Hi[II(sat2,opt)]=-didxj;
                    }
                }
                /* double-differenced tropospheric delay term */
                if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                    if (v) v[nv]-=(tropu[pi]-tropu[pj])-(tropr[pi]-tropr[pj]);
                    for (k=0;k<(opt->tropopt<TROPOPT_ESTG?1:3);k++) {

                        if (!H) continue;
                        /* design matrix of double-differenced tropospheric delay term */
                        Hi[IT(0,opt)+k]= (dtdxu[k+pi*3]-dtdxu[k+pj*3]);
                        Hi[IT(1,opt)+k]=-(dtdxr[k+pi*3]-dtdxr[k+pj*3]);
                    }
                }
                /* double-differenced phase-bias term */
                if (f<nf) {
                    ib=pamb->amb[i].id+rtk->na;/* index of double-difference ambiguity states */
                    if (v) v[nv]-=lami*x[ib];
                    if (H) Hi[ib]=lami;
                }
                /* test innovation for carrier-phase and pseudorange */
                if (v) {
                    if (f<nf&&opt->maxinnoc>0.0
                        &&fabs(v[nv])>opt->maxinnoc) { /* for carrier-phase */

                        rtk->ssat[sat1-1].rejc[f]++;
                        rtk->ssat[sat2-1].rejc[f]++;

                        /* no updates this ambiguity */
                        pamb->amb[i].update=0;

                        arc_log(ARC_WARNING,"arc_ddres_ddamb: outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                                sat1,sat2,f<nf?"L":"P",f%nf+1,v[nv]);
                        continue;
                    }
                    else if (f==nf&&opt->maxinno>0.0&&
                             fabs(v[nv])>opt->maxinno) { /* for pseudorange */

                        arc_log(ARC_WARNING,"arc_ddres_ddamb: outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                                sat1,sat2,f<nf?"L":"P",f%nf+1,v[nv]);
                        continue;
                    }
                }
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
                    rtk->ceres_active_x[ib]=1;
                    /* hold phase observation numbers */
                    rtk->nc++;
                }
                else rtk->np++;  /* hold code observation numbers */
                if (v) {
                    if (f<nf) rtk->ssat[sat1-1].resc[f   ]=v[nv];
                    else      rtk->ssat[sat2-1].resp[f-nf]=v[nv];
                }
                /* single-differenced measurement error variances */
                if (R) {
                    Rj[nv]=arc_varerr(sat2,sysj,azel[1+iu[pj]*2],bl,dt,f,opt)+
                           arc_snr_varerr(rtk->ssat[sat2-1].snr[0]*0.25,f,nf,opt);
                    Ri[nv]=arc_varerr(sat1,sysi,azel[1+iu[pi]*2],bl,dt,f,opt)+
                           arc_snr_varerr(rtk->ssat[sat1-1].snr[0]*0.25,f,nf,opt);
                }
                /* set valid data flags */
                if (opt->mode>PMODE_DGPS) {
                    if (f<nf) rtk->ssat[sat1-1].vsat[f]=rtk->ssat[sat2-1].vsat[f]=1;
                }
                else {
                    rtk->ssat[sat1-1].vsat[f-nf]=rtk->ssat[sat2-1].vsat[f-nf]=1;
                }
                arc_log(ARC_INFO,"arc_ddres_ddamb: sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n",
                        sat1,sat2,f<nf?"L":"P",f%nf+1,v[nv],
                        R==NULL?-999.0:Ri[nv],R==NULL?-999.0:Rj[nv]);
                if (vflg) vflg[nv]=(sat1<<16)|(sat2<<8)|((f<nf?0:1)<<4)|(f%nf);
                nb[b]++;
                nv++;  /* increase double-residual index */
            }
            b++;
        }
    } /* end of system loop */

    if (H&&nv) {
        arc_log(ARC_INFO,"arc_ddres_ddamb : H=\n");
        arc_tracemat(ARC_MATPRINTF,H,rtk->nx,nv,7,4);
    }
    /* double-differenced measurement error covariance */
    if (R&&nv) {
        arc_ddcov(nb,b,Ri,Rj,nv,R);
    }
    /* doppler partial derivatives */
    if (opt->est_doppler) {

        /* get doppler residuals */
        nd=arc_doppler_res(rtk,x,sat,ns,v,H,Ri,nv,rtk->nx,
                           rtk->opt.dynamics,rtk->opt.est_doppler,IC(opt));
        nv+=nd; /* numbers of residuals */

        /* rerset covariance of active states */
        if (R&&nv) {
            RR=arc_zeros(nv,nv);
            for (i=0;i<nv-nd;i++) {
                for (j=0;j<nv-nd;j++) RR[i+j*nv]=R[i+j*(nv-nd)];
            }
            for (i=0;i<nd;i++) RR[nv-nd+i+(nv-nd+i)*nv]=Ri[i];

            /* save measurements variance matrix */
            arc_matcpy(R,RR,nv,nv); free(RR);

            arc_log(ARC_INFO,"arc_ddres_ddamb : add doppler,R=\n");
            arc_tracemat(ARC_MATPRINTF,R,nv,nv,10,4);
        }
        if (H&&nv) {
            arc_log(ARC_INFO,"arc_ddres_ddamb : add doppler,H=\n");
            arc_tracemat(ARC_MATPRINTF,H,rtk->nx,nv,8,4);
        }
    }
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
    double *Ri,*Rj,lami,lamj,*Hi=NULL,*RR=NULL;
    int i,j,k,m,f,ff=0,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=1,nx=NXDC(&rtk->opt),nd=0;

    arc_log(ARC_INFO,"arc_init_dc_res   : dt=%.1f nx=%d ns=%d\n",dt,nx,ns);

    bl=arc_baseline(x,rtk->rb,dr);
    if (bl<=0.0) {
        arc_log(ARC_WARNING,"base station and rover station position is error\n");
        return 0;
    }
    ecef2pos(x,posu); ecef2pos(rtk->rb,posr);

    Ri=arc_mat(ns*nf*2+2,1); Rj=arc_mat(ns*nf*2+2,1);

    for (m=0;m<NUMOFSYS;m++) /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        for (f=nf;f<nf*2;f++) {

            /* search reference satellite with highest elevation */
            for (i=-1,j=0;j<ns;j++) {
                sysi=rtk->ssat[sat[j]-1].sys;
                if (!arc_test_sys(sysi,m)) continue;

                /* snr outlier */
                if (rtk->ssat[sat[j]-1].snrf[0]) continue;

                if (rtk->ssat[sat[j]-1].dcvl[0]) continue;

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

                if (y) if (!arc_validobs(iu[j],ir[j],f,nf,y)) continue;

                /* if detect snr outlier */
                if (rtk->ssat[sat[i]-1].snrf[0]||rtk->ssat[sat[j]-1].snrf[0]) continue;

                if (rtk->ssat[sat[j]-1].dcvl[0]) continue;

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
                if (v) if (opt->maxinno>0.0
                       &&fabs(v[nv])>opt->maxinno) {
                    arc_log(ARC_WARNING,"arc_init_dc_res : outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                            sat[i],sat[j],f<nf?"L":"P",f%nf+1,v[nv]);
                    continue;
                }
                /* single-differenced measurement error variances */
                if (R) {
                    Rj[nv]=arc_varerr(sat[j],sysj,azel[1+iu[j]*2],bl,dt,f,opt)+
                           arc_snr_varerr(rtk->ssat[sat[j]-1].snr[0]*0.25,f,nf,opt);
                    Ri[nv]=arc_varerr(sat[i],sysi,azel[1+iu[i]*2],bl,dt,f,opt)+
                           arc_snr_varerr(rtk->ssat[sat[i]-1].snr[0]*0.25,f,nf,opt);
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

    /* double-differenced measurement error covariance */
    if (R&&nv) arc_ddcov(nb,b,Ri,Rj,nv,R);
    if (H&&nv) {
        arc_log(ARC_INFO,"arc_init_dc_res : add doppler,H=\n");
        arc_tracemat(ARC_MATPRINTF,H,nx,nv,7,4);
    }
    /* doppler partial derivatives */
    if (opt->est_doppler) {

        nd=arc_doppler_res(rtk,x,sat,ns,v,H,Ri,nv,nx,
                           rtk->opt.dynamics_dc,rtk->opt.est_doppler,ICDC(opt));
        nv+=nd; /* numbers of residuals */

        /* rerset covariance of active states */
        if (R&&nv) {
            RR=arc_zeros(nv,nv);
            for (i=0;i<nv-nd;i++) {
                for (j=0;j<nv-nd;j++) RR[i+j*nv]=R[i+j*(nv-nd)];
            }
            for (i=0;i<nd;i++) RR[nv-nd+i+(nv-nd+i)*nv]=Ri[i];

            /* save measurements variance matrix */
            arc_matcpy(R,RR,nv,nv); free(RR);

            arc_log(ARC_INFO,"arc_init_dc_res : add doppler,R=\n");
            arc_tracemat(ARC_MATPRINTF,R,nv,nv,10,4);
        }
        if (H&&nv) {
            arc_log(ARC_INFO,"arc_init_dc_res : add doppler,H=\n");
            arc_tracemat(ARC_MATPRINTF,H,nx,nv,7,4);
        }
    }
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

    arc_log(ARC_INFO,"arc_intpres : n=%d tt=%.1f\n",n,tt);

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
    double el=-999.0;

    arc_log(ARC_INFO,"arc_ddmat   :\n");

    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].fix[0]=0; /* initial fix flag,this step is very importance */
        rtk->amb_nb=0; /* reset numbers of double difference ambiguity */
    }
    for (i=0;i<NUMOFSYS;i++) rtk->amb_refsat[i]=0; /* initial for reference satellite */

    /* initial single-difference to double-diffrence transformation matrix */
    if (D) for (i=0;i<na;i++) D[i+i*nx]=1.0;

    for (m=0;m<4;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

        nofix=(m==3&&rtk->opt.bdsmodear==0)
              ||(m==0&&rtk->opt.gpsmodear==0); /* this flag is importance */
        for (f=0,k=na;f<nf;f++,k+=MAXSAT) {  /* loop for frequency */

            for (i=k,el=-999.0,ref=-1;i<k+MAXSAT;i++) { /* loop for ambiguity list */

                /* select the reference single-difference ambiguity */
                if (rtk->x[i]==0.0||!rtk->ceres_active_x[i] /* importance */
                    ||!arc_test_sys(rtk->ssat[i-k].sys,m)
                    ||!rtk->ssat[i-k].vsat[f]||!rtk->ssat[i-k].half[f]) {
                    continue;
                }
                if (rtk->ssat[i-k].lock[f]>0&&!(rtk->ssat[i-k].slip[f]&2)&&
                    rtk->ssat[i-k].azel[1]>=el&&!nofix&&rtk->ceres_active_x[i]) {
                    el=rtk->ssat[i-k].azel[1]; ref=i; /* index of reference single-difference ambiguity */
                    continue;
                }
            }
            /* set the reference ambiguity index */
            if (ref<0) continue;
            rtk->ssat[ref-k].fix[f]=2;  /* fix */
            rtk->amb_refsat[m]=ref-k+1; /* sat no. */

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
                    rtk->amb_index[nb]   =j-k+1;    /* double-difference ambiguity index (sat no),no included reference satellite */
                    rtk->ssat[j-k].fix[f]=2;        /* fix this ambiguity*/
                    rtk->ddsat[2*nb  ]=rtk->amb_refsat[m];
                    rtk->ddsat[2*nb+1]=j-k+1;
                    nb++;                           /* numbers of double-difference ambiguity */
                }
                else rtk->ssat[j-k].fix[0]=1;
            }
        }
    }
    rtk->amb_nb=nb; /* save numbers of double-difference ambiguity */
    if (D) {
        arc_log(ARC_INFO,"D=\n"); arc_tracemat(ARC_MATPRINTF,D,nx,na+nb,2,0);
    }
    return nb; /* numbers of double-difference ambiguity */
}
/* adop---------------------------------------------------------------------- */
static double arc_amb_adop(const double *Qb,const double *xb,int nb)
{
    double det=0.0; arc_matdet(Qb,nb,&det); return pow(det,1.0/(2.0*nb));
}
/* compute the bootstrapped success rate ------------------------------------*/
extern double arc_amb_bs_success(const double *D,int n)
{
    double s=1.0;
    int i; for (i=0;i<n;i++) s*=(2.0*arc_norm_distri(0.5/SQRT(D[i]))-1.0);
    return s;
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

    if (rtk->opt.use_dd_sol) return; /* no hold double-difference ambiguity */

    v=arc_mat(nb,1); H=arc_zeros(nb,rtk->nx);

    for (j=1;j<(group==0?2:3);j++) /* two group double-difference ambiguity */

        for (m=0;m<NUMOFSYS;m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */

            /* reference single-difference satellite */
            refsat=group==0?rtk->amb_refsat[m]:
                   rtk->amb_group_refsat[m][j];

            /* no double-difference satellites */
            if (refsat<=0) continue;

            for (n=i=0;i<MAXSAT;i++) {
                if (!arc_test_sys(rtk->ssat[i].sys,m)
                    ||rtk->ssat[i].fix[0]!=2
                    ||rtk->ssat[i].azel[1]<rtk->opt.elmaskhold
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
                nv++; /* numbers of holding ambiguity */
            }
        }
    if (nv>0) {
        R=arc_zeros(nv,nv);
        for (i=0;i<nv;i++) R[i+i*nv]=VAR_HOLDAMB; /* measurement variance matrix */

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
    }
    free(db); free(QQ);
    return info;
}
/* boostraping resolve integer ambiguity--------------------------------------*/
static int arc_resamb_BOOST(rtk_t *rtk,double *bias,double *xa)
{
    prcopt_t *opt=&rtk->opt;
    int i,j,ny,nb,nx=rtk->nx,na=rtk->na;
    double *D,*y,*Qy,*b,*Qb,*Qab,Ps=0.0;

    arc_log(ARC_INFO,"arc_resamb_BOOST: nx=%d\n",nx);

    rtk->sol.ratio=0.0;
    rtk->sol.dop.dops[4]=-999.0;
    rtk->sol.p_ar=(float)-999.0;

    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
        rtk->opt.thresar[0]<1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D=arc_zeros(nx,nx);
    if ((nb=arc_ddmat(rtk,D))<=0) {
        arc_log(ARC_WARNING,"arc_resamb_BOOST: no valid double-difference\n");
        free(D); return 0;
    }
    ny=na+nb;
    y=arc_mat(ny,1); Qy=arc_mat(ny,ny);
    b=arc_mat(nb,2); Qb=arc_mat(nb,nb);Qab=arc_mat(na,nb);

    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    arc_amb_s2d(rtk,D,y,Qy,ny,nx);

    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    arc_amb_Qb(nb,na,ny,Qb,Qab,Qy);

    /* updates adop */
    if (rtk->opt.amb_adop) {
        rtk->sol.dop.dops[4]=arc_amb_adop(Qb,y+na,nb);
        arc_log(ARC_INFO,"ADOP=%8.4lf \n",rtk->sol.dop.dops[4]);
    }
    /* boostrap fix ambiguity */
    if (!arc_bootstrap(nb,y+na,Qb,b,&Ps)) {

        arc_log(ARC_INFO,"N=");
        arc_tracemat(ARC_MATPRINTF,b,1,nb,10,3);

        rtk->sol.p_ar=Ps; /* boostraping success rate */

        if (Ps>=opt->amb_boostps) {

            arc_log(ARC_INFO,"arc_resamb_BOOST: fix ambiguity ok,nb=%3d,Ps=%.3lf\n",nb,Ps);

            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i=0;i<na;i++) {
                rtk->xa[i]=rtk->x[i];
                for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
            }
            for (i=0;i<nb;i++) {
                bias[i]=b[i]; y[na+i]-=b[i];
            }
            /* fix solutions */
            if (!arc_amb_fix_sol(rtk,nb,na,Qb,y,Qab)) {
                /* restore single-differenced ambiguity */
                arc_restamb(rtk,bias,xa,0);
            }
            else nb=0; /* set numbers of double-difference ambiguity to zero */
        }
        else {
            arc_log(ARC_INFO,"arc_resamb_BOOST: fix ambiguity failed,nb=%3d,Ps=%.3lf\n",nb,Ps);
            nb=0; /* set numbers of double-difference ambiguity to zero */
        }
    }
    else nb=0; /* set numbers of double-difference ambiguity to zero */

    if (nb==0) {
        arc_log(ARC_INFO,"arc_resamb_BOOST: fix ambiguity failed,Ps=%.3lf\n",Ps);
    }
    free(D); free(y); free(Qy);
    free(b); free(Qb); free(Qab);
    return nb; /* number of ambiguities */
}
/* save double-difference ambiguity-------------------------------------------*/
static void arc_save_ddamb(rtk_t *rtk,double *bias,int nb,double *xa)
{
    arc_log(ARC_INFO,"arc_save_ddamb: \n");

    int i,j,index=-1;
    amb_t *pamb=&rtk->sol.bias;

    for (i=0;i<rtk->nx;i++) xa[i]=rtk->x [i];  /* ambiguity */
    for (i=0;i<rtk->na;i++) xa[i]=rtk->xa[i];  /* station position/trp/iono/clock and so on */

    for (i=0,j=0;i<pamb->nb;i++) {
        if (!pamb->amb[i].update
            ||rtk->ssat[pamb->amb[i].sat2-1].azel[1]<rtk->opt.elmaskhold) continue;
        pamb->amb[i].pt   =pamb->amb[i].t;      /* precious epoch time */
        pamb->amb[i].t    =rtk->sol.time;       /* current epoch time */
        pamb->amb[i].pratio=pamb->amb[i].ratio; /* precious epoch ratio */
        pamb->amb[i].ratio=rtk->sol.ratio;      /* lambda ratio */
        pamb->amb[i].dv   =bias[j]-pamb->amb[i].b; /* difference of current and precious epoch */
        pamb->amb[i].pb   =pamb->amb[i].b;         /* precious epoch ambiguity value */
        pamb->amb[i].b    =bias[j];                /* current epoch double-difference */
        pamb->amb[i].flag =AMB_FIX;                /* fix flag */

        index=rtk->na+pamb->amb[i].id; /* index of ambiguity in states list */
        xa[index]=pamb->amb[i].b; /* double-differene ambiguity */

        j++; /* numbers of having fixed-ambiguity */

        if      (fabs(pamb->amb[i].dv)<=1E-5) pamb->amb[i].c++; /* fix counts */
        else if (fabs(pamb->amb[i].dv)>=1E-4) pamb->amb[i].c=0; /* reset fix counts */
    }
    arc_log(ARC_INFO,"xa=\n");
    arc_tracemat(ARC_MATPRINTF,xa,1,rtk->nx,10,4);
}
/* if lambda failed,use last time fix-ambiguity to resolve ambiguity----------*/
static int arc_resamb_by_inherit(const rtk_t *rtk,const int* index,int nb,
                                 const double *yb,const double *b,int *ix,double *bias)
{
    arc_log(ARC_INFO,"arc_resamb_by_inherit :\n");

    int i,j,na=rtk->na;
    ddamb_t *pamb=NULL;

    for (i=0;i<nb;i++) {
        pamb=&rtk->sol.bias.amb[index[i+na]-na];
        pamb->update=0; /* initial double-difference ambiguity updates flag */
    }
    for (i=0,j=0;i<nb;i++) { /* select all meet double-difference ambiguity */
        if ((pamb=&rtk->sol.bias.amb[index[i+na]-na])
            &&pamb->c>=rtk->opt.inherit_fixc /* check fix counts */
            &&fabs(timediff(pamb->t,rtk->sol.time))
              <=(fabs(rtk->tt)*rtk->opt.inherit_age+2*DTTOL) /* check inherit time age */
            &&fabs(yb[i]-pamb->b)<=rtk->opt.inherit_thres) { /* inherit thres of ambiguity */
            bias[j]=pamb->b; /* double-difference ambiguity */
            ix[j++]=i; /* index of ambiguity */
            pamb->update=1; /* updates this ambiguity */
        }
    }
    for (i=0;i<nb;i++) {
        pamb=&rtk->sol.bias.amb[index[i+na]-na];
        if (!pamb->update) continue;
        arc_log(ARC_INFO,"%s: %3d-%3d, %6.4lf \n",time_str(pamb->t,2),pamb->sat1,
                pamb->sat2,pamb->b);
    }
    return j;
}
/* difference test for lambda-------------------------------------------------*/
static int arc_lambda_diff_test(const rtk_t *rtk,const double r1,const double r2)
{
    arc_log(ARC_INFO,"arc_lambda_diff_test: \n");
    return fabs(r1-r2)>=rtk->opt.lambda_diff;
}
/* lambda project test--------------------------------------------------------*/
static int arc_lambda_project_test(const rtk_t *rtk,const double *b,const double *y,
                                   const double *Qb,int na,int nb)
{
    arc_log(ARC_INFO,"arc_lmabda_project_test: \n");

    int i;
    double *db=arc_mat(1,nb),*T=arc_mat(nb,nb),*d=arc_mat(1,nb),r1,r2;

    for (i=0;i<nb;i++) db[i]=b[nb+i]-y[na+i];
    arc_matmul("TN",1,nb,nb,1.0,db,Qb,0.0,T);
    arc_matmul("NN",1,1, nb,1.0,T,db,0.0,&r1);

    for (i=0;i<nb;i++) d [i]=b[i   ]-y[na+i];
    for (i=0;i<nb;i++) db[i]=b[nb+i]-y[na+i];

    arc_tracemat(ARC_MATPRINTF,d, 1, nb,10,4);
    arc_tracemat(ARC_MATPRINTF,db,1, nb,10,4);
    arc_tracemat(ARC_MATPRINTF,Qb,nb,nb,10,4);

    arc_matmul("TN",1,nb,nb,1.0,db,Qb,0.0,T);
    arc_matmul("NN",1,1,nb,1.0,T,d,0.0,&r2);

    arc_log(ARC_INFO,"r1=%8.4lf, r2=%8.4lf, r2/r1=%8.4lf \n",r1,r2,r2/r1);

    free(db); free(T); free(d);

    return fabs(r2/r1)<=rtk->opt.lambda_project_thres; /* test lambda project */
}
/* resolve integer ambiguity for direct double-difference ambiguity-----------*/
static int arc_resamb_DirectDD_LAMBDA(rtk_t *rtk,double *bias,double *xa)
{
    arc_log(ARC_INFO,"arc_resamb_DirectDD_LAMBDA :\n");

    amb_t *pamb=&rtk->sol.bias;
    prcopt_t *opt=&rtk->opt;
    int i,j,index[MAXSAT],ny,nb,na=rtk->na,nx=rtk->nx,info=-1,nofix;
    double *y,*Qy,*b,*Qb,*Qab,s[2];
    double varf=0.0,var0=0.0;

    rtk->sol.ratio=0.0; /* initial lambda ratio */
    for (i=0;i<pamb->nb;i++) pamb->amb[i].flag=AMB_FLOAT;

    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
        rtk->opt.thresar[0]<1.0) {
        return 0;
    }
    /* exctract double-difference ambiguity from state list */
    for (i=0,j=0;i<pamb->nb;i++) {

        nofix=((satsys(pamb->amb[i].sat1,NULL)==SYS_GPS
                &&satsys(pamb->amb[i].sat2,NULL)==SYS_GPS&&opt->gpsmodear==0))
              ||((satsys(pamb->amb[i].sat1,NULL)==SYS_CMP
                  &&satsys(pamb->amb[i].sat2,NULL)==SYS_CMP&&opt->bdsmodear==0));

        if (nofix) continue; /* no fix this ambiguity */

        if (pamb->amb[i].update
            &&rtk->ssat[pamb->amb[i].sat2-1].azel[1]>=opt->elmaskhold)
            index[na+(j++)]=pamb->amb[i].id+na; /* fix this ambiguity */
    }
    for (i=0;i<na;i++) index[i]=i; /* index of double-difference ambiguity */

    ny=na+j; nb=j;
    y=arc_mat(1,ny); b=arc_mat(2,nb);
    Qy=arc_mat(ny,ny); Qb=arc_mat(nb,nb); Qab=arc_mat(na,nb);

    for (i=0;i<ny;i++) y[i]=rtk->x[index[i]];
    for (i=0;i<ny;i++) for (j=0;j<ny;j++) Qy [i+j*ny]=rtk->P[index[i   ]+index[j   ]*nx];
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]=rtk->P[index[na+i]+index[na+j]*nx];
    for (i=0;i<na;i++) for (j=0;j<nb;j++) Qab[i+j*na]=rtk->P[i          +index[na+j]*nx];

    arc_log(ARC_INFO,"y=\n"); arc_tracemat(ARC_MATPRINTF,y,1,ny,10,4);

    arc_log(ARC_INFO,"Qy=\n");
    arc_tracemat(ARC_MATPRINTF,Qy,ny,ny,10,4);
    arc_log(ARC_INFO,"Qb=\n");
    arc_tracemat(ARC_MATPRINTF,Qb,nb,nb,10,4);
    arc_log(ARC_INFO,"Qab=\n");
    arc_tracemat(ARC_MATPRINTF,Qab,na,nb,10,4);

    /* updates adop */
    if (rtk->opt.amb_adop) {
        rtk->sol.dop.dops[4]=arc_amb_adop(Qb,y+na,nb);
        arc_log(ARC_INFO,"ADOP=%8.4lf \n",rtk->sol.dop.dops[4]);
    }
    /* lambda/mlambda integer least-square estimation */
    if (!(info=arc_lambda(nb,2,y+na,Qb,b,s,NULL,NULL))) {

        arc_log(ARC_INFO,"N(1)="); arc_tracemat(ARC_MATPRINTF,b   ,1,nb,10,4);
        arc_log(ARC_INFO,"N(2)="); arc_tracemat(ARC_MATPRINTF,b+nb,1,nb,10,4);

        rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;

        /* validation by popular ratio-test */
        if (s[0]<=0.0||s[1]/s[0]>=opt->thresar[0]
            ||arc_lambda_diff_test(rtk,s[1],s[0])
            ||arc_lambda_project_test(rtk,b,y,Qb,na,nb)) {

            arc_log(ARC_INFO,"arc_resamb_DirectDD_LAMBDA: "
                            "validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                    nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);

            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i=0;i<na;i++) {
                rtk->xa[i]=rtk->x[i];
                for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
            }
            for (i=0;i<nb;i++) {
                bias[i]=b[i]; y[na+i]-=b[i];
            }
            /* fix solutions */
            if (!arc_amb_fix_sol(rtk,nb,na,Qb,y,Qab)) {
                /* restore double-difference ambiguity */
                arc_save_ddamb(rtk,bias,nb,xa);
            }
            else nb=0; /* set numbers of double-difference ambiguity to zero */
        }
        else if (opt->amb_inherit) {

            arc_log(ARC_INFO,"double-difference inherit when lambda failed\n");

            int *ix=arc_imat(1,nb),n;
            double *Qb_=arc_mat(nb,nb),*Qab_=arc_mat(na,nb),*b_=arc_mat(1,nb);

            arc_matcpy(b_,b,1,nb);
            if ((n=arc_resamb_by_inherit(rtk,index,nb,y+na,b_,ix,b))) {

                for (i=0;i<n;i++) y[na+i]=y[na+ix[i]];
                for (i=0;i<n;i++) for (j=0;j<n; j++) Qb_ [i+j*n ]=Qb [ix[i]+ix[j]*nb];
                for (i=0;i<n;i++) for (j=0;j<na;j++) Qab_[j+i*na]=Qab[j+ix[i]*na];

                nb=n; /* reset size of double-difference ambiguity */

                arc_log(ARC_INFO,"yb=\n");
                arc_tracemat(ARC_MATPRINTF,y+na,1,nb,10,4);
                arc_log(ARC_INFO,"b=\n");
                arc_tracemat(ARC_MATPRINTF,b,1,nb,10,4);
                arc_log(ARC_INFO,"Qb_=\n");
                arc_tracemat(ARC_MATPRINTF,Qb_,nb,nb,10,4);
                arc_log(ARC_INFO,"Qab_=\n");
                arc_tracemat(ARC_MATPRINTF,Qab_,na,nb,10,4);

                /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
                for (i=0;i<na;i++) {
                    rtk->xa[i]=rtk->x[i];
                    for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
                }
                for (i=0;i<nb;i++) {
                    bias[i]=b[i];   /* fixed solutions */
                    y[na+i]-=b[i]; /* b0-b */
                }
                /* fix solutions */
                if (!arc_amb_fix_sol(rtk,nb,na,Qb_,y,Qab_)) {

                    /* check inherit fix solution */
                    for (i=0;i<na;i++) var0+=rtk->P[i+i*nx],varf+=rtk->Pa[i+i*na];
                    
                    if (varf<var0) {
                        /* set inherit double-difference ambiguity flag */
                        rtk->inherit_fix=AMB_INHERIT_FIX;
                        /* restore double-difference ambiguity */
                        arc_save_ddamb(rtk,bias,nb,xa);
                    }
                    else {
                        arc_log(ARC_WARNING,"arc_resamb_DirectDD_LAMBDA: "
                                "double-difference ambiguity inherit fix is failed \n");
                        nb=0; /* set numbers of double-difference ambiguity to zero */
                    }
                }
                else nb=0;
            }
            else {
                arc_log(ARC_WARNING,"no fixed double difference ambiguity \n");
                nb=0; /* set numbers of double-difference ambiguity to zero */
            }
            free(ix); free(Qb_); free(Qab_);
        }
        else {
            arc_log(ARC_WARNING,"arc_resamb_DirectDD_LAMBDA : ambiguity validation "
                    "failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",nb,s[1]/s[0],s[0],s[1]);
            nb=0; /* set numbers of double-difference ambiguity to zero */
        }
    }
    else {
        arc_log(ARC_WARNING,"lambda error (info=%d)\n",info);
        nb=0; /* set numbers of double-difference ambiguity to zero */
    }
    free(y); free(Qy);
    free(b); free(Qb); free(Qab);

    return nb; /* number of ambiguities */
}
/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
static int arc_resamb_LAMBDA(rtk_t *rtk,double *bias,double *xa)
{
    prcopt_t *opt=&rtk->opt;
    int i,j,ny,nb,nx=rtk->nx,na=rtk->na,k,ok=0,index[MAXSAT];
    double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*QQ,s[2],*yb,*DB=NULL,*DD=NULL;
    double varf=0.0,var0=0.0;
    ddamb_t *pamb=NULL;

    arc_log(ARC_INFO,"arc_resamb_LAMBDA : nx=%d\n",nx);

    rtk->sol.ratio=0.0;
    rtk->sol.dop.dops[4]=-999.0;
    rtk->sol.p_ar=(float)-999.0; /* initial ambiguity solution success rate */
    rtk->inherit_fix=AMB_INHERIT_FLOAT; /* initial ambiguity inherit fix flag */

    if (opt->use_dd_sol) return arc_resamb_DirectDD_LAMBDA(rtk,bias,xa);

    if (rtk->opt.mode<=PMODE_DGPS||rtk->opt.modear==ARMODE_OFF||
        rtk->opt.thresar[0]<1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D=arc_zeros(nx,nx);
    if ((nb=arc_ddmat(rtk,D))<=0) {
        arc_log(ARC_WARNING,"arc_resamb_LAMBDA: no valid double-difference\n");
        free(D);
        return 0;
    }
    ny=na+nb;
    y=arc_mat(ny,1); Qy=arc_mat(ny,ny);DP=arc_mat(ny,nx);
    b=arc_mat(nb,2); db=arc_mat(nb,1); Qb=arc_mat(nb,nb);
    Qab=arc_mat(na,nb); QQ=arc_mat(na,nb);
    DB=arc_mat(1,nb);

    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    arc_matmul("TN",ny,1,nx,1.0, D,rtk->x,0.0,y);
    arc_matmul("TN",ny,nx,nx,1.0,D,rtk->P,0.0,DP);
    arc_matmul("NN",ny,ny,nx,1.0,DP,D,0.0,Qy);

    arc_log(ARC_INFO,"arc_resamb_LAMBDA: Qy=\n");
    arc_tracemat(ARC_MATPRINTF,Qy,na+nb,na+nb,10,4);

    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]=Qy[na+i+(na+j)*ny];
    for (i=0;i<na;i++) for (j=0;j<nb;j++) Qab[i+j*na]=Qy[   i+(na+j)*ny];

    /* add process noice to Qb */
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb[i+j*nb]+=SQR(rtk->opt.fixn);

    arc_log(ARC_INFO,"arc_resamb_LAMBDA: N(0)=\n");
    arc_tracemat(ARC_MATPRINTF,y+na,1,nb,10,3);

    arc_log(ARC_INFO,"arc_resamb_LAMBDA: Qb=\n");
    arc_tracemat(ARC_MATPRINTF,Qb,nb,nb,10,4);
    arc_log(ARC_INFO,"arc_resamb_LAMBDA: Qab=\n");
    arc_tracemat(ARC_MATPRINTF,Qab,na,nb,10,4);

    /* updates adop */
    if (rtk->opt.amb_adop) {
        rtk->sol.dop.dops[4]=arc_amb_adop(Qb,y+na,nb);
        arc_log(ARC_INFO,"ADOP=%8.4lf \n",rtk->sol.dop.dops[4]);
    }

    /* lambda/mlambda integer least-square estimation */
    if (!(arc_lambda(nb,2,y+na,Qb,b,s,DB,NULL))) {

        arc_log(ARC_INFO, "N(1)=");
        arc_tracemat(ARC_MATPRINTF,b,1,nb,10,3);
        arc_log(ARC_INFO, "N(2)=");
        arc_tracemat(ARC_MATPRINTF,b+nb,1,nb,10,3);

        if (DD==NULL) {
            DD=arc_mat(1,nb); arc_matcpy(DD,DB,nb,1);
        }
        /* updates ar sucess probability value */
        rtk->sol.p_ar=arc_amb_bs_success(DB,nb);

        arc_log(ARC_INFO,"AR sucess probability=%8.4f \n",rtk->sol.p_ar);

        rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;

        /* validation by popular ratio-test, difference test and project test */
        ok=((s[0]<=0.0||s[1]/s[0]>=(opt->thresar[0]))&&(s[0]<MAXAMBSQ))
           ||arc_lambda_diff_test(rtk,s[1],s[0])
           ||arc_lambda_project_test(rtk,b,y,Qb,na,nb);

        if (ok) {

            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i=0;i<na;i++) {
                rtk->xa[i]=rtk->x[i];
                for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
            }
            for (i=0;i<nb;i++) {
                bias[i]=b[i];   /* fixed solutions */
                y[na+i]-=b[i]; /* b0-b */
            }
            if (!arc_matinv(Qb,nb)) {
                arc_matmul("NN",nb,1,nb,1.0,Qb,y+na,0.0,db);
                arc_matmul("NN",na,1,nb,-1.0,Qab,db,1.0,rtk->xa);

                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                arc_matmul("NN",na,nb,nb,1.0,Qab,Qb,0.0,QQ);
                arc_matmul("NT",na,na,nb,-1.0,QQ,Qab,1.0,rtk->Pa);

                arc_log(ARC_INFO,"arc_resamb_LAMBDA: validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                        nb,s[0]==0.0?0.0:s[1]/s[0],s[0],s[1]);

                /* restore single-differenced ambiguity */
                arc_restamb(rtk,bias,xa,0); /* bias is the fixed solutions */
            }
            else nb=0; /* set numbers of double-difference ambiguity to zero */
        }
        else if (opt->amb_inherit) {
            arc_log(ARC_INFO,"double-difference inherit when lambda failed\n");
            int *ixb=arc_imat(1,nb);
            double *Qb_=arc_zeros(nb,nb),*Qab_=arc_zeros(na,nb);
            for (i=0,j=0;i<nb;i++) {

                rtk->ssat[rtk->amb_index[i]-1].fix[0]=0; /* reset fix flag */

                if ((pamb=arc_get_ddamb(&rtk->bias,
                                        rtk->ddsat[2*i],rtk->ddsat[2*i+1]))) {

                    if (fabs(timediff(pamb->t,rtk->sol.time))
                        <=(fabs(rtk->tt)*opt->inherit_age+2*DTTOL) /* check age of inherit ambiguity */
                        &&pamb->c>=opt->inherit_fixc /* check fix counts */
                        &&fabs(pamb->b-y[na+i])<=rtk->opt.inherit_thres) { /* check thres of inherit ambiguity */
                        b[j]=pamb->b; /* set double-difference ambiguity */
                        ixb[j++]=i; /* index of double-difference ambiguity */
                        rtk->ssat[rtk->amb_index[i]-1].fix[0]=2; /* fix this ambiguity */
                    }
                }
            }
            /* re-extract double-difference ambiguity covariance matrix */
            for (i=0;i<j;i++) for (k=0;k< j;k++) Qb_ [i+k*j ]=Qb [ixb[i]+ixb[k]*nb];
            for (i=0;i<j;i++) for (k=0;k<na;k++) Qab_[k+i*na]=Qab[k+ixb[i]*na];

            /* adjust for covanriance of Qab */
            for (i=0;i<j;i++) for (k=0;k<na;k++) Qab_[k+i*na]*=(1.0+opt->inherit_IF);

            for (i=0;i<j;i++) y[na+i]=y[na+ixb[i]];

            /* new size of double-difference ambiguity */
            nb=j;

            if (nb>0) {
                arc_log(ARC_INFO,"yb=\n");
                arc_tracemat(ARC_MATPRINTF,y+na,1,nb,10,4);

                arc_log(ARC_INFO,"b=\n");
                arc_tracemat(ARC_MATPRINTF,b,1,nb,10,4);

                arc_log(ARC_INFO,"Qb_=\n");
                arc_tracemat(ARC_MATPRINTF,Qb_,nb,nb,10,4);

                arc_log(ARC_INFO,"Qab_=\n");
                arc_tracemat(ARC_MATPRINTF,Qab_,na,nb,10,4);

                /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
                for (i=0;i<na;i++) {
                    rtk->xa[i]=rtk->x[i];
                    for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];
                }
                for (i=0;i<nb;i++) {
                    bias[i]=b[i];  /* fixed solutions */
                    y[na+i]-=b[i]; /* b0-b */
                }
                arc_log(ARC_INFO,"xf=\n"); arc_tracemat(ARC_MATPRINTF,rtk->x,1,na,10,4);

                /* fix solutions */
                if (!arc_amb_fix_sol(rtk,nb,na,Qb_,y,Qab_)) {

                    arc_log(ARC_INFO,"xa=\n"); arc_tracemat(ARC_MATPRINTF,rtk->xa,1, na,10,4);
                    arc_log(ARC_INFO,"Pa=\n"); arc_tracemat(ARC_MATPRINTF,rtk->Pa,na,na,10,4);

                    /* check inherit fix solution */
                    for (i=0;i<na;i++) var0+=rtk->P[i+i*nx],varf+=rtk->Pa[i+i*na];

                    if (varf<var0) { /* todo:may be have better way to check that */
                        /* restore single-differenced ambiguity */
                        arc_restamb(rtk,bias,xa,0); /* bias is the fixed solutions */
                        /* save solution status */
                        rtk->inherit_fix=AMB_INHERIT_FIX;
                    }
                    else {
                        arc_log(ARC_WARNING,"arc_resamb_LAMBDA: "
                                "double-difference ambiguity inherit fix is failed \n");
                        nb=0; /* set numbers of double-difference ambiguity to zero */
                    }
                }
                else nb=0; /* set numbers of double-difference ambiguity to zero */
            }
            free(ixb); free(Qb_); free(Qab_);
        }
        else nb=0; /* set numbers of double-difference ambiguity to zero */
    }
    else {
        arc_log(ARC_WARNING,"lambda error \n");
        nb=0; /* set numbers of double-difference ambiguity to zero */
    }
    if (rtk->sol.ratio<opt->thresar[0]) {
        arc_log(ARC_WARNING,"arc_resamb_LAMBDA : ambiguity validation "
                "failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",nb,s[1]/s[0],s[0],s[1]);
    }
    if (DD) free(DD);
    free(D); free(y); free(Qy); free(DP); free(DB);
    free(b); free(db); free(Qb); free(Qab); free(QQ);
    return nb; /* number of ambiguities */
}
/* resolve integer ambiguity by group-LAMBDA ---------------------------------*/
static int arc_resamb_group_LAMBDA(rtk_t *rtk,double *bias,double *xa)
{
    return arc_resamb_LAMBDA(rtk,bias,xa);
}
/* partial resolve integer ambiguity------------------------------------------*/
static int arc_resamb_PART(rtk_t *rtk,double *bias,double *xa)
{
    return arc_resamb_LAMBDA(rtk,bias,xa);
}
/* ff-ratio resolve integer ambiguity-----------------------------------------*/
static int arc_resamb_FFRATIO(rtk_t *rtk,double *bias,double *xa)
{
    return arc_resamb_LAMBDA(rtk,bias,xa);;
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

        sat1=(vflg[i]>>16)&0xFF;
        sat2=(vflg[i]>> 8)&0xFF;
        type=(vflg[i]>> 4)&0xF;
        freq=vflg[i]&0xF;
        strcpy(stype,type==0?"L":(type==1?"L":"C"));

        if (v[i]*v[i]<=fact*R[i+i*nv]) {
            rtk->ssat[sat2-1].large_resc[0]=
                    rtk->ssat[sat2-1].large_resc[0]<=0?0:
                    rtk->ssat[sat2-1].large_resc[0]--;
            continue;
        }
        rtk->ssat[sat2-1].large_resc[0]=
                rtk->ssat[sat2-1].large_resc[0]>=MAXLARGERESC?MAXLARGERESC:
                rtk->ssat[sat2-1].large_resc[0]++; /* larger double-difference counter */

        if (strcmp(stype,"C")==0) rtk->ssat[sat2-1].dcvl[0]=1;

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
/* N matrix for adaptive kalman filter-----------------------------------------*/
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
/* states updates of difference pseudorange positioning-----------------------*/
static void arc_diff_pr_update(const rtk_t* rtk,double *xp,double *Pp,double tt)
{
    int i,j,nx=NXDC(&rtk->opt),icdc=ICDC(&rtk->opt);
    double var=0.0;

    arc_log(ARC_INFO,"arc_diff_pr_update   : tt=%.3f\n",tt);

    /* fixed mode */
    if (rtk->opt.mode==PMODE_FIXED) {
        for (i=0;i<3;i++) { /* no update other states */
            arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],1E-6,i,nx);
        }
        return;
    }
    /* initialize position for first epoch */
    if (arc_norm(rtk->xd,3)<=0.0) {
        for (i=0;i<3;i++) {
            arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_POS,i,nx);
        }
        /* initial rover station velecity */
        if (rtk->opt.dynamics_dc) for (i=3;i<nx;i++) {
                arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_VEL,i,nx);
            }
        /* initial rover station clock drift */
        if (rtk->opt.est_doppler) arc_diff_pr_initx(xp,Pp,rtk->sol.clk_dri,
                                                    VAR_CLKDRI,ICDC(&rtk->opt),nx);
    }
    /* add some system noice to states,here only velecity */
    if (rtk->opt.dynamics_dc) {
        for (i=3;i<6;i++) {
            if (xp[i]==0.0) xp[i]=rtk->opt.vel_prn[i-3];
        }
    }
    /* add rover station clock drift process noice */
    if (rtk->opt.est_doppler) {
        if (xp[icdc]==0.0) xp[icdc]=rtk->opt.clk_dri_prn;
    }
    /* static mode */
    if (rtk->opt.mode==PMODE_STATIC) return;

    /* todo:using standard positioning to initial ukf prior states and its covariacne matrix,
     * todo:and may have more better methids to do this */

    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics_dc) {
        if (rtk->opt.init_pnt) {
            for (i=0;i<3;i++) xp[i]=rtk->sol.rr[i];
            if (rtk->opt.est_doppler) {
                arc_diff_pr_initx(xp,Pp,rtk->sol.clk_dri,VAR_CLKDRI,icdc,nx);
            }
            for (i=0;i<3;i++) Pp[i+i*nx]=rtk->sol.qr[i];
            Pp[1+nx*0]=Pp[0+nx*1]=rtk->sol.qr[3];
            Pp[2+nx*0]=Pp[0+nx*2]=rtk->sol.qr[5];
            Pp[1+nx*2]=Pp[2+nx*1]=rtk->sol.qr[4];
            return;
        }
        else {
            if (rtk->opt.est_doppler) {
                arc_diff_pr_initx(xp,Pp,rtk->sol.clk_dri,VAR_CLKDRI,icdc,nx);
            }
            if (rtk->sol.ratio>=rtk->opt.thresar[0]) {
                for (i=0;i<3;i++) arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_POS_AMB,i,nx);
                return;
            }
            for (i=0;i<3;i++) arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_POS,i,nx);
            return;
        }
    }
    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=Pp[i+i*nx]; var/=3.0;

    if (var>VAR_POS||var<0.0) {
        /* reset position with large variance */
        for (i=0;i<3;i++) arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_POS,i,nx);
        for (i=3;i<6;i++) arc_diff_pr_initx(xp,Pp,rtk->sol.rr[i],VAR_VEL,i,nx);
        if (rtk->opt.est_doppler) {
            arc_diff_pr_initx(xp,Pp,rtk->sol.clk_dri,VAR_CLKDRI,icdc,nx);
        }
        arc_log(ARC_INFO,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
    /* state transition of position/velocity/acceleration for dynamic mode */
    double *F=arc_eye(nx),*FP=arc_mat(nx,nx),
            *xpp=arc_mat(nx,1),pos[3]={0},Q[9]={0},Qv[9]={0};

    /* compute the state transition matrix */
    for (i=0;i<3;i++) F[i+(i+3)*nx]=tt;

    arc_log(ARC_INFO,"arc_diff_pr_update : state transition matrix=\n");
    arc_tracemat(ARC_MATPRINTF,F,nx,nx,10,4);

    arc_log(ARC_INFO,"arc_diff_pr_update : before states transition x=\n");
    arc_tracemat(ARC_MATPRINTF,xp,nx,1,10,4);

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
    /* process noice of clock drift */
    if (rtk->opt.est_doppler) Pp[icdc+icdc*nx]
                                      +=SQR(rtk->opt.clk_dri_prn);
    free(F); free(FP); free(xpp);
}
/* difference pseudorange positioning-----------------------------------------*/
static int arc_diff_pr_relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                              const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,*v,*H,*R,*xp,*Pp,*bias,dt;
    int i,j,n=nu+nr,ns,ny,nv,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],niter,nx=NXDC(&rtk->opt);
    int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
    int nf=1,stat=PMODE_DGPS;

    arc_log(ARC_INFO,"arc_diff_pr_relpos  : nx=%d nu=%d nr=%d\n",nx,nu,nr);

    dt=timediff(time,obs[nu].time);

    rs=arc_mat(6,n); dts=arc_mat(2,n);
    var=arc_mat(1,n); y=arc_mat(nf*2,n); e=arc_mat(3,n);
    azel=arc_zeros(2,n);

    /* initial satellite valid flag and snr */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        rtk->ssat[i].vsat[0]=0;
        rtk->ssat[i].snrf[0]=0;
        rtk->ssat[i].dcvl[0]=0; /* valid flag of dc-solution */
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
        stat=SOLQ_NONE;
        return 0;
    }
    /* time-interpolation of residuals (for post-processing) */
    if (opt->intpref) {
        dt=arc_intpres(time,obs+nu,nr,nav,rtk,y+nu*nf*2);
        arc_log(ARC_INFO,"arc_diff_pr_relpos： "
                "time-interpolation of residuals (for post-processing)\n");
    }
    /* select common satellites between rover and base-station */
    if ((ns=arc_selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
        arc_log(ARC_WARNING,"arc_diff_pr_relpos : no common satellite\n");
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        stat=SOLQ_NONE;
        return 0;
    }
    /* initial states for kalman filter */
    xp=arc_mat(nx,1); Pp=arc_zeros(nx,nx);
    arc_matcpy(xp,rtk->xd,nx,1);
    for (i=0;i<nx;i++) for (j=0;j<nx;j++) Pp[i+j*nx]=rtk->Pd[i+nx*j];

    /* updates states */
    arc_diff_pr_update(rtk,xp,Pp,rtk->tt);

    arc_log(ARC_INFO,"arc_diff_pr_relpos: xp=\n");
    arc_tracemat(ARC_MATPRINTF,xp,1,nx,13,4);

    arc_log(ARC_INFO,"arc_diff_pr_relpos : Pp=\n");
    arc_tracemat(ARC_MATPRINTF,Pp,nx,nx,10,4);

    ny=ns*nf*2+2+ns;
    v=arc_mat(ny,1); H=arc_zeros(nx,ny);
    R=arc_zeros(ny,ny); bias=arc_mat(nx,1);

    /* add 2 iterations for baseline-constraint moving-base */
    niter=opt->niter+(opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0?2:0);

    for (i=0;i<niter;i++) {  /* iterations compute */

        /* undifferenced residuals for rover */
        if (!arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk,NULL)) {
            arc_log(ARC_WARNING,"arc_diff_pr_relpos : rover initial position error\n");
            stat=SOLQ_NONE;
            break;
        }
        /* double-differenced residuals and partial derivatives */
        if ((nv=arc_diff_pr_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg))<1) {
            arc_log(ARC_WARNING,"arc_diff_pr_relpos : no double-differenced residual\n");
            stat=SOLQ_NONE;
            break;
        }
        arc_log(ARC_INFO,"arc_diff_pr_relpos ：double-differenced residual vector : \n");
        arc_tracemat(ARC_MATPRINTF,v,nv,1,10,4);

        /* adaptive kaman filter */
        if (opt->adapt_filter) {
            if (!adap_kaman_filter(rtk,xp,Pp,H,v,R,nx,nv)) {
                arc_log(ARC_WARNING,"arc_diff_pr_relpos : adaptive filter error (info=%d)\n",info);
                stat=SOLQ_NONE;
                break;
            }
        }
        /* kalman filter measurement update */
        else {
            if ((info=arc_filter(xp,Pp,H,v,R,nx,nv,NULL))) {
                arc_log(ARC_WARNING,"arc_diff_pr_relpos : filter error (info=%d)\n",info);
                stat=SOLQ_NONE;
                break;
            }
            arc_log(ARC_INFO,"arc_diff_pr_relpos : x(%d)=",i+1);
            arc_tracemat(ARC_MATPRINTF,xp,nx,1,10,4);

            arc_log(ARC_INFO,"arc_diff_pr_relpos : P(%d)=",i+1);
            arc_tracemat(ARC_MATPRINTF,Pp,nx,nx,10,4);
        }
    }
    /* post-fit residuals for float solution */
    if (stat!=SOLQ_NONE&&arc_zdres(0,obs,nu,rs,dts,svh,
                                   nav,xp,opt,0,y,e,azel,rtk,NULL)) {

        /* double-differecen residuals */
        nv=arc_diff_pr_ddres(rtk,nav,dt,xp,Pp,sat,
                             y,e,azel,iu,ir,ns,v,NULL,R,vflg);

        /* validation of float solution */
        if (arc_valpos(rtk,v,R,vflg,nv,ARC_SOLVALTHRES)) {
            /* save solution status */
            arc_matcpy(rtk->sol.prsol.rr,xp,3,1); /* position/velecity */

            for (j=0;j<3;j++) rtk->sol.prsol.qr[j]=Pp[j+j*nx];

            rtk->sol.prsol.qr[3]=Pp[1];    /* cov xy */
            rtk->sol.prsol.qr[4]=Pp[2+nx]; /* cov yz */
            rtk->sol.prsol.qr[5]=Pp[2];    /* cov zx */

            arc_matcpy(rtk->xd,xp,nx,1);
            arc_matcpy(rtk->Pd,Pp,nx,nx);
        }
    }
    free(rs); free(dts); free(var); free(y); free(e); free(azel);
    free(xp); free(Pp);  free(v); free(H); free(R); free(bias);
    return stat;
}
/* another version for extract double-difference ambiguity--------------------*/
static int arc_extract_ddamb_1(rtk_t *rtk,const double *bias)
{
    arc_log(ARC_INFO,"arc_extract_ddamb_1 :\n");

    int i,j;
    amb_t *pamb=&rtk->sol.bias;
    ddamb_t *amb=NULL;

    for (i=0,j=0;i<pamb->nb;i++) {

        if (!pamb->amb[i].update) continue;

        /* get double-difference ambiguity */
        if ((amb=arc_get_ddamb(&rtk->bias,
                               pamb->amb[i].sat1,
                               pamb->amb[i].sat2))==NULL) {
            arc_add_ddamb(&rtk->bias);
            amb=&rtk->bias.amb[rtk->bias.nb++]; /* new double-difference ambiguity */
        }
        *amb=pamb->amb[i]; /* updates double-difference ambiguity */
        j++; /* numbers of double-difference ambiguity */
    }
    pamb=&rtk->bias;
    for (i=0;i<pamb->nb;i++) {
        arc_log(ARC_INFO,"sat=%3d - %3d : %s  %8.3lf  %3d  %8.3lf  %8.3lf \n",
                pamb->amb[i].sat1,pamb->amb[i].sat2,
                time_str(pamb->amb[i].t,2),pamb->amb[i].b,pamb->amb[i].c,
                pamb->amb[i].ratio,pamb->amb[i].dv);
    }
    return j;
}
/* extract double-difference ambiguity----------------------------------------*/
static int arc_extract_ddamb(rtk_t *rtk,const double *bias)
{
    arc_log(ARC_INFO,"arc_extract_ddamb: \n");

    int i,j,k,sys,sat1,sat2,*ix;
    ddamb_t *amb=NULL;
    double p=0.0;

    if (rtk->opt.use_dd_sol) return arc_extract_ddamb_1(rtk,bias);

    ix=arc_imat(rtk->amb_nb,1);

    /* get double-difference ambiguity list */
    for (j=0,i=0;i<rtk->amb_nb;i++) {
        if (rtk->ssat[rtk->amb_index[i]-1].fix[0]==2) ix[j++]=rtk->amb_index[i];
    }
    for (k=0,i=0;i<j;i++) {
        sat2=ix[i]; sys=satsys(sat2,NULL);
        sat1=(sys==SYS_GPS?rtk->amb_refsat[0]
                          :sys==SYS_CMP?rtk->amb_refsat[3]:-1);

        if ((amb=arc_get_ddamb(&rtk->bias,sat1,sat2))==NULL) {
            arc_add_ddamb(&rtk->bias);
            amb=&rtk->bias.amb[rtk->bias.nb++]; /* new double-difference ambiguity */
        }
        /* set observation time and numbers of ambiguity */
        amb->pt=amb->t; amb->t=rtk->sol.time;
        p=amb->b;

        /* updates double-difference ambiguity */
        amb->pb=amb->b; amb->b=bias[i]; /* double-difference ambiguity */
        amb->sat1=sat1; amb->sat2=sat2; /* double-difference satellite */
        amb->pratio=amb->ratio; /* precious epoch lambda ratio */
        amb->ratio =rtk->sol.ratio; /* current epoch lambda ratio */
        amb->dv=amb->b-p; /* delta value of precious and current epoch */

        if      (fabs(amb->dv)<=1E-5) amb->c++; /* fix counts */
        else if (fabs(amb->dv)>=1E-4) amb->c=0; /* reset fix counts */
    }
    for (i=0;i<rtk->bias.nb;i++) {
        arc_log(ARC_INFO,"sat=%3d - %3d : %s  %8.3lf  %3d  %8.3lf  %8.3lf \n",
                rtk->bias.amb[i].sat1,rtk->bias.amb[i].sat2,
                time_str(rtk->bias.amb[i].t,2),rtk->bias.amb[i].b,rtk->bias.amb[i].c,
                rtk->bias.amb[i].ratio,rtk->bias.amb[i].dv);
    }
    free(ix);
    return j; /* numbers of double-difference ambiguity */
}
/* ceres active to filter active index----------------------------------------*/
static int arc_filter_index(const rtk_t *rtk,int *index)
{
    arc_log(ARC_INFO,"arc_filter_index :\n");

    int i,j; for (i=0,j=0;i<rtk->nx;i++) if (rtk->ceres_active_x[i]) index[j++]=i;

    arc_log(ARC_INFO,"filter index=\n");
    arc_tracemati(ARC_MATPRINTF,index,1,j,4,2);

    return j;
}
/* relative positioning ------------------------------------------------------*/
static int arc_relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                      const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,*v,*H,*R,*xp,*Pp,*xa,*bias,dt,*D=NULL,*Qv=NULL,*ve=NULL;
    int i,j,f,n=nu+nr,ns,ny,nv,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT],niter;
    int info,vflg[MAXOBS*NFREQ*2+1],svh[MAXOBS*2];
    int stat=rtk->opt.mode<=PMODE_DGPS?SOLQ_DGPS:SOLQ_FLOAT;
    int nf=1,nk=0,pnx=rtk->nx,index[MAXSAT],ni=0;

    arc_log(ARC_INFO,"arc_relpos  : nx=%d nu=%d nr=%d\n",rtk->nx,nu,nr);

    dt=timediff(time,obs[nu].time);

    rs=arc_mat(6,n); dts=arc_mat(2,n);
    var=arc_mat(1,n); y=arc_mat(nf*2,n); e=arc_mat(3,n);
    azel=arc_zeros(2,n);

    /* initial satellite valid flag and snr */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        rtk->ssat[i].vsat[0]=0;
        rtk->ssat[i].snrf[0]=0;
    }
    if (opt->use_dd_sol) for (i=0;i<rtk->sol.bias.nb;i++) {
        rtk->sol.bias.amb[i].update=0; /* initial ambiguity updates flag */
    }
    /* reset ceres problem solver active states index list */
    for (i=0;i<rtk->nx;i++) rtk->ceres_active_x[i]=0;

    /* satellite positions/clocks */
    arc_satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    /* exclude measurements of eclipsing satellite (block IIA) */
    if (rtk->opt.posopt[3]) {
        arc_testeclipse(obs,n,nav,rs);
    }
    /* undifferenced residuals for base station */
    if (!arc_zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,svh+nu,nav,rtk->rb,opt,1,
                   y+nu*nf*2,e+nu*3,azel+nu*2,rtk,NULL)) {
        arc_log(ARC_WARNING,"arc_relpos : initial base station position error\n");
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    /* time-interpolation of residuals (for post-processing) */
    if (opt->intpref) {
        dt=arc_intpres(time,obs+nu,nr,nav,rtk,y+nu*nf*2);
    }
    /* select common satellites between rover and base-station */
    if ((ns=arc_selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
        arc_log(ARC_WARNING,"arc_relpos : no common satellite\n");
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    /* get the middle point of satellites-elvations */
    if (opt->half_fix) {
        if (opt->use_dd_sol) { /* double-difference ambiguity estimate */
            opt->elmaskhold=arc_amb_adjust_el(rtk,sat,ns);
        }
        else { /* single-difference ambiguity estamate */
            opt->elmaskar=arc_amb_adjust_el(rtk,sat,ns);
        }
    }
    /* temporal update of states */
    arc_udstate(rtk,obs,sat,iu,ir,ns,nav,rs,y,azel);

    xp=arc_mat(rtk->nx,1); Pp=arc_zeros(rtk->nx,rtk->nx);
    xa=arc_mat(rtk->nx,1);
    arc_matcpy(xp,rtk->x,rtk->nx,1);

    ny=ns*nf*2+2+ns;
    v=arc_mat(ny,1); H=arc_zeros(rtk->nx,ny);
    R=arc_zeros(ny,ny); bias=arc_zeros(rtk->nx,1);

    /* add 2 iterations for baseline-constraint moving-base */
    niter=opt->niter+(opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0?2:0);

    if (opt->ceres==0) {

        if (opt->ukf==0) {

            for (i=0;i<niter;i++) {  /* iterations compute */

                /* undifferenced residuals for rover */
                if (!arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk,NULL)) {
                    arc_log(ARC_WARNING,"arc_relpos : rover initial position error\n");
                    stat=SOLQ_NONE;
                    break;
                }
                /* direct double-difference ambiguity solutions */
                if (opt->use_dd_sol) {

                    /* updates double-difference ambiguity */
                    arc_update_ddamb(rtk,obs,sat,iu,ir,ns,nav,rs,y,azel);

                    if (rtk->nx!=pnx&&i==0) { /* states size changed */
                        free(xp); xp=arc_mat(1,      rtk->nx); /* new size */
                        free(Pp); Pp=arc_mat(rtk->nx,rtk->nx); /* new size */
                        free(xa); xa=arc_mat(1,      rtk->nx); /* new size */
                    }
                    arc_matcpy(xp,rtk->x,rtk->nx,1);

                    arc_log(ARC_INFO,"xp=\n");
                    arc_tracemat(ARC_MATPRINTF,xp,1,rtk->nx,10,4);

                    /* double-differenced residuals and partial derivatives */
                    if ((nv=arc_ddres_ddamb(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg))<1) {
                        arc_log(ARC_WARNING,"arc_relpos : no double-differenced residual\n");
                        stat=SOLQ_NONE;
                        break;
                    }
                }
                /* estimate single-difference ambiguity(rtklib mode) */
                else {
                    /* double-differenced residuals and partial derivatives */
                    if ((nv=arc_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,H,R,vflg))<1) {
                        arc_log(ARC_WARNING,"arc_relpos : no double-differenced residual\n");
                        stat=SOLQ_NONE;
                        break;
                    }
                }
                arc_log(ARC_INFO,"arc_relpos ： double-differenced residual vector : \n");
                arc_tracemat(ARC_MATPRINTF,v,nv,1,10,4);

                /* robust kalman filter */
                if (opt->kalman_robust) {
                    D=arc_eye(nv); Qv=arc_mat(nv,nv); ve=arc_mat(nv,1);
                    arc_kalman_norm_inno(rtk,v,nv,H,R,ve); arc_kalman_norm_Qino(ve,nv,Qv);
                    arc_kalman_robust_phi(rtk,ve,nv,nk,D,Qv);
                }
                arc_matcpy(Pp,rtk->P,rtk->nx,rtk->nx);

                arc_log(ARC_INFO,"Pp=\n");
                arc_tracemat(ARC_MATPRINTF,Pp,rtk->nx,rtk->nx,10,4);

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
                    if (opt->use_dd_sol) {
                        ni=arc_filter_index(rtk,index);
                        if ((info=arc_filter_active(xp,Pp,H,v,R,rtk->nx,nv,
                                                    opt->kalman_robust?D:NULL,index,ni))) {
                            arc_log(ARC_WARNING,"arc_relpos : filter error (info=%d)\n",info);
                            stat=SOLQ_NONE;
                            break;
                        }
                    }
                    else {
                        if ((info=arc_filter(xp,Pp,H,v,R,rtk->nx,nv,
                                             opt->kalman_robust?D:NULL))) {
                            arc_log(ARC_WARNING,"arc_relpos : filter error (info=%d)\n",info);
                            stat=SOLQ_NONE;
                            break;
                        }
                    }
                    arc_log(ARC_INFO,"arc_relpos : x(%d)=",i+1);
                    arc_tracemat(ARC_MATPRINTF,xp,1,rtk->nx,10,4);

                    arc_log(ARC_INFO,"arc_relpos : P(%d)=",i+1);
                    arc_tracemat(ARC_MATPRINTF,Pp,rtk->nx,rtk->nx,10,4);
                }
            }
        }
        else if (opt->ukf) {
            fprintf(stderr,"ukf filter haven't been open-source\n");
        }
    }
    if (stat!=SOLQ_NONE&&arc_zdres(0,obs,nu,rs,dts,svh,nav,xp,opt,0,y,e,azel,rtk,NULL)) {

        /* post-fit residuals for float solution */
        if (opt->use_dd_sol) { /* double-difference ambiguity solutions */
            nv=arc_ddres_ddamb(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,NULL,NULL,vflg);
        }
            /* single-difference ambiguity solutions */
        else nv=arc_ddres(rtk,nav,dt,xp,Pp,sat,y,e,azel,iu,ir,ns,v,NULL,NULL,vflg);

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
                         :opt->amb_fix_mode==AMBFIX_LAMBDA?arc_resamb_LAMBDA(rtk,bias,xa)
                         :opt->amb_fix_mode==AMBFIX_BOOTS?arc_resamb_BOOST(rtk,bias,xa)
                         :opt->amb_fix_mode==AMBFIX_PART?arc_resamb_PART(rtk,bias,xa)
                         :arc_resamb_LAMBDA(rtk,bias,xa))) {

        /* extract double-difference ambiguity */
        if (!arc_extract_ddamb(rtk,bias)) {
            arc_log(ARC_WARNING,"arc_relpos: extract double-difference ambiguity,but no ambiguitys \n");
        }
        /* resolve integer ambiguity by LAMBDA success */
        if (arc_zdres(0,obs,nu,rs,dts,svh,nav,xa,opt,0,y,e,azel,rtk,NULL)) {

            /* post-fit reisiduals for fixed solution */
            if (opt->use_dd_sol) {
                nv=arc_ddres_ddamb(rtk,nav,dt,xa,NULL,sat,y,e,azel,iu,ir,ns,v,NULL,NULL,vflg);
            }
            else nv=arc_ddres(rtk,nav,dt,xa,NULL,sat,y,e,azel,iu,ir,ns,v,NULL,NULL,vflg);

            /* validation of fixed solution */
            if (arc_valpos(rtk,v,R,vflg,nv,ARC_SOLVALTHRES)) {

                /* hold integer ambiguity */
                if (++rtk->nfix>=rtk->opt.minfix&&
                    (rtk->opt.modear==ARMODE_FIXHOLD)) {
                    arc_holdamb(rtk,xa,opt->amb_group);
                    rtk->nfix=0; /* reset fix counts */
                }
                else if (opt->inst_amb) { /* hold all single-difference ambiguity */
                    arc_matcpy(rtk->x+rtk->na,
                               xa+rtk->na,rtk->nx-rtk->na,1);
                }
                stat=SOLQ_FIX; /* fix flag */
            }
        }
    }
    else if (stat!=SOLQ_NONE&&opt->no_amb_sol&&!opt->use_dd_sol) { /* lambda is failed */

        /* no-ambiguity-double-difference phase residuals */
        if ((nv=arc_ddres_noamb(rtk,nav,dt,xp,sat,y,
                                e,azel,ns,v,H,R,vflg,ir,iu))) {

            /* kalman filter measurement update */
            if ((info=arc_filter(xp,Pp,H,v,R,rtk->nx,nv,NULL))) {
                arc_log(ARC_WARNING,
                        "no-ambiguity-double-difference : filter error (info=%d)\n",info);
                stat=SOLQ_FLOAT;
            }
            stat=SOLQ_HALFFIX;

            arc_log(ARC_INFO,"no-ambiguity-double-difference : x(%d)=",i+1);
            arc_tracemat(ARC_MATPRINTF,xp,NP(opt),1,10,4);

            arc_log(ARC_INFO,"no-ambiguity-double-difference : P(%d)=",i+1);
            arc_tracemat(ARC_MATPRINTF,Pp,rtk->nx,rtk->nx,10,4);
        }
        else {
            arc_log(ARC_WARNING,"resolve rover station position by no-ambiguity-double-difference"
                    " when lambda is fialed \n");
            stat=SOLQ_FLOAT;
        }
        if (stat==SOLQ_HALFFIX) {
            /* update state and covariance matrix */
            arc_matcpy(rtk->x,xp,rtk->nx,1);
            arc_matcpy(rtk->P,Pp,rtk->nx,rtk->nx);
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
        rtk->fixc++;
    }
    else if (stat==SOLQ_FLOAT) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->x[i];
            rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)rtk->P[1];
        rtk->sol.qr[4]=(float)rtk->P[1+2*rtk->nx];
        rtk->sol.qr[5]=(float)rtk->P[2];
        rtk->halffix++;
    }
    else if (stat==SOLQ_HALFFIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=xp[i];
            rtk->sol.qr[i]=Pp[i+i*rtk->nx];
        }
        rtk->sol.qr[3]=(float)Pp[1];
        rtk->sol.qr[4]=(float)Pp[1+2*rtk->nx];
        rtk->sol.qr[5]=(float)Pp[2];
        rtk->nfix=0;
    }
    else if (stat==SOLQ_INHERITFIX) {
        for (i=0;i<3;i++) {
            rtk->sol.rr[i]=rtk->xa[i];  /* fix solutions */
            rtk->sol.qr[i]=(float)rtk->Pa[i+i*rtk->na];
        }
        rtk->sol.qr[3]=(float)rtk->Pa[1];
        rtk->sol.qr[4]=(float)rtk->Pa[1+2*rtk->na];
        rtk->sol.qr[5]=(float)rtk->Pa[2];
        rtk->inherix_fixc++;
    }
    for (i=0;i<n;i++) for (j=0;j<nf;j++) {
            if (obs[i].L[j]==0.0) continue;
            rtk->ssat[obs[i].sat-1].pt[obs[i].rcv-1][j]=obs[i].time;
            rtk->ssat[obs[i].sat-1].ph[obs[i].rcv-1][j]=obs[i].L[j];
        }
    for (i=0;i<n;i++) {  /* save rover and base station to satellite's distance */
        if (obs[i].rcv-1==0)
            rtk->ssat[obs[i].sat-1].r0[obs[i].rcv-1]=arc_geodist(rs+i*6,rtk->x,e+i*3);
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
    if (D) free(D); if (Qv) free(Qv); if (ve) free(ve);
    free(rs); free(dts); free(var);
    free(y); free(e); free(azel);
    free(xp); free(Pp); free(xa);
    free(v); free(H); free(R); free(bias);
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
    ddamb_t amb0={0};
    int i;

    arc_log(ARC_INFO,"rtkinit :\n");

    rtk->sol=sol0;
    rtk->opt=*opt;

    for (i=0;i<6;i++) rtk->rb[i]=0.0;
    rtk->nx=opt->mode<=PMODE_FIXED?NX(opt):arc_pppnx(opt);
    rtk->na=opt->mode<=PMODE_FIXED?NR(opt):arc_pppnx(opt);
    rtk->tt=0.0;
    rtk->x=arc_zeros(rtk->nx,1);
    rtk->P=arc_zeros(rtk->nx,rtk->nx);
    rtk->xa=arc_zeros(rtk->na,1);
    rtk->Pa=arc_zeros(rtk->na,rtk->na);
    rtk->xd=arc_zeros(NXDC(&rtk->opt),1);
    rtk->Pd=arc_zeros(NXDC(&rtk->opt),NXDC(&rtk->opt));
    rtk->nfix=rtk->neb=rtk->fixc=rtk->halffix=0;
    for (i=0;i<MAXSAT;i++) {
        rtk->ambc[i]=ambc0;
        rtk->ssat[i]=ssat0;
    }
    for (i=0;i<MAXERRMSG;i++) rtk->errbuf[i]=0;

    /* ceres solver options */
    rtk->ceres_active_x=arc_imat(rtk->nx,1);

    /* ambiguity solver options */
    for (i=0;i<MAXSAT;i++) rtk->amb_index[i]=0;

    for (i=0;i<NUMOFSYS;i++) rtk->prefsat[i]  =rtk->refsat[i]=0;
    for (i=0;i<NUMOFSYS;i++) rtk->ref_delay[i]=0;

    rtk->bias.amb=(ddamb_t*)malloc(sizeof(ddamb_t)*MAXSAT);
    rtk->bias.nmax=MAXSAT;
    rtk->bias.nb=0;
    for (i=0;i<MAXSAT;i++) rtk->bias.amb[i]=amb0;

    rtk->sol.bias.amb=(ddamb_t*)malloc(sizeof(ddamb_t)*MAXSAT);
    rtk->sol.bias.nmax=MAXSAT; rtk->sol.bias.nb=0;
    for (i=0;i<MAXSAT;i++) rtk->sol.bias.amb[i]=amb0;
}
/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void arc_rtkfree(rtk_t *rtk)
{
    arc_log(ARC_INFO,"rtkfree :\n");

    rtk->nx=rtk->na=0;
    if (rtk->x)  free(rtk->x ); rtk->x =NULL;
    if (rtk->P)  free(rtk->P ); rtk->P =NULL;
    if (rtk->xa) free(rtk->xa); rtk->xa=NULL;
    if (rtk->Pa) free(rtk->Pa); rtk->Pa=NULL;
    if (rtk->Pd) free(rtk->Pd); rtk->Pd=NULL;
    if (rtk->xd) free(rtk->xd); rtk->xd=NULL;
    if (rtk->ceres_active_x) {
        free(rtk->ceres_active_x); rtk->ceres_active_x=NULL;
    }
    if (rtk->bias.amb) {
        free(rtk->bias.amb); rtk->bias.nb=rtk->bias.nmax=0;
    }
    if (rtk->sol.bias.amb) {
        free(rtk->sol.bias.amb); rtk->sol.bias.nb=rtk->sol.bias.nmax=0;
    }
}
/* arc single rtk precise positioning ---------------------------------------*/
extern int arc_srtkpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    prcopt_t *opt=&rtk->opt;
    sol_t solb={{0}};
    gtime_t time;
    int i,nu,nr,stat=SOLQ_NONE;
    char msg[128]="";

#ifdef ARC_TEST
    arc_log(ARC_INFO,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@==>%d\n",EPOCH++);
#endif

    arc_log(ARC_INFO,"arc_srtkpos  : time=%s n=%d\n",time_str(obs[0].time,3),n);
    arc_log(ARC_WARNING,"arc_srtkpos : obs=\n"); arc_traceobs(ARC_INFO,obs,n);

    /* set base staion position */
    if (opt->refpos<=POSOPT_RINEX&&opt->mode!=PMODE_SINGLE&&
        opt->mode!=PMODE_MOVEB) {
        for (i=0;i<6;i++) rtk->rb[i]=i<3?opt->rb[i]:0.0;
    }
    /* count rover/base station observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ; /* rover */
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ; /* base */

    /* time difference of base and rover */
    rtk->sol.age=(float)timediff(obs[0].time,obs[nu].time);

    time=rtk->sol.time; /* previous epoch */

    /* rover position by single point positioning */
    if (!arc_pntpos(obs,nu,nav,&rtk->opt,
                    &rtk->sol,NULL,rtk->ssat,msg)) {
        arc_log(ARC_WARNING,"arc_srtkpos : point pos error (%s)\n",msg);
        if (opt->mode==PMODE_STATIC) {
            outsolstat(rtk);
            return 0;
        }
        /* use velecity to estimate current epoch rover station,only use doppler */
        else if (opt->mode==PMODE_KINEMA&&rtk->sol.doppler) {

            for (i=0;i<3;i++) rtk->sol.rr[i]+=(rtk->sol.rr[i+3]
                                               *timediff(rtk->sol.time,time));
            arc_log(ARC_INFO,"arc_srtkpos: rover station dr\n");
            arc_tracemat(ARC_MATPRINTF,rtk->sol.rr,NP(opt),1,10,4);
            rtk->sol.stat=SOLQ_DR;
            outsolstat(rtk);
            return 1;
        }
    }
    arc_matcpy(rtk->sol.prsol.rr,rtk->sol.rr,1,6);

    /* when large time difference of base and rover no positioning */
    if (fabs(rtk->sol.age)>=opt->maxage) {
        arc_log(ARC_WARNING,"large time difference of base and rover and no positioning \n");
    }
    if (time.time!=0) rtk->tt=timediff(rtk->sol.time,time);

    /* single point positioning */
    if (opt->mode==PMODE_SINGLE) {outsolstat(rtk);return 1;}

    /* suppress output of single solution */
    if (!opt->outsingle) {
        rtk->sol.stat=SOLQ_NONE;
    }
    /* check number of data of base station and age of differential */
    if (nr==0) {
        arc_log(ARC_ERROR,"arc_srtkpos : no base station observation data for rtk\n");
        outsolstat(rtk);
        return 1;
    }
    if (opt->mode==PMODE_MOVEB) { /* moving baseline */

        /* estimate position/velocity of base station */
        if (!arc_pntpos(obs+nu,nr,nav,&rtk->opt,&solb,NULL,NULL,msg)) {
            arc_log(ARC_WARNING,"arc_srtkpos : base station position error (%s)\n",msg);
            return 0;
        }
        rtk->sol.age=(float)timediff(rtk->sol.time,solb.time);

        if (fabs(rtk->sol.age)>TTOL_MOVEB) {
            arc_log(ARC_WARNING,"arc_srtkpos : time sync error for "
                    "moving-base (age=%.1f)\n",rtk->sol.age);
            return 0;
        }
        for (i=0;i<6;i++) rtk->rb[i]=solb.rr[i];

        /* time-synchronized position of base station */
        for (i=0;i<3;i++) rtk->rb[i]+=rtk->rb[i+3]*rtk->sol.age;
    }
    else {
        if (fabs(rtk->sol.age)>opt->maxtdiff) {
            arc_log(ARC_WARNING,"arc_srtkpos : age of differential "
                    "error (age=%.1f)\n", rtk->sol.age);
            outsolstat(rtk);
            return 1;
        }
    }
    /* difference pseudorange positioning for initialing */
    if (opt->init_dc&&!opt->dynamics) {

        stat=arc_diff_pr_relpos(rtk,obs,nu,nr,nav);

        if (stat==SOLQ_NONE) {

            arc_log(ARC_WARNING,"arc_srtkpos: "
                    "difference pseudorange positioning for initialing is failed\n");

            if (opt->mode==PMODE_KINEMA&&rtk->sol.doppler) {
                /* use velecity to estimate current epoch rover station,only use doppler */
                for (i=0;i<3;i++) {
                    rtk->sol.prsol.rr[i]+=(rtk->sol.prsol.rr[i+3]*rtk->tt);
                }
                arc_matcpy(rtk->sol.rr,rtk->sol.prsol.rr,3,1);
                arc_log(ARC_INFO,"arc_srtkpos: rover station dr\n");
                arc_tracemat(ARC_MATPRINTF,rtk->sol.rr,NP(opt),1,10,4);
                rtk->sol.stat=SOLQ_DR;
            }
            return 0;
        }
    }
    /* relative potitioning */
    arc_relpos(rtk,obs,nu,nr,nav);
    outsolstat(rtk);
    return 1;
}
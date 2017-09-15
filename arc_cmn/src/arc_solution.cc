
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
 *      Author: SuJingLan
 *********************************************************************************/

/**
 * @file arc_solution.cc
 * @brief Source file for the ARC-SRTK Library
 */

#include "arc.h"

/* constants and macros ------------------------------------------------------*/

#define SQR(x)     ((x)<0.0?-(x)*(x):(x)*(x))
#define SQRT(x)    ((x)<0.0?0.0:sqrt(x))

#define MAXFIELD   64           /* max number of fields in a record */

static const int solq_nmea[]={  /* nmea quality flags to rtklib sol quality */
        /* nmea 0183 v.2.3 quality flags: */
        /*  0=invalid, 1=gps fix (sps), 2=dgps fix, 3=pps fix, 4=rtk, 5=float rtk */
        /*  6=estimated (dead reckoning), 7=manual input, 8=simulation */

        SOLQ_NONE ,SOLQ_SINGLE, SOLQ_DGPS, SOLQ_PPP , SOLQ_FIX,
        SOLQ_FLOAT,SOLQ_DR    , SOLQ_NONE, SOLQ_NONE, SOLQ_NONE
};

static const char *c1[2]={"\033[30m","\033[0m"}; /* black */
static const char *c2[2]={"\033[31m","\033[0m"}; /* red */
static const char *c3[2]={"\033[32m","\033[0m"}; /* green */
static const char *c4[2]={"\033[33m","\033[0m"}; /* yellow */

/* solution option to field separator ----------------------------------------*/
static const char *opt2sep(const solopt_t *opt)
{
    if (!*opt->sep) return " ";
    else if (!strcmp(opt->sep,"\\t")) return "\t";
    return opt->sep;
}
/* separate fields -----------------------------------------------------------*/
static int tonum(char *buff, const char *sep, double *v)
{
    int n,len=(int)strlen(sep);
    char *p,*q;

    for (p=buff,n=0;n<MAXFIELD;p=q+len) {
        if ((q=strstr(p,sep))) *q='\0';
        if (*p) v[n++]=atof(p);
        if (!q) break;
    }
    return n;
}
/* sqrt of covariance --------------------------------------------------------*/
static double sqvar(double covar)
{
    return covar<0.0?-sqrt(-covar):sqrt(covar);
}
/* convert ddmm.mm in nmea format to deg -------------------------------------*/
static double dmm2deg(double dmm)
{
    return floor(dmm/100.0)+fmod(dmm,100.0)/60.0;
}
/* convert time in nmea format to time ---------------------------------------*/
static void septime(double t, double *t1, double *t2, double *t3)
{
    *t1=floor(t/10000.0);
    t-=*t1*10000.0;
    *t2=floor(t/100.0);
    *t3=t-*t2*100.0;
}
/* solution to covariance ----------------------------------------------------*/
static void soltocov(const sol_t *sol, double *P)
{
    P[0]     =sol->qr[0]; /* xx or ee */
    P[4]     =sol->qr[1]; /* yy or nn */
    P[8]     =sol->qr[2]; /* zz or uu */
    P[1]=P[3]=sol->qr[3]; /* xy or en */
    P[5]=P[7]=sol->qr[4]; /* yz or nu */
    P[2]=P[6]=sol->qr[5]; /* zx or ue */
}
/* covariance to solution ----------------------------------------------------*/
static void covtosol(const double *P, sol_t *sol)
{
    sol->qr[0]=(float)P[0]; /* xx or ee */
    sol->qr[1]=(float)P[4]; /* yy or nn */
    sol->qr[2]=(float)P[8]; /* zz or uu */
    sol->qr[3]=(float)P[1]; /* xy or en */
    sol->qr[4]=(float)P[5]; /* yz or nu */
    sol->qr[5]=(float)P[2]; /* zx or ue */
}
/* test solution status ------------------------------------------------------*/
static int test_solstat(const char *buff)
{
    if (strlen(buff)<7||buff[0]!='$') return 0;
    return !strncmp(buff+1,"POS" ,3)||!strncmp(buff+1,"VELACC",6)||
           !strncmp(buff+1,"CLK" ,3)||!strncmp(buff+1,"ION"   ,3)||
           !strncmp(buff+1,"TROP",4)||!strncmp(buff+1,"HWBIAS",6)||
           !strncmp(buff+1,"TRPG",4)||!strncmp(buff+1,"AMB"   ,3)||
           !strncmp(buff+1,"SAT" ,3);
}
/* output solution as the form of x/y/z-ecef ---------------------------------*/
static int outecef(unsigned char *buff, const char *s, const sol_t *sol,
                   const solopt_t *opt)
{
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;

    arc_log(ARC_INFO,"outecef:\n");

    p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f\n",
               s,sep,sol->rr[0],sep,sol->rr[1],sep,sol->rr[2],sep,sol->stat,sep,
               sol->ns,sep,SQRT(sol->qr[0]),sep,SQRT(sol->qr[1]),sep,SQRT(sol->qr[2]),
               sep,sqvar(sol->qr[3]),sep,sqvar(sol->qr[4]),sep,sqvar(sol->qr[5]),
               sep,sol->age,sep,sol->ratio);
    return p-(char *)buff;
}
/* output solution as the form of lat/lon/height -----------------------------*/
static int outpos(unsigned char *buff, const char *s, const sol_t *sol,
                  const solopt_t *opt)
{
    double pos[3],dms1[3],dms2[3],P[9],Q[9];
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;

    arc_log(ARC_INFO,"outpos  :\n");

    ecef2pos(sol->rr,pos);
    soltocov(sol,P);
    covenu(pos,P,Q);
    if (opt->height==1) { /* geodetic height */

    }
    if (opt->degf) {
        deg2dms(pos[0]*R2D,dms1,5);
        deg2dms(pos[1]*R2D,dms2,5);
        p+=sprintf(p,"%s%s%4.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f",s,sep,
                   dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],sep,
                   dms2[2]);
    }
    else p+=sprintf(p,"%s%s%14.9f%s%14.9f",s,sep,pos[0]*R2D,sep,pos[1]*R2D);
    p+=sprintf(p,"%s%10.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f\n",
               sep,pos[2],sep,sol->stat,sep,sol->ns,sep,SQRT(Q[4]),sep,
               SQRT(Q[0]),sep,SQRT(Q[8]),sep,sqvar(Q[1]),sep,sqvar(Q[2]),
               sep,sqvar(Q[5]),sep,sol->age,sep,sol->ratio);
    return p-(char *)buff;
}
/* output solution as the form of e/n/u-baseline -----------------------------*/
static int outenu(unsigned char *buff, const char *s, const sol_t *sol,
                  const double *rb, const solopt_t *opt)
{
    double pos[3],rr[3],enu[3],P[9],Q[9];
    int i;
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;

    arc_log(ARC_INFO,"outenu  :\n");

    for (i=0;i<3;i++) rr[i]=sol->rr[i]-rb[i];
    ecef2pos(rb,pos);
    soltocov(sol,P);
    covenu(pos,P,Q);
    ecef2enu(pos,rr,enu);
    p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f\n",
               s,sep,enu[0],sep,enu[1],sep,enu[2],sep,sol->stat,sep,sol->ns,sep,
               SQRT(Q[0]),sep,SQRT(Q[4]),sep,SQRT(Q[8]),sep,sqvar(Q[1]),
               sep,sqvar(Q[5]),sep,sqvar(Q[2]),sep,sol->age,sep,sol->ratio);
    return p-(char *)buff;
}
/* output processing options ---------------------------------------------------
* output processing options to buffer
* args   : unsigned char *buff IO output buffer
*          prcopt_t *opt    I   processign options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outprcopts(unsigned char *buff, const prcopt_t *opt)
{
    const int sys[]={SYS_GPS,SYS_GLO,SYS_GAL,SYS_QZS,SYS_SBS,0};
    const char *s1[]={"single","dgps","kinematic","static","moving-base","fixed",
                      "ppp-kinematic","ppp-static","ppp-fixed",""};
    const char *s2[]={"L1","L1+L2","L1+L2+L5","L1+L2+L5+L6","L1+L2+L5+L6+L7",
                      "L1+L2+L5+L6+L7+L8",""};
    const char *s3[]={"forward","backward","combined"};
    const char *s4[]={"off","broadcast","sbas","iono-free","estimation",
                      "ionex tec","qzs","lex","vtec_sf","vtec_ef","gtec",""};
    const char *s5[]={"off","saastamoinen","sbas","est ztd","est ztd+grad",""};
    const char *s6[]={"broadcast","precise","broadcast+sbas","broadcast+ssr apc",
                      "broadcast+ssr com","qzss lex",""};
    const char *s7[]={"gps","glonass","galileo","qzss","sbas",""};
    const char *s8[]={"off","continuous","instantaneous","fix and hold",""};
    const char *s9[]={"off","on","auto calib","external calib",""};
    int i;
    char *p=(char *)buff;

    arc_log(ARC_INFO,"outprcopts:\n");

    p+=sprintf(p,"%s pos mode  : %s\n",COMMENTH,s1[opt->mode]);

    if (PMODE_DGPS<=opt->mode&&opt->mode<=PMODE_FIXED) {
        p+=sprintf(p,"%s freqs     : %s\n",COMMENTH,s2[opt->nf-1]);
    }
    if (opt->mode>PMODE_SINGLE) {
        p+=sprintf(p,"%s solution  : %s\n",COMMENTH,s3[opt->soltype]);
    }
    p+=sprintf(p,"%s elev mask : %.1f deg\n",COMMENTH,opt->elmin*R2D);
    if (opt->mode>PMODE_SINGLE) {
        p+=sprintf(p,"%s dynamics  : %s\n",COMMENTH,opt->dynamics?"on":"off");
        p+=sprintf(p,"%s tidecorr  : %s\n",COMMENTH,opt->tidecorr?"on":"off");
    }
    if (opt->mode<=PMODE_FIXED) {
        p+=sprintf(p,"%s ionos opt : %s\n",COMMENTH,s4[opt->ionoopt]);
    }
    p+=sprintf(p,"%s tropo opt : %s\n",COMMENTH,s5[opt->tropopt]);
    p+=sprintf(p,"%s ephemeris : %s\n",COMMENTH,s6[opt->sateph]);
    if (opt->navsys!=SYS_GPS) {
        p+=sprintf(p,"%s navi sys  :",COMMENTH);
        for (i=0;sys[i];i++) {
            if (opt->navsys&sys[i]) p+=sprintf(p," %s",s7[i]);
        }
        p+=sprintf(p,"\n");
    }
    if (PMODE_KINEMA<=opt->mode&&opt->mode<=PMODE_FIXED) {
        p+=sprintf(p,"%s amb res   : %s\n",COMMENTH,s8[opt->modear]);
        if (opt->navsys&SYS_GLO) {
            p+=sprintf(p,"%s amb glo   : %s\n",COMMENTH,s9[opt->glomodear]);
        }
        if (opt->thresar[0]>0.0) {
            p+=sprintf(p,"%s val thres : %.1f\n",COMMENTH,opt->thresar[0]);
        }
    }
    if (opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0) {
        p+=sprintf(p,"%s baseline  : %.4f %.4f m\n",COMMENTH,
                   opt->baseline[0],opt->baseline[1]);
    }
    for (i=0;i<2;i++) {
        if (opt->mode==PMODE_SINGLE||(i>=1&&opt->mode>PMODE_FIXED)) continue;
        p+=sprintf(p,"%s antenna%d  : %-21s (%7.4f %7.4f %7.4f)\n",COMMENTH,
                   i+1,opt->anttype[i],opt->antdel[i][0],opt->antdel[i][1],
                   opt->antdel[i][2]);
    }
    return p-(char *)buff;
}
/* output solution header ------------------------------------------------------
* output solution header to buffer
* args   : unsigned char *buff IO output buffer
*          solopt_t *opt    I   solution options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outsolheads(unsigned char *buff, const solopt_t *opt)
{
    const char *s1[]={"WGS84","Tokyo"},*s2[]={"ellipsoidal","geodetic"};
    const char *s3[]={"GPST","UTC ","JST "},*sep=opt2sep(opt);
    char *p=(char *)buff;
    int timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);

    arc_log(ARC_INFO,"outsolheads:\n");

    if (opt->posf==SOLF_NMEA||opt->posf==SOLF_STAT||opt->posf==SOLF_GSIF) {
        return 0;
    }
    if (opt->outhead) {
        p+=sprintf(p,"%s (",COMMENTH);
        if      (opt->posf==SOLF_XYZ) p+=sprintf(p,"x/y/z-ecef=WGS84");
        else if (opt->posf==SOLF_ENU) p+=sprintf(p,"e/n/u-baseline=WGS84");
        else p+=sprintf(p,"lat/lon/height=%s/%s",s1[opt->datum],s2[opt->height]);
        p+=sprintf(p,",Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)\n");
    }
    p+=sprintf(p,"%s  %-*s%s",COMMENTH,(opt->timef?16:8)+timeu+1,s3[opt->times],sep);

    if (opt->posf==SOLF_LLH) { /* lat/lon/hgt */
        if (opt->degf) {
            p+=sprintf(p,"%16s%s%16s%s%10s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s\n",
                       "latitude(d'\")",sep,"longitude(d'\")",sep,"height(m)",sep,
                       "Q",sep,"ns",sep,"sdn(m)",sep,"sde(m)",sep,"sdu(m)",sep,
                       "sdne(m)",sep,"sdeu(m)",sep,"sdue(m)",sep,"age(s)",sep,"ratio");
        }
        else {
            p+=sprintf(p,"%14s%s%14s%s%10s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s\n",
                       "latitude(deg)",sep,"longitude(deg)",sep,"height(m)",sep,
                       "Q",sep,"ns",sep,"sdn(m)",sep,"sde(m)",sep,"sdu(m)",sep,
                       "sdne(m)",sep,"sdeu(m)",sep,"sdun(m)",sep,"age(s)",sep,"ratio");
        }
    }
    else if (opt->posf==SOLF_XYZ) { /* x/y/z-ecef */
        p+=sprintf(p,"%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s\n",
                   "x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,"Q",sep,"ns",sep,
                   "sdx(m)",sep,"sdy(m)",sep,"sdz(m)",sep,"sdxy(m)",sep,
                   "sdyz(m)",sep,"sdzx(m)",sep,"age(s)",sep,"ratio");
    }
    else if (opt->posf==SOLF_ENU) { /* e/n/u-baseline */
        p+=sprintf(p,"%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s\n",
                   "e-baseline(m)",sep,"n-baseline(m)",sep,"u-baseline(m)",sep,
                   "Q",sep,"ns",sep,"sde(m)",sep,"sdn(m)",sep,"sdu(m)",sep,
                   "sden(m)",sep,"sdnu(m)",sep,"sdue(m)",sep,"age(s)",sep,"ratio");
    }
    return p-(char *)buff;
}
/* std-dev of soltuion -------------------------------------------------------*/
static double sol_std(const sol_t *sol)
{
    /* approximate as max std-dev of 3-axis std-devs */
    if (sol->qr[0]>sol->qr[1]&&sol->qr[0]>sol->qr[2]) return SQRT(sol->qr[0]);
    if (sol->qr[1]>sol->qr[2]) return SQRT(sol->qr[1]);
    return SQRT(sol->qr[2]);
}
/* output solution body --------------------------------------------------------
* output solution body to buffer
* args   : unsigned char *buff IO output buffer
*          sol_t  *sol      I   solution
*          double *rb       I   base station position {x,y,z} (ecef) (m)
*          solopt_t *opt    I   solution options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outsols(unsigned char *buff, const sol_t *sol, const double *rb,
                   const solopt_t *opt)
{
    gtime_t time,ts={0};
    double gpst;
    int week,timeu;
    const char *sep=opt2sep(opt);
    char s[64],str[126];
    unsigned char *p=buff;

    arc_log(ARC_INFO,"outsols :\n");

    /* suppress output if std is over opt->maxsolstd */
    if (opt->maxsolstd>0.0&&sol_std(sol)>opt->maxsolstd) {
        return 0;
    }
    if (opt->posf==SOLF_NMEA) {
        if (opt->nmeaintv[0]<0.0) return 0;
        if (!screent(sol->time,ts,ts,opt->nmeaintv[0])) return 0;
    }
    if (sol->stat<=SOLQ_NONE||(opt->posf==SOLF_ENU&&arc_norm(rb,3)<=0.0)) {
        return 0;
    }
    timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);

    time=sol->time;
    if (opt->times>=TIMES_UTC) time=gpst2utc(time);
    if (opt->times==TIMES_JST) time=timeadd(time,9*3600.0);

    if (opt->timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d%s%*.*f",week,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }

    switch (opt->posf) {
        case SOLF_LLH:  p+=outpos (p,s,sol,opt);   break;
        case SOLF_XYZ:  p+=outecef(p,s,sol,opt);   break;
        case SOLF_ENU:  p+=outenu(p,s,sol,rb,opt); break;
    }
    return p-buff;
}
/* output solution extended ----------------------------------------------------
* output solution exteneded infomation
* args   : unsigned char *buff IO output buffer
*          sol_t  *sol      I   solution
*          ssat_t *ssat     I   satellite status
*          solopt_t *opt    I   solution options
* return : number of output bytes
* notes  : only support nmea
*-----------------------------------------------------------------------------*/
extern int outsolexs(unsigned char *buff, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt)
{
    gtime_t ts={0};
    unsigned char *p=buff;

    arc_log(ARC_INFO,"outsolexs:\n");

    /* suppress output if std is over opt->maxsolstd */
    if (opt->maxsolstd>0.0&&sol_std(sol)>opt->maxsolstd) {
        return 0;
    }
    if (opt->posf==SOLF_NMEA) {
        if (opt->nmeaintv[1]<0.0) return 0;
        if (!screent(sol->time,ts,ts,opt->nmeaintv[1])) return 0;
    }
    return p-buff;
}
/* output processing option ----------------------------------------------------
* output processing option to file
* args   : FILE   *fp       I   output file pointer
*          prcopt_t *opt    I   processing options
* return : none
*-----------------------------------------------------------------------------*/
extern void outprcopt(FILE *fp, const prcopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;

    arc_log(ARC_INFO,"outprcopt:\n");

    if ((n=outprcopts(buff,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output solution header ------------------------------------------------------
* args   : FILE   *fp       I   output file pointer
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
extern void outsolhead(FILE *fp, const solopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;

    arc_log(ARC_INFO,"outsolhead:\n");

    if ((n=outsolheads(buff,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output solution body --------------------------------------------------------
* output solution body to file
* args   : FILE   *fp       I   output file pointer
*          sol_t  *sol      I   solution
*          double *rb       I   base station position {x,y,z} (ecef) (m)
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
extern void outsol(FILE *fp, const sol_t *sol, const double *rb,
                   const solopt_t *opt)
{
#ifdef ARC_NOSTATICS
    return;
#endif
    unsigned char buff[MAXSOLMSG+1];
    int n;

    arc_log(ARC_INFO,"outsol  :\n");

    if ((n=outsols(buff,sol,rb,opt))>0) {
        if (fp==stderr) {
            fprintf(fp,"%s",buff);
        }
        else fwrite(buff,n,1,fp);
    }
}
/* output solution extended ----------------------------------------------------
* output solution exteneded infomation to file
* args   : FILE   *fp       I   output file pointer
*          sol_t  *sol      I   solution
*          ssat_t *ssat     I   satellite status
*          solopt_t *opt    I   solution options
* return : output size (bytes)
* notes  : only support nmea
*-----------------------------------------------------------------------------*/
extern void outsolex(FILE *fp, const sol_t *sol, const ssat_t *ssat,
                     const solopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;

    arc_log(ARC_INFO,"outsolex:\n");

    if ((n=outsolexs(buff,sol,ssat,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
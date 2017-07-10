
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
 * @file arc_core.cc
 * @brief Source file for the ARC-SRTK Library
 * @author SuJingLan
 */

#include "arc.h"
#include "glog/logging.h"
#include "libPF/ParticleFilter.h"
#include "arc_ObservationModel.h"
#include "arc_MovementModel.h"
#include <iomanip>
using namespace std;
/* constants/global variables ------------------------------------------------*/
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))

#define MAXINFILE   1000        /* max number of input files */
#define USEPNTINI   1           /* using the standard position to inital solution*/

static pcvs_t pcvss={0};        /* receiver antenna parameters */
static pcvs_t pcvsr={0};        /* satellite antenna parameters */
static obs_t obss={0};          /* observation data */
static nav_t navs={0};          /* navigation data */

static sta_t stas[MAXRCV];      /* station infomation */
static int nepoch=0;            /* number of observation epochs */
static int iobsu =0;            /* current rover observation data index */
static int iobsr =0;            /* current reference observation data index */
static int isbs  =0;            /* current sbas message index */
static int revs  =0;            /* analysis direction (0:forward,1:backward) */
static int aborts=0;            /* abort status */

static char proc_rov [64]="";   /* rover for current processing */
static char proc_base[64]="";   /* base station for current processing */

typedef libPF::ParticleFilter<ARC::ARC_States> ParticleFilterType;
                                /* particle filter type */

/* show message and check break ----------------------------------------------*/
static int checkbrk(const char *format, ...)
{
    va_list arg;
    char buff[1024],*p=buff;
    if (!*format) return showmsg("");
    va_start(arg,format);
    p+=vsprintf(p,format,arg);
    va_end(arg);
    if (*proc_rov&&*proc_base) sprintf(p," (%s-%s)",proc_rov,proc_base);
    else if (*proc_rov ) sprintf(p," (%s)",proc_rov );
    else if (*proc_base) sprintf(p," (%s)",proc_base);
    return showmsg(buff);
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
static int nextobsb(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i>=0;(*i)--) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i-n>=0;n++) {
        tt=timediff(obs->data[*i-n].time,obs->data[*i].time);
        if (obs->data[*i-n].rcv!=rcv||tt<-DTTOL) break;
    }
    return n;
}
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs, int solq, const prcopt_t *popt,int *nu,int *nr)
{
    gtime_t time={0};
    int i,n=0;

    trace(ARC_INFO,"infunc  : revs=%d iobsu=%d iobsr=%d isbs=%d\n",revs,iobsu,iobsr,isbs);

    if (0<=iobsu&&iobsu<obss.n) {
        settime((time=obss.data[iobsu].time));
        if (checkbrk("processing : %s Q=%d",time_str(time,0),solq)) {
            aborts=1; showmsg("aborted"); return -1;
        }
    }
    if (!revs) { /* input forward data */
        if ((*nu=nextobsf(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(*nr=nextobsf(&obss,&iobsr,2))>0;iobsr+=*nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)>-DTTOL) break;
        }
        else {
            for (i=iobsr;(*nr=nextobsf(&obss,&i,2))>0;iobsr=i,i+=*nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)>DTTOL) break;
        }
        *nr=nextobsf(&obss,&iobsr,2);
        for (i=0;i<*nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu+i];
        for (i=0;i<*nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr+i];
        iobsu+=*nu;
    }
    else { /* input backward data */
        if ((*nu=nextobsb(&obss,&iobsu,1))<=0) return -1;
        if (popt->intpref) {
            for (;(*nr=nextobsb(&obss,&iobsr,2))>0;iobsr-=*nr)
                if (timediff(obss.data[iobsr].time,obss.data[iobsu].time)<DTTOL) break;
        }
        else {
            for (i=iobsr;(*nr=nextobsb(&obss,&i,2))>0;iobsr=i,i-=*nr)
                if (timediff(obss.data[i].time,obss.data[iobsu].time)<-DTTOL) break;
        }
        *nr=nextobsb(&obss,&iobsr,2);
        for (i=0;i<*nu&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsu-*nu+1+i];
        for (i=0;i<*nr&&n<MAXOBS*2;i++) obs[n++]=obss.data[iobsr-*nr+1+i];
        iobsu-=*nu;
    }
    return n;
}
/* carrier-phase bias correction by fcb --------------------------------------*/
static void corr_phase_bias_fcb(obsd_t *obs, int n, const nav_t *nav)
{
    int i,j,k;

    for (i=0;i<nav->nf;i++) {
        if (timediff(nav->fcb[i].te,obs[0].time)<-1E-3) continue;
        if (timediff(nav->fcb[i].ts,obs[0].time)> 1E-3) break;
        for (j=0;j<n;j++) {
            for (k=0;k<NFREQ;k++) {
                if (obs[j].L[k]==0.0) continue;
                obs[j].L[k]-=nav->fcb[i].bias[obs[j].sat-1][k];
            }
        }
        return;
    }
}
/* select common satellites between rover and reference station --------------*/
static int selcomsat(const obsd_t *obs, const rtk_t* rtk,int nu, int nr,
                     const prcopt_t *opt, int *sat, int *iu, int *ir)
{
    int i,j,k=0;

    trace(ARC_INFO,"selsat  : nu=%d nr=%d\n",nu,nr);

    for (i=0,j=nu;i<nu&&j<nu+nr;i++,j++) {
        if      (obs[i].sat<obs[j].sat) j--;
        else if (obs[i].sat>obs[j].sat) i--;
        else if (rtk->ssat[obs[j].sat-1].azel[1]>=opt->elmin) {
            sat[k]=obs[i].sat; iu[k]=i; ir[k++]=j;
            trace(ARC_INFO,"(%2d) sat=%3d iu=%2d ir=%2d\n",k-1,obs[i].sat,i,j);
        }
    }
    return k;
}
/* process positioning -------------------------------------------------------*/
static void procpos(const prcopt_t *popt, const solopt_t *sopt,
                    int mode)
{
    rtk_t rtk;
    obsd_t obs[MAXOBS*2];
    gtime_t time={0};
    double rb[3]={0};
    int i,nobs,n,ns,rsat[MAXSAT],usat[MAXSAT],nu=0,nr=0,sat[MAXSAT],first=1;
    char msg[126];

    trace(ARC_INFO,"procpos : mode=%d\n",mode);

    /* rtk struct date type initial */
    rtkinit(&rtk,popt);

    /* arc-srtk observation model and states movement model  */
    ARC::ARC_ObservationModel ObsModel;
    ARC::ARC_MovementModel    MoveModel(popt,&rtk);
    ARC::ARC_States           States(popt);
    ParticleFilterType        PF(500,&ObsModel,&MoveModel);
    /* set the observation settings and solution data */
    ObsModel.SetOpt(popt);
    ObsModel.SetNav(&navs);
    ObsModel.SetSRTK(&rtk);
    MoveModel.SetNav(&navs);
    MoveModel.SetSRTK(&rtk);

    while ((nobs=inputobs(obs,rtk.sol.stat,popt,&nu,&nr))>=0) {
        /*abort */
        if (aborts) break;

        /* precious epoch time */
        time=rtk.sol.time;

        /* exclude satellites */
        for (i=n=0;i<nobs;i++) {
            if ((satsys(obs[i].sat,NULL)&popt->navsys)&&
                popt->exsats[obs[i].sat-1]!=1) obs[n++]=obs[i];
        }
        if (n<=0) continue;
        /* count rover/base station observations */
        for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;   /* rover */
        for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;   /* base */
        
        /* carrier-phase bias correction */
        if (navs.nf>0) {
            corr_phase_bias_fcb(obs,n,&navs);
        }
#if USEPNTINI
        /* rover position by single point positioning */
        if (!pntpos(obs,nu,&navs,&rtk.opt,&rtk.sol,NULL,rtk.ssat,msg)) {
            trace(ARC_WARNING, "arc-srtk point pos error (%s)\n",msg);
            continue;
        }
        /* inital base station and rover station observation time difference */
        if (time.time!=0) rtk.tt=timediff(rtk.sol.time,time);

        /* to integrate a known prior state use */
        for (i=0;i<rtk.nx;i++) States.SetStatesValue(rtk.x[i],i);
        PF.setPriorState(States);

        /* inital the states standard deviation */
        for (i=0;i<rtk.nx;i++) MoveModel.SetStdX(rtk.P[i],i);
#endif
        /* select common satellites between rover and reference station */
        if ((ns=selcomsat(obs,&rtk,nu,nr,popt,sat,usat,rsat))<=0) continue;

        /* asign common satellites to movement model */
        MoveModel.SetComSatList(rsat,usat,sat,ns);
        /* set gps and bds observation of states movement */
        MoveModel.SetObs(obs,nu+nr);
        /* states movement model for predicting the states */
        MoveModel.drift(States,rtk.tt);

        /* set the states of observation model */
        ObsModel.setStates(States);
        /* set the observations of gps and bds */
        ObsModel.SetObs(obs,n);
        
        /* particle filter */
        PF.filter();
        
        LOG(WARNING)<<"Particle filter position is : "
                    << setiosflags(ios::fixed) << setprecision(10)
                    <<PF.getBestState().getStateValue(0)<<" , "
                    <<PF.getBestState().getStateValue(1)<<" , "
                    <<PF.getBestState().getStateValue(2);
    }
    rtkfree(&rtk);
}
/* read prec ephemeris, sbas data, lex data, tec grid and open rtcm ----------*/
static void readpreceph(char **infile, int n, const prcopt_t *prcopt,
                        nav_t *nav)
{
    seph_t seph0={0};
    int i;
    char *ext;

    trace(ARC_INFO,"readpreceph: n=%d\n",n);

    nav->ne=nav->nemax=0;
    nav->nc=nav->ncmax=0;
    nav->nf=nav->nfmax=0;

    /* read precise ephemeris files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readsp3(infile[i],nav,0);
    }
    /* read precise clock files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        readrnxc(infile[i],nav);
    }
    /* read satellite fcb files */
    for (i=0;i<n;i++) {
        if (strstr(infile[i],"%r")||strstr(infile[i],"%b")) continue;
        if ((ext=strrchr(infile[i],'.'))&&
            (!strcmp(ext,".fcb")||!strcmp(ext,".FCB"))) {
            readfcb(infile[i],nav);
        }
    }
    /* allocate sbas ephemeris */
    nav->ns=nav->nsmax=NSATSBS*2;
    if (!(nav->seph=(seph_t *)malloc(sizeof(seph_t)*nav->ns))) {
        showmsg("error : sbas ephem memory allocation");
        trace(ARC_ERROR,"error : sbas ephem memory allocation");
        return;
    }
    for (i=0;i<nav->ns;i++) nav->seph[i]=seph0;
}
/* free prec ephemeris and sbas data -----------------------------------------*/
static void freepreceph(nav_t *nav)
{
    int i;

    trace(ARC_INFO,"freepreceph:\n");

    free(nav->peph); nav->peph=NULL; nav->ne=nav->nemax=0;
    free(nav->pclk); nav->pclk=NULL; nav->nc=nav->ncmax=0;
    free(nav->fcb ); nav->fcb =NULL; nav->nf=nav->nfmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}
/* read obs and nav data -----------------------------------------------------*/
static int readobsnav(gtime_t ts, gtime_t te, double ti, char **infile,
                      const int *index, int n, const prcopt_t *prcopt,
                      obs_t *obs, nav_t *nav, sta_t *sta)
{
    int i,j,ind=0,nobs=0,rcv=1;

    trace(ARC_INFO,"readobsnav: ts=%s n=%d\n",time_str(ts,0),n);

    obs->data=NULL; obs->n =obs->nmax =0;
    nav->eph =NULL; nav->n =nav->nmax =0;
    nav->geph=NULL; nav->ng=nav->ngmax=0;
    nepoch=0;

    for (i=0;i<n;i++) {
        if (checkbrk("")) return 0;

        if (index[i]!=ind) {
            if (obs->n>nobs) rcv++;
            ind=index[i]; nobs=obs->n;
        }
        /* read rinex obs and nav file */
        if (readrnxt(infile[i],rcv,ts,te,ti,prcopt->rnxopt[rcv<=1?0:1],obs,nav,
                     rcv<=2?sta+rcv-1:NULL)<0) {
            checkbrk("error : insufficient memory");
            trace(ARC_WARNING,"insufficient memory\n");
            return 0;
        }
    }
    if (obs->n<=0) {
        checkbrk("readobsnav : error , no obs data");
        return 0;
    }
    if (nav->n<=0&&nav->ng<=0&&nav->ns<=0) {
        checkbrk("readobsnav : error , no nav data");
        trace(ARC_WARNING,"\n");
        return 0;
    }
    /* sort observation data */
    nepoch=sortobs(obs);

    /* delete duplicated ephemeris */
    uniqnav(nav);
    return 1;
}
/* free obs and nav data -----------------------------------------------------*/
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    trace(ARC_INFO,"freeobsnav:\n");

    free(obs->data); obs->data=NULL; obs->n =obs->nmax =0;
    free(nav->eph ); nav->eph =NULL; nav->n =nav->nmax =0;
    free(nav->geph); nav->geph=NULL; nav->ng=nav->ngmax=0;
    free(nav->seph); nav->seph=NULL; nav->ns=nav->nsmax=0;
}
/* average of single position ------------------------------------------------*/
static int avepos(double *ra, int rcv, const obs_t *obs, const nav_t *nav,
                  const prcopt_t *opt)
{
    obsd_t data[MAXOBS];
    gtime_t ts={0};
    sol_t sol={{0}};
    int i,j,n=0,m,iobs;
    char msg[128];

    trace(ARC_INFO,"avepos: rcv=%d obs.n=%d\n",rcv,obs->n);

    for (i=0;i<3;i++) ra[i]=0.0;

    for (iobs=0;(m=nextobsf(obs,&iobs,rcv))>0;iobs+=m) {

        for (i=j=0;i<m&&i<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j<=0||!screent(data[0].time,ts,ts,1.0)) continue; /* only 1 hz */

        if (!pntpos(data,j,nav,opt,&sol,NULL,NULL,msg)) continue;

        for (i=0;i<3;i++) ra[i]+=sol.rr[i];
        n++;
    }
    if (n<=0) {
        trace(ARC_WARNING,"avepos : no average of base station position\n");
        return 0;
    }
    for (i=0;i<3;i++) ra[i]/=n;
    return 1;
}
/* antenna phase center position ---------------------------------------------*/
static int antpos(prcopt_t *opt, int rcvno, const obs_t *obs, const nav_t *nav,
                  const sta_t *sta, const char *posfile)
{
    double *rr=rcvno==1?opt->ru:opt->rb,del[3],pos[3],dr[3]={0};
    int i,postype=rcvno==1?opt->rovpos:opt->refpos;
    char *name;

    trace(ARC_INFO,"antpos  : rcvno=%d\n",rcvno);

    if (postype==POSOPT_SINGLE) { /* average of single position */
        if (!avepos(rr,rcvno,obs,nav,opt)) {
            trace(ARC_ERROR,"error : station pos computation");
            return 0;
        }
    }
    else if (postype==POSOPT_RINEX) { /* get from rinex header */
        if (norm(stas[rcvno==1?0:1].pos,3)<=0.0) {
            showmsg("error : no position in rinex header");
            trace(ARC_WARNING,"no position position in rinex header\n");
            return 0;
        }
        /* antenna delta */
        if (stas[rcvno==1?0:1].deltype==0) { /* enu */
            for (i=0;i<3;i++) del[i]=stas[rcvno==1?0:1].del[i];
            del[2]+=stas[rcvno==1?0:1].hgt;
            ecef2pos(stas[rcvno==1?0:1].pos,pos);
            enu2ecef(pos,del,dr);
        }
        else { /* xyz */
            for (i=0;i<3;i++) dr[i]=stas[rcvno==1?0:1].del[i];
        }
        for (i=0;i<3;i++) rr[i]=stas[rcvno==1?0:1].pos[i]+dr[i];
    }
    return 1;
}
/* open procssing session ----------------------------------------------------*/
static int openses(const prcopt_t *popt, const solopt_t *sopt,
                   const filopt_t *fopt, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    int i;

    trace(ARC_INFO,"openses :\n");

    /* read satellite antenna parameters */
    if (*fopt->satantp&&!(readpcv(fopt->satantp,pcvs))) {
        showmsg("error : no sat ant pcv in %s",fopt->satantp);
        trace(ARC_WARNING,"sat antenna pcv read error: %s\n",fopt->satantp);
        return 0;
    }
    /* read receiver antenna parameters */
    if (*fopt->rcvantp&&!(readpcv(fopt->rcvantp,pcvr))) {
        showmsg("error : no rec ant pcv in %s",fopt->rcvantp);
        trace(ARC_WARNING,"rec antenna pcv read error: %s\n",fopt->rcvantp);
        return 0;
    }
    /* use satellite L2 offset if L5 offset does not exists */
    for (i=0;i<pcvs->n;i++) {
        if (norm(pcvs->pcv[i].off[2],3)>0.0) continue;
        matcpy(pcvs->pcv[i].off[2],pcvs->pcv[i].off[1], 3,1);
        matcpy(pcvs->pcv[i].var[2],pcvs->pcv[i].var[1],19,1);
    }
    for (i=0;i<pcvr->n;i++) {
        if (norm(pcvr->pcv[i].off[2],3)>0.0) continue;
        matcpy(pcvr->pcv[i].off[2],pcvr->pcv[i].off[1], 3,1);
        matcpy(pcvr->pcv[i].var[2],pcvr->pcv[i].var[1],19,1);
    }
    return 1;
}
/* close procssing session ---------------------------------------------------*/
static void closeses(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr)
{
    trace(ARC_INFO,"closeses:\n");

    /* free antenna parameters */
    free(pcvs->pcv); pcvs->pcv=NULL; pcvs->n=pcvs->nmax=0;
    free(pcvr->pcv); pcvr->pcv=NULL; pcvr->n=pcvr->nmax=0;

    /* free erp data */
    free(nav->erp.data); nav->erp.data=NULL; nav->erp.n=nav->erp.nmax=0;
}
/* set antenna parameters ----------------------------------------------------*/
static void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs,
                   const pcvs_t *pcvr, const sta_t *sta)
{
    pcv_t *pcv;
    double pos[3],del[3];
    int i,j,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
    char id[64];

    /* set satellite antenna parameters */
    for (i=0;i<MAXSAT;i++) {
        if (!(satsys(i+1,NULL)&popt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,pcvs))) {
            satno2id(i+1,id);
            trace(ARC_WARNING,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    for (i=0;i<(mode?2:1);i++) {
        if (!strcmp(popt->anttype[i],"*")) { /* set by station parameters */
            strcpy(popt->anttype[i],sta[i].antdes);
            if (sta[i].deltype==1) { /* xyz */
                if (norm(sta[i].pos,3)>0.0) {
                    ecef2pos(sta[i].pos,pos);
                    ecef2enu(pos,sta[i].del,del);
                    for (j=0;j<3;j++) popt->antdel[i][j]=del[j];
                }
            }
            else { /* enu */
                for (j=0;j<3;j++) popt->antdel[i][j]=stas[i].del[j];
            }
        }
        if (!(pcv=searchpcv(0,popt->anttype[i],time,pcvr))) {
            trace(ARC_ERROR,"no receiver antenna pcv: %s\n",popt->anttype[i]);
            *popt->anttype[i]='\0';
            continue;
        }
        strcpy(popt->anttype[i],pcv->type);
        popt->pcvr[i]=*pcv;
    }
}
/* read ocean tide loading parameters ----------------------------------------*/
static void readotl(prcopt_t *popt, const char *file, const sta_t *sta)
{
    int i,mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;

    for (i=0;i<(mode?2:1);i++) {
        readblq(file,sta[i].name,popt->odisp[i]);
    }
}
/* execute processing session ------------------------------------------------*/
static int execses(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                   const solopt_t *sopt, const filopt_t *fopt, int flag,
                   char **infile, const int *index, int n, char *outfile)
{
    prcopt_t popt_=*popt;
    char path[1024];

    trace(ARC_INFO,"execses : n=%d outfile=%s\n",n,outfile);

    /* read erp data */
    trace(ARC_INFO,"read erp data : %s \n",fopt->eop);
    if (*fopt->eop) {
        if (navs.erp.data) free(navs.erp.data);
        navs.erp.data=NULL; navs.erp.n=navs.erp.nmax=0;
        reppath(fopt->eop,path,ts,"","");
        if (!readerp(path,&navs.erp)) {
            showmsg("error : no erp data %s",path);
            trace(ARC_WARNING,"no erp data %s\n",path);
        }
    }
    /* read obs and nav data */
    trace(ARC_INFO,"read obs and nav data \n");
    if (!readobsnav(ts,te,ti,infile,index,n,
                    &popt_,&obss,&navs,stas)) return 0;

    /* read dcb parameters */
    trace(ARC_INFO,"read dcb parameters : %s \n",fopt->dcb);
    if (*fopt->dcb) {
        reppath(fopt->dcb,path,ts,"","");
        readdcb(path,&navs,stas);
    }
    /* set antenna paramters */
    trace(ARC_INFO,"set antenna paramters \n");
    if (popt_.mode!=PMODE_SINGLE) {
        setpcv(obss.n>0?obss.data[0].time:timeget(),
               &popt_,&navs,&pcvss,&pcvsr,stas);
    }
    /* read ocean tide loading parameters */
    trace(ARC_INFO,"read ocean tide loading parameters \n");
    if (popt_.mode>PMODE_SINGLE&&*fopt->blq) {
        readotl(&popt_,fopt->blq,stas);
    }
    if (PMODE_DGPS<=popt_.mode&&popt_.mode<=PMODE_STATIC) {
        if (!antpos(&popt_,2,&obss,&navs,stas,fopt->stapos)) {
            freeobsnav(&obss,&navs);
            return 0;
        }
    }
    iobsu=iobsr=isbs=revs=aborts=0;

    if (popt_.mode==PMODE_SINGLE||popt_.soltype==0) {
        procpos(&popt_,sopt,0); /* forward */
    }
    else if (popt_.soltype==1) {
        revs=1; iobsu=iobsr=obss.n-1;
        procpos(&popt_,sopt,0); /* backward */
    }
    /* free obs and nav data */
    freeobsnav(&obss,&navs);

    return aborts?1:0;
}
/* execute processing session for each rover ---------------------------------*/
static int execses_r(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov)
{
    int stat=0;

    trace(ARC_INFO,"execses_r: n=%d outfile=%s\n",n,outfile);
    
    /* execute processing session */
    stat=execses(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile);
    return stat;
}
/* execute processing session for each base station --------------------------*/
static int execses_b(gtime_t ts, gtime_t te, double ti, const prcopt_t *popt,
                     const solopt_t *sopt, const filopt_t *fopt, int flag,
                     char **infile, const int *index, int n, char *outfile,
                     const char *rov, const char *base)
{
    int i,stat=0;

    trace(ARC_INFO,"execses_b: n=%d outfile=%s\n",n,outfile);

    /* read prec ephemeris and sbas data */
    readpreceph(infile,n,popt,&navs);

    /* execute processing session */
    stat=execses_r(ts,te,ti,popt,sopt,fopt,flag,infile,index,n,outfile,rov);

    /* free prec ephemeris and sbas data */
    freepreceph(&navs);

    return stat;
}
/* post-processing positioning -------------------------------------------------
* post-processing positioning
* args   : gtime_t ts       I   processing start time (ts.time==0: no limit)
*        : gtime_t te       I   processing end time   (te.time==0: no limit)
*          double ti        I   processing interval  (s) (0:all)
*          double tu        I   processing unit time (s) (0:all)
*          prcopt_t *popt   I   processing options
*          solopt_t *sopt   I   solution options
*          filopt_t *fopt   I   file options
*          char   **infile  I   input files (see below)
*          int    n         I   number of input files
*          char   *outfile  I   output file ("":stdout, see below)
*          char   *rov      I   rover id list        (separated by " ")
*          char   *base     I   base station id list (separated by " ")
* return : status (0:ok,0>:error,1:aborted)
* notes  : input files should contain observation data, navigation data, precise
*          ephemeris/clock (optional), sbas log file (optional), ssr message
*          log file (optional) and tec grid file (optional). only the first
*          observation data file in the input files is recognized as the rover
*          data.
*
*          the type of an input file is recognized by the file extention as ]
*          follows:
*              .sp3,.SP3,.eph*,.EPH*: precise ephemeris (sp3c)
*              .sbs,.SBS,.ems,.EMS  : sbas message log files (rtklib or ems)
*              .lex,.LEX            : qzss lex message log files
*              .rtcm3,.RTCM3        : ssr message log files (rtcm3)
*              .*i,.*I              : tec grid files (ionex)
*              .fcb,.FCB            : satellite fcb
*              others               : rinex obs, nav, gnav, hnav, qnav or clock
*
*          inputs files can include wild-cards (*). if an file includes
*          wild-cards, the wild-card expanded multiple files are used.
*
*          inputs files can include keywords. if an file includes keywords,
*          the keywords are replaced by date, time, rover id and base station
*          id and multiple session analyses run. refer reppath() for the
*          keywords.
*
*          the output file can also include keywords. if the output file does
*          not include keywords. the results of all multiple session analyses
*          are output to a single output file.
*
*          ssr corrections are valid only for forward estimation.
*-----------------------------------------------------------------------------*/
extern int arc_srtk(gtime_t ts, gtime_t te, double ti, double tu,
                    const prcopt_t *popt, const solopt_t *sopt,
                    const filopt_t *fopt, char **infile, int n, char *outfile,
                    const char *rov, const char *base)
{
    int i,stat=0,index[MAXINFILE]={0};

    trace(ARC_INFO,"postpos : ti=%.0f tu=%.0f n=%d outfile=%s\n",ti,tu,n,outfile);

    /* open processing session */
    if (!openses(popt,sopt,fopt,&navs,&pcvss,&pcvsr)) return -1;

    for (i=0;i<n;i++) index[i]=i;

    /* execute processing session */
    stat=execses_b(ts,te,ti,popt,sopt,fopt,1,infile,index,n,outfile,rov,base);

    /* close processing session */
    closeses(&navs,&pcvss,&pcvsr);
    return stat;
}


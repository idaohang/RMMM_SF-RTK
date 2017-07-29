
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

/* constants -----------------------------------------------------------------*/
#define SQR(x)      ((x)*(x))

#define NX          (4+3)       /* # of estimated parameters */

#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define ERR_CBIAS   0.3         /* code bias error std (m) */
#define REL_HUMI    0.7         /* relative humidity for saastamoinen model */

/* pseudorange measurement error variance ------------------------------------*/
static double arc_varerr(const prcopt_t *opt, double el, int sys)
{
    double fact,varr;
    fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
    varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
    return SQR(fact)*varr;
}
/* get tgd parameter (m) -----------------------------------------------------*/
static double arc_gettgd(int sat, const nav_t *nav)
{
    int i;
    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;
        return CLIGHT*nav->eph[i].tgd[0];
    }
    return 0.0;
}
/* psendorange with code bias correction -------------------------------------*/
static double arc_prange(const obsd_t *obs, const nav_t *nav, const double *azel,
                         int iter, const prcopt_t *opt, double *var)
{
    const double *lam=nav->lam[obs->sat-1];
    double PC,P1,P1_P2,P1_C1,gamma;
    int i=0,j=1;
    
    *var=0.0;

    /* test snr mask */
    if (iter>0) {
        if (testsnr(0,i,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) {
            arc_log(ARC_WARNING, "snr mask: %s sat=%2d el=%.1f snr=%.1f\n",
                    time_str(obs->time, 0), obs->sat, azel[1] * R2D, obs->SNR[i] * 0.25);
            return 0.0;
        }
    }
    gamma=SQR(lam[j])/SQR(lam[i]); /* f1^2/f2^2 */
    P1=obs->P[i];
    P1_P2=nav->cbias[obs->sat-1][0];
    P1_C1=nav->cbias[obs->sat-1][1];

    /* single-frequency */
    if (P1==0.0) return 0.0;
    if (obs->code[i]==CODE_L1C) P1+=P1_C1; /* C1->P1 */
    PC=P1-P1_P2/(1.0-gamma);

    *var=SQR(ERR_CBIAS);
    return PC;
}
/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int arc_ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                        const double *azel, int ionoopt, double *ion, double *var)
{
    arc_log(ARC_INFO, "ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
            time_str(time, 3), ionoopt, sat, pos[0] * R2D, pos[1] * R2D, azel[0] * R2D,
            azel[1] * R2D);
    /* broadcast model */
    if (ionoopt==IONOOPT_BRDC) {
        *ion= arc_ionmodel(time, nav->ion_gps, pos, azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* qzss broadcast model */
    if (ionoopt==IONOOPT_QZS&& arc_norm(nav->ion_qzs, 8)>0.0) {
        *ion= arc_ionmodel(time, nav->ion_qzs, pos, azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int arc_tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                        const double *azel, int tropopt, double *trp, double *var)
{
    arc_log(ARC_INFO, "tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
            time_str(time, 3), tropopt, pos[0] * R2D, pos[1] * R2D, azel[0] * R2D,
            azel[1] * R2D);
    
    /* saastamoinen model */
    if (tropopt==TROPOPT_SAAS
        ||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp= arc_tropmodel(time, pos, azel, REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* hopf model */
    else if (tropopt==TROPOPT_HOPF
             ||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=arc_tropmodel_hopf(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* unb3 model */
    else if (tropopt==TROPOPT_UNB3
             ||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=arc_tropmodel_unb3(time,pos,azel,REL_HUMI,NULL,NULL);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* mops model */
    else if (tropopt==TROPOPT_MOPS
             ||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=arc_tropmodel_unb3(time,pos,azel,REL_HUMI,NULL,NULL);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* gcat model */
    else if (tropopt==TROPOPT_GCAT
             ||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=arc_tropmodel_gcat(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* black model */
    else if (tropopt==TROPOPT_BALCK
             ||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=arc_tropmodel_black(time,pos,azel,REL_HUMI,NULL,NULL);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* black model */
    else if (tropopt==TROPOPT_WAAS
             ||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=arc_tropmodel_waas(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}
/* pseudorange residuals -----------------------------------------------------*/
static int arc_rescode(int iter, const obsd_t *obs, int n, const double *rs,
                       const double *dts, const double *vare, const int *svh,
                       const nav_t *nav, const double *x, const prcopt_t *opt,
                       double *v, double *H, double *var, double *azel, int *vsat,
                       double *resp, int *ns)
{
    double r,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P,lam_L1;
    int i,j,nv=0,sys,mask[4]={0};

    arc_log(ARC_INFO, "resprng : n=%d\n", n);
    
    for (i=0;i<3;i++) rr[i]=x[i]; dtr=x[3];
    
    ecef2pos(rr,pos);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        
        if (!(sys=satsys(obs[i].sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
            arc_log(ARC_WARNING, "duplicated observation data %s sat=%2d\n",
                    time_str(obs[i].time, 3), obs[i].sat);
            i++;
            continue;
        }
        /* geometric distance/azimuth/elevation angle */
        if ((r= arc_geodist(rs + i * 6, rr, e))<=0.0||
                arc_satazel(pos, e, azel + i * 2)<opt->elmin) continue;
        
        /* psudorange with code bias correction */
        if ((P=arc_prange(obs+i,nav,azel+i*2,iter,opt,&vmeas))==0.0) continue;
        
        /* excluded satellite? */
        if (satexclude(obs[i].sat,svh[i],opt)) continue;
        
        /* ionospheric corrections */
        if (!arc_ionocorr(obs[i].time, nav, obs[i].sat, pos, azel + i * 2,
                          iter > 0 ? opt->ionoopt : IONOOPT_BRDC, &dion, &vion)) continue;
        
        /* GPS-L1 -> L1/B1 */
        if ((lam_L1=nav->lam[obs[i].sat-1][0])>0.0) {
            dion*=SQR(lam_L1/lam_carr[0]);
        }
        /* tropospheric corrections */
        if (!arc_tropcorr(obs[i].time, nav, pos, azel + i * 2,
                          iter > 0 ? opt->tropopt : TROPOPT_SAAS, &dtrp, &vtrp)) {
            continue;
        }
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        
        /* design matrix */
        for (j=0;j<NX;j++) H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        
        /* time system and receiver bias offset correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else mask[0]=1;
        
        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        /* error variance */
        var[nv++]=arc_varerr(opt,azel[1+i*2],sys)+vare[i]+vmeas+vion+vtrp;

        arc_log(ARC_INFO, "sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n", obs[i].sat,
                azel[i * 2] * R2D, azel[1 + i * 2] * R2D, resp[i], sqrt(var[nv - 1]));
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<4;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int arc_valsol(const double *azel, const int *vsat, int n,
                      const prcopt_t *opt, const double *v, int nv, int nx,
                      char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;

    arc_log(ARC_INFO, "valsol  : n=%d nv=%d\n", n, nv);
    
    /* chi-square validation of residuals */
    vv= arc_dot(v, v, nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        sprintf(msg,"chi-square error nv=%d vv=%.1f"
                " cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        return 0;
    }
    /* large gdop check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    arc_dops(ns, azels, opt->elmin, dop);
    if (dop[0]<=0.0||dop[0]>opt->maxgdop) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
    return 1;
}
/* estimate receiver position ------------------------------------------------*/
static int arc_estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                      const double *vare, const int *svh, const nav_t *nav,
                      const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,
                      double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;

    arc_log(ARC_INFO, "estpos  : n=%d\n", n);
    
    v= arc_mat(n + 4, 1); H= arc_mat(NX, n + 4); var= arc_mat(n + 4, 1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    
    for (i=0;i<MAXITR;i++) {
        
        /* pseudorange residuals */
        nv=arc_rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,v,H,var,azel,vsat,resp,
                       &ns);
        if (nv<NX) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
        /* weight by variance */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* least square estimation */
        if ((info= arc_lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) x[j]+=dx[j];
        
        if (arc_norm(dx, NX)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* glo-gps time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* gal-gps time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* bds-gps time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(unsigned char)ns;
            sol->age=0.0;
            
            /* validate solution */
            if ((stat=arc_valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);
    
    free(v); free(H); free(var);
    return 0;
}
/* raim fde (failure detection and exclution) -------------------------------*/
static int arc_raim_fde(const obsd_t *obs, int n, const double *rs,
                        const double *dts, const double *vare, const int *svh,
                        const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                        double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;

    arc_log(ARC_INFO, "raim_fde: %s n=%2d\n", time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e=arc_mat(6,n); dts_e=arc_mat(2,n);vare_e=arc_mat(1,n); azel_e=arc_zeros(2,n);
    svh_e=arc_imat(1,n); vsat_e= arc_imat(1,n); resp_e= arc_mat(1,n);
    
    for (i=0;i<n;i++) {
        
        /* satellite exclution */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            arc_matcpy(rs_e+6*k,rs+6*j,6,1);
            arc_matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!arc_estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,&sol_e,azel_e,
                        vsat_e,resp_e,msg_e)) {
            arc_log(ARC_ERROR, "raim_fde: exsat=%2d (%s)\n", obs[i].sat, msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            arc_log(ARC_ERROR, "raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                    obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);

        arc_log(ARC_INFO, "raim_fde: exsat=%2d rms=%8.3f\n", obs[i].sat,rms_e);
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            arc_matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        arc_log(ARC_WARNING, "%s: %s excluded by raim\n", tstr+11,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* doppler residuals ---------------------------------------------------------*/
static int arc_resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                      const nav_t *nav, const double *rr, const double *x,
                      const double *azel, const int *vsat, double *v, double *H)
{
    double lam,rate,pos[3],E[9],a[3],e[3],vs[3],cosel;
    int i,j,nv=0;

    arc_log(ARC_INFO, "resdop  : n=%d\n", n);
    
    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        lam=nav->lam[obs[i].sat-1][0];
        
        if (obs[i].D[0]==0.0||lam==0.0
            ||!vsat[i]|| arc_norm(rs + 3 + i * 6, 3)<=0.0) {
            continue;
        }
        /* line-of-sight vector in ecef */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        arc_matmul("TN", 3, 1, 3, 1.0, E, a, 0.0, e);
        
        /* satellite velocity relative to receiver in ecef */
        for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j];
        
        /* range rate with earth rotation correction */
        rate= arc_dot(vs, e, 3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        /* doppler residual */
        v[nv]=-lam*obs[i].D[0]-(rate+x[3]-CLIGHT*dts[1+i*2]);
        
        /* design matrix */
        for (j=0;j<4;j++) H[j+nv*4]=j<3?-e[j]:1.0;
        nv++;
    }
    return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
static void arc_estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                       const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                       const double *azel, const int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H;
    int i,j,nv;

    arc_log(ARC_INFO, "estvel  : n=%d\n", n);
    
    v=arc_mat(n,1); H=arc_mat(4,n);
    
    for (i=0;i<MAXITR;i++) {
        
        /* doppler residuals */   
        if ((nv=arc_resdop(obs,n,rs,dts,nav,sol->rr,
                           x,azel,vsat,v,H))<4) {
            break;
        }
        /* least square estimation */
        if (arc_lsq(H,v,4,nv,dx,Q)) {
			break;
		}
        for (j=0;j<4;j++) x[j]+=dx[j];
        if (arc_norm(dx,4)<1E-6) {
            for (i=0;i<3;i++) sol->rr[i+3]=x[i];
            break;
        }
    }
    free(v); free(H);
}
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
* notes  : assuming sbas-gps, galileo-gps, qzss-gps, compass-gps time offset and
*          receiver bias are negligible (only involving glonass-gps time offset
*          and receiver bias)
*-----------------------------------------------------------------------------*/
extern int arc_pntpos(const obsd_t *obs, int n, const nav_t *nav,
                      const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                      char *msg)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    
    sol->stat=SOLQ_NONE;
    
    if (n<=0) {strcpy(msg,"no observation data"); return 0;}

    arc_log(ARC_INFO, "pntpos  : tobs=%s n=%d\n", time_str(obs[0].time, 3), n);
    
    sol->time=obs[0].time; if (msg) msg[0]='\0';
    
    rs=arc_mat(6,n); dts=arc_mat(2,n); var=arc_mat(1,n);
    azel_=arc_zeros(2,n); resp=arc_mat(1,n);
    
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
#if 0
        opt_.sateph =EPHOPT_BRDC;
#endif
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positons, velocities and clocks */
    arc_satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position with pseudorange */
    stat=arc_estpos(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    
    /* raim fde */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=arc_raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,sol,azel_,vsat,resp,msg);
    }
    /* estimate receiver velocity with doppler */
    if (stat) arc_estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
            ssat[i].snr[0]=0;
        }
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            ssat[obs[i].sat-1].snr[0]=obs[i].SNR[0];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}

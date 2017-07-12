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
 *  Created on: July 08, 2017
 *      Author: SuJingLan
 *********************************************************************************/

/**
 * @class ARC_ObservationMoel
 * @file arc_Obsservation.cpp
 * @brief SRTK Observation Model source file
 * @author SuJingLan
 */

#include "arc.h"
#include "arc_ObservationModel.h"
#include "glog/logging.h"

///////////////////////////////////////////////////////////////////////////////////
/* single-differenced observable ---------------------------------------------*/
static double arc_sdobs(const obsd_t *obs, int i, int j, int f)
{
    double pi=f<NFREQ?obs[i].L[f]:obs[i].P[f-NFREQ];
    double pj=f<NFREQ?obs[j].L[f]:obs[j].P[f-NFREQ];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return norm(dr,3);
}
/* select common satellites between rover and reference station --------------*/
static int arc_selsat(const obsd_t *obs, double *azel, int nu, int nr,
                      const prcopt_t *opt, int *sat, int *iu, int *ir)
{
    int i,j,k=0;
    for (i=0,j=nu;i<nu&&j<nu+nr;i++,j++) {
        if      (obs[i].sat<obs[j].sat) j--;
        else if (obs[i].sat>obs[j].sat) i--;
        else if (azel[1+j*2]>=opt->elmin) { 
            sat[k]=obs[i].sat; iu[k]=i; ir[k++]=j;
        }
    }
    return k;
}
/* undifferenced phase/code residual for satellite ---------------------------*/
static void arc_zdres_sat(int base, double r, const obsd_t *obs,const nav_t *nav,
                          const double *azel, const double *dant,const prcopt_t *opt,
                          double *y,double dion)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=0,nf=1;

    if (lam[i]==0.0) return;

    /* check snr mask */
    if (testsnr(base,i,azel[1],obs->SNR[i]*0.25,
                &opt->snrmask)) return;

    /* residuals = observable - pseudorange */
    if (obs->L[i]!=0.0) y[i   ]=obs->L[i]*lam[i]-r-dant[i]+dion;
    if (obs->P[i]!=0.0) y[i+nf]=obs->P[i]       -r-dant[i]-dion;
}
/* undifferenced phase/code residuals ----------------------------------------*/
static int arc_zdres(int base, const obsd_t *obs, int n, const double *rs,
                     const double *dts, const int *svh, const nav_t *nav,
                     const double *rr, const prcopt_t *opt, int index, double *y,
                     double *e, double *azel)
{
    double r,rr_[3],pos[3],dant[NFREQ]={0},disp[3];
    double zhd,zazel[]={0.0,90.0*D2R},dion=0.0,vion=0.0;
    int i=0,nf=1;

    for (i=0;i<n*nf*2;i++) y[i]=0.0;

    if (norm(rr,3)<=0.0) return 0; /* no receiver position */

    for (i=0;i<3;i++) rr_[i]=rr[i];

    /* earth tide correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),
                 rr_,opt->tidecorr,&nav->erp,
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
        antmodel(opt->pcvr+index,opt->antdel[index],azel+i*2,opt->posopt[1],
                 dant);
        /* undifferenced phase/code residual for satellite */
        arc_zdres_sat(base,r,obs+i,nav,azel+i*2,dant,opt,y+i*nf*2,dion);
    }
    return 1;
}
/* test valid observation data -----------------------------------------------*/
static int validobs(int i, int j, int f, int nf, double *y)
{
    /* if no phase observable, psudorange is also unusable */
    return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
           (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
}
/* test navi system (m=0:gps/qzs/sbs,1:glo,2:gal,3:bds) ----------------------*/
static int test_sys(int sys, int m)
{
    switch (sys) {
        case SYS_GPS: return m==0;
        case SYS_QZS: return m==0;
        case SYS_SBS: return m==0;
        case SYS_GLO: return m==1;
        case SYS_GAL: return m==2;
        case SYS_CMP: return m==3;
    }
    return 0;
}
/* double-differenced phase/code residuals -----------------------------------*/
static int arc_ddres(rtk_t *rtk, const nav_t *nav, double dt, const double *x,
                     const int *sat, double *y,double *azel, const int *iu,
                     const int *ir, int ns, double *v)
{
    prcopt_t *opt=&rtk->opt;
    double posu[3],posr[3];
    double lami,lamj;
    int i,j,m,f,ff,nv=0,sysi,sysj,nf=1;

    ecef2pos(x,posu); ecef2pos(opt->rb,posr);

    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].resp[j]=rtk->ssat[i].resc[j]=0.0;
    }
    for (m=0;m<4;m++)

        for (f=0;f<2*nf;f++) {

            /* search reference satellite with highest elevation */
            for (i=-1,j=0;j<ns;j++) {
                sysi=rtk->ssat[sat[j]-1].sys;
                if (!test_sys(sysi,m)) continue;
                if (!validobs(iu[j],ir[j],f,nf,y)) continue;
                if (i<0||azel[1+iu[j]*2]>=azel[1+iu[i]*2]) i=j;
            }
            if (i<0) continue;

            /* make double difference */
            for (j=0;j<ns;j++) {
                if (i==j) continue;
                sysj=rtk->ssat[sat[j]-1].sys;
                if (!test_sys(sysj,m)) continue;
                if (!validobs(iu[j],ir[j],f,nf,y)) continue;

                ff=f%nf;
                lami=nav->lam[sat[i]-1][ff];
                lamj=nav->lam[sat[j]-1][ff];
                if (lami<=0.0||lamj<=0.0) continue;

                /* double-differenced residual */
                v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])-(y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);

                /* double-differenced phase-bias term */
                if (f<nf) {
                    v[nv]-=lami*x[IB(sat[i],f,opt)]-lamj*x[IB(sat[j],f,opt)];
                }
                nv++;
            }
        }
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

    if (nb==0||fabs(tt)<DTTOL) {
        nb=n; for (i=0;i<n;i++) obsb[i]=obs[i];
        return tt;
    }
    ttb=timediff(time,obsb[0].time);
    if (fabs(ttb)>opt->maxtdiff*2.0||ttb==tt) return tt;

    satposs(time,obsb,nb,nav,opt->sateph,rs,dts,var,svh);

    if (!arc_zdres(1,obsb,nb,rs,dts,svh,nav,rtk->rb,opt,1,yb,e,azel)) {
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
/* single-difference residual ------------------------------------------------*/
static int arc_sdres(rtk_t *rtk, const nav_t *nav, double dt, const double *x,
                     const int *sat, double *y,double *azel, const int *iu,
                     const int *ir, int ns, double *v)
{
    prcopt_t *opt=&rtk->opt;
    double posu[3],posr[3];
    double lamj;
    int i,j,m,f,ff,nv=0,sysj,nf=1;

    ecef2pos(x,posu); ecef2pos(opt->rb,posr);

    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].resp[0]=rtk->ssat[i].resc[0]=0.0;
    }
    for (m=0;m<4;m++)

        for (f=0;f<2*nf;f++) {

            /* make double difference */
            for (j=0;j<ns;j++) {
                sysj=rtk->ssat[sat[j]-1].sys;
                if (!test_sys(sysj,m)) continue;
                if (!validobs(iu[j],ir[j],f,nf,y)) continue;

                ff=f%nf;
                lamj=nav->lam[sat[j]-1][ff];
                if (lamj<=0.0) continue;

                /* single-differenced residual */
                v[nv]=y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2];

                /* double-differenced phase-bias term */
                if (f<nf) {
                    v[nv]-=lamj*x[IB(sat[j],f,opt)];
                }
                nv++;
            }
        }
    return nv;
}
/* relative positioning ------------------------------------------------------*/
static int arc_relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                      const nav_t *nav,const double *rr,double *ddy,int *nddy,
                      double *sdy,int* nsdy)
{
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,dt;
    int i=0,n=nu+nr,ns,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT];
    int svh[MAXOBS*2];
    int nf=1;

    dt=timediff(time,obs[nu].time);
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n);
    y=mat(nf*2,n); e=mat(3,n);
    azel=zeros(2,n);

    /* initial satellite status informations */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        rtk->ssat[i].vsat[0]=0;
        rtk->ssat[i].snr [0]=0;
    }
    /* compute the satellite position and velecitys */
    satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    if (!arc_zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,svh+nu,nav,opt->rb,opt,1,
                  y+nu*nf*2,e+nu*3,azel+nu*2)) {
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    if (opt->intpref) {
        dt=arc_intpres(time,obs+nu,nr,nav,rtk,y+nu*nf*2);
    }
    if ((ns=arc_selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    arc_zdres(0,obs,nu,rs,dts,svh,nav,rr,opt,0,y,e,azel);

    *nddy=arc_ddres(rtk,nav,dt,rr,sat,y,azel,iu,ir,ns,ddy);
    *nsdy=arc_sdres(rtk,nav,dt,rr,sat,y,azel,iu,ir,ns,sdy);
    free(rs); free(dts); free(var);
    free(y); free(e); free(azel);
    return 1;
}
static int arc_measure(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                       const double *rr,double *ddy,int *nddy, double* sdy,int *nsdy)
{
    int nu,nr;

    /* count rover/base station observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;   /* rover */
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;   /* base */

    /* compute relative positioning double-difference residuals,modified from rtklib */
    arc_relpos(rtk,obs,nu,nr,nav,rr,ddy,nddy,sdy,nsdy);
    return 1;
}
///////////////////////////////////////////////////////////////////////////////
namespace ARC {
    ARC_ObservationModel::ARC_ObservationModel() :
            libPF::ObservationModel<ARC_States>() {
    }

    ARC_ObservationModel::ARC_ObservationModel(const ARC_OPT *OPT,
                                               const ARC_OBSD *OBS,
                                               const ARC_NAV *NAV,
                                               ARC_RTK* SRTK,
                                               int NObs)
            : libPF::ObservationModel<ARC_States>(){
        m_OPT=OPT;
        m_RTK=SRTK;
        m_OBS=OBS;
        m_NAV=NAV;
        m_NObs=NObs;
        USE_SD_OR_DD=1;
    }
    ARC_ObservationModel::~ARC_ObservationModel() {
    }
    double ARC_ObservationModel::measure(const ARC_States &state) const {

        int NDDY,NSDY;
        static double DDY[MAXSAT];
        static double SDY[MAXSAT];
        static double Xp[MAXPFSTETAS];
        double sum=0.0;

        ARC_ASSERT_TRUE(Exception,m_RTK->nx > 0,"Particle Filter States is Zero");

        for (int i=0;i<m_RTK->nx;i++) Xp[i]=state.getStateValue(i);
        for (int i=0;i<m_RTK->nx;i++) {
            m_RTK->x[i]=state.getStateValue(i);
        }
        arc_measure(m_RTK,m_OBS,m_NObs,m_NAV,Xp,DDY,&NDDY,SDY,&NSDY);

        if (USE_SD_OR_DD==1) {
            for (int i=0; i<NDDY;i++) sum+=SQR(DDY[i]);
            return 1.0/SQRT(sum/double(NDDY));
        }
        else if (USE_SD_OR_DD==0) {
            for (int i=0; i<NSDY;i++) sum+=SQR(SDY[i]);
            return 1.0/SQRT(sum/double(NDDY));
        }
    }
}


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
static double sdobs(const obsd_t *obs, int i, int j, int f)
{
    double pi=f<NFREQ?obs[i].L[f]:obs[i].P[f-NFREQ];
    double pj=f<NFREQ?obs[j].L[f]:obs[j].P[f-NFREQ];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* single-differenced measurement error variance -----------------------------*/
static double varerr(int sat, int sys, double el, double bl, double dt, int f,
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
    return 2.0*(opt->ionoopt==IONOOPT_IFLC?3.0:1.0)*
                   (a*a+b*b/sinel/sinel+c*c)+d*d;
}
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return norm(dr,3);
}
/* select common satellites between rover and reference station --------------*/
static int selsat(const obsd_t *obs, double *azel, int nu, int nr,
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
static void zdres_sat(int base, double r, const obsd_t *obs,
                      const nav_t *nav,
                      const double *azel, const double *dant,
                      const prcopt_t *opt, double *y)
{
    const double *lam=nav->lam[obs->sat-1];
    int i,nf=NF(opt);

    for (i=0;i<nf;i++) {
        if (lam[i]==0.0) continue;

        /* check snr mask */
        if (testsnr(base,i,azel[1],obs->SNR[i]*0.25,
                    &opt->snrmask)) {
            continue;
        }
        /* residuals = observable - pseudorange */
        if (obs->L[i]!=0.0) y[i   ]=obs->L[i]*lam[i]-r-dant[i];
        if (obs->P[i]!=0.0) y[i+nf]=obs->P[i]       -r-dant[i];
    }
}
/* undifferenced phase/code residuals ----------------------------------------*/
static int zdres(int base, const obsd_t *obs, int n, const double *rs,
                 const double *dts, const int *svh, const nav_t *nav,
                 const double *rr, const prcopt_t *opt, int index, double *y,
                 double *e, double *azel)
{
    double r,rr_[3],pos[3],dant[NFREQ]={0},disp[3];
    double zhd,zazel[]={0.0,90.0*D2R};
    int i,nf=NF(opt);

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

        /* receiver antenna phase center correction */
        antmodel(opt->pcvr+index,opt->antdel[index],azel+i*2,opt->posopt[1],
                 dant);
        /* undifferenced phase/code residual for satellite */
        zdres_sat(base,r,obs+i,nav,azel+i*2,dant,opt,y+i*nf*2);
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
/* double-differenced measurement error covariance ---------------------------*/
static void ddcov(const int *nb, int n, const double *Ri, const double *Rj,
                  int nv, double *R)
{
    int i,j,k=0,b;

    for (i=0;i<nv*nv;i++) R[i]=0.0;
    for (b=0;b<n;k+=nb[b++]) {
        for (i=0;i<nb[b];i++) for (j=0;j<nb[b];j++) {
            R[k+i+(k+j)*nv]=Ri[k+i]+(i==j?Rj[k+i]:0.0);
        }
    }
}
/* precise tropspheric model -------------------------------------------------*/
 static double prectrop(gtime_t time, const double *pos, int r,
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
static int ddres(rtk_t *rtk, const nav_t *nav, double dt, const double *x,
                 const int *sat, double *y,
                 double *azel, const int *iu, const int *ir, int ns, double *v,
                 double *R)
{
    prcopt_t *opt=&rtk->opt;
    double bl,dr[3],posu[3],posr[3],didxi=0.0,didxj=0.0,*im;
    double *tropr,*tropu,*dtdxr,*dtdxu,*Ri,*Rj,lami,lamj,fi,fj;
    int i,j,k,m,f,ff,nv=0,nb[NFREQ*4*2+2]={0},b=0,sysi,sysj,nf=NF(opt);

    bl=baseline(x,rtk->rb,dr);
    ecef2pos(x,posu); ecef2pos(opt->rb,posr);

    Ri=mat(ns*nf*2+2,1); Rj=mat(ns*nf*2+2,1); im=mat(ns,1);
    tropu=mat(ns,1); tropr=mat(ns,1); dtdxu=mat(ns,3); dtdxr=mat(ns,3);

    for (i=0;i<MAXSAT;i++) for (j=0;j<NFREQ;j++) {
        rtk->ssat[i].resp[j]=rtk->ssat[i].resc[j]=0.0;
    }
    /* compute factors of ionospheric and tropospheric delay */
    for (i=0;i<ns;i++) {
        if (opt->ionoopt>=IONOOPT_EST) {
            im[i]=(ionmapf(posu,azel+iu[i]*2)+ionmapf(posr,azel+ir[i]*2))/2.0;
        }
        if (opt->tropopt>=TROPOPT_EST) {
            tropu[i]=prectrop(rtk->sol.time,posu,0,azel+iu[i]*2,opt,x,dtdxu+i*3);
            tropr[i]=prectrop(rtk->sol.time,posr,1,azel+ir[i]*2,opt,x,dtdxr+i*3);
        }
    }
    for (m=0;m<4;m++)

        for (f=0;f<nf*2;f++) {

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
                sysi=rtk->ssat[sat[i]-1].sys;
                sysj=rtk->ssat[sat[j]-1].sys;
                if (!test_sys(sysj,m)) continue;
                if (!validobs(iu[j],ir[j],f,nf,y)) continue;

                ff=f%nf;
                lami=nav->lam[sat[i]-1][ff];
                lamj=nav->lam[sat[j]-1][ff];
                if (lami<=0.0||lamj<=0.0) continue;

                /* double-differenced residual */
                v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])-
                      (y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);
                /* double-differenced ionospheric delay term */
                if (opt->ionoopt==IONOOPT_EST) {
                    fi=lami/lam_carr[0]; fj=lamj/lam_carr[0];
                    didxi=(f<nf?-1.0:1.0)*fi*fi*im[i];
                    didxj=(f<nf?-1.0:1.0)*fj*fj*im[j];
                    v[nv]-=didxi*x[II(sat[i],opt)]-didxj*x[II(sat[j],opt)];
                }
                /* double-differenced tropospheric delay term */
                if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                    v[nv]-=(tropu[i]-tropu[j])-(tropr[i]-tropr[j]);
                }
                /* double-differenced phase-bias term */
                if (f<nf) {
                    v[nv]-=lami*x[IB(sat[i],f,opt)]-lamj*x[IB(sat[j],f,opt)];
                }
                /* single-differenced measurement error variances */
                Ri[nv]=varerr(sat[i],sysi,azel[1+iu[i]*2],bl,dt,f,opt);
                Rj[nv]=varerr(sat[j],sysj,azel[1+iu[j]*2],bl,dt,f,opt);
                nv++;
                nb[b]++;
            }
            b++;
        }
    ddcov(nb,b,Ri,Rj,nv,R);
    free(Ri); free(Rj); free(im);
    free(tropu); free(tropr); free(dtdxu); free(dtdxr);
    return nv;
}
/* time-interpolation of residuals (for post-mission) ------------------------*/
static double intpres(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                      rtk_t *rtk, double *y)
{
    static obsd_t obsb[MAXOBS];
    static double yb[MAXOBS*NFREQ*2],rs[MAXOBS*6],dts[MAXOBS*2],var[MAXOBS];
    static double e[MAXOBS*3],azel[MAXOBS*2];
    static int nb=0,svh[MAXOBS*2];
    prcopt_t *opt=&rtk->opt;
    double tt=timediff(time,obs[0].time),ttb,*p,*q;
    int i,j,k,nf=NF(opt);

    if (nb==0||fabs(tt)<DTTOL) {
        nb=n; for (i=0;i<n;i++) obsb[i]=obs[i];
        return tt;
    }
    ttb=timediff(time,obsb[0].time);
    if (fabs(ttb)>opt->maxtdiff*2.0||ttb==tt) return tt;

    satposs(time,obsb,nb,nav,opt->sateph,rs,dts,var,svh);

    if (!zdres(1,obsb,nb,rs,dts,svh,nav,rtk->rb,opt,1,yb,e,azel)) {
        return tt;
    }
    for (i=0;i<n;i++) {
        for (j=0;j<nb;j++) if (obsb[j].sat==obs[i].sat) break;
        if (j>=nb) continue;
        for (k=0,p=y+i*nf*2,q=yb+j*nf*2;k<nf*2;k++,p++,q++) {
            if (*p==0.0||*q==0.0) *p=0.0; else *p=(ttb*(*p)-tt*(*q))/(ttb-tt);
        }
    }
    return fabs(ttb)>fabs(tt)?ttb:tt;
}
/* relative positioning ------------------------------------------------------*/
static int relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                  const nav_t *nav,double *sdy,double *ddy,int *nsdy,int *nddy,
                  double *DDR)
{
    prcopt_t *opt=&rtk->opt;
    gtime_t time=obs[0].time;
    double *rs,*dts,*var,*y,*e,*azel,*v,*R,dt;
    int i,j,n=nu+nr,ns,ny,nv,sat[MAXSAT],iu[MAXSAT],ir[MAXSAT];
    int svh[MAXOBS*2];
    int nf=1;

    dt=timediff(time,obs[nu].time);
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n);
    y=mat(nf*2,n); e=mat(3,n);
    azel=zeros(2,n);

    /* initial satellite status informations */
    for (i=0;i<MAXSAT;i++) {
        rtk->ssat[i].sys=satsys(i+1,NULL);
        for (j=0;j<NFREQ;j++) rtk->ssat[i].vsat[j]=0;
        for (j=1;j<NFREQ;j++) rtk->ssat[i].snr [j]=0;
    }
    /* compute the satellite position and velecitys */
    satposs(time,obs,n,nav,opt->sateph,rs,dts,var,svh);

    if (!zdres(1,obs+nu,nr,rs+nu*6,dts+nu*2,svh+nu,nav,opt->rb,opt,1,
               y+nu*nf*2,e+nu*3,azel+nu*2)) {
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    if (opt->intpref) {
        dt=intpres(time,obs+nu,nr,nav,rtk,y+nu*nf*2);
    }
    if ((ns=selsat(obs,azel,nu,nr,opt,sat,iu,ir))<=0) {
        free(rs); free(dts); free(var); free(y); free(e); free(azel);
        return 0;
    }
    ny=ns*nf*2+2;
    v=mat(ny,1); R=mat(ny,ny);

    zdres(0,obs,nu,rs,dts,svh,nav,rtk->x,opt,0,y,e,azel);
    nv=ddres(rtk,nav,dt,rtk->x,sat,y,azel,iu,ir,ns,v,R);
    for (i=0;i<nv;i++) ddy[i]=v[i];
    for (i=0;i<n;i++) {
        for (j=0;j<2;j++) sdy[2*i+j]=y[i*2+j];
    }
    for (i=0;i<nv;i++) DDR[i]=R[i+nv*i];
    *nsdy=n; *nddy=nv;
    free(rs); free(dts); free(var);
    free(y); free(e); free(azel); free(v); free(R);
    return 1;
}
static int arc_measure(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
                       double *sdy,double *ddy,int *nsdy,int *nddy,double *DDR)
{
    int nu,nr;

    /* count rover/base station observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++) ;   /* rover */
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++) ;   /* base */

    /* compute relative positioning double-difference residuals,modified from rtklib */
    relpos(rtk,obs,nu,nr,nav,sdy,ddy,nsdy,nddy,DDR);
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
    }
    ARC_ObservationModel::~ARC_ObservationModel() {
    }
    double ARC_ObservationModel::measure(const ARC_States &state) const {

        double Y[MAXSAT*2];
        double DDY[MAXSAT];
        double DDR[MAXSAT];

        ARC_ASSERT_TRUE(Exception,m_RTK->nx>0,"the states'numbers is zero");
        
        double *xb=mat(1,m_RTK->nx);
        int NY,NDDY;
        matcpy(xb,m_RTK->x,m_RTK->nx,1);
        for (int i=0;i<state.getStatesNum();i++) {
            m_RTK->x[i]=state.getStateValue(i);
        }
        arc_measure(m_RTK,m_OBS,m_NObs,m_NAV,Y,DDY,&NY,&NDDY,DDR);
        matcpy(m_RTK->x,xb,m_RTK->nx,1);
        double sum=0.0;
        for (int i=0;i<NDDY;i++) sum+=SQR(DDY[i])/DDR[i];
        ARC_ASSERT_TRUE(Exception,sum>0.0,"Warning : division by zero");
        free(xb);
        return 1.0/SQRT(sum);
    }
}


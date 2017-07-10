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
 * @file arc_MovementModel.cc
 * @brief Source file for the ARC_MovementModel class.
 * @author SuJingLan
 */

#include "arc.h"
#include "arc_MovementModel.h"

/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
    int j;
    rtk->x[i]=xi;
    for (j=0;j<rtk->nx;j++) {
        rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
    }
}
/* single-differenced observable ---------------------------------------------*/
static double sdobs(const obsd_t *obs, int i, int j, int f)
{
    double pi=f<NFREQ?obs[i].L[f]:obs[i].P[f-NFREQ];
    double pj=f<NFREQ?obs[j].L[f]:obs[j].P[f-NFREQ];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* temporal update of position/velocity/acceleration -------------------------*/
static void arc_udpos(rtk_t *rtk, double tt)
{
    double *F,*FP,*xp,pos[3],Q[9]={0},Qv[9],var=0.0;
    int i,j;

    /* initialize position for first epoch */
    if (norm(rtk->x,3)<=0.0) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        if (rtk->opt.dynamics) {
            for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
            for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        }
    }
    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        return;
    }
    /* check variance of estimated postion */
    for (i=0;i<3;i++) var+=rtk->P[i+i*rtk->nx]; var/=3.0;

    if (var>VAR_POS) {
        /* reset position with large variance */
        for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i);
        for (i=3;i<6;i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i);
        for (i=6;i<9;i++) initx(rtk,1E-6,VAR_ACC,i);
        return;
    }
    /* state transition of position/velocity/acceleration */
    F=eye(rtk->nx); FP=mat(rtk->nx,rtk->nx); xp=mat(rtk->nx,1);

    for (i=0;i<6;i++) {
        F[i+(i+3)*rtk->nx]=tt;
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN",rtk->nx,1,rtk->nx,1.0,F,rtk->x,0.0,xp);
    matcpy(rtk->x,xp,rtk->nx,1);
    matmul("NN",rtk->nx,rtk->nx,rtk->nx,1.0,F,rtk->P,0.0,FP);
    matmul("NT",rtk->nx,rtk->nx,rtk->nx,1.0,FP,F,0.0,rtk->P);

    /* process noise added to only acceleration */
    Q[0]=Q[4]=SQR(rtk->opt.prn[3]); Q[8]=SQR(rtk->opt.prn[4]);
    ecef2pos(rtk->x,pos);
    covecef(pos,Q,Qv);
    for (i=0;i<3;i++) for (j=0;j<3;j++) {
        rtk->P[i+6+(j+6)*rtk->nx]+=Qv[i+j*3];
    }
    free(F); free(FP); free(xp);
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void arc_udion(rtk_t *rtk, double tt, double bl, const int *sat, int ns)
{
    double el,fact;
    int i,j;

    for (i=1;i<=MAXSAT;i++) {
        j=II(i,&rtk->opt);
        if (rtk->x[j]!=0.0&&
            rtk->ssat[i-1].outc[0]>GAP_RESION&&rtk->ssat[i-1].outc[1]>GAP_RESION)
            rtk->x[j]=0.0;
    }
    for (i=0;i<ns;i++) {
        j=II(sat[i],&rtk->opt);

        if (rtk->x[j]==0.0) {
            initx(rtk,1E-6,SQR(rtk->opt.std[1]*bl/1E4),j);
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

    for (i=0;i<2;i++) {
        j=IT(i,&rtk->opt);

        if (rtk->x[j]==0.0) {
            initx(rtk,INIT_ZWD,SQR(rtk->opt.std[2]),j); /* initial zwd */

            if (rtk->opt.tropopt>=TROPOPT_ESTG) {
                for (k=0;k<2;k++) initx(rtk,1E-6,VAR_GRA,++j);
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
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv)
{
    unsigned int slip,LLI;
    int f,sat=obs[i].sat;

    for (f=0;f<rtk->opt.nf;f++) {

        if (obs[i].L[f]==0.0) continue;

        /* restore previous LLI */
        if (rcv==1) LLI=getbitu(&rtk->ssat[sat-1].slip[f],0,2); /* rover */
        else        LLI=getbitu(&rtk->ssat[sat-1].slip[f],2,2); /* base  */

        /* detect slip by cycle slip flag in LLI */
        if (rtk->tt>=0.0) { /* forward */
            slip=obs[i].LLI[f];
        }
        else { /* backward */
            slip=LLI;
        }
        /* detect slip by parity unknown flag transition in LLI */
        if (((LLI&2)&&!(obs[i].LLI[f]&2))||(!(LLI&2)&&(obs[i].LLI[f]&2))) {
            slip|=1;
        }
        /* save current LLI */
        if (rcv==1) setbitu(&rtk->ssat[sat-1].slip[f],0,2,obs[i].LLI[f]);
        else        setbitu(&rtk->ssat[sat-1].slip[f],2,2,obs[i].LLI[f]);

        /* save slip and half-cycle valid flag */
        rtk->ssat[sat-1].slip[f]|=(unsigned char)slip;
        rtk->ssat[sat-1].half[f]=(obs[i].LLI[f]&2)?0:1;
    }
}
/* temporal update of phase biases -------------------------------------------*/
static void arc_udbias(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav)
{
    double cp,pr,*bias,offset,lami;
    int i,j,f,slip,reset,nf=NF(&rtk->opt);

    for (i=0;i<ns;i++) {

        /* detect cycle slip by LLI */
        for (f=0;f<rtk->opt.nf;f++) rtk->ssat[sat[i]-1].slip[f]&=0xFC;
        detslp_ll(rtk,obs,iu[i],1);
        detslp_ll(rtk,obs,ir[i],2);

        /* update half-cycle valid flag */
        for (f=0;f<nf;f++) {
            rtk->ssat[sat[i]-1].half[f]=
                    !((obs[iu[i]].LLI[f]&2)||(obs[ir[i]].LLI[f]&2));
        }
    }
    for (f=0;f<nf;f++) {
        /* reset phase-bias if instantaneous AR or expire obs outage counter */
        for (i=1;i<=MAXSAT;i++) {

            reset=++rtk->ssat[i-1].outc[f]>(unsigned int)rtk->opt.maxout;

            if (rtk->opt.modear==ARMODE_INST&&rtk->x[IB(i,f,&rtk->opt)]!=0.0) {
                initx(rtk,0.0,0.0,IB(i,f,&rtk->opt));
            }
            else if (reset&&rtk->x[IB(i,f,&rtk->opt)]!=0.0) {
                initx(rtk,0.0,0.0,IB(i,f,&rtk->opt));
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
        bias=zeros(ns,1);

        /* estimate approximate phase-bias by phase - code */
        for (i=j=0,offset=0.0;i<ns;i++) {
            if (rtk->opt.ionoopt!=IONOOPT_IFLC) {
                cp=sdobs(obs,iu[i],ir[i],f); /* cycle */
                pr=sdobs(obs,iu[i],ir[i],f+NFREQ);
                lami=nav->lam[sat[i]-1][f];
                if (cp==0.0||pr==0.0||lami<=0.0) continue;
                bias[i]=cp-pr/lami;
            }
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
            initx(rtk,bias[i],SQR(rtk->opt.std[0]),IB(sat[i],f,&rtk->opt));
        }
        free(bias);
    }
}
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
    int i;
    for (i=0;i<3;i++) dr[i]=ru[i]-rb[i];
    return norm(dr,3);
}
/* temporal update of states --------------------------------------------------*/
static void arc_udstate(rtk_t *rtk, const obsd_t *obs, const int *sat,
                        const int *iu, const int *ir, int ns, const nav_t *nav,
                        double tt)
{
    double dr[3]={0},bl;
    /* temporal update of position/velocity/acceleration */
    arc_udpos(rtk,tt);

    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt>=IONOOPT_EST) {
        bl=baseline(rtk->x,rtk->rb,dr);
        arc_udion(rtk,tt,bl,sat,ns);
    }
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt>=TROPOPT_EST) {
        arc_udtrop(rtk,tt,bl);
    }
    arc_udbias(rtk,tt,obs,sat,iu,ir,ns,nav);
}
/////////////////////////////////////////////////////////////////////////////////
/// \brief ARC Main namespace of this package.
namespace ARC {
    
    ARC_MovementModel::ARC_MovementModel() :
            libPF::MovementModel<ARC_States>() {
        m_RNG = new libPF::CRandomNumberGenerator();
    }
    ARC_MovementModel::~ARC_MovementModel() {
        if (m_RNG) delete m_RNG;
        if (StdX) delete StdX;
    }
    ARC_MovementModel::ARC_MovementModel(const ARC_OPT *OPT,ARC_RTK* SRTK):
            libPF::MovementModel<ARC_States>() {
        NX=SRTK->nx;
        StdX=mat(1,NX);
        m_RNG = new libPF::CRandomNumberGenerator();
    }
    void ARC_MovementModel::drift(ARC_States &state, double dt) const
    {
        ARC_ASSERT_TRUE(Exception,m_SRTK->nx>0,"states numbers is less "
                                               "zero,can check function:rtkinit()");
        if (Ns<=0) return;
        arc_udstate(m_SRTK,m_OBS,SatList,m_RoverSat,m_BaseSat,Ns,m_NAV,dt);
        for (int i=0;i<state.getStatesNum();i++) {
            state.SetStatesValue(m_SRTK->x[i],i);
        }
        for (int i=0;i<m_SRTK->nx;i++) StdX[i]=SQRT(m_SRTK->P[i+i*m_SRTK->nx]);
    }
    void ARC_MovementModel::diffuse(ARC_States &state, double dt) const
    {
        ARC_ASSERT_TRUE(Exception,m_SRTK->nx>0,"states numbers is less zero")
        for (int i=0;i<m_SRTK->nx;i++) {
            state.SetStatesValue(state.getStateValue(i)
                                 +m_RNG->getGaussian(StdX[i])*dt,i);
        }
    }
}
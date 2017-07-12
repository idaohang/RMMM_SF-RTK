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

/* single-differenced observable ---------------------------------------------*/
static double sdobs(const obsd_t *obs, int i, int j, int f)
{
    double pi=f<NFREQ?obs[i].L[f]:obs[i].P[f-NFREQ];
    double pj=f<NFREQ?obs[j].L[f]:obs[j].P[f-NFREQ];
    return pi==0.0||pj==0.0?0.0:pi-pj;
}
/* temporal update of position/velocity/acceleration -------------------------*/
static void arc_udpos(int nx,double *X, double tt,const rtk_t *rtk)
{
    int i;

    /* initial states value from standard positioning*/
    for (i=0;i<3;i++) rtk->sol.rr[i];
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
                       const int *iu, const int *ir, int ns, const nav_t *nav,
                       double *X)
{
    double cp,pr,*bias,offset,lami;
    int i,j,f,slip,reset,nf=1;

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
                X[IB(i,f,&rtk->opt)]=0.0;
            }
            else if (reset&&rtk->x[IB(i,f,&rtk->opt)]!=0.0) {
                X[IB(i,f,&rtk->opt)]=0.0;
            }
            if (rtk->opt.modear!=ARMODE_INST&&reset) {
                rtk->ssat[i-1].lock[f]=-rtk->opt.minlock;
            }
        }
        /* reset phase-bias if detecting cycle slip */
        for (i=0;i<ns;i++) {
            j=IB(sat[i],f,&rtk->opt);
            slip=rtk->ssat[sat[i]-1].slip[f];
            if (rtk->opt.modear==ARMODE_INST||!(slip&1)) continue;
            X[j]=0.0;
            rtk->ssat[sat[i]-1].lock[f]=-rtk->opt.minlock;
        }
        bias=zeros(ns,1);

        /* estimate approximate phase-bias by phase - code */
        for (i=j=0,offset=0.0;i<ns;i++) {
            cp=sdobs(obs,iu[i],ir[i],f); /* cycle */
            pr=sdobs(obs,iu[i],ir[i],f+NFREQ);
            lami=nav->lam[sat[i]-1][f];
            if (cp==0.0||pr==0.0||lami<=0.0) continue;
            bias[i]=cp-pr/lami;

            if (X[IB(sat[i],f,&rtk->opt)]!=0.0) {
                offset+=bias[i]-X[IB(sat[i],f,&rtk->opt)];
                j++;
            }
        }
        /* correct phase-bias offset to enssure phase-code coherency */
        if (j>0) {
            for (i=1;i<=MAXSAT;i++) {
                if (X[IB(i,f,&rtk->opt)]!=0.0) X[IB(i,f,&rtk->opt)]+=offset/j;
            }
        }
        /* set initial states of phase-bias */
        for (i=0;i<ns;i++) {
            if (bias[i]==0.0||X[IB(sat[i],f,&rtk->opt)]!=0.0) continue;
            X[IB(sat[i],f,&rtk->opt)]=bias[i];
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
                        double tt,double *X)
{
    /* temporal update of position/velocity/acceleration */
    arc_udpos(rtk->nx,X,tt,rtk);

    /* temporal update of phase-bias */
    arc_udbias(rtk,tt,obs,sat,iu,ir,ns,nav,X);
}
/////////////////////////////////////////////////////////////////////////////////
/// \brief ARC Main namespace of this package.
namespace ARC {
    
    ARC_MovementModel::ARC_MovementModel() :
            libPF::MovementModel<ARC_States>() {
        PF_ROVERPOS_STD=ARC_PF_ROVERPOS_STD;
        PF_AMB_MIN=ARC_PF_AMB_MIN;
        PF_AMB_MAX=ARC_PF_AMB_MAX;
        m_RNG = new libPF::CRandomNumberGenerator();
    }
    ARC_MovementModel::~ARC_MovementModel() {
        if (m_RNG) delete m_RNG;
    }
    ARC_MovementModel::ARC_MovementModel(const ARC_OPT *OPT,ARC_RTK* SRTK):
            libPF::MovementModel<ARC_States>() {
        m_RNG = new libPF::CRandomNumberGenerator();
    }
    void ARC_MovementModel::drift(ARC_States &state, double dt) const
    {
        if (Ns<=0) return;
        arc_udstate(m_SRTK,m_OBS,SatList,m_RoverSat,m_BaseSat,Ns,m_NAV,dt,state.getStatesVal());
    }
    void ARC_MovementModel::diffuse(ARC_States &state, double dt) const
    {
        for (int i=0;i<m_SRTK->nx;i++) {
            if (state.getStateValue(i)==0.0)
                continue;
            if (i<3) state.SetStatesValue(state.getStateValue(i)+m_RNG->getGaussian(PF_ROVERPOS_STD),i);
            else     state.SetStatesValue(double(int(state.getStateValue(i)+m_RNG->getUniform(PF_AMB_MIN,PF_AMB_MAX))),i);
        }
    }
}
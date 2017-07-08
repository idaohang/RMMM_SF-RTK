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
#include "arc_ObservationModel.h"

namespace ARC {
    ARC_ObservationModel::ARC_ObservationMoel() : libPF::ObservationModel<ARC_States>() {
    }

    ARC_ObservationModel::~ARC_ObservationModel() {
    }
    void ARC_ObservationModel::setTrueCarState(const ARC_States& state)
    {
        m_TrueCarState = state;
    }
    double ARC_ObservationModel::measure(const ARC_States &state) const
    {
        double v[MAXSAT]={0.0};
        for (size_t i=0;i<MAXSAT;i++) {
            
        }
    }
    void ARC_ObservationModel::ComputeZd()
    {
        
    }
    void ARC_ObservationModel::ComputeSatPos()
    {
        double *SatPos=mat(6,Nu+Nb),*SatClk=mat(2,Nu+Nb),*Var=mat(1,Nu+Nb);
        int Svh[MAXSAT*2];
        satposs(m_RoverTime,m_Obs,Nu+Nb,&m_Nav,m_Opt.sateph,SatPos,SatClk,Var,Svh);
        for (int i=0;i<Nu+Nb;i++) {
            if (SatList[i]) {
                for (int j=0;j<6;i++) {
                    m_SatPos[(SatList[i]-1)*6+j]=SatPos[i*6+j];
                }
                for (int j=0;j<2;j++) {
                    m_SatClk[(SatList[i]-1)*2+j]=SatClk[i*2+j];
                }
                m_SatPosVar[SatList[i]-1]=Var[i];
                m_SVH[SatList[i]-1]=Svh[i];
            }
        }
        delete SatPos; delete SatClk; delete Var;
    }
}


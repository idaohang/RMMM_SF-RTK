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

#include "arc_MovementModel.h"
#include "arc_PF.h"
#include <libPF/CRandomNumberGenerator.h>
/// \brief ARC Main namespace of this package.

namespace ARC {
    ARC_MovementModel::ARC_MovementModel() : libPF::MovementModel<ARC_States>(),
                                             m_Std_BaseX(ARC_PF_BASEPOS_STD),
                                             m_Std_BaseY(ARC_PF_BASEPOS_STD),
                                             m_Std_BaseZ(ARC_PF_BASEPOS_STD),
                                             m_Std_RoverX(ARC_PF_ROVERPOS_STD),
                                             m_Std_RoverY(ARC_PF_ROVERPOS_STD),
                                             m_Std_RoverZ(ARC_PF_ROVERPOS_STD),
                                             m_Std_Rover_Clk(ARC_PF_ROVERCLK_STD),
                                             m_Std_Base_Clk(ARC_PF_BASECLK_STD),
                                             m_Std_Base_Trp(ARC_PF_BASETROP_STD),
                                             m_Std_Rover_Trp(ARC_PF_ROVERTROP_STD) {
        m_RNG = new libPF::CRandomNumberGenerator();
        for (size_t i = 0; i < MAXSAT; i++) m_Std_N[i] = ARC_PF_AMB_STD;
    }

    ARC_MovementModel::~ARC_MovementModel() {
        if (m_RNG) delete m_RNG;
    }

    void ARC_MovementModel::drift(ARC_States &state, double dt) const
    {
        state.setRoverPosX(state.getRoverPosX()+state.getRoverVelX()*dt);
        state.setRoverPosY(state.getRoverPosY()+state.getRoverVelY()*dt);
        state.setRoverPosZ(state.getRoverPosZ()+state.getRoverVelZ()*dt);
    }
    void ARC_MovementModel::diffuse(ARC_States &state, double dt) const
    {
        state.setBasePosX(state.getBasePosX()+m_RNG->getGaussian(m_Std_BaseX)*dt);
        state.setBasePosY(state.getBasePosY()+m_RNG->getGaussian(m_Std_BaseY)*dt);
        state.setBasePosZ(state.getBasePosZ()+m_RNG->getGaussian(m_Std_BaseZ)*dt);
        state.setRoverPosX(state.getRoverPosX()+m_RNG->getGaussian(m_Std_RoverX)*dt);
        state.setRoverPosY(state.getRoverPosY()+m_RNG->getGaussian(m_Std_RoverY)*dt);
        state.setRoverPosZ(state.getRoverPosZ()+m_RNG->getGaussian(m_Std_RoverZ)*dt);

        state.setBaseVelX(state.getBaseVelX()+m_RNG->getGaussian(m_Std_Base_VelX)*dt);
        state.setBaseVelY(state.getBaseVelY()+m_RNG->getGaussian(m_Std_Base_VelY)*dt);
        state.setBaseVelZ(state.getBaseVelZ()+m_RNG->getGaussian(m_Std_Base_VelZ)*dt);
        state.setRoverVelX(state.getRoverVelX()+m_RNG->getGaussian(m_Std_Rover_VelX)*dt);
        state.setRoverVelX(state.getRoverVelX()+m_RNG->getGaussian(m_Std_Rover_VelY)*dt);
        state.setRoverVelX(state.getRoverVelX()+m_RNG->getGaussian(m_Std_Rover_VelZ)*dt);

        state.setRoverClk(state.getRoverClk()+m_RNG->getGaussian(m_Std_Rover_Clk)*dt);
        state.setBaseClk(state.getRoverClk()+m_RNG->getGaussian(m_Std_Base_Clk)*dt);

        state.setRoverTrop(state.getRoverTrop()+m_RNG->getGaussian(m_Std_Rover_Trp)*dt);
        state.setBaseTrop(state.getBaseTrop()+m_RNG->getGaussian(m_Std_Base_Trp)*dt);
        state.setRoverIon(state.getRoverIon()+m_RNG->getGaussian(m_Std_Rover_Ion)*dt);
        state.setBaseIon(state.getBaseIon()+m_RNG->getGaussian(m_Std_Base_Ion)*dt);

        for (size_t i=0;i<MAXSAT;i++) {
            state.setAmb(state.getAmb(i+1)+m_RNG->getGaussian(m_Std_N[i])*dt,i+1);
        }

    }
}
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

ARC_MovementModel::ARC_MovementModel() : libPF::MovementModel<ARC_States>() ,
                                         m_Std_BaseX(ARC_PF_BASEPOS_STD),
                                         m_Std_BaseY(ARC_PF_BASEPOS_STD),
                                         m_Std_BaseZ(ARC_PF_BASEPOS_STD),
                                         m_Std_RoverX(ARC_PF_ROVERPOS_STD),
                                         m_Std_RoverY(ARC_PF_ROVERPOS_STD),
                                         m_Std_RoverZ(ARC_PF_ROVERPOS_STD),
                                         m_Std_Rover_Clk(ARC_PF_ROVERCLK_STD),
                                         m_Std_Base_Clk(ARC_PF_BASECLK_STD),
                                         m_Std_Base_Trp(ARC_PF_BASETROP_STD),
                                         m_Std_Rover_Trp(ARC_PF_ROVERTROP_STD)
{
    m_RNG = new libPF::CRandomNumberGenerator();
    for (size_t i=0;i<MAXSAT;i++) m_Std_N[i]=ARC_PF_AMB_STD;
}

ARC_MovementModel::~ARC_MovementModel() {
    if (m_RNG) delete m_RNG;
}
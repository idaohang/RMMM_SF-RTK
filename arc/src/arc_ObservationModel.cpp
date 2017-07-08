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
    int ARC_ObservationModel::ComputeZd()
    {
        
    }
    int ARC_ObservationModel::ComputeSatPos()
    {
        int svh[MAXOBS*2];
        double *var=mat(1,Nu+Nb);
        
    }
}


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
 * @file arc_Obsservation.h
 * @brief SRTK Observation Model
 * @author SuJingLan
 */

#ifndef ARC_ARC_OBSERVATIONMODEL_H
#define ARC_ARC_OBSERVATIONMODEL_H

#include "libPF/ObservationModel.h"
#include "arc_States.h"
#include "arc_PF.h"

namespace ARC {

    class ARC_ObservationModel:
            public libPF::ObservationModel<ARC_States> {
        /**
         * Constructor
         */
    public:
        ARC_ObservationModel();
        ARC_ObservationModel(const ARC_OPT *OPT,const ARC_OBSD* OBS,
                             const ARC_NAV *NAV,ARC_RTK * SRTK,
                             int Nobs);
        /**
         * Destructor
         */
        ~ARC_ObservationModel();
        /**
         * @param state Reference to the state that has to be weightened.
         * @return weight for the given state.
         */
        double measure(const ARC_States& state) const;
        /**
         * \brief set the true states of ARC-SRTK Observation Model
         * @param state the true state of the ARC-SRTK Observation Model
         */
        void setStates(const ARC_States& state){
            m_TrueArcState=state;
        }
        /// \brief get the states value given index
        double getStateVal(int index) {
            return m_TrueArcState.getStateValue(index);
        }
        /// \brief set the settings of arc-srtk
        void SetOpt(const ARC_OPT* OPT) {
            m_OPT=OPT;
        }
        /// \brief set the navigation data
        void SetNav(const ARC_NAV* NAV) {
            m_NAV=NAV;
        }
        /// \brief set the observation data
        void SetObs(const ARC_OBSD* OBS,int Nobs) {
            m_OBS=OBS; m_NObs=Nobs;
        }
        /// \brief set the arc-srtk solution data
        void SetSRTK(ARC_RTK* SRTK) {
            m_RTK=SRTK;
        }
    private:
        /// \brief the true States of Observation Model
        ARC_States m_TrueArcState;
        /// \brief the numbers of observations
        int m_NObs;
        /// \brief arc-srtk solution type
        ARC_RTK* m_RTK;
        /// \brief arc-srtk settings
        const ARC_OPT*  m_OPT;
        /// \brief arc-srtk observations
        const ARC_OBSD* m_OBS;
        /// \brief arc-srtk navigation
        const ARC_NAV*  m_NAV;
    };
}
#endif //ARC_ARC_OBSERVATIONMODEL_H

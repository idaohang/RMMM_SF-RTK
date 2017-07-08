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

#include <libPF/ObservationModel.h>
#include "arc_States.h"
#include "arc_PF.h"

namespace ARC {

    class ARC_ObservationModel: public libPF::ObservationModel<ARC_States> {

        /**
         * Constructor
         */
        ARC_ObservationModel();

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
        void setTrueCarState(const ARC_States& state);
    private:
        /// \brief the true States of Observation Model
        ARC_States m_TrueCarState;
   
        /// \brief the navigation data,(x,y,z,vx,vy,vz){m,m/s}
        double m_SatPos[MAXSAT*6];
        /// \brief save all satelite information
        ARC_SAT m_SatInfo[MAXSAT];
        /// \brief navigation data
        ARC_NAV m_Nav;
        /// \brief observation data
        ARC_OBSD m_Obs[MAXSAT*2];
        /// \brief processing options type
        ARC_OPT m_Opt;
        /// \brief hydro‐static mapping function of rover and base station
        double m_THMap[MAXSAT][2];
        /// \brief wet mapping function of rover and base station
        double m_TWMap[MAXSAT][2];
        /// \brief tropospheric zenith total delay (m) of base and rover station
        double m_TH[MAXSAT][2];
        /// \brief tropospheric zenith hydro‐static delay (m) of base and rover station
        double m_TW[MAXSAT][2];
        /// \brief the vertical ionospheric delay in L1 frequency (m) of base and rover station
        double m_Ion[MAXSAT][2];
        /// \brief the ionosphere map function
        double m_IonMap[MAXSAT][2];
        /// \brief the ambguity of gps and bds satelites
        double m_N[MAXSAT];
        /// \brief the measurement variance matrix
        double m_R[MAXSAT];

        /// \brief the number of rover and base station obsservations
        int Nb,Nu;

        /// \brief the flag satellite heathy
        int m_SVH[MAXSAT];

        /// \brief the satellite clock
        double m_SatClk[MAXSAT*2];

        /// \brief variance of satellite position (ECEF)
        double m_SatPosVar[MAXSAT];

        /// \brief the observation time of base station
        ARC_Time m_BaseTime;
        /// \brief the observation time of rover station
        ARC_Time m_RoverTime;

        /// \brief this epoch satellite no. list
        int SatList[MAXSAT];

        /**
         * Method
         */
        /**
         * \brief compute the undifferenced phase/code residuals
         */
        void ComputeZd();
        void ComputeSatPos();


    };
}
#endif //ARC_ARC_OBSERVATIONMODEL_H

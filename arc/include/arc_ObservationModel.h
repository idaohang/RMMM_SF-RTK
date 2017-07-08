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

namespace ARC {
    class ARC_ObservationMoel: public libPF::ObservationModel<ARC_States> {
        /**
             * Constructor
             */
        ARC_ObservationMoel();

        /**
         * Destructor
         */
        ~ARC_ObservationMoel();

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
    };
}
#endif //ARC_ARC_OBSERVATIONMODEL_H

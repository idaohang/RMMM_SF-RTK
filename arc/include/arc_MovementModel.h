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

#ifndef ARC_ARC_MOVEMENTMODEL_H
#define ARC_ARC_MOVEMENTMODEL_H
#include <libPF/MovementModel.h>
#include <libPF/RandomNumberGenerationStrategy.h>
#include "arc_States.h"
#include "arc_PF.h"
/**
 * @class ARC_MovementModel
 *
 * @brief single frequency rtk position movement model propagates class
 *
 * This movement model propagates a car state according to the translation and
 * rotation speed.
 *
 * @author sujinglan
 */
namespace ARC {
    class ARC_MovementModel : public libPF::MovementModel<ARC_States> {
        /**
         * Constructor
         */
        ARC_MovementModel();
        /**
         * Destructor
         */
        ~ARC_MovementModel();
        /**
         * The drift method propagates the car using its speed.
         * @param state Pointer to the state that has to be manipulated.
         */
        void drift(ARC_States &state, double dt) const;

        /**
         * The diffusion consists of a very small gaussian jitter on the
         * state's variable.
         * @param state Pointer to the state that has to be manipulated.
         */
        void diffuse(ARC_States &state, double dt) const;

        /// \brief set the standard deviation of index-th states
        void SetStdX(double Std,int Index){
            StdX[Index]=Std;
        }
        /// \brief get the standard deviation of index-th states
        double getStdX(int Index) const {
            return StdX[Index];
        }
    private:
        /// \brief standard deviation for the states
        double *StdX;
        /// \brief the numbers of estimate states
        int NX;
        /// \brief Stores the random number generator
        libPF::RandomNumberGenerationStrategy* m_RNG;
        /// \brief arc-srtk solution type
        ARC_RTK* m_SRTK;
        /// \brief observation data
        ARC_OBSD* m_OBS;
        /// \brief navigation data type
        ARC_NAV* m_NAV;
        /// \brief common sat numbers
        int Ns;
        /// \brief base and rover station satellite list
        int m_RoverSat[MAXSAT];
        int m_BaseSat[MAXSAT];
        /// \brief All satellite list
        int SatList[MAXSAT];
    };
}
#endif //ARC_ARC_MOVEMENTMODEL_H

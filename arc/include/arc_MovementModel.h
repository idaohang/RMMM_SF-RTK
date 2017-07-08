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

        /**
         * @param d new standard deviation for the diffusion of x,y,z
         *          of base and rover station position
         */
        inline void setBaesXStdDev(double d) { m_Std_BaseX = d; }

        inline void setBaseYstdDev(double d) { m_Std_BaseY = d; }

        inline void setBaseZStdDev(double d) { m_Std_BaseZ = d; }

        inline void setRoverXStdDev(double d) { m_Std_RoverX = d; }

        inline void setRoverYStdDev(double d) { m_Std_RoverY = d; }

        inline void setRoverZStdDev(double d) { m_Std_RoverZ = d; }

        /**
         * @param d new standard deviation for the diffusion of clock drift
         *          of base and rover station
         */
        inline void setBaseClkStdDev(double d) { m_Std_Base_Clk = d; }

        inline void setRoverClkStdDev(double d) { m_Std_Rover_Clk = d; }

        /**
         * @param d new standard deviation for the diffusion of trosphere delay
         *          of base and rover station
         */
        inline void setBaseTrpStdDev(double d) { m_Std_Base_Trp = d; }

        inline void setRoverTrpStdDev(double d) { m_Std_Rover_Trp = d; }

        /**
         * @param d new standard deviation for the diffusion of ionsphere delay
         *          of base and rover station
         */
        inline void setBaseIonStdDev(double d) {m_Std_Base_Ion=d;}
        inline void setRoverIonStDev(double d) {m_Std_Rover_Ion=d;}

        /**
         * @param d new standard deviation for the diffusion of ambguity
         *          of single difference ambguity on base and rover station
         */
        inline void setNStdDev(double *d) {
            for (size_t i = 0; i < MAXSAT; i++) m_Std_N[i] = d[i];
        }

        /**
         * @param d new standard deviation for the diffusion of base and rover station velecity
         */
        inline void setBaseStdDevX(double d) {m_Std_Base_VelX=d;}
        inline void setBaseStdDevY(double d) {m_Std_Base_VelY=d;}
        inline void setBaseStdDevZ(double d) {m_Std_Base_VelZ=d;}

        inline void setRoverStdDevX(double d) {m_Std_Rover_VelX=d;}
        inline void setRoverStdDevY(double d) {m_Std_Rover_VelY=d;}
        inline void setRoverStdDevZ(double d) {m_Std_Rover_VelZ=d;}

    private:
        /// \brief standard deviation for the base and rover position
        double m_Std_BaseX, m_Std_BaseY, m_Std_BaseZ;
        double m_Std_RoverX, m_Std_RoverY, m_Std_RoverZ;
        /// \brief standard deviation for the base and rover clock drift
        double m_Std_Base_Clk, m_Std_Rover_Clk;
        /// \brief standard deviation for the base and rover trosphere delay
        double m_Std_Base_Trp, m_Std_Rover_Trp;
        /// \brief standard deviation for ambguity
        double m_Std_N[MAXSAT];

        /// \brief standard deviation of base station and rover station ionsphere delay (m)
        double m_Std_Base_Ion;
        double m_Std_Rover_Ion;

        /// \brief standard deviation of the velecity of base station and rover station (m/s)
        double m_Std_Base_VelX,m_Std_Base_VelY,m_Std_Base_VelZ;
        double m_Std_Rover_VelX,m_Std_Rover_VelY,m_Std_Rover_VelZ;

        /// \brief Stores the random number generator
        libPF::RandomNumberGenerationStrategy* m_RNG;

        /**
         * Private Method
         */
        double SQR(double x) { return x * x; }
        double SQRT(double x) { return x<0.0?0.0:sqrt(x);}
    };
}
#endif //ARC_ARC_MOVEMENTMODEL_H

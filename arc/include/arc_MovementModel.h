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
#include <libPF/CRandomNumberGenerator.h>
#include <libPF/RandomNumberGenerationStrategy.h>
#include "arc_States.h"
#include "arc_PF.h"
#include "arc_assert_macros.hpp"
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
    class ARC_MovementModel :
            public libPF::MovementModel<ARC_States> {
    public:
        ARC_DEFINE_EXCEPTION(Exception, std::runtime_error);
        /**
         * Constructor
         */
        ARC_MovementModel();
        ARC_MovementModel(const ARC_OPT* OPT,ARC_RTK* SRTK);
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
        inline void SetStdX(double Std,int Index){
            StdX[Index]=Std;
        }
        /// \brief get the standard deviation of index-th states
        inline double getStdX(int Index) const {
            return StdX[Index];
        }
        /// \brief set the navigation data
        inline void SetNav(const ARC_NAV* NAV){
            m_NAV=NAV;
        }
        /// \brief set the observation data
        inline void SetObs(const ARC_OBSD* OBS,int Nobs){
            m_OBS=OBS; m_Nobs=Nobs;
        }
        /// \brief set the arc-srtk solution data
        inline void SetSRTK(ARC_RTK* SRTK) {
            m_SRTK=SRTK;
        }
        /// \brief set the numbers of common satellite
        inline void SetComNum(int Num) {
            Ns=Num;
        }
        /// \brief set the base and rover station satellite list
        inline void SetComSatList(const int* rSat,const int* uSat,
                             const int* Sat,int Num) {
            Ns=Num;
            for (int i=0;i<Num;i++) m_RoverSat[i]=uSat[i],
                                    m_BaseSat[i]=rSat[i],SatList[i]=Sat[i];
        }
        inline void SetRoverStd(double std) {
            PF_ROVERPOS_STD=std;
        }
        inline double getRoverStd() {
            return PF_ROVERPOS_STD;
        }
        inline double getAmbMin() {
            return PF_AMB_MIN;
        }
        inline double getAmbMax() {
            return PF_AMB_MAX;
        }
        inline void SetAmbStd(double Min,double Max) {
            PF_AMB_MIN=Min;
            PF_AMB_MAX=Max;
        }
    private:
        /// \brief standard deviation for the states
        double StdX[MAXPFSTETAS];
        /// \brief Stores the random number generator
        libPF::RandomNumberGenerationStrategy* m_RNG;
        /// \brief arc-srtk solution type
        ARC_RTK* m_SRTK;
        /// \brief observation data
        const ARC_OBSD* m_OBS;
        /// \brief navigation data type
        const ARC_NAV* m_NAV;
        /// \brief observation data numbers
        int m_Nobs;
        /// \brief common sat numbers
        int Ns;
        /// \brief base and rover station satellite list
        int m_RoverSat[MAXSAT];
        int m_BaseSat[MAXSAT];
        /// \brief All satellite list
        int SatList[MAXSAT];

        /// \brief particle filter rover position std
        double PF_ROVERPOS_STD;
        /// \brief particle fiter ambguity min and max search
        double PF_AMB_MIN,PF_AMB_MAX;
    };
}
#endif //ARC_ARC_MOVEMENTMODEL_H

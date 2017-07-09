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
        ARC_States m_TrueArcState;
   
        /// \brief the navigation data,(x,y,z,vx,vy,vz){m,m/s}
        double m_SatPos[MAXSAT*6];
        /// \brief save all satelite information
        ARC_SAT m_SatInfo[MAXSAT];
        /// \brief navigation data
        const ARC_NAV* m_Nav;
        /// \brief observation data
        ARC_OBSD m_Obs[MAXSAT*2];
        /// \brief processing options type
        const ARC_OPT* m_Opt;
        /// \brief arc-srtk solution type
        ARC_RTK m_ARC;
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

        /// \brief undifferenced phase/code residuals
        double m_ZDRes[MAXSAT*2];

        /**
         * Methods
         */
        /**
         * \brief compute the undifferenced phase/code residuals
         */
        /// \brief undifferenced phase/code residuals
        int ComputeZd() const;
        /// \brief compute the satellite position
        void ComputeSatPos() const;
        /// \brief set the gnss observation and navgation data
        void SetObsData(const ARC_OBSD *Obs,int Nu,int Nb) ;
        void SetNavData(const ARC_NAV* Nav) ;
        /// \brief set the process settings
        void SetOpt(const ARC_OPT* Opt);
        /// \brief set the satellite information
        void SetSatInfo(const ARC_SAT* SatInfo);
        /// \brief time-interpolation of residuals (for post-mission)
        double IntPres(ARC_Time time, const ARC_OBSD *obs, int n, const ARC_NAV *nav,
                       double *y) const;
        /// \brief select common satellites between rover and reference station
        int SelectCommonSat(const ARC_OBSD *obs, double *azel, int nu, int nr,
                            const ARC_OPT *opt, int *sat, int *iu, int *ir) const;
        /// \brief precise tropspheric model
        double PrecTrop(ARC_Time time, const double *pos, int r,
                        const double *azel, const ARC_OPT *opt,
                        double *dtdx) const;
        /// \brief test navi system (m=0:gps/qzs/sbs,1:glo,2:gal,3:bds)
        inline int test_sys(int sys, int m) const {
            switch (sys) {
                case SYS_GPS: return m==0;
                case SYS_QZS: return m==0;
                case SYS_SBS: return m==0;
                case SYS_GLO: return m==1;
                case SYS_GAL: return m==2;
                case SYS_CMP: return m==3;
            }
            return 0;
        };
        /// \brief test valid observation data
        inline  int validobs(int i, int j, int f, int nf, double *y) const {
            return y[f+i*nf*2]!=0.0&&y[f+j*nf*2]!=0.0&&
                (f<nf||(y[f-nf+i*nf*2]!=0.0&&y[f-nf+j*nf*2]!=0.0));
        }
    };
}
#endif //ARC_ARC_OBSERVATIONMODEL_H

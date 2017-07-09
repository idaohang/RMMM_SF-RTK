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

#ifndef ARC_ARC_STATES_H
#define ARC_ARC_STATES_H
#include "arc.h"

/**
 * @class ARC_States
 * @brief the SRTK states for particel filter
 *
 * This state has the following parameters:
 * @li <b>xb</b>    the X-position of the base station (ECEF/m)
 * @li <b>yb</b>    the Y-Position of the base station (ECEF/m)
 * @li <b>zb</b>    the Z-Position of the base station (ECEF/m)
 * @li <b>xu</b>    the X-position of the rover station (ECEF/m)
 * @li <b>yu</b>    the Y-position of the rover station (ECEF/m)
 * @li <b>zu</b>    the Z-position of the rover station (ECEF/m)
 * @li <b>dtb</b>   the clock drift of the base station (m)
 * @li <b>dtu</b>   the clock drift of the rover station (m)
 * @li <b>Ni</b>    the ambiguity of GPS and BDS ,every satelite have a ambiguity (cycle)
 * @li <b>Tb</b>    the trosphere delay of the base station (m)
 * @li <b>Tu</b>    the trosphere delay of the rover station (m)
 * @li <b>Ib</b>    the ionosphere delay of base station (m)
 * @li <b>Iu</b>    the ionosphere delay of rover station (m)
 * @author sujinglan
 */
namespace ARC {
    class ARC_States {
    public:
        /**
         * Constructor
         */
        ARC_States();
        /**
         * Destructor
         */
        virtual ~ARC_States();
        /**
         * Opereator
         */
        ARC_States operator*(double factor) const;

        ARC_States &operator+=(const ARC_States &other);

        /// \brief set the position of base station in ECEF
        inline void setBasePosX(double BasePosX);

        inline void setBasePosY(double BasePosY);

        inline void setBasePosZ(double BasePosZ);

        /// \brief set the position of rover station in ECEF
        inline void setRoverPosX(double RoverX);

        inline void setRoverPosY(double RoverY);

        inline void setRoverPosZ(double RoverZ);

        /// \brief get the postion of base station in ECEF
        inline double getBasePosX() const;

        inline double getBasePosY() const;

        inline double getBasePosZ() const;

        /// \brief get the postion of rover station in ECEF
        inline double getRoverPosX() const;

        inline double getRoverPosY() const;

        inline double getRoverPosZ() const;

        /// \brief set the clock drift value of base and rover station (m)
        inline void setBaseClk(double BaseClk);

        inline void setRoverClk(double RoverClk);

        /// \brief get the clock drift value of base and rover station (m)
        inline double getRoverClk() const;

        inline double getBaseClk() const;

        /// \brief set the trosphere delay of base and rover station
        inline void setBaseTrop(double BaseTrop);

        inline void setRoverTrop(double RoverTrop);

        /// \brief get the trosphere delay of base and rover station
        inline double getBaseTrop() const;

        inline double getRoverTrop() const;

        /// \brief set the ambguity of every satelite for single-difference on base-rover station
        /// @param N[in] is the ambguity of satelite of single difference on base and rover station
        inline void setAmb(double *N);
        /// @param N[in] the ith satelite single difference ambguity
        inline void setAmb(double N,int i);

        /// \brief get the ambguity
        inline double *getAmb() const;
        /// \brief get the ith satelite single difference ambguity
        inline double getAmb(int sat) const { return m_N[sat-1];}

        /// \brief set the ionosphere delay of base station
        inline void setBaseIon(double BaseIon) {m_BaseIon=BaseIon;}
        /// \brief set the ionosphere delay of rover station
        inline void setRoverIon(double RoverIon) {m_RoverIon=RoverIon;}

        /// \brief get the ionosphere delay of base and rover station
        inline double getBaseIon(int sat) const { return m_BaseIon[sat-1];}
        inline double getRoverIon(int sat) const { return m_RoverIon[sat-1];}

        /// \brief get and set the velecity of base station and rover station
        inline void setBaseVelX(double Vel) {m_Base_VelX=Vel;}
        inline void setRoverVelX(double Vel) {m_Rover_VelX=Vel;}
        inline double getBaseVelX() const{ return m_Base_VelX;}
        inline double getRoverVelX() const { return m_Rover_VelX;}

        inline void setBaseVelY(double Vel) {m_Base_VelY=Vel;}
        inline void setRoverVelY(double Vel) {m_Rover_VelY=Vel;}
        inline double getBaseVelY() const{ return m_Base_VelY;}
        inline double getRoverVelY() const { return m_Rover_VelY;}

        inline void setBaseVelZ(double Vel) {m_Base_VelZ=Vel;}
        inline void setRoverVelZ(double Vel) {m_Rover_VelZ=Vel;}
        inline double getBaseVelZ() const{ return m_Base_VelZ;}
        inline double getRoverVelZ() const { return m_Rover_VelZ;}

        inline double getBaseTrpGN() const { return m_BaseTrpGN;}
        inline double getBaseTrpGE() const { return m_BaseTrpGE;}

        inline double getRoverTrpGN() const { return m_RoverTrpGN;}
        inline double getRoverTrpGE() const { return m_RoverTrpGE;}

        inline void setBaseTrpGN(double GN) {m_BaseTrpGN=GN;}
        inline void setBaseTrpGE(double GE) {m_BaseTrpGE=GE;}

        inline void setRoverTrpGN(double GN) {m_RoverTrpGN=GN;}
        inline void setRoverTrpGE(double GE) {m_RoverTrpGE=GE;}


    private:
        /// \brief the states for base station position
        double m_BaseX;
        double m_BaseY;
        double m_BaseZ;
        /// \brief the states for rover station postion
        double m_RoverX;
        double m_RoverY;
        double m_RoverZ;
        /// \brief the clock drift of base and rover station
        double m_BaseClk;
        double m_RoverClk;
        /// \brief the troposphere delay of base station
        double m_BaseTrop;
        double m_RoverTrop;
        /// \brief the ambiguity value of every satelite
        double m_N[MAXSAT];
        /// \brief the ionosphere delay of base station
        double m_BaseIon[MAXSAT];
        /// \brief the ionosphere delay of rover station
        double m_RoverIon[MAXSAT];

        /// \brief the velecity of base station
        double m_Base_VelX,m_Base_VelY,m_Base_VelZ;
        /// \brief the velecity of rover station
        double m_Rover_VelX,m_Rover_VelY,m_Rover_VelZ;

        /// \brief the gradient parameters of troposphere
        double m_BaseTrpGN,m_BaseTrpGE;
        double m_RoverTrpGn,m_RoverTrpGE;
    };
}
#endif //ARC_ARC_STATES_H

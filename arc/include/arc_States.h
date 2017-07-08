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
        void setBasePosX(double BasePosX);

        void setBasePosY(double BasePosY);

        void setBasePosZ(double BasePosZ);

        /// \brief set the position of rover station in ECEF
        void setRoverPosX(double RoverX);

        void setRoverPosY(double RoverY);

        void setRoverPosZ(double RoverZ);

        /// \brief get the postion of base station in ECEF
        double getBasePosX() const;

        double getBasePosY() const;

        double getBasePosZ() const;

        /// \brief get the postion of rover station in ECEF
        double getRoverPosX() const;

        double getRoverPosY() const;

        double getRoverPosZ() const;

        /// \brief set the clock drift value of base and rover station (m)
        void setBaseClk(double BaseClk);

        void setRoverClk(double RoverClk);

        /// \brief get the clock drift value of base and rover station (m)
        double getRoverClk() const;

        double getBaseClk() const;

        /// \brief set the trosphere delay of base and rover station
        void setBaseTrop(double BaseTrop);

        void setRoverTrop(double RoverTrop);

        /// \brief get the trosphere delay of base and rover station
        double getBaseTrop() const;

        double getRoverTrop() const;

        /// \brief set the ambguity of every satelite for single-difference on base-rover station
        void setAmb(double *N);

        /// \brief get the ambguity
        double *getAmb() const;

        /// \brief set the ionosphere delay of base station
        void setBaseIon(double BaseIon) {m_BaseIon=BaseIon;}
        /// \brief set the ionosphere delay of rover station
        void setRoverIon(double RoverIon) {m_RoverIon=RoverIon;}

        /// \brief get the ionosphere delay of base station
        double getBaseIon() const { return m_BaseIon;}
        double getRoverIon() const { return m_RoverIon;}

        /// \brief get and set the velecity of base station and rover station
        void setBaseVelX(double Vel) {m_Base_VelX=Vel;}
        void setRoverVelX(double Vel) {m_Rover_VelX=Vel;}
        double getBaseVelX() const{ return m_Base_VelX;}
        double getRoverVelX() const { return m_Rover_VelX;}

        void setBaseVelY(double Vel) {m_Base_VelY=Vel;}
        void setRoverVelY(double Vel) {m_Rover_VelY=Vel;}
        double getBaseVelY() const{ return m_Base_VelY;}
        double getRoverVelY() const { return m_Rover_VelY;}

        void setBaseVelZ(double Vel) {m_Base_VelZ=Vel;}
        void setRoverVelZ(double Vel) {m_Rover_VelZ=Vel;}
        double getBaseVelZ() const{ return m_Base_VelZ;}
        double getRoverVelZ() const { return m_Rover_VelZ;}

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
        /// \brief the trosphere delay of base station
        double m_BaseTrop;
        double m_RoverTrop;
        /// \brief the ambiguity value of every satelite
        double m_N[MAXSAT];
        /// \brief the ionosphere delay of base station
        double m_BaseIon;
        /// \brief the ionosphere delay of rover station
        double m_RoverIon;

        /// \brief the velecity of base station
        double m_Base_VelX,m_Base_VelY,m_Base_VelZ;
        /// \brief the velecity of rover station
        double m_Rover_Vel,m_Rover_VelY,m_Rover_VelZ;
    };
}
#endif //ARC_ARC_STATES_H

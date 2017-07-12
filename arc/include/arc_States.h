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
#include "arc_PF.h"

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
        ARC_States(const ARC_OPT* OPT);
        ARC_States(int Num);
        /**
         * Destructor
         */
        virtual ~ARC_States();
        /**
         * Opereator
         */
        ARC_States operator*(double factor) const{
            ARC_States State;
            for (int i=0;i<MAXPFSTETAS;i++) State.SetStatesValue(X[i]*factor,i);
            return State;
        };
        ARC_States &operator+=(const ARC_States &other){
            for (int i=0;i<MAXPFSTETAS;i++) X[i]=X[i]+other.getStateValue(i);
        };
        /// \brief set the value of states
        inline void SetStatesValue(const double Val,int Index) {
            X[Index]=Val;
        };
        /// \brief get the value of states
        inline double getStateValue(int Index) const{
            return X[Index];
        };
        /// \brief get the numbers of states
        inline int getStatesNum() const {
            return MAXPFSTETAS;
        };
        /// \brief get the states value
        inline double *getStatesVal() {
            return X;
        }
    private:
        /// \brief need to estimate states,this is same to rtklib
        mutable double X[MAXPFSTETAS];
    };
}
#endif //ARC_ARC_STATES_H

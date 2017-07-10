/*********************************************************************************
 *  ARC-SRTK - Single Frequency RTK Pisotioning Library
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
 *  Created on: July 09, 2017
 *      Author: SuJingLan
 *********************************************************************************/

/**
 * @file arc_StateDistribution.h
 * @brief This file contains states distribution.
 * @author SuJingLan
 */

#ifndef ARC_ARC_STATEDISTRIBUTION_H
#define ARC_ARC_STATEDISTRIBUTION_H

#include "arc_States.h"
#include <libPF/StateDistribution.h>

namespace libPF {
    class RandomNumberGenerationStrategy;
}

namespace ARC {
    class ARC_StateDistribution :
            public libPF::StateDistribution<ARC_States>
    {
    public:
        ARC_StateDistribution();
        ARC_StateDistribution(const ARC_OPT *OPT,const ARC_RTK* SRTK);
        ~ARC_StateDistribution();

        const ARC_States draw() const;
        void SetUniformofStates(double min,double max,int index);
        /// \brief setting of arc-srtk
        void SetOpt(const ARC_OPT* OPT);
        /// \brief set the arc-srtk solution data
        void SetSRTK(const ARC_RTK* SRTK) {
            m_SRTK=SRTK;
        }
    private:
        libPF::RandomNumberGenerationStrategy* m_RNG;
        /// \brief the arc-srtk solution settings
        const ARC_OPT* m_OPT;
        /// \brief the arc-srtk soliution type
        const ARC_RTK* m_SRTK;
        /// \brief the uniform distribution of states
        double *m_UniState;
    };
}
#endif //ARC_ARC_STATEDISTRIBUTION_H

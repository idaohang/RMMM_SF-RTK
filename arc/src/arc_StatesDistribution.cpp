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
 * @file arc_StateDistribution.cpp
 * @brief This source file of states distribution
 * @author SuJingLan
 */

#include <cmath>
#include <libPF/CRandomNumberGenerator.h>
#include <arc_StateDistribution.h>
namespace ARC {
    ARC_StateDistribution::ARC_StateDistribution() {
        m_RNG = new libPF::CRandomNumberGenerator();
    }
    ARC_StateDistribution::ARC_StateDistribution(const ARC_OPT *OPT,
                                                 const ARC_RTK *SRTK) {
        m_OPT=OPT; m_SRTK=SRTK;
        m_UniState=mat(2,m_SRTK->nx);
    }
    ARC_StateDistribution::~ARC_StateDistribution() {
        delete m_RNG;
        delete m_UniState;
    }

    const ARC_States ARC_StateDistribution::draw() const {
        ARC_States state;
        for (int i=0;i<m_SRTK->nx;i++) {
            state.SetStatesValue(m_RNG->getUniform(m_UniState[2*i],
                                                   m_UniState[2*i+1]),i);
        }
        return state;
    }
    void ARC_StateDistribution::SetUniformofStates(double min,
                                                   double max,
                                                   int index) {
        ARC_ASSERT_TRUE(Exception,m_SRTK->nx > 0,"Particle Filter States is Zero");
        if (index>=m_SRTK->nx)
            return;
        m_UniState[2*index  ]=min;
        m_UniState[2*index+1]=max;
    }
    void ARC_StateDistribution::SetOpt(const ARC_OPT *OPT) {
        if (OPT) m_OPT=OPT;
    }
}

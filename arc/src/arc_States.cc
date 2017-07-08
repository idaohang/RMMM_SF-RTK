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

/**
 * @file arc_States.cc
 * @brief Source file for the ARC_States class.
 * @author sujinglan
 */

#include "arc_States.h"
#include "arc_assert_macros.hpp"

/// \brief ARC Main namespace of this package.
namespace ARC {
    /// The default constructor.
    ARC_States::ARC_States() : m_BaseX(0.0), m_BaseY(0.0), m_BaseZ(0.0),
                               m_RoverX(0.0), m_RoverY(0.0), m_RoverZ(0.0),
                               m_RoverClk(0.0), m_BaseClk(0.0),
                               m_RoverTrop(0.0), m_BaseTrop(0.0) {
        for (size_t i = 0; i < MAXSAT; i++) m_N[i] = 0.0;
    }

    ARC_States::~ARC_States() {}

    /// set the base position in ECEF (m)
    void ARC_States::setBasePosX(double BasePosX) { m_BaseX = BasePosX; }

    void ARC_States::setBasePosY(double BasePosY) { m_BaseY = BasePosY; }

    void ARC_States::setBasePosZ(double BasePosZ) { m_BaseZ = BasePosZ; }

    /// set the rover position in ECEF (m)
    void ARC_States::setRoverPosX(double RoverPosX) { m_RoverX = RoverPosX; }

    void ARC_States::setRoverPosY(double RoverPosY) { m_RoverY = RoverPosY; }

    void ARC_States::setRoverPosZ(double RoverPosZ) { m_RoverZ = RoverPosZ; }

    /// set the base and rover trosphere delay (m)
    void ARC_States::setBaseTrop(double BaseTrop) { m_BaseTrop = BaseTrop; }

    void ARC_States::setRoverTrop(double RoverTrop) { m_RoverTrop = RoverTrop; }

    /// set the clock drift (m) of rover station and base station
    void ARC_States::setRoverClk(double RoverClk) { m_RoverClk = RoverClk; }

    void ARC_States::setBaseClk(double BaseTrop) { m_BaseTrop = BaseTrop; }

    /// set the ambguity of every satelite (cycle)
    void ARC_States::setAmb(double *N) {
        for (size_t i = 0; i < MAXSAT; i++) m_N[i] = N[i];
    }

    /// get the base postion in ECEF (m)
    double ARC_States::getBasePosX() const { return m_BaseX; }

    double ARC_States::getBasePosY() const { return m_BaseY; }

    double ARC_States::getBasePosZ() const { return m_BaseZ; }

    /// get the rover postion in ECEF (m)
    double ARC_States::getRoverPosX() const { return m_RoverX; }

    double ARC_States::getRoverPosY() const { return m_RoverY; }

    double ARC_States::getRoverPosZ() const { return m_RoverZ; }

    /// get the base and rover clock drift (m)
    double ARC_States::getBaseClk() const { return m_BaseClk; }

    double ARC_States::getRoverClk() const { return m_RoverClk; }

    /// get the single-difference ambguity of every satelite
    double *ARC_States::getAmb() const {
        double *Amb = new double[MAXSAT];
        for (size_t i = 0; i < MAXSAT; i++) Amb[i] = m_N[i];
        return Amb;
    }

    /// get the troephere delay of base and rover station
    double ARC_States::getBaseTrop() const { return m_RoverTrop; }

    double ARC_States::getRoverTrop() const { return m_BaseTrop; }

    ARC_States ARC_States::operator*(double factor) const {
        ARC_States states;
        states.m_RoverX = m_RoverX * factor;
        states.m_RoverY = m_RoverY * factor;
        states.m_RoverZ = m_RoverZ * factor;
        states.m_BaseX = m_BaseX * factor;
        states.m_BaseY = m_BaseY * factor;
        states.m_BaseZ = m_BaseZ * factor;
        states.m_RoverClk = m_RoverClk * factor;
        states.m_BaseClk = m_BaseClk * factor;
        states.m_RoverTrop = m_RoverTrop * factor;
        states.m_BaseTrop = m_BaseTrop * factor;

        for (size_t i = 0; i < MAXSAT; i++) states.m_N[i] = m_N[i] * factor;

        return states;
    }

    ARC_States &ARC_States::operator+=(const ARC_States &other) {
        m_BaseX = m_BaseX + other.m_BaseX;
        m_BaseY = m_BaseY + other.m_BaseY;
        m_BaseZ = m_BaseZ + other.m_BaseZ;
        m_RoverX = m_RoverX + other.m_RoverX;
        m_RoverY = m_RoverY + other.m_RoverY;
        m_RoverZ = m_RoverZ + other.m_RoverZ;
        m_RoverClk = m_RoverClk + other.m_RoverClk;
        m_BaseClk = m_BaseClk + other.m_BaseClk;
        m_BaseTrop = m_BaseTrop + other.m_BaseTrop;
        m_RoverTrop = m_RoverTrop + other.m_RoverTrop;
        for (size_t i = 0; i < MAXSAT; i++) m_N[i] = m_N[i] + other.m_N[i];
        return *this;
    }
}



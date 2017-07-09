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
 * @brief ARC-SRTK Particle Filter Some Const-Variable and Marco
 * @author sujinglan
 */
#ifndef ARC_ARC_PF_H
#define ARC_ARC_PF_H

#define ARC_PF_BASEPOS_STD                 (1.0)           /// initial standard deviation of base station position (m)
#define ARC_PF_ROVERPOS_STD                (100.0)         /// initial standard deviation of rover station position (m)
#define ARC_PF_BASECLK_STD                 (200.0)         /// initial standard deviation of base station clock drift (m)
#define ARC_PF_ROVERCLK_STD                (200.0)         /// initial standard deviation of rover station clock drift (m)
#define ARC_PF_BASETROP_STD                (5.0)           /// initial standard deviation of base station trosphere delay (m)
#define ARC_PF_ROVERTROP_STD               (5.0)           /// initial standard deviation of rover station trosphere delay (m)
#define ARC_PF_AMB_STD                     (10.0)          /// initial standard deviation of single-difference ambguity (cycle)
typedef nav_t ARC_NAV;                                     /// arc-srtk navigation data type
typedef obs_t ARC_OBS;                                     /// arc-srtk gnss observation data type
typedef ssat_t ARC_SAT;                                    /// arc-srtk satellite status type
typedef obsd_t ARC_OBSD;                                   /// arc-srtk gnss observation data record
typedef prcopt_t ARC_OPT;                                  /// arc-srtk processing options type
typedef gtime_t  ARC_Time;                                 /// arc-srtk observation time
typedef rtk_t ARC_RTK;                                     /// arc-srtk solution data type

#define NF(opt) ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf) /// how many frequency of gnss data



#endif //ARC_ARC_PF_H

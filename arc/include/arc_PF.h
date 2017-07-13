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

#define USE_PF_FILTER                                      /// the flag of using particle filter

#define ARC_PF_BASEPOS_STD                 (1.0)           /// initial standard deviation of base station position (m)
#define ARC_PF_ROVERPOS_STD                (2.0)           /// initial standard deviation of rover station position (m)
#define ARC_PF_BASECLK_STD                 (200.0)         /// initial standard deviation of base station clock drift (m)
#define ARC_PF_ROVERCLK_STD                (200.0)         /// initial standard deviation of rover station clock drift (m)
#define ARC_PF_BASETROP_STD                (5.0)           /// initial standard deviation of base station trosphere delay (m)
#define ARC_PF_ROVERTROP_STD               (5.0)           /// initial standard deviation of rover station trosphere delay (m)
#define ARC_PF_AMB_STD                     (0.5)           /// initial standard deviation of single-difference ambguity (cycle)
#define NF(opt) ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf) /// how many frequency of gnss data
#define ARC_PF_AMB_MIN                     (-2)
#define ARC_PF_AMB_MAX                     (+2)
#define ARC_PF_ITERS                       1
#define ARC_PF_NUM                         5000            /// particle numbers
#define ARC_PF_USE_SD                      0
#define ARC_PF_USE_DD                      1

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define VAR_POS     SQR(30.0)                              /// initial variance of receiver pos (m^2)
#define VAR_VEL     SQR(10.0)                              /// initial variance of receiver vel ((m/s)^2)
#define VAR_ACC     SQR(10.0)                              /// initial variance of receiver acc ((m/ss)^2)
#define VAR_HWBIAS  SQR(1.0)                               /// initial variance of h/w bias ((m/MHz)^2)
#define VAR_GRA     SQR(0.001)                             /// initial variance of gradient (m^2)
#define INIT_ZWD    0.15                                   /// initial zwd (m)

#define PRN_HWBIAS  1E-6                                   /// process noise of h/w bias (m/MHz/sqrt(s))
#define GAP_RESION  120                                    /// gap to reset ionosphere parameters (epochs)

#define VAR_HOLDAMB 0.001                                  /// constraint to hold ambiguity (cycle^2)

#define TTOL_MOVEB  (1.0+2*DTTOL)                          /// time sync tolerance for moving-baseline (s)

                                                           /// number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated)
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))
                                                           /// state variable index
#define II(s,opt)   (NP(opt)+(s)-1)                        /// ionos (s:satellite no)
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r))        /// tropos (r:0=rov,1:ref)
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))          /// receiver h/w bias
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)             /// phase bias (s:satno,f:freq)
#define MAXPFSTETAS                         80             /// max particle filter states
typedef nav_t                               ARC_NAV;       /// arc-srtk navigation data type
typedef obs_t                               ARC_OBS;       /// arc-srtk gnss observation data type
typedef ssat_t                              ARC_SAT;       /// arc-srtk satellite status type
typedef obsd_t                              ARC_OBSD;      /// arc-srtk gnss observation data record
typedef prcopt_t                            ARC_OPT;       /// arc-srtk processing options type
typedef gtime_t                             ARC_Time;      /// arc-srtk observation time
typedef rtk_t                               ARC_RTK;       /// arc-srtk solution data type

#endif //ARC_ARC_PF_H

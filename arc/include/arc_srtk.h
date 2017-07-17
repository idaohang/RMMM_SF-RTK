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
 *  Created on: July 13, 2017
 *      Author: SuJingLan
 *********************************************************************************/

/**
 * @file arc.h
 * @brief Header file for the ARC-SRTK Library
 * @author SuJingLan
 */

#ifndef ARC_SRTK_H
#define ARC_SRTK_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif
/* constants/global variables -------------------------------------------------*/
//#define ARC_TRACE_MAT                     /* matrix printf */
#define GLOG            1                 /* google log for debug */
#define ARC_NOLOG       -1                /* disable log informations */
#define ARC_INFO        0				  /* google information log */
#define ARC_WARNING 	1				  /* google warnings */
#define ARC_ERROR       2				  /* google errors */
#define ARC_FATAL       3				  /* google fatals */
#define ARC_LOGFILE     4 				  /* recore log information to file */
#define ARC_MATPRINTF   5                 /* matrix printf flag */
#define ARC_SOLVALTHRES 4.0               /* validation of solution thres */
/* single frequency rtk position post-processing ------------------------------*/
extern int  arc_postpos(gtime_t ts, gtime_t te, double ti, double tu,
                        const prcopt_t *popt, const solopt_t *sopt,
                        const filopt_t *fopt, char **infile, int n, char *outfile,
                        const char *rov, const char *base);
/* using particle filter solve single frequency position ---------------------*/
extern int  arc_pf_srtk(gtime_t ts, gtime_t te, double ti, double tu,
                        const prcopt_t *popt, const solopt_t *sopt,
                        const filopt_t *fopt, char **infile, int n, char *outfile);
/* single frequency rtk positioning ------------------------------------------*/
extern int  arc_srtkpos(rtk_t *rtk, const obsd_t *obs, int nobs, const nav_t *nav);
/* arc trace log functions ---------------------------------------------------*/
extern void arc_traceopen(const char *file);
extern void arc_traceclose(void);
extern void arc_tracelevel(int level);
extern void arc_log(int level, const char *format, ...);
extern void arc_tracet(int level, const char *format, ...);
extern void arc_tracemat(int level, const double *A, int n, int m, int p, int q);
extern void arc_traceobs(int level, const obsd_t *obs, int n);
extern void arc_tracenav(int level, const nav_t *nav);
extern void arc_tracegnav(int level, const nav_t *nav);
extern void arc_tracehnav(int level, const nav_t *nav);
extern void arc_tracepeph(int level, const nav_t *nav);
extern void arc_tracepclk(int level, const nav_t *nav);
extern void arc_traceb(int level, const unsigned char *p, int n);
extern void arc_tracebuf(int buffcount);
/* arc cholesky functions ------------------------------------------------------*/
extern double *arc_cholesky(double *A,int n);

#ifdef __cplusplus
}
#endif
#endif
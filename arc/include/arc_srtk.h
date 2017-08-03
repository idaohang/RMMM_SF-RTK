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
 * @file arc_srtk.h
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
#define ARC_TRACE_MAT                     /* matrix printf */
#define GLOG            1                 /* google log for debug */
#define ARC_LOGTOFILE   1                 /* google log output file */
#define ARC_NOLOG       -1                /* disable log informations */
#define ARC_INFO        0				  /* google information log */
#define ARC_WARNING 	1				  /* google warnings */
#define ARC_ERROR       2				  /* google errors */
#define ARC_FATAL       3				  /* google fatals */
#define ARC_LOGFILE     4 				  /* recore log information to file */
#define ARC_MATPRINTF   5                 /* matrix printf flag */
#define ARC_SOLVALTHRES 4.0               /* validation of solution thres */
#define ARC_CERES_SINGLE      1           /* ceres solver single epoch solution */
#define ARC_CERES_WINDOWS     2           /* ceres solver windows solution */
#define ARC_TEST                          /* arc-srtk test define */
#define ARC_UKF_USEPNT_INIT               /* using standard position to initial ukf prior states and its covariance matrix */
#define ARC_AMB_PART    1
#define ARC_AMB_NOPART  0
#define ARC_UKF         1
#define ARC_NOUKF       0
/* single frequency rtk position post-processing ------------------------------*/
extern int  arc_postpos(gtime_t ts, gtime_t te, double ti, double tu,
                        prcopt_t *popt, const solopt_t *sopt,
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
extern void arc_tracemat(int level,const double *A,int n,int m,int p,int q);
extern void arc_tracemati(int level,const int *A,int n,int m,int p,int q);
extern void arc_traceobs(int level,const obsd_t *obs,int n);
extern void arc_tracenav(int level,const nav_t *nav);
extern void arc_tracegnav(int level,const nav_t *nav);
extern void arc_tracehnav(int level,const nav_t *nav);
extern void arc_tracepeph(int level,const nav_t *nav);
extern void arc_tracepclk(int level,const nav_t *nav);
extern void arc_traceb(int level,const unsigned char *p,int n);
extern void arc_tracebuf(int buffcount);
extern void arc_set_glog_tofile(int opt);
/* arc cholesky functions ------------------------------------------------------*/
extern double *arc_cholesky(double *A,int n);
/* computethe trace of matrix---------------------------------------------------*/
extern double arc_mattrace(double *A,int n);
/* adaptive Kaman filter -------------------------------------------------------*/
extern int adap_kaman_filter(rtk_t* rtk,double *x, double *P, const double *H,
                             const double *v,const double *R,int n,int m);
/* tropmodel model--------------------------------------------------------------*/
extern double arc_tropmodel_hopf(gtime_t time, const double *pos, const double *azel,
                                 double humi);
extern double arc_tropmodel_unb3(gtime_t time, const double *pos, const double *azel,
                                 double humi,double *zhd,double *zwd);
extern double arc_tropmodel_mops(gtime_t time, const double *pos, const double *azel,
                                 double humi,double *zhd,double *zwd);
extern double arc_tropmodel_gcat(gtime_t time, const double *pos, const double *azel,
                                 double humi);
extern double arc_tropmodel_black(gtime_t time, const double *pos, const double *azel,
                                  double humi,double *zhd,double *zwd);
extern double arc_tropmodel_waas(gtime_t time, const double *pos, const double *azel,
                                 double humi);
/* troposphere mapping function-------------------------------------------------*/
extern double arc_tropmapf_UNSW931(gtime_t time, const double *pos, const double *azel,
                                   double *mapfw);
extern double arc_tropmapf_exp(gtime_t time, const double *pos, const double *azel,
                               double *mapfw);
extern double arc_tropmapf_ifadis(gtime_t time, const double *pos, const double *azel,
                                  double *mapfw);
extern double arc_tropmapf_ma_mu(gtime_t time, const double *pos, const double *azel,
                                 double *mapfw);
extern double arc_tropmapf_mtt(gtime_t time, const double *pos, const double *azel,
                               double *mapfw);
extern double arc_tropmapf_chao(gtime_t time, const double *pos, const double *azel,
                                double *mapfw);
extern double arc_tropmapf_cfa2_2(gtime_t time, const double *pos, const double *azel,
                                  double *mapfw);
/* ukf filter -------------------------------------------------------------------*/
extern ukf_t* arc_ukf_filter_new(unsigned int state_dim,
                                 unsigned int measure_dim,
                                 double *Q,
                                 double *R,
                                 filter_function ffun,
                                 measure_function mfun);
extern void arc_ukf_filter_delete(ukf_t *filter);
extern void arc_ukf_filter_compute_weights(ukf_t *filter,double alpha,double k,double beta);
extern void arc_ukf_filter_reset(ukf_t *filter, double *x0,double *PO);
extern void arc_ukf_filter_get_state(ukf_t *filter, double *x, double *P);
extern int  arc_ukf_filter_update(ukf_t *filter, double *y, double *u,double*F,
                                  double *G);
extern void arc_ukf_free_problem();

/* satellites, systems, codes functions --------------------------------------*/
extern int  satno   (int sys, int prn);
extern int  satsys  (int sat, int *prn);
extern int  satid2no(const char *id);
extern void satno2id(int sat, char *id);
extern unsigned char obs2code(const char *obs, int *freq);
extern char *code2obs(unsigned char code, int *freq);
extern int  satexclude(int sat, int svh, const prcopt_t *opt);
extern int  testsnr(int base, int freq, double el, double snr,
                    const snrmask_t *mask);
extern void setcodepri(int sys, int freq, const char *pri);
extern int  getcodepri(int sys, unsigned char code, const char *opt);
extern unsigned int getbitu(const unsigned char *buff, int pos, int len);
extern int execcmd(const char *cmd);
extern void setbitu(unsigned char *buff, int pos, int len, unsigned int data);
extern double arc_snr_varerr(const double dbhz,int f,int nf,const prcopt_t* opt);

/* matrix and vector functions -----------------------------------------------*/
extern double *arc_mat(int n, int m);
extern int    *arc_imat(int n, int m);
extern double *arc_zeros(int n, int m);
extern double *arc_eye(int n);
extern double arc_dot(const double *a, const double *b, int n);
extern double arc_norm(const double *a, int n);
extern void arc_cross3(const double *a, const double *b, double *c);
extern int  arc_normv3(const double *a, double *b);
extern void arc_matcpy(double *A, const double *B, int n, int m);
extern void arc_matmul(const char *tr, int n, int k, int m, double alpha,
                       const double *A, const double *B, double beta, double *C);
extern int  arc_matinv(double *A, int n);
extern int  arc_solve(const char *tr, const double *A, const double *Y, int n,
                      int m, double *X);
extern int  arc_lsq(const double *A, const double *y, int n, int m, double *x,
                    double *Q);
extern int  arc_filter(double *x, double *P, const double *H, const double *v,
                       const double *R, int n, int m,double * D);
extern int  arc_smoother(const double *xf, const double *Qf, const double *xb,
                         const double *Qb, int n, double *xs, double *Qs);
extern void arc_add_fatal(fatalfunc_t *func);

/* time and string functions -------------------------------------------------*/
extern double  str2num(const char *s, int i, int n);
extern int     str2time(const char *s, int i, int n, gtime_t *t);
extern void    time2str(gtime_t t, char *str, int n);
extern gtime_t epoch2time(const double *ep);
extern void    time2epoch(gtime_t t, double *ep);
extern gtime_t gpst2time(int week, double sec);
extern double  time2gpst(gtime_t t, int *week);
extern gtime_t gst2time(int week, double sec);
extern double  time2gst(gtime_t t, int *week);
extern gtime_t bdt2time(int week, double sec);
extern double  time2bdt(gtime_t t, int *week);
extern char    *time_str(gtime_t t, int n);

extern gtime_t timeadd  (gtime_t t, double sec);
extern double  timediff (gtime_t t1, gtime_t t2);
extern gtime_t gpst2utc (gtime_t t);
extern gtime_t utc2gpst (gtime_t t);
extern gtime_t gpst2bdt (gtime_t t);
extern gtime_t bdt2gpst (gtime_t t);
extern gtime_t timeget  (void);
extern void    timeset  (gtime_t t);
extern double  time2doy (gtime_t t);
extern double  utc2gmst (gtime_t t, double ut1_utc);
extern int read_leaps(const char *file);

extern int adjgpsweek(int week);
extern unsigned int tickget(void);
extern void sleepms(int ms);

extern int reppath(const char *path, char *rpath, gtime_t time, const char *rov,
                   const char *base);
extern int reppaths(const char *path, char *rpaths[], int nmax, gtime_t ts,
                    gtime_t te, const char *rov, const char *base);

/* coordinates transformation ------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos);
extern void pos2ecef(const double *pos, double *r);
extern void ecef2enu(const double *pos, const double *r, double *e);
extern void enu2ecef(const double *pos, const double *e, double *r);
extern void covenu  (const double *pos, const double *P, double *Q);
extern void covecef (const double *pos, const double *Q, double *P);
extern void xyz2enu (const double *pos, double *E);
extern void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);
extern void deg2dms (double deg, double *dms, int ndec);
extern double dms2deg(const double *dms);

/* input and output functions ------------------------------------------------*/
extern void readpos(const char *file, const char *rcv, double *pos);
extern int  sortobs(obs_t *obs);
extern void uniqnav(nav_t *nav);
extern int  screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
extern int  readnav(const char *file, nav_t *nav);
extern int  savenav(const char *file, const nav_t *nav);
extern void freeobs(obs_t *obs);
extern void freenav(nav_t *nav, int opt);
extern int  readblq(const char *file, const char *sta, double *odisp);
extern int  readerp(const char *file, erp_t *erp);
extern int  geterp (const erp_t *erp, gtime_t time, double *val);

/* platform dependent functions ----------------------------------------------*/
extern int expath (const char *path, char *paths[], int nmax);

/* positioning models --------------------------------------------------------*/
extern double arc_satwavelen(int sat, int frq, const nav_t *nav);
extern double arc_satazel(const double *pos, const double *e, double *azel);
extern double arc_geodist(const double *rs, const double *rr, double *e);
extern void arc_dops(int ns, const double *azel, double elmin, double *dop);
extern void arc_csmooth(obs_t *obs, int ns);

/* atmosphere models ---------------------------------------------------------*/
extern double arc_ionmodel(gtime_t t, const double *ion, const double *pos,
                           const double *azel);
extern double arc_ionmapf(const double *pos, const double *azel);
extern double arc_ionppp(const double *pos, const double *azel, double re,
                         double hion, double *pppos);
extern double arc_tropmodel(gtime_t time, const double *pos, const double *azel,
                            double humi);
extern double arc_tropmapf(gtime_t time, const double *pos, const double *azel,
                           double *mapfw);
extern int arc_ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                        const double *azel, int ionoopt, double *ion, double *var);
extern int arc_tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                        const double *azel, int tropopt, double *trp, double *var);

/* antenna models ------------------------------------------------------------*/
extern int  arc_readpcv(const char *file, pcvs_t *pcvs);
extern pcv_t *arc_searchpcv(int sat, const char *type, gtime_t time,
                            const pcvs_t *pcvs);
extern void arc_antmodel(const pcv_t *pcv, const double *del, const double *azel,
                         int opt, double *dant);
extern void arc_antmodel_s(const pcv_t *pcv, double nadir, double *dant);

/* earth tide models ---------------------------------------------------------*/
extern void arc_sunmoonpos(gtime_t tutc, const double *erpv, double *rsun,
                           double *rmoon, double *gmst);
extern void arc_tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
                         const double *odisp, double *dr);

/* rinex functions -----------------------------------------------------------*/
extern int arc_readrnx(const char *file, int rcv, const char *opt, obs_t *obs,
                       nav_t *nav, sta_t *sta);
extern int arc_readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te,
                        double tint, const char *opt, obs_t *obs, nav_t *nav,
                        sta_t *sta);
extern int arc_readrnxc(const char *file, nav_t *nav);
extern int arc_rtk_uncompress(const char *file, char *uncfile);

/* ephemeris and clock functions ---------------------------------------------*/
extern double arc_eph2clk(gtime_t time, const eph_t *eph);
extern double arc_geph2clk(gtime_t time, const geph_t *geph);
extern double arc_seph2clk(gtime_t time, const seph_t *seph);
extern void arc_eph2pos(gtime_t time, const eph_t *eph, double *rs, double *dts,
                        double *var);
extern void arc_geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts,
                         double *var);
extern void arc_seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts,
                         double *var);
extern int  arc_peph2pos(gtime_t time, int sat, const nav_t *nav, int opt,
                         double *rs, double *dts, double *var);
extern void arc_satantoff(gtime_t time, const double *rs, int sat, const nav_t *nav,
                          double *dant);
extern int  arc_satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                       const nav_t *nav, double *rs, double *dts, double *var,
                       int *svh);
extern void arc_satposs(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                        int sateph, double *rs, double *dts, double *var, int *svh);
extern void arc_readsp3(const char *file, nav_t *nav, int opt);
extern int  arc_readsap(const char *file, gtime_t time, nav_t *nav);
extern int  arc_readdcb(const char *file, nav_t *nav, const sta_t *sta);
extern int  arc_readfcb(const char *file, nav_t *nav);
extern void arc_alm2pos(gtime_t time, const alm_t *alm, double *rs, double *dts);
extern int arc_pppnx(const prcopt_t *opt);

/* integer ambiguity resolution ----------------------------------------------*/
extern int arc_lambda(int n, int m, const double *a, const double *Q, double *F,
                      double *s,double *D,double *L);
extern int lambda_reduction(int n, const double *Q, double *Z);
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s);

extern int arc_par_lambda(const double *a,const double *Qa,int n,int m,double *F,
                          double *s,double p0);

extern int arc_bootstrap(int n,const double *a, const double *Q, double *F,double *Ps);
extern double arc_invRatio(int namb,double ffailure,double pf);
extern double arc_amb_bs_success(const double *D,int n);

/* standard positioning ------------------------------------------------------*/
extern int arc_pntpos(const obsd_t *obs, int n, const nav_t *nav,
                      const prcopt_t *opt, sol_t *sol, double *azel,
                      ssat_t *ssat, char *msg);
/* precise positioning -------------------------------------------------------*/
extern void arc_rtkinit(rtk_t *rtk, const prcopt_t *opt);
extern void arc_rtkfree(rtk_t *rtk);
extern double arc_conffunc(int N, double B, double sig);

/* application defined functions ---------------------------------------------*/
extern int arc_showmsg(char *format, ...);
extern void arc_settspan(gtime_t ts, gtime_t te);
extern void arc_settime(gtime_t time);
extern void arc_exclude_bds_geo(prcopt_t *opt);
extern int arc_is_bds_geo(int sat);
extern double arc_chi2(int na,double x,double *f);
extern double arc_re_chi2(int na,double p);
extern double arc_norm_distri(const double u);
extern double arc_re_norm(double p);
extern int arc_mateigenvalue(const double* A,int n,double *u,double *v);
extern int arc_matdet(const double*A,int n,double*det);
extern double arc_normcdf(double X);

#ifdef __cplusplus
}
#endif
#endif
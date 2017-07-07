//
// Created by sujinglan on 7/7/17.
//

#include "arc.h"
#if GLOG
#include "glog/logging.h"
#endif

static int arc_raim_fde(const obsd_t *obs, int n, int iter,const double *rs,
                 const double *dts, const double *vare, const int *svh,
                 const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                 double *azel, int *vsat, double *resp,char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e = { { 0 } };
    char tstr[32], name[16], msg_e[128]="";
    double *rs_e, *dts_e, *vare_e, *azel_e, *resp_e, rms_e, rms = 0.0;
    int i, j, k, nvsat, stat = 0, *svh_e, *vsat_e, sat = 0;

    trace(ARC_INFO, "raim_fde: %s n=%2d\n", time_str(obs[0].time, 0), n);

    if (!(obs_e = (obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6, n); dts_e = mat(2, n); vare_e = mat(1, n); azel_e = zeros(2, n);
    svh_e = imat(1, n); vsat_e = imat(1, n); resp_e = mat(1, n);

    for (i = nvsat = 0; i < n; i++)
    {
        if (!vsat[i])continue;
        rms += resp[i] * resp[i];
        nvsat++;
    }
    rms = sqrt(rms / nvsat - 1);

    double maxresp = 0;
    int m = 0;
    for (i = 0; i < n ; i++)
    {
        if (fabs(resp[i])>fabs(maxresp) && vsat[i] == 1)
        {
            maxresp = resp[i];
            m = i;
        }
    }
    vsat[m] = 0;//invalid satellite

    for (j = k = 0; j < n - iter; j++)
    {
        if (j == m) continue;
        obs_e[k] = obs[j];
        matcpy(rs_e + 6 * k, rs + 6 * j, 6, 1);
        matcpy(dts_e + 2 * k, dts + 2 * j, 2, 1);
        vare_e[k] = vare[j];
        vsat_e[k] = vsat[j];
        svh_e[k++] = svh[j];
    }
    /* estimate receiver position without a satellite */
    if (stat = estpos(obs_e, n - iter, rs_e, dts_e,
                      vare_e, svh_e, nav, opt, &sol_e, azel_e,
                      vsat_e, resp_e, msg_e)) {
        trace(ARC_WARNING, "raim_fde: exsat=%2d\n", obs[m].sat);
        //return 0;
    }
    for (j = nvsat = 0, rms_e = 0.0; j < n - iter; j++) {
        if (!vsat_e[j]) continue;
        rms_e += resp_e[j] * resp_e[j];
        nvsat++;
    }
    if (nvsat<5) {
        trace(ARC_ERROR, "raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
              obs[m].sat, nvsat);
        return 0;
    }
    rms_e = sqrt(rms_e / nvsat - 1);

    trace(ARC_WARNING, "raim_fde: exsat=%2d rms=%8.3f\n", obs[m].sat, rms_e);

    if (rms_e>rms) return 0;

    /* save result */
    for (j = k = 0; j < n - iter; j++) {
        if (j == m) continue;
        matcpy(azel + 2 * j, azel_e + 2 * k, 2, 1);
        vsat[j] = vsat_e[k];
        resp[j] = resp_e[k++];
    }
    *sol = sol_e;
    sat = obs[m].sat;
    rms = rms_e;
    strcpy(msg, msg_e);

    time2str(obs[0].time, tstr, 2); satno2id(sat, name);
    trace(ARC_WARNING, "%s: %s excluded by raim\n", tstr + 11, name);

    free(obs_e);
    free(rs_e); free(dts_e); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* check distance on pseduorange -------------------------------------------*/
static int arc_checkdist(obsd_t *obs,int n)
{
    double maxdist = 0, mindist = 1.5e7;
    for (int i = 0; i < n; i++)
    {
        if (satsys(obs[i].sat, NULL) == SYS_GPS) maxdist = 3.5e7;
        else if (satsys(obs[i].sat, NULL) == SYS_CMP) maxdist = 5.5e7;

        if (obs[i].P[0]<mindist || obs[i].P[0]>maxdist)
        {
            obs[i].Pv[0] = 0;
            obs[i].P[0] = 0.0;
            trace(ARC_INFO, "pseduorange invalid %s sat=%2d\n",
                  time_str(obs[i].time, 0),obs[i].sat);
        }
    }
    return 1;
}
/* detect pseduorange ------------------------------------------------------*/
static int arc_detepseduo(rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav)
{
    prcopt_t opt = rtk->opt;
    opt.ionoopt = IONOOPT_BRDC;
    opt.tropopt = TROPOPT_SAAS;
    sol_t sol = rtk->sol;
    double *rs, *dts, *var, *azel_, *resp;
    int i, stat, vsat[MAXOBS], svh[MAXOBS];
    char msg[128] = "";
    for (int j = 0; j < MAXOBS; j++){ vsat[j] = 1; }
    sol.stat = SOLQ_NONE;
    trace(ARC_INFO, "arc_detepseduo : tobs=%s n=%d\n",
          time_str(obs[0].time, 3), n);
    sol.time = obs[0].time;
    rs = mat(6, n);
    dts = mat(2, n); var = mat(1, n);
    azel_ = zeros(2, n); resp = mat(1, n);
    satposs(sol.time, obs, n, nav, opt.sateph, rs, dts, var, svh);

    stat = estpos(obs, n, rs, dts, var, svh,
                  nav, &opt, &sol, azel_, vsat, resp, msg);

    for (i = 1; i < 10; i++) {
        if (!stat && (n - i) >= 5) {
            stat = arc_raim_fde(obs, n,i, rs, dts,
                                var, svh, nav, &opt, &sol,
                                azel_, vsat, resp,msg);
        }
        if (stat != 0)break;
    }
    for (i = 0; i < n; i++) {
        obs[i].Pv[0] = vsat[i];
        if (obs[i].Pv[0] == 0)  {
            obs[i].P[0] = 0.0;
            trace(ARC_INFO, "pseduorange invalid %s sat=%2d\n",
                  time_str(obs[i].time, 0), obs[i].sat);
        }

    }
    return 1;
}
/* gnss(gps and bds) rinex format data preprocess ----------------------------
 * gps and bds rinex format data preprocess
 * args   : rtk_t   *rtk       I  rtk_t strcut data
 *          obsd_t  *obs       IO input gnss data ,and save the preprocess gnss data
 *          int      n         I  numbers of obsdatas ( include rovers and base)
 *          nav_t   *nav       I  navigation struct (const)
 * return : status (0 : ok , others : errors)
 * -------------------------------------------------------------------------*/
extern int arc_predata(rtk_t *rtk,obsd_t *obs, int n, const nav_t *nav)
{
    int nu, nr;
    int i;
    if (n<=0) return 0;
    for (nu = 0; nu < n&&obs[nu].rcv == 1; nu++);
    for (nr = 0; nu + nr < n&&obs[nu + nr].rcv == 2; nr++);
    for (i = 0; i < n; i++) {
        obs[i].Pv[0] = 1;
        obs[i].Lv[0] = 1;
    }
    arc_checkdist(obs, n);
    arc_detepseduo(rtk, obs, nu, nav);
    arc_detepseduo(rtk, obs + nu, nr, nav);
    return 1;
}



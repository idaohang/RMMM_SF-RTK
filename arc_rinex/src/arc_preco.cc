//
// Created by sujinglan on 7/7/17.
//

#include "arc.h"
#if GLOG
#include "glog/logging.h"
#endif

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



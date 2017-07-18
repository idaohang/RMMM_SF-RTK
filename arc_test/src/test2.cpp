
#include "arc.h"
#define NUMINFILE 10

int main()
{
    gtime_t ts = { 0 };
    gtime_t te = { 0 };
    double ti = 0.0, tu = 0.0;
    char *infile[NUMINFILE];
    char *outfile = "";
    char *base = "", *rover = "";

    prcopt_t arc_opt = prcopt_default;
    solopt_t arc_solopt = solopt_default;

    arc_opt.mode = PMODE_STATIC;
    arc_opt.ionoopt = IONOOPT_SBAS;
    arc_opt.tropopt = TROPOPT_SAAS;
    arc_opt.modear = ARMODE_FIXHOLD;                    //AR mode(0:off, 1 : continuous, 2 : instantaneous（瞬时的）, 3 : fix and hold, 4 : ppp - ar) * /
    arc_opt.bdsmodear = 1;
    arc_opt.elmaskar = 15.0* D2R;                       // elevation mask of AR for rising satellite (deg)
    arc_opt.elmaskhold = 20.0* D2R;                     // elevation mask to hold ambiguity (deg)
    arc_opt.niter = 1;                                  // 滤波迭代次数
    arc_opt.refpos = 0;                                 // 0:pos in prcopt
                                                        // (0:pos in prcopt,  1:average of single pos,
                                                        //  2:read from file, 3:rinex header, 4:rtcm pos)
    arc_opt.rb[0] =  -2267819.9911;
    arc_opt.rb[1] = 5009410.5991;
    arc_opt.rb[2] = 3220957.0217;
    arc_opt.refpos=0;

    arc_opt.nf = 1;                                     // 解算频率
    arc_opt.elmin = 15.0 * D2R;
    arc_opt.navsys = SYS_CMP|SYS_GPS;
    arc_opt.dynamics=0;
    arc_opt.ceres=0;
    arc_opt.ceres_cholesky=1;

    arc_solopt.posf = SOLF_XYZ;
    arc_solopt.outopt = 1;                              // output processing options (0:no,1:yes) */
    arc_solopt.sstat = 0;                               // solution statistics level (0:off,1:states,2:residuals) */
    arc_solopt.trace = 1;

    filopt_t arc_fopt = fileopt_default;

    const char *roverobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/03/GPSC1Log201608230510_3.16o";
    const char *baseobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/base/5202K81320201608230000B.16o";
    const char *navf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/03/GPSC1Log201608230510_3.16n";
    const char *navf1 = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/03/GPSC1Log201608230510_3.16c";
    const char *navf2 = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/base/5202K81320201608230000B.16c";
    const char *navf3 = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/base/5202K81320201608230000B.16n";

    for (int i = 0; i < NUMINFILE; i++) {
        if (!(infile[i] = (char *)malloc(1024))) {
            for (; i >= 0; i--) free(infile[i]);
            return -1;
        }
    }
    strcpy(infile[0], roverobsf);
    strcpy(infile[1], baseobsf);
    strcpy(infile[2], navf);
    strcpy(infile[3],navf1);
    strcpy(infile[4],navf2);
    strcpy(infile[5],navf3);

    arc_tracelevel(ARC_NOLOG);
    arc_tracebuf(10);

    arc_postpos(ts, te, ti, tu, &arc_opt, &arc_solopt, &arc_fopt, infile, 6, outfile,rover,base);

    return 0;
}


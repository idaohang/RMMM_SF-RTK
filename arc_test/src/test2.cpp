
#include "arc.h"

using namespace std;

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

    arc_opt.mode    = PMODE_STATIC;                       // -静态模式或者动态模式
    arc_opt.ionoopt = IONOOPT_OFF;                        // -电离层延迟改正
    arc_opt.tropopt = TROPOPT_SAAS;                       // -对流层延迟改正
    arc_opt.modear  = ARMODE_FIXHOLD;                     // -模糊度固定后的处理策略：
                                                          //  ARMODE_FIXHOLD表示继承上一个历元的权值信息，ARMODE_INST表示每一个历元都初始化模糊度
    arc_opt.bdsmodear = 1;
    arc_opt.elmaskar = 10.0* D2R;
    arc_opt.elmaskhold = 20.0* D2R;
    arc_opt.niter = 1;
    arc_opt.refpos = 0;
    // (0:pos in prcopt,  1:average of single pos,
    //  2:read from file, 3:rinex header, 4:rtcm pos)
    arc_opt.rb[0] = -2405143.8990;
    arc_opt.rb[1] = 5385195.2549;
    arc_opt.rb[2] = 2420032.5440;

    //arc_opt.rb[0] = -2392060.4126;
    //arc_opt.rb[1] = 5383990.1051;
    //arc_opt.rb[2] = 2435512.0331;

    //arc_opt.rb[0] =  -2267819.9911;
    //arc_opt.rb[1] = 5009410.5991;
    //arc_opt.rb[2] = 3220957.0217;

    arc_opt.refpos=0;

    arc_opt.nf=1;
    arc_opt.elmin=15.0*D2R;
    arc_opt.navsys=SYS_GPS;
    arc_opt.dynamics=0;                                   // 1表示用动力学方程，0表示不使用动力学方程
    arc_opt.ceres=0;                                      // 线性优化开关
    arc_opt.ceres_cholesky=0;
    arc_opt.ukf=0;                                        // ukf非线性滤波开关
    arc_opt.ukf_alpha=1E-4;
    arc_opt.ukf_beta=4;
    arc_opt.ukf_ZCount=0;
    arc_opt.reset_amb_all=0;                              // 1表示每一个历元重置所有模糊度
    arc_opt.amb_part_var=0;                               // 根据模糊度的误差值，进行部分模糊度固定
    arc_opt.amb_iter=3;                                   // 部分模糊度的迭代次数
    arc_opt.amb_ref_thres=0.99;                           // 参考模糊度取整的判断条件
    arc_opt.exclude_bds_geo=0;                            // 是否去除GEO卫星
    arc_opt.amb_group=0;                                  // 分组模糊度固定开关，效果不好，已经删除该函数
    arc_opt.amb_el_group=45.0*D2R;
    arc_opt.init_dc=0;                                    // 用伪距差分定位对每一个历元初始化

    arc_opt.amb_fix_mode=AMBFIX_LAMBDA;                   // 模糊度固定策略
    arc_opt.amb_boostps=0.99;                             // boostrap模糊度固定时的成功率限值
    arc_opt.amb_partps=0.95;                              // 基于成功率的部分模糊度固定的参数

    arc_opt.kalman_robust=0;                              // 抗野值卡儿曼滤波开关
    arc_opt.kalman_robust_alpha=0.001;

    arc_opt.snr_det=0;                                    // 信噪比异常值探测
    arc_opt.snr_alpha=0.01;

    arc_opt.amb_part_D=0;                                 // 基于LDL分解的部分模糊度固定

    arc_opt.amb_delay=1;                                  // 变换参考卫星时的延迟历元
    arc_opt.amb_ref_delayc=5;

    arc_opt.dynamics_dc=1;                                // 伪距差分动力学模型开关

    arc_opt.inst_amb=0;

    arc_opt.no_amb_sol=0;

    arc_solopt.posf = SOLF_XYZ;
    arc_solopt.outopt = 1;
    arc_solopt.sstat = 0;
    arc_solopt.trace = 1;

    arc_opt.vel_prn[0]=arc_opt.vel_prn[1]=arc_opt.vel_prn[2]=0.001;
    arc_opt.est_doppler=0;
    arc_opt.clk_dri_prn=0.5;

    filopt_t arc_fopt = fileopt_default;

    const char *roverobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hklt/hklt1910.17o";
    const char *baseobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hkkt/hkkt1910.17o";
    const char *navf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hklt/hklt1910.17n";
    const char *navf1 = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hkkt/hkkt1910.17n";
    const char *navf2="/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hkkt/hkkt1910.17f";

    //const char *roverobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/static/CLP10160.16o";
    //const char *baseobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/static/JZD10160.16o";
    //const char *navf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/static/brdm0160.16p";

    //const char *roverobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/03/GPSC1Log201608230510_3.16o";
    //const char *baseobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/base/5202K81320201608230000B.16o";
    //const char *navf  = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/03/GPSC1Log201608230510_3.16n";
    //const char *navf1 = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/03/GPSC1Log201608230510_3.16c";
    //const char *navf2 = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/base/5202K81320201608230000B.16c";
    //const char *navf3 = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/base/5202K81320201608230000B.16n";

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
    //strcpy(infile[5],navf3);

    arc_tracelevel(ARC_INFO);
    arc_tracebuf(0);
    arc_set_glog_tofile(0);

    arc_postpos(ts,te,ti,tu,&arc_opt,&arc_solopt,&arc_fopt,infile,5,outfile,rover,base);

    return 0;
}


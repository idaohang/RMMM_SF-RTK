
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
 * @file arc_core.cc
 * @brief Source file for the ARC-SRTK Library
 * @author SuJingLan
 */

#include "arc.h"

/* constants/global variables --------------------------------------------------*/
#define MAXFILEPATH 1024                    /* max length of file path */
#define NUMINFILE   20                      /* max numbers of input file */

#define ARC_SRTK_INFO                                                              "\n\n\
======================================ARC-SRTK======================================\n\n\
    -r   rover station observation file \n\
    -b   base station observation file \n\
    -n   navigation file \n\
    -d   bds navigation file \n\
    -g   gps navigation file\n\
    -c   base station position {x,y,z} \n\
    -m   position mode {static/single/dgps/kinema/move} \n\
    -i   option file \n\
    -t   trace level {info/warnings/errors/fatals} \n\
    -y   dynamic mode {off/on} \n\
    -u   ukf filter {off/on} \n\
    -s   robust kalman filter {off/on} \n\
    -e   snr outlier detect {off/on} \n\
    -p   using difference-pseudorange to initial if it sets to 1, \n\
         ohterwise using standard positioning (only active in none-dynamic mode) \n\
    -h   help information \n\
    -a   ambiguity fix partial {off/on} \n\
    -l   instantaneous ar mode,different from rtklib {off/on} \n\
    -q   estimate rover station clock drift rate {off/on} \n\n\
====================================================================================\n\n"

static char roevr_obsfile[MAXFILEPATH]=""; /* rover station observation file path */
static char base_obsfile[MAXFILEPATH]="";  /* base station observation file path */
static char navfile[MAXFILEPATH]="";       /* navigation file path */
static char outfile[MAXFILEPATH]="";       /* arc-srtk solution output file path */
static double basepos[3]={0};              /* base station position */
static char bdsnavfile[MAXFILEPATH]="";    /* bds navigation file */
static char gpsnavfile[MAXFILEPATH]="";    /* gps navigation file */
static char inifile[MAXFILEPATH]="";       /* option file */
static char staticfile[MAXFILEPATH]="";    /* solution static file */

static int dynamic=0;                       /* dynamic mode */
static int ukf=0;                           /* ukf filter mode */
static int mode=PMODE_STATIC;               /* position mode */
static int tracelevel=ARC_NOLOG;              /* trace level */
static int robust=0;                        /* robust kalman fiter */
static int snr_detect=0;                    /* snr outlier detect */
static int init_dc=0;                       /* using difference-pseudorange to initial */
static int dynamics_dc=0;                   /* difference-pseudorange position on dynamic mode */
static int ambpart=0;                       /* fix partial ambiguity */
static int ambinst=0;                       /* instantaneous ar mode,different from rtklib */
static int estclk=0;                        /* estimate rover station clock drift rate */

int ExecCmd(const QString &cmd, int show)
{
    Q_UNUSED(show);
    return QProcess::startDetached(cmd);    /* FIXME: show option not yet supported */
}

/* ARC-SRTK option----------------------------------------------------------- */
static const char *arc_opt="r:b:n:d:g:c:m:i:t:y:u:s:e:p:ha:l:q:";
/* options for arc-srtk-------------------------------------------------------

   -r   rover station observation file
   -b   base station observation file
   -n   navigation file
   -d   bds navigation file
   -g   gps navigation file
   -c   base station position {x,y,z}
   -m   position mode {static/single/dgps/kinema/move}
   -i   option file
   -t   trace level {off/info/warnings/errors/fatals}
   -y   dynamic mode {off/on}
   -u   ukf filter {off/on}
   -s   robust kalman filter {off/on}
   -e   snr outlier detect {off/on}
   -p   using difference-pseudorange to initial if it sets to 1,
        ohterwise using standard positioning (only active in none-dynamic mode)
   -h   help information
   -a   ambiguity fix partial {off/on}
   -l   instantaneous ar mode,different from rtklib {off/on}
   -q   estimate rover station clock drift rate {off/on}
--------------------------------------------------------------------------------*/

/* ARC-SRTK MAIN FUNCTION-------------------------------------------------------*/
int main(int argc,char *argv[])
{
    int i,j=0,ch;

    gtime_t ts={0};
    gtime_t te={0};
    double ti=0.0,tu=0.0;
    char *infile[MAXFILEPATH],*p=NULL;

    opterr=1;

    while((ch=getopt(argc,argv,arc_opt))!= -1) {
        switch(ch) {
            case 'r' :
                if (optarg) strcpy(roevr_obsfile,optarg);
                break;
            case 'b' :
                if (optarg) strcpy(base_obsfile,optarg);
                break;
            case 'n' :
                if (optarg) strcpy(navfile,optarg);
                break;
            case 'd' :
                if (optarg) strcpy(bdsnavfile,optarg);
                break;
            case 'c':
                if (optarg) {
                    p=strtok(optarg,","); basepos[0]=atof(p);
                    p=strtok(NULL,",");   basepos[1]=atof(p);
                    p=strtok(NULL,",");   basepos[2]=atof(p);
                }
                break;
            case 'g':
                if (optarg) strcpy(gpsnavfile,optarg);
                break;
            case 'o' :
                if (optarg) strcpy(navfile,outfile);
            case 'm' :
                if (optarg) {
                    if      (strcmp(optarg,"static")==0)  mode=PMODE_STATIC;
                    else if (strcmp(optarg,"kinamic")==0) mode=PMODE_KINEMA;
                    else if (strcmp(optarg,"dgps")==0)    mode=PMODE_DGPS;
                    else if (strcmp(optarg,"move")==0)    mode=PMODE_MOVEB;
                    else if (strcmp(optarg,"single")==0)  mode=PMODE_SINGLE;
                    else                                  mode=PMODE_STATIC;
                }
            case 'i':
                if (optarg) {
                    strcpy(inifile,optarg);
                }
                break;
            case 't':
                if      (!strcmp(optarg,"info"))     tracelevel=ARC_INFO;
                else if (!strcmp(optarg,"warnings")) tracelevel=ARC_WARNING;
                else if (!strcmp(optarg,"errors"))   tracelevel=ARC_WARNING;
                else if (!strcmp(optarg,"fatals"))   tracelevel=ARC_FATAL;
                else if (!strcmp(optarg,"fatals"))   tracelevel=ARC_NOLOG;
                break;
            case 'y':
                if      (!strcmp(optarg,"off")) dynamic=0;
                else if (!strcmp(optarg,"on"))  dynamic=1;
                break;
            case 'u':
                if      (!strcmp(optarg,"off")) ukf=0;
                else if (!strcmp(optarg,"on"))  ukf=1;
                break;
            case 's':
                if      (!strcmp(optarg,"off")) robust=0;
                else if (!strcmp(optarg,"on"))  robust=1;
                break;
            case 'e':
                if      (!strcmp(optarg,"off")) snr_detect=0;
                else if (!strcmp(optarg,"on"))  snr_detect=1;
                break;
            case 'p':
                if      (!strcmp(optarg,"off")) init_dc=0;
                else if (!strcmp(optarg,"on"))  init_dc=1;
                break;
            case 'h':
                fprintf(stderr,"%s \n",ARC_SRTK_INFO);
                return 0;
            case 'a':
                if      (!strcmp(optarg,"off")) ambpart=0;
                else if (!strcmp(optarg,"on"))  ambpart=1;
                break;
            case 'l':
                if      (!strcmp(optarg,"off")) ambinst=0;
                else if (!strcmp(optarg,"on"))  ambinst=1;
                break;
            case 'q':
                if      (!strcmp(optarg,"off")) estclk=0;
                else if (!strcmp(optarg,"on"))  estclk=1;
                break;
            default:
                break;
        }
    }

    prcopt_t arc_opt   =prcopt_default;
    solopt_t arc_solopt=solopt_default;
    filopt_t arc_fopt  =fileopt_default;

    arc_opt.mode   =mode<0?PMODE_STATIC:mode;
    arc_opt.ionoopt=IONOOPT_OFF;
    arc_opt.tropopt=TROPOPT_SAAS;
    arc_opt.modear =ARMODE_FIXHOLD;

    arc_matcpy(arc_opt.rb,basepos,1,3);

    arc_opt.bdsmodear =1;
    arc_opt.elmaskar  =10.0*D2R;
    arc_opt.elmaskhold=20.0*D2R;
    arc_opt.niter =1;
    arc_opt.refpos=0;
    arc_opt.nf=1;
    arc_opt.elmin=15.0*D2R;
    arc_opt.navsys=SYS_GPS;
    arc_opt.dynamics=dynamic;
    arc_opt.kalman_robust=robust;

    arc_opt.ukf=ukf;
    arc_opt.snr_det=snr_detect;

    arc_opt.init_dc=init_dc;
    arc_opt.inst_amb=ambinst;
    arc_opt.amb_part_var=ambpart;
    arc_opt.est_doppler=estclk;

    arc_tracelevel(tracelevel);
    arc_tracebuf(0);
    arc_set_glog_tofile(0);

    if (strcmp(inifile,"")!=0) {
        resetsysopts();
        if (!loadopts(inifile,sysopts)) {
            arc_info(0,2,"ARC-SRTK: No Options File: %s. Defaults Used");
        }
        getsysopts(&arc_opt,&arc_solopt,&arc_fopt);
        strcpy(roevr_obsfile,arc_fopt.roverobs);
        strcpy(base_obsfile, arc_fopt.baseobs);
        strcpy(navfile,      arc_fopt.navfile);
        strcpy(bdsnavfile,   arc_fopt.bdsnav);
        strcpy(gpsnavfile,   arc_fopt.gpsnav);
    }
    if (!strcmp(roevr_obsfile,"")||!strcmp(base_obsfile,"")||!strcmp(navfile,"")) {
        arc_info(0,2,"ARC-SRTK: Pleas Input Rover and Base Station Observation Files");
        return 0;
    }
    for (i=0;i<NUMINFILE;i++) {
        if (!(infile[i]=(char *)malloc(1024))) {
            for (;i>=0;i--) free(infile[i]); return -1;
        }
    }
    strcpy(infile[0],roevr_obsfile);
    strcpy(infile[1],base_obsfile);
    strcpy(infile[2],navfile);
    strcpy(outfile,arc_fopt.output);
    strcpy(staticfile,outfile);

    if (strcmp(bdsnavfile,"")!=0) strcpy(infile[3+(j++)],bdsnavfile);
    if (strcmp(gpsnavfile,"")!=0) strcpy(infile[3+(j++)],gpsnavfile);

    arc_info(0,3,"ARC-SRTK Start Processing");

    arc_postpos(ts,te,ti,tu,&arc_opt,
                &arc_solopt,&arc_fopt,infile,3+j,outfile,"","");

    for (i=0;i<NUMINFILE;i++) free(infile[i]);

    if (!strcmp(bdsnavfile,"")) strcpy(bdsnavfile,"arc-srtk");
    if (!strcmp(gpsnavfile,"")) strcpy(gpsnavfile,"arc-srtk");
    if (!strcmp(outfile,""))    strcpy(outfile,"arc-srtk");
    if (!strcmp(roevr_obsfile,"")) strcpy(roevr_obsfile,"arc-srtk");
    if (!strcmp(navfile,""))       strcpy(navfile,"arc-srtk");

    QString QRoverObs=QString(QLatin1String(roevr_obsfile));
    QString QNavp=QString(QLatin1String(navfile));
    QString QNavn=QString(QLatin1String(gpsnavfile));
    QString QNavc=QString(QLatin1String(bdsnavfile));
    QString QSol=QString(QLatin1String(outfile));

    QString cmd="../arc_plot/bin/plot",opts="";
    QString cmd1="./plot";
    opts=" -r \""+QRoverObs+"\" \""+QNavp+"\" \""+QNavn+"\" \""+ QNavc+"\" \""+QSol;

    if (!ExecCmd(cmd+opts,1)&&!ExecCmd(cmd1+opts,1)) {
        arc_info(100,2,"ARC-SRTK: Plot Error");
    }
    return 0;
}































































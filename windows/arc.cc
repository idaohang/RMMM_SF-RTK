// arc.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "arc.h"

#define NUMINFILE 10

int main()
{
	gtime_t ts={0};
	gtime_t te={0};
	double ti=0.0,tu=0.0;
	char *infile[NUMINFILE];
	char *outfile="";
	char *base="",*rover="";

	prcopt_t arc_opt=prcopt_default;
	solopt_t arc_solopt=solopt_default;

	filopt_t arc_fopt=fileopt_default;

	const char *roverobsf="C:\\Users\\sujinglan\\Desktop\\ins_gps\\imu_data_pos320\\03\\GPSC1Log201608230510_3.16o";
	const char *baseobsf="C:\\Users\\sujinglan\\Desktop\\ins_gps\\imu_data_pos320\\base\\5202K81320201608230000B.16o";

	const char *rover_gps_navf="C:\\Users\\sujinglan\\Desktop\\ins_gps\\imu_data_pos320\\03\\GPSC1Log201608230510_3.16n";
	const char *base_gps_navf="C:\\Users\\sujinglan\\Desktop\\ins_gps\\imu_data_pos320\\base\\5202K81320201608230000B.16n";

	const char *rover_bds_navf="C:\\Users\\sujinglan\\Desktop\\ins_gps\\imu_data_pos320\\03\\GPSC1Log201608230510_3.16c";
	const char *base_bds_navf="C:\\Users\\sujinglan\\Desktop\\ins_gps\\imu_data_pos320\\base\\5202K81320201608230000B.16c";

	for (int i=0;i<NUMINFILE;i++) {
	    if (!(infile[i]=(char *)malloc(1024))) {
            for (;i>=0;i--) free(infile[i]);
            return -1;
        }    
	}

	strcpy(infile[0],roverobsf);
	strcpy(infile[1],baseobsf);
	strcpy(infile[2],rover_gps_navf);
	strcpy(infile[3],base_gps_navf);
	strcpy(infile[4],rover_bds_navf);
	strcpy(infile[5],base_bds_navf);

	tracelevel(ARC_NOLOG);
	tracebuf(10);

	postpos_(ts,te,ti,tu,&arc_opt,&arc_solopt,&arc_fopt,infile,6,outfile,base,rover);

    return 0;
}


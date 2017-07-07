/***************************************************************************
 * full_test.cpp is part of Math Graphic Library
 * Copyright (C) 2007-2016 Alexey Balakin <mathgl.abalakin@gmail.ru>       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include <time.h>
#include <locale.h>
#include <time.h>
#include <getopt.h>
#include <unistd.h>

#include "mgl2/mgl.h"
#include "mgl2/font.h"
#include "mgl2/eval.h"
//-----------------------------------------------------------------------------
void MGL_EXPORT mgl_create_cpp_font(HMGL gr, const wchar_t *how);
long MGL_EXPORT mgl_check_tex_table();
//-----------------------------------------------------------------------------
struct mglSample	/// Structure for list of samples
{
	const char *name;
	void (*func)(mglGraph*);
	const char *mgl;
};
extern mglSample samp[];
extern const char *mmgl_dat_prepare;
//-----------------------------------------------------------------------------
int MGL_LOCAL_PURE mgl_cmd_smp(const void *a, const void *b)
{
	const mglSample *aa = (const mglSample *)a;
	const mglSample *bb = (const mglSample *)b;
	return strcmp(aa->name, bb->name);
}
//-----------------------------------------------------------------------------
int type = 0;
int dotest  = 0;
int width  = 800;
int height = 600;
int big  = 0;
int srnd = 0;
int use_mgl = 0;
int verbose = 0;
int quality  = MGL_DRAW_NORM;
//-----------------------------------------------------------------------------
void mgls_prepare1d(mglData *y, mglData *y1=0, mglData *y2=0, mglData *x1=0, mglData *x2=0);
void mgls_prepare2d(mglData *a, mglData *b=0, mglData *v=0);
void mgls_prepare3d(mglData *a, mglData *b=0);
void mgls_prepare2v(mglData *a, mglData *b);
void mgls_prepare3v(mglData *ex, mglData *ey, mglData *ez);
//-----------------------------------------------------------------------------
void save(mglGraph *gr,const char *name,const char *suf);
void test(mglGraph *gr)
{
	mglParse par;
	par.Execute(gr,"text 0 -0.2 'abcde'\ntext 0 0 'abcde'[2]\n"
	"text 0 0.2 'abcde'+2\ntext 0 0.4 'abcde',2\ntext 0 0.6 'abcde',2,'k'"
	"text 0 0.8 'abcde',2,'k'+1\ntext 0 1 'abcde',2,'k'[5]\ntext 0 1 'abcde''k'[5]");
	return;

	gr->SubPlot(2,2,0);	gr->Axis();
	gr->SubPlot(2,2,1);	gr->Rotate(40,60);	gr->Axis();
	gr->SetRotatedText(false);	gr->SubPlot(2,2,2);	gr->Axis();
	gr->SubPlot(2,2,3);	gr->Rotate(40,60);	gr->Axis();	return;
	mglData a;	a.SetList(5,0.,1.,0.,1.,-1.,2.);
	gr->Plot(a);
	return;
	par.Execute(gr,"call 'test' -1\n func 'test' 1\nline $1 0 1 1 'b'\nreturn\n");
//	par.Execute(gr,"load '/home/balakin/mathgl-code/mathgl-2x/build/examples/libmgl_module.so':baxis\n");
//	par.Execute(gr,"subplot 1 1 0:#rotate 40 60\nperspective 1.22:box:axis\n");
	return;
}
//-----------------------------------------------------------------------------
static struct option longopts[] =
{
	{ "big",	no_argument,	&big,		1 },
	{ "web",	no_argument,	&big,		2 },
	{ "mini",	no_argument,	&big,		3 },
	{ "help",	no_argument,	NULL,		'?' },
	{ "height",	required_argument,	NULL,	'h' },
	{ "kind",	required_argument,	NULL,	'k' },
	{ "list",	no_argument,	NULL,		'l' },
	{ "mgl",	no_argument,	&use_mgl,	1 },
	{ "srnd",	no_argument,	&srnd,		1 },

	{ "png",	no_argument,	&type,		0 },
	{ "eps",	no_argument,	&type,		1 },
	{ "svg",	no_argument,	&type,		2 },
	{ "solid",	no_argument,	&type,		3 },
	{ "jpeg",	no_argument,	&type,		4 },
	{ "prc",	no_argument,	&type,		5 },
	{ "gif",	no_argument,	&type,		6 },
	{ "none",	no_argument,	&type,		7 },
	{ "bps",	no_argument,	&type,		8 },
	{ "pdf",	no_argument,	&type,		9 },
	{ "obj_old",no_argument,	&type,		10 },
	{ "obj",	no_argument,	&type,		11 },
	{ "off",	no_argument,	&type,		12 },
	{ "stl",	no_argument,	&type,		13 },
	{ "tex",	no_argument,	&type,		14 },
	{ "json",	no_argument,	&type,		15 },
	{ "jsonz",	no_argument,	&type,		16 },
	{ "docs",	no_argument,	&type,		17 },

	{ "test",	no_argument,	&dotest,	1 },
	{ "font",	no_argument,	&dotest,	2 },
	{ "time",	no_argument,	&dotest,	3 },
	{ "fexport",no_argument,	&dotest,	4 },
	{ "textbl",	no_argument,	&dotest,	5 },

	{ "thread",	required_argument,	NULL,	't' },
	{ "verbose",no_argument,	&verbose,	1 },
	{ "width",	required_argument,	NULL,	'w' },
	{ "quality",required_argument,	NULL,	'q' },
	{ NULL,		0,				NULL,		0 }
};
//-----------------------------------------------------------------------------
void usage()
{
	puts (
//		"--png		- output png\n"
		"--width=num	- png picture width\n"
		"--height=num	- png picture height\n"
		"--mini		- png picture is 200x150\n"
		"--big		- png picture is 1920x1440\n"
		"--web		- png picture is 640x480\n"
		"--prc		- output prc\n"
		"--pdf		- output pdf\n"
		"--eps		- output EPS\n"
		"--tex		- output LaTeX\n"
		"--jpeg		- output JPEG\n"
		"--json		- output JSON\n"
		"--jsonz		- output JSONz\n"
		"--solid	- output solid PNG\n"
		"--svg		- output SVG\n"
		"--obj		- output obj/mtl\n"
		"--obj_old	- output obj/mtl in old way\n"
		"--off		- output off\n"
		"--stl		- output stl\n"
		"--docs		- output in png, prc/pdf and json\n"
		"--none		- none output\n"
		"--srnd		- use the same random numbers in any run\n"
		"--list		- print list of sample names\n"
		"--kind=name	- produce only this sample\n"
		"--thread=num	- number of threads used\n"
		"--mgl		- use MGL scripts for samples\n"
		"--test		- run in test mode\n"
		"--time		- measure execution time for all samples\n"
		"--font		- write current font as C++ file\n"
		"--quality=val	- use specified quality for plot(s)\n"
		"--fexport	- test most of output formats\n"
	);
}
//-----------------------------------------------------------------------------
void save(mglGraph *gr,const char *name,const char *suf="")
{
	//	return;
	char buf[128];
	switch(type)
	{
		case 1:	// EPS
			snprintf(buf,128,"%s%s.eps",name,suf);
			gr->WriteEPS(buf);
			break;
		case 2:	// SVG
			snprintf(buf,128,"%s%s.svg",name,suf);
			gr->WriteSVG(buf);	break;
		case 3:	// PNG
			snprintf(buf,128,"%s%s.png",name,suf);
			gr->WritePNG(buf,0,true);	break;
		case 4:	// JPEG
			snprintf(buf,128,"%s%s.jpg",name,suf);
			gr->WriteJPEG(buf);	break;
		case 5:	// PRC
			snprintf(buf,128,"%s%s.prc",name,suf);
			gr->WritePRC(buf,"",false);	break;
		case 6:	// GIF
			snprintf(buf,128,"%s%s.gif",name,suf);
			gr->WriteGIF(buf);	break;
		case 7:	gr->Finish();	// none
			break;
		case 8:	// EPS to PNG
			snprintf(buf,128,"%s%s.png",name,suf);
			gr->WritePNG(buf,0,false);
			break;
 		case 9:	// PDF
			snprintf(buf,128,"%s%s.prc",name,suf);
			gr->WritePRC(buf);	remove(buf);	break;
		case 10:	// old OBJ
			snprintf(buf,128,"%s%s.obj",name,suf);
			gr->WriteOBJold(buf);	break;
		case 11:	// OBJ
			snprintf(buf,128,"%s%s.obj",name,suf);
			gr->WriteOBJ(buf);	break;
		case 12:	// OFF
			snprintf(buf,128,"%s%s.off",name,suf);
			gr->WriteOFF(buf);	break;
		case 13:	// STL
			snprintf(buf,128,"%s%s.stl",name,suf);
			gr->WriteSTL(buf);	break;
		case 14:	// TeX
			snprintf(buf,128,"%s%s.tex",name,suf);
			gr->WriteTEX(buf);	break;
		case 15:	// JSON
			snprintf(buf,128,"%s%s.json",name,suf);
			gr->WriteJSON(buf);	break;
		case 16:	// JSON
			snprintf(buf,128,"%s%s.jsonz",name,suf);
			gr->WriteJSON(buf,"",true);	break;
		case 17:	// PNG + JSON + PDF
			snprintf(buf,128,"%s%s.png",name,suf);
			gr->WritePNG(buf,0,true);
			snprintf(buf,128,"%s%s.json",name,suf);
			gr->WriteJSON(buf);
			gr->SetSize(height,height,false);
			snprintf(buf,128,"%s%s.prc",name,suf);
			gr->WritePRC(buf);	remove(buf);	break;
		default:// PNG (no alpha)
#if MGL_HAVE_PNG
			snprintf(buf,128,"%s%s.png",name,suf);
			gr->WritePNG(buf,0,false);	break;
#else
			snprintf(buf,128,"%s%s.bmp",name,suf);
			gr->WriteBMP(buf);	break;
#endif
	}
}
//-----------------------------------------------------------------------------
void smgl_fexport(mglGraph *gr)	// test file export
{
	gr->SubPlot(3,2,0);
	double d,x1,x2,x0,y=0.95;
	d=0.3, x0=0.2, x1=0.5, x2=0.6;
	gr->Line(mglPoint(x0,1-0*d),mglPoint(x1,1-0*d),"k-");	gr->Puts(mglPoint(x2,y-0*d),"Solid '-'",":rL");
	gr->Line(mglPoint(x0,1-1*d),mglPoint(x1,1-1*d),"k|");	gr->Puts(mglPoint(x2,y-1*d),"Long Dash '|'",":rL");
	gr->Line(mglPoint(x0,1-2*d),mglPoint(x1,1-2*d),"k;");	gr->Puts(mglPoint(x2,y-2*d),"Dash ';'",":rL");
	gr->Line(mglPoint(x0,1-3*d),mglPoint(x1,1-3*d),"k=");	gr->Puts(mglPoint(x2,y-3*d),"Small dash '='",":rL");
	gr->Line(mglPoint(x0,1-4*d),mglPoint(x1,1-4*d),"kj");	gr->Puts(mglPoint(x2,y-4*d),"Dash-dot 'j'",":rL");
	gr->Line(mglPoint(x0,1-5*d),mglPoint(x1,1-5*d),"ki");	gr->Puts(mglPoint(x2,y-5*d),"Small dash-dot 'i'",":rL");
	gr->Line(mglPoint(x0,1-6*d),mglPoint(x1,1-6*d),"k:");	gr->Puts(mglPoint(x2,y-6*d),"Dots ':'",":rL");
	gr->Line(mglPoint(x0,1-7*d),mglPoint(x1,1-7*d),"k ");	gr->Puts(mglPoint(x2,y-7*d),"None ' '",":rL");

	d=0.25; x1=-1; x0=-0.8;	y = -0.05;
	gr->Mark(mglPoint(x1,5*d),"k.");	gr->Puts(mglPoint(x0,y+5*d),"'.'",":rL");
	gr->Mark(mglPoint(x1,4*d),"k+");	gr->Puts(mglPoint(x0,y+4*d),"'+'",":rL");
	gr->Mark(mglPoint(x1,3*d),"kx");	gr->Puts(mglPoint(x0,y+3*d),"'x'",":rL");
	gr->Mark(mglPoint(x1,2*d),"k*");	gr->Puts(mglPoint(x0,y+2*d),"'*'",":rL");
	gr->Mark(mglPoint(x1,d),"ks");		gr->Puts(mglPoint(x0,y+d),"'s'",":rL");
	gr->Mark(mglPoint(x1,0),"kd");		gr->Puts(mglPoint(x0,y),"'d'",":rL");
	gr->Mark(mglPoint(x1,-d,0),"ko");	gr->Puts(mglPoint(x0,y-d),"'o'",":rL");
	gr->Mark(mglPoint(x1,-2*d,0),"k^");	gr->Puts(mglPoint(x0,y-2*d),"'\\^'",":rL");
	gr->Mark(mglPoint(x1,-3*d,0),"kv");	gr->Puts(mglPoint(x0,y-3*d),"'v'",":rL");
	gr->Mark(mglPoint(x1,-4*d,0),"k<");	gr->Puts(mglPoint(x0,y-4*d),"'<'",":rL");
	gr->Mark(mglPoint(x1,-5*d,0),"k>");	gr->Puts(mglPoint(x0,y-5*d),"'>'",":rL");

	d=0.25; x1=-0.5; x0=-0.3;	y = -0.05;
	gr->Mark(mglPoint(x1,5*d),"k#.");	gr->Puts(mglPoint(x0,y+5*d),"'\\#.'",":rL");
	gr->Mark(mglPoint(x1,4*d),"k#+");	gr->Puts(mglPoint(x0,y+4*d),"'\\#+'",":rL");
	gr->Mark(mglPoint(x1,3*d),"k#x");	gr->Puts(mglPoint(x0,y+3*d),"'\\#x'",":rL");
	gr->Mark(mglPoint(x1,2*d),"k#*");	gr->Puts(mglPoint(x0,y+2*d),"'\\#*'",":rL");
	gr->Mark(mglPoint(x1,d),"k#s");		gr->Puts(mglPoint(x0,y+d),"'\\#s'",":rL");
	gr->Mark(mglPoint(x1,0),"k#d");		gr->Puts(mglPoint(x0,y),"'\\#d'",":rL");
	gr->Mark(mglPoint(x1,-d,0),"k#o");	gr->Puts(mglPoint(x0,y-d),"'\\#o'",":rL");
	gr->Mark(mglPoint(x1,-2*d,0),"k#^");	gr->Puts(mglPoint(x0,y-2*d),"'\\#\\^'",":rL");
	gr->Mark(mglPoint(x1,-3*d,0),"k#v");	gr->Puts(mglPoint(x0,y-3*d),"'\\#v'",":rL");
	gr->Mark(mglPoint(x1,-4*d,0),"k#<");	gr->Puts(mglPoint(x0,y-4*d),"'\\#<'",":rL");
	gr->Mark(mglPoint(x1,-5*d,0),"k#>");	gr->Puts(mglPoint(x0,y-5*d),"'\\#>'",":rL");

	gr->SubPlot(3,2,1);
	double a=0.1,b=0.4,c=0.5;
	gr->Line(mglPoint(a,1),mglPoint(b,1),"k-A");		gr->Puts(mglPoint(c,1),"Style 'A' or 'A\\_'",":rL");
	gr->Line(mglPoint(a,0.8),mglPoint(b,0.8),"k-V");	gr->Puts(mglPoint(c,0.8),"Style 'V' or 'V\\_'",":rL");
	gr->Line(mglPoint(a,0.6),mglPoint(b,0.6),"k-K");	gr->Puts(mglPoint(c,0.6),"Style 'K' or 'K\\_'",":rL");
	gr->Line(mglPoint(a,0.4),mglPoint(b,0.4),"k-I");	gr->Puts(mglPoint(c,0.4),"Style 'I' or 'I\\_'",":rL");
	gr->Line(mglPoint(a,0.2),mglPoint(b,0.2),"k-D");	gr->Puts(mglPoint(c,0.2),"Style 'D' or 'D\\_'",":rL");
	gr->Line(mglPoint(a,0),mglPoint(b,0),"k-S");		gr->Puts(mglPoint(c,0),"Style 'S' or 'S\\_'",":rL");
	gr->Line(mglPoint(a,-0.2),mglPoint(b,-0.2),"k-O");	gr->Puts(mglPoint(c,-0.2),"Style 'O' or 'O\\_'",":rL");
	gr->Line(mglPoint(a,-0.4),mglPoint(b,-0.4),"k-T");	gr->Puts(mglPoint(c,-0.4),"Style 'T' or 'T\\_'",":rL");
	gr->Line(mglPoint(a,-0.6),mglPoint(b,-0.6),"k-_");	gr->Puts(mglPoint(c,-0.6),"Style '\\_' or none",":rL");
	gr->Line(mglPoint(a,-0.8),mglPoint(b,-0.8),"k-AS");	gr->Puts(mglPoint(c,-0.8),"Style 'AS'",":rL");
	gr->Line(mglPoint(a,-1),mglPoint(b,-1),"k-_A");		gr->Puts(mglPoint(c,-1),"Style '\\_A'",":rL");

	a=-1;	b=-0.7;	c=-0.6;
	gr->Line(mglPoint(a,1),mglPoint(b,1),"kAA");		gr->Puts(mglPoint(c,1),"Style 'AA'",":rL");
	gr->Line(mglPoint(a,0.8),mglPoint(b,0.8),"kVV");	gr->Puts(mglPoint(c,0.8),"Style 'VV'",":rL");
	gr->Line(mglPoint(a,0.6),mglPoint(b,0.6),"kKK");	gr->Puts(mglPoint(c,0.6),"Style 'KK'",":rL");
	gr->Line(mglPoint(a,0.4),mglPoint(b,0.4),"kII");	gr->Puts(mglPoint(c,0.4),"Style 'II'",":rL");
	gr->Line(mglPoint(a,0.2),mglPoint(b,0.2),"kDD");	gr->Puts(mglPoint(c,0.2),"Style 'DD'",":rL");
	gr->Line(mglPoint(a,0),mglPoint(b,0),"kSS");		gr->Puts(mglPoint(c,0),"Style 'SS'",":rL");
	gr->Line(mglPoint(a,-0.2),mglPoint(b,-0.2),"kOO");	gr->Puts(mglPoint(c,-0.2),"Style 'OO'",":rL");
	gr->Line(mglPoint(a,-0.4),mglPoint(b,-0.4),"kTT");	gr->Puts(mglPoint(c,-0.4),"Style 'TT'",":rL");
	gr->Line(mglPoint(a,-0.6),mglPoint(b,-0.6),"k-__");	gr->Puts(mglPoint(c,-0.6),"Style '\\_\\_'",":rL");
	gr->Line(mglPoint(a,-0.8),mglPoint(b,-0.8),"k-VA");	gr->Puts(mglPoint(c,-0.8),"Style 'VA'",":rL");
	gr->Line(mglPoint(a,-1),mglPoint(b,-1),"k-AV");		gr->Puts(mglPoint(c,-1),"Style 'AV'",":rL");

	gr->SubPlot(3,2,2);
	//#LENUQ
	gr->FaceZ(mglPoint(-1,	-1), 0.4, 0.3, "L#");	gr->Puts(mglPoint(-0.8,-0.9), "L", "w:C", -1.4);
	gr->FaceZ(mglPoint(-0.6,-1), 0.4, 0.3, "E#");	gr->Puts(mglPoint(-0.4,-0.9), "E", "w:C", -1.4);
	gr->FaceZ(mglPoint(-0.2,-1), 0.4, 0.3, "N#");	gr->Puts(mglPoint(0,  -0.9), "N", "w:C", -1.4);
	gr->FaceZ(mglPoint(0.2,	-1), 0.4, 0.3, "U#");	gr->Puts(mglPoint(0.4,-0.9), "U", "w:C", -1.4);
	gr->FaceZ(mglPoint(0.6,	-1), 0.4, 0.3, "Q#");	gr->Puts(mglPoint(0.8,-0.9), "Q", "w:C", -1.4);
	//#lenuq
	gr->FaceZ(mglPoint(-1,	-0.7), 0.4, 0.3, "l#");	gr->Puts(mglPoint(-0.8,-0.6), "l", "k:C", -1.4);
	gr->FaceZ(mglPoint(-0.6,-0.7), 0.4, 0.3, "e#");	gr->Puts(mglPoint(-0.4,-0.6), "e", "k:C", -1.4);
	gr->FaceZ(mglPoint(-0.2,-0.7), 0.4, 0.3, "n#");	gr->Puts(mglPoint(0,  -0.6), "n", "k:C", -1.4);
	gr->FaceZ(mglPoint(0.2,	-0.7), 0.4, 0.3, "u#");	gr->Puts(mglPoint(0.4,-0.6), "u", "k:C", -1.4);
	gr->FaceZ(mglPoint(0.6,	-0.7), 0.4, 0.3, "q#");	gr->Puts(mglPoint(0.8,-0.6), "q", "k:C", -1.4);
	//#CMYkP
	gr->FaceZ(mglPoint(-1,	-0.4), 0.4, 0.3, "C#");	gr->Puts(mglPoint(-0.8,-0.3), "C", "w:C", -1.4);
	gr->FaceZ(mglPoint(-0.6,-0.4), 0.4, 0.3, "M#");	gr->Puts(mglPoint(-0.4,-0.3), "M", "w:C", -1.4);
	gr->FaceZ(mglPoint(-0.2,-0.4), 0.4, 0.3, "Y#");	gr->Puts(mglPoint(0,  -0.3), "Y", "w:C", -1.4);
	gr->FaceZ(mglPoint(0.2,	-0.4), 0.4, 0.3, "k#");	gr->Puts(mglPoint(0.4,-0.3), "k", "w:C", -1.4);
	gr->FaceZ(mglPoint(0.6,	-0.4), 0.4, 0.3, "P#");	gr->Puts(mglPoint(0.8,-0.3), "P", "w:C", -1.4);
	//#cmywp
	gr->FaceZ(mglPoint(-1,	-0.1), 0.4, 0.3, "c#");	gr->Puts(mglPoint(-0.8, 0), "c", "k:C", -1.4);
	gr->FaceZ(mglPoint(-0.6,-0.1), 0.4, 0.3, "m#");	gr->Puts(mglPoint(-0.4, 0), "m", "k:C", -1.4);
	gr->FaceZ(mglPoint(-0.2,-0.1), 0.4, 0.3, "y#");	gr->Puts(mglPoint(0,   0), "y", "k:C", -1.4);
	gr->FaceZ(mglPoint(0.2,	-0.1), 0.4, 0.3, "w#");	gr->Puts(mglPoint(0.4, 0), "w", "k:C", -1.4);
	gr->FaceZ(mglPoint(0.6,	-0.1), 0.4, 0.3, "p#");	gr->Puts(mglPoint(0.8, 0), "p", "k:C", -1.4);
	//#BGRHW
	gr->FaceZ(mglPoint(-1,	0.2), 0.4, 0.3, "B#");	gr->Puts(mglPoint(-0.8, 0.3), "B", "w:C", -1.4);
	gr->FaceZ(mglPoint(-0.6,0.2), 0.4, 0.3, "G#");	gr->Puts(mglPoint(-0.4, 0.3), "G", "w:C", -1.4);
	gr->FaceZ(mglPoint(-0.2,0.2), 0.4, 0.3, "R#");	gr->Puts(mglPoint(0,   0.3), "R", "w:C", -1.4);
	gr->FaceZ(mglPoint(0.2,	0.2), 0.4, 0.3, "H#");	gr->Puts(mglPoint(0.4, 0.3), "H", "w:C", -1.4);
	gr->FaceZ(mglPoint(0.6,	0.2), 0.4, 0.3, "W#");	gr->Puts(mglPoint(0.8, 0.3), "W", "w:C", -1.4);
	//#bgrhw
	gr->FaceZ(mglPoint(-1,	0.5), 0.4, 0.3, "b#");	gr->Puts(mglPoint(-0.8, 0.6), "b", "k:C", -1.4);
	gr->FaceZ(mglPoint(-0.6,0.5), 0.4, 0.3, "g#");	gr->Puts(mglPoint(-0.4, 0.6), "g", "k:C", -1.4);
	gr->FaceZ(mglPoint(-0.2,0.5), 0.4, 0.3, "r#");	gr->Puts(mglPoint(0,   0.6), "r", "k:C", -1.4);
	gr->FaceZ(mglPoint(0.2,	0.5), 0.4, 0.3, "h#");	gr->Puts(mglPoint(0.4, 0.6), "h", "k:C", -1.4);
	gr->FaceZ(mglPoint(0.6,	0.5), 0.4, 0.3, "w#");	gr->Puts(mglPoint(0.8, 0.6), "w", "k:C", -1.4);
	//#brighted
	gr->FaceZ(mglPoint(-1,	0.8), 0.4, 0.3, "{r1}#");	gr->Puts(mglPoint(-0.8, 0.9), "\\{r1\\}", "w:C", -1.4);
	gr->FaceZ(mglPoint(-0.6,0.8), 0.4, 0.3, "{r3}#");	gr->Puts(mglPoint(-0.4, 0.9), "\\{r3\\}", "w:C", -1.4);
	gr->FaceZ(mglPoint(-0.2,0.8), 0.4, 0.3, "{r5}#");	gr->Puts(mglPoint(0,   0.9), "\\{r5\\}", "k:C", -1.4);
	gr->FaceZ(mglPoint(0.2,	0.8), 0.4, 0.3, "{r7}#");	gr->Puts(mglPoint(0.4, 0.9), "\\{r7\\}", "k:C", -1.4);
	gr->FaceZ(mglPoint(0.6,	0.8), 0.4, 0.3, "{r9}#");	gr->Puts(mglPoint(0.8, 0.9), "\\{r9\\}", "k:C", -1.4);
	// HEX
	gr->FaceZ(mglPoint(-1, -1.3), 1, 0.3, "{xff9966}#");	gr->Puts(mglPoint(-0.5,-1.2), "\\{xff9966\\}", "k:C", -1.4);
	gr->FaceZ(mglPoint(0,  -1.3), 1, 0.3, "{x83CAFF}#");	gr->Puts(mglPoint( 0.5,-1.2), "\\{x83CAFF\\}", "k:C", -1.4);

	gr->SubPlot(3,2,3);
	char stl[3]="r1", txt[4]="'1'";
	for(int i=0;i<10;i++)
	{
		txt[1]=stl[1]='0'+i;
		gr->Line(mglPoint(-1,0.2*i-1),mglPoint(1,0.2*i-1),stl);
		gr->Puts(mglPoint(1.05,0.2*i-1),txt,":L");
	}

	gr->SubPlot(3,2,4);	gr->Title("TriPlot sample");	gr->Rotate(50,60);
	double t[] = {0,1,2, 0,1,3, 0,2,3, 1,2,3};
	double xt[] = {-1,1,0,0}, yt[] = {-1,-1,1,0}, zt[] = {-1,-1,-1,1};
	mglData tt(4,3,t), uu(4,xt), vv(4,yt), ww(4,zt);
	gr->TriPlot(tt,uu,vv,ww,"b");
	gr->TriPlot(tt,uu,vv,ww,"k#");

	gr->SubPlot(3,2,5);
	mglData r(4);	r.Fill(1,4);
	gr->SetRanges(1,4,1,4);	gr->Axis();
	gr->Mark(r,r,"s");
	gr->Plot(r,"b");

	gr->WriteJPEG("fexport.jpg");
//	gr->WritePNG("fexport.png");
	gr->WriteBMP("fexport.bmp");
	gr->WriteTGA("fexport.tga");
	gr->WriteEPS("fexport.eps");
	gr->WriteSVG("fexport.svg");
	gr->WriteGIF("fexport.gif");

	gr->WriteXYZ("fexport.xyz");
	gr->WriteSTL("fexport.stl");
	gr->WriteOFF("fexport.off");
	gr->WriteTEX("fexport.tex");
	gr->WriteOBJ("fexport.obj");
	gr->WritePRC("fexport.prc");
	gr->WriteJSON("fexport.json");

	gr->ExportMGLD("fexport.mgld");
	gr->Clf();
	gr->ImportMGLD("fexport.mgld");
}
//-----------------------------------------------------------------------------
int main(int argc,char **argv)
{
// const char *f = strrchr(argv[0],'/');
// std::string p(argv[0],f-argv[0]);
// printf("getcwd = '%s', argv = '%s', path = '%s', inst = '%s'\n", getcwd(NULL,0), argv[0], p.c_str(), MGL_INSTALL_DIR);
// fflush(stdout);

	mgl_textdomain(argv?argv[0]:NULL,"");
	mgl_suppress_warn(true);
	const char *suf = "";
	char name[256]="", *tmp;
	int ch;
	time_t st,en;	time(&st);
	mglGraph *gr = NULL;
	mglSample *s=samp;
	while(( ch = getopt_long_only(argc, argv, "", longopts, NULL)) != -1)
		switch(ch)
		{
			case 0:		break;
			case 'w':	width =atoi(optarg);	break;
			case 'h':	height=atoi(optarg);	break;
			case 'q':	quality =atoi(optarg);	break;
			case 'k':	strncpy(name, optarg,256);
						tmp=strchr(name,'.');	if(tmp)	*tmp=0;
						tmp=strchr(name,'-');	if(tmp)	*tmp=0;
						break;
			case 't':	mgl_set_num_thr(atoi(optarg));	break;
			case 'l':
				while(s->name[0])	{	printf("%s ",s->name);	s++;	}
				printf("\n");	return 0;
			case '?':
			default:	usage();	return 0;
		}

	if(dotest==1)	printf("Global (before):%s\n",mglGlobalMess.c_str());
	gr = new mglGraph;
	if(	type==11|| type==12|| type==5 || type==9)	width=height;
	switch(big)
	{
	case 1:	gr->SetSize(1920,1440);	suf = "-lg";	break;
	case 2:	gr->SetSize(640,480);	break;
	case 3:	gr->SetSize(192,144);	suf = "-sm";	break;
	default:	gr->SetSize(width,height);
	}
	gr->SetQuality(quality);

	if(dotest==1)
	{
		mgl_set_test_mode(true);	test(gr);
		time(&en);	printf("time is %g sec\n",difftime(en,st));
#if MGL_HAVE_PNG
		gr->WritePNG("test.png","",false);
#else
		gr->WriteBMP("test.bmp");
#endif
		gr->WriteSVG("test.svg");
		gr->WriteEPS("test.eps");
		printf("Messages:%s\n",gr->Message());
		printf("Global:%s\n",mglGlobalMess.c_str());
		delete gr;	return 0;
	}
	else if(dotest==2)	// NOTE mgl_gen_fnt[###][6] have to be updated if new glyphs will be added to built-in font
	{	mgl_create_cpp_font(gr->Self(), L"!-~,¡-ÿ,̀-̏,Α-ω,ϑ,ϕ,ϖ,ϰ,ϱ,ϵ,А-я,ℏ,ℑ,ℓ,ℜ,←-↙,∀-∯,≠-≯,⟂");
		delete gr;	return 0;	}
	else if(dotest==3)
	{
		int qual[7]={0,1,2,4,5,6,8};
		size_t ll=strlen(mmgl_dat_prepare)+1;
		mglParse par;
		par.AllowSetSize(true);
		FILE *fp = fopen(big?"time_big.texi":"time.texi","w");
		fprintf(fp,"@multitable @columnfractions .16 .12 .12 .12 .12 .12 .12 .12\n");
		fprintf(fp,"@headitem Name");
		for(int i=0;i<7;i++)	fprintf(fp," @tab q=%d",qual[i]);
		clock_t beg,end;
		while(s->name[0])	// all samples
		{
			char *buf = new char[strlen(s->mgl)+ll];
			strcpy(buf,s->mgl);	strcat(buf,mmgl_dat_prepare);
			fprintf(fp,"\n@item %s",s->name);

			printf("%s",s->name);
			for(int i=0;i<7;i++)
			{
				gr->DefaultPlotParam();
				gr->SetQuality(qual[i]);	gr->Clf();
				beg = clock();
				if(!use_mgl)	s->func(gr);
				else 	par.Execute(gr,buf);
				gr->Finish();
				end = clock();
				fprintf(fp," @tab %.3g",double(end-beg)/CLOCKS_PER_SEC);
				printf("\t%d->%g",qual[i],double(end-beg)/CLOCKS_PER_SEC);
				fflush(fp);	fflush(stdout);
			}
			printf("\n");	delete []buf;	s++;
		}
		fprintf(fp,"\n@end multitable\n");	fclose(fp);
	}
	else if(dotest==4)
	{	smgl_fexport(gr);	delete gr;	return 0;	}
	else if(dotest==5)
	{
		mgl_check_tex_table();
		delete gr;	return 0;
	}

	if(type==15 || type==16)	big=3;	// save mini version for json

	if(srnd)	mgl_srnd(1);
	gr->VertexColor(false);	gr->Compression(false);
	if(name[0]==0)
	{
		while(s->name[0])	// all samples
		{
			gr->DefaultPlotParam();	gr->Clf();
			if(use_mgl)
			{
				mglParse par;
				par.AllowSetSize(true);
				char *buf = new char[strlen(s->mgl)+strlen(mmgl_dat_prepare)+1];
				strcpy(buf,s->mgl);		strcat(buf,mmgl_dat_prepare);
				if(type!=7)	printf("\n-------\n%s\n-------\n",verbose?buf:s->mgl);
				par.Execute(gr,buf);	delete []buf;
				const char *mess = gr->Message();
				if(*mess)	printf("Warnings: %s\n-------\n",mess);
			}
			else	s->func(gr);
			save(gr, s->name, suf);
			printf("%s ",s->name);	fflush(stdout);	s++;
			gr->SetQuality(quality);
		}
		printf("\n");
	}
	else	// manual sample
	{
		mglSample tst;	tst.name=name;
		int i=0;
		for(i=0;samp[i].name[0];i++);	// determine the number of samples
		s = (mglSample *) bsearch(&tst, samp, i, sizeof(mglSample), mgl_cmd_smp);
		if(s)
		{
			gr->DefaultPlotParam();	gr->Clf();
			if(use_mgl)
			{
				mglParse par;
				par.AllowSetSize(true);
				char *buf = new char[strlen(s->mgl)+strlen(mmgl_dat_prepare)+1];
				strcpy(buf,s->mgl);		strcat(buf,mmgl_dat_prepare);
				if(type!=7)	printf("\n-------\n%s\n-------\n",verbose?buf:s->mgl);
				par.Execute(gr,buf);	delete []buf;
				const char *mess = gr->Message();
				if(*mess)	printf("Warnings: %s\n-------\n",mess);
			}
			else	s->func(gr);
			save(gr, s->name, suf);
		}
		else	printf("no sample %s\n",name);
	}
	delete gr;	return 0;
}
//-----------------------------------------------------------------------------

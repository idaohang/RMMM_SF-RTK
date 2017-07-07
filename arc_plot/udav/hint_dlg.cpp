/***************************************************************************
 *   Copyright (C) 2008 by Alexey Balakin                                  *
 *   mathgl.abalakin@gmail.com                                             *
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
#include <QSettings>
#include <QCheckBox>
#include <QLayout>
#include <QPushButton>
#include "hint_dlg.h"
#include "mgl2/data_cf.h"
//-----------------------------------------------------------------------------
//
//	Hint dialog
//
//-----------------------------------------------------------------------------
HintDialog::HintDialog(QWidget *parent) : QDialog(parent)
{
	hints.append(_("You can shift axis range by pressing middle button and moving mouse. Also, you can zoom in/out axis range by using mouse wheel."));
	hints.append(_("You can rotate/shift/zoom whole plot by mouse. Just press 'Rotate' toolbutton, click image and hold a mouse button: left button for rotation, right button for zoom/perspective, middle button for shift."));
	hints.append(_("You may quickly draw the data from file. Just use: udav 'filename.dat' in command line."));
	hints.append(_("You can copy the current image to clipboard by pressing Ctrl-Shift-C. Later you can paste it directly into yours document or presentation."));
	hints.append(_("You can export image into a set of format (EPS, SVG, PNG, JPEG) by pressing right mouse button inside image and selecting 'Export as ...'."));
	hints.append(_("You can setup colors for script highlighting in Property dialog. Just select menu item 'Settings/Properties'."));
	hints.append(_("You can save the parameter of animation inside MGL script by using comment started from '##a ' or '##c ' for loops."));
	hints.append(_("New drawing never clears things drawn already. For example, you can make a surface with contour lines by calling commands 'surf' and 'cont' one after another (in any order). "));
	hints.append(_("You can put several plots in the same image by help of commands 'subplot' or 'inplot'."));
	hints.append(_("All indexes (of data arrays, subplots and so on) are always start from 0."));
	hints.append(_("You can edit MGL file in any text editor. Also you can run it in console by help of commands: mglconv, mglview."));
	hints.append(_("You can use command 'once on|off' for marking the block which should be executed only once. For example, this can be the block of large data reading/creating/handling. Press F9 (or menu item 'Graphics/Reload') to re-execute this block."));
	hints.append(_("You can use command 'stop' for terminating script parsing. It is useful if you don't want to execute a part of script."));
	hints.append(_("You can type arbitrary expression as input argument for data or number. In last case (for numbers), the first value of data array is used."));
	hints.append(_("There is powerful calculator with a lot of special functions. You can use buttons or keyboard to type the expression. Also you can use existed variables in the expression."));
	hints.append(_("The calculator can help you to put complex expression in the script. Just type the expression (which may depend on coordinates x,y,z and so on) and put it into the script."));
	hints.append(_("You can easily insert file or folder names, last fitted formula or numerical value of selection by using menu Edit|Insert."));
	hints.append(_("The special dialog (Edit|Insert|New Command) help you select the command, fill its arguments and put it into the script."));
	hints.append(_("You can put several plotting commands in the same line or in separate function, for highlighting all of them simultaneously."));

	numHints=hints.size();
	cur = int(mgl_rnd()*numHints);
	setWindowTitle(_("UDAV - Hint"));
	QHBoxLayout *a;
	QPushButton *b;
	QVBoxLayout *o = new QVBoxLayout(this);
	text = new QTextEdit(this);	o->addWidget(text);
	text->setReadOnly(true);	text->setText(hints[cur]);

	start = new QCheckBox(_("Show at startup"), this);	o->addWidget(start);
	start->setChecked(true);

	a = new QHBoxLayout;	o->addLayout(a);
	b = new QPushButton(_("Prev"), this);		a->addWidget(b);
	connect(b, SIGNAL(clicked()), this, SLOT(prevClicked()));
	b = new QPushButton(_("Next"), this);		a->addWidget(b);
	connect(b, SIGNAL(clicked()), this, SLOT(nextClicked()));
	b = new QPushButton(_("Close"), this);	a->addWidget(b);
	connect(b, SIGNAL(clicked()),this, SLOT(close()));
}
//-----------------------------------------------------------------------------
void HintDialog::closeEvent(QCloseEvent *)
{
	QSettings settings("udav","UDAV");
	settings.setPath(QSettings::IniFormat, QSettings::UserScope, "UDAV");
	settings.beginGroup("/UDAV");
	settings.setValue("/showHint", start->isChecked());
	settings.endGroup();
}
//-----------------------------------------------------------------------------
void udavShowHint(QWidget *p)
{	HintDialog *hd = new HintDialog(p);	hd->exec();	}
//-----------------------------------------------------------------------------

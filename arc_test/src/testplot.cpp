
#include "arc.h"

using namespace std;

extern int arc_plot(int argc,char *argv[])
{
#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
    QApplication::setGraphicsSystem("raster");
#endif
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
int ExecCmd(const QString &cmd, int show)
{
    Q_UNUSED(show);
    return QProcess::startDetached(cmd);  /* FIXME: show option not yet supported */
}
#define NUMINFILE 10

int main(int argc, char *argv[])
{
    const char *roverobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hklt/hklt1910.17o";
    const char *baseobsf = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hkkt/hkkt1910.17o";
    const char *navf  = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hklt/hklt1910.17n";
    const char *navf1 = "/home/sujinglan/arc_rtk/arc_test/data/gps_bds/hkcors/hkkt/hkkt1910.17f";

    QString QRoverObs=QString(QLatin1String(roverobsf));
    QString QBaseObs=QString(QLatin1String(baseobsf));
    QString QNav1=QString(QLatin1String(navf));
    QString QNav2=QString(QLatin1String(navf1));

    QString cmd="/home/sujinglan/arc_rtk/arc_plot/bin/rtkplot_qt",opts="";

    opts=" -r \""+QRoverObs+"\" \""+QBaseObs+"\" \""+QNav1+"\" \""+ QNav2;

    if (!ExecCmd(cmd+opts,1)) {
        QMessageBox *msgBox;
        msgBox = new QMessageBox("ARC-SRTK",
                                 "Plot Error",
                                 QMessageBox::Critical,
                                 QMessageBox::Ok|QMessageBox::Default,
                                 QMessageBox::Cancel|QMessageBox::Escape,
                                 0);
        msgBox->show();
    }
    return 0;
}


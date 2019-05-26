#include "dep/qcustomplot.h"
#include "view/widget.h"
#include <cmath>

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Widget w;
    w.show();

    return a.exec();
}

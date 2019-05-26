#include "widget.h"
#include "ui_widget.h"
#include <iostream>
#include <cmath>

// --------------------------
// Set stop time here
// --------------------------
#define ENDOFTIME 50
// --------------------------
// Set stop time here
// --------------------------

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    ui->setupUi(this);

    // Set window size
    this->setFixedSize(1400,700);

    // Add main layout with two plots
    mainlayout = new QGridLayout(this);
    inputPlot = new QCustomPlot(this);
    outputPlot = new QCustomPlot(this);
    mainlayout->addWidget(inputPlot,0,0);
    mainlayout->addWidget(outputPlot,0,1);
    inputPlot->setFixedSize(this->width()/3,this->height());
    outputPlot->setFixedSize(this->width()/3,this->height());

    // Give the axes some labels:
    inputPlot->xAxis->setLabel("t");
    inputPlot->yAxis->setLabel("input");
    outputPlot->xAxis->setLabel("t");
    outputPlot->yAxis->setLabel("output");

    // --------------------------
    // Change ranges if you need
    // --------------------------
    // Set axes ranges so see all data:
    inputPlot->xAxis->setRange(0, ENDOFTIME);
    inputPlot->yAxis->setRange(-3, 3);
    outputPlot->xAxis->setRange(0, ENDOFTIME);
    outputPlot->yAxis->setRange(-3, 3);

    // --------------------------
    // Create the object here
    // --------------------------
    object = new Model();   // <=
    // --------------------------
    // Create the object here
    // --------------------------

    // Get time in msec
    // --------------------------
    // Google for MacOS timings
    // --------------------------
#ifdef __linux__
    struct timeval tmpStruct;
    gettimeofday(&tmpStruct, NULL);
    startTime = tmpStruct.tv_sec * 1000 + tmpStruct.tv_usec / 1000 + 0.5;
#endif
#ifdef _WIN32
    SYSTEMTIME tmpStruct;
    GetSystemTime(&tmpStruct);
    startTime = tmpStruct.wSecond * 1000 + tmpStruct.wMilliseconds + 0.5;
#endif

    makePlot();
    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(makePlot()));

    // --------------------------
    // Set sampling time here
    // --------------------------
    timer->start(100);
    // --------------------------
    // Set sampling time here
    // --------------------------
}

Widget::~Widget()
{
    delete ui;
    delete inputPlot;
    delete outputPlot;
    delete timer;
    delete mainlayout;

    // --------------------------
    // Delete the object here
    // --------------------------
    delete object;
    // --------------------------
    // Delete the object here
    // --------------------------
}

void Widget::makePlot() {

    // generate some data:


#ifdef __linux__
    struct timeval tmpTime;
    gettimeofday(&tmpTime, NULL);
    double tmp = (tmpTime.tv_sec * 1000 + tmpTime.tv_usec / 1000 + 0.5)-startTime;
#endif
#ifdef _WIN32
    SYSTEMTIME tmpTime;
    GetSystemTime(&tmpTime);
    double tmp = tmpTime.wSecond * 1000 + tmpTime.wMilliseconds + 0.5 - startTime;
#endif

    // --------------------------
    // Replace input signal with ours
    // --------------------------
    // double signal = std::sin(tmp/1000);
    double signal = 1;
    // --------------------------
    // Replace input signal with ours
    // --------------------------

    // Update input array to plot
    input.append(signal);

    // Get elapsed time
    if (time.empty()) {
        dt = 0;
    } else {
        dt = tmp / 1000.0 - time.last();
    }
    // Update time array to plot
    time.append(tmp/1000);

    // --------------------------
    // Update the object here
    // --------------------------
    output.append(object->update(signal));
    // --------------------------
    // Update the object here
    // --------------------------

    inputPlot->addGraph();
    inputPlot->graph(0)->setData(time, input);

    outputPlot->addGraph();
    outputPlot->graph(0)->setData(time, output);

    inputPlot->replot();
    outputPlot->replot();

    if (tmp/1000 > ENDOFTIME) {timer->stop();}
}

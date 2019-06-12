#include "widget.h"
#include "ui_widget.h"
#include <iostream>
#include <cmath>

Widget::Widget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Widget)
{
    /// Create the object
    ///////////////////////////////////////////////////////////
    object = new Model();
    // object->send(1);
    ///////////////////////////////////////////////////////////

    ui->setupUi(this);

    /// Set window size
    this->setFixedSize(1400,700);

    /// Add main layout with two plots
    mainlayout = new QGridLayout(this);
    inputPlot = new QCustomPlot(this);
    outputPlot = new QCustomPlot(this);
    mainlayout->addWidget(inputPlot,0,0);
    mainlayout->addWidget(outputPlot,0,1);
    inputPlot->setFixedSize(this->width()/2,this->height());
    outputPlot->setFixedSize(this->width()/2,this->height());

    /// Give the axes some labels:
    inputPlot->xAxis->setLabel("t");
    inputPlot->yAxis->setLabel("input");
    outputPlot->xAxis->setLabel("t");
    outputPlot->yAxis->setLabel("output");

    /// Change ranges if you need
    ///////////////////////////////////////////////////////////
    /// Set axes ranges so see all data:
    inputPlot->xAxis->setRange(0, object->SIMULATION_TIME);
    inputPlot->yAxis->setRange(-3, 3);
    outputPlot->xAxis->setRange(0, object->SIMULATION_TIME);
    outputPlot->yAxis->setRange(-1, 1);

    /// Get time in msec
    ///////////////////////////////////////////////////////////
    #ifdef __linux__
        struct timeval tmpStruct;
        gettimeofday(&tmpStruct, nullptr);
        startTime = tmpStruct.tv_sec * 1000 + tmpStruct.tv_usec / 1000 + 0.5;
    #endif
    #ifdef _WIN32
        SYSTEMTIME tmpStruct;
        GetSystemTime(&tmpStruct);
        startTime = tmpStruct.wSecond * 1000 + tmpStruct.wMilliseconds + 0.5;
    #endif
    ///////////////////////////////////////////////////////////

    makePlot();
    timer = new QTimer(this);
    connect(timer, SIGNAL(timeout()), this, SLOT(makePlot()));

    /// Set sampling time
    ///////////////////////////////////////////////////////////
    timer->start(object->TIME_STEP);
    ///////////////////////////////////////////////////////////
}

Widget::~Widget()
{
    delete ui;
    delete inputPlot;
    delete outputPlot;
    delete timer;
    delete mainlayout;

    /// Delete the object
    ///////////////////////////////////////////////////////////
    delete object;
    ///////////////////////////////////////////////////////////
}

void Widget::makePlot() {
// generate some data:
#ifdef __linux__
    struct timeval tmpTime;
    gettimeofday(&tmpTime, nullptr);
    double tmp = (tmpTime.tv_sec * 1000 + tmpTime.tv_usec / 1000 + 0.5)-startTime;
#endif
#ifdef _WIN32
    SYSTEMTIME tmpTime;
    GetSystemTime(&tmpTime);
    double tmp = tmpTime.wSecond * 1000 + tmpTime.wMilliseconds + 0.5 - startTime;
#endif

    /// Replace input signal with ours
    ///////////////////////////////////////////////////////////
    // double signal = std::sin(tmp/1000);
    double signal = object->control();
    // double signal = 1;
    // object->send(signal);
    ///////////////////////////////////////////////////////////

    /// Update input array to plot
    input.append(signal);

    /// Get elapsed time
    if (time.empty()) {
        dt = 0;
    } else {
        dt = tmp / 1000.0 - time.last();
    }

    object->time_now += object->TIME_STEP;

    /// Update time array to plot
    time.append(object->time_now);

    /// Update the object here
    ///////////////////////////////////////////////////////////
//    output.append(object->update_discrete(signal));
    output.append(object->update_continuous(signal, object->TIME_STEP));
    ///////////////////////////////////////////////////////////

    inputPlot->addGraph();
    inputPlot->graph(0)->setData(time, input);

    outputPlot->addGraph();
    outputPlot->graph(0)->setData(time, output);

    inputPlot->replot();
    outputPlot->replot();

    if (object->time_now > object->SIMULATION_TIME) {timer->stop();}
}

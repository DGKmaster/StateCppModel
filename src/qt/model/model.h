#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <inttypes.h>
#include <QSerialPort>
#include <sstream>

#include "matrix.h"
#include "integrator.h"

class Model
{
/// Attributes
///////////////////////////////////////////////////////////
private:
    const double OMEGA = 0.1;

    const QString SERIAL_PORT_NAME = "/dev/pts/11";
    const int SERIAL_PORT_BAUD_RATE = QSerialPort::Baud115200;
    const QSerialPort::DataBits DATA_BITS = QSerialPort::Data8;
    const QSerialPort::Parity PARITY = QSerialPort::NoParity;
    const QSerialPort::StopBits STOP_BITS = QSerialPort::OneStop;

    Matrix Ad = Matrix(3, 3);
    Matrix Bd = Matrix(3, 1);
    Matrix Cd = Matrix(1, 3);
    Matrix x = Matrix(3, 1);
    Matrix y = Matrix(1, 1);
	Matrix x_prev = Matrix(3, 1);
    Matrix u_state = Matrix(2, 1);

    Integrator _du_state = Integrator(-3 * OMEGA * 0.84147, 0);
    Integrator _u_state = Integrator(3 * 0.54030, 0);
	double feedback;

    QSerialPort serial_port;

public:
    /// 0.2 -> 5 Hz
    /// 0.02 -> 50 Hz
    /// 0.01 -> 100 Hz
    const double TIME_STEP = 0.02;
    const double SIMULATION_TIME = 100;
    double time_now = 0;

///////////////////////////////////////////////////////////

/// Methods
///////////////////////////////////////////////////////////
public:
    /// Constructors and destructors
    ///////////////////////////////////////////////////////////
    Model();

    ~Model();
    ///////////////////////////////////////////////////////////

    /// 
    ///////////////////////////////////////////////////////////
    void send(const double& value);

    double update(const double& input);

    double control();
    ///////////////////////////////////////////////////////////
};
#endif // MODEL_H

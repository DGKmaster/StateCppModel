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

    Matrix Ad = Matrix(3, 3);
    Matrix Bd = Matrix(3, 1);
    Matrix Cd = Matrix(1, 3);
    Matrix x = Matrix(3, 1);
    Matrix y = Matrix(1, 1);
	Matrix x_prev = Matrix(3, 1);
    Matrix u_state = Matrix(2, 1);

    Integrator _du_state = Integrator(-3*OMEGA*std::sin(1), 0);
    Integrator _u_state = Integrator(3 * std::cos(1), 0);
	double feedback;

    QSerialPort serial_port;

public:
    const double TIME_STEP = 0.1;
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

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
    Model() {
        ///
        ///////////////////////////////////////////////////////////
        double _In[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        Matrix In(3, 3, _In);

        double _A[] = {0, 1, 0, 0, 0, 1, -1.5, -5, -2};
        Matrix A = Matrix(3, 3, _A);

		this->Ad = A.matExp(TIME_STEP);

        double _B[] = {0, 0, 1};
        Matrix B = Matrix(3, 1, _B);
        this->Bd = A.inverse()*(Ad - In)*B;

        double _C[] = {0.5, 0, 0};
        Matrix C = Matrix(1, 3, _C);
        this->Cd = C;

        double _x[] = {0, 0, 0};
        this->x = Matrix(3, 1, _x);

        this->y = Matrix(1, 1);

        /// cos(1)*amplitude
        /// -sin(1rad)*omega*amplitude
        double _u_state[] = {0.54*3, 0.84147*-0.3};
        this->u_state = Matrix(2, 1, _u_state);
        ///////////////////////////////////////////////////////////
        
        ///
        ///////////////////////////////////////////////////////////
        this->serial_port.setPortName(this->SERIAL_PORT_NAME);
        this->serial_port.setBaudRate(this->SERIAL_PORT_BAUD_RATE);

        // if (!serial_port.open(QIODevice::WriteOnly)) {
        //     std::cout << "Can not open port" << std::endl;
        // }
        ///////////////////////////////////////////////////////////
    }

    ~Model() {
        serial_port.close();
    }
    ///////////////////////////////////////////////////////////

    /// 
    ///////////////////////////////////////////////////////////
    void send(const double& value) {
        std::stringstream os;
        os << value;
        QByteArray write_data(os.str().c_str());
        qint64 bytes_written = serial_port.write(write_data);

        if (bytes_written == -1) {
            std::cout << "Failed to write the data" << std::endl;
        } else if (bytes_written != write_data.size()) {
            std::cout << "Failed to write all the data" << std::endl;
        } else if (!serial_port.waitForBytesWritten(5000)) {
            std::cout << "Operation timed out or an error" << std::endl;  
        }
    }

    double update(const double& input) {
		x_prev = x;
        x = Ad*x + Bd*input;
        y = Cd*x_prev;

        double output = this->y.getM();

        return output;
    }

    double control() {
        /// Constants
        ///////////////////////////////////////////////////////////
        const double T = this->TIME_STEP;

        const double _In[] = {1, 0, 0, 1};
        const Matrix In(2, 2, _In);

        const double _u_A[] = {0, -1, OMEGA*OMEGA, 0};
        const Matrix u_A = Matrix(2, 2, _u_A);
        ///////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////
        Matrix mat_exp = Matrix(2, 1);
        Matrix u = Matrix(2, 1);
        
		double signal = _u_state.getState();
	
        this->_du_state.update(feedback, T);
		this->_u_state.update(_du_state.getState(), T);
		
		feedback = -OMEGA*OMEGA*_u_state.getState();
        ///////////////////////////////////////////////////////////

        return signal;
    }
    ///////////////////////////////////////////////////////////
};
#endif // MODEL_H

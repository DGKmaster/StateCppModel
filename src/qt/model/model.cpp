#include "model.h"

Model::Model() {
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
    this->serial_port.setDataBits(this->DATA_BITS);
    this->serial_port.setParity(this->PARITY);
    this->serial_port.setStopBits(this->STOP_BITS);
    // if (!serial_port.open(QIODevice::WriteOnly)) {
    //     std::cout << "Can not open port" << std::endl;
    // }
    ///////////////////////////////////////////////////////////
}

Model::~Model() {
    serial_port.close();
}

void Model::send(const double& value) {
    const uint LENGTH = 20;
    char buffer_char[LENGTH];
    uint8_t buffer_uint[LENGTH];

    snprintf(buffer_char, LENGTH, "%A", value);

    sscanf(buffer_char, "%p", &buffer_uint);

	cobs::encode(buffer_uint, LENGTH, 1);

    std::stringstream os;
    os << buffer_uint;
    
    // std::cout << std::hex << buffer_uint << std::endl;
    
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

double Model::update(const double& input) {
    this->x_prev = this->x;
    this->x = this->Ad*x + this->Bd*input;
    this->y = this->Cd*this->x_prev;

    double output = this->y.getM();

    return output;
}

double Model::control() {
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

#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <inttypes.h>
#include <QSerialPort>
#include <sstream>

class Matrix {
    private:
        /// Number of rows
        const uint16_t num_rows;
        /// Number of columns
        const uint16_t num_columns;
        /// Pointer to matrix elements in 1D
        double* mx;

    public:
        /// Constructors and destructors
        ///////////////////////////////////////////////////////////
        /**
        * @brief Default constructor
        * @return Nothing
        */
        Matrix():
            num_rows(0),
            num_columns(0) {
            this->mx = nullptr;
        }

        /**
         * @brief Default destructor
         * @return Nothing
         */
        ~Matrix() {
            delete[] mx;
        }

        /**
         * @brief Initialise matrix RxC with specified array values
         * @param r Number of rows
         * @param c Number of columns
         * @param m Matrix elements
         */
        Matrix(const uint16_t& r, const uint16_t& c, const double m[]):
            num_rows(r),
            num_columns(c),
            mx(new double[num_rows * num_columns]) {
            for(uint16_t i = 0; i < num_rows; i++) {
                for(uint16_t j = 0; j < num_columns; j++) {
                    this->mx[i*num_columns + j] = m[i*num_columns + j];
                }
            }
        }

        /**
         * @brief Initialise matrix RxC with zero values
         * @param r Number of rows
         * @param c Number of columns
         */
        Matrix(const uint16_t& r, const uint16_t& c):
            num_rows(r),
            num_columns(c),
            mx(new double[num_rows * num_columns]) {
            for(uint16_t i = 0; i < num_rows; i++) {
                for(uint16_t j = 0; j < num_columns; j++) {
                    this->mx[i*num_columns + j] = 0;
                }
            }
        }
        ///////////////////////////////////////////////////////////
        
        /// Accessors
        ///////////////////////////////////////////////////////////
        double getM() {
            double output = this->mx[0];

            return output;
        }
        ///////////////////////////////////////////////////////////

        /// Operators overload
        ///////////////////////////////////////////////////////////
        /**
         * @brief Element-wise multiplication
         * @param matrix Right operand
         * @return Result
         */
        Matrix operator*(const Matrix &matrix) const {
            if(this->num_columns != matrix.num_rows) {
                std::cout << "Bad dimensions =" << std::endl;
            }

            Matrix matrix_out(this->num_rows, matrix.num_columns);

            /// Multiplying matrix a and b and storing in array
            for(uint16_t i = 0; i < this->num_rows; i++) {
                for(uint16_t j = 0; j < matrix.num_columns; j++) {
                    for(uint16_t k = 0; k < this->num_columns; k++) {
                        matrix_out.mx[i*matrix_out.num_columns + j] += this->mx[i*this->num_columns + k] * matrix.mx[k*matrix.num_columns + j];
                    }
                }
            }

            return matrix_out;
        }

        Matrix operator*(const double &value) const {
            const uint16_t N = this->num_rows;
            const uint16_t M = this->num_columns;

            Matrix matrix_out(N, M);

            for(uint16_t i = 0; i < N; i++) {
                for(uint16_t j = 0; j < M; j++) {
                    matrix_out.mx[i*M + j] = this->mx[i*M + j] * value;
                }
            }

            return matrix_out;
        }

        Matrix operator/(const double &value) const {
            const uint16_t N = this->num_rows;
            const uint16_t M = this->num_columns;

            Matrix matrix_out(N, M);

            for(uint16_t i = 0; i < N; i++) {
                for(uint16_t j = 0; j < M; j++) {
                    matrix_out.mx[i*M + j] = this->mx[i*M + j] / value;
                }
            }

            return matrix_out;
        }

        /**
         * @brief Element-wise sum
         * @param matrix Right operand
         * @return Result
         */
        Matrix operator+(const Matrix &matrix) const {
            const uint16_t N = this->num_rows;
            const uint16_t M = this->num_columns;

            if(N != matrix.num_rows && M != matrix.num_columns) {
                std::cout << "Bad dimensions +" << std::endl;
            }

            Matrix matrix_out(N, M);

            for(uint16_t i = 0; i < N; i++) {
                for(uint16_t j = 0; j < M; j++) {
                    matrix_out.mx[i*M + j] = this->mx[i*M + j] + matrix.mx[i*M + j];
                }
            }

            return matrix_out;
        }

        /**
         * @brief Element-wise subtraction
         * @param matrix Right operand
         * @return Result
         */
        Matrix operator-(const Matrix &matrix) const {
            const uint16_t N = this->num_rows;
            const uint16_t M = this->num_columns;

            if(N != matrix.num_rows && M != matrix.num_columns) {
                std::cout << "Bad dimensions -" << std::endl;
            }

            Matrix matrix_out(N, M);

            /// Multiplying matrix a and b and storing in array
            for(uint16_t i = 0; i < N; i++) {
                for(uint16_t j = 0; j < M; j++) {
                    matrix_out.mx[i*M + j] = this->mx[i*M + j] - matrix.mx[i*M + j];
                }
            }

            return matrix_out;
        }

        /**
         * @brief Element-wise assignment
         * @param matrix rvalue
         */
        void operator=(const Matrix &matrix) {
            const uint16_t N = this->num_rows;
            const uint16_t M = this->num_columns;

            for(uint16_t i = 0; i < N; i++) {
                for(uint16_t j = 0; j < M; j++) {
                    this->mx[i*M + j] = matrix.mx[i*M + j];
                }
            }
        }
        ///////////////////////////////////////////////////////////

        /// Other functions
        ///////////////////////////////////////////////////////////
        /**
         * @brief Display matrix using stdout
         */
        void show(const std::string& text_out = "Output matrix: ") const {
            const uint16_t N = this->num_rows;
            const uint16_t M = this->num_columns;

            std::cout << text_out << std::endl;

            for(uint16_t i = 0; i < N; ++i) {
                for(int16_t j = 0; j < M; ++j) {
                    std::cout << " " << this->mx[i*M + j];
                    if(j == M - 1) {
                        std::cout << std::endl;
                    }
                }
            }
        }

        /**
        * @brief Function to get cofactor of A[p][q] in temp[][]. n is current dimension of A[][]
        * @param A
        * @param temp
        * @param p
        * @param q
        * @param n
         */
        Matrix getCofactor(const uint16_t& p, const uint16_t& q, const int16_t& n) {
            const uint16_t N = this->num_rows;

            Matrix matrix_out(N, N);

            uint16_t i = 0;
            uint16_t j = 0;

            /// Looping for each element of the matrix
            for (uint16_t row = 0; row < n; row++) {
                for (uint16_t col = 0; col < n; col++) {
                    ///  Copying into temporary matrix only those element
                    ///  which are not in given row and column
                    if (row != p && col != q) {
                        matrix_out.mx[i*N + j] = this->mx[row*N + col];
                        j++;

                        /// Row is filled, so increase row index and
                        /// reset col index
                        if (j == n - 1) {
                            j = 0;
                            i++;
                        }
                    }
                }
            }

            return matrix_out;
        }

        /**
        * @brief Recursive function for finding determinant of matrix. n is current dimension of A[][].
        * @param A
        * @param n
        * @return
        */
        double determinant(const uint16_t& n) {
            const uint16_t N = this->num_rows;

            /// Initialize result
            double D = 0;

            /// Base case : if matrix contains single element
            if (n == 1) {
                return this->mx[0];
            }

            /// To store cofactors
            Matrix _temp(N, N);

            /// To store sign multiplier
            int16_t _sign = 1;

            /// Iterate for each element of first row
            for (uint16_t i = 0; i < n; i++) {
                /// Getting Cofactor
                _temp = this->getCofactor(0, i, n);

                D += _sign * this->mx[0*N + i] * _temp.determinant(n - 1);

                /// terms are to be added with alternate sign
                _sign = -_sign;
            }

            return D;
        }

        /**
         * @brief Function to get adjoint of A[N][N] in adj[N][N].
         * @param A
         * @param adj
         */
        Matrix adjoint() {
            const uint16_t N = this->num_rows;

            //if (N == 1) {
            //    double _out[] = {this->mx[0]};
            //    Matrix out(1, 1, _out);
            //    return out;
            //}

            Matrix matrix_out(N, N);

            /// temp is used to store cofactors of A[][]
            int16_t _sign = 1;
            Matrix _temp(N, N);

            for (uint16_t i = 0; i < N; i++) {
                for (uint16_t j = 0; j < N; j++) {
                    /// Get cofactor of A[i][j]
                    _temp = this->getCofactor(i, j, N);

                    /// sign of adj[j][i] positive if sum of row
                    /// and column indexes is even.
                    _sign = ((i + j) % 2 == 0)? 1: -1;

                    /// Interchanging rows and columns to get the
                    /// transpose of the cofactor matrix
                    matrix_out.mx[j*N + i] = _sign*_temp.determinant(N - 1);
                }
            }

            return matrix_out;
        }

        /**
        * @brief Function to calculate and store inverse, returns false if matrix is singular
        */
        Matrix inverse() {
            const uint16_t N = this->num_rows;

            /// Find determinant of A[][]
            double det = this->determinant(N);

            if (det == 0.0) {
                std::cout << "Singular matrix, can't find its inverse" << std::endl;
            }
            Matrix matrix_out(N, N);

            /// Find adjoint
            Matrix adj(N, N);
            adj = this->adjoint();

            /// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
            for (uint16_t i = 0; i < N; i++) {
                for (uint16_t j = 0; j < N; j++) {
                    matrix_out.mx[i*N + j] = adj.mx[i*N + j] / det;
                }
            }

            return matrix_out;
        }

		Matrix MatExp(double dt) {

			const uint16_t N = this->num_rows;
			const uint16_t M = this->num_columns;
			
			double random_init[] = {1, 1, 1};
			Matrix x_internal_init(N, 1, random_init);
			Matrix x_internal_init_transpose(1, N, random_init);
			Matrix x_internal(N, 1, random_init);

			Matrix x_der1_internal_init(N, 1);
			x_der1_internal_init = *this * x_internal;

			Matrix x_der1_internal(N, 1);
			x_der1_internal = x_der1_internal_init;

			Matrix matrix_out(N, N);

			for (uint16_t i = 0; i < N; i++) {
				x_internal.mx[i] = x_internal.mx[i] + x_der1_internal.mx[i] * dt;
			}

			Matrix inverse_matrix_mult = x_internal_init * x_internal_init_transpose;
			matrix_out = x_internal * x_internal_init_transpose * inverse_matrix_mult.inverse();
			
			return matrix_out;
		}


		Matrix MatExp_v2(double dt) {
			const uint16_t N = this->num_rows;

			double _matrix_out[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
			Matrix matrix_out(N, N, _matrix_out);

			Matrix A_power(N, N);
			A_power = *this;
			double dt_power = dt;
			double factorial = 1;
            
            /// First step - just matrix_out
            /// Second step - matrix_out + ...
			matrix_out = matrix_out + (A_power * dt_power) / factorial;

            /// Third step and following
			for (uint16_t k = 2; k < 4; k++) {
				A_power = A_power * *this;
				dt_power = dt_power * dt;
                factorial = factorial * k;

				matrix_out = matrix_out + (A_power * dt_power) / factorial;
			}

			return matrix_out;
		}
        ///////////////////////////////////////////////////////////
};


class Integrator {
public:
	Integrator(double init_state, double init_in) {state = init_state;  prev_in = init_in; }

	void update(double input, double dt) {
        this->state_prev = this->state;
        this->state = state_prev + input*dt;
        this->prev_in = input;
	}

    double getState_prev() {
		return this->state_prev;
	}

	double getState() {
		return this->state;
	}

private:
	double state;
    double state_prev;
	double prev_in = 0;
};

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
        ///////////////////////////////////////////////////////////
        double _In[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        Matrix In(3, 3, _In);

        // double _A[] = {0, 1, 0, 0, 0, 1, -0.2, -1, -1};
        double _A[] = {0, 1, 0, 0, 0, 1, -1.5, -5, -2};
        Matrix A = Matrix(3, 3, _A);
		// matexp for TIME_STEP = 0.1
		double _Ad_test[] = { 0.999762603520798,0.099202680178427,0.004663356472870, -0.006995034709305,0.976445821156447,0.089875967232687, -0.134813950849031,-0.456374870872741,0.796693886691073 };
		Matrix Ad_test = Matrix(3, 3, _Ad_test);
        
        // this->Ad = Ad_test;
        // this->Ad = A.MatExp(TIME_STEP);
		this->Ad = A.MatExp_v2(TIME_STEP);
        (this->Ad - Ad_test).show();

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

        Matrix mat_exp = Matrix(2, 1);
        Matrix u = Matrix(2, 1);
        
        /// New method
        ///////////////////////////////////////////////////////////
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

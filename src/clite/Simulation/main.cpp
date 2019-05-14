#include <iostream>

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
            for(uint8_t i = 0; i < num_rows; i++) {
                for(uint8_t j = 0; j < num_columns; j++) {
                    this->mx[i*num_columns + j] =  m[i*num_columns + j];
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
            for(uint8_t i = 0; i < num_rows; i++) {
                for(uint8_t j = 0; j < num_columns; j++) {
                    this->mx[i*num_columns + j] = 0;
                }
            }
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
                std::cout << "Bad dimensions" << std::endl;
            }

            Matrix matrix_out(this->num_rows, matrix.num_columns);

            /// Multiplying matrix a and b and storing in array
            for(uint8_t i = 0; i < this->num_rows; ++i) {
                for(uint8_t j = 0; j < matrix.num_columns; ++j) {
                    for(uint8_t k = 0; k < this->num_columns; ++k) {
                        matrix_out.mx[i*matrix_out.num_columns + j] += matrix.mx[i*matrix.num_columns + k] * this->mx[k*this->num_columns + j];
                    }
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
            if(this->num_columns != matrix.num_columns && this->num_rows != matrix.num_rows) {
                std::cout << "Bad dimensions" << std::endl;
            }

            Matrix matrix_out(this->num_rows, this->num_columns);

            /// Multiplying matrix a and b and storing in array
            for(uint8_t i = 0; i < matrix_out.num_rows; ++i) {
                for(uint8_t j = 0; j < matrix_out.num_columns; ++j) {
                    matrix_out.mx[i*matrix_out.num_columns + j] = this->mx[i*matrix.num_columns + j] + matrix.mx[i*matrix.num_columns + j];
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
            if(this->num_columns != matrix.num_columns && this->num_rows != matrix.num_rows) {
                std::cout << "Bad dimensions" << std::endl;
            }

            Matrix matrix_out(this->num_rows, this->num_columns);

            /// Multiplying matrix a and b and storing in array
            for(uint8_t i = 0; i < matrix_out.num_rows; ++i) {
                for(uint8_t j = 0; j < matrix_out.num_columns; ++j) {
                    matrix_out.mx[i*matrix_out.num_columns + j] = this->mx[i*matrix.num_columns + j] - matrix.mx[i*matrix.num_columns + j];
                }
            }

            return matrix_out;
        }

        /**
         * @brief Element-wise assignment
         * @param matrix rvalue
         */
        void operator=(const Matrix &matrix) {
            for(uint8_t i = 0; i < this->num_rows; ++i) {
                for(uint8_t j = 0; j < this->num_columns; ++j) {
                    this->mx[i*this->num_columns + j] = matrix.mx[i*matrix.num_columns + j];
                }
            }
        }
        ///////////////////////////////////////////////////////////

        /**
         * @brief Display matrix using stdout
         */
        void show() {
            /// Displaying the multiplication of two matrix.
            std::cout << "Output Matrix: " << std::endl;
            for(uint8_t i = 0; i < this->num_rows; ++i) {
                for(int8_t j = 0; j < this->num_columns; ++j) {
                    std::cout << " " << this->mx[i*this->num_columns + j];
                    if(j == this->num_columns - 1) {
                        std::cout << std::endl;
                    }
                }
            }
        }
};

//Matrix control_signal(const double& time) {
//    return Matrix(1, 1);
//}

/// Task data
///////////////////////////////////////////////////////////
/// u(t) = 3*cos(0.1*t + 1)

/// A = [0 1 0; 0 0 1; -1.5 -5 -2]
/// B = [0 0 1]
/// C = [0.5 0 0]

/// x' = A*x + B*u;
/// y = C*x
///////////////////////////////////////////////////////////

int main() {
    /// Simulation parameters
    ///////////////////////////////////////////////////////////
    double time = 0;
    double time_step = 0.01;
    double simulation_time = 5;
    ///////////////////////////////////////////////////////////

    /// Initialise all matrices
    ///////////////////////////////////////////////////////////
    double _A[] = {0, 1, 0, 0, 0, 1, -1.5, -5, -2};
    Matrix A = Matrix(3, 3, _A);

    double _B[] = {0, 0, 1};
    Matrix B = Matrix(3, 1, _B);

    double _C[] = {0.5, 0, 0};
    Matrix C = Matrix(1, 3, _C);


    double _x[] = {1, 0, 0};
    Matrix x = Matrix(3, 1, _x);

    Matrix u = Matrix(1, 1);
    Matrix y = Matrix(3, 1);
    ///////////////////////////////////////////////////////////

    /// Simulation
    ///////////////////////////////////////////////////////////
    x = A*x + B*u;
    y = C*x;

    x.show();
    y.show();

    //while (time < simulation_time) {
    //    Matrix u = control_signal(time);

    //    time += time_step;
    //    std::cout<<"Time: "<< time << std::endl;
    //}
    ///////////////////////////////////////////////////////////

    return 0;
}

#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <inttypes.h>

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

            if(N != matrix.num_columns && M != matrix.num_rows) {
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
        ///////////////////////////////////////////////////////////
};

class Model
{
/// Attributes
///////////////////////////////////////////////////////////
private:
    Matrix Ad = Matrix(3, 3);
    Matrix Bd = Matrix(3, 1);
    Matrix Cd = Matrix(1, 3);
    Matrix x = Matrix(3, 1);
    Matrix y = Matrix(1, 1);

    const double OMEGA = 10;
    double x_prev_down = 0.84147 * -0.3; //-sin(1rad)*omega*amplitude
    double x_up = 0;
    double x_down = 0;
    double x_prev_up = 0.54*3; //cos(1)*amplitude
public:
    const double TIME_STEP = 0.01;
    const double SIMULATION_TIME = 20;
///////////////////////////////////////////////////////////

/// Methods
///////////////////////////////////////////////////////////
public:
    ///////////////////////////////////////////////////////////
    Model() {
        double _In[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        Matrix In(3, 3, _In);

        // double _A[] = {0, 1, 0, 0, 0, 1, -0.2, -1, -1};
        double _A[] = {0, 1, 0, 0, 0, 1, -1.5, -5, -2};
        Matrix A = Matrix(3, 3, _A);
        this->Ad = A;
        // this->Ad = A.expm();

        double _B[] = {0, 0, 1};
        Matrix B = Matrix(3, 1, _B);
        this->Bd = A.inverse()*(Ad - In)*B;

        // double _C[] = {2, 1, 0};
        double _C[] = {0.5, 0, 0};
        Matrix C = Matrix(1, 3, _C);
        this->Cd = C;

        double _x[] = {0, 0, 0};
        this->x = Matrix(3, 1, _x);

        this->y = Matrix(1, 1);
    }

    ~Model() {

    }
///////////////////////////////////////////////////////////
    
    double update(const double& input) {
        x = Ad*x + Bd*input;
        y = Cd*x;

        double output = this->y.getM();

        return output;
    }

    double control() {
        double T = this->TIME_STEP;

        x_up = x_prev_up + T * (x_prev_down / (T * T * OMEGA * OMEGA + 1) - T * OMEGA * OMEGA * x_prev_up / (T * T * OMEGA * OMEGA + 1));
        x_down = x_prev_down - T * OMEGA * OMEGA * (x_prev_up / (T * T * OMEGA * OMEGA + 1) + T * x_prev_down / (T * T * OMEGA * OMEGA + 1));
    
        x_prev_down = x_down;
        x_prev_up = x_up;

        double signal = x_up;

        return signal;
    }
///////////////////////////////////////////////////////////
};
#endif // MODEL_H

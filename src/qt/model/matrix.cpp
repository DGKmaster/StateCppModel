#include "matrix.h"

Matrix::Matrix(): 
    num_rows(0),
    num_columns(0) {
    this->mx = nullptr;
}

Matrix::~Matrix() {
    delete[] this->mx;
}

Matrix::Matrix(
    const uint16_t& r, 
    const uint16_t& c, 
    const double m[]):
    num_rows(r),
    num_columns(c),
    mx(new double[num_rows * num_columns]) {
    for(uint16_t i = 0; i < num_rows; i++) {
        for(uint16_t j = 0; j < num_columns; j++) {
            this->mx[i*num_columns + j] = 
            m[i*num_columns + j];
        }
    }
}

Matrix::Matrix(const uint16_t& r, const uint16_t& c):
    num_rows(r),
    num_columns(c),
    mx(new double[num_rows * num_columns]) {
    for(uint16_t i = 0; i < num_rows; i++) {
        for(uint16_t j = 0; j < num_columns; j++) {
            this->mx[i*num_columns + j] = 0;
        }
    }
}

double Matrix::getElement(const uint16_t& i) const {
    double output = this->mx[i];

    return output;
}

void Matrix::setElement(
    const uint16_t i, 
    const double& element) {
    this->mx[i] = element;
}

Matrix Matrix::operator*(const Matrix &matrix) const {
    if(this->num_columns != matrix.num_rows) {
        std::cout << "Bad dimensions =" << std::endl;
    }

    Matrix matrix_out(
        this->num_rows, 
        matrix.num_columns);

    /// Multiplying matrix a and b and storing in array
    for(uint16_t i = 0; i < this->num_rows; i++) {
        for(uint16_t j = 0; 
        j < matrix.num_columns; 
        j++) {
            for(uint16_t k = 0; 
            k < this->num_columns; 
            k++) {
                matrix_out.mx[
                    i*matrix_out.num_columns + j] 
                += this->mx[i*this->num_columns + k] 
                * matrix.mx[k*matrix.num_columns + j];
            }
        }
    }

    return matrix_out;
}

Matrix Matrix::operator*(const double &value) const {
    const uint16_t N = this->num_rows;
    const uint16_t M = this->num_columns;

    Matrix matrix_out(N, M);

    for(uint16_t i = 0; i < N; i++) {
        for(uint16_t j = 0; 
        j < M; 
        j++) {
            matrix_out.mx[i*M + j] = 
            this->mx[i*M + j] * value;
        }
    }

    return matrix_out;
}

Matrix Matrix::operator/(const double &value) const {
    const uint16_t N = this->num_rows;
    const uint16_t M = this->num_columns;

    Matrix matrix_out(N, M);

    for(uint16_t i = 0; 
    i < N; i++) {
        for(uint16_t j = 0; 
        j < M; j++) {
            matrix_out.mx[i*M + j] = 
            this->mx[i*M + j] / value;
        }
    }

    return matrix_out;
}

Matrix Matrix::operator+(const Matrix &matrix) const {
    const uint16_t N = this->num_rows;
    const uint16_t M = this->num_columns;

    if(N != matrix.num_rows && 
    M != matrix.num_columns) {
        std::cout << "Bad dimensions +" << std::endl;
    }

    Matrix matrix_out(N, M);

    for(uint16_t i = 0; i < N; i++) {
        for(uint16_t j = 0; j < M; j++) {
            matrix_out.mx[i*M + j] = 
            this->mx[i*M + j] + matrix.mx[i*M + j];
        }
    }

    return matrix_out;
}


Matrix Matrix::operator-(const Matrix &matrix) const {
    const uint16_t N = this->num_rows;
    const uint16_t M = this->num_columns;

    if(N != matrix.num_rows && 
    M != matrix.num_columns) {
        std::cout << "Bad dimensions -" << std::endl;
    }

    Matrix matrix_out(N, M);

    /// Multiplying matrix a and b and storing in array
    for(uint16_t i = 0; i < N; i++) {
        for(uint16_t j = 0; j < M; j++) {
            matrix_out.mx[i*M + j] = 
            this->mx[i*M + j] - matrix.mx[i*M + j];
        }
    }

    return matrix_out;
}

void Matrix::operator=(const Matrix &matrix) {
    const uint16_t N = this->num_rows;
    const uint16_t M = this->num_columns;

    for(uint16_t i = 0; i < N; i++) {
        for(uint16_t j = 0; j < M; j++) {
            this->mx[i*M + j] = matrix.mx[i*M + j];
        }
    }
}

void Matrix::show(const std::string& text_out) const {
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


Matrix Matrix::getCofactor(
    const uint16_t& p, 
    const uint16_t& q, 
    const int16_t& n) {
    const uint16_t N = this->num_rows;

    Matrix matrix_out(N, N);

    uint16_t i = 0;
    uint16_t j = 0;

    /// Looping for each element of the matrix
    for (uint16_t row = 0; row < n; row++) {
        for (uint16_t col = 0; 
        col < n; col++) {
            /// Copying into temporary 
            /// matrix only those element
            /// which are not in given row and column
            if (row != p && col != q) {
                matrix_out.mx[i*N + j] = 
                this->mx[row*N + col];
                j++;

                /// Row is filled, so 
                /// increase row index and
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

double Matrix::determinant(const uint16_t& n) {
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

        D += _sign * 
        this->mx[0*N + i] * _temp.determinant(n - 1);

        /// terms are to be added with alternate sign
        _sign = -_sign;
    }

    return D;
}

Matrix Matrix::adjoint() {
    const uint16_t N = this->num_rows;

    Matrix matrix_out(N, N);

    /// temp is used to store cofactors of A[][]
    int16_t _sign = 1;
    Matrix _temp(N, N);

    for (uint16_t i = 0; i < N; i++) {
        for (uint16_t j = 0; j < N; j++) {
            /// Get cofactor of A[i][j]
            _temp = this->getCofactor(i, j, N);

            /// sign of adj[j][i] 
            /// positive if sum of row
            /// and column indexes is even.
            _sign = ((i + j) % 2 == 0)? 1: -1;

            /// Interchanging rows 
            /// and columns to get the
            /// transpose of the cofactor matrix
            matrix_out.mx[j*N + i] = 
            _sign*_temp.determinant(N - 1);
        }
    }

    return matrix_out;
}


Matrix Matrix::inverse() {
    const uint16_t N = this->num_rows;

    /// Find determinant of A[][]
    double det = this->determinant(N);

    if (det == 0.0) {
        std::cout << 
        "Singular matrix, can't find its inverse" 
        << std::endl;
    }
    Matrix matrix_out(N, N);

    /// Find adjoint
    Matrix adj(N, N);
    adj = this->adjoint();

    /// Find Inverse using formula 
    /// "inverse(A) = adj(A)/det(A)"
    for (uint16_t i = 0; i < N; i++) {
        for (uint16_t j = 0; j < N; j++) {
            matrix_out.mx[i*N + j] = 
            adj.mx[i*N + j] / det;
        }
    }

    return matrix_out;
}

Matrix Matrix::matExp(double dt) {
    const uint16_t N = this->num_rows;

    double _matrix_out[] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    Matrix matrix_out(N, N, _matrix_out);

    Matrix A_power(N, N);
    A_power = *this;
    double dt_power = dt;
    double factorial = 1;
    
    /// First step - just matrix_out
    /// Second step - matrix_out + ...
    matrix_out = matrix_out + 
    (A_power * dt_power) / factorial;

    /// Third step and following
    for (uint16_t k = 2; k < 4; k++) {
        A_power = A_power * *this;
        dt_power = dt_power * dt;
        factorial = factorial * k;

        matrix_out = matrix_out 
        + (A_power * dt_power) / factorial;
    }

    return matrix_out;
}

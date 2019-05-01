#include <iostream>

/// u(t) = 3cos(0.1t + 1)

/// A = [0 1 0; 0 0 1; -1.5 -5 -2]
/// B = [0 0 1]
/// C = [0.5 0 0]

/// x' = A*x + B*u;
/// y = C*x

double control_signal(const double& time) {

}

int main() {
    double time = 0;
    double time_step = 0.01;
    double simulation_time = 5;

    while (time < simulation_time) {
        double u = control_signal(time);

        time += time_step;
        std::cout<<"Time: "<< time << std::endl;
    }

    return 0;
}

class Matrix {
    private:
        const uint16_t num_rows;
        const uint16_t num_columns;

        double* mx[][];

    public:
        Matrix(uint16_t num_rows, uint16_t num_columns):
            num_rows(num_rows),
            num_columns(num_columns),
            mx(new double[num_rows][num_columns]) {}

        Matrix matrix_multiply(Matrix matrix_1, Matrix matrix_2) {
            Matrix matrix_out(matrix_1.num_rows, matrix_2.num_columns);

            /// Multiplying matrix a and b and storing in array
            for(size_t i = 0; i < matrix_1.num_rows; ++i) {
                for(size_t j = 0; j < matrix_2.num_columns; ++j) {
                    for(size_t k = 0; k < matrix_1.num_columns; ++k) {
                        matrix_out.mx[i][j] += matrix_1.mx[i][k] * matrix_2.mx[k][j];
                    }
                }
            }

            /// Displaying the multiplication of two matrix.
            std::cout << "Output Matrix: " << std::endl;
            for(size_t i = 0; i < matrix_out.num_rows; ++i) {
                for(size_t j = 0; j < matrix_out.num_columns; ++j) {
                    std::cout << " " << matrix_out.mx[i][j];
                    if(j == matrix_out.num_columns - 1) {
                        std::cout << std::endl;
                    }
                }
            }
            
            return matrix_out;
        }
};

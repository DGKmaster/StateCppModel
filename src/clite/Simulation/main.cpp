#include <iostream>

int main() {
    std::cout<<"Hello World!"<<std::endl;
    
    
    return 0;
}

/// u(t) = 3cos(0.1t + 1)
/// A = [0 1 0; 0 0 1; -1.5 -5 -2]
/// B = [0 0 1]
/// C = [0.5 0 0]

double matrix_multiply(double arg_1, double arg_2);

double control_signal(const double& time);
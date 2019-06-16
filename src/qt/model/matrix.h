#include <inttypes.h>
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
    //////////////////////////////////////
    /**
    * @brief Default constructor
    * @return Nothing
    */
    Matrix();

    /**
     * @brief Default destructor
     * @return Nothing
     */
    ~Matrix();

    /**
     * @brief Initialise matrix RxC 
     * with specified array values
     * @param r Number of rows
     * @param c Number of columns
     * @param m Matrix elements
     */
    Matrix(
        const uint16_t& r, 
        const uint16_t& c, 
        const double m[]);

    /**
     * @brief Initialise matrix RxC with zero values
     * @param r Number of rows
     * @param c Number of columns
     */
    Matrix(const uint16_t& r, const uint16_t& c);
    //////////////////////////////////////
    
    /// Accessors
    //////////////////////////////////////
    double getElement(const uint16_t& i) const;

    void setElement(
        const uint16_t i, 
        const double& element);
    //////////////////////////////////////

    /// Operators overload
    //////////////////////////////////////
    /**
     * @brief Element-wise multiplication
     * @param matrix Right operand
     * @return Result
     */
    Matrix operator*(const Matrix &matrix) const;

    Matrix operator*(const double &value) const;

    Matrix operator/(const double &value) const;

    /**
     * @brief Element-wise sum
     * @param matrix Right operand
     * @return Result
     */
    Matrix operator+(const Matrix &matrix) const;

    /**
     * @brief Element-wise subtraction
     * @param matrix Right operand
     * @return Result
     */
    Matrix operator-(const Matrix &matrix) const;

    /**
     * @brief Element-wise assignment
     * @param matrix rvalue
     */
    void operator=(const Matrix &matrix);
    //////////////////////////////////////

    /// Other functions
    //////////////////////////////////////
    /**
     * @brief Display matrix using stdout
     */
    void show(
        const std::string& text_out = 
        "Output matrix: ") const;

    /**
    * @brief Function to get cofactor of 
    * A[p][q] in temp[][]. 
    * n is current dimension of A[][]
    * @param A
    * @param temp
    * @param p
    * @param q
    * @param n
        */
    Matrix getCofactor(
        const uint16_t& p, 
        const uint16_t& q, 
        const int16_t& n);

    /**
    * @brief Recursive function for 
    * finding determinant of matrix. 
    * n is current dimension of A[][].
    * @param A
    * @param n
    * @return
    */
    double determinant(const uint16_t& n);

    /**
     * @brief Function to get adjoint
     * of A[N][N] in adj[N][N].
     * @param A
     * @param adj
     */
    Matrix adjoint();

    /**
    * @brief Function to calculate and store inverse,
    * returns false if matrix is singular
    */
    Matrix inverse();

    Matrix matExp(double dt);
    //////////////////////////////////////
};

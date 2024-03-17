/*
ICS Minor Project

Matrix Algebra Utility Program

1. Matrix Operations
    a. Addition/Subtraction
    b. Multiplication
    c. Determinant
    d. Adjoint
    e. Transpose
    f. Inverse
    g. Determine whether matrix is
            Symmetric/Skew Sym
            Row/Column
            Diagonol
            Upper/Lower Triangular
            Singular/Non-Singular
            Idempotent
            Nilpotent
            Orthogonal
            Involutory

2. Single Variable Polynomial Operations
    a. Addition
    b. Multiplication
    c. Estimate Real Roots (based on Newton's Method)
    d. Derivative
    e. Integration

3. Solve Linear Equations (Row-Echelon or Gaussian Elimination + Back Substitution Method)

*/

#include <stdio.h>
#include <math.h>

// Matrix Start

// Function to input matrix
void inpMat(int rows, int cols, double mat[rows][cols])
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            scanf("%lf", &mat[i][j]);
}

// Function to print matrix on terminal
void outMat(int rows, int cols, double mat[rows][cols])
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
            printf("%lf ", mat[i][j]);
        printf("\n");
    }
}

// Matrix Operation Calculations Functions Start

// Make Matrix 2 = Matrix 1
void equateMat(int rows, int cols, double mat1[rows][cols], double mat2[rows][cols])
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            mat2[i][j] = mat1[i][j];
}

// Add two matrices
void addMat(int rows, int cols, double mat1[rows][cols], double mat2[rows][cols], double add[rows][cols])
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            add[i][j] = mat1[i][j] + mat2[i][j];
}

// Multiply two matrices
void mulMat(int rows1, int cols1, double mat1[rows1][cols1], int rows2, int cols2, double mat2[rows2][cols2], double mul[rows1][cols2])
{
    for (int i = 0; i < rows1; i++)
    {
        for (int j = 0; j < cols2; j++)
        {
            mul[i][j] = 0;
            for (int k = 0; k < cols1; k++)
                mul[i][j] += mat1[i][k] * mat2[k][j];
        }
    }
}

// Calculate Determininant of a Square Matrix
double determinant(int n, double mat[n][n])
{
    double det = 0;
    if (n == 1)
        return mat[0][0];
    else if (n == 2)
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    else
    {
        for (int k = 0; k < n; k++)
        {
            double submatrix[n - 1][n - 1];
            for (int i = 1; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    if (j < k)
                        submatrix[i - 1][j] = mat[i][j];
                    else if (j > k)
                        submatrix[i - 1][j - 1] = mat[i][j];
                }
            det += (k % 2 == 0 ? 1 : -1) * mat[0][k] * determinant(n - 1, submatrix);
        }
    }
    return det;
}

// Find transpose of a matrix
void transpose(int rows, int cols, double mat[rows][cols], double trans[cols][rows])
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            trans[j][i] = mat[i][j];
}

// Find Adjoint of a matrix
void adjoint(int n, double mat[n][n], double adj[n][n])
{
    double temp[n][n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double submatrix[n - 1][n - 1];
            double O = (i + j) % 2 == 0 ? 1 : -1;
            for (int row = 0; row < n; row++)
            {
                for (int col = 0; col < n; col++)
                {
                    if (row != i && col != j)
                    {
                        int newRow = row < i ? row : row - 1;
                        int newCol = col < j ? col : col - 1;
                        submatrix[newRow][newCol] = mat[row][col];
                    }
                }
            }
            temp[i][j] = O * determinant(n - 1, submatrix);
        }
    }
    transpose(n, n, temp, adj);
}

// Invert a matrix
void inverse(int n, double mat[n][n], double inv[n][n])
{ // Inverse(A) = Adjoint(A)/determinant(A)
    double adj[n][n];
    adjoint(n, mat, adj);
    double det = determinant(n, mat);
    double inv_det = 1.0 / det;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adj[i][j] * inv_det;
}

// Matrix Calculations End
// Matrix Properties Start

// Find if a matrix is Symmetric, i.e. A = A^
int isSym(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;
    double trans[cols][rows];
    transpose(rows, cols, mat, trans);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (trans[i][j] != mat[i][j])
                return 0;

    return 1;
}
// Find if a matrix is Skew-Symmetric, i.e. A+A^ = O

int isSkewSym(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;
    double trans[cols][rows];
    transpose(rows, cols, mat, trans);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (trans[i][j] != -mat[i][j])
                return 0;

    return 1;
}

// Find if a matrix is Diagonal, i.e. Aij = 0 for j!=i
int isDiag(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (i != j && mat[i][j])
                return 0;

    return 1;
}

// Find if a matrix is Upper Triangular, i.e. Aij = 0 for j<i
int isUpTr(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (j < i && mat[i][j])
                return 0;

    return 1;
}

// Find if a matrix is Lower Triangular, i.e. Aij = 0 for j>i
int isLoTr(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (i < j && mat[i][j])
                return 0;

    return 1;
}

// Find if a matrix is Idempotent, i.e. A^2 = A
int isIdem(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;

    double sq[rows][cols];
    mulMat(rows, cols, mat, rows, cols, mat, sq);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (sq[i][j] != mat[i][j])
                return 0;

    return 1;
}

// Find if a matrix is Zero, i.e. A = O
int isZero(int rows, int cols, double mat[rows][cols])
{
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            if (mat[i][j] != 0)
                return 0;
    return 1;
}

// Find if a matrix is Nilpotent, i.e. A^k = O for k<n
int isNil(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;

    int k = 1;
    double temp[rows][cols], temp2[rows][cols];
    mulMat(rows, cols, mat, rows, cols, mat, temp);

    while (isZero(rows, cols, temp) == 0 && k <= rows)
    {
        mulMat(rows, cols, mat, rows, cols, temp, temp2);
        equateMat(rows, cols, temp, temp2);
        k++;
    }

    if (k >= rows)
        return 0;
    return 1;
}

// Find if a matrix is Orthogonal, i.e. A*A^ = I
int isOrt(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;

    double trans[cols][rows];
    transpose(rows, cols, mat, trans);

    double mul[rows][cols];
    mulMat(rows, cols, mat, cols, rows, trans, mul);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            if (i == j && mul[i][j] != 1)
                return 0;
            if (i != j && mul[i][j])
                return 0;
        }

    return 1;
}

// Find if a matrix is Involutory, i.e. A^2 = I
int isInv(int rows, int cols, double mat[rows][cols])
{
    if (rows != cols)
        return 0;

    double sq[rows][cols];
    mulMat(rows, cols, mat, rows, cols, mat, sq);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
        {
            if (i == j && sq[i][j] != 1)
                return 0;
            if (i != j && sq[i][j])
                return 0;
        }

    return 1;
}

// Matrix Properties End
// Matrix Operation Utility Functions Start

// Operation to call addition
void OPaddMat()
{
    system("cls");
    printf("Enter number of Rows and Columns of Matrices:\n");
    int rows, cols;
    scanf("%d %d", &rows, &cols);
    double mat1[rows][cols], mat2[rows][cols], add[rows][cols];
    printf("Enter Matrix 1\n");
    inpMat(rows, cols, mat1);
    printf("Enter Matrix 2\n");
    inpMat(rows, cols, mat2);

    addMat(rows, cols, mat1, mat2, add);

    printf("Added Matrix\n");
    outMat(rows, cols, add);
}

// Operation to call subtraction
void OPsubMat()
{
    system("cls");
    printf("Enter number of Rows and Columns of Matrices:\n");
    int rows, cols;
    scanf("%d %d", &rows, &cols);
    double mat1[rows][cols], mat2[rows][cols], add[rows][cols];
    printf("Enter Matrix 1\n");
    inpMat(rows, cols, mat1);
    printf("Enter Matrix 2\n");
    inpMat(rows, cols, mat2);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            mat2[i][j] *= -1;

    addMat(rows, cols, mat1, mat2, add);

    printf("Subtracted Matrix\n");
    outMat(rows, cols, add);
}

// Operation to call multiplication
void OPmulMat()
{
    system("cls");
    printf("Enter number of Rows and Columns of Matrix 1:\n");
    int rows1, cols1;
    scanf("%d %d", &rows1, &cols1);
    printf("Enter number of Rows and Columns of Matrix 2:\n");
    int rows2, cols2;
    scanf("%d %d", &rows2, &cols2);

    if (cols1 != rows2)
    {
        printf("Cannot multiply matrices: invalid dimensions\n");
        return;
    }

    double mat1[rows1][cols1], mat2[rows2][cols2], mul[rows1][cols2];
    printf("Enter Matrix 1\n");
    inpMat(rows1, cols1, mat1);
    printf("Enter Matrix 2\n");
    inpMat(rows2, cols2, mat2);

    mulMat(rows1, cols1, mat1, rows2, cols2, mat2, mul);

    printf("Multiplied Matrix\n");
    outMat(rows1, cols2, mul);
}

// Operation to call determinant
void OPdetMat()
{
    system("cls");
    printf("Enter number of Rows and Columns of the Square Matrix:\n");
    int n;
    scanf("%d", &n);
    double mat[n][n];
    printf("Enter the Square Matrix\n");
    inpMat(n, n, mat);

    double det = determinant(n, mat);

    printf("Determinant of the Matrix: %lf\n", det);
}

// Operation to call adjoint
void OPadjMat()
{
    system("cls");
    printf("Enter number of Rows and Columns of the Square Matrix:\n");
    int n;
    scanf("%d", &n);
    double mat[n][n];
    printf("Enter the Square Matrix\n");
    inpMat(n, n, mat);

    double adj[n][n];
    adjoint(n, mat, adj);

    printf("Adjoint of the Matrix\n");
    outMat(n, n, adj);
}

// Operation to call transpose
void OPtransMat()
{
    system("cls");
    printf("Enter number of Rows and Columns of the Matrix:\n");
    int rows, cols;
    scanf("%d %d", &rows, &cols);
    double mat[rows][cols];
    printf("Enter the Matrix\n");
    inpMat(rows, cols, mat);

    double trans[cols][rows];
    transpose(rows, cols, mat, trans);

    printf("Transpose of the Matrix\n");
    outMat(cols, rows, trans);
}

// Operation to call inversion
void OPinvMat()
{
    system("cls");
    printf("Enter the size of the Square Matrix:\n");
    int n;
    scanf("%d", &n);
    double mat[n][n];
    printf("Enter the Square Matrix\n");
    inpMat(n, n, mat);

    double inv[n][n];
    inverse(n, mat, inv);

    printf("Inverse of the Matrix\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%lf ", inv[i][j]);
        printf("\n");
    }
}

// Operation to find properties
void OPpropertiesMat()
{
    system("cls");
    printf("Enter number of Rows and Columns of the Matrix:\n");
    int rows, cols;
    scanf("%d %d", &rows, &cols);
    double mat[rows][cols];
    printf("Enter the Matrix\n");
    inpMat(rows, cols, mat);

    if (isSym(rows, cols, mat))
        printf("Matrix is Symmetric\n"); // A=A^
    if (isSkewSym(rows, cols, mat))
        printf("Matrix is Skew-Symmetric\n"); // A=-A^
    if (isDiag(rows, cols, mat))
        printf("Matrix is Diagnol\n"); // Aij=0; i!=j;
    if (isUpTr(rows, cols, mat))
        printf("Matrix is Upper Triangular\n"); // Aij=0; j<i;
    if (isLoTr(rows, cols, mat))
        printf("Matrix is Lower Triangular\n"); // Aij=0; i<j
    if (rows == cols && determinant(rows, mat) == 0)
        printf("Matrix is Singular\n"); // detA=0
    else if (rows == cols && determinant(rows, mat))
        printf("Matrix is Non-Singular\n"); // detA!=0
    if (isIdem(rows, cols, mat))
        printf("Matrix is Idempotent\n"); // A^2=A
    if (isNil(rows, cols, mat))
        printf("Matrix is Nilpotent\n"); // A^k=0, k<n
    if (isOrt(rows, cols, mat))
        printf("Matrix is Orthogonal\n"); // A*A^=I
    if (isInv(rows, cols, mat))
        printf("Matrix is Involutory\n"); // A^2 = I
}

// Matrix OP end
// Matrix end

// Polynomial Start

// Function to input Polynomial
void inpPol(int d, double pol[d + 1])
{
    for (int i = d; i >= 0; i--)
        scanf("%lf", &pol[i]);
}

// Function to print Polynomial on terminal
void outPol(int d, double pol[d + 1])
{
    for (int j = 0; j < d + 1; j++)
    {
        if (j < d)
            printf("%lfx^%d + ", pol[j], j);
        else
            printf("%lfx^%d\n", pol[j], d);
    }
}

// Polynomial Operation Calculations Start

// Function to make an array zero
void zeroArr(int rows, double mat[rows])
{
    for (int i = 0; i < rows; i++)
        mat[i] = 0;
}

// Function to add two polynomials
void addPol(int d1, double pol1[d1 + 1], int d2, double pol2[d2 + 1], int d_res, double res_pol[d_res + 1])
{
    for (int i = 0; i < d_res; i++)
    {
        if (d_res > d2)
        {
            for (int i = 0; i < d2 + 1; i++)
            {
                res_pol[i] = pol1[i] + pol2[i];
            }
            for (int i = d2 + 1; i < d_res + 1; i++)
            {
                res_pol[i] = pol1[i];
            }
        }
        else
        {
            for (int i = 0; i < d1 + 1; i++)
            {
                res_pol[i] = pol1[i] + pol2[i];
            }
            for (int i = d1 + 1; i < d_res + 1; i++)
            {
                res_pol[i] = pol2[i];
            }
        }
    }
}

// Function to divide polynomial by linear function
void divPol(int d, double poly[d + 1], int r, double q[d])
{
    q[0] = poly[0];
    for (int i = 1; i <= d; i++)
    {
        q[i] = (q[i - 1] * r) + poly[i];
    }
}

// Function to multiply two polynomials
void mulPol(int deg1, int deg2, double coef1[deg1 + 1], double coef2[deg2 + 2], double mul[deg1 + deg2 + 1])
{

    zeroArr(deg1 + deg2 + 1, mul);
    for (int i = 0; i < deg1 + 1; i++)
    {
        for (int j = 0; j < deg2 + 1; j++)
        {
            mul[i + j] += coef1[i] * coef2[j];
        }
    }
}

// Function to differentiate a polynomial
void derPol(int d, double pol[d + 1], double der[d])
{
    for (int i = 0; i < d; i++)
    {
        der[i] = (i + 1) * pol[i + 1];
    }
}

// Function to indefintely integrate a polynomial
void intPol(int d, double pol[d + 1], double inte[d + 2])
{
    for (int i = 0; i < d + 1; i++) // Fix: start from 0
    {
        inte[i + 1] = pol[i] / (i + 1); // Fix: divide by i instead of i-1
    }
}

// Function to reverse a matrix
void reverse(int n, double arr[n])
{
    int left = 0;
    int right = n - 1;

    while (left < right)
    {
        double temp = arr[left];
        arr[left] = arr[right];
        arr[right] = temp;
        left++;
        right--;
    }
}

// Second Function to Print a polynomial
void print_polynomial(int degree, double coefficients[])
{
    printf("Polynomial: ");
    for (int i = degree; i >= 0; i--)
    {
        if (coefficients[i] != 0)
        {
            if (i == 0)
            {
                printf("%lfx^%d ", coefficients[i], i);
            }
            else
            {
                printf("%lfx^%d + ", coefficients[i], i);
            }
        }
    }
    printf("\n");
}

// Function to evaluate the value of the polynomial at a given point x
double evaluate_polynomial(int degree, double coefficients[], double x)
{
    double result = 0.0;
    for (int i = degree; i >= 0; i--)
    {
        result += coefficients[i] * pow(x, i);
    }
    return result;
}

// Function to evaluate the derivative of the polynomial at a given point x
double evaluate_derivative(int degree, double coefficients[], double x)
{
    double result = 0.0;
    for (int i = degree; i > 0; i--)
    {
        result += i * coefficients[i] * pow(x, i - 1);
    }
    return result;
}

// Function to perform synthetic division of a polynomial by a linear polynomial (x - r)
void synthetic_division(int degree, double coefficients[], double r, double q[])
{
    double temp[100 + 1]; // Temporary array to hold intermediate results
    for (int i = 0; i <= degree; i++)
    {
        temp[i] = coefficients[i];
    }
    q[0] = temp[degree];
    for (int i = degree - 1; i >= 0; i--)
    {
        q[i] = temp[i + 1] + q[i + 1] * r;
    }
    // print_polynomial(degree - 1, q);
}

// Function to estimate the root of a polynomial using Newton's method
double estimate_root(int degree, double coefficients[], double epsilon)
{
    // Automatically taking initial guess as 0
    int it = 0;
    double x = 3.141592653598979;
    double f_x, f_prime_x;
    do
    {
        f_x = evaluate_polynomial(degree, coefficients, x);
        f_prime_x = evaluate_derivative(degree, coefficients, x);
        x = x - f_x / f_prime_x; // Newton's Theorem
        it++;
    } while (fabs(f_x) > epsilon && it < 5000 && fabs(f_prime_x) > epsilon);

    if (fabs(f_x) < epsilon)
        return x;
    else
        return 3.141592653598979;
}

// Operation to call for solving
void OPsolvePol()
{
    system("cls");

    int degree, D, it = 0;

    printf("Enter the degree of the polynomial: ");
    scanf("%d", &degree);
    D = degree;
    double coefficients[degree + 1];

    for (int i = degree; i >= 0; i--)
    {
        printf("Coefficient of x^%d\n = ", i);
        scanf("%lf", &coefficients[i]);
    }

    print_polynomial(degree, coefficients);

    double roots[degree]; // Array to store roots
    int num_roots = 0;

    while (num_roots < D)
    {
        double root = estimate_root(degree, coefficients, 0.000000000001);

        if (root != 3.141592653598979)
            roots[num_roots++] = root;
        else if (root == 3.141592653598979)
        {
            // printf("No More Real Root Found\n");
            break;
        }
        // Divide the polynomial by (x - root)
        double q[degree];
        synthetic_division(degree, coefficients, root, q);

        for (int i = 0; i < degree; i++)
        {
            coefficients[i] = q[i];
        }
        // print_polynomial(degree, coefficients);
        degree--;
    }

    if (num_roots == 0)
    {
        printf("No real roots of the polynomial");
    }
    else
    {
        printf("Roots of the polynomial are: ");
        for (int i = 0; i < num_roots; i++)
        {
            printf("%lf ", roots[i]);
        }
        printf("\n");
    }
}

// Poly Op. Calc End
// Polynomial Operation Utility Functions Start

// Function to call integration
void OPintPol()
{
    system("cls");
    printf("Enter degree of Polynomial\n");
    int d;
    scanf("%d", &d);
    printf("Enter the Polynomial\n");
    double pol[d + 1];
    inpPol(d, pol);

    double inte[d + 2];

    intPol(d, pol, inte);

    printf("C");
    for (int i = 1; i < d + 2; i++)
    {
        printf("+%lfx^%d", inte[i], i);
    }
}

// Function to call Differentiation
void OPderPol()
{
    system("cls");
    printf("Enter degree of Polynomial\n");
    int d;
    scanf("%d", &d);
    printf("Enter the Polynomial\n");
    double pol[d + 1];
    inpPol(d, pol);

    double der[d];
    derPol(d, pol, der);

    outPol(d - 1, der);
}

// Function to call division
void OPdivPol()
{
    system("cls");
    printf("Enter degree of divident Polyonmial\n");
    int d;
    scanf("%d", &d);
    printf("Enter Divident Polynomiial\n");
    double poly[d + 1];
    for (int i = 0; i <= d; i++)
    {
        printf("Coefficient of x^%d\n = ", d - i);
        scanf("%lf", &poly[i]);
    }
    int r;
    printf("Enter value of constt 'r' in divisor (x-r)\n");
    scanf("%d", &r);
    double q[d], rem;

    divPol(d, poly, r, q);
    rem = q[d - 1] * r + poly[d];
    outPol(d - 1, q);
}

// Function to call Polynomial Multiplication
void OPmulPol()
{
    system("cls");
    printf("Enter degree of Polynomial 1\n");
    int d1;
    scanf("%d", &d1);
    printf("Enter Polynomial 1\n");
    double pol1[d1 + 1];
    inpPol(d1, pol1);

    printf("Enter degree of Polynomial 2\n");
    int d2;
    scanf("%d", &d2);
    printf("Enter Polynomial 2\n");
    double pol2[d2 + 1];
    inpPol(d2, pol2);

    double mul[d1 + d2 + 1];
    mulPol(d1, d2, pol1, pol2, mul);

    outPol(d1 + d2, mul);
}

// Function to call Polynomial Addition
void OPaddPol()
{
    system("cls");
    printf("Enter degree of Polynomial 1\n");
    int d1;
    scanf(" %d", &d1);
    printf("Enter Polynomial 1\n");
    double pol1[d1 + 1];
    inpPol(d1, pol1);
    outPol(d1, pol1);
    printf("Enter degree of Polynomial 2\n");
    int d2;
    scanf("%d", &d2);
    printf("Enter Polynomial 2\n");
    double pol2[d2 + 1];
    inpPol(d2, pol2);
    outPol(d2, pol2);
    int d_res;
    if (d1 > d2)
        d_res = d1;
    else
        d_res = d2;
    double res_pol[d_res + 1];
    addPol(d1, pol1, d2, pol2, d_res, res_pol);

    outPol(d_res, res_pol);

    return;
}
// Polynomial OP End
// Polynomial End

// Linear Equations Start

void back_subs(double eq[][101], double x[], int n)
{
    int i, j;

    for (i = n - 1; i >= 0; i--)
    {
        x[i] = eq[i][n];
        for (j = i + 1; j < n; j++)
            x[i] -= eq[i][j] * x[j];
        x[i] /= eq[i][i];
    }
}

void gauss_elim(double eq[][101], int n, int *solvable)
{
    int i, j, k;
    double factor;

    *solvable = 1;

    for (i = 0; i < n; i++)
    {
        if (eq[i][i] == 0)
        {
            *solvable = 0;
            return;
        }
        for (j = i + 1; j < n; j++)
        {
            factor = eq[j][i] / eq[i][i];
            for (k = i; k <= n; k++)
                eq[j][k] -= factor * eq[i][k];
        }
    }
}

void OPlin()
{
    system("cls");
    int n, solvable;
    double eq[100][101];

    printf("Enter the number of equations (maximum 100): ");
    scanf("%d", &n);

    if (n > 100)
    {
        printf("Number of equations exceeds the maximum allowed limit.\n");
        return;
    }

    printf("Enter the coefficients of the equations:\n");
    for (int i = 0; i < n; i++)
    {
        printf("Equation %d:\n", i + 1);
        for (int j = 0; j < n + 1; j++)
        {
            printf("Coefficient %d: ", j + 1);
            scanf("%lf", &eq[i][j]);
        }
    }

    double soln[100] = {0};

    gauss_elim(eq, n, &solvable);

    if (solvable == 0)
    {
        printf("The system of equations is singular or inconsistent.\n");
        return;
    }

    back_subs(eq, soln, n);

    printf("Solution:\n");
    for (int i = 0; i < n; i++)
        printf("x[%d] = %lf\n", i + 1, soln[i]);

    return;
}

// Linear Equatons End

// Main Functions Start

void OPmat()
{
    system("cls");
    printf("Choose Matrix Operation\n");
    printf("For Addition press A\n");
    printf("For Subtraction press S\n");
    printf("For Multiplication press M\n");
    printf("For Determinant press D\n");
    printf("For Adjoint press F\n");
    printf("For Transpose Matrix press T\n");
    printf("For Inverting Matrix press I\n");
    printf("For Finding Matrix Properties press P\n");

    char op; // User input of Matrix Operation
    scanf(" %c", &op);

    switch (op)
    {
    case 'A':
        OPaddMat();
        break;

    case 'S':
        OPsubMat();
        break;

    case 'M':
        OPmulMat();
        break;

    case 'D':
        OPdetMat();
        break;

    case 'F':
        OPadjMat();
        break;

    case 'T':
        OPtransMat();
        break;

    case 'I':
        OPinvMat();
        break;

    case 'P':
        OPpropertiesMat();
        break;

    default:
        printf("Invalid operation\n");
        break;
    }

    return;
}

void OPpol()
{
    system("cls");
    printf("Choose Polynomial Operation\n");
    printf("For adding two polynomials, press A\n");
    printf("For multiplying two polynomials, press M\n");
    printf("For Synthetic Division, press S\n");
    printf("For differentiating a polynomial, press D\n");
    printf("For integrating a polynomial, press I\n");
    printf("For finding roots, press R\n");

    char op_pol;
    scanf(" %c", &op_pol);

    switch (op_pol)
    {
    case 'A':
        OPaddPol();
        break;

    case 'M':
        OPmulPol();
        break;

    case 'S':
        OPdivPol();
        break;

    case 'D':
        OPderPol();
        break;

    case 'I':
        OPintPol();
        break;

    case 'R':
        OPsolvePol();
        break;

    default:
        break;
    }
}

int main()
{
    system("cls");
    printf("Welcome to Linear Algebraic Utility tool\n");
    printf("For Matrix Operations/Calculations, press M\n");
    printf("For Polynomial Operations/Calculations, press P\n");
    printf("For Solving Linear Equation, press L\n");
    char op_main;
    scanf("%c", &op_main);

    switch (op_main)
    {
    case 'M':
        OPmat();
        break;

    case 'P':
        OPpol();
        break;

    case 'L':
        OPlin();
        break;

    default:
        break;
    }
    printf("Press a key to End Program\n");
    scanf("%d", &op_main);
    return 0;
}
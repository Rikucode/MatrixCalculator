#include <vector>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <cmath>

#ifndef MATRIXCALCULATOR_H
#define MATRIXCALCULATOR_H

#endif

/* Определитель матрицы +

Нахождение обратной матрицы

Транспонирование матрицы*/

class Vector {
public:
    int dim;
    bool isRow;
    std::vector<double> coordinates;
    Vector(double x, double y, double z) {
        isRow = true;
        dim = 3;
        coordinates.resize(dim);
        this->coordinates[0] = x;
        this->coordinates[1] = y;
        this->coordinates[2] = z;
    }

    Vector(double x, double y) {
        isRow = true;
        dim = 2;
        coordinates.resize(dim);
        this->coordinates[0] = x;
        this->coordinates[1] = y;
    }

    Vector (std::vector<double> coordinates) {
        isRow = true;
        dim = coordinates.size();
        for (int i = 0; i < dim; ++i) {
            this->coordinates.push_back(coordinates[i]);
        }
    }

    Vector operator() (Vector vector1, Vector vector2) {
        if (vector1.dim == 3 && vector2.dim == 3) {
            coordinates[0] = vector1.coordinates[1] * vector2.coordinates[2] - vector2.coordinates[1] * vector1.coordinates[2];
            coordinates[1] = vector2.coordinates[0] * vector1.coordinates[2] - vector1.coordinates[0] * vector2.coordinates[2];
            coordinates[2] = vector1.coordinates[0] * vector2.coordinates[1] - vector2.coordinates[0] * vector1.coordinates[1];
        }
        else {
            std::cerr << "Vectors' dimensions must be equal 3!";
            std::exit(-1);
        }
        Vector vector(coordinates[0],coordinates[1],coordinates[2]);
        return vector;
    }
};


Vector operator+ (Vector vector1, Vector vector2) {
    if (vector1.dim != vector2.dim) {
        std::cerr << "Different dimensions!";
       std::exit(-1);
    } else {
        for (int i = 0; i < vector1.dim; ++i) {
            vector1.coordinates[i] += vector2.coordinates[i];
        }
    }
    return vector1;
}

Vector operator- (Vector vector1, Vector vector2) {
    if (vector1.dim != vector2.dim) {
        std::cerr << "Different dimensions!";
        std::exit(-1);
    } else {
        for (int i = 0; i < vector1.dim; ++i) {
            vector1.coordinates[i] -= vector2.coordinates[i];
        }
    }
    return vector1;
}

Vector operator* (double scalar, Vector vector) {
    for (int i = 0; i < vector.dim; ++i) {
        vector.coordinates[i] *= scalar;
    }
    return vector;
}

Vector operator* (Vector vector, double scalar) {
    for (int i = 0; i < vector.dim; ++i) {
        vector.coordinates[i] *= scalar;
    }
    return vector;
}

double operator* (Vector vector1, Vector vector2) {
    double sum = 0;
    if (vector1.dim != vector2.dim) {
        std::cerr << "Different dimensions!";
        std::exit(-1);
    }  else for (int i = 0; i < vector1.dim; ++i) {
        sum += vector1.coordinates[i] * vector2.coordinates[i];
    }
    return sum;
}


class Matrix {
public:
    int n, m;
    std::vector<std::vector<double>> values;
    Matrix(std::vector<std::vector<double>> values) {
        this->n = values.size();
        this->m = values[0].size();
        std::vector<double> temp_values;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                temp_values.push_back(values[i][j]);
            }
            this->values.push_back(temp_values);
            temp_values.clear();
        }
    }
};

Matrix operator + (Matrix matrix1,Matrix matrix2) {
    if (matrix1.n != matrix2.n || matrix1.m != matrix2.m) {
        std::cerr << "Different dimensions!";
        std::exit(-1);
    }
    for (int i = 0; i < matrix1.n; ++i) {
        for (int j = 0; j < matrix1.m; ++j) {
            matrix1.values[i][j]+=matrix2.values[i][j];
        }
    }
    return matrix1;
}

Matrix operator - (Matrix matrix1,Matrix matrix2) {
    if (matrix1.n != matrix2.n || matrix1.m != matrix2.m) {
        std::cerr << "Different dimensions!";
        std::exit(-1);
    }
    for (int i = 0; i < matrix1.n; ++i) {
        for (int j = 0; j < matrix1.m; ++j) {
            matrix1.values[i][j]-=matrix2.values[i][j];
        }
    }
    return matrix1;
}

Matrix operator * (Matrix matrix1, Matrix matrix2) {
    if (matrix1.m != matrix2.n) {
        std::cerr << "Different dimensions!";
        std::exit(-1);
    }
    std::vector<std::vector<double>> matrix_values;
    std::vector<double> tempValues;
    matrix_values.resize(matrix1.n);
    for (int i = 0; i < matrix1.n; ++i) {
        matrix_values[i].resize(matrix2.m);
    }
    for (int i = 0; i < matrix1.n; ++i) {
        for (int j = 0; j < matrix2.m; ++j) {
            matrix_values [i][j] = 0;
            for (int k = 0; k < matrix1.m; ++k) {
                matrix_values[i][j] += matrix1.values[i][k] * matrix2.values[k][j];
            }
        }
    }
    Matrix matrix(matrix_values);
    return matrix;
}

Matrix operator * (Matrix matrix, double scalar){
    for (int i = 0; i < matrix.n; ++i) {
        for (int j = 0; j < matrix.m; ++j) {
            matrix.values[i][j] *= scalar;
        }
    }
    return matrix;
}

Matrix operator * (double scalar, Matrix matrix){
    for (int i = 0; i < matrix.n; ++i) {
        for (int j = 0; j < matrix.m; ++j) {
            matrix.values[i][j] *= scalar;
        }
    }
    return matrix;
}

void transpMatrix (Matrix& matrix) {
    std::vector<std::vector<double>> tmpValues;
    tmpValues.resize(matrix.m);
    for (int i = 0; i < matrix.m; ++i) {
        tmpValues[i].resize(matrix.n);
    }
    for (int i = 0; i < matrix.n; ++i) {
        for (int j = 0; j < matrix.m; ++j) {
            tmpValues[j][i] = matrix.values[i][j];
        }
    }
    Matrix matrix1(tmpValues);
    matrix = matrix1;
}

Matrix vectorToMatrix (Vector vector) {
    std::vector<std::vector<double>> tmpVector;
    if (vector.isRow) {
        tmpVector.push_back(vector.coordinates);
        Matrix matrix(tmpVector);
        return matrix;
    } else {
        tmpVector.resize(vector.coordinates.size());
        for (int i = 0; i < vector.coordinates.size(); i++) {
            tmpVector[i].resize(1);
        }
        for (int i = 0; i < vector.coordinates.size(); i++){
            tmpVector[i][0] = vector.coordinates[i];
        }
        Matrix matrix(tmpVector);
        return matrix;
    }
}

void transpVector (Vector& vector) {
    if (vector.isRow){
        vector.isRow = false;
    } else {
        vector.isRow = true;
    }
}

Matrix operator* (Vector vector, Matrix matrix) {
    Matrix tmpMatrix = vectorToMatrix(vector);
    matrix = tmpMatrix * matrix;
    return matrix;
}

Matrix operator* (Matrix matrix, Vector vector) {
    Matrix tmpMatrix = vectorToMatrix(vector);
    matrix = matrix * tmpMatrix;
    return matrix;
}

int search(Matrix matrix, double target,
           bool match, unsigned int& uI, unsigned int& uJ, unsigned int starti, unsigned int startj) {
    if ((!matrix.m) || (!matrix.n)) return 0;
    if ((starti >= matrix.n) || (startj >= matrix.m)) return 0;
    for (unsigned int i = starti; i < matrix.n; i++)
        for (unsigned int j = startj; j < matrix.m; j++) {
            if (match) {
                if (matrix.values[i][i] == target) {
                    uI = i; uJ = j; return 1;
                }
            }
            else if (matrix.values[i][j] != target) {
                uI = i; uJ = j; return 1;
            }
        }
    return 0;
}

void swaprows(Matrix matrix, int x1, int x2) {
    if ((!matrix.n) || (!matrix.m)) return;
    if ((x1 >= matrix.n) || (x2 >= matrix.n) || (x1 == x2)) return;
    double tmp;
    for (unsigned int x = 0; x < matrix.m; x++) {
        tmp = matrix.values[x1][x];
        matrix.values[x1][x] = matrix.values[x2][x];
        matrix.values[x2][x] = tmp;
    }
    return;
};

void swapcolumns(Matrix matrix, int x1, int x2) {
    if ((!matrix.n) || (!matrix.m)) return;
    if ((x1 >= matrix.m) || (x2 >= matrix.m) || (x1 == x2)) return;
    double tmp;
    for (unsigned int x = 0; x < matrix.n; x++) {
        tmp = matrix.values[x][x1];
        matrix.values[x][x1] = matrix.values[x][x2];
        matrix.values[x][x2] = tmp;
    }
    return;
};

double determinant(Matrix matrix) {
    if (matrix.m != matrix.n){
        std::cerr << "Matrix should be square";
        std::exit(-1);
    }
    unsigned int m = matrix.n;
    if (m == 0) return 0;
    if (m == 1) return matrix.values[0][0];
    if (m == 2) return (matrix.values[0][0] * matrix.values[1][1] - matrix.values[1][0] * matrix.values[0][1]);
    bool sign = false;
    double det = 1;
    double tmp;
    unsigned int x, y;
    for (unsigned int i = 0; i < m; i++) {
        if (matrix.values[i][i] == 0) {
            if (!search(matrix.values, 0, false, y, x, i, i)) return 0;
            if (i != y) {
                swaprows(matrix, i, y);
                sign = !sign;
            }
            if (i != x) {
                swapcolumns(matrix, i, x);
                sign = !sign;
            }
        }
        det *= matrix.values[i][i];
        tmp = matrix.values[i][i];
        for (x = i; x < m; x++) {
            matrix.values[i][x] = matrix.values[i][x] / tmp;
        }
        for (y = i + 1; y < matrix.n; y++) {
            tmp = matrix.values[y][i];
            for (x = i; x < m; x++)
                matrix.values[y][x] -= (matrix.values[i][x] * tmp);
        }
    }
    if (sign) return det * (-1);
    return det;
};

Matrix vectorsImp(Vector vector1, Vector vector2) {
    Matrix matrix1 = vectorToMatrix(vector1);
    Matrix matrix2 = vectorToMatrix(vector2);
    matrix1 = matrix1 * matrix2;
    return matrix1;
}

void inversion(Matrix& matrix)
{
    if (matrix.m != matrix.n){
        std::cerr << "Matrix should be square";
        std::exit(-1);
    }
    if (!determinant(matrix)){
        std::cerr << "Determinant mustn't be equal 0";
        std::exit(-1);
    }
    double temp;
    int N = matrix.n;

    std::vector<std::vector<double>> E;
    E.resize(N);

    for (int i = 0; i < N; i++)
        E[i].resize(N);

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            E[i][j] = 0.0;

            if (i == j)
                E[i][j] = 1.0;
        }

    for (int k = 0; k < N; k++)
    {
        temp = matrix.values[k][k];

        for (int j = 0; j < N; j++)
        {
            matrix.values[k][j] /= temp;
            E[k][j] /= temp;
        }

        for (int i = k + 1; i < N; i++)
        {
            temp = matrix.values[i][k];

            for (int j = 0; j < N; j++)
            {
                matrix.values[i][j] -= matrix.values[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = matrix.values[i][k];

            for (int j = 0; j < N; j++)
            {
                matrix.values[i][j] -= matrix.values[k][j] * temp;
                E[i][j] -= E[k][j] * temp;
            }
        }
    }

    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            matrix.values[i][j] = E[i][j];
        }
    }
}
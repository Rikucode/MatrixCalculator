#include <vector>
#include <stdlib.h>
#include <iostream>
#ifndef MATRIXCALCULATOR_H
#define MATRIXCALCULATOR_H

#endif


class Vector {
public:
    int dim;
    double x, y, z;

    Vector(double x, double y, double z) {
        dim = 3;
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Vector(double x, double y) {
        dim = 2;
        this->x = x;
        this->y = y;
    }

    Vector operator() (Vector vector1, Vector vector2) {
        if (vector1.dim == 3 && vector2.dim == 3) {
            x = vector1.y * vector2.z - vector2.y * vector1.z;
            y = vector2.x * vector1.z - vector1.x * vector2.z;
            z = vector1.x * vector2.y - vector2.x * vector1.y;
        }
        else {
            std::cerr << "Vectors' dimensions must be equal 3!";
            std::exit(-1);
        }
        Vector vector(x,y,z);
        return vector;
    }
};


Vector operator+ (Vector vector1, Vector vector2) {
    if (vector1.dim != vector2.dim) {
        std::cerr << "Different dimensions!";
       std::exit(-1);
    } else if (vector1.dim == 3){
        vector1.x += vector2.x;
        vector1.y += vector2.y;
        vector1.z += vector2.z;
    } else if (vector1.dim == 2){
        vector1.x += vector2.x;
        vector1.y += vector2.y;
    }
    return vector1;
}

Vector operator- (Vector vector1, Vector vector2) {
    if (vector1.dim != vector2.dim) {
        std::cerr << "Different dimensions!";
        std::exit(-1);
    } else if (vector1.dim == 3){
        vector1.x -= vector2.x;
        vector1.y -= vector2.y;
        vector1.z -= vector2.z;
    } else if (vector1.dim == 2){
        vector1.x -= vector2.x;
        vector1.y -= vector2.y;
    }
    return vector1;
}

Vector operator* (double scalar, Vector vector) {
    if (vector.dim == 3){
        vector.x *= scalar;
        vector.y *= scalar;
        vector.z *= scalar;
    } else if (vector.dim == 2){
        vector.x *= scalar;
        vector.y *= scalar;
    }
    return vector;
}

Vector operator* (Vector vector, double scalar) {
    if (vector.dim == 3){
        vector.x *= scalar;
        vector.y *= scalar;
        vector.z *= scalar;
    } else if (vector.dim == 2){
        vector.x *= scalar;
        vector.y *= scalar;
    }
    return vector;
}

double operator* (Vector vector1, Vector vector2) {
    double sum = 0;
    if (vector1.dim != vector2.dim) {
        std::cerr << "Different dimensions!";
        std::exit(-1);
    } else if (vector1.dim == 2){
        sum = vector1.x * vector2.x + vector1.y * vector2.y;
    } else if (vector1.dim == 3){
        sum = vector1.x * vector2.x + vector1.y * vector2.y + vector1.z * vector2.z;
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
                matrix_values[i][j] += matrix1.values[i][k]*matrix2.values[k][j];
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


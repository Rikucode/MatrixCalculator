#include "MatrixCalculator.h"
#include <iostream>

using namespace std;

void v_out(Vector vector) {
    cout << "x: " << vector.x << "; y: " << vector.y << "; z: " << vector.z << endl;
}

void m_out(Matrix matrix) {
    for (int i = 0; i < matrix.n; ++i) {
        for (int j = 0; j < matrix.m; ++j) {
            cout << matrix.values[i][j] << ' ';
        }
        cout << endl;
    }
}

int main(){
//    vector<vector<double>> matrixV = {{1,2,2},{3,1,1}};
//    vector<vector<double>> matrix1 = {{4,2},{3,1},{1,5}};
//    vector<vector<double>> matrix2 = {{1,4,5},{32,5,0},{3,2,9}};
//    Matrix matrix(matrixV);
//    Matrix matrixA(matrix1);
//    matrix = matrixA * matrix;
//    m_out(matrix);
//    Vector vector1(2,2,3);
//    Vector vector2(3,4,5);
//    Vector vector(0,0,0);
//    vector (vector1,vector2);
//    v_out(vector);
    return 0;
}

#include "MatrixCalculator.h"
#include <iostream>
#include <iomanip>

using namespace std;

void v_out(Vector vector) {
    if (vector.isRow) {
        for (int i = 0; i < vector.dim; ++i) {
            cout << vector.coordinates[i] << ' ';
        }
        cout << endl;
    } else {
        for (int i = 0; i < vector.dim; ++i) {
            cout << vector.coordinates[i] << endl;
        }
    }
}

void m_out(Matrix matrix) {
    for (int i = 0; i < matrix.n; ++i) {
        for (int j = 0; j < matrix.m; ++j) {
            if (abs(matrix.values[i][j]) < 1.0e-10){
                cout << "0 ";
            } else {
                cout << matrix.values[i][j] << setprecision(10) << ' ';
            }
        }
        cout << endl;
    }
}

int main(){
    //Поле для экспериментов
    Vector vector1(2,3,4);
    Vector vector2(4,2,9);
    cout << "Vector's sum: ";
    v_out(vector1 + vector2);
    cout << "Vector's dif: ";
    v_out(vector1 - vector2);
    cout << "Vector and scalar product: ";
    v_out(vector1 * 3);
    transpVector(vector2);
    cout << "Vectors' scalar product (row on column): " << vector1*vector2 << endl;
    transpVector(vector2);
    cout << "Vectors' cross product: ";
    v_out(vector1(vector1,vector2));

    return 0;
}

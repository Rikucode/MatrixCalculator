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
    return 0;
}

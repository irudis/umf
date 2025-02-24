#include <iostream>
#include <iomanip>
#include <algorithm>

#include "data.h"

using namespace std;

int main() {
    setlocale(LC_ALL, "");

    diff_data data;

    data.input();
    data.points();
    data.apply_conditions();
    data.build_matrix();
    data.gauss_zeidel(1e-016, 1000, 1.1);
    cout << "Вектор:" << endl;
    for (int i = 0; i < data.nodes.size(); i++)
        if (!binary_search(data.fict.begin(), data.fict.end(), i))
            cout << scientific << setprecision(15) << data.u[i] << endl;
    cout << "Истинное значение:" << endl;
    for (int i = 0; i < data.nodes.size(); i++)
        if (!binary_search(data.fict.begin(), data.fict.end(), i))
            cout << scientific << setprecision(15) << data.U(data.nodes[i].first, data.nodes[i].second) << endl;
    cout << "Погрешность" << endl;
    for (int i = 0; i < data.nodes.size(); i++)
        if (!binary_search(data.fict.begin(), data.fict.end(), i))
            cout << scientific << setprecision(15) << abs(data.u[i] - data.U(data.nodes[i].first, data.nodes[i].second)) << endl;
    return 0;
}

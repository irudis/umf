#include <fstream>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <iomanip>

#include "data.h"

using namespace std;

double diff_data::func(double x, double y)  {
    return x * x + y * y - 4;
}

double diff_data::U(double x, double y) {
    return x * x + y * y;
}

void diff_data::input() {
    ifstream in("input/area.txt");
    in >> N;
    blocks.resize(N);
    blocks[0].resize(4);
    in >> blocks[0][0] >> blocks[0][1] >> blocks[0][2] >> blocks[0][3];
    X_min = blocks[0][0];
    X_max = blocks[0][1];
    Y_min = blocks[0][2];
    Y_max = blocks[0][3];
    for (int i = 1; i < N; i++) {
        blocks[i].resize(4);
        in >> blocks[i][0] >> blocks[i][1] >> blocks[i][2] >> blocks[i][3];
        if (X_min > blocks[i][0])
            X_min = blocks[i][0];
        if (X_max < blocks[i][1])
            X_max = blocks[i][1];
        if (Y_min > blocks[i][2])
            Y_min = blocks[i][2];
        if (Y_max < blocks[i][3])
            Y_max = blocks[i][3];
    }
    in.close();
    in.open("input/coords.txt");
    in >> nx >> ny >> kx >> ky >> lambda >> gamma;
    in.close();
}

void diff_data::points() {
    if (kx == 1)
        hx = (X_max - X_min) / nx;
    else {
        double sum = 0;
        for (int i = 0; i < nx; i++)
            sum += pow(kx, i);
        hx = (X_max - X_min) / sum;
    }

    if (ky == 1)
        hy = (Y_max - Y_min) / ny;
    else {
        double sum = 0;
        for (int i = 0; i < ny; i++)
            sum += pow(ky, i);
        hy = (Y_max - Y_min) / sum;
    }
    nodes.resize((ny + 1) * (nx + 1));
    for (int k = 0; k < (ny + 1) * (nx + 1); k++) {
        int row = k / (nx + 1);
        int col = k - row * (nx + 1);
        double sx = 0, sy = 0;
        for (int i = 0; i < col; i++)
            sx += hx * pow(kx, i);
        for (int i = 0; i < row; i++)
            sy += hy * pow(ky, i);
        if (col == nx)
            nodes[k].first = X_max;
        else
            nodes[k].first = X_min + sx;
        if (row == ny)
            nodes[k].second = Y_max;
        else
            nodes[k].second = Y_min + sy;
    }
    for (int i = 0; i < nodes.size(); i++) {
        int f = 0;
        for (int j = 0; j < N; j++)
            if (nodes[i].first <= blocks[j][1] && nodes[i].first >= blocks[j][0] && nodes[i].second <= blocks[j][3] && nodes[i].second >= blocks[j][2]) {
                f = 1;
                break;
            }
        if (!f) {
            fict.resize(fict.size() + 1);
            fict.back() = i;
        }
    }
}

// void diff_data::apply_conditions() {
//     int num;
//     ifstream in("bc2.txt");
//     in >> num;
//     for (int i = 0; i < num; i++) {
//         double f, s;
//         int bl, gr;
//         in >> bl >> gr;
//         switch (gr) {
//         case 1:
//             for (int j = 0; j < nodes.size(); j++)
//                 if (nodes[j].first >= blocks[bl - 1][0] && nodes[j].first <= blocks[bl - 1][1] && nodes[j].second == blocks[bl - 1][2]) {
//                     bc2.resize(bc2.size() + 1);
//                     bc2.back() = j;
//                 }
//             break;
//         case 2:
//             for (int j = 0; j < nodes.size(); j++)
//                 if (nodes[j].second >= blocks[bl - 1][2] && nodes[j].second <= blocks[bl - 1][3] && nodes[j].first == blocks[bl - 1][1]) {
//                     bc2.resize(bc2.size() + 1);
//                     bc2.back() = j;
//                 }
//             break;
//         case 3:
//             for (int j = 0; j < nodes.size(); j++)
//                 if (nodes[j].first >= blocks[bl - 1][0] && nodes[j].first <= blocks[bl - 1][1] && nodes[j].second == blocks[bl - 1][3]) {
//                     bc2.resize(bc2.size() + 1);
//                     bc2.back() = j;
//                 }
//             break;
//         case 4:
//             for (int j = 0; j < nodes.size(); j++)
//                 if (nodes[j].second >= blocks[bl - 1][2] && nodes[j].second <= blocks[bl - 1][3] && nodes[j].first == blocks[bl - 1][0]) {
//                     bc2.resize(bc2.size() + 1);
//                     bc2.back() = j;
//                 }
//             break;
//         }
//     }
//     sort(bc2.begin(), bc2.end());
//     auto last = unique(bc2.begin(), bc2.end());
//     bc2.erase(last, bc2.end());
//     in.close();
//     in.open("input/conditions_1.txt");
//     in >> num;
//     for (int i = 0; i < num; i++) {
//         double f, s;
//         int bl, gr;
//         in >> bl >> gr;
//         switch (gr) {
//         case 1:
//             for (int j = 0; j < nodes.size(); j++)
//                 if (nodes[j].first >= blocks[bl - 1][0] && nodes[j].first <= blocks[bl - 1][1] && nodes[j].second == blocks[bl - 1][2]) {
//                     bc1.resize(bc1.size() + 1);
//                     bc1.back() = j;
//                 }
//             break;
//         case 2:
//             for (int j = 0; j < nodes.size(); j++)
//                 if (nodes[j].second >= blocks[bl - 1][2] && nodes[j].second <= blocks[bl - 1][3] && nodes[j].first == blocks[bl - 1][1]) {
//                     bc1.resize(bc1.size() + 1);
//                     bc1.back() = j;
//                 }
//             break;
//         case 3:
//             for (int j = 0; j < nodes.size(); j++)
//                 if (nodes[j].first >= blocks[bl - 1][0] && nodes[j].first <= blocks[bl - 1][1] && nodes[j].second == blocks[bl - 1][3]) {
//                     bc1.resize(bc1.size() + 1);
//                     bc1.back() = j;
//                 }
//             break;
//         case 4:
//             for (int j = 0; j < nodes.size(); j++)
//                 if (nodes[j].second >= blocks[bl - 1][2] && nodes[j].second <= blocks[bl - 1][3] && nodes[j].first == blocks[bl - 1][0]) {
//                     bc1.resize(bc1.size() + 1);
//                     bc1.back() = j;
//                 }
//             break;
//         }
//     }
//     sort(bc1.begin(), bc1.end());
//     last = unique(bc1.begin(), bc1.end());
//     bc1.erase(last, bc1.end());
//     in.close();
//     for (int i = 0; i < bc2.size(); i++)
//         if (binary_search(bc1.begin(), bc1.end(), bc2[i])) {
//             bc2.erase(bc2.begin() + i);
//             i--;
//         }
// }

void diff_data::apply_conditions() {
    auto process = [this](const string& file, vector<double>& bc) {
        ifstream in(file);
        for (int num; in >> num;)
            while (num--) {
                int bl, gr; in >> bl >> gr;
                const auto& b = blocks[bl - 1];
                bool is_x; double val, lo, hi;
                switch (gr) {
                    case 1: is_x = false; val = b[2]; lo = b[0]; hi = b[1]; break;
                    case 2: is_x = true;  val = b[1]; lo = b[2]; hi = b[3]; break;
                    case 3: is_x = false; val = b[3]; lo = b[0]; hi = b[1]; break;
                    case 4: is_x = true;  val = b[0]; lo = b[2]; hi = b[3]; break;
                    default: continue;
                }
                for (int j = 0; j < nodes.size(); ++j) {
                    const auto& n = nodes[j];
                    if ((is_x  && n.first == val && n.second >= lo && n.second <= hi) ||
                        (!is_x && n.second == val && n.first >= lo && n.first <= hi))
                        bc.push_back(j);
                }
            }
        sort(bc.begin(), bc.end());
        bc.erase(unique(bc.begin(), bc.end()), bc.end());
    };

    process("bc2.txt", bc2);
    process("input/conditions_1.txt", bc1);

    bc2.erase(remove_if(bc2.begin(), bc2.end(), [&](int i) {
        return binary_search(bc1.begin(), bc1.end(), i);
    }), bc2.end());
}

void diff_data::build_matrix() {
    m = nx + 1;
    matrix.resize(5);
    matrix[2].resize(nodes.size());
    matrix[1].resize(nodes.size() - 1);
    matrix[3].resize(nodes.size() - 1);
    matrix[0].resize(nodes.size() - m);
    matrix[4].resize(nodes.size() - m);
    vec.resize(nodes.size());
    u.resize(nodes.size());
    vector<double> Hx;
    Hx.resize(nx);
    Hx.insert(Hx.begin(), hx);
    vector<double> Hy;
    Hy.resize(nx);
    Hy.insert(Hy.begin(), hy);
    for (int i = 0; i <= nx; i++)
        Hx[i] = hx * pow(kx, i);
    for (int i = 0; i <= ny; i++)
        Hy[i] = hy * pow(ky, i);
    for (int i = 0; i <= ny; i++)
        for (int j = 0; j <= nx; j++) {
            int k = i * m + j;
            if (binary_search(fict.begin(), fict.end(), k)) {
                matrix[2][k] = 1;
                vec[k] = 0;
            }
            else {
                if (i != 0 && i != ny && j != 0 && j != nx && !binary_search(fict.begin(), fict.end(), k + 1) && !binary_search(fict.begin(), fict.end(), k - 1) && !binary_search(fict.begin(), fict.end(), k + m) && !binary_search(fict.begin(), fict.end(), k - m) && !binary_search(fict.begin(), fict.end(), k + m + 1) && !binary_search(fict.begin(), fict.end(), k + m - 1) && !binary_search(fict.begin(), fict.end(), k - m + 1) && !binary_search(fict.begin(), fict.end(), k - m - 1)) {
                    matrix[2][k] = lambda * (2 / (Hx[j] * Hx[j - 1]) + 2 / (Hy[i] * Hy[i - 1])) + gamma;
                    matrix[1][k - 1] = -2 * lambda / (Hx[j - 1] * (Hx[j] + Hx[j - 1]));
                    matrix[3][k] = -2 * lambda / (Hx[j] * (Hx[j] + Hx[j - 1]));
                    matrix[0][k - m] = -2 * lambda / (Hy[i - 1] * (Hy[i] + Hy[i - 1]));
                    matrix[4][k] = -2 * lambda / (Hy[i] * (Hy[i] + Hy[i - 1]));
                    vec[k] = func(nodes[k].first, nodes[k].second);
                }
                else {
                    if (binary_search(bc1.begin(), bc1.end(), k)) {
                        matrix[2][k] = 1;
                        vec[k] = U(nodes[k].first, nodes[k].second);
                    }
                    else if (binary_search(bc2.begin(), bc2.end(), k)) {
                        if (i == 0 || binary_search(fict.begin(), fict.end(), k - m)) {
                            matrix[2][k] = -lambda / Hy[i];
                            matrix[3][k] = lambda / Hy[i];
                            vec[k] = -tetta(nodes[k].first, nodes[k].second);
                        }
                        else if (i == ny || binary_search(fict.begin(), fict.end(), k + m)) {
                            matrix[2][k] = lambda / Hy[i - 1];
                            matrix[1][k - 1] = -lambda / Hy[i - 1];
                            vec[k] = tetta(nodes[k].first, nodes[k].second);
                        }
                        else if (j == 0 || binary_search(fict.begin(), fict.end(), k - 1)) {
                            matrix[2][k] = -lambda / Hx[j];
                            matrix[3][k] = lambda / Hx[j];
                            vec[k] = -tetta(nodes[k].first, nodes[k].second);
                        }
                        else if (j == nx || binary_search(fict.begin(), fict.end(), k + 1)) {
                            matrix[2][k] = lambda / Hx[j - 1];
                            matrix[1][k - 1] = -lambda / Hx[j - 1];
                            vec[k] = tetta(nodes[k].first, nodes[k].second);
                        }
                    }
                }
            }
        }
}

double diff_data::tetta(double x, double y) {
    return 4 * y;
}

double diff_data::norm(int n, vector<double>& vec) {
    double sum = 0;
    for (int i = 0; i < n; i++)
        sum += vec[i] * vec[i];
    return sqrt(sum);
}

double diff_data::iteration(int i, int n, int m, double w, vector<vector<double>>& matrix, vector<double>& vec, vector<double>& x, vector<double>& xnext) {
    double sum = 0, a = 0;
    for (int j = 0; j < n; j++) {
        int r = j - i;
        if (r == 0) a = matrix[2][i]; 
        else if (r == 1) a = matrix[3][i];				
        else if (r == m) a = matrix[4][i];
        else if (r == -1) a = matrix[1][j];			
        else if (r == -m) a = matrix[0][j];
        else  a = 0;
        sum += a * x[j];
    }
    xnext[i] = x[i] + w / matrix[2][i] * (vec[i] - sum);
    return sum;
}

void diff_data::gauss_zeidel(double e, int maxit, double w) {
    int n = nodes.size();
    vector<double> Ax;
    int k = 0;
    Ax.resize(n);
    do {
        for (int i = 0; i < n; i++)
            Ax[i] = vec[i] - iteration(i, n, m, w, matrix, vec, u, u);
        k++;
        cout << "Итерация: " << k << '\t' << "Относительная невязка: " << scientific << setprecision(15) << norm(n, Ax) / norm(n, vec) << endl;
    } while (k < maxit && norm(n, Ax) / norm(n, vec) >= e);
}
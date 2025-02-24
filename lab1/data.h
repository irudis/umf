#ifndef DATA_DATA
#define DATA_DATA

#include <vector>

using namespace std;

class diff_data {
public:
    void input();
    void points();
    void apply_conditions();
    void build_matrix();
    double func(double x, double y);
    double U(double x, double y);
    double tetta(double x, double y);
    double norm(int n, vector<double>& vec);
    double iteration(int i, int n, int m, double w, vector<vector<double>>& matrix, vector<double>& vec, vector<double>& x, vector<double>& xnext);
    void gauss_zeidel(double e, int maxit, double w);

    vector<pair<double, double>> nodes;
    vector<double> fict;
    vector<double> u;
private:
    vector<vector<double>> blocks;
    vector<double> bc1;
    vector<double> bc2;
    vector<double> vec;
    vector<vector<double>> matrix;

    int N;
    int nx;
    int ny;
    double hx;
    double hy;
    double kx;
    double ky;
    double X_min, X_max, Y_min, Y_max;
    double lambda, gamma;
    int m; 
};

#endif
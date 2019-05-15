#include <iostream>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "CoDecBCH.h"

using namespace std;
using namespace Eigen;

bool iterative_decoder(const vector<vector<int>> &G,
        vector<vector<int>> &Ge,
        CoDecBCH &C) {
    /*
     * для с-регулярного двудольного графа
     */
    int n = G.size();
    vector<vector<int> > raws;
    vector<vector<int> > columns;
    for (int i = 0; i != n; ++i) {
        raws.push_back(vector<int>());
        for (int j = 0; j != n; ++j) {
            if (G[i][j])
                raws.back().push_back(j);
        }
    }
    for (int j = 0; j != n; ++j) {
        columns.push_back(vector<int>());
        for (int i = 0; i != n; ++i) {
            if (G[i][j])
                columns.back().push_back(i);
        }
    }
    int delta = raws[0].size();
    vector<int> recd(delta);
    bool end = true;
    for (int i = 0; i != n; ++i) {
        for (int j = 0; j != delta; ++j) {
            recd[j] = Ge[i][raws[i][j]];
        }
        C.DecoderBCH(recd);
        for (int j = 0; j != delta; ++j) {
            if (end && (Ge[i][raws[i][j]] != recd[j]))
                end = false;
            Ge[i][raws[i][j]] = recd[j];
        }
    }
    if (end)
        return end;
    end = true;
    for (int j = 0; j != n; ++j) {
        for (int i = 0; i != delta; ++i) {
            recd[i] = Ge[columns[j][i]][i];
        }
        C.DecoderBCH(recd);
        for (int i = 0; i != delta; ++i) {
            if (end && (Ge[columns[j][i]][i] != recd[i]))
                end = false;
            Ge[columns[j][i]][i] = recd[i];
        }
    }
    return end;
}

void forward(vector<vector<int>> &Gc,
        vector<vector<int>> &C1,
        vector<vector<int>> &C2,
        int &m, int &n, int &c, int &d) {
    int s = 0, t = 0;
    for (int i = 0; i != m; ++i) {
        for (int j = 0; j != n; ++j) {
            if (Gc[i][j] != -1) {
                C1[s][t] = Gc[i][j];
                if (t == c - 1)
                    ++s;
                else
                    ++t;
            }
        }
    }
    s = 0, t = 0;
    for (int i = 0; i != n; ++i) {
        for (int j = 0; j != m; ++j) {
            if (Gc[j][i] != -1) {
                C2[s][t] = Gc[j][i];
                if (s == d - 1)
                    ++s;
                else
                    ++t;
            }
        }
    }
}


void backward_c1(vector<vector<int>> &Gc,
             vector<vector<int>> &C1,
             int &m, int &n, int &c) {
    int s = 0, t = 0;
    for (int i = 0; i != m; ++i) {
        for (int j = 0; j != n; ++j) {
            if (Gc[i][j] != -1) {
                Gc[i][j] = C1[s][t];
                if (s == c - 1)
                    ++s;
                else
                    ++t;
            }
        }
    }
}

bool backward_c2(vector<vector<int>> &Gc,
                 vector<vector<int>> &C2,
                 int &m, int &n, int &d) {
    bool counter = true;
    int s = 0, t = 0;
    for (int i = 0; i != n; ++i) {
        for (int j = 0; j != m; ++j) {
            if (Gc[j][i] != -1) {
                if (Gc[j][i] == C2[s][t])
                    counter = false;
                Gc[j][i] = C2[s][t];
                if (s == d - 1)
                    ++s;
                else
                    ++t;
            }
        }
    }
    return counter;
}
// #include <typeinfo>
double second_eigen(vector<vector<int>> &G_source) {
    int m = G_source.size();
    int n = G_source[0].size();
    int k = n + m;
    MatrixXd G = MatrixXd::Zero(k, k);
    for (int i = 0; i != m; ++i) {
        for (int j = 0; j != n; ++j) {
            G(i, m + j) = static_cast<double> (G_source[i][j]);
            G(m + j, i) = static_cast<double> (G_source[i][j]);
        }
    }
    SelfAdjointEigenSolver<MatrixXd> es(G);
    // cout << typeid(es.eigenvalues()).name() << endl;
    VectorXd eigenvalues = es.eigenvalues();
    double tmp, mu = 0;
    k = eigenvalues.size();
    tmp = eigenvalues(0);
    mu = tmp;
    for (int i = 0; i != k; ++i) {
        if (eigenvalues(i) > tmp) {
            mu = tmp;
            tmp = eigenvalues(i);
        }

    }
    return mu;
}

int main() {
    /* Полный двудольный граф
     */
    int N, n = 20, t = 5;
    // cin >> n >> t;
    CoDecBCH C(n, t);
    N = n * n;
    vector<vector<int>> G(n, vector<int> (n, 1));
    double  delta = C.get_delta();
    double mu = second_eigen(G);
    int errors_max;
    errors_max = static_cast<int> ((delta / 2) * (delta / 2 - mu / n) * n * n) - 1;
    cout << "errors_were: " << errors_max << "\n";
//    if (C.get_d() >= 3 * mu) {
//        errors_max = static_cast<int> ((delta / 2) * (delta / 2 - mu / n) * n * n) - 1;
//    } else {
//        cout << "d0 < 3 * mu\n";
//        return 0;
//    }
    vector<int> zero(N, 0);
    std::vector<int> indices(N);
    for (int i = 0; i != N; ++i)
        indices[i] = i;
    random_shuffle(indices.begin(), indices.end());
    for (int i = 0; i != errors_max; ++i)
        zero[indices[i]] = 1;
    vector<vector<int>> Ge(n, vector<int> (n, 0));
    int ii = 0;
    for (int i = 0; i != n; ++i) {
        for (int j = 0; j != n; ++j) {
            Ge[i][j] = zero[ii + j];
        }
        ii += n;
    }
    int iter = 0, max_iter = 1000;
    while (!iterative_decoder(G, Ge, C) || iter > max_iter)
        ++max_iter;
    int errors = 0;
    for (int i = 0; i != n; ++i) {
        for (int j = 0; j != n; ++j) {
            if (Ge[i][j])
                ++errors;
        }
    }
    cout << "errors_are: " << errors << "\n";
}

/*int main_() {
    int m, n, c, d;
    cin >> m >> n >> c >> d;
    int             m_1, N_1, n_1, k_1, t_1, d_1;
    int             m_2, N_2, n_2, k_2, t_2, d_2;
    int             p1[21], p2[21];
    int             alpha_to_1[1048576], index_of_1[1048576], g_1[548576];
    int             alpha_to_2[1048576], index_of_2[1048576], g_2[548576];
    int             recd_1[1048576], data_1[1048576], bb_1[548576];
    int             recd_2[1048576], data_2[1048576], bb_2[548576];
    GetPrimePoly(p1, c);
    GetPrimePoly(p2, d);
    KURSACH_CODECBCH_H::GetGF(alpha_to_1, index_of_1, m_1, p1, N_1);
    KURSACH_CODECBCH_H::GetGF(alpha_to_2, index_of_2, m_2, p2, N_2);
    KURSACH_CODECBCH_H::GetGenPoly(alpha_to_1, index_of_1, g_1, k_1, n_1, d_1, t_1, N_1);
    KURSACH_CODECBCH_H::GetGenPoly(alpha_to_2, index_of_2, g_2, k_2, n_2, d_2, t_2, N_2);
    vector<vector<int>> C1(m, vector<int>(c));
    vector<vector<int>> C2(n, vector<int>(d));
    vector<vector<int>> G(m, vector<int>(n));
    vector<vector<int>> Gc(m, vector<int>(n));
    for (int i = 0; i != m; ++i) {
        for (int j = 0; j != n; ++j) {
            cin >> G[i][j];
            if (G[i][j] != 0) {
                Gc[i][j] = (rand() % static_cast<int>(2));
            } else {
                Gc[i][j] = -1;
            }
        }
    }
    double mu = second_eigen(G, m, n);
    bool flag = false;
    int iter = 0;
    forward(Gc, C1, C2, m, n, c, d);
    while (!flag || iter < 1000) {
        for (int i = 0; i != m; ++i) {
            KURSACH_CODECBCH_H::DecoderBCH(C1[i], alpha_to_1, index_of_1, g_1, k_1, n_1, t_1, N_1);
        }
        backward_c1(Gc, C1, m, n, c);
        for (int i = 0; i != n; ++i) {
            KURSACH_CODECBCH_H::DecoderBCH(C2[i], alpha_to_2, index_of_2, g_2, k_2, n_2, t_2, N_2);
        }
        flag = backward_c2(Gc, C2, m, n, d);
    }
}*/
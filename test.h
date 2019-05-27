//
// Created by Руслан Сафронов on 2019-05-27.
//

#include <iostream>
#include <chrono>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "CoDecBCH.h"
#include <cmath>
#include "matplotlibcpp.h"
#include "tanner.h"

class Test{
private:
    int max_iter, k;
    int m, tm, n, tn;
    int c, d;
public:
    Test(int max_iter_, int k_, int m_, int tm_, int n_, int tn_) :
            max_iter(max_iter_),
            k(k_),
            m(m_),
            tm(tm_),
            n(n_),
            tn(tn_)
    {};

    ~Test(){};

    void test() {
        CoDecBCH Cm(m, tm);
        CoDecBCH Cn(n, tn);
        int N = n * m;
        vector<vector<int>> G(m, vector<int>(n, 1));
        vector<vector<int> > raws;
        vector<vector<int> > columns;
        for (int i = 0; i != m; ++i) {
            raws.push_back(vector<int>());
            for (int j = 0; j != n; ++j) {
                if (G[i][j])
                    raws.back().push_back(j);
            }
        }
        for (int j = 0; j != n; ++j) {
            columns.push_back(vector<int>());
            for (int i = 0; i != m; ++i) {
                if (G[i][j])
                    columns.back().push_back(i);
            }
        }
        c = raws[0].size(), d = columns[0].size();
        double boundZ, boundReal = 0.;
        double beta_0 = static_cast<double > (min(tm * m, tn * n)) / N;
        double beta_1 = static_cast<double > (max(tm * m, tn * n)) / N;
        double  deltm = Cm.get_delta(), deltn = Cn.get_delta(), mu = second_eigen(G);
        if (c == d) {
            boundZ = 0.99 * (deltm / 2) * (deltm / 2 - mu / c);
        } else {
            double gamma = mu / sqrt(static_cast<double >(c * d));
            double alpha = 2 * gamma / sqrt(deltm * deltn);
            if (alpha < 1)
                alpha = 1.0 - 1e-6;
            boundZ = (alpha * deltm * deltn
                      - 2 * gamma * sqrt(deltm * deltn))
                     / (4 * (1 - gamma));
        }
        int BoundZ = static_cast<int >(boundZ * N), BoundReal = 0;
        int Beta_0 = static_cast<int >(beta_0 * N);
        int Beta_1 = static_cast<int >(beta_1 * N);
        VectorXd beta = VectorXd::LinSpaced(k, boundZ, 1.0);
        std::vector<double> x, y, iters, time, Z, R, v(k, Beta_0), z(k, Beta_1);
        int iterG = 0, elapsed_secondsG = 0;
        for (int t = 0; t != k; ++t) {
            int errors_max = static_cast<int> (beta(t) * N);
            vector<int> zero(N, 0);
            std::vector<int> indices(N);
            for (int i = 0; i != N; ++i)
                indices[i] = i;
            random_shuffle(indices.begin(), indices.end());
            for (int i = 0; i != errors_max; ++i)
                zero[indices[i]] = 1;
            vector<vector<int>> Ge(m, vector<int>(n, 0));
            int ii = 0;
            for (int i = 0; i != m; ++i) {
                for (int j = 0; j != n; ++j) {
                    Ge[i][j] = zero[ii + j];
                }
                ii += n;
            }
            int iter = 0, errors = 0;
            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();
            while (!iterative_decoder_mn(G, raws, columns, Ge, c, d, Cm, Cn)
            && iter < max_iter)
                ++iter;
            iterG += iter;
            end = std::chrono::system_clock::now();
            int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>
                    (end-start).count();
            elapsed_secondsG += elapsed_seconds;
            for (int i = 0; i != m; ++i) {
                for (int j = 0; j != n; ++j) {
                    if (Ge[i][j])
                        ++errors;
                }
            }
            if (!errors || !BoundReal) {
                BoundReal = errors_max;
                boundReal = beta(t);
            }

            iters.push_back(iter);
            time.push_back(elapsed_seconds);
            x.push_back(errors_max);
            y.push_back(errors);
        }
        for (int t = 0; t != k; ++t) {
            Z.push_back(BoundZ);
            R.push_back(BoundReal);
        }
        char title [50], path[50], boundZplt[50], boundRplt[50], timePlt[50], itersPlt[50];
        sprintf(title,
                "(%d, %d; %d, %d) - BCH Code | (%d, %d) - Degrees",
                m, tm, n, tn, c, d);
        sprintf(path,
                "../graphics/%d_%d_%d_%d_BCH_%d_%d.png",
                m, tm, n, tn, c, d);
        sprintf(boundZplt,
                "Estimated bound of errors = %d (%.2f %%)",
                BoundZ,
                boundZ * 100);
        sprintf(boundRplt,
                "Real bound of errors = %d (%.2f %%)",
                BoundReal,
                boundReal * 100);
        sprintf(timePlt,
                "Mean time = %.2f s",
                static_cast<double >(elapsed_secondsG) / k);
        sprintf(itersPlt,
                "Mean iters = %.2f",
                static_cast<double >(iterG) / k);

        plt::figure_size(780, 1200);
        plt::suptitle(title);

        plt::subplot(3, 1, 1);
        plt::ylabel("Errors remaining");
        plt::plot(x, y);
        plt::named_plot(boundZplt, Z, y, "r--");
        plt::named_plot(boundRplt, R, y, "c--");
        plt::legend();
        plt::grid(true);

        plt::subplot(3, 1, 2);
        plt::ylabel("Time(s)");
        plt::named_plot(timePlt, x, time, "y--");
        plt::legend();
        plt::grid(true);

        plt::subplot(3, 1, 3);
        plt::ylabel("Iterations");
        plt::named_plot(itersPlt, x, iters, "go--");
        plt::legend();
        plt::grid(true);

        plt::save(path);
    }
};

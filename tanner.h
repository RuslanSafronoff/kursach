//
// Created by Руслан Сафронов on 2019-05-27.
//

#pragma once
#include <iostream>
#include <chrono>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include "CoDecBCH.h"
#include <cmath>
#include "matplotlibcpp.h"


using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;

bool iterative_decoder(const vector<vector<int>> &G,
                       vector<vector<int>> &raws,
                       vector<vector<int>> &columns,
                       vector<vector<int>> &Ge,
                       const int &c,
                       const int &d,
                       const CoDecBCH &Cm,
                       const CoDecBCH &Cn) {
    int m = G.size();
    int n = G[0].size();
    vector<int> recdc(c);
    vector<int> recdd(d);
    bool end = true;
    for (int i = 0; i != m; ++i) {
        for (int j = 0; j != c; ++j) {
            recdc[j] = Ge[i][raws[i][j]];
        }
        Cm.DecoderBCH(recdc);
        for (int j = 0; j != c; ++j) {
            if (end && (Ge[i][raws[i][j]] != recdc[j]))
                // end = false;
                Ge[i][raws[i][j]] = recdc[j];
        }
    }
    end = true;
    for (int j = 0; j != n; ++j) {
        for (int i = 0; i != d; ++i) {
            recdd[i] = Ge[columns[j][i]][j];
        }
        Cn.DecoderBCH(recdd);
        for (int i = 0; i != d; ++i) {
            if (end && (Ge[columns[j][i]][j] != recdd[i]))
                end = false;
            Ge[columns[j][i]][j] = recdd[i];
        }
    }
    return end;
}

double second_eigen(const vector<vector<int>> &G_source) {
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
    VectorXd eigenvalues = es.eigenvalues();
    for (int i = 1; i != k + 1; ++i) {
        if (eigenvalues(k - i) != eigenvalues(k - 1)) {
            return eigenvalues(k - i);
        }
    }
}
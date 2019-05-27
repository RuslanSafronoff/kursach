//
// Created by Руслан Сафронов on 2019-04-18.
//
#pragma once
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

class CoDecBCH {
private:
    int             m, N, n, k, t, d;
    int             p[21];
    vector<int>     alpha_to, index_of, g; //(1048576), index_of(1048576), g(548576);
public:
    CoDecBCH(int n_, int t_) : n(n_), t(t_) {
        g.resize(548576);
        GetPrimePoly();
        GetGF();
        GetGenPoly();
    }

    double get_delta() {
        return static_cast<double> (d) / static_cast<double> (n);
    }

    double get_d() {
        return static_cast<double> (d);
    }

    void
    GetPrimePoly()
    /*
     *  Read m, the degree of a primitive irreducible polynomial p(x).
     *  Get precomputed coefficients p[] of p(x).
     *  Read the code length.
     */
    {
        int nMin;
        m = 0;
        N = 1;
        do  {
            nMin = N - 1;
            ++m;
            N *= 2;
        } while ( !((n < N) && (n > nMin)) );
        alpha_to.resize(N);
        index_of.resize(N);
        --N;
        for (int i = 1; i < m; i++)
            p[i] = 0;
        p[0] = p[m] = 1;
        if (m == 2) p[1] = 1;
        else if (m == 3) p[1] = 1;
        else if (m == 4) p[1] = 1;
        else if (m == 5) p[2] = 1;
        else if (m == 6) p[1] = 1;
        else if (m == 7) p[1] = 1;
        else if (m == 8) p[4] = p[5] = p[6] = 1;
        else if (m == 9) p[4] = 1;
        else if (m == 10) p[3] = 1;
        else if (m == 11) p[2] = 1;
        else if (m == 12) p[3] = p[4] = p[7] = 1;
        else if (m == 13) p[1] = p[3] = p[4] = 1;
        else if (m == 14) p[1] = p[11] = p[12] = 1;
        else if (m == 15) p[1] = 1;
        else if (m == 16) p[2] = p[3] = p[5] = 1;
        else if (m == 17) p[3] = 1;
        else if (m == 18) p[7] = 1;
        else if (m == 19) p[1] = p[5] = p[6] = 1;
        else if (m == 20) p[3] = 1;
    }

    void
    GetGF()
    /*
     * Получение поля GF(2**m) из примитивного неприводимого многочлена p(x)
     * с коэффициентами p[0]..p[m].
     *
     * Получение таблицы:
     *   индексная форма -> полиномиальная форма: alpha_to[] содержит j=alpha^i;
     *   полиномиальная форма -> индексная форма:	index_of[j=alpha^i] = i
     *
     * alpha = 2 примитивный элемент GF(2**m)
     */
    {
        int i, mask;

        mask = 1;
        alpha_to[m] = 0;
        for (i = 0; i < m; i++) {
            alpha_to[i] = mask;
            index_of[alpha_to[i]] = i;
            if (p[i] != 0)
                alpha_to[m] ^= mask;
            mask <<= 1;
        }
        index_of[alpha_to[m]] = m;
        mask >>= 1;
        for (i = m + 1; i < N; i++) {
            if (alpha_to[i - 1] >= mask)
                alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
            else
                alpha_to[i] = alpha_to[i - 1] << 1;
            index_of[alpha_to[i]] = i;
        }
        index_of[0] = -1;
    }


    void
    GetGenPoly()
    /*
     * Вычисление порождающего многочлена бинарного БЧХ-кода. Сначала генерируются
     * циклические множества modulo 2**m - 1, cycle[][] =  (i, 2*i, 4*i, ..., 2^l*i). Затем
     * определяются множества, которые содержат целые числа из {1..(d-1)}.
     * Порождающий многочлен вычисляется
     * как продукт линейных факторов формы (x+alpha^i), для каждого i в
     * циклических множествах.
     */
    {
        int ii, jj, ll, kaux;
        int test, aux, nocycles, root, noterms, rdncy;
        int cycle[1024][21], size[1024], min[1024], zeros[1024];

        /* Генерируются циклические множества modulo N, N = 2**m - 1*/
        cycle[0][0] = 0;
        size[0] = 1;
        cycle[1][0] = 1;
        size[1] = 1;
        jj = 1;            /* индекс циклического множества */
        do {
            /* Генерируется jj-ое килическое множество */
            ii = 0;
            do {
                ii++;
                cycle[jj][ii] = (cycle[jj][ii - 1] * 2) % N;
                size[jj]++;
                aux = (cycle[jj][ii] * 2) % N;
            } while (aux != cycle[jj][0]);
            /* Следующий представитель циклического множества */
            ll = 0;
            do {
                ll++;
                test = 0;
                for (ii = 1; ((ii <= jj) && (!test)); ii++)
                    /* Проверяется предыдущее циклическое множество */
                    for (kaux = 0; ((kaux < size[ii]) && (!test)); kaux++)
                        if (ll == cycle[ii][kaux])
                            test = 1;
            } while ((test) && (ll < (N - 1)));
            if (!(test)) {
                jj++;    /* следующий индекс циклического множества */
                cycle[jj][0] = ll;
                size[jj] = 1;
            }
        } while (ll < (N - 1));
        nocycles = jj;        /* число циклических множеств modulo N */

        d = 2 * t + 1;

        /* Поиск корней 1, 2, ..., d-1 в циклических множествах */
        kaux = 0;
        rdncy = 0;
        for (ii = 1; ii <= nocycles; ii++) {
            min[kaux] = 0;
            test = 0;
            for (jj = 0; ((jj < size[ii]) && (!test)); jj++)
                for (root = 1; ((root < d) && (!test)); root++)
                    if (root == cycle[ii][jj]) {
                        test = 1;
                        min[kaux] = ii;
                    }
            if (min[kaux]) {
                rdncy += size[min[kaux]];
                kaux++;
            }
        }
        noterms = kaux;
        kaux = 1;
        for (ii = 0; ii < noterms; ii++)
            for (jj = 0; jj < size[min[ii]]; jj++) {
                zeros[kaux] = cycle[min[ii]][jj];
                kaux++;
            }

        k = n - rdncy;

        if (k < 0) {
            printf("==Недопустимые параметры!==\n");
            exit(0);
        }

        /* Вычисление порождающего многочлена */
        g[0] = alpha_to[zeros[1]];
        g[1] = 1;        /* g(x) = (X + zeros[1]) изначально */
        for (ii = 2; ii <= rdncy; ii++) {
            g[ii] = 1;
            for (jj = ii - 1; jj > 0; jj--)
                if (g[jj] != 0)
                    g[jj] = g[jj - 1] ^ alpha_to[(index_of[g[jj]] + zeros[ii]) % N];
                else
                    g[jj] = g[jj - 1];
            g[0] = alpha_to[(index_of[g[0]] + zeros[ii]) % N];
        }
    }

    void
    DecoderBCH(vector<int> &recd) const
    /*
     * Мы получили биты в recd[i], i=0..(N-1).
     *
     * Вычисление 2*t синдромов подстановкой alpha^i в rec(X) и
     * оценкой, синдромы хранятся в s[i], i=1..2t (s[0] = 0).
     * Затем мы используем алгоритм Берлекэмпа, чтобы найти
     * многочлен  локаторов ошибки elp[i].
     *
     * Если степень elp > t, мы не сможем исправить все ошибки, и
     * и мы обнаружили неисправляемый шаблон ошибок. Мы выводим
     * информационные биты неисправленными
     *
     * Если степень elp is <= t, мы последовательно подставляем alpha^i , i=1..N в elp,
     * чтобы найти корни, затем инвертированные корни, номера локации ошибок.
     * Эта так называемая процедура Ченя.
     *
     * Если число ошибок не равно степени elp,
     * декодер предполагает, что ошибок больше чем t и не может исправить
     * их, а только обнаружить. Мы выводим
     * информационные биты неисправленными.
     */
    {
        int i, j, u, q, t2, count = 0, syn_error = 0;
        int elp[1026][1024], d[1026], l[1026], u_lu[1026], s[1025];
        int root[200], loc[200], err[1024], reg[201];

        t2 = 2 * t;

        /* Сначала формируются синдромы */
        // printf("S(x) = ");
        for (i = 1; i <= t2; i++) {
            s[i] = 0;
            for (j = 0; j < n; j++)
                if (recd[j] != 0)
                    s[i] ^= alpha_to[(i * j) % N];
            if (s[i] != 0)
                syn_error = 1; /* установить флаг ошибки, если синдром ненулевой */
            /*
             * Заметьте: Выйдите из программы здесь,
             *           если используете код только для индикации ошибок
             */
            /* конвертирование синдрома из полиномиальной формы в индексную  */
            s[i] = index_of[s[i]];
        }

        if (syn_error) {    /* если есть ошибки, пытаемся исправить их */
            /*
             * Вычисление многочлена локаторов ошибки с помощью
             * итеративного алгоритма Берлекэмпа-Мэсси.
             * d[u] – mu-ое расхождение, где
             * u='mu'+1 and 'mu' (греческая буква) – номер шага
             * от -1 до 2*t, l[u] – степень
             * elp на том шаге, и u_l[u] – различие между
             * номером шага и степенью elp.
             */
            /* инициализация табличных записей */
            d[0] = 0;            /* индексная форма */
            d[1] = s[1];        /* индексная форма */
            elp[0][0] = 0;        /* индексная форма */
            elp[1][0] = 1;        /* полиномиальная форма */
            for (i = 1; i < t2; i++) {
                elp[0][i] = -1;    /* индексная форма */
                elp[1][i] = 0;    /* полиномиальная форма */
            }
            l[0] = 0;
            l[1] = 0;
            u_lu[0] = -1;
            u_lu[1] = 0;
            u = 0;

            do {
                u++;
                if (d[u] == -1) {
                    l[u + 1] = l[u];
                    for (i = 0; i <= l[u]; i++) {
                        elp[u + 1][i] = elp[u][i];
                        elp[u][i] = index_of[elp[u][i]];
                    }
                } else
                    /*
                     * поиск слов с наибольшим u_lu[q], для которых d[q]!=0
                     */
                {
                    q = u - 1;
                    while ((d[q] == -1) && (q > 0))
                        q--;
                    /* нашли первый ненулевой d[q]  */
                    if (q > 0) {
                        j = q;
                        do {
                            j--;
                            if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
                                q = j;
                        } while (j > 0);
                    }

                    /*
                     * нашли такое q, что d[u]!=0 и u_lu[q] максимальное
                     */
                    /* сохранение степени ноаого многочлена elp */
                    if (l[u] > l[q] + u - q)
                        l[u + 1] = l[u];
                    else
                        l[u + 1] = l[q] + u - q;

                    /* формируем новый elp(x) */
                    for (i = 0; i < t2; i++)
                        elp[u + 1][i] = 0;
                    for (i = 0; i <= l[q]; i++)
                        if (elp[q][i] != -1)
                            elp[u + 1][i + u - q] =
                                    alpha_to[(d[u] + N - d[q] + elp[q][i]) % N];
                    for (i = 0; i <= l[u]; i++) {
                        elp[u + 1][i] ^= elp[u][i];
                        elp[u][i] = index_of[elp[u][i]];
                    }
                }
                u_lu[u + 1] = u - l[u + 1];

                /* формируем (u+1)-ое различие */
                if (u < t2) {
                    /* не вычислено различий на последней итерации */
                    if (s[u + 1] != -1)
                        d[u + 1] = alpha_to[s[u + 1]];
                    else
                        d[u + 1] = 0;
                    for (i = 1; i <= l[u + 1]; i++)
                        if ((s[u + 1 - i] != -1) && (elp[u + 1][i] != 0))
                            d[u + 1] ^= alpha_to[(s[u + 1 - i]
                                                  + index_of[elp[u + 1][i]]) % N];
                    /* кладём d[u+1] в индексную форму */
                    d[u + 1] = index_of[d[u + 1]];
                }
            } while ((u < t2) && (l[u + 1] <= t));

            u++;
            if (l[u] <= t) {/* Можем исправить ошибки */
                /* кладём elp в индексную форму */
                for (i = 0; i <= l[u]; i++)
                    elp[u][i] = index_of[elp[u][i]];

                /* Процедура Ченя: находим корни многочлена локаторов ошибки */
                for (i = 1; i <= l[u]; i++)
                    reg[i] = elp[u][i];
                count = 0;
                for (i = 1; i <= N; i++) {
                    q = 1;
                    for (j = 1; j <= l[u]; j++)
                        if (reg[j] != -1) {
                            reg[j] = (reg[j] + j) % N;
                            q ^= alpha_to[reg[j]];
                        }
                    if (!q) {    /* сохраняем корень и индексы
						 * номеров локации ошибки */
                        root[count] = i;
                        loc[count] = N - i;
                        count++;
                    }
                }
                if (count == l[u])
                    /* число корней = степени elp, значит <= t ошибок */
                    for (i = 0; i < l[u]; i++)
                        recd[loc[i]] ^= 1;
                else {
                    // printf("==Незавершённое декодирование==\n");
                }    /* elp имеет степень >t, значит неразрешимо */
            } else;
                // printf("==Незавершённое декодирование==\n");
        }
    }
};

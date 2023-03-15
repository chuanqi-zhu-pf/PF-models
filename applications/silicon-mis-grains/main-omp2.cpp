
// Phase-field model of multirystalline silicon growth
// Implemented by openmp, which allow the moving frame follwing the motion of the solid-liquid growth front

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define NTH 1
#define NDX 200 // 差分計算における計算領域一辺の分割数
#define NDY 150 // 差分計算における計算領域一辺の分割数
#define NDZ 1
#define FN 4
#define N 4

int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int nm = N - 1;
double PI = 3.1415926;
double RR = 8.3145;

double aij[N][N], wij[N][N], mij[N][N], fij[N][N], anij[N][N];
double th00ij[N][N], thij[N][N][3];
double ax00ij[N][N][3], axij[N][N][3][3];
double a100, a200, a300, a1, a2, a3;
double th00, th, vp, eta;
double thii, vpii, etaii;

double vpij[N][N], etaij[N][N];

double face[FN][3], ang[N][3];
double face111[4][3], face110[6][3], face100[3][3];

int ni, nj, ix, iy, iz;

double ****phi, ****phi2;
int ***phiNum, ****phiIdx;
double ***temp, ***mob, ***intphi;

double dtime, L, dx;
double M0, W0, A0;
double W1, A1;
double zeta1, ww1, coef, rp0, rp1, del, al0;

double Tm, Ti, Tg, Tr;
double gt0, gt1; // 粒界エネルギ密度
double delta;    // 粒界幅（差分ブロック数にて表現）
double mobi;     // 粒界の易動度

int allL, allS;
int fronti, endi, intpos;

void datasave2d(int step);

int main()
{
    dx = 2.0e-5;
    dtime = 1.0e-2;
    delta = 7.0;

    Tm = 1687.0;
    Ti = 1686.0;
    Tg = 900.0;
    Tr = 0.009;

    gt0 = 10.0;
    gt1 = 14.0;
    mobi = 1.0;

    A0 = 8.0 * delta * gt0 / PI / PI;
    W0 = 4.0 * gt0 / delta;
    M0 = mobi * PI * PI / (8.0 * delta);

    A1 = 8.0 * delta * gt1 / PI / PI;
    W1 = 4.0 * gt1 / delta;

    del = 0.36;
    coef = 1.0 / 1.62454431;
    al0 = 54.7 / 180.0 * PI;
    zeta1 = 0.8;
    ww1 = 90.0 / 10.0;
    rp0 = 0.05;
    rp1 = 0.05;

    for (ni = 0; ni <= nm; ni++)
    {
        for (nj = 0; nj <= nm; nj++)
        {
            wij[ni][nj] = W0;
            aij[ni][nj] = A0;
            mij[ni][nj] = M0;
            anij[ni][nj] = 0.0;
            if ((ni == 0) || (nj == 0))
            {
                anij[ni][nj] = 1.0;
            }
            if (ni > 0 && nj > 0)
            {
                wij[ni][nj] = W1;
                aij[ni][nj] = A1;
                mij[ni][nj] = mij[ni][nj] * 0.00;
            }
            if (ni == nj)
            {
                wij[ni][nj] = 0.0;
                aij[ni][nj] = 0.0;
                mij[ni][nj] = 0.0;
                anij[ni][nj] = 0.0;
            }
        }
    }

    // grain 1
    // [grain][axis]
    ang[1][0] = 45.0 / 180.0 * PI;
    ang[1][1] = 54.7 / 180.0 * PI;

    // [phase][phase][axis]
    thij[1][0][0] = thij[0][1][0] = ang[1][0];
    thij[1][0][1] = thij[0][1][1] = ang[1][1];

    // [phase][phase][axis][component]
    axij[1][0][0][0] = axij[0][1][0][0] = 1.0;
    axij[1][0][0][1] = axij[0][1][0][1] = 0.0;
    axij[1][0][0][2] = axij[0][1][0][2] = 0.0;

    axij[1][0][1][0] = axij[0][1][1][0] = 0.0;
    axij[1][0][1][1] = axij[0][1][1][1] = 0.0;
    axij[1][0][1][2] = axij[0][1][1][2] = 1.0;

    // grain 2
    // [grain][axis]
    ang[2][0] = 45.0 / 180.0 * PI;

    // [phase][phase][axis]
    thij[2][0][0] = thij[0][2][0] = ang[2][0];

    // [phase][phase][axis][component]
    axij[2][0][0][0] = axij[0][2][0][0] = 1.0;
    axij[2][0][0][1] = axij[0][2][0][1] = 0.0;
    axij[2][0][0][2] = axij[0][2][0][2] = 0.0;

    phi = new double ***[N];
    phi2 = new double ***[N];
    for (ni = 0; ni <= nm; ni++)
    {
        phi[ni] = new double **[NDX];
        phi2[ni] = new double **[NDX];
        for (ix = 0; ix <= ndmx; ix++)
        {
            phi[ni][ix] = new double *[NDY];
            phi2[ni][ix] = new double *[NDY];
            for (iy = 0; iy <= ndmy; iy++)
            {
                phi[ni][ix][iy] = new double[NDZ];
                phi2[ni][ix][iy] = new double[NDZ];
            }
        }
    }

    phiIdx = new int ***[N + 1];
    for (ni = 0; ni <= N; ni++)
    {
        phiIdx[ni] = new int **[NDX];
        for (ix = 0; ix <= ndmx; ix++)
        {
            phiIdx[ni][ix] = new int *[NDY];
            for (iy = 0; iy <= ndmy; iy++)
            {
                phiIdx[ni][ix][iy] = new int[NDZ];
            }
        }
    }

    intphi = new double **[NDX];
    temp = new double **[NDX];
    mob = new double **[NDX];
    for (ix = 0; ix <= ndmx; ix++)
    {
        intphi[ix] = new double *[NDY];
        temp[ix] = new double *[NDY];
        mob[ix] = new double *[NDY];
        for (iy = 0; iy <= ndmy; iy++)
        {
            intphi[ix][iy] = new double[NDZ];
            temp[ix][iy] = new double[NDZ];
            mob[ix][iy] = new double[NDZ];
        }
    }

    phiNum = new int **[NDX];
    for (ix = 0; ix <= ndmx; ix++)
    {
        phiNum[ix] = new int *[NDY];

        for (iy = 0; iy <= ndmy; iy++)
        {
            phiNum[ix][iy] = new int[NDZ];
        }
    }

    face[0][0] = 1.0;
    face[0][1] = 1.0;
    face[0][2] = 1.0;

    face[1][0] = -1.0;
    face[1][1] = 1.0;
    face[1][2] = 1.0;

    face[2][0] = 1.0;
    face[2][1] = -1.0;
    face[2][2] = 1.0;

    face[3][0] = 1.0;
    face[3][1] = 1.0;
    face[3][2] = -1.0;

    // 111
    face111[0][0] = 1.0;
    face111[0][1] = 1.0;
    face111[0][2] = 1.0;

    face111[1][0] = -1.0;
    face111[1][1] = 1.0;
    face111[1][2] = 1.0;

    face111[2][0] = 1.0;
    face111[2][1] = -1.0;
    face111[2][2] = 1.0;

    face111[3][0] = 1.0;
    face111[3][1] = 1.0;
    face111[3][2] = -1.0;

    // 110
    face110[0][0] = 1.0;
    face110[0][1] = 1.0;
    face110[0][2] = 0.0;

    face110[1][0] = -1.0;
    face110[1][1] = 1.0;
    face110[1][2] = 0.0;

    face110[2][0] = 0.0;
    face110[2][1] = 1.0;
    face110[2][2] = 1.0;

    face110[3][0] = 0.0;
    face110[3][1] = -1.0;
    face110[3][2] = 1.0;

    face110[4][0] = 1.0;
    face110[4][1] = 0.0;
    face110[4][2] = 1.0;

    face110[5][0] = -1.0;
    face110[5][1] = 0.0;
    face110[5][2] = 1.0;

    // 100
    face100[0][0] = 1.0;
    face100[0][1] = 0.0;
    face100[0][2] = 0.0;

    face100[1][0] = 0.0;
    face100[1][1] = 1.0;
    face100[1][2] = 0.0;

    face100[2][0] = 0.0;
    face100[2][1] = 0.0;
    face100[2][2] = 1.0;

    printf("Diffuson stability: %e\n", A0 * M0 * dtime);
    printf("Reaction stability: %e\n", (Tm - Ti) * delta / gt0 / PI);

    for (ix = 0; ix <= ndmx; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                for (ni = 1; ni <= nm; ni++)
                {
                    phi[ni][ix][iy][iz] = 0.0;
                }
                phi[0][ix][iy][iz] = 1.0; // nm番目のフェーズフィールドを１に初期化
            }
        }
    }

    for (ix = 0; ix <= ndmx; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                if (ix < NDX / 8 && iy < NDY / 2)
                // if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) < NDX / 6 * NDX / 6)
                {
                    phi[1][ix][iy][iz] = 1.0;
                    phi[2][ix][iy][iz] = 0.0;
                    phi[3][ix][iy][iz] = 0.0;
                    phi[0][ix][iy][iz] = 0.0;
                }
                else if (ix < NDX / 8 && iy >= NDY / 2) // && j < NDY * 2 / 3)
                {
                    phi[1][ix][iy][iz] = 0.0;
                    phi[2][ix][iy][iz] = 1.0;
                    phi[3][ix][iy][iz] = 0.0;
                    phi[0][ix][iy][iz] = 0.0;
                }
                // else if (i < NDX / 8 && j >= NDY * 2 / 3)
                // {
                //     phi[1][i][j][k] = 0.0;
                //     phi[2][i][j][k] = 0.0;
                //     phi[3][i][j][k] = 1.0;
                //     phi[0][i][j][k] = 0.0;
                // }
                else
                {
                    phi[1][ix][iy][iz] = 0.0;
                    phi[2][ix][iy][iz] = 0.0;
                    phi[3][ix][iy][iz] = 0.0;
                    phi[0][ix][iy][iz] = 1.0;
                }
                mob[ix][iy][iz] = 0.0;
                temp[ix][iy][iz] = Ti - NDX / 8 * dx * Tg + ix * dx * Tg;
            }
        }
    }

#pragma omp parallel num_threads(NTH)
    {
        //=========================================================
        double F0;                                      // 粒界移動の駆動力
        int i, j, k, l, ii, jj, kk, ll, it, rn, phinum; // 整数
        int ip, im, jp, jm, kp, km;                     // 整数
        int n1, n2, n3;                                 // 整数
        int istep, nstep, p, pstep;                     // 計算カウント数の最大値（計算終了カウント）
        double sum1, sum2, sum3;                        // 各種の和の作業変数
        double pddtt;                                   // フェーズフィールドの時間変化率

        double miijj, fiijj;
        double min_val;
        int min_idx;

        double phidx, phidy, phidz, phiabs, phiabs2;
        double phidxx, phidxy, phidxz;
        double phidyx, phidyy, phidyz;
        double phidzx, phidzy, phidzz;

        double xxp, xyp, xzp;
        double yxp, yyp, yzp;
        double zxp, zyp, zzp;

        double xxpii, xypii, xzpii;
        double yxpii, yypii, yzpii;
        double zxpii, zypii, zzpii;

        double dphiabs2dx, dphiabs2dy, dphiabs2dz;
        double dphiabs2dphix, dphiabs2dphiy, dphiabs2dphiz;
        double dphiabs2dphixdx, dphiabs2dphiydy, dphiabs2dphizdz;

        double al, alm, am;

        double nx, ny, nz;
        double nxx, nxy, nxz;
        double nyx, nyy, nyz;
        double nzx, nzy, nzz;
        double ux00, uy00, uz00;
        double ux0, uy0, uz0;
        double ux, uy, uz, uu;

        double phidxii, phidyii, phidzii, phiabs2ii;
        double nxii, nyii, nzii, alii;

        double nxphix, nyphix, nzphix;
        double nxphiy, nyphiy, nzphiy;
        double nxphiz, nyphiz, nzphiz;

        double nxphixdx, nyphixdx, nzphixdx;
        double nxphiydy, nyphiydy, nzphiydy;
        double nxphizdz, nyphizdz, nzphizdz;

        double PP, QQ, CC, SS;
        double dPdx, dPdy, dPdz;
        double dQdx, dQdy, dQdz;
        double dPdphix, dPdphiy, dPdphiz;
        double dQdphix, dQdphiy, dQdphiz;
        double dPdphixdx, dPdphiydy, dPdphizdz;
        double dQdphixdx, dQdphiydy, dQdphizdz;

        double epsilon0, ep;
        double epdx, epdy, epdz;
        double epdphix, epdphiy, epdphiz;
        double epdphixdx, epdphiydy, epdphizdz;

        double termx, termy, termz;
        double termiikk, termjjkk;

        int start, end, rows, offset, th_id;

        istep = 0;
        nstep = 15001;
        pstep = 500;
        th_id = omp_get_thread_num();
        rows = NDX / NTH;
        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;

    start:;

        //--------------------------------------------
        if (th_id == 0)
        {
            for (i = 0; i <= ndmx; i++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    for (k = 0; k <= ndmz; k++)
                    {
                        temp[i][j][k] -= Tr * dtime;
                    }
                }
            }

            allL = 1;
            // search interface front
            for (i = ndmx; i <= ndmx; i--)
            {
                if (allL == 0)
                {
                    fronti = i;
                    intpos = i;
                    break;
                }
                for (j = 0; j <= ndmy; j++)
                {
                    for (k = 0; k <= ndmz; k++)
                    {
                        if (phi[0][i][j][k] == 0.0)
                        {
                            allL = 0;
                            break;
                        }
                    }
                    if (allL == 0)
                    {
                        break;
                    }
                }
            }

            fronti = fronti + int(delta * 1.5);

            allS = 1;
            // search interface end
            for (i = 0; i <= ndmx; i++)
            {
                if (allS == 0)
                {
                    endi = i;
                    break;
                }
                for (j = 0; j <= ndmy; j++)
                {
                    for (k = 0; k <= ndmz; k++)
                    {
                        if (phi[0][i][j][k] == 1.0)
                        {
                            allS = 0;
                            break;
                        }
                    }
                    if (allS == 0)
                    {
                        break;
                    }
                }
            }

            endi = endi - int(delta * 1.5);

            if (istep % pstep == 0)
            {
                printf("--------------------\n");
                printf("- current step is: %f (s)\n", istep * dtime);
                printf("- interface is at: %d (grid)\n", intpos);
                printf("- undercooling of interface is: %f (K)\n", (Tm - temp[intpos][0][0]));

                datasave2d(istep);
            }
        }

        //--------------------------------------------
#pragma omp barrier

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    kp = k + 1;
                    km = k - 1;
                    if (i == ndmx)
                    {
                        ip = ndmx;
                    }
                    if (i == 0)
                    {
                        im = 0;
                    }
                    if (j == ndmy)
                    {
                        jp = ndmy;
                    }
                    if (j == 0)
                    {
                        jm = 0;
                    }
                    if (k == ndmz)
                    {
                        kp = 0;
                    }
                    if (k == 0)
                    {
                        km = ndmz;
                    }

                    //--- 位置(i,j)およびその周囲(i±1,j±1)において、pが０ではない方位の個数---
                    phinum = 0;
                    for (ii = 0; ii <= nm; ii++)
                    {
                        if ((phi[ii][i][j][k] > 0.0) ||
                            ((phi[ii][i][j][k] == 0.0) && (phi[ii][ip][j][k] > 0.0) ||
                             (phi[ii][im][j][k] > 0.0) ||
                             (phi[ii][i][jp][k] > 0.0) ||
                             (phi[ii][i][jm][k] > 0.0) ||
                             (phi[ii][i][j][kp] > 0.0) ||
                             (phi[ii][i][j][km] > 0.0)))
                        {
                            phinum++;
                            phiIdx[phinum][i][j][k] = ii;
                            // printf("%d  ", n00);
                        }
                    }
                    phiNum[i][j][k] = phinum;
                    mob[i][j][k] = 0.0;
                }
            }
        }

#pragma omp barrier

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    ip = i + 1;
                    im = i - 1;
                    jp = j + 1;
                    jm = j - 1;
                    kp = k + 1;
                    km = k - 1;
                    if (i == ndmx)
                    {
                        ip = ndmx;
                    }
                    if (i == 0)
                    {
                        im = 0;
                    }
                    if (j == ndmy)
                    {
                        jp = ndmy;
                    }
                    if (j == 0)
                    {
                        jm = 0;
                    }
                    if (k == ndmz)
                    {
                        kp = 0;
                    }
                    if (k == 0)
                    {
                        km = ndmz;
                    }

                    for (n1 = 1; n1 <= phiNum[i][j][k]; n1++)
                    {
                        ii = phiIdx[n1][i][j][k];
                        pddtt = 0.0;

                        for (n2 = 1; n2 <= phiNum[i][j][k]; n2++)
                        {
                            jj = phiIdx[n2][i][j][k];
                            sum1 = 0.0;
                            for (n3 = 1; n3 <= phiNum[i][j][k]; n3++)
                            {
                                kk = phiIdx[n3][i][j][k];

                                // calculate the interface normal and deirivatives of the phase field
                                phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]);
                                phidxy = (phi[kk][ip][jp][k] + phi[kk][im][jm][k] - phi[kk][im][jp][k] - phi[kk][ip][jm][k]) / 4.0;
                                phidxz = (phi[kk][ip][j][kp] + phi[kk][im][j][km] - phi[kk][im][j][kp] - phi[kk][ip][j][km]) / 4.0;

                                phidyx = phidxy;
                                phidyz = (phi[kk][i][jp][kp] + phi[kk][i][jm][km] - phi[kk][i][jm][kp] - phi[kk][i][jp][km]) / 4.0;
                                phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]);

                                phidzx = phidxz;
                                phidzy = phidyz;
                                phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]);

                                phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0;
                                phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0;
                                phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0;

                                phiabs2 = phidx * phidx + phidy * phidy + phidz * phidz;
                                phiabs = sqrt(phiabs2);

                                dphiabs2dx = 2.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz);
                                dphiabs2dy = 2.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz);
                                dphiabs2dz = 2.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz);

                                dphiabs2dphix = 2.0 * phidx;
                                dphiabs2dphiy = 2.0 * phidy;
                                dphiabs2dphiz = 2.0 * phidz;

                                dphiabs2dphixdx = 2.0 * phidxx;
                                dphiabs2dphiydy = 2.0 * phidyy;
                                dphiabs2dphizdz = 2.0 * phidzz;

                                if (anij[ii][kk] == 1.0 && phiabs2 != 0.0)
                                {
                                    nx = phidx / phiabs;
                                    ny = phidy / phiabs;
                                    nz = phidz / phiabs;

                                    min_val = 0.0;
                                    min_idx = 0;
                                    for (l = 0; l < FN; l++)
                                    {
                                        ux0 = face[l][0];
                                        uy0 = face[l][1];
                                        uz0 = face[l][2];

                                        for (rn = 0; rn <= 2; rn++)
                                        {
                                            th = thij[ii][kk][rn];
                                            a1 = axij[ii][kk][rn][0];
                                            a2 = axij[ii][kk][rn][1];
                                            a3 = axij[ii][kk][rn][2];

                                            ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                            uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                            uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

                                            ux0 = ux;
                                            uy0 = uy;
                                            uz0 = uz;
                                        }

                                        uu = ux * ux + uy * uy + uz * uz;

                                        al = acos(fabs(nx * ux + ny * uy + nz * uz) / sqrt(uu));

                                        if (l == 0)
                                        {
                                            min_val = al;
                                            min_idx = l;
                                        }
                                        else
                                        {
                                            if (min_val > al)
                                            {
                                                min_val = al;
                                                min_idx = l;
                                            }
                                        }
                                    }

                                    ux0 = face[min_idx][0];
                                    uy0 = face[min_idx][1];
                                    uz0 = face[min_idx][2];

                                    for (rn = 0; rn <= 2; rn++)
                                    {
                                        th = thij[ii][kk][rn];
                                        a1 = axij[ii][kk][rn][0];
                                        a2 = axij[ii][kk][rn][1];
                                        a3 = axij[ii][kk][rn][2];

                                        ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                        uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                        uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

                                        ux0 = ux;
                                        uy0 = uy;
                                        uz0 = uz;
                                    }

                                    uu = ux * ux + uy * uy + uz * uz;

                                    al = acos((nx * ux + ny * uy + nz * uz) / sqrt(uu));

                                    epsilon0 = sqrt(aij[ii][kk]);

                                    PP = nx * ny * (ux * uy) + ny * nz * (uy * uz) + nx * nz * (ux * uz);
                                    QQ = pow(nx * ux, 2.0) + pow(ny * uy, 2.0) + pow(nz * uz, 2.0);

                                    CC = sqrt((2.0 * PP + QQ) / uu + rp0 * rp0);
                                    SS = sqrt((uu - QQ - 2.0 * PP) / uu + rp0 * rp0);

                                    ep = epsilon0 * (1.0 + del * (CC + tan(al0) * SS)) * coef;

                                    nxx = -phidx / phiabs / phiabs2 * dphiabs2dx / 2.0 + phidxx / phiabs;
                                    nyx = -phidy / phiabs / phiabs2 * dphiabs2dx / 2.0 + phidxy / phiabs;
                                    nzx = -phidz / phiabs / phiabs2 * dphiabs2dx / 2.0 + phidxz / phiabs;

                                    nxy = -phidx / phiabs / phiabs2 * dphiabs2dy / 2.0 + phidxy / phiabs;
                                    nyy = -phidy / phiabs / phiabs2 * dphiabs2dy / 2.0 + phidyy / phiabs;
                                    nzy = -phidz / phiabs / phiabs2 * dphiabs2dy / 2.0 + phidyz / phiabs;

                                    nxz = -phidx / phiabs / phiabs2 * dphiabs2dz / 2.0 + phidxz / phiabs;
                                    nyz = -phidy / phiabs / phiabs2 * dphiabs2dz / 2.0 + phidyz / phiabs;
                                    nzz = -phidz / phiabs / phiabs2 * dphiabs2dz / 2.0 + phidzz / phiabs;

                                    nxphix = -phidx * dphiabs2dphix / phiabs / phiabs2 / 2.0 + 1.0 / phiabs;
                                    nyphix = -phidy * dphiabs2dphix / phiabs / phiabs2 / 2.0;
                                    nzphix = -phidz * dphiabs2dphix / phiabs / phiabs2 / 2.0;

                                    nxphiy = -phidx * dphiabs2dphiy / phiabs / phiabs2 / 2.0;
                                    nyphiy = -phidy * dphiabs2dphiy / phiabs / phiabs2 / 2.0 + 1.0 / phiabs;
                                    nzphiy = -phidz * dphiabs2dphiy / phiabs / phiabs2 / 2.0;

                                    nxphiz = -phidx * dphiabs2dphiz / phiabs / phiabs2 / 2.0;
                                    nyphiz = -phidy * dphiabs2dphiz / phiabs / phiabs2 / 2.0;
                                    nzphiz = -phidz * dphiabs2dphiz / phiabs / phiabs2 / 2.0 + 1.0 / phiabs;

                                    nxphixdx = -0.5 / pow(phiabs, 6.0) * ((phidxx * dphiabs2dphix + phidx * dphiabs2dphixdx) * pow(phiabs, 3.0) - 1.5 * phidx * dphiabs2dphix * phiabs * dphiabs2dx) - 0.5 / phiabs / phiabs2 * dphiabs2dx;
                                    nyphixdx = -0.5 / pow(phiabs, 6.0) * ((phidxy * dphiabs2dphix + phidy * dphiabs2dphixdx) * pow(phiabs, 3.0) - 1.5 * phidy * dphiabs2dphix * phiabs * dphiabs2dx);
                                    nzphixdx = -0.5 / pow(phiabs, 6.0) * ((phidxz * dphiabs2dphix + phidz * dphiabs2dphixdx) * pow(phiabs, 3.0) - 1.5 * phidz * dphiabs2dphix * phiabs * dphiabs2dx);

                                    nxphiydy = -0.5 / pow(phiabs, 6.0) * ((phidxy * dphiabs2dphiy + phidx * dphiabs2dphiydy) * pow(phiabs, 3.0) - 1.5 * phidx * dphiabs2dphiy * phiabs * dphiabs2dy);
                                    nyphiydy = -0.5 / pow(phiabs, 6.0) * ((phidyy * dphiabs2dphiy + phidy * dphiabs2dphiydy) * pow(phiabs, 3.0) - 1.5 * phidy * dphiabs2dphiy * phiabs * dphiabs2dy) - 0.5 / phiabs / phiabs2 * dphiabs2dy;
                                    nzphiydy = -0.5 / pow(phiabs, 6.0) * ((phidyz * dphiabs2dphiy + phidz * dphiabs2dphiydy) * pow(phiabs, 3.0) - 1.5 * phidz * dphiabs2dphiy * phiabs * dphiabs2dy);

                                    nxphizdz = -0.5 / pow(phiabs, 6.0) * ((phidxz * dphiabs2dphiz + phidx * dphiabs2dphizdz) * pow(phiabs, 3.0) - 1.5 * phidx * dphiabs2dphiz * phiabs * dphiabs2dz);
                                    nyphizdz = -0.5 / pow(phiabs, 6.0) * ((phidyz * dphiabs2dphiz + phidy * dphiabs2dphizdz) * pow(phiabs, 3.0) - 1.5 * phidy * dphiabs2dphiz * phiabs * dphiabs2dz);
                                    nzphizdz = -0.5 / pow(phiabs, 6.0) * ((phidzz * dphiabs2dphiz + phidz * dphiabs2dphizdz) * pow(phiabs, 3.0) - 1.5 * phidz * dphiabs2dphiz * phiabs * dphiabs2dz) - 0.5 / phiabs / phiabs2 * dphiabs2dz;

                                    dPdx = nxx * ny * (ux * uy) + nx * nyx * (ux * uy) + nyx * nz * (uy * uz) + ny * nzx * (uy * uz) + nxx * nz * (ux * uz) + nx * nzx * (ux * uz);
                                    dPdy = nxy * ny * (ux * uy) + nx * nyy * (ux * uy) + nyy * nz * (uy * uz) + ny * nzy * (uy * uz) + nxy * nz * (ux * uz) + nx * nzy * (ux * uz);
                                    dPdz = nxz * ny * (ux * uy) + nx * nyz * (ux * uy) + nyz * nz * (uy * uz) + ny * nzz * (uy * uz) + nxz * nz * (ux * uz) + nx * nzz * (ux * uz);

                                    dPdphix = nxphix * ny * (ux * uy) + nx * nyphix * (ux * uy) + nyphix * nz * (uy * uz) + ny * nzphix * (uy * uz) + nxphix * nz * (ux * uz) + nx * nzphix * (ux * uz);
                                    dPdphiy = nxphiy * ny * (ux * uy) + nx * nyphiy * (ux * uy) + nyphiy * nz * (uy * uz) + ny * nzphiy * (uy * uz) + nxphiy * nz * (ux * uz) + nx * nzphiy * (ux * uz);
                                    dPdphiz = nxphiz * ny * (ux * uy) + nx * nyphiz * (ux * uy) + nyphiz * nz * (uy * uz) + ny * nzphiz * (uy * uz) + nxphiz * nz * (ux * uz) + nx * nzphiz * (ux * uz);

                                    dPdphixdx = nxphixdx * ny * (ux * uy) + nxphix * nyx * (ux * uy) + nxx * nyphix * (ux * uy) + nx * nyphixdx * (ux * uy) + nyphixdx * nz * (uy * uz) + nyphix * nzx * (uy * uz) + nyx * nzphix * (uy * uz) + ny * nzphixdx * (uy * uz) + nxphixdx * nz * (ux * uz) + nxphix * nzx * (ux * uz) + nxx * nzphix * (ux * uz) + nx * nzphixdx * (ux * uz);
                                    dPdphiydy = nxphiydy * ny * (ux * uy) + nxphiy * nyy * (ux * uy) + nxy * nyphiy * (ux * uy) + nx * nyphiydy * (ux * uy) + nyphiydy * nz * (uy * uz) + nyphiy * nzy * (uy * uz) + nyy * nzphiy * (uy * uz) + ny * nzphiydy * (uy * uz) + nxphiydy * nz * (ux * uz) + nxphiy * nzy * (ux * uz) + nxy * nzphiy * (ux * uz) + nx * nzphiydy * (ux * uz);
                                    dPdphizdz = nxphizdz * ny * (ux * uy) + nxphiz * nyz * (ux * uy) + nxz * nyphiz * (ux * uy) + nx * nyphizdz * (ux * uy) + nyphizdz * nz * (uy * uz) + nyphiz * nzz * (uy * uz) + nyz * nzphiz * (uy * uz) + ny * nzphizdz * (uy * uz) + nxphizdz * nz * (ux * uz) + nxphiz * nzz * (ux * uz) + nxz * nzphiz * (ux * uz) + nx * nzphizdz * (ux * uz);

                                    dQdx = 2.0 * (ux * ux * nx * nxx + uy * uy * ny * nyx + uz * uz * nz * nzx);
                                    dQdy = 2.0 * (ux * ux * nx * nxy + uy * uy * ny * nyy + uz * uz * nz * nzy);
                                    dQdz = 2.0 * (ux * ux * nx * nxz + uy * uy * ny * nyz + uz * uz * nz * nzz);

                                    dQdphix = 2.0 * (ux * ux * nx * nxphix + uy * uy * ny * nyphix + uz * uz * nz * nzphix);
                                    dQdphiy = 2.0 * (ux * ux * nx * nxphiy + uy * uy * ny * nyphiy + uz * uz * nz * nzphiy);
                                    dQdphiz = 2.0 * (ux * ux * nx * nxphiz + uy * uy * ny * nyphiz + uz * uz * nz * nzphiz);

                                    dQdphixdx = 2.0 * (ux * ux * (nxx * nxphix + nx * nxphixdx) + uy * uy * (nyx * nyphix + ny * nyphixdx) + uz * uz * (nzx * nzphix + nz * nzphixdx));
                                    dQdphiydy = 2.0 * (ux * ux * (nxy * nxphiy + nx * nxphiydy) + uy * uy * (nyy * nyphiy + ny * nyphiydy) + uz * uz * (nzy * nzphiy + nz * nzphiydy));
                                    dQdphizdz = 2.0 * (ux * ux * (nxz * nxphiz + nx * nxphizdz) + uy * uy * (nyz * nyphiz + ny * nyphizdz) + uz * uz * (nzz * nzphiz + nz * nzphizdz));

                                    epdx = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdx + dQdx) * coef;
                                    epdy = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdy + dQdy) * coef;
                                    epdz = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdz + dQdz) * coef;

                                    epdphix = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphix + dQdphix) * coef;
                                    epdphiy = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphiy + dQdphiy) * coef;
                                    epdphiz = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphiz + dQdphiz) * coef;

                                    epdphixdx = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdx + dQdx) * (2.0 * dPdphix + dQdphix) + (2.0 * dPdphixdx + dQdphixdx) * (1.0 / CC - tan(al0) / SS)) * coef;
                                    epdphiydy = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdy + dQdy) * (2.0 * dPdphiy + dQdphiy) + (2.0 * dPdphiydy + dQdphiydy) * (1.0 / CC - tan(al0) / SS)) * coef;
                                    epdphizdz = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdz + dQdz) * (2.0 * dPdphiz + dQdphiz) + (2.0 * dPdphizdz + dQdphizdz) * (1.0 / CC - tan(al0) / SS)) * coef;

                                    termx = ep * ep * phidxx + 2.0 * ep * epdx * phidx + ep * epdx * dphiabs2dx + ep * epdphixdx * phiabs2 + epdx * epdphix * phiabs2;
                                    termy = ep * ep * phidyy + 2.0 * ep * epdy * phidy + ep * epdy * dphiabs2dy + ep * epdphiydy * phiabs2 + epdy * epdphiy * phiabs2;
                                    termz = ep * ep * phidzz + 2.0 * ep * epdz * phidz + ep * epdz * dphiabs2dz + ep * epdphizdz * phiabs2 + epdz * epdphiz * phiabs2;

                                    termiikk = termx + termy + termz;
                                }
                                else
                                {
                                    termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                                }

                                if (anij[jj][kk] == 1.0 && phiabs2 != 0.0)
                                {
                                    nx = phidx / phiabs;
                                    ny = phidy / phiabs;
                                    nz = phidz / phiabs;

                                    min_val = 0.0;
                                    min_idx = 0;
                                    for (l = 0; l < FN; l++)
                                    {

                                        ux0 = face[l][0];
                                        uy0 = face[l][1];
                                        uz0 = face[l][2];

                                        for (rn = 0; rn <= 2; rn++)
                                        {
                                            th = thij[jj][kk][rn];
                                            a1 = axij[jj][kk][rn][0];
                                            a2 = axij[jj][kk][rn][1];
                                            a3 = axij[jj][kk][rn][2];

                                            ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                            uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                            uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

                                            ux0 = ux;
                                            uy0 = uy;
                                            uz0 = uz;
                                        }

                                        uu = ux * ux + uy * uy + uz * uz;

                                        al = acos(fabs(nx * ux + ny * uy + nz * uz) / sqrt(uu));
                                        if (l == 0)
                                        {
                                            min_val = al;
                                            min_idx = l;
                                        }
                                        else
                                        {
                                            if (min_val > al)
                                            {
                                                min_val = al;
                                                min_idx = l;
                                            }
                                        }
                                    }

                                    ux0 = face[min_idx][0];
                                    uy0 = face[min_idx][1];
                                    uz0 = face[min_idx][2];

                                    for (rn = 0; rn <= 2; rn++)
                                    {
                                        th = thij[jj][kk][rn];
                                        a1 = axij[jj][kk][rn][0];
                                        a2 = axij[jj][kk][rn][1];
                                        a3 = axij[jj][kk][rn][2];

                                        ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                        uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                        uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

                                        ux0 = ux;
                                        uy0 = uy;
                                        uz0 = uz;
                                    }

                                    uu = ux * ux + uy * uy + uz * uz;

                                    al = acos((nx * ux + ny * uy + nz * uz) / sqrt(uu));

                                    epsilon0 = sqrt(aij[jj][kk]);

                                    PP = nx * ny * (ux * uy) + ny * nz * (uy * uz) + nx * nz * (ux * uz);
                                    QQ = pow(nx * ux, 2.0) + pow(ny * uy, 2.0) + pow(nz * uz, 2.0);

                                    CC = sqrt((2.0 * PP + QQ) / uu + rp0 * rp0);
                                    SS = sqrt((uu - QQ - 2.0 * PP) / uu + rp0 * rp0);

                                    ep = epsilon0 * (1.0 + del * (CC + tan(al0) * SS)) * coef;

                                    nxx = -phidx / phiabs / phiabs2 * dphiabs2dx / 2.0 + phidxx / phiabs;
                                    nyx = -phidy / phiabs / phiabs2 * dphiabs2dx / 2.0 + phidxy / phiabs;
                                    nzx = -phidz / phiabs / phiabs2 * dphiabs2dx / 2.0 + phidxz / phiabs;

                                    nxy = -phidx / phiabs / phiabs2 * dphiabs2dy / 2.0 + phidxy / phiabs;
                                    nyy = -phidy / phiabs / phiabs2 * dphiabs2dy / 2.0 + phidyy / phiabs;
                                    nzy = -phidz / phiabs / phiabs2 * dphiabs2dy / 2.0 + phidyz / phiabs;

                                    nxz = -phidx / phiabs / phiabs2 * dphiabs2dz / 2.0 + phidxz / phiabs;
                                    nyz = -phidy / phiabs / phiabs2 * dphiabs2dz / 2.0 + phidyz / phiabs;
                                    nzz = -phidz / phiabs / phiabs2 * dphiabs2dz / 2.0 + phidzz / phiabs;

                                    nxphix = -phidx * dphiabs2dphix / phiabs / phiabs2 / 2.0 + 1.0 / phiabs;
                                    nyphix = -phidy * dphiabs2dphix / phiabs / phiabs2 / 2.0;
                                    nzphix = -phidz * dphiabs2dphix / phiabs / phiabs2 / 2.0;

                                    nxphiy = -phidx * dphiabs2dphiy / phiabs / phiabs2 / 2.0;
                                    nyphiy = -phidy * dphiabs2dphiy / phiabs / phiabs2 / 2.0 + 1.0 / phiabs;
                                    nzphiy = -phidz * dphiabs2dphiy / phiabs / phiabs2 / 2.0;

                                    nxphiz = -phidx * dphiabs2dphiz / phiabs / phiabs2 / 2.0;
                                    nyphiz = -phidy * dphiabs2dphiz / phiabs / phiabs2 / 2.0;
                                    nzphiz = -phidz * dphiabs2dphiz / phiabs / phiabs2 / 2.0 + 1.0 / phiabs;

                                    nxphixdx = -0.5 / pow(phiabs, 6.0) * ((phidxx * dphiabs2dphix + phidx * dphiabs2dphixdx) * pow(phiabs, 3.0) - 1.5 * phidx * dphiabs2dphix * phiabs * dphiabs2dx) - 0.5 / phiabs / phiabs2 * dphiabs2dx;
                                    nyphixdx = -0.5 / pow(phiabs, 6.0) * ((phidxy * dphiabs2dphix + phidy * dphiabs2dphixdx) * pow(phiabs, 3.0) - 1.5 * phidy * dphiabs2dphix * phiabs * dphiabs2dx);
                                    nzphixdx = -0.5 / pow(phiabs, 6.0) * ((phidxz * dphiabs2dphix + phidz * dphiabs2dphixdx) * pow(phiabs, 3.0) - 1.5 * phidz * dphiabs2dphix * phiabs * dphiabs2dx);
                                    nxphiydy = -0.5 / pow(phiabs, 6.0) * ((phidxy * dphiabs2dphiy + phidx * dphiabs2dphiydy) * pow(phiabs, 3.0) - 1.5 * phidx * dphiabs2dphiy * phiabs * dphiabs2dy);
                                    nyphiydy = -0.5 / pow(phiabs, 6.0) * ((phidyy * dphiabs2dphiy + phidy * dphiabs2dphiydy) * pow(phiabs, 3.0) - 1.5 * phidy * dphiabs2dphiy * phiabs * dphiabs2dy) - 0.5 / phiabs / phiabs2 * dphiabs2dy;
                                    nzphiydy = -0.5 / pow(phiabs, 6.0) * ((phidyz * dphiabs2dphiy + phidz * dphiabs2dphiydy) * pow(phiabs, 3.0) - 1.5 * phidz * dphiabs2dphiy * phiabs * dphiabs2dy);

                                    nxphizdz = -0.5 / pow(phiabs, 6.0) * ((phidxz * dphiabs2dphiz + phidx * dphiabs2dphizdz) * pow(phiabs, 3.0) - 1.5 * phidx * dphiabs2dphiz * phiabs * dphiabs2dz);
                                    nyphizdz = -0.5 / pow(phiabs, 6.0) * ((phidyz * dphiabs2dphiz + phidy * dphiabs2dphizdz) * pow(phiabs, 3.0) - 1.5 * phidy * dphiabs2dphiz * phiabs * dphiabs2dz);
                                    nzphizdz = -0.5 / pow(phiabs, 6.0) * ((phidzz * dphiabs2dphiz + phidz * dphiabs2dphizdz) * pow(phiabs, 3.0) - 1.5 * phidz * dphiabs2dphiz * phiabs * dphiabs2dz) - 0.5 / phiabs / phiabs2 * dphiabs2dz;

                                    dPdx = nxx * ny * (ux * uy) + nx * nyx * (ux * uy) + nyx * nz * (uy * uz) + ny * nzx * (uy * uz) + nxx * nz * (ux * uz) + nx * nzx * (ux * uz);
                                    dPdy = nxy * ny * (ux * uy) + nx * nyy * (ux * uy) + nyy * nz * (uy * uz) + ny * nzy * (uy * uz) + nxy * nz * (ux * uz) + nx * nzy * (ux * uz);
                                    dPdz = nxz * ny * (ux * uy) + nx * nyz * (ux * uy) + nyz * nz * (uy * uz) + ny * nzz * (uy * uz) + nxz * nz * (ux * uz) + nx * nzz * (ux * uz);

                                    dPdphix = nxphix * ny * (ux * uy) + nx * nyphix * (ux * uy) + nyphix * nz * (uy * uz) + ny * nzphix * (uy * uz) + nxphix * nz * (ux * uz) + nx * nzphix * (ux * uz);
                                    dPdphiy = nxphiy * ny * (ux * uy) + nx * nyphiy * (ux * uy) + nyphiy * nz * (uy * uz) + ny * nzphiy * (uy * uz) + nxphiy * nz * (ux * uz) + nx * nzphiy * (ux * uz);
                                    dPdphiz = nxphiz * ny * (ux * uy) + nx * nyphiz * (ux * uy) + nyphiz * nz * (uy * uz) + ny * nzphiz * (uy * uz) + nxphiz * nz * (ux * uz) + nx * nzphiz * (ux * uz);

                                    dPdphixdx = nxphixdx * ny * (ux * uy) + nxphix * nyx * (ux * uy) + nxx * nyphix * (ux * uy) + nx * nyphixdx * (ux * uy) + nyphixdx * nz * (uy * uz) + nyphix * nzx * (uy * uz) + nyx * nzphix * (uy * uz) + ny * nzphixdx * (uy * uz) + nxphixdx * nz * (ux * uz) + nxphix * nzx * (ux * uz) + nxx * nzphix * (ux * uz) + nx * nzphixdx * (ux * uz);
                                    dPdphiydy = nxphiydy * ny * (ux * uy) + nxphiy * nyy * (ux * uy) + nxy * nyphiy * (ux * uy) + nx * nyphiydy * (ux * uy) + nyphiydy * nz * (uy * uz) + nyphiy * nzy * (uy * uz) + nyy * nzphiy * (uy * uz) + ny * nzphiydy * (uy * uz) + nxphiydy * nz * (ux * uz) + nxphiy * nzy * (ux * uz) + nxy * nzphiy * (ux * uz) + nx * nzphiydy * (ux * uz);
                                    dPdphizdz = nxphizdz * ny * (ux * uy) + nxphiz * nyz * (ux * uy) + nxz * nyphiz * (ux * uy) + nx * nyphizdz * (ux * uy) + nyphizdz * nz * (uy * uz) + nyphiz * nzz * (uy * uz) + nyz * nzphiz * (uy * uz) + ny * nzphizdz * (uy * uz) + nxphizdz * nz * (ux * uz) + nxphiz * nzz * (ux * uz) + nxz * nzphiz * (ux * uz) + nx * nzphizdz * (ux * uz);

                                    dQdx = 2.0 * (ux * ux * nx * nxx + uy * uy * ny * nyx + uz * uz * nz * nzx);
                                    dQdy = 2.0 * (ux * ux * nx * nxy + uy * uy * ny * nyy + uz * uz * nz * nzy);
                                    dQdz = 2.0 * (ux * ux * nx * nxz + uy * uy * ny * nyz + uz * uz * nz * nzz);

                                    dQdphix = 2.0 * (ux * ux * nx * nxphix + uy * uy * ny * nyphix + uz * uz * nz * nzphix);
                                    dQdphiy = 2.0 * (ux * ux * nx * nxphiy + uy * uy * ny * nyphiy + uz * uz * nz * nzphiy);
                                    dQdphiz = 2.0 * (ux * ux * nx * nxphiz + uy * uy * ny * nyphiz + uz * uz * nz * nzphiz);

                                    dQdphixdx = 2.0 * (ux * ux * (nxx * nxphix + nx * nxphixdx) + uy * uy * (nyx * nyphix + ny * nyphixdx) + uz * uz * (nzx * nzphix + nz * nzphixdx));
                                    dQdphiydy = 2.0 * (ux * ux * (nxy * nxphiy + nx * nxphiydy) + uy * uy * (nyy * nyphiy + ny * nyphiydy) + uz * uz * (nzy * nzphiy + nz * nzphiydy));
                                    dQdphizdz = 2.0 * (ux * ux * (nxz * nxphiz + nx * nxphizdz) + uy * uy * (nyz * nyphiz + ny * nyphizdz) + uz * uz * (nzz * nzphiz + nz * nzphizdz));

                                    epdx = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdx + dQdx) * coef;
                                    epdy = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdy + dQdy) * coef;
                                    epdz = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdz + dQdz) * coef;

                                    epdphix = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphix + dQdphix) * coef;
                                    epdphiy = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphiy + dQdphiy) * coef;
                                    epdphiz = epsilon0 * del * (1.0 / (2.0 * uu * CC) - tan(al0) / (2.0 * uu * SS)) * (2.0 * dPdphiz + dQdphiz) * coef;

                                    epdphixdx = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdx + dQdx) * (2.0 * dPdphix + dQdphix) + (2.0 * dPdphixdx + dQdphixdx) * (1.0 / CC - tan(al0) / SS)) * coef;
                                    epdphiydy = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdy + dQdy) * (2.0 * dPdphiy + dQdphiy) + (2.0 * dPdphiydy + dQdphiydy) * (1.0 / CC - tan(al0) / SS)) * coef;
                                    epdphizdz = 1.0 / 2.0 / uu * epsilon0 * del * ((-1.0 / 2.0 / uu / pow(CC, 3.0) - tan(al0) / (2.0 * uu * pow(SS, 3.0))) * (2.0 * dPdz + dQdz) * (2.0 * dPdphiz + dQdphiz) + (2.0 * dPdphizdz + dQdphizdz) * (1.0 / CC - tan(al0) / SS)) * coef;

                                    termx = ep * ep * phidxx + 2.0 * ep * epdx * phidx + ep * epdx * dphiabs2dx + ep * epdphixdx * phiabs2 + epdx * epdphix * phiabs2;
                                    termy = ep * ep * phidyy + 2.0 * ep * epdy * phidy + ep * epdy * dphiabs2dy + ep * epdphiydy * phiabs2 + epdy * epdphiy * phiabs2;
                                    termz = ep * ep * phidzz + 2.0 * ep * epdz * phidz + ep * epdz * dphiabs2dz + ep * epdphizdz * phiabs2 + epdz * epdphiz * phiabs2;

                                    termjjkk = termx + termy + termz;
                                }
                                else
                                {
                                    termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                                }

                                sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
                            }

                            phidxii = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0;
                            phidyii = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0;
                            phidzii = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0;
                            phiabs2ii = phidxii * phidxii + phidyii * phidyii + phidzii * phidzii;
                            if (anij[ii][jj] == 1.0 && phiabs2ii != 0.0)
                            {

                                nxii = phidxii / sqrt(phiabs2ii);
                                nyii = phidyii / sqrt(phiabs2ii);
                                nzii = phidzii / sqrt(phiabs2ii);

                                min_val = 0.0;
                                min_idx = 0;
                                for (l = 0; l < FN; l++)
                                {
                                    ux0 = face[l][0];
                                    uy0 = face[l][1];
                                    uz0 = face[l][2];

                                    for (rn = 0; rn <= 2; rn++)
                                    {
                                        th = thij[ii][jj][rn];
                                        a1 = axij[ii][jj][rn][0];
                                        a2 = axij[ii][jj][rn][1];
                                        a3 = axij[ii][jj][rn][2];

                                        ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                        uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                                        uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

                                        ux0 = ux;
                                        uy0 = uy;
                                        uz0 = uz;
                                    }

                                    uu = ux * ux + uy * uy + uz * uz;
                                    al = acos(fabs(nxii * ux + nyii * uy + nzii * uz) / sqrt(uu));

                                    if (l == 0)
                                    {
                                        min_val = al;
                                        min_idx = l;
                                    }
                                    else
                                    {
                                        if (min_val > al)
                                        {
                                            min_val = al;
                                            min_idx = l;
                                        }
                                    }
                                }

                                if (min_idx <= 3 && min_val < 90.0 / ww1 / 180.0 * PI)
                                {
                                    miijj = mij[ii][jj] * (zeta1 + (1 - zeta1) * sqrt(pow(tan(ww1 * min_val), 2.0) + rp1 * rp1) * tanh(1.0 / sqrt(pow(tan(ww1 * min_val), 2.0) + rp1 * rp1)));
                                }
                                else
                                {
                                    miijj = mij[ii][jj];
                                }

                                // add interface mobility to the interface field
                                // if (phiNum[i][j][k] == 2)
                                if (phi[0][i][j][k] > 0.2 && phi[0][i][j][k] < 0.8 && phiNum[i][j][k] == 2)
                                // if (((*phi)[0][i][j][k] * (*phi)[1][i][j][k] > 0.0 || (*phi)[0][i][j][k] * (*phi)[2][i][j][k] > 0.0) && (*phi)[1][i][j][k] * (*phi)[2][i][j][k] == 0.0)
                                {
                                    mob[i][j][k] = miijj / mij[ii][jj]; // min_val / PI * 180.0; //
                                }
                            }
                            else
                            {
                                miijj = mij[ii][jj];
                            }

                            if (ii > 0 && jj == 0)
                            {
                                F0 = -(temp[i][j][k] - Tm);
                            }
                            else if (ii == 0 && jj > 0)
                            {
                                F0 = (temp[i][j][k] - Tm);
                            }
                            else
                            {
                                F0 = 0.0;
                            }
                            pddtt += -2.0 * miijj / (double)(phiNum[i][j][k]) * (sum1 - 8.0 / PI * F0 * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
                        }
                        phi2[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime;
                        if (phi2[ii][i][j][k] >= 1.0)
                        {
                            phi2[ii][i][j][k] = 1.0;
                        }
                        if (phi2[ii][i][j][k] <= 0.0)
                        {
                            phi2[ii][i][j][k] = 0.0;
                        }
                    }
                } // k
            }     // j
        }         // i

#pragma omp barrier

        for (ii = 0; ii <= nm; ii++)
        {
            for (i = start; i <= end; i++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    for (k = 0; k <= ndmz; k++)
                    {
                        phi[ii][i][j][k] = phi2[ii][i][j][k];
                    }
                }
            }
        }

        for (i = start; i <= end; i++)
        {
            for (j = 0; j <= ndmy; j++)
            {
                for (k = 0; k <= ndmz; k++)
                {
                    sum1 = 0.0;
                    for (ii = 0; ii <= nm; ii++)
                    {
                        sum1 += phi[ii][i][j][k];
                    }
                    for (ii = 0; ii <= nm; ii++)
                    {
                        phi[ii][i][j][k] = phi[ii][i][j][k] / sum1;
                    }
                }
            }
        }

#pragma omp barrier

        istep = istep + 1;
        if (istep < nstep)
        {
            goto start;
        }
    end:;
    }
    return 0;
}

void datasave2d(int step)
{
    int is, js;
    FILE *stream; // ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/phi/2d%d.csv", step);
    stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン

    for (is = 0; is <= ndmx; is++)
    {
        for (js = 0; js <= ndmy; js++)
        {
            fprintf(stream, "%lf   ", phi[0][is][js][0] * (1.0 - phi[0][is][js][0]));
            fprintf(stream, "\n");
        }
    }
    fclose(stream); // ファイルをクローズ

    FILE *streamm; // ストリームのポインタ設定
    char bufferm[30];
    sprintf(bufferm, "data/mob/2d%d.csv", step);
    streamm = fopen(bufferm, "a"); // 書き込む先のファイルを追記方式でオープン

    for (is = 0; is <= ndmx; is++)
    {
        for (js = 0; js <= ndmy; js++)
        {
            fprintf(streamm, "%lf   ", mob[is][js][0]);
            fprintf(streamm, "\n");
        }
    }
    fclose(streamm); // ファイルをクローズ
}
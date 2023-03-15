
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <omp.h>

using namespace std;

#define N 3
#define NDX 128
#define NDY 128
#define NDZ 1
#define PI 3.14159
#define RR 8.314
#define NTH 8

int nm = N - 1;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int rows = NDX / NTH;

int nstep = 6001;
int pstep = 2000;

double dx = 1.0e-7;
double gamma0 = 0.245;
double astre = 0.0192;
double astrem = 0.1;
double delta = 5.0 * dx;
double mobi = 0.1e-8;

double A0 = 8.0 * delta * gamma0 / PI / PI;
double W0 = 4.0 * gamma0 / delta;
double M0 = mobi * PI * PI / (8.0 * delta);
double S0 = 1.475e6;

double gradT = 0.0e5;
double rateT = 0.0e-3;
double temp0 = 922.5 - NDX / 16 * dx * gradT;
double cl = 0.01;

double Dl = 1.34e-7 * exp(-3.0e4 / RR / temp0);
double Ds = 0.89e-4 * exp(-1.36e5 / RR / temp0);

double dtime = 0.2 * dx * dx / Dl;

// Linear phae diagram
double Te = 859.6;
double ce = 0.122;
double ml1 = -605.3;
double kap1 = 0.105;
double kap2 = 1.8;
double ml2 = 1302.0;

double mij[N][N], aij[N][N], wij[N][N], anij[N][N];
double thij[N][N], vpij[N][N], etaij[N][N];

double
    ****phi,
    ****phi2,
    ****conp,
    ****conp2,
    ****dphi0,
    ***con,
    ***temp,
    ***intphi;
int
    ***phiNum,
    ****phiIdx;

int ni, nj, ix, iy, iz;

double cc, c0, c00, dc0;
double cm0, dcm0, sumc, suma, suma0, sumb, sumil;
double tpv, tpdphidx, tpdphidt, tpn;
int tpx, tpy, tpxp, tpxm;

void fracsave(int step);
void datasave(int step);
double calC01e(double temp);
double calC1e(double temp);
double calC02e(double temp);
double calC2e(double temp);

int main(void)
{
    cout << "----------------------------------------------" << endl;
    cout << "Computation Started!" << endl;
    cout << "concenration field stablity number is: " << dtime * Dl / dx / dx << endl;
    cout << "phase field stability number is: " << dtime / dx / dx * M0 * A0 << endl;
    cout << "input mobility is " << mobi << endl;

    phi = new double ***[N];
    phi2 = new double ***[N];
    conp = new double ***[N];
    conp2 = new double ***[N];
    for (ni = 0; ni <= nm; ni++)
    {
        phi[ni] = new double **[NDX];
        phi2[ni] = new double **[NDX];
        conp[ni] = new double **[NDX];
        conp2[ni] = new double **[NDX];
        for (ix = 0; ix <= ndmx; ix++)
        {
            phi[ni][ix] = new double *[NDY];
            phi2[ni][ix] = new double *[NDY];
            conp[ni][ix] = new double *[NDY];
            conp2[ni][ix] = new double *[NDY];
            for (iy = 0; iy <= ndmy; iy++)
            {
                phi[ni][ix][iy] = new double[NDZ];
                phi2[ni][ix][iy] = new double[NDZ];
                conp[ni][ix][iy] = new double[NDZ];
                conp2[ni][ix][iy] = new double[NDZ];
            }
        }
    }

    dphi0 = new double ***[3];
    for (ni = 0; ni <= nm; ni++)
    {
        dphi0[ni] = new double **[NDX];
        for (ix = 0; ix <= ndmx; ix++)
        {
            dphi0[ni][ix] = new double *[NDY];
            for (iy = 0; iy <= ndmy; iy++)
            {
                dphi0[ni][ix][iy] = new double[NDZ];
            }
        }
    }

    con = new double **[NDX];
    temp = new double **[NDX];
    intphi = new double **[NDX];
    for (ix = 0; ix <= ndmx; ix++)
    {
        con[ix] = new double *[NDY];
        temp[ix] = new double *[NDY];
        intphi[ix] = new double *[NDY];
        for (iy = 0; iy <= ndmy; iy++)
        {
            con[ix][iy] = new double[NDZ];
            temp[ix][iy] = new double[NDZ];
            intphi[ix][iy] = new double[NDZ];
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

    for (ni = 0; ni <= nm; ni++)
    {
        for (nj = 0; nj <= nm; nj++)
        {
            wij[ni][nj] = W0;
            aij[ni][nj] = A0;
            mij[ni][nj] = M0;
            thij[ni][nj] = 0.0;
            vpij[ni][nj] = 0.0;
            etaij[ni][nj] = 0.0;
            if (ni == nj)
            {
                wij[ni][nj] = 0.0;
                aij[ni][nj] = 0.0;
                mij[ni][nj] = 0.0;
            }
        }
    }

    anij[1][0] = 1.0;
    anij[0][1] = 1.0;
    // thij[1][0] = PI / 4;
    // thij[0][1] = PI / 4;

    sumc = 0.0;
    for (ix = 0; ix <= ndmx; ix++)
    {
        for (iy = 0; iy <= ndmy; iy++)
        {
            for (iz = 0; iz <= ndmz; iz++)
            {
                temp[ix][iy][iz] = temp0 + gradT * ix * dx;
                // if (ix <= NDX / 16) // && iy <= NDY * (1.0 - cl))
                if ((iy) * (iy) + (ix) * (ix) <= NDX / 64 * NDX / 64)
                {
                    phi[1][ix][iy][iz] = 1.0;
                    conp[1][ix][iy][iz] = calC1e(temp[ix][iy][iz]);
                    phi[2][ix][iy][iz] = 0.0;
                    conp[2][ix][iy][iz] = calC2e(temp[ix][iy][iz]);
                    phi[0][ix][iy][iz] = 0.0;
                    conp[0][ix][iy][iz] = calC01e(temp[ix][iy][iz]);
                }
                // else if (ix <= NDX / 16 && iy > NDY * (1.0 - cl))
                // else if ((j - NDY / 4) * (j - NDY / 4) + (i - NDX / 4) * (i - NDX / 4) <= NDX / 16 * NDX / 16)
                // {
                //     phi[1][ix][iy][iz] = 0.0;
                //     conp[1][ix][iy][iz] = calC1e(temp[ix][iy][iz]);
                //     phi[2][ix][iy][iz] = 1.0;
                //     conp[2][ix][iy][iz] = calC2e(temp[ix][iy][iz]);
                //     phi[0][ix][iy][iz] = 0.0;
                //     conp[0][ix][iy][iz] = calC02e(temp[ix][iy][iz]);
                // }
                else
                {
                    phi[1][ix][iy][iz] = 0.0;
                    conp[1][ix][iy][iz] = calC1e(temp[ix][iy][iz]);
                    phi[2][ix][iy][iz] = 0.0;
                    conp[2][ix][iy][iz] = calC2e(temp[ix][iy][iz]);
                    phi[0][ix][iy][iz] = 1.0;
                    conp[0][ix][iy][iz] = cl;
                }
                con[ix][iy][iz] = conp[1][ix][iy][iz] * phi[1][ix][iy][iz] + conp[2][ix][iy][iz] * phi[2][ix][iy][iz] + conp[0][ix][iy][iz] * phi[0][ix][iy][iz];
                sumc += con[ix][iy][iz];
            }
        }
    }
    cc = sumc / NDX / NDY / NDZ;
    cm0 = sumc;

#pragma omp parallel num_threads(NTH)
    {
        int th_id, istep;
        int start, end, offset;

        th_id = omp_get_thread_num();

        istep = 0;
        offset = th_id * rows;
        start = offset;
        end = offset + rows - 1;

        int phinum;

        int i, j, im, ip, jp, jm, kp, km, k;
        int ii, jj, kk, di, dim;
        int n1, n2, n3;

        double dF, pddtt, sum1;
        double termiikk, termjjkk;
        double miijj;
        double gc0i, c0i;
        double phidxi, phidyi;
        double nxi, nyi;
        int xdi, xdip, xdim, ydi, ydip, ydim;
        int ixdi, ixdip, jydi, jydip;
        double maxphi01;
        int phix, phiy;

        double phidx, phidy, phidz;
        double phidxx, phidyy, phidzz;
        double phidxy, phidxz, phidyz;
        double phiabs2;

        double th, vp, eta;
        double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;
        double phidxp, phidyp, phidzp;
        double phidxpx, phidypx, phidzpx;
        double phidxpy, phidypy, phidzpy;
        double phidxpz, phidypz, phidzpz;
        double ep, epdx, epdy, epdz;
        double epsilon0;
        double term0;
        double termx, termx0, termx1, termx0dx, termx1dx;
        double termy, termy0, termy1, termy0dy, termy1dy;
        double termz, termz0, termz1, termz0dz, termz1dz;

        double phidxm, phidym, phidzm;
        double phidxpm, phidypm, phidzpm, phiabsm2;

        double dphi01, dphi02;
        double Dl0, Ds0;
        double flxci;
        double cddttl1, cddttl2, cddtts1, cddtts2;
        double fel1, fwl1, fnl1, fsl1, fil1, fol1;
        double fel2, fwl2, fnl2, fsl2, fil2, fol2;
        double fes1, fws1, fns1, fss1, fis1, fos1;
        double fes2, fws2, fns2, fss2, fis2, fos2;
        double lfrac1, lfrac2;

    start:;

        if (istep % pstep == 0 && th_id == 0)
        {
            datasave(istep);
            cout << istep << " steps(" << istep * dtime << " seconds) has done!" << endl;
            cout << "The actual concnetration is                  " << cc << endl;
            cout << "The current concnetration is                 " << c0 << endl;
            cout << "The concnetration deviation is              " << dcm0 / NDX / NDY / NDZ << endl;
            cout << " Al " << suma / NDX / NDY / NDZ * 100 << "%" << endl;
            cout << " Si " << sumb / NDX / NDY / NDZ * 100 << "%" << endl;
            // cout << " Tip postion: (" << tpx << ", " << tpy << ")" << endl;
        }

        if (istep % 200 == 0 && th_id == 0)
        {
            fracsave(istep);
        }

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
                        }
                    }
                    phiNum[i][j][k] = phinum;
                    dphi0[2][i][j][k] = 0.0;
                    dphi0[1][i][j][k] = 0.0;
                }
            }
        }
#pragma omp barrier
        // Evolution Equation of Phase Fields
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
                    intphi[i][j][k] = 1.0;
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
                                phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 / dx;
                                phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 / dx;
                                phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 / dx;

                                phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx; // フェーズフィールドの空間２階微分
                                phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
                                phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

                                phidxy = (phi[kk][ip][jp][k] + phi[kk][im][jm][k] - phi[kk][im][jp][k] - phi[kk][ip][jm][k]) / 4.0 / dx / dx;
                                phidxz = (phi[kk][ip][j][kp] + phi[kk][im][j][km] - phi[kk][im][j][kp] - phi[kk][ip][j][km]) / 4.0 / dx / dx;
                                phidyz = (phi[kk][i][jp][kp] + phi[kk][i][jm][km] - phi[kk][i][jm][kp] - phi[kk][i][jp][km]) / 4.0 / dx / dx;

                                phiabs2 = phidx * phidx + phidy * phidy + phidz * phidz;

                                if (anij[ii][kk] == 1.0 && phiabs2 != 0.0)
                                {
                                    epsilon0 = sqrt(aij[ii][kk]);

                                    th = thij[ii][kk];
                                    vp = vpij[ii][kk];
                                    eta = etaij[ii][kk];

                                    xxp = cos(th) * cos(vp);
                                    yxp = sin(th) * cos(vp);
                                    zxp = sin(vp);
                                    xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                    yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                    zyp = cos(vp) * sin(eta);
                                    xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                    yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                    zzp = cos(eta) * cos(vp);

                                    phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                    phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                    phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                    phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                    phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                    phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                    phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                    phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                    phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                    phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                    phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                    phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                    ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 2.0));

                                    epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs2, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0));
                                    epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs2, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0));
                                    epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs2, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0));

                                    term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                    termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs2;
                                    termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs2;
                                    termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs2;

                                    termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs2, 2.0);
                                    termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs2, 2.0);
                                    termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs2, 2.0);

                                    termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs2 - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs2, 2.0);
                                    termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs2 - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs2, 2.0);
                                    termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs2 - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs2, 2.0);

                                    termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs2, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0);
                                    termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs2, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0);
                                    termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs2, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0);

                                    termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                    termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                    termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                    termiikk = term0 + termx + termy + termz;
                                }
                                else
                                {
                                    termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
                                }

                                if (anij[jj][kk] == 1.0 && phiabs2 != 0.0)
                                {
                                    epsilon0 = sqrt(aij[jj][kk]);

                                    th = thij[jj][kk];
                                    vp = vpij[jj][kk];
                                    eta = etaij[jj][kk];

                                    xxp = cos(th) * cos(vp);
                                    yxp = sin(th) * cos(vp);
                                    zxp = sin(vp);
                                    xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                    yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                    zyp = cos(vp) * sin(eta);
                                    xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                    yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                    zzp = cos(eta) * cos(vp);

                                    phidxp = phidx * xxp + phidy * yxp + phidz * zxp;
                                    phidyp = phidx * xyp + phidy * yyp + phidz * zyp;
                                    phidzp = phidx * xzp + phidy * yzp + phidz * zzp;

                                    phidxpx = phidxx * xxp + phidxy * yxp + phidxz * zxp;
                                    phidypx = phidxx * xyp + phidxy * yyp + phidxz * zyp;
                                    phidzpx = phidxx * xzp + phidxy * yzp + phidxz * zzp;

                                    phidxpy = phidxy * xxp + phidyy * yxp + phidyz * zxp;
                                    phidypy = phidxy * xyp + phidyy * yyp + phidyz * zyp;
                                    phidzpy = phidxy * xzp + phidyy * yzp + phidyz * zzp;

                                    phidxpz = phidxz * xxp + phidyz * yxp + phidzz * zxp;
                                    phidypz = phidxz * xyp + phidyz * yyp + phidzz * zyp;
                                    phidzpz = phidxz * xzp + phidyz * yzp + phidzz * zzp;

                                    ep = epsilon0 * (1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 2.0));

                                    epdx = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpx + pow(phidyp, 3.0) * phidypx + pow(phidzp, 3.0) * phidzpx) / pow(phiabs2, 2.0) - (phidx * phidxx + phidy * phidxy + phidz * phidxz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0));
                                    epdy = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpy + pow(phidyp, 3.0) * phidypy + pow(phidzp, 3.0) * phidzpy) / pow(phiabs2, 2.0) - (phidx * phidxy + phidy * phidyy + phidz * phidyz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0));
                                    epdz = 16.0 * epsilon0 * astre * ((pow(phidxp, 3.0) * phidxpz + pow(phidyp, 3.0) * phidypz + pow(phidzp, 3.0) * phidzpz) / pow(phiabs2, 2.0) - (phidx * phidxz + phidy * phidyz + phidz * phidzz) * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0));

                                    term0 = 2.0 * ep * epdx * phidx + phidxx * ep * ep + 2.0 * ep * epdy * phidy + phidyy * ep * ep + 2.0 * ep * epdz * phidz + phidzz * ep * ep;

                                    termx0 = (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / phiabs2;
                                    termy0 = (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / phiabs2;
                                    termz0 = (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / phiabs2;

                                    termx1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidx / pow(phiabs2, 2.0);
                                    termy1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidy / pow(phiabs2, 2.0);
                                    termz1 = (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) * phidz / pow(phiabs2, 2.0);

                                    termx0dx = (3.0 * pow(phidxp, 2.0) * phidxpx * xxp + 3.0 * pow(phidyp, 2.0) * phidypx * xyp + 3.0 * pow(phidzp, 2.0) * phidzpx * xzp) / phiabs2 - (2.0 * phidx * phidxx + 2.0 * phidy * phidxy + 2.0 * phidz * phidxz) * (pow(phidxp, 3.0) * xxp + pow(phidyp, 3.0) * xyp + pow(phidzp, 3.0) * xzp) / pow(phiabs2, 2.0);
                                    termy0dy = (3.0 * pow(phidxp, 2.0) * phidxpy * yxp + 3.0 * pow(phidyp, 2.0) * phidypy * yyp + 3.0 * pow(phidzp, 2.0) * phidzpy * yzp) / phiabs2 - (2.0 * phidx * phidxy + 2.0 * phidy * phidyy + 2.0 * phidz * phidyz) * (pow(phidxp, 3.0) * yxp + pow(phidyp, 3.0) * yyp + pow(phidzp, 3.0) * yzp) / pow(phiabs2, 2.0);
                                    termz0dz = (3.0 * pow(phidxp, 2.0) * phidxpz * zxp + 3.0 * pow(phidyp, 2.0) * phidypz * zyp + 3.0 * pow(phidzp, 2.0) * phidzpz * zzp) / phiabs2 - (2.0 * phidx * phidxz + 2.0 * phidy * phidyz + 2.0 * phidz * phidzz) * (pow(phidxp, 3.0) * zxp + pow(phidyp, 3.0) * zyp + pow(phidzp, 3.0) * zzp) / pow(phiabs2, 2.0);

                                    termx1dx = ((phidxx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidx * (4.0 * pow(phidxp, 3.0) * phidxpx + 4.0 * pow(phidyp, 3.0) * phidypx + 4.0 * pow(phidzp, 3.0) * phidzpx))) / pow(phiabs2, 2.0) - 4.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz) * phidx * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0);
                                    termy1dy = ((phidyy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidy * (4.0 * pow(phidxp, 3.0) * phidxpy + 4.0 * pow(phidyp, 3.0) * phidypy + 4.0 * pow(phidzp, 3.0) * phidzpy))) / pow(phiabs2, 2.0) - 4.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz) * phidy * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0);
                                    termz1dz = ((phidzz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) + phidz * (4.0 * pow(phidxp, 3.0) * phidxpz + 4.0 * pow(phidyp, 3.0) * phidypz + 4.0 * pow(phidzp, 3.0) * phidzpz))) / pow(phiabs2, 2.0) - 4.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz) * phidz * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 3.0);

                                    termx = 16.0 * epsilon0 * astre * (epdx * (termx0 - termx1) + ep * (termx0dx - termx1dx));
                                    termy = 16.0 * epsilon0 * astre * (epdy * (termy0 - termy1) + ep * (termy0dy - termy1dy));
                                    termz = 16.0 * epsilon0 * astre * (epdz * (termz0 - termz1) + ep * (termz0dz - termz1dz));

                                    termjjkk = term0 + termx + termy + termz;
                                }
                                else
                                {
                                    termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
                                }

                                sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
                            }
                            if (ii == 1 && jj == 0)
                            {
                                dF = ((conp[0][i][j][k] - ce) * ml1 + Te - temp[i][j][k]) * S0;
                            }
                            else if (ii == 0 && jj == 1)
                            {
                                dF = -(((conp[0][i][j][k] - ce) * ml1 + Te - temp[i][j][k])) * S0;
                            }
                            else if (ii == 2 && jj == 0)
                            {
                                dF = ((conp[0][i][j][k] - ce) * ml2 + Te - temp[i][j][k]) * S0;
                            }
                            else if (ii == 0 && jj == 2)
                            {
                                dF = -(((conp[0][i][j][k] - ce) * ml2 + Te - temp[i][j][k])) * S0;
                            }
                            else
                            {
                                dF = 0.0;
                            }
                            miijj = mij[ii][jj];
                            if (ii + jj == 3 && phi[0][i][j][k] == 0.0)
                            {
                                miijj = 0.01 * mij[ii][jj];
                            }
                            if ((ii == 1 && jj == 0) || (ii == 0 && jj == 1))
                            {
                                th = thij[ii][jj];
                                vp = vpij[ii][jj];
                                eta = etaij[ii][jj];

                                xxp = cos(th) * cos(vp);
                                yxp = sin(th) * cos(vp);
                                zxp = sin(vp);
                                xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
                                yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
                                zyp = cos(vp) * sin(eta);
                                xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
                                yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
                                zzp = cos(eta) * cos(vp);

                                phidxm = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0 / dx;
                                phidym = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0 / dx;
                                phidzm = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0 / dx;

                                phidxpm = phidxm * xxp + phidym * yxp + phidzm * zxp;
                                phidypm = phidxm * xyp + phidym * yyp + phidzm * zyp;
                                phidzpm = phidxm * xzp + phidym * yzp + phidzm * zzp;

                                phiabsm2 = phidxpm * phidxpm + phidypm * phidypm + phidzpm * phidzpm;

                                miijj = miijj * (1.0 - 3.0 * astrem + 4.0 * astrem * (pow(phidxpm, 4.0) + pow(phidypm, 4.0) + pow(phidzpm, 4.0)) / pow(phiabsm2, 2.0));
                            }
                            pddtt += -2.0 * miijj / double(phiNum[i][j][k]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));

                            if (ii == 0 && jj == 1)
                            {
                                dphi01 = -2.0 * miijj / double(phiNum[i][j][k]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
                            }
                            if (ii == 0 && jj == 2)
                            {
                                dphi02 = -2.0 * miijj / double(phiNum[i][j][k]) * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
                            }
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
                        if (ii == 0)
                        {
                            if (dphi01 != 0.0 && dphi02 != 0.0)
                            {
                                dphi0[1][i][j][k] = (phi2[0][i][j][k] - phi[0][i][j][k]) * dphi01 / (dphi01 + dphi02);
                                dphi0[2][i][j][k] = (phi2[0][i][j][k] - phi[0][i][j][k]) * dphi02 / (dphi01 + dphi02);
                            }
                            else if (dphi01 != 0.0 && dphi02 == 0.0)
                            {
                                dphi0[1][i][j][k] = (phi2[0][i][j][k] - phi[0][i][j][k]);
                            }
                            else if (dphi01 == 0.0 && dphi02 != 0.0)
                            {
                                dphi0[2][i][j][k] = (phi2[0][i][j][k] - phi[0][i][j][k]);
                            }
                            dphi01 = 0.0;
                            dphi02 = 0.0;
                        }
                    }
                }
            } // j
        }     // i

#pragma omp barrier
        // Partition and Ejection
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
                    if (phi[0][i][j][k] > 0.0 && phi[0][i][j][k] < 1.0)
                    {
                        phidxi = (phi[0][ip][j][k] - phi[0][im][j][k]) / 2.0 / dx;
                        phidyi = (phi[0][i][jp][k] - phi[0][i][jm][k]) / 2.0 / dx;
                        nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi);
                        nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi);

                        flxci = dphi0[1][i][j][k] * (calC1e(temp[i][j][k]) - conp[0][i][j][k]) + dphi0[2][i][j][k] * (calC2e(temp[i][j][k]) - conp[0][i][j][k]);
                        con[i][j][k] -= flxci;

                        di = 0;
                        do
                        {
                            di++;
                            xdi = int(round(nxi * di));
                            xdip = int(round(nxi * (di + 1)));
                            ydi = int(round(nyi * di));
                            ydip = int(round(nyi * (di + 1)));

                            ixdi = i + xdi;
                            if (ixdi > ndmx)
                            {
                                ixdi = ixdi - ndmx - 1;
                            }
                            if (ixdi < 0)
                            {
                                ixdi = ndmx + ixdi + 1;
                            }
                            ixdip = i + xdip;
                            if (ixdip > ndmx)
                            {
                                ixdip = ixdip - ndmx - 1;
                            }
                            if (ixdip < 0)
                            {
                                ixdip = ndmx + ixdip + 1;
                            }

                            jydi = j + ydi;
                            if (jydi > ndmy)
                            {
                                jydi = jydi - ndmy - 1;
                            }
                            if (jydi < 0)
                            {
                                jydi = ndmy + jydi + 1;
                            }

                            jydip = j + ydip;
                            if (jydip > ndmy)
                            {
                                jydip = jydip - ndmy - 1;
                            }
                            if (jydip < 0)
                            {
                                jydip = ndmy + jydip + 1;
                            }

                            if (phi[0][ixdi][jydi][k] == 1.0)
                            {
                                con[ixdi][jydi][k] += flxci;
                                break;
                            }
                        } while (phi[0][ixdi][jydi][k] < 1.0);
                    }
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
                    for (ii = 0; ii <= nm; ii++)
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

        // Calculate the local concentration
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
                    if (phi[0][i][j][k] == 0.0)
                    {
                        conp[1][i][j][k] = con[i][j][k] * phi[1][i][j][k];
                        conp[2][i][j][k] = con[i][j][k] * phi[2][i][j][k];
                        if (phi[1][i][j][k] < 1.0 && phi[2][i][j][k] < 1.0)
                        {
                            conp[1][i][j][k] = calC1e(temp[i][j][k]);
                            conp[2][i][j][k] = calC2e(temp[i][j][k]);
                        }
                        conp[0][i][j][k] = calC01e(temp[i][j][k]) * phi[1][i][j][k] + calC02e(temp[i][j][k]) * phi[2][i][j][k];
                    }
                    else if (phi[0][i][j][k] > 0.0 && phi[0][i][j][k] < 1.0)
                    {
                        conp[1][i][j][k] = calC1e(temp[i][j][k]);
                        conp[2][i][j][k] = calC2e(temp[i][j][k]);
                        conp[0][i][j][k] = (con[i][j][k] - phi[1][i][j][k] * conp[1][i][j][k] - phi[2][i][j][k] * conp[2][i][j][k]) / phi[0][i][j][k];
                    }
                    else if (phi[0][i][j][k] == 1.0)
                    {
                        conp[1][i][j][k] = calC1e(temp[i][j][k]);
                        conp[2][i][j][k] = calC2e(temp[i][j][k]);
                        conp[0][i][j][k] = con[i][j][k];
                    }
                    con[i][j][k] = conp[1][i][j][k] * phi[1][i][j][k] + conp[2][i][j][k] * phi[2][i][j][k] + conp[0][i][j][k] * phi[0][i][j][k];
                    if (con[i][j][k] > 1.0)
                    {
                        con[i][j][k] = 1.0;
                    }
                    if (con[i][j][k] < 0.0)
                    {
                        con[i][j][k] = 0.0;
                    }
                }
            }
        }

#pragma omp barrier
        // Evolution Equation of Concentration field
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

                    fes1 = Ds * (phi[1][i][j][k] * (conp[1][ip][j][k] - conp[1][i][j][k]) / dx);
                    fws1 = Ds * (phi[1][i][j][k] * (conp[1][i][j][k] - conp[1][im][j][k]) / dx);
                    if (phi[1][ip][j][k] == 0.0)
                    {
                        fes1 = 0.0;
                    }
                    if (phi[1][im][j][k] == 0.0)
                    {
                        fws1 = 0.0;
                    }
                    fns1 = Ds * (phi[1][i][j][k] * (conp[1][i][jp][k] - conp[1][i][j][k]) / dx);
                    fss1 = Ds * (phi[1][i][j][k] * (conp[1][i][j][k] - conp[1][i][jm][k]) / dx);
                    if (phi[1][i][jp][k] == 0.0)
                    {
                        fns1 = 0.0;
                    }
                    if (phi[1][i][jm][k] == 0.0)
                    {
                        fss1 = 0.0;
                    }
                    fis1 = Ds * (phi[1][i][j][k] * (conp[1][i][j][kp] - conp[1][i][j][k]) / dx);
                    fos1 = Ds * (phi[1][i][j][k] * (conp[1][i][j][k] - conp[1][i][j][km]) / dx);
                    if (phi[1][i][j][kp] == 0.0)
                    {
                        fis1 = 0.0;
                    }
                    if (phi[1][i][j][km] == 0.0)
                    {
                        fos1 = 0.0;
                    }

                    fes2 = Ds * (phi[2][i][j][k] * (conp[2][ip][j][k] - conp[2][i][j][k]) / dx);
                    fws2 = Ds * (phi[2][i][j][k] * (conp[2][i][j][k] - conp[2][im][j][k]) / dx);
                    if (phi[2][ip][j][k] == 0.0)
                    {
                        fes2 = 0.0;
                    }
                    if (phi[2][im][j][k] == 0.0)
                    {
                        fws2 = 0.0;
                    }
                    fns2 = Ds * (phi[2][i][j][k] * (conp[2][i][jp][k] - conp[2][i][j][k]) / dx);
                    fss2 = Ds * (phi[2][i][j][k] * (conp[2][i][j][k] - conp[2][i][jm][k]) / dx);
                    if (phi[2][i][jp][k] == 0.0)
                    {
                        fns2 = 0.0;
                    }
                    if (phi[2][i][jm][k] == 0.0)
                    {
                        fss2 = 0.0;
                    }
                    fis2 = Ds * (phi[2][i][j][k] * (conp[2][i][j][kp] - conp[2][i][j][k]) / dx);
                    fos2 = Ds * (phi[2][i][j][k] * (conp[2][i][j][k] - conp[2][i][j][km]) / dx);
                    if (phi[2][i][j][kp] == 0.0)
                    {
                        fis2 = 0.0;
                    }
                    if (phi[2][i][j][km] == 0.0)
                    {
                        fos2 = 0.0;
                    }

                    // diffsion in the liquid of phase 1
                    if (phi[1][i][j][k] > 0.0 || phi[2][i][j][k] > 0.0)
                    {
                        lfrac1 = (phi[1][i][j][k] / (phi[1][i][j][k] + phi[2][i][j][k]));
                    }
                    else
                    {
                        lfrac1 = 1.0;
                    }

                    fel1 = Dl * phi[0][i][j][k] * lfrac1 * (conp[0][ip][j][k] - conp[0][i][j][k]) / dx;
                    fwl1 = Dl * phi[0][i][j][k] * lfrac1 * (conp[0][i][j][k] - conp[0][im][j][k]) / dx;
                    if (phi[0][ip][j][k] == 0.0 || (phi[1][ip][j][k] == 0.0 && phi[2][ip][j][k] > 0.0))
                    {
                        fel1 = 0.0;
                    }
                    if (phi[0][im][j][k] == 0.0 || (phi[1][im][j][k] == 0.0 && phi[2][im][j][k] > 0.0))
                    {
                        fwl1 = 0.0;
                    }
                    fnl1 = Dl * phi[0][i][j][k] * lfrac1 * (conp[0][i][jp][k] - conp[0][i][j][k]) / dx;
                    fsl1 = Dl * phi[0][i][j][k] * lfrac1 * (conp[0][i][j][k] - conp[0][i][jm][k]) / dx;
                    if (phi[0][i][jp][k] == 0.0 || (phi[1][i][jp][k] == 0.0 && phi[2][i][jp][k] > 0.0))
                    {
                        fnl1 = 0.0;
                    }
                    if (phi[0][i][jm][k] == 0.0 || (phi[1][i][jm][k] == 0.0 && phi[2][i][jm][k] > 0.0))
                    {
                        fsl1 = 0.0;
                    }
                    fil1 = Dl * phi[0][i][j][k] * lfrac1 * (conp[0][i][j][kp] - conp[0][i][j][k]) / dx;
                    fol1 = Dl * phi[0][i][j][k] * lfrac1 * (conp[0][i][j][k] - conp[0][i][j][km]) / dx;
                    if (phi[0][i][j][kp] == 0.0 || (phi[1][i][j][kp] == 0.0 && phi[2][i][j][kp] > 0.0))
                    {
                        fil1 = 0.0;
                    }
                    if (phi[0][i][j][km] == 0.0 || (phi[1][i][j][km] == 0.0 && phi[2][i][j][km] > 0.0))
                    {
                        fol1 = 0.0;
                    }

                    cddttl1 = (fel1 - fwl1 + fnl1 - fsl1 + fil1 - fol1) / dx;
                    cddtts1 = (fes1 - fws1 + fns1 - fss1 + fis1 - fos1) / dx;
                    cddtts2 = (fes2 - fws2 + fns2 - fss2 + fis2 - fos2) / dx;

                    // liquid phase
                    if (phi[0][i][j][k] > 0.0)
                    {
                        conp2[0][i][j][k] = conp[0][i][j][k] + cddttl1 * dtime / phi[0][i][j][k];
                    }
                    else
                    {
                        conp2[0][i][j][k] = conp[0][i][j][k];
                    }
                    // phase 1
                    if (phi[1][i][j][k] > 0.0)
                    {
                        conp2[1][i][j][k] = conp[1][i][j][k] + cddtts1 * dtime / phi[1][i][j][k];
                    }
                    else
                    {
                        conp2[1][i][j][k] = conp[1][i][j][k];
                    }
                    // phase 2
                    if (phi[2][i][j][k] > 0.0)
                    {
                        conp2[2][i][j][k] = conp[2][i][j][k] + cddtts2 * dtime / phi[2][i][j][k];
                    }
                    else
                    {
                        conp2[2][i][j][k] = conp[2][i][j][k];
                    }
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
                    con[i][j][k] = conp2[1][i][j][k] * phi[1][i][j][k] + conp2[2][i][j][k] * phi[2][i][j][k] + conp2[0][i][j][k] * phi[0][i][j][k];
                    if (con[i][j][k] > 1.0)
                    {
                        con[i][j][k] = 1.0;
                    }
                    if (con[i][j][k] < 0.0)
                    {
                        con[i][j][k] = 0.0;
                    }
                    conp[0][i][j][k] = conp2[0][i][j][k];
                    conp[1][i][j][k] = conp2[1][i][j][k];
                    conp[2][i][j][k] = conp2[2][i][j][k];
                    temp[i][j][k] -= rateT * dtime;
                }
            }
        }

#pragma omp barrier

        if (th_id == 0)
        {
            sumc = 0.0;
            suma = 0.0;
            sumb = 0.0;
            sumil = 0.0;
            maxphi01 = 0.0;
            for (i = 0; i <= ndmx; i++)
            {
                for (j = 0; j <= ndmy; j++)
                {
                    for (k = 0; k <= ndmz; k++)
                    {
                        suma += phi[1][i][j][k];
                        sumb += phi[2][i][j][k];
                        sumc += con[i][j][k];
                    }
                }
            }

            for (i = 0; i <= ndmx; i++)
            {
                if (maxphi01 < phi[1][i][0][0] * (1.0 - phi[1][i][0][0]))
                {
                    maxphi01 = phi[1][i][0][0] * (1.0 - phi[1][i][0][0]);
                    tpx = i;
                }
            }

            tpxp = tpx + 1;
            tpxm = tpx - 1;
            tpdphidt = dphi0[1][tpx][0][0] / dtime;
            tpdphidx = (phi[0][tpxp][0][0] - phi[0][tpxm][0][0]) / 2.0 / dx;
            tpv = -tpdphidt / tpdphidx;

            c0 = sumc / NDX / NDY / NDZ;
            dcm0 = sumc - cm0;
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

double calC01e(double temp)
{
    return (ce + (temp - Te) / ml1);
}

double calC1e(double temp)
{
    return (ce + (temp - Te) / ml1) * kap1;
}

double calC02e(double temp)
{
    return (ce + (temp - Te) / ml2);
}

double calC2e(double temp)
{
    return (ce + (temp - Te) / ml2) * kap2;
}

void datasave(int step)
{
    int is, js, ks;

    // FILE *streamc; //ストリームのポインタ設定
    // char bufferc[30];
    // sprintf(bufferc, "data/con/1d%d.csv", step);
    // streamc = fopen(bufferc, "a"); //書き込む先のファイルを追記方式でオープン

    // for (is = 0; is <= ndmx; is++)
    // {
    //     fprintf(streamc, "%lf   ", con[is][0][0]);
    //     fprintf(streamc, "\n");
    // }
    // fclose(streamc); //ファイルをクローズ

    // FILE *streamp; //ストリームのポインタ設定
    // char bufferp[30];
    // sprintf(bufferp, "data/phi/1d%d.csv", step);
    // streamp = fopen(bufferp, "a"); //書き込む先のファイルを追記方式でオープン

    // for (is = 0; is <= ndmx; is++)
    // {
    //     fprintf(streamp, "%lf   ", conp[0][is][0][0]);
    //     fprintf(streamp, "\n");
    // }
    // fclose(streamp); //ファイルをクローズ

    FILE *streamc; // ストリームのポインタ設定
    char bufferc[30];
    sprintf(bufferc, "data/con/2d%d.csv", step);
    streamc = fopen(bufferc, "a"); // 書き込む先のファイルを追記方式でオープン

    for (is = 0; is <= ndmx; is++)
    {
        for (js = 0; js <= ndmy; js++)
        {
            fprintf(streamc, "%lf   ", con[is][js][0]);
            fprintf(streamc, "\n");
        }
    }
    fclose(streamc); // ファイルをクローズ

    FILE *streamp; // ストリームのポインタ設定
    char bufferp[30];
    sprintf(bufferp, "data/phi/2d%d.csv", step);
    streamp = fopen(bufferp, "a"); // 書き込む先のファイルを追記方式でオープン

    for (is = 0; is <= ndmx; is++)
    {
        for (js = 0; js <= ndmy; js++)
        {
            fprintf(streamp, "%lf   ", conp[0][is][js][0]);
            fprintf(streamp, "\n");
        }
    }
    fclose(streamp); // ファイルをクローズ
}

void fracsave(int step)
{
    // FILE *streamd; //ストリームのポインタ設定
    // char bufferd[30];
    // sprintf(bufferd, "data/frac.csv", step);
    // streamd = fopen(bufferd, "a"); //書き込む先のファイルを追記方式でオープン
    // fprintf(streamd, "%lf   ", (suma - suma0) / NDX / NDY / NDZ);
    // fprintf(streamd, "\n");
    // fclose(streamd); //ファイルをクローズ

    FILE *streamv; // ストリームのポインタ設定
    char bufferv[30];
    sprintf(bufferv, "data/vel.csv", step);
    streamv = fopen(bufferv, "a"); // 書き込む先のファイルを追記方式でオープン
    fprintf(streamv, "%lf   ", tpv);
    fprintf(streamv, "\n");
    fclose(streamv); // ファイルをクローズ
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#define TM 1000.0 // melting point (K)
#define DT 0.1    // undercooling (K)
#define TG 1000.0 // temperature gradient (K/m)
#define TR 0.01   // cooling rate (K/s)
#define VV 1.0e-5 // steady-state velocity (m/s)
#define LL 1.0e-3 // domain size (m)
#define TT 1.0e2  // total time (s)
#define GA 0.5    // surface energy (J/m2)
#define S0 1.0e6  // fusion entropy (J/m3K)
#define HC 100.0  // heat conductivity (W/mK)
#define SH 5.0e6  // specific heat (J/Km3)
#define LH 1.0e9  // latent heat (J/m3)
#define PI 3.141592

int NDX, NDY, NDZ;
int ndmx, ndmy, ndmz;

int NDXT, NDYT, NDZT;
int ndmxt, ndmyt, ndmzt;

int i, j, k, ip, im, jp, jm, kp, km;
double phidxx, phidyy, phidzz;
double delphi, deltemp;

double Tm, temp0;
double gamma0, eta, mu, s0, gt, Dt;

double dx, dtime;
int istep, nstep, pstep;

double dxt, dtimet, calctime;
int id, jd, idp, jdp, is, js, isp, jsp;

double ***phi, ***phi2;
double ***temp, ***temp2;

void datasave1d(int step), datasave2d(int step);

int main()
{
    // determine the grid size
    gt = GA / S0;
    eta = 0.9 * PI * gt / DT;
    dx = eta / 5.0;

    // determine the time step
    mu = VV / DT;
    dtime = 0.2 * dx * dx / (mu * gt);
    Dt = HC / SH;
    // dtime = 0.2 * dx * dx / Dt;

    // calculate the the grid number
    NDX = int(LL / dx / 20) * 20;
    NDY = 1;
    NDZ = 1;

    ndmx = NDX - 1;
    ndmy = NDY - 1;
    ndmz = NDZ - 1;

    // calculate the the step number
    nstep = TT / dtime + 1;
    pstep = (nstep - 1) / 10;

    phi = new double **[NDX];
    phi2 = new double **[NDX];
    temp = new double **[NDX];
    temp2 = new double **[NDX];
    for (i = 0; i <= ndmx; i++)
    {
        phi[i] = new double *[NDY];
        phi2[i] = new double *[NDY];
        temp[i] = new double *[NDY];
        temp2[i] = new double *[NDY];
        for (j = 0; j <= ndmy; j++)
        {
            phi[i][j] = new double[NDZ];
            phi2[i][j] = new double[NDZ];
            temp[i][j] = new double[NDZ];
            temp2[i][j] = new double[NDZ];
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                if (i < NDX / 8)
                {
                    phi[i][j][k] = 0.0;
                }
                else
                {
                    phi[i][j][k] = 1.0;
                }
                temp[i][j][k] = (TM - DT) + (i - NDX / 8) * dx * TG;
            }
        }
    }

    FILE *stream; //ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/domain.csv");
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン
    fprintf(stream, "%d   ", NDX);
    fprintf(stream, "\n");
    fprintf(stream, "%d   ", pstep);
    fprintf(stream, "\n");
    fclose(stream); //ファイルをクローズ

    std::cout << "dx: " << dx << std::endl;
    std::cout << "NDX: " << NDX << std::endl;
    std::cout << "dtime: " << dtime << std::endl;
    std::cout << "nstep: " << nstep << std::endl;

    std::cout << "M: " << mu * gt << std::endl;
    std::cout << "D: " << Dt << std::endl;

start:;

    if (istep % pstep == 0)
    {
        datasave1d(istep);
    }

    for (i = 0; i <= ndmx; i++)
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
                    jp = 0;
                }
                if (j == 0)
                {
                    jm = ndmy;
                }
                if (k == ndmz)
                {
                    kp = 0;
                }
                if (k == 0)
                {
                    km = ndmz;
                }

                phidxx = (phi[ip][j][k] + phi[im][j][k] - 2.0 * phi[i][j][k]) / dx / dx; //フェーズフィールドの空間２階微分
                phidyy = (phi[i][jp][k] + phi[i][jm][k] - 2.0 * phi[i][j][k]) / dx / dx;
                phidzz = (phi[i][j][kp] + phi[i][j][km] - 2.0 * phi[i][j][k]) / dx / dx;

                delphi = mu * (gt * (phidxx + phidyy + phidzz + PI * PI / eta / eta * (phi[i][j][k] - 0.5)) - PI / eta * sqrt(phi[i][j][k] * (1.0 - phi[i][j][k])) * (DT));

                phi2[i][j][k] = phi[i][j][k] + delphi * dtime;

                if (phi2[i][j][k] >= 1.0)
                {
                    phi2[i][j][k] = 1.0;
                }
                if (phi2[i][j][k] <= 0.0)
                {
                    phi2[i][j][k] = 0.0;
                }

                deltemp = HC / SH * ((temp[ip][j][k] + temp[im][j][k] + temp[i][jp][k] + temp[i][jm][k] + temp[i][j][kp] + temp[i][j][km] - 6.0 * temp[i][j][k]) / dx / dx) - 30.0 * phi[i][j][k] * (1.0 - phi[i][j][k]) * phi[i][j][k] * (1.0 - phi[i][j][k]) * LH / SH * delphi;
                temp2[i][j][k] = temp[i][j][k] + deltemp * dtime;
            }
        }
    }

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                phi[i][j][k] = phi2[i][j][k];
                temp2[i][j][k] -= TR * dtime;
                temp[i][j][k] = temp2[i][j][k];
            }
        }
    }

    istep = istep + 1;
    if (istep < nstep)
    {
        goto start;
    }

end:;
}

void datasave1d(int step)
{
    FILE *stream; //ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/1d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndmx; i++)
    {
        fprintf(stream, "%lf   ", phi[i][0][0]);
        fprintf(stream, "\n");
    }
    fclose(stream); //ファイルをクローズ
}

void datasave2d(int step)
{
    FILE *stream; //ストリームのポインタ設定
    char buffer[30];
    sprintf(buffer, "data/2d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            fprintf(stream, "%lf   ", phi[i][j][0]);
            fprintf(stream, "\n");
        }
    }
    fclose(stream); //ファイルをクローズ
}
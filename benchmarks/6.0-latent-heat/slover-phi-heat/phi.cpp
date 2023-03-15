#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#define LL 1.0e-6
#define PI 3.141592

int NDX, NDY, NDZ;
int ndx, ndy, ndz;
int ndmx, ndmy, ndmz;

int i, j, k, ip, im, jp, jm, kp, km;
double phidxx, phidyy, phidzz;
double delphi;
double sumie, sumis, phidx, phidy;

double Tm, temp0;
double gamma0, eta, mu, s0;
double cndct, speht, rlate;

double deltemp;

double dx, dtime, dtimem, dtimet;
int istep, nstep, pstep;

double ***phi, ***phi2;
double ***temp, ***temp2;

void datasave1d(int step), datasave2d(int step);

int main()
{
    cndct = 84.01;     //熱伝導率[W/mK]
    speht = 5.42e+06;  //比熱[J/Km^3]
    rlate = 2.350e+09; //潜熱[J/m^3]

    Tm = 1728.0;
    temp0 = 1511.2;
    s0 = 1.35e6;
    gamma0 = 0.37;
    eta = 0.9 * PI * gamma0 / (Tm - temp0) / s0;
    dx = eta / 5.0;

    mu = 1.0e-5;
    dtimem = 0.2 * dx * dx / mu / gamma0;
    dtimet = 0.2 * dx * dx / cndct * speht;
    if (dtimem > dtimet)
    {
        dtime = dtimet;
    }
    if (dtimem < dtimet)
    {
        dtime = dtimem;
    }

    NDX = int(LL / dx);
    NDY = 1;
    NDZ = 1;

    ndmx = NDX - 1;
    ndmy = NDY - 1;
    ndmz = NDZ - 1;

    nstep = 100001;
    pstep = 10000;

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
                if (i < NDX / 8) // && j < NDY / 8)
                // if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) < NDX/8  * NDX/8)
                {
                    phi[i][j][k] = 0.0;
                }
                else
                {
                    phi[i][j][k] = 1.0;
                }
                temp[i][j][k] = temp0;
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
    std::cout << "dtimem: " << dtimem << std::endl;
    std::cout << "dtimet: " << dtimet << std::endl;

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

                delphi = mu * (gamma0 * (phidxx + phidyy + phidzz + PI * PI / eta / eta * (phi[i][j][k] - 0.5)) - PI / eta * sqrt(phi[i][j][k] * (1.0 - phi[i][j][k])) * s0 * (Tm - temp[i][j][k]));

                phi2[i][j][k] = phi[i][j][k] + delphi * dtime;

                if (phi2[i][j][k] >= 1.0)
                {
                    phi2[i][j][k] = 1.0;
                }
                if (phi2[i][j][k] <= 0.0)
                {
                    phi2[i][j][k] = 0.0;
                }

                deltemp = cndct / speht * ((temp[ip][j][k] + temp[im][j][k] + temp[i][jp][k] + temp[i][jm][k] + temp[i][j][kp] + temp[i][j][km] - 6.0 * temp[i][j][k]) / dx / dx) - 30.0 * phi[i][j][k] * (1.0 - phi[i][j][k]) * phi[i][j][k] * (1.0 - phi[i][j][k]) * rlate / speht * delphi;
                temp2[i][j][k] = temp[i][j][k] + deltemp * dtime; //温度場の時間発展(陽解法)
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
                temp[i][j][k] = temp2[i][j][k];
            }
        }
    }

    sumie = 0.0;
    sumis = 0.0;
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
                if (phi[i][j][k] > 0.0 && phi[i][j][k] < 1.0)
                {
                    sumis += dx * dx;
                }
                phidx = (phi[ip][j][k] - phi[im][j][k]) / 2.0 / dx;
                phidy = (phi[i][jp][k] - phi[i][jm][k]) / 2.0 / dx;
                sumie += ((phidx * phidx + phidy * phidy) * 4.0 * eta * eta / PI / PI + 4.0 * phi[i][j][k] * (1.0 - phi[i][j][k])) * gamma0 / eta * dx * dx;
            }
        }
    }

    sumis /= eta;

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
    sprintf(buffer, "data/phi/1d%d.csv", step);
    stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndmx; i++)
    {
        fprintf(stream, "%lf   ", phi[i][0][0]);
        fprintf(stream, "\n");
    }
    fclose(stream); //ファイルをクローズ

    FILE *streamt; //ストリームのポインタ設定
    char buffert[30];
    sprintf(buffert, "data/temp/1d%d.csv", step);
    streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

    for (i = 0; i <= ndmx; i++)
    {
        fprintf(streamt, "%lf   ", temp[i][0][0]);
        fprintf(streamt, "\n");
    }
    fclose(streamt); //ファイルをクローズ
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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#define PI 3.141592
#define NDX 128
#define NDY 1
#define NDZ 1

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;

double phi[NDX][NDY][NDZ], phi2[NDX][NDY][NDZ];

int i, j, k, ip, im, jp, jm, kp, km;
double phidxx, phidyy, phidzz;
double delphi;
double sumie, sumis, phidx, phidy;

double LL;
double Tm, temp0;
double gamma0, eta, mu, s0;

double dx, dtime;
int istep, nstep, pstep;

void datasave1d(int step), datasave2d(int step);

int main()
{
    LL = 1.28e-3;
    dx = LL / NDX;
    Tm = 1000.0;
    temp0 = 999.99;
    gamma0 = 0.5;
    eta = 8.0 * dx;
    mu = 1.0;
    dtime = 0.2 * dx * dx / mu / gamma0;
    s0 = 1.0e6;

    nstep = 1001;
    pstep = 100;

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                if (i < NDX / 4) // && j < NDY / 8)
                // if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) < NDX / 8 * NDX / 8)
                {
                    phi[i][j][k] = 0.0;
                }

                else
                {
                    phi[i][j][k] = 1.0;
                }
            }
        }
    }

start:;

    if (istep % pstep == 0)
    {
        datasave1d(istep);
        std::cout << "The excessive energy of the diffuse inteface is: " << sumie / sumis << std::endl;
        std::cout << "stability number of reaction-diffusion:" << eta * (Tm - temp0) * s0 / PI / gamma0 << std::endl;
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

                delphi = mu * (gamma0 * (phidxx + phidyy + phidzz + PI * PI / eta / eta * (phi[i][j][k] - 0.5)) - PI / eta * sqrt(phi[i][j][k] * (1.0 - phi[i][j][k])) * s0 * (Tm - temp0));

                phi2[i][j][k] = phi[i][j][k] + delphi * dtime;

                if (phi2[i][j][k] >= 1.0)
                {
                    phi2[i][j][k] = 1.0;
                }
                if (phi2[i][j][k] <= 0.0)
                {
                    phi2[i][j][k] = 0.0;
                }
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
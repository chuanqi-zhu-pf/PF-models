#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#define LL 1.0e-6
#define PI 3.141592
#define NDX 1280
#define NDY 1
#define NDZ 1

int ndmx = NDX - 1;
int ndmy = NDY - 1;
int  ndmz = NDZ - 1;
int ndx, ndy, ndz;

int i, j, k, ip, im, jp, jm, kp, km;
double phidxx, phidyy, phidzz;
double dpptt;
double phidx, phidy;

double Tm, temp0;
double gt0, eta, mu, s0;
double Dt, Lh, Cp;

double dtttt;

double dx, dtime, dtimem, dtimet;
int istep, nstep, pstep;

double phi[NDX][NDY][NDZ], phi2[NDX][NDY][NDZ];
double temp[NDX][NDY][NDZ],temp2[NDX][NDY][NDZ];

void datasave1d(int step), datasave2d(int step);

int main()
{
    Dt = 1.6e-5;  
    Lh = 2.35e9;  //潜熱[J/m^3]
    Cp = 5.42e6;  //比熱[J/Km^3]

    Tm = 1728.0;
    temp0 = 1511.2;

    mu = 13.5;
    gt0 = 2.74e-7;

    eta = 3.5e-09;
    dx = eta / 5.0;
    dtime = 0.2 * dx * dx / Dt;
   
    std::cout << "temp stability" <<  Dt  * dtime / dx / dx << std::endl;
    std::cout << "phi stability" << gt0 * mu * dtime / dx / dx << std::endl;
	std::cout << "dT stability" << eta * (Tm - temp0) / PI / gt0 << std::endl;

    nstep = 100001;
    pstep = 10000;

    for (i = 0; i <= ndmx; i++)
    {
        for (j = 0; j <= ndmy; j++)
        {
            for (k = 0; k <= ndmz; k++)
            {
                if (i < NDX / 8) // && j < NDY / 8)
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

                // change of solid
                dpptt = - mu * ( gt0  * (phidxx + phidyy + phidzz + PI * PI / eta / eta * (phi[i][j][k] - 0.5)) - PI / eta * sqrt(phi[i][j][k] * (1.0 - phi[i][j][k]))  * (Tm - temp[i][j][k]));

                phi2[i][j][k] = phi[i][j][k] - dpptt * dtime;

                if (phi2[i][j][k] >= 1.0)
                {
                    phi2[i][j][k] = 1.0;
                }
                if (phi2[i][j][k] <= 0.0)
                {
                    phi2[i][j][k] = 0.0;
                }

                dtttt = Dt* ((temp[ip][j][k] + temp[im][j][k] + temp[i][jp][k] + temp[i][jm][k] + temp[i][j][kp] + temp[i][j][km] - 6.0 * temp[i][j][k]) / dx / dx) +  8.0 / PI*sqrt(phi[i][j][k] * (1.0 - phi[i][j][k])) * Lh / Cp * dpptt;
                temp2[i][j][k] = temp[i][j][k] + dtttt* dtime; //温度場の時間発展(陽解法)
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
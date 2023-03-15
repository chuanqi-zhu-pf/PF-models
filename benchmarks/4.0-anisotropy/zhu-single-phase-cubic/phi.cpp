#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#define TM 1000.0 // melting point (K)
#define DT 0.05   // undercooling (K)
#define TG 0.0    // temperature gradient (K/m)
#define TR 0.01   // cooling rate (K/s)
#define VV 1.0e-5 // steady-state velocity (m/s)
#define LL 1.0e-3 // domain size (m)
#define TT 1.5e1  // total time (s)
#define GA 0.5    // surface energy (J/m2)
#define S0 1.0e6  // fusion entropy (J/(m3*K))
#define PI 3.141592

int NDX, NDY, NDZ;
int ndx, ndy, ndz;
int ndmx, ndmy, ndmz;

int i, j, k, ip, im, jp, jm, kp, km;
double phidxx, phidyy, phidzz;
double delphi;

double Tm, temp0;
double gamma0, eta, mu, s0, gt;

double dx, dtime;
int istep, nstep, pstep;

double nx, ny, nz;
double astre, af;
double phidx, phidy, phidz, phiabs2, phiabs;
double dphiabs2dx, dphiabs2dy, dphiabs2dz;
double dphiabs2dphix, dphiabs2dphiy, dphiabs2dphiz;
double dphiabs2dphixdx, dphiabs2dphiydy, dphiabs2dphizdz;
double phidxy, phidyz, phidxz;
double phidyx, phidzy, phidzx;
double nxx, nxy, nxz;
double nyx, nyy, nyz;
double nzx, nzy, nzz;
double nxphix, nyphix, nzphix;
double nxphiy, nyphiy, nzphiy;
double nxphiz, nyphiz, nzphiz;
double nxphixdx, nyphixdx, nzphixdx;
double nxphiydy, nyphiydy, nzphiydy;
double nxphizdz, nyphizdz, nzphizdz;

double AA, II;
double IIdx, IIdphix, IIdphixdx;
double IIdy, IIdphiy, IIdphiydy;
double IIdz, IIdphiz, IIdphizdz;

double AAdx, AAdphix, AAdphixdx;
double AAdy, AAdphiy, AAdphiydy;
double AAdz, AAdphiz, AAdphizdz;

double termx, termy, termz;

double ***phi, ***phi2;
double ***temp, ***temp2;

void datasave1d(int step), datasave2d(int step);

int main()
{
    astre = 0.10;
    // determine the grid size and time step
    gt = GA / S0;
    eta = 0.5 * PI * gt / DT;
    dx = eta / 5.0;
    mu = VV / DT;
    dtime = 0.2 * dx * dx / (mu * gt);

    // calculate the the grid number
    NDX = int(LL / dx);
    NDY = int(LL / dx);
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
                // if (i < NDX / 8)
                if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) < NDX / 64 * NDX / 64)
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

start:;

    if (istep % pstep == 0)
    {
        datasave2d(istep);
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

                phidx = (phi[ip][j][k] - phi[im][j][k]) / 2.0 / dx;
                phidy = (phi[i][jp][k] - phi[i][jm][k]) / 2.0 / dx;
                phidz = (phi[i][j][kp] - phi[i][j][km]) / 2.0 / dx;

                phiabs2 = phidx * phidx + phidy * phidy + phidz * phidz;
                phiabs = sqrt(phiabs2);

                phidxx = (phi[ip][j][k] + phi[im][j][k] - 2.0 * phi[i][j][k]) / dx / dx; //フェーズフィールドの空間２階微分
                phidyy = (phi[i][jp][k] + phi[i][jm][k] - 2.0 * phi[i][j][k]) / dx / dx;
                phidzz = (phi[i][j][kp] + phi[i][j][km] - 2.0 * phi[i][j][k]) / dx / dx;

                if (phiabs2 != 0.0)
                {
                    phidxy = (phi[ip][jp][k] + phi[im][jm][k] - phi[im][jp][k] - phi[ip][jm][k]) / 4.0 / dx / dx;
                    phidxz = (phi[ip][j][kp] + phi[im][j][km] - phi[im][j][kp] - phi[ip][j][km]) / 4.0 / dx / dx;
                    phidyz = (phi[i][jp][kp] + phi[i][jm][km] - phi[i][jm][kp] - phi[i][jp][km]) / 4.0 / dx / dx;

                    phidyx = phidxy;
                    phidzx = phidxz;
                    phidzy = phidyz;

                    nx = phidx / phiabs;
                    ny = phidy / phiabs;
                    nz = phidz / phiabs;

                    dphiabs2dx = 2.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz);
                    dphiabs2dy = 2.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz);
                    dphiabs2dz = 2.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz);

                    dphiabs2dphix = 2.0 * phidx;
                    dphiabs2dphiy = 2.0 * phidy;
                    dphiabs2dphiz = 2.0 * phidz;

                    dphiabs2dphixdx = 2.0 * phidxx;
                    dphiabs2dphiydy = 2.0 * phidyy;
                    dphiabs2dphizdz = 2.0 * phidzz;

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

                    AA = 1.0 - 3.0 * astre + 4.0 * astre * (pow(nx, 4.0) + pow(ny, 4.0) + pow(nz, 4.0));
                    II = 4.0 * eta / PI / PI * phiabs2 + 4.0 * phi[i][j][k] * (1.0 - phi[i][j][k]);

                    IIdx = 4.0 * eta / PI / PI * dphiabs2dx + 4.0 * phidx * (1.0 - 2.0 * phi[i][j][k]);
                    IIdy = 4.0 * eta / PI / PI * dphiabs2dy + 4.0 * phidy * (1.0 - 2.0 * phi[i][j][k]);
                    IIdz = 4.0 * eta / PI / PI * dphiabs2dz + 4.0 * phidz * (1.0 - 2.0 * phi[i][j][k]);

                    IIdphix = 4.0 * eta / PI / PI * 2.0 * phidx;
                    IIdphiy = 4.0 * eta / PI / PI * 2.0 * phidy;
                    IIdphiz = 4.0 * eta / PI / PI * 2.0 * phidz;

                    IIdphixdx = 4.0 * eta / PI / PI * 2.0 * phidxx;
                    IIdphiydy = 4.0 * eta / PI / PI * 2.0 * phidyy;
                    IIdphizdz = 4.0 * eta / PI / PI * 2.0 * phidzz;

                    AAdx = 4.0 * astre * (4.0 * pow(nx, 3.0) * nxx + 4.0 * pow(ny, 3.0) * nyx + 4.0 * pow(nz, 3.0) * nzx);
                    AAdy = 4.0 * astre * (4.0 * pow(nx, 3.0) * nxy + 4.0 * pow(ny, 3.0) * nyy + 4.0 * pow(nz, 3.0) * nzy);
                    AAdz = 4.0 * astre * (4.0 * pow(nx, 3.0) * nxz + 4.0 * pow(ny, 3.0) * nyz + 4.0 * pow(nz, 3.0) * nzz);

                    AAdphix = 4.0 * astre * (4.0 * pow(nx, 3.0) * nxphix + 4.0 * pow(ny, 3.0) * nyphix + 4.0 * pow(nz, 3.0) * nzphix);
                    AAdphiy = 4.0 * astre * (4.0 * pow(nx, 3.0) * nxphiy + 4.0 * pow(ny, 3.0) * nyphiy + 4.0 * pow(nz, 3.0) * nzphiy);
                    AAdphiz = 4.0 * astre * (4.0 * pow(nx, 3.0) * nxphiz + 4.0 * pow(ny, 3.0) * nyphiz + 4.0 * pow(nz, 3.0) * nzphiz);

                    AAdphixdx = 16.0 * astre * (3.0 * pow(nx, 2.0) * nxx * nxphix + pow(nx, 3.0) * nxphixdx + 3.0 * pow(ny, 2.0) * nyx * nyphix + pow(ny, 3.0) * nyphixdx + 3.0 * pow(nz, 2.0) * nzx * nzphix + pow(nz, 3.0) * nzphixdx);
                    AAdphiydy = 16.0 * astre * (3.0 * pow(nx, 2.0) * nxy * nxphiy + pow(nx, 3.0) * nxphiydy + 3.0 * pow(ny, 2.0) * nyy * nyphiy + pow(ny, 3.0) * nyphiydy + 3.0 * pow(nz, 2.0) * nzy * nzphiy + pow(nz, 3.0) * nzphiydy);
                    AAdphizdz = 16.0 * astre * (3.0 * pow(nx, 2.0) * nxz * nxphiz + pow(nx, 3.0) * nxphizdz + 3.0 * pow(ny, 2.0) * nyz * nyphiz + pow(ny, 3.0) * nyphizdz + 3.0 * pow(nz, 2.0) * nzz * nzphiz + pow(nz, 3.0) * nzphizdz);

                    termx = (IIdx * AAdphix + II * AAdphixdx + AAdx * IIdphix + AA * IIdphixdx) * PI * PI / 8.0 / eta;
                    termy = (IIdy * AAdphiy + II * AAdphiydy + AAdy * IIdphiy + AA * IIdphiydy) * PI * PI / 8.0 / eta;
                    termz = (IIdz * AAdphiz + II * AAdphizdz + AAdz * IIdphiz + AA * IIdphizdz) * PI * PI / 8.0 / eta;

                    delphi = mu * (gt * (termx + termy + termz + PI * PI / eta / eta * (phi[i][j][k] - 0.5)) - PI / eta * sqrt(phi[i][j][k] * (1.0 - phi[i][j][k])) * (TM - temp[i][j][k]));
                }
                else
                {
                    delphi = mu * (gt * (phidxx + phidyy + phidzz + PI * PI / eta / eta * (phi[i][j][k] - 0.5)) - PI / eta * sqrt(phi[i][j][k] * (1.0 - phi[i][j][k])) * (TM - temp[i][j][k]));
                }

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
                temp[i][j][k] -= TR * dtime;
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
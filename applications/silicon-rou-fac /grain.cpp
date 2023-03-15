// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/11/11

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#define LL 1.4e-3 // domain size (m)
#define TT 1.2e1  // total time (s)

#define TG 900.0  // temperature gradient (K/m)
#define TR 9.0e-3 // cooling rate (K/s)
#define VV 1.0e-5 // steady-state velocity (m/s)
// #define DT 8.3e-2	  // lowest undercooling (K) (real lowest value = 8.3e-5)
#define DT 0.2	  // largest undercooling (K)


#define TM 1000.0
#define GA 0.5	 // surface energy (J/m2)
#define S0 2.4e4 // fusion entropy (J/(m3*K)) (real value = 2.4e6)

#define N 3

int NDX, NDY, NDZ;
int ndmx, ndmy, ndmz;
int nm = N - 1;
double PI = 3.141592;
double RR = 8.3145;

double ****phi, ****phi2;
int ****phiIdx, ***phiNum;
double ***temp, ***mob;

double aij[N][N],wij[N][N], mij[N][N], sij[N][N];
double  thij[N][N][3], axij[N][N][3][3];

int phinum;
int i, j, k, ip, im, jp, jm, kp, km;
int ii, jj, kk;
int n1, n2, n3;

int istep, nstep, pstep;
double dtime, L, dx;
double M0;				 // 粒界の易動度
double W0;				 // ペナルティー項の係数
double A0;				 // 勾配エネルギー係数
double sum1, sum2, sum3; // 各種の和の作業変数
double pddtt;			 // フェーズフィールドの時間変化率

double gt0;	  // 粒界エネルギ密度
double delta; // 粒界幅（差分ブロック数にて表現）
double mu0, Dtl;	  // 粒界の易動度
double staD, staG;

double astre;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidz;
double phidxx, phidyy, phidzz;

int intpos, allL;
double inttemp;
int rn;

double rofr;
double k_ro, k_sm;

double miijj;
double zeta1, ww1, rp0, rp1;
double th;
double nxii, nyii, nzii;
double phidxii, phidyii, phidzii, phiabs2ii;
double min_val;
int min_idx, l;
double a1, a2, a3;
double ux0, uy0, uz0;
double ux, uy, uz, uu;
double al;
#define FN 4
double face[FN][3],ang[N][3];

int fn, fi, idx;
double max;
double sums, sumr, sumg, sumb;
double colii[N][3];
double FaceVec[N][3][3];
double FaceIdx[N][3], FaceAng[N][3];
double xx1, yy1, zz1;
double xx2, yy2, zz2;
double xx3, yy3, zz3;
double aa, bb;
double rr0, gg0, bb0;
double rr00, gg00, bb00;
double face111[4][3], face110[6][3], face100[3][3];

void datasave1d(int step), datasave2d(int step), datasave3d(int step);

//******* メインプログラム ******************************************
int main()
{
	// k_ro = 1.2e-4;
	// k_sm = 0.6e-4;
	mu0 = 1.2e-4;
	gt0 = GA / S0;

	dx = 7.0e-6;
	dtime = 4.0e-3;
	delta = 7.0 * dx;

	nstep = int(TT / dtime) + 1;
	pstep = nstep / 10;

	NDX = int(LL / dx);
	NDY = NDX / 4 * 3;
	NDZ = 1;

	ndmx = NDX - 1;
	ndmy = NDY - 1;
	ndmz = NDZ - 1;

	zeta1 = 0.1;
	ww1 = 90.0 / 10.0;
	rp0 = 0.05;
	rp1 = 0.05;

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

	phi = new double ***[N];
	phi2 = new double ***[N];
	for (ii = 0; ii <= nm; ii++)
	{
		phi[ii] = new double **[NDX];
		phi2[ii] = new double **[NDX];
		for (i = 0; i <= ndmx; i++)
		{
			phi[ii][i] = new double *[NDY];
			phi2[ii][i] = new double *[NDY];
			for (j = 0; j <= ndmy; j++)
			{
				phi[ii][i][j] = new double[NDZ];
				phi2[ii][i][j] = new double[NDZ];
			}
		}
	}

	temp = new double **[NDX];
	mob = new double **[NDX];
	for (i = 0; i <= ndmx; i++)
	{
		temp[i] = new double *[NDY];
		mob[i] = new double *[NDY];
		for (j = 0; j <= ndmy; j++)
		{
			temp[i][j] = new double[NDZ];
			mob[i][j] = new double[NDZ];
		}
	}

	phiNum = new int **[NDX];
	for (i = 0; i <= ndmx; i++)
	{
		phiNum[i] = new int *[NDY];

		for (j = 0; j <= ndmy; j++)
		{
			phiNum[i][j] = new int[NDZ];
		}
	}

	phiIdx = new int ***[N + 1];
	for (ii = 0; ii <= N; ii++)
	{
		phiIdx[ii] = new int **[NDX];
		for (i = 0; i <= ndmx; i++)
		{
			phiIdx[ii][i] = new int *[NDY];
			for (j = 0; j <= ndmy; j++)
			{
				phiIdx[ii][i][j] = new int[NDZ];
			}
		}
	}

	FILE *stream; // ストリームのポインタ設定
	char buffer[30];
	sprintf(buffer, "data/domain.csv");
	stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン
	fprintf(stream, "%d   ", NDX);
	fprintf(stream, "\n");
	fprintf(stream, "%d   ", pstep);
	fprintf(stream, "\n");
	fclose(stream); // ファイルをクローズ

	std::cout << "kinetics of rough surface: " << k_ro << std::endl;
	std::cout << "dtime: " << dtime << std::endl;
	std::cout << "nstep: " << nstep << std::endl;
	std::cout << "diffusion stability for phi: " << dtime / dx / dx * k_ro * gt0 << std::endl;

	A0 = 8.0 * delta * gt0 / PI / PI;	// 勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gt0 / delta;				// ペナルティー項の係数[式(4.40)]
	M0 = mu0 * PI * PI / (8.0 * delta); // 粒界の易動度[式(4.40)]

	for (ii = 0; ii <= nm; ii++)
	{
		for (jj = 0; jj <= nm; jj++)
		{
			wij[ii][jj] = W0;
			aij[ii][jj] = A0;
			mij[ii][jj] = M0;
			sij[ii][jj] = 0.0;
			if ((ii == 0) || (jj == 0))
			{
				sij[ii][jj] = 1.0;
			}
			if (ii > 0 && jj > 0)
            {
                mij[ii][jj] = mij[ii][jj] * 0.05;
            }
			if (ii < jj)
			{
				sij[ii][jj] = -sij[ii][jj];
			}
			if (ii == jj)
			{
				wij[ii][jj] = 0.0;
				aij[ii][jj] = 0.0;
				mij[ii][jj] = 0.0;
				sij[ii][jj] = 0.0;
			}
		}
	}

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				// if ((i - NDX / 4) * (i - NDX / 4) + (j - NDY / 2) * (j - NDY / 2) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 8 * NDX / 8)
				if (i < NDX / 32 && j < NDY /2)
				{
					phi[1][i][j][k] = 1.0;
					phi[2][i][j][k] = 0.0;
					phi[0][i][j][k] = 0.0;
				}
				// else if ((i - NDX *3/ 4) * (i - NDX *3/ 4) + (j - NDY / 2) * (j - NDY / 2) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 8 * NDX / 8)
				else if (i < NDX / 32 && j >= NDY/2)
				{
					phi[1][i][j][k] = 0.0;
					phi[2][i][j][k] = 1.0;
					phi[0][i][j][k] = 0.0;
				}
				else
				{
					phi[1][i][j][k] = 0.0;
					phi[2][i][j][k] = 0.0;
					phi[0][i][j][k] = 1.0;
				}
				temp[i][j][k] = (TM - DT) + (i - NDX / 32) * dx * TG;
			}
		}
	}

start:;

	if (istep % pstep == 0)
	{
		datasave3d(istep);
		printf("----------------------------------\n");
        printf("%lf secs (%d steps) has passed!\n", istep * dtime, istep);
        printf("Interface position: %d\n", intpos);
        printf("Interface temperature: %f K\n", inttemp);
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
				mob[i][j][k] = 0.0;
			}
		}
	}

	// Evolution Equations
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
							phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 / dx;
							phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 / dx;
							phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 / dx;

							phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx; // フェーズフィールドの空間２階微分
							phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
							phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

							termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
							termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);

							sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
						}

						// miijj = mij[ii][jj];
						// phidxii = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0;
						// phidyii = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0;
						// phidzii = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0;
						// phiabs2ii = phidxii * phidxii + phidyii * phidyii + phidzii * phidzii;
						// if (phiabs2ii != 0.0)
						// {

						// 	nxii = phidxii / sqrt(phiabs2ii);
						// 	nyii = phidyii / sqrt(phiabs2ii);
						// 	nzii = phidzii / sqrt(phiabs2ii);

						// 	min_val = 0.0;
						// 	min_idx = 0;
						// 	for (l = 0; l < FN; l++)
						// 	{
						// 		ux0 = face[l][0];
						// 		uy0 = face[l][1];
						// 		uz0 = face[l][2];

						// 		for (rn = 0; rn <= 2; rn++)
                        //             {
                        //                 th = thij[ii][jj][rn];
                        //                 a1 = axij[ii][jj][rn][0];
                        //                 a2 = axij[ii][jj][rn][1];
                        //                 a3 = axij[ii][jj][rn][2];

                        //                 ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                        //                 uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
                        //                 uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

                        //                 ux0 = ux;
                        //                 uy0 = uy;
                        //                 uz0 = uz;
                        //             }

						// 		uu = ux * ux + uy * uy + uz * uz;
						// 		al = acos(fabs(nxii * ux + nyii * uy + nzii * uz) / sqrt(uu));

						// 		if (l == 0)
						// 		{
						// 			min_val = al;
						// 			min_idx = l;
						// 		}
						// 		else
						// 		{
						// 			if (min_val > al)
						// 			{
						// 				min_val = al;
						// 				min_idx = l;
						// 			}
						// 		}
						// 	}
						// 	if (min_val < 90.0 / ww1 / 180.0 * PI)
						// 	{
						// 		rofr = zeta1 + (1 - zeta1) * sqrt(pow(tan(ww1 * min_val), 2.0) + rp1 * rp1) * tanh(1.0 / sqrt(pow(tan(ww1 * min_val), 2.0) + rp1 * rp1));
						// 	    mob[i][j][k] = rofr;
						// 	}
						// 	else
						// 	{
						// 		rofr = 1.0;
						// 		mob[i][j][k] = rofr;
						// 	}
						// 	// miijj = (k_sm*(TM-temp[i][j][k])* (1-rofr) + k_ro*rofr)*mij[ii][jj];     
						// 	if (min_val < 90.0 / ww1 / 180.0 * PI)
						// 	{
						// 		miijj = mij[ii][jj] * (zeta1 + (1 - zeta1) * sqrt(pow(tan(ww1 * min_val), 2.0) + rp1 * rp1) * tanh(1.0 / sqrt(pow(tan(ww1 * min_val), 2.0) + rp1 * rp1)));
						// 	}
						// 	else
						// 	{
						// 		miijj = mij[ii][jj];
						// 	} 
							//  // add interface mobility to the interface field
							// if (phi[0][i][j][k] > 0.2 && phi[0][i][j][k] < 0.8)
							// // if (((*phi)[0][i][j][k] * (*phi)[1][i][j][k] > 0.0 || (*phi)[0][i][j][k] * (*phi)[2][i][j][k] > 0.0) && (*phi)[1][i][j][k] * (*phi)[2][i][j][k] == 0.0)
							// {
							// 	mob[i][j][k] = rofr;
							// } 
						// }

					        phidxii = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0;
						    phidyii = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0;
						    phidzii = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0;
                            phiabs2ii = phidxii * phidxii + phidyii * phidyii + phidzii * phidzii;
                            if ( phiabs2ii != 0.0)
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
                                if (phi[0][i][j][k] > 0.2 && phi[0][i][j][k] < 0.8 && phiNum[i][j][k] == 2)
                                // if (((*phi)[0][i][j][k] * (*phi)[1][i][j][k] > 0.0 || (*phi)[0][i][j][k] * (*phi)[2][i][j][k] > 0.0) && (*phi)[1][i][j][k] * (*phi)[2][i][j][k] == 0.0)
                                {
                                    mob[i][j][k] = miijj/mij[ii][jj];
                                }
                            }
                            else
                            {
                                miijj = mij[ii][jj];
                            }
						pddtt += -2.0 * miijj / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * sij[ii][jj] * (TM - temp[i][j][k]) * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
					}
					phi2[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime; // フェーズフィールドの時間発展（陽解法）
					if (phi2[ii][i][j][k] >= 1.0)
					{
						phi2[ii][i][j][k] = 1.0;
					} // フェーズフィールドの変域補正
					if (phi2[ii][i][j][k] <= 0.0)
					{
						phi2[ii][i][j][k] = 0.0;
					}
				}
			} // j
		}	  // i
	}

	for (ii = 0; ii <= nm; ii++)
	{
		for (i = 0; i <= ndmx; i++)
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

	for (i = 0; i <= ndmx; i++)
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
				temp[i][j][k] -= TR * dtime;
			}
		}
	}

	allL = 1;
	// search interface front
	for (i = ndmx; i <= ndmx; i--)
	{
		if (allL == 0)
		{
			intpos = i;
			inttemp = temp[intpos][0][0];
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

	istep = istep + 1;
	if (istep < nstep)
	{
		goto start;
	}

end:;
	return 0;
}

void datasave1d(int step)
{
	FILE *stream; // ストリームのポインタ設定
	char buffer[30];
	sprintf(buffer, "data/1d%d.csv", step);
	stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		fprintf(stream, "%lf   ", phi[0][i][0][0]);
		fprintf(stream, "\n");
	}
	fclose(stream); // ファイルをクローズ
}

void datasave2d(int step)
{
	FILE *stream; // ストリームのポインタ設定
	char buffer[30];
	sprintf(buffer, "data/2d%d.csv", step);
	stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(stream, "%lf   ", mob[i][j][0]);
			fprintf(stream, "\n");
		}
	}
	fclose(stream); // ファイルをクローズ

	
}

void datasave3d(int step)
{

FILE *streami;
char bufferi[30];
sprintf(bufferi, "data/mob/3d%d.vtk", step);
streami = fopen(bufferi, "a");

fprintf(streami, "# vtk DataFile Version 1.0\n");
fprintf(streami, "phi_%d.vtk\n", step);
fprintf(streami, "ASCII\n");
fprintf(streami, "DATASET STRUCTURED_POINTS\n");
fprintf(streami, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
fprintf(streami, "ORIGIN 0.0 0.0 0.0\n");
fprintf(streami, "ASPECT_RATIO 1.0 1.0 1.0\n");
fprintf(streami, "\n");
fprintf(streami, "POINT_DATA %d\n", NDX * NDY * NDZ);
fprintf(streami, "SCALARS scalars float\n");
fprintf(streami, "LOOKUP_TABLE default\n");

for (k = 0; k <= ndmz; k++)
{
	for (j = 0; j <= ndmy; j++)
	{
		for (i = 0; i <= ndmx; i++)
		{
			fprintf(streami, "%e\n", mob[i][j][k]);
		}
	}
}
fclose(streami);

//  for (ii = 1; ii <= nm; ii++)
//         {
//             max = 0.0;
//             idx = 0;
//             for (fn = 0; fn < 3; fn++)
//             {
//                 ux0 = face100[fn][0];
//                 uy0 = face100[fn][1];
//                 uz0 = face100[fn][2];

//                 for (rn = 0; rn <= 2; rn++)
//                 {
//                     th = thij[ii][0][rn];
//                     a1 = axij[ii][0][rn][0];
//                     a2 = axij[ii][0][rn][1];
//                     a3 = axij[ii][0][rn][2];

//                     ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                     uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                     uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

//                     ux0 = ux;
//                     uy0 = uy;
//                     uz0 = uz;
//                 }

//                 if (max < fabs(ux))
//                 {
//                     max = fabs(ux);
//                     idx = fn;
//                 }
//             }
//             FaceAng[ii][0] = acos(max);
//             FaceIdx[ii][0] = idx;

//             max = 0.0;
//             idx = 0;
//             for (fn = 0; fn < 6; fn++)
//             {
//                 ux0 = face110[fn][0] / sqrt(2.0);
//                 uy0 = face110[fn][1] / sqrt(2.0);
//                 uz0 = face110[fn][2] / sqrt(2.0);

//                 for (rn = 0; rn <= 2; rn++)
//                 {
//                     th = thij[ii][0][rn];
//                     a1 = axij[ii][0][rn][0];
//                     a2 = axij[ii][0][rn][1];
//                     a3 = axij[ii][0][rn][2];

//                     ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                     uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                     uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

//                     ux0 = ux;
//                     uy0 = uy;
//                     uz0 = uz;
//                 }

//                 if (max < fabs(ux))
//                 {
//                     max = fabs(ux);
//                     idx = fn;
//                 }
//             }
//             FaceAng[ii][1] = acos(max);
//             FaceIdx[ii][1] = idx;

//             max = 0.0;
//             idx = 0;
//             for (fn = 0; fn < 4; fn++)
//             {
//                 ux0 = face111[fn][0] / sqrt(3.0);
//                 uy0 = face111[fn][1] / sqrt(3.0);
//                 uz0 = face111[fn][2] / sqrt(3.0);

//                 for (rn = 0; rn <= 2; rn++)
//                 {
//                     th = thij[ii][0][rn];
//                     a1 = axij[ii][0][rn][0];
//                     a2 = axij[ii][0][rn][1];
//                     a3 = axij[ii][0][rn][2];

//                     ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                     uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                     uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

//                     ux0 = ux;
//                     uy0 = uy;
//                     uz0 = uz;
//                 }

//                 if (max < fabs(ux))
//                 {
//                     max = fabs(ux);
//                     idx = fn;
//                 }
//             }
//             FaceAng[ii][2] = acos(max);
//             FaceIdx[ii][2] = idx;
//         }

//         for (ii = 1; ii <= nm; ii++)
//         {
//             fn = FaceIdx[ii][0];

//             ux0 = face100[fn][0];
//             uy0 = face100[fn][1];
//             uz0 = face100[fn][2];

//             for (rn = 0; rn <= 2; rn++)
//             {
//                 th = thij[ii][0][rn];
//                 a1 = axij[ii][0][rn][0];
//                 a2 = axij[ii][0][rn][1];
//                 a3 = axij[ii][0][rn][2];

//                 ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                 uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                 uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

//                 ux0 = ux;
//                 uy0 = uy;
//                 uz0 = uz;
//             }
//             FaceVec[ii][0][0] = ux;
//             FaceVec[ii][0][1] = uy;
//             FaceVec[ii][0][2] = uz;

//             fn = FaceIdx[ii][1];

//             ux0 = face110[fn][0] / sqrt(2.0);
//             uy0 = face110[fn][1] / sqrt(2.0);
//             uz0 = face110[fn][2] / sqrt(2.0);

//             for (rn = 0; rn <= 2; rn++)
//             {
//                 th = thij[ii][0][rn];
//                 a1 = axij[ii][0][rn][0];
//                 a2 = axij[ii][0][rn][1];
//                 a3 = axij[ii][0][rn][2];

//                 ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                 uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                 uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

//                 ux0 = ux;
//                 uy0 = uy;
//                 uz0 = uz;
//             }
//             FaceVec[ii][1][0] = ux;
//             FaceVec[ii][1][1] = uy;
//             FaceVec[ii][1][2] = uz;

//             fn = FaceIdx[ii][2];

//             ux0 = face111[fn][0] / sqrt(3.0);
//             uy0 = face111[fn][1] / sqrt(3.0);
//             uz0 = face111[fn][2] / sqrt(3.0);

//             for (rn = 0; rn <= 2; rn++)
//             {
//                 th = thij[ii][0][rn];
//                 a1 = axij[ii][0][rn][0];
//                 a2 = axij[ii][0][rn][1];
//                 a3 = axij[ii][0][rn][2];

//                 ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                 uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
//                 uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

//                 ux0 = ux;
//                 uy0 = uy;
//                 uz0 = uz;
//             }
//             FaceVec[ii][2][0] = ux;
//             FaceVec[ii][2][1] = uy;
//             FaceVec[ii][2][2] = uz;
//         }

//         for (ii = 1; ii <= nm; ii++)
//         {
//             // 001
//             xx1 = FaceVec[ii][0][0];
//             yy1 = FaceVec[ii][0][1] + 1.0e-10;
//             zz1 = FaceVec[ii][0][2] + 1.0e-10;

//             // 011
//             xx2 = FaceVec[ii][1][0];
//             yy2 = FaceVec[ii][1][1] + 1.0e-10;
//             zz2 = FaceVec[ii][1][2] + 1.0e-10;

//             // 111
//             xx3 = FaceVec[ii][2][0];
//             yy3 = FaceVec[ii][2][1];
//             zz3 = FaceVec[ii][2][2];

//             aa = (yy1 * zz3 - yy3 * zz1) / (yy2 * zz1 - yy1 * zz2);
//             bb = -(aa * zz2 + zz3) / zz1;

//             rr0 = bb / (bb * xx1 + aa * xx2 + xx3);
//             gg0 = aa / (bb * xx1 + aa * xx2 + xx3);
//             bb0 = 1.0 / (bb * xx1 + aa * xx2 + xx3);

//             rr00 = rr0 / sqrt(rr0 * rr0 + gg0 * gg0 + bb0 * bb0);
//             gg00 = gg0 / sqrt(rr0 * rr0 + gg0 * gg0 + bb0 * bb0);
//             bb00 = bb0 / sqrt(rr0 * rr0 + gg0 * gg0 + bb0 * bb0);

//             colii[ii][0] = fabs(rr00);
//             colii[ii][1] = fabs(gg00);
//             colii[ii][2] = fabs(bb00);
//         }

//         FILE *stream;
//         char buffer[30];
//         sprintf(buffer, "data/phi/3d%d.vtk", step);
//         stream = fopen(buffer, "a");

//         fprintf(stream, "# vtk DataFile Version 3.0\n");
//         fprintf(stream, "phi_%d.vtk\n", step);
//         fprintf(stream, "ASCII\n");
//         fprintf(stream, "DATASET STRUCTURED_POINTS\n");
//         fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
//         fprintf(stream, "ORIGIN 0 0 0\n");
//         fprintf(stream, "SPACING 1 1 1\n");
//         fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
//         fprintf(stream, "SCALARS solid float 1\n");
//         fprintf(stream, "LOOKUP_TABLE default\n");

//         for (k = 0; k <= ndmz; k++)
//         {
//             for (j = 0; j <= ndmy; j++)
//             {
//                 for (i = 0; i <= ndmx; i++)
//                 {
//                     fprintf(stream, "%e\n", phi[0][i][j][k]);
//                 }
//             }
//         }

//         fprintf(stream, "VECTORS vectors float\n");

//         for (k = 0; k <= ndmz; k++)
//         {
//             for (j = 0; j <= ndmy; j++)
//             {
//                 for (i = 0; i <= ndmx; i++)
//                 {
//                     sumr = 0.0;
//                     sumg = 0.0;
//                     sumb = 0.0;
//                     for (ii = 1; ii <= nm; ii++)
//                     {
//                         if (phi[ii][i][j][k] > 0.0)
//                         {
//                             sumr += colii[ii][0] * pow(phi[ii][i][j][k], 1.3);
//                             sumg += colii[ii][1] * pow(phi[ii][i][j][k], 1.3);
//                             sumb += colii[ii][2] * pow(phi[ii][i][j][k], 1.3);
//                         }
//                     }
//                     fprintf(stream, "%e %e %e\n", sumr, sumg, sumb);
//                 }
//             }
//         }

//         fclose(stream);
}

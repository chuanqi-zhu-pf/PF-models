// Adapted from Prof. Koyama's MPF code in his textbook
// phase field coupled with concentration field
// Author: Chuanqi Zhu
// Created on: 2022/12/27

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <omp.h>

#define NDX 256
#define NDY 256
#define NDZ 1
#define NTH 8

#define N 2

// public varibales

int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int nm = N - 1;
int rows = NDX / NTH;
double PI = 3.141592;
double RR = 8.3145;

double phi[N][NDX][NDY][NDZ], phi2[N][NDX][NDY][NDZ];
int phiIdx[N + 1][NDX][NDY][NDZ], phiNum[NDX][NDY][NDZ];
double conp[N][NDX][NDY][NDZ], conp2[N][NDX][NDY][NDZ];
double con[NDX][NDY][NDZ], coni[NDX][NDY][NDZ];
double aij[N][N], wij[N][N], mij[N][N], fij[N][N];
double anij[N][N], thij[N][N], vpij[N][N], etaij[N][N];

int nstep, pstep;
double dtime, dx;

double temp0, gamma0, S0, Tm, ml, delta, mobi;
double M0, W0, A0;

int ni, nj, ix, iy, iz;
double sumc, c0, dc0;
double cl, dT0, dTi, dd0, ll0, ml1, kap;
double c1e, c10e;
double Dl, Ds;
double astre, astrem;

double velt;
int istepp, intpos, intposp, intdis;
int allL;

void datasave1d(int step), datasave2d(int step), datasave3d(int step);
double calKv(double veli);

// double calCs(double conp1, double conp0, double phi1, double 0);

//******* メインプログラム ******************************************
int main()
{
	dx = 1.0e-7;
	delta = 4.0 * dx;

	nstep = 80001;
	pstep = 4000;

	mobi = 2.0e-10;
	gamma0 = 0.245;
	S0 = 1.475e6;

	temp0 = 922.5;
	cl = 0.01;

	Tm = 933.47;
	ml = -605.3;
	kap = 0.105;
	Dl = 1.34e-7 * exp(-3.0e4 / RR / temp0);
	Ds = 0.89e-4 * exp(-1.36e5 / RR / temp0);
	dtime = 0.2 * dx * dx / Dl;

	c10e = -(Tm - temp0) / ml;
	c1e = c10e * kap;

	A0 = 8.0 * delta * gamma0 / PI / PI; // 勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gamma0 / delta;			 // ペナルティー項の係数[式(4.40)]
	M0 = mobi * PI * PI / (8.0 * delta); // 粒界の易動度[式(4.40)]

	for (ni = 0; ni <= nm; ni++)
	{
		for (nj = 0; nj <= nm; nj++)
		{
			wij[ni][nj] = W0;
			aij[ni][nj] = A0;
			mij[ni][nj] = M0;
			if (ni == nj)
			{
				wij[ni][nj] = 0.0;
				aij[ni][nj] = 0.0;
				mij[ni][nj] = 0.0;
			}
		}
	}

	//----------- variables for gradeint term with cubic anisotropy ------
	astre = 0.0192;
	astrem = 0.1;
	for (ni = 0; ni <= nm; ni++)
	{
		for (nj = 0; nj <= nm; nj++)
		{
			anij[ni][nj] = 0.0;
			thij[ni][nj] = 0.0;
			vpij[ni][nj] = 0.0;
			etaij[ni][nj] = 0.0;
			if ((ni == 0) || (nj == 0))
			{
				anij[ni][nj] = 1.0;
			}
			if (ni == nj)
			{
				anij[ni][nj] = 0.0;
			}
		}
	}

	// Initialization of phase and concentration fields
	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			for (iz = 0; iz <= ndmz; iz++)
			{
				// if (i < NDX / 8)
				// if ((ix) * (ix) + (iy) * (iy) < NDX / 32 * NDX / 32)
				if (ix + iy < NDX / 16)
				{
					phi[1][ix][iy][iz] = 1.0;
					conp[1][ix][iy][iz] = c1e;
					phi[0][ix][iy][iz] = 0.0;
					conp[0][ix][iy][iz] = c10e;
				}
				else
				{
					phi[1][ix][iy][iz] = 0.0;
					conp[1][ix][iy][iz] = c1e;
					phi[0][ix][iy][iz] = 1.0;
					conp[0][ix][iy][iz] = cl;
				}
				con[ix][iy][iz] = conp[0][ix][iy][iz] * phi[0][ix][iy][iz] + conp[1][ix][iy][iz] * phi[1][ix][iy][iz];
			}
		}
	}

	std::cout << "W/d0: " << delta / dd0 << std::endl;
	std::cout << "diffusion stability for phi: " << dtime * mobi * gamma0 / dx / dx << std::endl;
	std::cout << "diffusion stability for con: " << dtime / dx / dx * Dl << std::endl;

#pragma omp parallel num_threads(NTH)
	{
		// private varibales
		int th_id, istep;
		int start, end, offset;

		th_id = omp_get_thread_num();
		istep = 0;
		offset = th_id * rows;
		start = offset;
		end = offset + rows - 1;

		int phinum;
		int i, j, k, ip, im, jp, jm, kp, km;
		int ii, jj, kk;
		int n1, n2, n3;

		double dG;
		double kvi, veli;

		double flxci;
		double phidxi, phidyi;
		double nxi, nyi;
		int di, xdi, xdip, xdim, ydi, ydip, ydim;
		int ixdi, ixdip, jydi, jydip;

		double sum1, pddtt;
		double wwiikk, wwjjkk;

		double epsilon0;
		double termiikk, termjjkk;

		double phidx, phidy, phidz;
		double phidxx, phidyy, phidzz;
		double phiabs2;
		double phidxy, phidyz, phidxz;
		double phidyx, phidzy, phidzx;

		double th, vp, eta;
		double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;
		double phidxp, phidyp, phidzp;
		double phidxpx, phidypx, phidzpx;
		double phidxpy, phidypy, phidzpy;
		double phidxpz, phidypz, phidzpz;
		double ep, epdx, epdy, epdz;
		double term0;
		double termx, termx0, termx1, termx0dx, termx1dx;
		double termy, termy0, termy1, termy0dy, termy1dy;
		double termz, termz0, termz1, termz0dz, termz1dz;
		double miijj, phidxm, phidym, phidzm, phidxpm, phidypm, phidzpm, phiabsm2;

		double fes1, fns1, fos1, fws1, fss1, fis1;
		double fel, fnl, fol, fwl, fsl, fil;
		double cddttl, cddtts1;
		double conp0ip, conp0im, conp0jp, conp0jm;

	start:;

		if (th_id == 0)
		{
			if (istep % pstep == 0)
			{
				datasave2d(istep);
				std::cout << "the average concentration is: " << sumc / NDX / NDY << std::endl;
				std::cout << "time passed: " << dtime * istep << " s" << std::endl;
			}
			// calculate the average concnetraion
			sumc = 0.0;
			for (i = 0; i <= ndmx; i++)
			{
				for (j = 0; j <= ndmy; j++)
				{
					for (k = 0; k <= ndmz; k++)
					{
						sumc += con[i][j][k];
					}
				}
			}

			// search the tip position
			allL = 1;
			for (i = ndmx; i <= ndmx; i--)
			{
				if (allL == 0)
				{
					intpos = i;
					break;
				}
				for (j = 0; j <= ndmy; j++)
				{
					for (k = 0; k <= ndmz; k++)
					{
						if (phi[0][i][j][k] < 1.0)
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

			if (intposp == 0)
			{
				intposp = intpos;
				istepp = istep;
			}
			else
			{
				if ((intpos - intposp) > 1)
				{
					FILE *stream; // ストリームのポインタ設定
					char buffer[30];
					sprintf(buffer, "data/tip_pos.csv", istep);
					stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン
					fprintf(stream, "%e   ", istep * dtime);
					fprintf(stream, "%e   ", (intpos - intposp) * dx);
					fprintf(stream, "%e   ", (istep - istepp) * dtime);
					fprintf(stream, "\n");
					fclose(stream);
					intposp = intpos;
					istepp = istep;
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
				}
			}
		}

#pragma omp barrier
		// Evolution Equations
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
								phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 / dx;
								phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 / dx;
								phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 / dx;

								phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx; // フェーズフィールドの空間２階微分
								phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
								phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

								//------------------------------ Implementation of gradeint term with cubic anisotropy ---------------------

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

									wwiikk = wij[ii][kk];
									termiikk = term0 + termx + termy + termz;
								}
								else
								{
									wwiikk = wij[ii][kk];
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

									wwjjkk = wij[jj][kk];
									termjjkk = term0 + termx + termy + termz;
								}
								else
								{
									wwjjkk = wij[jj][kk];
									termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
								}

								sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
							}

							if (ii == 1 && jj == 0)
							{
								dG = (conp[0][i][j][k] * ml + Tm - temp0) * S0;
							}
							else if (ii == 0 && jj == 1)
							{
								dG = -(conp[0][i][j][k] * ml + Tm - temp0) * S0;
							}
							miijj = mij[ii][jj];
							if ((ii > 0 && jj == 0) || (ii == 0 && jj > 0))
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

								phidxm = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0;
								phidym = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0;
								phidzm = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0;

								phidxpm = phidxm * xxp + phidym * yxp + phidzm * zxp;
								phidypm = phidxm * xyp + phidym * yyp + phidzm * zyp;
								phidzpm = phidxm * xzp + phidym * yzp + phidzm * zzp;

								phiabsm2 = phidxpm * phidxpm + phidypm * phidypm + phidzpm * phidzpm;

								miijj = miijj * (1.0 - 3.0 * astrem + 4.0 * astrem * (pow(phidxpm, 4.0) + pow(phidypm, 4.0) + pow(phidzpm, 4.0)) / pow(phiabsm2, 2.0));
							}
							pddtt += -2.0 * miijj / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * dG * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
							// フェーズフィールドの発展方程式[式(4.31)]
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

#pragma omp barrier

		for (i = start; i <= end; i++)
		{
			for (j = 0; j <= ndmy; j++)
			{
				for (k = 0; k <= ndmz; k++)
				{
					sum1 = 0.0;
					for (ii = 0; ii <= nm; ii++)
					{
						sum1 += phi2[ii][i][j][k];
					}
					for (ii = 0; ii <= nm; ii++)
					{
						phi2[ii][i][j][k] = phi2[ii][i][j][k] / sum1;
					}
				}
			}
		}

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
					if (phi[0][i][j][k] > 0.0 && phi2[0][i][j][k] < 1.0)
					{

						conp[1][i][j][k] = (conp[1][i][j][k] * phi[1][i][j][k] + (phi2[1][i][j][k] - phi[1][i][j][k]) * (kap * conp[0][i][j][k])) / phi2[1][i][j][k];
						flxci = (phi2[1][i][j][k] - phi[1][i][j][k]) * conp[0][i][j][k] * (1.0 - kap);

						if (phi2[0][i][j][k] == 0.0)
						{
							phidxi = -(phi[1][ip][j][k] - phi[1][im][j][k]) / 2.0;
							phidyi = -(phi[1][i][jp][k] - phi[1][i][jm][k]) / 2.0;
						}
						else
						{
							phidxi = -(phi2[1][ip][j][k] - phi2[1][im][j][k]) / 2.0;
							phidyi = -(phi2[1][i][jp][k] - phi2[1][i][jm][k]) / 2.0;
						}

						nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi);
						nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi);

						di = 0;
						do
						{
							di++;
							xdi = int(round(nxi * di));
							ydi = int(round(nyi * di));

							ixdi = i + xdi;
							jydi = j + ydi;

							if (phi2[0][ixdi][jydi][k] == 1.0)
							{
								con[i][j][k] -= flxci;
								con[ixdi][jydi][k] += flxci;
								break;
							}
						} while (phi2[0][ixdi][jydi][k] < 1.0);
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

		// redistribution of solute to each concentraton field
		for (i = start; i <= end; i++)
		{
			for (j = 0; j <= ndmy; j++)
			{
				for (k = 0; k <= ndmz; k++)
				{
					if (phi[0][i][j][k] == 0.0)
					{
						conp[0][i][j][k] = c10e;
					}
					else if (phi[0][i][j][k] == 1.0)
					{
						conp[1][i][j][k] = c1e;
						conp[0][i][j][k] = con[i][j][k];
					}
				}
			}
		}

#pragma omp barrier
		// Diffusion equation
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
					fos1 = Ds * (phi[1][i][j][k] * (conp[1][i][j][kp] - conp[1][i][j][k]) / dx);
					fis1 = Ds * (phi[1][i][j][k] * (conp[1][i][j][k] - conp[1][i][j][km]) / dx);
					if (phi[1][i][j][kp] == 0.0)
					{
						fos1 = 0.0;
					}
					if (phi[1][i][j][km] == 0.0)
					{
						fis1 = 0.0;
					}

					fel = Dl * phi[0][i][j][k] * (conp[0][ip][j][k] - conp[0][i][j][k]) / dx;
					fwl = Dl * phi[0][i][j][k] * (conp[0][i][j][k] - conp[0][im][j][k]) / dx;
					// back obstacle
					if (phi[0][ip][j][k] == 0.0 && phi[0][i][j][k] > 0.0)
					{
						fel = 0.0;
					}
					if (phi[0][im][j][k] == 0.0 && phi[0][i][j][k] > 0.0)
					{
						fwl = 0.0;
					}
					// front obstacle
					if (phi[0][ip][j][k] == 1.0 && phi[0][i][j][k] < 1.0)
					{
						fel = 0.0;
					}
					if (phi[0][im][j][k] == 1.0 && phi[0][i][j][k] < 1.0)
					{
						fwl = 0.0;
					}
					if (phi[0][ip][j][k] < 1.0 && phi[0][i][j][k] == 1.0)
					{
						fel = 0.0;
					}
					if (phi[0][im][j][k] < 1.0 && phi[0][i][j][k] == 1.0)
					{
						fwl = 0.0;
					}

					fnl = Dl * phi[0][i][j][k] * (conp[0][i][jp][k] - conp[0][i][j][k]) / dx;
					fsl = Dl * phi[0][i][j][k] * (conp[0][i][j][k] - conp[0][i][jm][k]) / dx;
					// back obstacle
					if (phi[0][i][jp][k] == 0.0)
					{
						fnl = 0.0;
					}
					if (phi[0][i][jm][k] == 0.0)
					{
						fsl = 0.0;
					}
					// front obstacle
					if (phi[0][i][jp][k] == 1.0 && phi[0][i][j][k] < 1.0)
					{
						fnl = 0.0;
					}
					if (phi[0][i][jm][k] == 1.0 && phi[0][i][j][k] < 1.0)
					{
						fsl = 0.0;
					}
					if (phi[0][i][jp][k] < 1.0 && phi[0][i][j][k] == 1.0)
					{
						fnl = 0.0;
					}
					if (phi[0][i][jm][k] < 1.0 && phi[0][i][j][k] == 1.0)
					{
						fsl = 0.0;
					}
					fol = Dl * phi[0][i][j][k] * (conp[0][i][j][kp] - conp[0][i][j][k]) / dx;
					fil = Dl * phi[0][i][j][k] * (conp[0][i][j][k] - conp[0][i][j][km]) / dx;
					if (phi[0][i][j][kp] == 0.0)
					{
						fol = 0.0;
					}
					if (phi[0][i][j][km] == 0.0)
					{
						fil = 0.0;
					}

					cddttl = (fel - fwl + fnl - fsl + fol - fil) / dx;
					cddtts1 = (fes1 - fws1 + fns1 - fss1 + fos1 - fis1) / dx;

					if (phi[0][i][j][k] > 0.0)
					{
						conp2[0][i][j][k] = conp[0][i][j][k] + cddttl * dtime / phi[0][i][j][k];
					}
					else
					{
						conp2[0][i][j][k] = conp[0][i][j][k];
					}
					if (phi[1][i][j][k] > 0.0)
					{
						conp2[1][i][j][k] = conp[1][i][j][k] + cddtts1 * dtime / phi[1][i][j][k];
					}
					else
					{
						conp2[1][i][j][k] = conp[1][i][j][k];
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
					conp[0][i][j][k] = conp2[0][i][j][k];
					conp[1][i][j][k] = conp2[1][i][j][k];
					con[i][j][k] = conp[0][i][j][k] * phi[0][i][j][k] + conp[1][i][j][k] * phi[1][i][j][k];
					coni[i][j][k] = 0.0;
				}
			}
		}

#pragma omp barrier
		// front concnetration for driving force and solute trapping
		for (i = start; i <= end; i++)
		{
			for (j = 0; j <= ndmy; j++)
			{
				for (k = 0; k <= ndmz; k++)
				{
					if (phi[0][i][j][k] < 1.0 && phi[0][i][j][k] > 0.0)
					{

						coni[i][j][k] = conp[0][i][j][k];
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

void datasave1d(int step)
{
	FILE *stream; // ストリームのポインタ設定
	char buffer[30];
	sprintf(buffer, "data/phi/1d%d.csv", step);
	stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		fprintf(stream, "%lf   ", conp[0][ix][0][0]);
		fprintf(stream, "\n");
	}
	fclose(stream); // ファイルをクローズ

	FILE *streamc; // ストリームのポインタ設定
	char bufferc[30];
	sprintf(bufferc, "data/con/1d%d.csv", step);
	streamc = fopen(bufferc, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		fprintf(streamc, "%lf   ", con[ix][0][0]);
		fprintf(streamc, "\n");
	}
	fclose(streamc);
}

void datasave2d(int step)
{
	FILE *streamc; // ストリームのポインタ設定
	char bufferc[30];
	sprintf(bufferc, "data/con/2d%d.csv", step);
	streamc = fopen(bufferc, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			fprintf(streamc, "%lf   ", con[ix][iy][0]);
			fprintf(streamc, "\n");
		}
	}
	fclose(streamc); // ファイルをクローズ

	FILE *streamcp; // ストリームのポインタ設定
	char buffercp[30];
	sprintf(buffercp, "data/coni/2d%d.csv", step);
	streamcp = fopen(buffercp, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			fprintf(streamcp, "%lf   ", coni[ix][iy][0]);
			fprintf(streamcp, "\n");
		}
	}
	fclose(streamcp); // ファイルをクローズ

	FILE *streamp; // ストリームのポインタ設定
	char bufferp[30];
	sprintf(bufferp, "data/phi/2d%d.csv", step);
	streamp = fopen(bufferp, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			// fprintf(streamp, "%lf   ", [0][ix][iy][0]);
			fprintf(streamp, "%lf   ", pow(4.0 * phi[1][ix][iy][0] * (1 - phi[1][ix][iy][0]), 5.0));
			fprintf(streamp, "\n");
		}
	}
	fclose(streamp); // ファイルをクローズ
}

double calKv(double veli)
{
	double kv;
	double vd = 0.68;
	double vdb = 2.6;

	if (veli > vdb)
	{
		kv = 0.99;
	}
	else
	{
		kv = (kap * (1.0 - pow((veli / vdb), 2.0)) + veli / vd) / (1.0 - pow((veli / vdb), 2.0) + veli / vd);
	}

	// return kv;
	return kap;
}

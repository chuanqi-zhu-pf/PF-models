// Adapted from Prof. Koyama's MPF code in his textbook
// phase field coupled with concentration field
// Author: Chuanqi Zhu
// Created on: 2022/11/11

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>

#define NDX 128
#define NDY 64
#define NDZ 1

#define N 2

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int nm = N - 1;
double PI = 3.141592;
double RR = 8.3145;

double phi[N][NDX][NDY][NDZ], phi2[N][NDX][NDY][NDZ];
int phiIdx[N + 1][NDX][NDY][NDZ], phiNum[NDX][NDY][NDZ];
double aij[N][N], wij[N][N], mij[N][N], fij[N][N];
double anij[N][N], thij[N][N], vpij[N][N], etaij[N][N];

int phinum;
int i, j, k, ip, im, jp, jm, kp, km;
int ii, jj, kk;
int n1, n2, n3;

int istep, nstep, pstep;
double dtime, L, dx;
double M0;				 // 粒界の易動度
double W0;				 // ペナルティー項の係数
double A0;				 // 勾配エネルギー係数
double F0;				 // 粒界移動の駆動力
double temp0;			 // 温度
double sum1, sum2, sum3; // 各種の和の作業変数
double pddtt;			 // フェーズフィールドの時間変化率

double gamma0; // 粒界エネルギ密度
double delta;  // 粒界幅（差分ブロック数にて表現）
double mobi;   // 粒界の易動度
double vm0;	   // モル体積

double astre, astrem;
// double epsilon0;
// double termiikk, termjjkk;

// double phidx, phidy, phidz;
// double phidxx, phidyy, phidzz;

//--------------------------- variables and functions for concentration and temperature field --------------------------------
double conp[N][NDX][NDY][NDZ], conp2[N][NDX][NDY][NDZ];
double con[NDX][NDY][NDZ], coni[NDX][NDY][NDZ], intvel[NDX][NDY][NDZ];
double sumc, c0, dc0;
double cl, Tm, ce, ml1, kap1;
double c1e, c10e;
double dT, S0, Dl, Ds;
double dF;
double kvi, veli;

double flxci;
int dd, ddx, ddxp, ddy, ddyp;
int ddi, ddj, ddip, ddjp;

double phidxi, phidyi;
double nxi, nyi;
int di, xdi, xdip, xdim, ydi, ydip, ydim;
int ixdi, ixdip, jydi, jydip;

double fes1, fns1, fos1, fws1, fss1, fis1;
double fel, fnl, fol, fwl, fsl, fil;
double cddttl, cddtts1;

int intpos, intdis;
double velt;
int allL, intposp;

//----------------------------------------------------------------------------------------------

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

void datasave1d(int step),
	datasave2d(int step), datasave3d(int step);

double calKv(double veli);

// double calCs(double conp1, double conp0, double phi1, double 0);

//******* メインプログラム ******************************************
int main()
{
	nstep = 20001;
	pstep = 1000;

	dtime = 1.0e-12;
	dx = 1.0e-10;

	delta = 7.0;
	mobi = 2.6e8;
	gamma0 = 800.0;

	Dl = 1.5e-9;
	Ds = 3.0e-13;

	temp0 = 1603.0;
	Tm = 1685.0;
	ml1 = -400.0;
	kap1 = 0.3;
	cl = 0.09;
	c10e = -(Tm - temp0) / ml1;
	c1e = c10e * kap1;

	A0 = 8.0 * delta * gamma0 / PI / PI; // 勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gamma0 / delta;			 // ペナルティー項の係数[式(4.40)]
	M0 = mobi * PI * PI / (8.0 * delta); // 粒界の易動度[式(4.40)]

	//----------- variables for gradeint term with cubic anisotropy ------
	astre = 0.02;
	astrem = 0.1;
	for (ii = 0; ii <= nm; ii++)
	{
		for (jj = 0; jj <= nm; jj++)
		{
			anij[ii][jj] = 0.0;
			thij[ii][jj] = 0.0;
			vpij[ii][jj] = 0.0;
			etaij[ii][jj] = 0.0;
			if ((ii == 0) || (jj == 0))
			{
				anij[ii][jj] = 1.0;
			}
			if (ii == jj)
			{
				anij[ii][jj] = 0.0;
			}
		}
	}

	for (ii = 0; ii <= nm; ii++)
	{
		for (jj = 0; jj <= nm; jj++)
		{
			wij[ii][jj] = W0;
			aij[ii][jj] = A0;
			mij[ii][jj] = M0;
			fij[ii][jj] = 0.0;
			if ((ii == 0) || (jj == 0))
			{
				fij[ii][jj] = F0;
			}
			if (ii < jj)
			{
				fij[ii][jj] = -fij[ii][jj];
			}
			if (ii == jj)
			{
				wij[ii][jj] = 0.0;
				aij[ii][jj] = 0.0;
				mij[ii][jj] = 0.0;
				fij[ii][jj] = 0.0;
			}
		}
	}

	// Initialization of phase and concentration field
	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				// if (i < NDX / 32)
				if ((i) * (i) + (j) * (j) < NDX / 3 * NDX / 3)
				{
					phi[1][i][j][k] = 1.0;
					conp[1][i][j][k] = c1e;
					phi[0][i][j][k] = 0.0;
					conp[0][i][j][k] = c10e;
				}
				else
				{
					phi[1][i][j][k] = 0.0;
					conp[1][i][j][k] = c1e;
					phi[0][i][j][k] = 1.0;
					conp[0][i][j][k] = cl;
				}
				con[i][j][k] = conp[0][i][j][k] * phi[0][i][j][k] + conp[1][i][j][k] * phi[1][i][j][k];
			}
		}
	}

	std::cout << "diffusion stability for phi: " << dtime * mobi * gamma0 << std::endl;
	std::cout << "maximum permitted undercooling: " << gamma0 / delta << std::endl;
	std::cout << "diffusion stability for con: " << dtime / dx / dx * Dl << std::endl;

start:;

	if (istep % pstep == 0)
	{
		datasave2d(istep);
		std::cout << "the average concentration is: " << sumc / NDX / NDY << std::endl;
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
				intvel[i][j][k] = 0.0;
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
							phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 ;
							phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 ;
							phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 ;

							phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]); // フェーズフィールドの空間２階微分
							phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) ;
							phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) ;

							//------------------------------ Implementation of gradeint term with cubic anisotropy ---------------------

							phidxy = (phi[kk][ip][jp][k] + phi[kk][im][jm][k] - phi[kk][im][jp][k] - phi[kk][ip][jm][k]) / 4.0 ;
							phidxz = (phi[kk][ip][j][kp] + phi[kk][im][j][km] - phi[kk][im][j][kp] - phi[kk][ip][j][km]) / 4.0 ;
							phidyz = (phi[kk][i][jp][kp] + phi[kk][i][jm][km] - phi[kk][i][jm][kp] - phi[kk][i][jp][km]) / 4.0 ;

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
							dF = Tm + conp[0][i][j][k] * ml1 - temp0;
						}
						else if (ii == 0 && jj == 1)
						{
							dF = -(Tm + conp[0][i][j][k] * ml1 - temp0);
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
						pddtt += -2.0 * miijj  / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
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

	for (i = 0; i <= ndmx; i++)
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

	// Partition and Ejection
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
				if (phi[0][i][j][k] > 0.0 && phi2[0][i][j][k] < 1.0)
				{
					if (phi2[0][i][j][k] == 0.0)
					{
						phidxi = -(phi[1][ip][j][k] - phi[1][im][j][k]) / 2.0;
						phidyi = -(phi[1][i][jp][k] - phi[1][i][jm][k]) / 2.0;
						veli = (phi2[1][i][j][k] - phi[1][i][j][k]) / dtime / sqrt(phi[1][i][j][k]*(1.0-phi[1][i][j][k]))/PI * delta*dx;
					}
					else
					{
						phidxi = -(phi2[1][ip][j][k] - phi2[1][im][j][k]) / 2.0;
						phidyi = -(phi2[1][i][jp][k] - phi2[1][i][jm][k]) / 2.0;
						veli= (phi2[1][i][j][k] - phi[1][i][j][k]) / dtime / sqrt(phi2[1][i][j][k]*(1.0-phi2[1][i][j][k])) /PI * delta*dx;
					}

					nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi);
					nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi);

					kvi = calKv(veli);
					conp[1][i][j][k] = (conp[1][i][j][k] * phi[1][i][j][k] + (phi2[1][i][j][k] - phi[1][i][j][k]) * (kvi * conp[0][i][j][k])) / phi2[1][i][j][k];
					flxci = (phi2[1][i][j][k] - phi[1][i][j][k]) * conp[0][i][j][k] * (1.0 - kvi);

					di = 0;
					do
					{
						di++;
						xdi = int(round(nxi * di));
						ydi = int(round(nyi * di));

						ixdi = i + xdi;
						jydi = j + ydi;
						
						if (phi2[0][ixdi][jydi][k] >= 1.0)
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

	// redistribution of solute to each concentraton field
	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				if (phi[0][i][j][k] == 0.0)
				{
					conp[1][i][j][k] = con[i][j][k];
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

	// Diffusion equation
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
				if (phi[0][ip][j][k] == 0.0)
				{
					fel = 0.0;
				}
				if (phi[0][im][j][k] == 0.0)
				{
					fwl = 0.0;
				}
				// front obstacle
				if (phi[0][ip][j][k] == 1.0 && phi[0][i][j][k]<1.0)
				{
					fel = 0.0;
				}
				if (phi[0][im][j][k] == 1.0 && phi[0][i][j][k]<1.0)
				{
					fwl = 0.0;
				}
				if (phi[0][ip][j][k] < 1.0 && phi[0][i][j][k]==1.0)
				{
					fel = 0.0;
				}
				if (phi[0][im][j][k] < 1.0 && phi[0][i][j][k]==1.0)
				{
					fwl = 0.0;
				}
				
				fnl = Dl * phi[0][i][j][k] * (conp[0][i][jp][k] - conp[0][i][j][k]) / dx;
				fsl = Dl * phi[0][i][j][k] * (conp[0][i][j][k] - conp[0][i][jm][k]) / dx;
				if (phi[0][i][jp][k] == 0.0)
				{
					fnl = 0.0;
				}
				if (phi[0][i][jm][k] == 0.0)
				{
					fsl = 0.0;
				}
				if (phi[0][i][jp][k] == 1.0 && phi[0][i][j][k]<1.0)
				{
					fnl = 0.0;
				}
				if (phi[0][i][jm][k] == 1.0 && phi[0][i][j][k]<1.0)
				{
					fsl = 0.0;
				}
				if (phi[0][i][jp][k] < 1.0 && phi[0][i][j][k]==1.0)
				{
					fnl = 0.0;
				}
				if (phi[0][i][jm][k] < 1.0 && phi[0][i][j][k]==1.0)
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

	for (i = 0; i <= ndmx; i++)
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

	// moving frame
	intdis = intpos - NDX / 2;
	if (intdis > 0)
	{
		intposp -= intdis;
		for (i = 0; i <= ndmx; i++)
		{
			for (j = 0; j <= ndmy; j++)
			{
				for (k = 0; k <= ndmz; k++)
				{
					if ((i + intdis) <= ndmx)
					{
						phi[0][i][j][k] = phi[0][i + intdis][j][k];
						phi[1][i][j][k] = phi[1][i + intdis][j][k];
						phi2[0][i][j][k] = phi2[0][i + intdis][j][k];
						phi2[1][i][j][k] = phi2[1][i + intdis][j][k];
						conp[0][i][j][k] = conp[0][i + intdis][j][k];
						conp[1][i][j][k] = conp[1][i + intdis][j][k];
						con[i][j][k] = con[i + intdis][j][k];
					}
					else
					{
						phi[0][i][j][k] = 1.0;
						phi[1][i][j][k] = 0.0;
						phi2[0][i][j][k] = 1.0;
						phi2[1][i][j][k] = 0.0;
						conp[0][i][j][k] = cl;
						conp[1][i][j][k] = c1e;
						con[i][j][k] = cl;
					}
				}
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
	sprintf(buffer, "data/phi/1d%d.csv", step);
	stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		fprintf(stream, "%lf   ", conp[0][i][0][0]);
		fprintf(stream, "\n");
	}
	fclose(stream); // ファイルをクローズ

	FILE *streamc; // ストリームのポインタ設定
	char bufferc[30];
	sprintf(bufferc, "data/con/1d%d.csv", step);
	streamc = fopen(bufferc, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		fprintf(streamc, "%lf   ", con[i][0][0]);
		fprintf(streamc, "\n");
	}
	fclose(streamc);
}

void datasave2d(int step)
{
	FILE *stream; // ストリームのポインタ設定
	char buffer[30];
	sprintf(buffer, "data/phi/2d%d.csv", step);
	stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(stream, "%lf   ", conp[0][i][j][0]);
			fprintf(stream, "\n");
		}
	}
	fclose(stream); // ファイルをクローズ

	FILE *streamc; // ストリームのポインタ設定
	char bufferc[30];
	sprintf(bufferc, "data/con/2d%d.csv", step);
	streamc = fopen(bufferc, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(streamc, "%lf   ", con[i][j][0]);
			fprintf(streamc, "\n");
		}
	}
	fclose(streamc); // ファイルをクローズ
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
		kv = (kap1 * (1.0 - pow((veli / vdb), 2.0)) + veli / vd) / (1.0 - pow((veli / vdb), 2.0) + veli / vd);
	}

	// return kv;
	return kap1;
}

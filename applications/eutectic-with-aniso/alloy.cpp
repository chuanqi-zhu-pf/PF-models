// Adapted from Prof. Koyama's MPF code in his textbook
// phase field coupled with concentration field
// Author: Chuanqi Zhu
// Created on: 2022/11/11

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>

#define NDX 128
#define NDY 64
#define NDZ 1

#define N 3

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

// double astre;
// double epsilon0;
// double termiikk, termjjkk;

// double phidx, phidy, phidz;
double phidxx, phidyy, phidzz;
double miijj;

//--------------------------- variables and functions for concentration and temperature field --------------------------------
double conp[N][NDX][NDY][NDZ], conp2[N][NDX][NDY][NDZ];
double con[NDX][NDY][NDZ], coni[NDX][NDY][NDZ];
double sumc, suml, c0, dc0;
double cl, Te, ce, ml1, ml2, kap1, kap2;
double c1e, c10e, c2e, c20e;
double dT, S0, Dl, Ds;
double dF;

double flxci;
int dd, ddx, ddxp, ddy, ddyp;
int ddi, ddj, ddip, ddjp;

double phidxi, phidyi;
double nxi, nyi;
int di, xdi, xdip, xdim, ydi, ydip, ydim;
int ixdi, ixdip, jydi, jydip;
double coni1, coni2;

double fes1, fns1, fos1, fws1, fss1, fis1;
double fes2, fns2, fos2, fws2, fss2, fis2;
double fel, fnl, fol, fwl, fsl, fil;
double cddttl, cddtts1, cddtts2;

//----------------------------------------------------------------------------------------------

//----------- variables for gradeint term with anisotropy ------

// common variables
double anij[N][N], thij[N][N], vpij[N][N], etaij[N][N];
double th, vp, eta;
double epsilon0;
double ep, epdx, epdy, epdz;
double termx, termy, termz;
double termiikk, termjjkk, wwiikk, wwjjkk;

double phidx, phidy, phidz, phiabs2, phiabs;
double phidxy, phidyz, phidxz;
double phidyx, phidzy, phidzx;

// Cubic anisotropy
double astre;
double xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp;
double phidxp, phidyp, phidzp;
double phidxpx, phidypx, phidzpx;
double phidxpy, phidypy, phidzpy;
double phidxpz, phidypz, phidzpz;

double term0;
double termx0, termx1, termx0dx, termx1dx;
double termy0, termy1, termy0dy, termy1dy;
double termz0, termz1, termz0dz, termz1dz;

double epdphix, epdphiy, epdphiz;
double epdphixdx, epdphiydy, epdphizdz;
double nxii, nyii, nzii;
double phidxii, phidyii, phidzii, phiabs2ii;
double astrem, phidxm, phidym, phidzm, phidxpm, phidypm, phidzpm, phiabsm2;

// Faceted anisotropy
#define FN 4
double face[FN][3];
double del, coef, al0, zeta1, rp0, rp1;

double dphiabs2dx, dphiabs2dy, dphiabs2dz;
double dphiabs2dphix, dphiabs2dphiy, dphiabs2dphiz;
double dphiabs2dphixdx, dphiabs2dphiydy, dphiabs2dphizdz;

double min_val;
int min_idx, l;
double a1, a2, a3;
double ux0, uy0, uz0;
double ux, uy, uz, uu;
double al;

double PP, QQ, CC, SS;

double nx, ny, nz;
double nxx, nxy, nxz;
double nyx, nyy, nyz;
double nzx, nzy, nzz;

double nxphix, nyphix, nzphix;
double nxphiy, nyphiy, nzphiy;
double nxphiz, nyphiz, nzphiz;

double nxphixdx, nyphixdx, nzphixdx;
double nxphiydy, nyphiydy, nzphiydy;
double nxphizdz, nyphizdz, nzphizdz;

double dPdx, dPdy, dPdz;
double dQdx, dQdy, dQdz;
double dPdphix, dPdphiy, dPdphiz;
double dQdphix, dQdphiy, dQdphiz;
double dPdphixdx, dPdphiydy, dPdphizdz;
double dQdphixdx, dQdphiydy, dQdphizdz;

//--------------------------------------------------------------------

void datasave1d(int step),
	datasave2d(int step), datasave3d(int step);

//******* メインプログラム ******************************************
int main()
{
	nstep = 40001;
	pstep = 1000;
	dtime = 5.0;
	temp0 = 1000.0;
	L = 2000.0;
	vm0 = 7.0e-6;
	delta = 5.0;
	mobi = 1.0;

	dx = L / double(NDX) * 1.0e-9;		  // 差分プロック１辺の長さ(m)
	gamma0 = 0.5 * vm0 / RR / temp0 / dx; // 粒界エネルギ密度（0.5J/m^2）を無次元化
	A0 = 8.0 * delta * gamma0 / PI / PI;  // 勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gamma0 / delta;			  // ペナルティー項の係数[式(4.40)]
	M0 = mobi * PI * PI / (8.0 * delta);  // 粒界の易動度[式(4.40)]
	F0 = 50.0 / RR / temp0;				  // 粒界移動の駆動力

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

	S0 = 4.5 / RR / temp0;
	Te = 1005.0;
	ce = 0.5;
	ml1 = -100.0;
	ml2 = 100.0;
	kap1 = 0.2;
	kap2 = 1.8;
	Dl = 0.1e-16;
	Ds = 0.2e-18;
	cl = 0.45;
	c10e = -(Te - temp0) / ml1 + ce;
	c20e = -(Te - temp0) / ml2 + ce;
	c1e = c10e * kap1;
	c2e = c20e * kap2;

	//----------- variables for gradeint term with cubic anisotropy ------

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
				anij[ii][jj] = 0.0;
			}
			if (ii == jj)
			{
				anij[ii][jj] = 0.0;
			}
		}
	}

	anij[1][0] = 1.0;
	anij[0][1] = 1.0;

	anij[2][0] = 2.0;
	anij[0][2] = 2.0;

	thij[1][0] = PI / 8.0;
	thij[0][1] = PI / 8.0;

	// cubic anisotropy
	astre = 0.03;
	astrem = 0.1;

	// faceted anisotropy
	del = 0.36;
	coef = 1.0 / 1.62454431;
	al0 = 54.7 / 180.0 * PI;
	zeta1 = 0.35;
	rp0 = 0.1;
	rp1 = 0.1;

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
	//--------------------------------------------------------------------

	// Initialization of phase and concentration field
	sumc = 0.0;
	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				// if (i < NDX / 16 && j < NDY / 2)
				if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY * 3 / 4) * (j - NDY * 3 / 4) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 8 * NDX / 8)
				{
					phi[1][i][j][k] = 1.0;
					conp[1][i][j][k] = c1e;
					phi[2][i][j][k] = 0.0;
					conp[2][i][j][k] = c2e;
					phi[0][i][j][k] = 0.0;
					conp[0][i][j][k] = c10e;
				}
				// else if (i < NDX / 16 && j >= NDY / 2)
				else if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 4) * (j - NDY / 4) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 8 * NDX / 8)
				{
					phi[1][i][j][k] = 0.0;
					conp[1][i][j][k] = c1e;
					phi[2][i][j][k] = 1.0;
					conp[2][i][j][k] = c2e;
					phi[0][i][j][k] = 0.0;
					conp[0][i][j][k] = c20e;
				}
				else
				{
					phi[1][i][j][k] = 0.0;
					conp[1][i][j][k] = c1e;
					phi[2][i][j][k] = 0.0;
					conp[2][i][j][k] = c2e;
					phi[0][i][j][k] = 1.0;
					conp[0][i][j][k] = cl;
				}
				con[i][j][k] = conp[0][i][j][k] * phi[0][i][j][k] + conp[1][i][j][k] * phi[1][i][j][k] + conp[2][i][j][k] * phi[2][i][j][k];
				sumc += con[i][j][k];
			}
		}
	}

start:;

	if (istep % pstep == 0)
	{
		datasave2d(istep);
		std::cout << sumc / NDX / NDY / NDZ << std::endl;
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
					ip = 0;
				}
				if (i == 0)
				{
					im = ndmx;
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
					ip = 0;
				}
				if (i == 0)
				{
					im = ndmx;
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
							phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]); // フェーズフィールドの空間２階微分
							phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]);
							phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]);

							//------------------------------ Implementation of gradeint term with cubic anisotropy ---------------------

							phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0;
							phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0;
							phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0;

							phidxy = (phi[kk][ip][jp][k] + phi[kk][im][jm][k] - phi[kk][im][jp][k] - phi[kk][ip][jm][k]) / 4.0;
							phidxz = (phi[kk][ip][j][kp] + phi[kk][im][j][km] - phi[kk][im][j][kp] - phi[kk][ip][j][km]) / 4.0;
							phidyz = (phi[kk][i][jp][kp] + phi[kk][i][jm][km] - phi[kk][i][jm][kp] - phi[kk][i][jp][km]) / 4.0;

							phiabs2 = phidx * phidx + phidy * phidy + phidz * phidz;
							phiabs = sqrt(phiabs2);

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

								wwiikk = wij[ii][kk] * pow((1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 2.0)), 2.0);
								termiikk = term0 + termx + termy + termz;
							}
							else if (anij[ii][kk] == 2.0 && phiabs2 != 0.0)
							{
								phidyx = phidxy;
								phidzx = phidxz;
								phidzy = phidyz;

								dphiabs2dx = 2.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz);
								dphiabs2dy = 2.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz);
								dphiabs2dz = 2.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz);

								dphiabs2dphix = 2.0 * phidx;
								dphiabs2dphiy = 2.0 * phidy;
								dphiabs2dphiz = 2.0 * phidz;

								dphiabs2dphixdx = 2.0 * phidxx;
								dphiabs2dphiydy = 2.0 * phidyy;
								dphiabs2dphizdz = 2.0 * phidzz;

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

									th = 45.0 / 180.0 * PI;
									a1 = 1.0;
									a2 = 0.0;
									a3 = 0.0;

									ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
									uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
									uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

									ux0 = ux;
									uy0 = uy;
									uz0 = uz;

									uu = ux * ux + uy * uy + uz * uz;
									al = acos(abs(nx * ux + ny * uy + nz * uz) / sqrt(uu));

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

								th = 45.0 / 180.0 * PI;
								a1 = 1.0;
								a2 = 0.0;
								a3 = 0.0;

								ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
								uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
								uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

								ux0 = ux;
								uy0 = uy;
								uz0 = uz;

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

								wwiikk = wij[ii][kk] * pow((1.0 + del * (CC + tan(al0) * SS)) * coef, 2.0);
								termiikk = termx + termy + termz;
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

								wwjjkk = wij[jj][kk] * pow((1.0 - 3.0 * astre + 4.0 * astre * (pow(phidxp, 4.0) + pow(phidyp, 4.0) + pow(phidzp, 4.0)) / pow(phiabs2, 2.0)), 2.0);
								termjjkk = term0 + termx + termy + termz;
							}
							else if (anij[jj][kk] == 2.0 && phiabs2 != 0.0)
							{
								phidyx = phidxy;
								phidzx = phidxz;
								phidzy = phidyz;

								dphiabs2dx = 2.0 * (phidx * phidxx + phidy * phidxy + phidz * phidxz);
								dphiabs2dy = 2.0 * (phidx * phidxy + phidy * phidyy + phidz * phidyz);
								dphiabs2dz = 2.0 * (phidx * phidxz + phidy * phidyz + phidz * phidzz);

								dphiabs2dphix = 2.0 * phidx;
								dphiabs2dphiy = 2.0 * phidy;
								dphiabs2dphiz = 2.0 * phidz;

								dphiabs2dphixdx = 2.0 * phidxx;
								dphiabs2dphiydy = 2.0 * phidyy;
								dphiabs2dphizdz = 2.0 * phidzz;

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

									th = 45.0 / 180.0 * PI;
									a1 = 1.0;
									a2 = 0.0;
									a3 = 0.0;

									ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
									uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
									uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

									ux0 = ux;
									uy0 = uy;
									uz0 = uz;

									uu = ux * ux + uy * uy + uz * uz;
									al = acos(abs(nx * ux + ny * uy + nz * uz) / sqrt(uu));

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

								th = 45.0 / 180.0 * PI;
								a1 = 1.0;
								a2 = 0.0;
								a3 = 0.0;

								ux = ux0 * cos(th) + (a2 * uz0 - a3 * uy0) * sin(th) + a1 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
								uy = uy0 * cos(th) + (a3 * ux0 - a1 * uz0) * sin(th) + a2 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));
								uz = uz0 * cos(th) + (a1 * uy0 - a2 * ux0) * sin(th) + a3 * (a1 * ux0 + a2 * uy0 + a3 * uz0) * (1.0 - cos(th));

								ux0 = ux;
								uy0 = uy;
								uz0 = uz;

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

								wwjjkk = wij[jj][kk] * pow((1.0 + del * (CC + tan(al0) * SS)) * coef, 2.0);
								termjjkk = termx + termy + termz;
							}
							else
							{
								wwjjkk = wij[jj][kk];
								termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
							}

							sum1 += 0.5 * (termiikk - termjjkk) + (wwiikk - wwjjkk) * phi[kk][i][j][k];
						}

						if (ii == 1 && jj == 0)
						{
							dF = ((coni[i][j][k] - ce) * ml1 + Te - temp0) * S0;
						}
						else if (ii == 0 && jj == 1)
						{
							dF = -(((coni[i][j][k] - ce) * ml1 + Te - temp0)) * S0;
						}
						else if (ii == 2 && jj == 0)
						{
							dF = ((coni[i][j][k] - ce) * ml2 + Te - temp0) * S0;
						}
						else if (ii == 0 && jj == 2)
						{
							dF = -(((coni[i][j][k] - ce) * ml2 + Te - temp0)) * S0;
						}
						else if (ii + jj == 3)
						{
							dF = 0.0;
						}
						miijj = mij[ii][jj];
						if (ii + jj == 3)
						{
							miijj = 0.000 * mij[ii][jj];
						}
						// if ((ii == 1 && jj == 0) || (ii == 0 && jj == 1))
						// {
						// 	th = thij[ii][jj];
						// 	vp = vpij[ii][jj];
						// 	eta = etaij[ii][jj];

						// 	xxp = cos(th) * cos(vp);
						// 	yxp = sin(th) * cos(vp);
						// 	zxp = sin(vp);
						// 	xyp = -sin(th) * cos(eta) - cos(th) * sin(vp) * sin(eta);
						// 	yyp = cos(th) * cos(eta) - sin(th) * sin(vp) * sin(eta);
						// 	zyp = cos(vp) * sin(eta);
						// 	xzp = sin(eta) * sin(th) - cos(eta) * cos(th) * sin(vp);
						// 	yzp = -sin(eta) * cos(th) - cos(eta) * sin(th) * sin(vp);
						// 	zzp = cos(eta) * cos(vp);

						// 	phidxm = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0;
						// 	phidym = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0;
						// 	phidzm = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0;

						// 	phidxpm = phidxm * xxp + phidym * yxp + phidzm * zxp;
						// 	phidypm = phidxm * xyp + phidym * yyp + phidzm * zyp;
						// 	phidzpm = phidxm * xzp + phidym * yzp + phidzm * zzp;

						// 	phiabsm2 = phidxpm * phidxpm + phidypm * phidypm + phidzpm * phidzpm;

						// 	miijj = mij[ii][jj] * (1.0 - 3.0 * astrem + 4.0 * astrem * (pow(phidxpm, 4.0) + pow(phidypm, 4.0) + pow(phidzpm, 4.0)) / pow(phiabsm2, 2.0));
						// }

						pddtt += -2.0 * miijj / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
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
				coni[i][j][k] = 0.0;
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
				if (phi[0][i][j][k] > 0.0 && phi2[0][i][j][k] < 1.0)
				{
					// calculate the inteface normal vector
					if (phi2[0][i][j][k] == 0.0)
					{
						phidxi = (phi[0][ip][j][k] - phi[0][im][j][k]) / 2.0 / dx;
						phidyi = (phi[0][i][jp][k] - phi[0][i][jm][k]) / 2.0 / dx;
					}
					else
					{
						phidxi = (phi2[0][ip][j][k] - phi2[0][im][j][k]) / 2.0 / dx;
						phidyi = (phi2[0][i][jp][k] - phi2[0][i][jm][k]) / 2.0 / dx;
					}

					nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi);
					nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi);

					// update the local solid concentration
					if (phi2[1][i][j][k] != 0.0)
					{
						conp[1][i][j][k] = (conp[1][i][j][k] * phi[1][i][j][k] + (phi2[1][i][j][k] - phi[1][i][j][k]) * (kap1 * conp[0][i][j][k])) / phi2[1][i][j][k];
					}
					else
					{
						conp[1][i][j][k] = c1e;
					}
					if (phi2[2][i][j][k] != 0.0)
					{
						conp[2][i][j][k] = (conp[2][i][j][k] * phi[2][i][j][k] + (phi2[2][i][j][k] - phi[2][i][j][k]) * (kap2 * conp[0][i][j][k])) / phi2[2][i][j][k];
					}
					else
					{
						conp[2][i][j][k] = c2e;
					}

					// calculate the ejection flux
					flxci = (phi2[1][i][j][k] - phi[1][i][j][k]) * conp[0][i][j][k] * (1.0 - kap1) + (phi2[2][i][j][k] - phi[2][i][j][k]) * conp[0][i][j][k] * (1.0 - kap2);

					// Search the front liquid and cast out the ejected solute
					di = 0;
					do
					{
						di++;
						xdi = int(round(nxi * di));
						ydi = int(round(nyi * di));

						ixdi = i + xdi;
						if (ixdi > ndmx)
						{
							ixdi = ixdi - ndmx - 1;
						}
						if (ixdi < 0)
						{
							ixdi = ndmx + ixdi + 1;
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
					if (phi[1][i][j][k] < 1.0 && phi[2][i][j][k] < 1.0)
					{
						conp[1][i][j][k] = c1e;
						conp[2][i][j][k] = c2e;
					}
					else
					{
						conp[1][i][j][k] = con[i][j][k] * phi[1][i][j][k];
						conp[2][i][j][k] = con[i][j][k] * phi[2][i][j][k];
					}
					conp[0][i][j][k] = c10e * phi[1][i][j][k] + c20e * phi[2][i][j][k];
				}
				else if (phi[0][i][j][k] == 1.0)
				{
					conp[1][i][j][k] = c1e;
					conp[2][i][j][k] = c2e;
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
					ip = 0;
				}
				if (i == 0)
				{
					im = ndmx;
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
				fos2 = Ds * (phi[2][i][j][k] * (conp[2][i][j][kp] - conp[2][i][j][k]) / dx);
				fis2 = Ds * (phi[2][i][j][k] * (conp[2][i][j][k] - conp[2][i][j][km]) / dx);
				if (phi[2][i][j][kp] == 0.0)
				{
					fos2 = 0.0;
				}
				if (phi[2][i][j][km] == 0.0)
				{
					fis2 = 0.0;
				}
				//* pow((1.0 - phi[0][i][j][k] * phi[1][i][j][k] * phi[2][i][j][k]), 10.0)
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
				cddtts2 = (fes2 - fws2 + fns2 - fss2 + fos2 - fis2) / dx;

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
	sumc = 0.0;
	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				conp[0][i][j][k] = conp2[0][i][j][k];
				conp[1][i][j][k] = conp2[1][i][j][k];
				conp[2][i][j][k] = conp2[2][i][j][k];
				con[i][j][k] = conp[0][i][j][k] * phi[0][i][j][k] + conp[1][i][j][k] * phi[1][i][j][k] + conp[2][i][j][k] * phi[2][i][j][k];
				if (con[i][j][k] > 1.0)
				{
					con[i][j][k] = 1.0;
				}
				if (con[i][j][k] < 0.0)
				{
					con[i][j][k] = 0.0;
				}
				sumc += con[i][j][k];
			}
		}
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
					ip = 0;
				}
				if (i == 0)
				{
					im = ndmx;
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
					phidxi = (phi[0][ip][j][k] - phi[0][im][j][k]) / 2.0;
					phidyi = (phi[0][i][jp][k] - phi[0][i][jm][k]) / 2.0;
					nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi);
					nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi);

					di = 0;
					do
					{
						di++;
						xdi = int(round(nxi * di));
						ydi = int(round(nyi * di));

						ixdi = i + xdi;
						if (ixdi > ndmx)
						{
							ixdi = ixdi - ndmx - 1;
						}
						if (ixdi < 0)
						{
							ixdi = ndmx + ixdi + 1;
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

						if (phi[0][ixdi][jydi][k] == 1.0)
						{
							coni[i][j][k] = conp[0][ixdi][jydi][k];
							break;
						}

					} while (di < int(delta));

					if (phi[0][ixdi][jydi][k] != 1.0)
					{
						coni[i][j][k] = conp[0][i][j][k];
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
		fprintf(stream, "%lf   ", phi[0][i][0][0]);
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
	fclose(streamc); // ファイルをクローズ
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
			fprintf(stream, "%lf   ", coni[i][j][0]);
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

	FILE *streamcl; // ストリームのポインタ設定
	char buffercl[30];
	sprintf(buffercl, "data/conl/2d%d.csv", step);
	streamcl = fopen(buffercl, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(streamcl, "%lf   ", conp[0][i][j][0]);
			fprintf(streamcl, "\n");
		}
	}
	fclose(streamcl); // ファイルをクローズ
}

void datasave3d(int step)
{
	FILE *stream;
	char buffer[30];
	sprintf(buffer, "data/3d%d.vtk", step);
	stream = fopen(buffer, "a");

	fprintf(stream, "# vtk DataFile Version 1.0\n");
	fprintf(stream, "phi_%d.vtk\n", step);
	fprintf(stream, "ASCII\n");
	fprintf(stream, "DATASET STRUCTURED_POINTS\n");
	fprintf(stream, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
	fprintf(stream, "ORIGIN 0.0 0.0 0.0\n");
	fprintf(stream, "ASPECT_RATIO 1.0 1.0 1.0\n");
	fprintf(stream, "\n");
	fprintf(stream, "POINT_DATA %d\n", NDX * NDY * NDZ);
	fprintf(stream, "SCALARS scalars float\n");
	fprintf(stream, "LOOKUP_TABLE default\n");

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				fprintf(stream, "%e\n", phi[0][i][j][k]);
			}
		}
	}
	fclose(stream);
}

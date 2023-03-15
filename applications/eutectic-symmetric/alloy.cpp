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

double astre;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidz;
double phidxx, phidyy, phidzz;

double miijj;

//--------------------------- variables and functions for concentration and temperature field --------------------------------
double conp[N][NDX][NDY][NDZ], conp2[N][NDX][NDY][NDZ];
double con[NDX][NDY][NDZ], coni[NDX][NDY][NDZ];
double sumc, suml, c00, c0, dc00;
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

double fes1, fns1, fos1, fws1, fss1, fis1;
double fes2, fns2, fos2, fws2, fss2, fis2;
double fel, fnl, fol, fwl, fsl, fil;
double cddttl, cddtts1, cddtts2;

//----------------------------------------------------------------------------------------------

void datasave1d(int step),
	datasave2d(int step), datasave3d(int step);

//******* メインプログラム ******************************************
int main()
{
	nstep = 12001;
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

	S0 = 10.5 / RR / temp0;
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

	// Initialization of phase and concentration field
	sumc = 0.0;
	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				if (i < NDX / 16 && j < NDY / 4)
				// if ((i - NDX / 4) * (i - NDX / 4) + (j - NDY / 4) * (j - NDY / 4) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 5 * NDX / 5)
				{
					phi[1][i][j][k] = 1.0;
					conp[1][i][j][k] = c1e;
					phi[2][i][j][k] = 0.0;
					conp[2][i][j][k] = c2e;
					phi[0][i][j][k] = 0.0;
					conp[0][i][j][k] = c10e;
				}
				else if (i < NDX / 16 && j >= NDY / 4)
				// else if ((i - NDX * 3 / 4) * (i - NDX * 3 / 4) + (j - NDY / 4) * (j - NDY / 4) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 5 * NDX / 5)
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
	c00 = sumc / NDX / NDY / NDZ;

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
							phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0;
							phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0;
							phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0;

							phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]); // フェーズフィールドの空間２階微分
							phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]);
							phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]);

							termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
							termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);

							sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
						}

						if (ii == 1 && jj == 0)
						{
							dF = ((conp[0][i][j][k] - ce) * ml1 + Te - temp0) * S0;
						}
						else if (ii == 0 && jj == 1)
						{
							dF = -(((conp[0][i][j][k] - ce) * ml1 + Te - temp0)) * S0;
						}
						else if (ii == 2 && jj == 0)
						{
							dF = ((conp[0][i][j][k] - ce) * ml2 + Te - temp0) * S0;
						}
						else if (ii == 0 && jj == 2)
						{
							dF = -(((conp[0][i][j][k] - ce) * ml2 + Te - temp0)) * S0;
						}
						else if (ii + jj == 3)
						{
							dF = 0.0;
						}
						miijj = mij[ii][jj];
						if (ii + jj == 3)
						{
							miijj = 0.01 * mij[ii][jj];
						}
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
				// if (phi[0][ip][j][k] == 1.0 && phi[0][i][j][k] < 1.0)
				// {
				// 	fel = 0.0;
				// }
				// if (phi[0][im][j][k] == 1.0 && phi[0][i][j][k] < 1.0)
				// {
				// 	fwl = 0.0;
				// }
				// if (phi[0][ip][j][k] < 1.0 && phi[0][i][j][k] == 1.0)
				// {
				// 	fel = 0.0;
				// }
				// if (phi[0][im][j][k] < 1.0 && phi[0][i][j][k] == 1.0)
				// {
				// 	fwl = 0.0;
				// }
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
				// front obstacle
				// if (phi[0][i][jp][k] == 1.0 && phi[0][i][j][k] < 1.0)
				// {
				// 	fnl = 0.0;
				// }
				// if (phi[0][i][jm][k] == 1.0 && phi[0][i][j][k] < 1.0)
				// {
				// 	fsl = 0.0;
				// }
				// if (phi[0][i][jp][k] < 1.0 && phi[0][i][j][k] == 1.0)
				// {
				// 	fnl = 0.0;
				// }
				// if (phi[0][i][jm][k] < 1.0 && phi[0][i][j][k] == 1.0)
				// {
				// 	fsl = 0.0;
				// }
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
	suml = 0.0;
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
				sumc += con[i][j][k];
				if (phi[0][i][j][k] == 1.0)
				{
					suml += phi[0][i][j][k];
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
			fprintf(stream, "%lf   ", conp[1][i][j][0]);
			fprintf(stream, "\n");
		}
	}
	fclose(stream); // ファイルをクローズ

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

	FILE *streamc0; // ストリームのポインタ設定
	char bufferc0[30];
	sprintf(bufferc0, "data/ave_con.csv", step);
	streamc0 = fopen(bufferc0, "a"); // 書き込む先のファイルを追記方式でオープン

	fprintf(streamc0, "%lf   ", sumc / NDX / NDY / NDZ);
	fprintf(streamc0, "\n");
	fclose(streamc); // ファイルをクローズ
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

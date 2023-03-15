// Adapted from Prof. Koyama's MPF code in his textbook
// phase field coupled with concentration field
// Author: Chuanqi Zhu
// Created on: 2022/11/11

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NDX 128
#define NDY 1
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

int phinum;
int i, j, k, ip, im, jp, jm, kp, km;
int ii, jj, kk;
int n1, n2, n3;

int istep, nstep, pstep;
double dtime, L, dx;
double M0;				 //粒界の易動度
double W0;				 //ペナルティー項の係数
double A0;				 //勾配エネルギー係数
double F0;				 //粒界移動の駆動力
double temp0;			 //温度
double sum1, sum2, sum3; //各種の和の作業変数
double pddtt;			 //フェーズフィールドの時間変化率

double gamma0; //粒界エネルギ密度
double delta;  //粒界幅（差分ブロック数にて表現）
double mobi;   //粒界の易動度
double vm0;	   //モル体積

double astre;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidz;
double phidxx, phidyy, phidzz;

int intpos, intdis;

//--------------------------- variables and functions for concentration and temperature field --------------------------------
double conp[N][NDX][NDY][NDZ], conp2[N][NDX][NDY][NDZ];
double con[NDX][NDY][NDZ], coni[NDX][NDY][NDZ], veli[NDX][NDY][NDZ];
double temp[NDX][NDY][NDZ];
double sumc, c0, dc0;
double cl, Te, Tm, ce, mle, kap1;
double c1e, c10e;
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
double fel, fnl, fol, fwl, fsl, fil;
double cddttl, cddtts1;

//----------------------------------------------------------------------------------------------

double vd, vdb, ke, kv, mv;
double Tr, Tg;

void datasave1d(int step),
	datasave2d(int step), datasave3d(int step);

double calC1e(double temp00), calC10e(double temp00);

//******* メインプログラム ******************************************
int main()
{
	nstep = 100001;
	pstep = 10000;
	dtime = 1.0e-1;
	temp0 = 1010.0;
	L = 16.0;
	vm0 = 7.0e-6;
	delta = 5.0;
	mobi = 1.0;

	dx = L * 1.0e-9;					  //差分プロック１辺の長さ(m)
	gamma0 = 0.5 * vm0 / RR / temp0 / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
	A0 = 8.0 * delta * gamma0 / PI / PI;  //勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gamma0 / delta;			  //ペナルティー項の係数[式(4.40)]
	M0 = mobi * PI * PI / (8.0 * delta);  //粒界の易動度[式(4.40)]
	F0 = 50.0 / RR / temp0;				  //粒界移動の駆動力

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

	Tm = 1050.0;
	ke = 0.5;
	mle = -(Tm - 1000.0) / 0.5;
	vd = 0.68;
	vdb = 2.6;

	S0 = 100.5 / RR / 1000.0;
	ce = 0.5;
	kap1 = 0.5;
	Dl = 0.1e-15;
	Ds = 0.2e-18;
	cl = 0.350;
	c10e = -(Tm - temp0) / mle;
	c1e = c10e * ke;

	Tr = 0.00000;
	Tg = 0.00;

	// Initialization of phase and concentration field
	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			for (k = 0; k <= ndmz; k++)
			{
				temp[i][j][k] = temp0 + Tg * (i - NDX / 16);
				if (i < NDX / 16)
				// if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 8 * NDX / 8)
				{
					phi[1][i][j][k] = 1.0;
					conp[1][i][j][k] = calC1e(temp[i][j][k]);
					phi[0][i][j][k] = 0.0;
					conp[0][i][j][k] = calC10e(temp[i][j][k]);
				}
				else if (i < (NDX * 4 / 64 + 1) && i >= NDX / 16)
				{
					phi[1][i][j][k] = 0.0;
					conp[1][i][j][k] = calC1e(temp[i][j][k]);
					phi[0][i][j][k] = 1.0;
					conp[0][i][j][k] = calC10e(temp[i][j][k]);
				}
				else
				{
					phi[1][i][j][k] = 0.0;
					conp[1][i][j][k] = calC1e(temp[i][j][k]);
					phi[0][i][j][k] = 1.0;
					conp[0][i][j][k] = cl;
				}
				con[i][j][k] = conp[0][i][j][k] * phi[0][i][j][k] + conp[1][i][j][k] * phi[1][i][j][k];
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

							phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]); //フェーズフィールドの空間２階微分
							phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]);
							phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]);

							termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
							termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);

							sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
						}

						if (ii == 1 && jj == 0)
						{
							dF = ((coni[i][j][k] - c10e) * mle) * S0;
						}
						else if (ii == 0 && jj == 1)
						{
							dF = -((coni[i][j][k] - c10e) * mle) * S0;
						}
						pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
						//フェーズフィールドの発展方程式[式(4.31)]
					}
					phi2[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）
					if (phi2[ii][i][j][k] >= 1.0)
					{
						phi2[ii][i][j][k] = 1.0;
					} //フェーズフィールドの変域補正
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
				veli[i][j][k] = 0.0;
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
					}
					else
					{
						phidxi = -(phi2[1][ip][j][k] - phi2[1][im][j][k]) / 2.0;
						phidyi = -(phi2[1][i][jp][k] - phi2[1][i][jm][k]) / 2.0;
					}

					nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi);
					nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi);

					veli[i][j][k] = (phi2[1][i][j][k] - phi[1][i][j][k]) / sqrt(phidxi * phidxi + phidyi * phidyi) / dtime;

					di = 0;
					do
					{
						di++;
						xdi = int(round(nxi * di));
						xdip = int(round(nxi * (di + 1)));
						ydi = int(round(nyi * di));
						ydip = int(round(nyi * (di + 1)));

						ixdi = i + xdi;
						if (ixdi > ndmx)
						{
							ixdi = ixdi - ndmx - 1;
						}
						if (ixdi < 0)
						{
							ixdi = ndmx + ixdi + 1;
						}
						ixdip = i + xdip;
						if (ixdip > ndmx)
						{
							ixdip = ixdip - ndmx - 1;
						}
						if (ixdip < 0)
						{
							ixdip = ndmx + ixdip + 1;
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

						jydip = j + ydip;
						if (jydip > ndmy)
						{
							jydip = jydip - ndmy - 1;
						}
						if (jydip < 0)
						{
							jydip = ndmy + jydip + 1;
						}

						if (phi2[0][ixdi][jydi][k] >= phi2[0][ixdip][jydip][k])
						{
							flxci = (phi[1][i][j][k] - phi2[1][i][j][k]) * (calC1e(temp[i][j][k]) - conp[0][i][j][k]);
							con[i][j][k] -= flxci;
							con[ixdip][jydip][k] += flxci;
							break;
						}
					} while (phi2[0][ixdi][jydi][k] < phi2[0][ixdip][jydip][k]);
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

				if (phi[0][i][j][k] == 0.0)
				{
					conp[1][i][j][k] = con[i][j][k];
					conp[0][i][j][k] = calC10e(temp[i][j][k]);
				}
				else if (phi[0][i][j][k] > 0.0 && phi[0][i][j][k] < 1.0)
				{
					conp[1][i][j][k] = calC1e(temp[i][j][k]);
					conp[0][i][j][k] = (con[i][j][k] - phi[1][i][j][k] * conp[1][i][j][k]) / phi[0][i][j][k];
				}
				else if (phi[0][i][j][k] == 1.0)
				{
					conp[1][i][j][k] = calC1e(temp[i][j][k]);
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
				if (con[i][j][k] > 1.0)
				{
					con[i][j][k] = 1.0;
				}
				if (con[i][j][k] < 0.0)
				{
					con[i][j][k] = 0.0;
				}
				temp[i][j][k] -= Tr * dtime;
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

				if (phi[0][i][j][k] < 1.0 && phi[0][i][j][k] > 0.0)
				{
					phidxi = (phi[0][ip][j][k] - phi[0][im][j][k]) / 2.0 / dx;
					phidyi = (phi[0][i][jp][k] - phi[0][i][jm][k]) / 2.0 / dx;
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

					} while (phi[0][ixdi][jydi][k] < 1.0);
				}
			}
		}
	}

	// for (i = 0; i <= ndmx; i++)
	// {
	// 	for (j = 0; j <= ndmy; j++)
	// 	{
	// 		for (k = 0; k <= ndmz; k++)
	// 		{
	// 			if (phi[0][i][j][k] < 1.0 && phi[0][i + 1][j][k] == 1.0)
	// 			{
	// 				intpos = i;
	// 				break;
	// 			}
	// 		}
	// 	}
	// }

	// intdis = intpos - NDX / 8;
	// if (intdis > 0)
	// {
	// 	for (i = 0; i <= ndmx; i++)
	// 	{
	// 		for (j = 0; j <= ndmy; j++)
	// 		{
	// 			for (k = 0; k <= ndmz; k++)
	// 			{
	// 				if ((i + intdis) <= ndmx)
	// 				{
	// 					phi[0][i][j][k] = phi[0][i + intdis][j][k];
	// 					phi[1][i][j][k] = phi[1][i + intdis][j][k];
	// 					conp[0][i][j][k] = conp[0][i + intdis][j][k];
	// 					conp[1][i][j][k] = conp[1][i + intdis][j][k];
	// 					con[i][j][k] = con[i + intdis][j][k];
	// 					temp[i][j][k] = temp[i + intdis][j][k];
	// 				}
	// 				else
	// 				{
	// 					phi[0][i][j][k] = 1.0;
	// 					phi[1][i][j][k] = 0.0;
	// 					conp[0][i][j][k] = cl;
	// 					conp[1][i][j][k] = calC1e(temp[i][j][k]);
	// 					con[i][j][k] = cl;
	// 					temp[i][j][k] = temp[i - 1][j][k] + Tg;
	// 				}
	// 			}
	// 		}
	// 	}
	// }

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
	FILE *stream; //ストリームのポインタ設定
	char buffer[30];
	sprintf(buffer, "data/phi/1d%d.csv", step);
	stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		fprintf(stream, "%lf   ", temp[i][0][0]);
		fprintf(stream, "\n");
	}
	fclose(stream); //ファイルをクローズ

	FILE *streamc; //ストリームのポインタ設定
	char bufferc[30];
	sprintf(bufferc, "data/con/1d%d.csv", step);
	streamc = fopen(bufferc, "a"); //書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		fprintf(streamc, "%lf   ", con[i][0][0]);
		fprintf(streamc, "\n");
	}
	fclose(streamc); //ファイルをクローズ
}

void datasave2d(int step)
{
	FILE *stream; //ストリームのポインタ設定
	char buffer[30];
	sprintf(buffer, "data/phi/2d%d.csv", step);
	stream = fopen(buffer, "a"); //書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(stream, "%lf   ", phi[0][i][j][0]);
			fprintf(stream, "\n");
		}
	}
	fclose(stream); //ファイルをクローズ

	FILE *streamc; //ストリームのポインタ設定
	char bufferc[30];
	sprintf(bufferc, "data/con/2d%d.csv", step);
	streamc = fopen(bufferc, "a"); //書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(streamc, "%lf   ", con[i][j][0]);
			fprintf(streamc, "\n");
		}
	}
	fclose(streamc); //ファイルをクローズ
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

double calC10e(double temp00)
{
	return 0.35;
	// return -(Tm - temp00) / mle;
}

double calC1e(double temp00)
{
	return 0.25;
	// return -(Tm - temp00) / mle * ke;
}

// kv = (ke * (1.0 - pow(veli[i][j][k] / vdb, 2.0)) + veli[i][j][k] / vd) / (1.0 - pow(veli[i][j][k] / vdb, 2.0) + veli[i][j][k] / vd);
// c10e = c1e / kv;
// mv = -(Tm - temp0) / c10e;
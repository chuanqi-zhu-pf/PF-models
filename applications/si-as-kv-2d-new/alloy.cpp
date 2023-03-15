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
#define NDY 128
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
double M0;				 // 粒界の易動度
double W0;				 // ペナルティー項の係数
double A0;				 // 勾配エネルギー係数
double F0;				 // 粒界移動の駆動力
double temp0, Tg, Tr;	 // 温度
double sum1, sum2, sum3; // 各種の和の作業変数
double pddtt;			 // フェーズフィールドの時間変化率

double gamma0; // 粒界エネルギ密度
double delta;  // 粒界幅（差分ブロック数にて表現）
double mobi;   // 粒界の易動度
double vm0;	   // モル体積

double pho;
double astre;
double epsilon0;
double termiikk, termjjkk;

double phidx, phidy, phidz;
double phidxx, phidyy, phidzz;

double conp[N][NDX][NDY][NDZ], conp2[N][NDX][NDY][NDZ];
double con[NDX][NDY][NDZ], coni[NDX][NDY][NDZ], temp[NDX][NDY][NDZ];
double sumc, c0, dc0;
double cl, Tm, ce, ml1, kap1;
double dT, S0, Dl, Ds;
double dF;
double kvi, veli;

double flxci;
double phidxi, phidyi, phidzi;
double nxi, nyi, nzi;
int di;
int xdi, ydi, zdi;
int ixdi, jydi, kzdi;
int xdip, ydip, zdip;
int ixdip, jydip, kzdip;

double fes1, fns1, fos1, fws1, fss1, fis1;
double fel, fnl, fol, fwl, fsl, fil;
double cddttl, cddtts1;

double conpii_ip, conpii_jp, conpii_kp, conpii_im, conpii_jm, conpii_km;
double phiii_ip, phiii_im;
double piiddc, dpiidc;
double Dii, miijj;

int intpos, intdis, allL, passdis;
int intposp, istepp;

//----------------------------------------------------------------------------------------------

void datasave1d(int step),
	datasave2d(int step), datasave3d(int step);

double calKv(double veli);
double calC0e(int pi, double tp), calC1e(int pi, double tp);

//******* メインプログラム ******************************************
int main()
{
	nstep = 20001;
	pstep = 1000;

	dtime = 1.0e-10;
	dx = 1.0e-9;

	delta = 6.0 * dx;
	mobi = 2.0e-2;
	gamma0 = 0.8e-7;

	Dl = 1.5e-9;
	Ds = 3.0e-11;

	// Tg = 5.0e7;
	// Tr = 2.5e7;

	temp0 = 1586.5;
	Tm = 1685.0;
	ml1 = -400.0;
	kap1 = 0.3;
	cl = 0.09;

	pho = 0.01;

	A0 = 8.0 * delta * gamma0 / PI / PI; // 勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gamma0 / delta;			 // ペナルティー項の係数[式(4.40)]
	M0 = mobi * PI * PI / (8.0 * delta); // 粒界の易動度[式(4.40)]

	for (ii = 0; ii <= nm; ii++)
	{
		for (jj = 0; jj <= nm; jj++)
		{
			wij[ii][jj] = W0;
			aij[ii][jj] = A0;
			mij[ii][jj] = M0;
			if (ii == jj)
			{
				wij[ii][jj] = 0.0;
				aij[ii][jj] = 0.0;
				mij[ii][jj] = 0.0;
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
				temp[i][j][k] = temp0 + Tg * i * dx;
				if ((i) * (i) + (j - NDY * 2 / 4) * (j - NDY * 2 / 4) < NDX / 32 * NDX / 32)
				{
					phi[1][i][j][k] = 1.0;
					conp[1][i][j][k] = calC1e(1, temp[i][j][k]);
					phi[0][i][j][k] = 0.0;
					conp[0][i][j][k] = calC0e(1, temp[i][j][k]);
				}
				else
				{
					phi[1][i][j][k] = 0.0;
					conp[1][i][j][k] = calC1e(1, temp[i][j][k]);
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

						if (ii == 1 && jj == 0)
						{
							dF = Tm + coni[i][j][k] * ml1 - temp[i][j][k];
						}
						else if (ii == 0 && jj == 1)
						{
							dF = -(Tm + coni[i][j][k] * ml1 - temp[i][j][k]);
						}
						pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * dF * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
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

					if (phi2[0][i][j][k] == 0.0)
					{
						phidxi = (phi[0][ip][j][k] - phi[0][im][j][k]) / 2.0;
						phidyi = (phi[0][i][jp][k] - phi[0][i][jm][k]) / 2.0;
					}
					else
					{
						phidxi = (phi2[0][ip][j][k] - phi2[0][im][j][k]) / 2.0;
						phidyi = (phi2[0][i][jp][k] - phi2[0][i][jm][k]) / 2.0;
					}

					nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi);
					nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi);

					if (phi2[0][i][j][k] < 0.5)
					{
						di = 0;
						do
						{
							di++;
							xdi = int(round(nxi * di));
							ydi = int(round(nyi * di));

							ixdi = i + xdi;

							jydi = j + ydi;
							if (jydi > ndmy)
							{
								jydi = jydi - ndmy - 1;
							}
							if (jydi < 0)
							{
								jydi = ndmy + jydi + 1;
							}

							if (phi2[0][ixdi][jydi][k] >= 0.5)
							{
								veli = (phi2[1][i][j][k] - phi[1][i][j][k]) / dtime / sqrt(phidxi * phidxi + phidyi * phidyi) * dx;
								kvi = calKv(veli);
								// average concentration of solid with trapping solute
								conp[1][i][j][k] = (conp[1][i][j][k] * phi[1][i][j][k] + (phi2[1][i][j][k] - phi[1][i][j][k]) * (kvi * conp[0][i][j][k])) / phi2[1][i][j][k];
								// solute that can not be trapped by newly-formed solid should be ejected outside
								flxci = (phi2[1][i][j][k] - phi[1][i][j][k]) * conp[0][i][j][k] * (1.0 - kvi);
								con[i][j][k] -= flxci;
								con[ixdi][jydi][k] += flxci;
								break;
							}
						} while (phi2[0][ixdi][jydi][k] < 0.5 && di < int(delta / dx / 2));
					}
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
					conp[0][i][j][k] = calC0e(1, temp[i][j][k]);
				}
				else if (phi[0][i][j][k] >= 0.5 && phi[0][i][j][k] < 1.0)
				{
					for (ii = 1; ii <= nm; ii++)
					{
						con[i][j][k] -= conp[ii][i][j][k] * phi[ii][i][j][k];
					}
					conp[0][i][j][k] = con[i][j][k] / phi[0][i][j][k];
				}
				else if (phi[0][i][j][k] == 1.0)
				{
					conp[1][i][j][k] = calC1e(1, temp[i][j][k]);
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

				for (ii = 0; ii <= nm; ii++)
				{
					if (ii == 0)
					{
						Dii = Dl;
					}
					else
					{
						Dii = Ds;
					}

					conpii_ip = conp[ii][ip][j][k];
					conpii_im = conp[ii][im][j][k];
					conpii_jp = conp[ii][i][jp][k];
					conpii_jm = conp[ii][i][jm][k];
					conpii_kp = conp[ii][i][j][kp];
					conpii_km = conp[ii][i][j][km];
					if (phi[ii][ip][j][k] <= pho && phi[ii][i][j][k] > pho)
					{
						conpii_ip = conp[ii][i][j][k];
					}
					if (phi[ii][im][j][k] <= pho && phi[ii][i][j][k] > pho)
					{
						conpii_im = conp[ii][i][j][k];
					}
					if (phi[ii][i][jp][k] <= pho && phi[ii][i][j][k] > pho)
					{
						conpii_jp = conp[ii][i][j][k];
					}
					if (phi[ii][i][jm][k] <= pho && phi[ii][i][j][k] > pho)
					{
						conpii_jm = conp[ii][i][j][k];
					}
					if (phi[ii][i][j][kp] <= pho && phi[ii][i][j][k] > pho)
					{
						conpii_kp = conp[ii][i][j][k];
					}
					if (phi[ii][i][j][km] <= pho && phi[ii][i][j][k] > pho)
					{
						conpii_km = conp[ii][i][j][k];
					}
					if (phi[ii][i][j][k] >= pho)
					{
						piiddc = Dii * phi[ii][i][j][k] * (conpii_ip + conpii_im + conpii_jp + conpii_jm + conpii_kp + conpii_km - 6.0 * conp[ii][i][j][k]) / dx / dx;
						dpiidc = Dii * ((phi[ii][ip][j][k] - phi[ii][im][j][k]) * (conpii_ip - conpii_im) + (phi[ii][i][jp][k] - phi[ii][i][jm][k]) * (conpii_jp - conpii_jm) + (phi[ii][i][j][kp] - phi[ii][i][j][km]) * (conpii_kp - conpii_km)) / 4.0 / dx / dx;
					}
					else
					{
						piiddc = 0.0;
						dpiidc = 0.0;
					}

					if (phi[ii][i][j][k] > 0.0)
					{
						conp2[ii][i][j][k] = conp[ii][i][j][k] + (piiddc + dpiidc) * dtime / phi[ii][i][j][k];
					}
					else
					{
						conp2[ii][i][j][k] = conp[ii][i][j][k];
					}
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

				if (phi[0][i][j][k] < 1.0 && phi[0][i][j][k] > 0.0)
				{
					coni[i][j][k] = conp[0][i][j][k];

					if (phi[0][i][j][k] >= 0.5 && conp[0][i][j][k] < calC0e(1, temp[i][j][k]))
					{
						phidxi = (phi[0][ip][j][k] - phi[0][im][j][k]) / 2.0;
						phidyi = (phi[0][i][jp][k] - phi[0][i][jm][k]) / 2.0;
						phidzi = (phi[0][i][j][kp] - phi[0][i][j][km]) / 2.0;
						nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi + phidzi * phidzi);
						nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi + phidzi * phidzi);
						nzi = phidzi / sqrt(phidxi * phidxi + phidyi * phidyi + phidzi * phidzi);

						di = 0;
						do
						{
							di--;
							xdi = int(round(nxi * di));
							ydi = int(round(nyi * di));
							zdi = int(round(nzi * di));

							ixdi = i + xdi;

							jydi = j + ydi;
							if (jydi > ndmy)
							{
								jydi = jydi - ndmy - 1;
							}
							if (jydi < 0)
							{
								jydi = ndmy + jydi + 1;
							}

							kzdi = k + zdi;
							if (kzdi > ndmz)
							{
								kzdi = kzdi - ndmz - 1;
							}
							if (kzdi < 0)
							{
								kzdi = ndmz + kzdi + 1;
							}

							if (phi[0][ixdi][jydi][kzdi] < 0.5 && phi[0][ixdi][jydi][kzdi] > 0.0)
							{
								coni[i][j][k] = conp[0][ixdi][jydi][kzdi];
								break;
							}

						} while (phi[0][ixdi][jydi][kzdi] >= 0.5 && abs(di) < (int(delta / dx / 2)));
					}
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
				sumc += con[i][j][k];
				temp[i][j][k] -= Tr * dtime;
			}
		}
	}

	// search the tip position
	allL = 1;
	for (i = ndmx; i >= 0; i--)
	{
		if (phi[0][i][NDY / 4][NDZ / 2] < 1.0)
		{
			intpos = i + passdis;
			allL = 0;
			break;
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

	intdis = intpos - passdis - NDX * 2 / 4;
	if (intdis > 0)
	{
		passdis++;
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
						conp[0][i][j][k] = conp[0][i + intdis][j][k];
						conp[1][i][j][k] = conp[1][i + intdis][j][k];
						con[i][j][k] = con[i + intdis][j][k];
						coni[i][j][k] = coni[i + intdis][j][k];
						temp[i][j][k] = temp[i + intdis][j][k];
					}
					else
					{
						phi[0][i][j][k] = 1.0;
						phi[1][i][j][k] = 0.0;
						conp[0][i][j][k] = cl;
						conp[1][i][j][k] = calC1e(1, temp[i][j][k]);
						con[i][j][k] = cl;
						temp[i][j][k] = temp[ndmx - intdis][j][k] + Tg * (i - ndmx + intdis) * dx;
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

	FILE *streaml; // ストリームのポインタ設定
	char bufferl[30];
	sprintf(bufferl, "data/conl/1d%d.csv", step);
	streaml = fopen(bufferl, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		fprintf(streaml, "%lf   ", conp[0][i][0][0]);
		fprintf(streaml, "\n");
	}
	fclose(streaml); // ファイルをクローズ

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

	FILE *streamt; // ストリームのポインタ設定
	char buffert[30];
	sprintf(buffert, "data/temp/1d%d.csv", step);
	streamt = fopen(buffert, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{

		fprintf(streamt, "%lf   ", temp[i][0][0]);
		fprintf(streamt, "\n");
	}
	fclose(streamt); // ファイルをクローズ
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
			fprintf(stream, "%lf   ", phi[0][i][j][0]);
			fprintf(stream, "\n");
		}
	}
	fclose(stream); // ファイルをクローズ

	FILE *streaml; // ストリームのポインタ設定
	char bufferl[30];
	sprintf(bufferl, "data/conl/2d%d.csv", step);
	streaml = fopen(bufferl, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(streaml, "%lf   ", conp[0][i][j][0]);
			fprintf(streaml, "\n");
		}
	}
	fclose(streaml); // ファイルをクローズ

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

	FILE *streamt; // ストリームのポインタ設定
	char buffert[30];
	sprintf(buffert, "data/temp/2d%d.csv", step);
	streamt = fopen(buffert, "a"); // 書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(streamt, "%lf   ", temp[i][j][0]);
			fprintf(streamt, "\n");
		}
	}
	fclose(streamt); // ファイルをクローズ
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

	return kv;
}

double calC0e(int pi, double tp)
{
	double cc;
	if (pi == 1)
	{
		cc = -(Tm - tp) / ml1;
	}
	return cc;
}

double calC1e(int pi, double tp)
{
	double cc;
	if (pi == 1)
	{
		cc = -(Tm - tp) / ml1 * kap1;
	}
	return cc;
}

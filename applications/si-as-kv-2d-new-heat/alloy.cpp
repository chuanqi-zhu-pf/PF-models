// Adapted from Prof. Koyama's MPF code in his textbook
// phase field coupled with concentration field
// Author: Chuanqi Zhu
// Created on: 2022/11/11

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <omp.h>

#define NDX 128
#define NDY 128
#define NDZ 1
#define NTH 8

#define N 2

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int rows = NDY / NTH;
int nm = N - 1;
double PI = 3.141592;
double RR = 8.3145;

double phi[N][NDX][NDY][NDZ], phi2[N][NDX][NDY][NDZ];
int phiIdx[N + 1][NDX][NDY][NDZ], phiNum[NDX][NDY][NDZ];
double aij[N][N], wij[N][N], mij[N][N], fij[N][N];
double anij[N][N], thij[N][N], vpij[N][N], etaij[N][N];
double conp[N][NDX][NDY][NDZ], conp2[N][NDX][NDY][NDZ];
double con[NDX][NDY][NDZ], coni[NDX][NDY][NDZ], temp[NDX][NDY][NDZ], temp2[NDX][NDY][NDZ];

int nstep, pstep;
double dtime, dx;
double M0, W0, A0;
double temp0, Tg, Tr;
double sumc, sump, sumpp;
double gamma0;
double delta;
double mobi;
double cl, Tm, ml1, kap1;
double Dl, Ds;
double pho;
double astre, astrem;
double dphidT;
int intpos, intdis, allL, passdis;
int intposp, istepp;
int ix, iy, iz, ni, nj;

//----------------------------------------------------------------------------------------------

void datasave1d(int step),
	datasave2d(int step), datasave3d(int step);

double calKv(double vi);
double calC0e(int pi, double tp), calC1e(int pi, double tp);

//******* メインプログラム ******************************************
int main()
{
	nstep = 300001;
	pstep = 3000;

	dtime = 0.25e-10;
	dx = 5.0e-10;

	delta = 6.0 * dx;
	mobi = 2.6e-2;
	gamma0 = 0.8e-7;

	Dl = 1.5e-9;
	Ds = 3.0e-11;

	// Tg = 5.0e6;
	// Tr = 2.0e6;

	temp0 = 1592.0;
	Tm = 1685.0;
	ml1 = -400.0;
	kap1 = 0.3;
	cl = 0.09;

	dphidT = 0.0045;

	pho = 0.01;
	astre = 0.02;
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

	// thij[0][1] = PI / 8.0;
	// thij[1][0] = PI / 8.0;

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

	// Initialization of fields
	sumc = 0.0;
	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			for (iz = 0; iz <= ndmz; iz++)
			{
				temp[ix][iy][iz] = temp0 + Tg * ix * dx;
				if ((ix) * (ix) + (iy - NDY * 2 / 4) * (iy - NDY * 2 / 4) < NDX / 32 * NDX / 32)
				// if (ix < NDX / 32 && iy < NDY / 2)
				{
					phi[1][ix][iy][iz] = 1.0;
					conp[1][ix][iy][iz] = calC1e(1, temp[ix][iy][iz]);
					phi[0][ix][iy][iz] = 0.0;
					conp[0][ix][iy][iz] = calC0e(1, temp[ix][iy][iz]);
				}
				// else if (ix < NDX / 25 && iy >= NDY / 2)
				// {
				// 	phi[1][ix][iy][iz] = 1.0;
				// 	conp[1][ix][iy][iz] = calC1e(1, temp[ix][iy][iz]);
				// 	phi[0][ix][iy][iz] = 0.0;
				// 	conp[0][ix][iy][iz] = calC0e(1, temp[ix][iy][iz]);
				// }
				else
				{
					phi[1][ix][iy][iz] = 0.0;
					conp[1][ix][iy][iz] = calC1e(1, temp[ix][iy][iz]);
					phi[0][ix][iy][iz] = 1.0;
					conp[0][ix][iy][iz] = cl;
				}
				con[ix][iy][iz] = 0.0;
				for (ni = 0; ni <= nm; ni++)
				{
					con[ix][iy][iz] += conp[ni][ix][iy][iz] * phi[ni][ix][iy][iz];
				}
				sumc += con[ix][iy][iz];
			}
		}
	}

	std::cout << "####### Computation Start! ########" << std::endl;
	std::cout << "phi stability is: " << mobi * gamma0 * dtime / dx / dx << std::endl;
	std::cout << "diffusion stability is: " << Dl * dtime / dx / dx << std::endl;
	std::cout << "dT stability 1: " << ((Tm - temp0) + ml1 * cl) * delta / PI / gamma0 << std::endl;
	std::cout << "###################################" << std::endl;

#pragma omp parallel num_threads(NTH)
	{
		int start, end, offset, th_id;

		th_id = omp_get_thread_num();
		offset = th_id * rows;
		start = offset;
		end = offset + rows - 1;

		int istep = 0;
		double sum1, pddtt, deltemp;
		int phinum;
		int i, j, k, ip, im, jp, jm, kp, km;
		int ii, jj, kk;
		int n1, n2, n3;
		double epsilon0;
		double termiikk, termjjkk;
		double phidx, phidy, phidz;
		double phidxx, phidyy, phidzz;
		double dF, kvi, veli;

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
		double phidxm, phidym, phidzm, phidxpm, phidypm, phidzpm, phiabsm2;

		double flxci;
		double phidxi, phidyi, phidzi;
		double nxi, nyi, nzi;
		int di;
		int xdi, ydi, zdi;
		int ixdi, jydi, kzdi;
		int xdip, ydip, zdip;
		int ixdip, jydip, kzdip;

		double conpii_ip, conpii_jp, conpii_kp, conpii_im, conpii_jm, conpii_km;
		double piiddc, dpiidc;
		double Dii, miijj;

	start:;

#pragma omp barrier

		if (th_id == 0)
		{
			sumc = 0.0;
			sump = 0.0;
			for (i = 0; i <= ndmx; i++)
			{
				for (j = 0; j <= ndmy; j++)
				{
					for (k = 0; k <= ndmz; k++)
					{
						sumc += con[i][j][k];
						sump += phi[1][i][j][k];
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
						temp[i][j][k] += (sump - sumpp) * dphidT;
					}
				}
			}

			sumpp = sump;

			// search the tip position
			allL = 1;
			for (i = ndmx; i >= 0; i--)
			{
				if (phi[0][i][NDY / 2][NDZ / 2] < 1.0 || phi[0][i][NDY / 4][NDZ / 2] < 1.0)
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

			if (istep % pstep == 0)
			{
				datasave2d(istep);
				std::cout << "the average concentration is: " << sumc / NDX / NDY / NDZ << std::endl;
			}
		}

#pragma omp barrier
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
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
#pragma omp barrier
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
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

									termiikk = term0 + termx + termy + termz;
								}
								else
								{
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

									termjjkk = term0 + termx + termy + termz;
								}
								else
								{
									termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);
								}

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

							miijj = mij[ii][jj];

							phidxm = (phi[ii][ip][j][k] - phi[ii][im][j][k]) / 2.0;
							phidym = (phi[ii][i][jp][k] - phi[ii][i][jm][k]) / 2.0;
							phidzm = (phi[ii][i][j][kp] - phi[ii][i][j][km]) / 2.0;
							phiabsm2 = phidxm * phidxm + phidym * phidym + phidzm * phidzm;
							if (anij[ii][jj] == 1.0 && phiabsm2 != 0.0)
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

								phidxpm = phidxm * xxp + phidym * yxp + phidzm * zxp;
								phidypm = phidxm * xyp + phidym * yyp + phidzm * zyp;
								phidzpm = phidxm * xzp + phidym * yzp + phidzm * zzp;

								miijj = miijj * (1.0 - 3.0 * astrem + 4.0 * astrem * (pow(phidxpm, 4.0) + pow(phidypm, 4.0) + pow(phidzpm, 4.0)) / pow(phiabsm2, 2.0));
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

#pragma omp barrier
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
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
#pragma omp barrier
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
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

#pragma omp barrier
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
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
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
			{
				for (k = 0; k <= ndmz; k++)
				{
					if (phi[0][i][j][k] == 0.0)
					{
						conp[1][i][j][k] = con[i][j][k];
						conp[0][i][j][k] = calC0e(1, temp[i][j][k]);
					}
					else if (phi[0][i][j][k] >= pho && phi[0][i][j][k] < 1.0)
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
#pragma omp barrier
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
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

#pragma omp barrier
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
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
		for (i = 0; i <= ndmx; i++)
		{
			for (j = start; j <= end; j++)
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

	FILE *streaml; // ストリームのポインタ設定
	char bufferl[30];
	sprintf(bufferl, "data/conl/1d%d.csv", step);
	streaml = fopen(bufferl, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		fprintf(streaml, "%lf   ", conp[0][ix][0][0]);
		fprintf(streaml, "\n");
	}
	fclose(streaml); // ファイルをクローズ

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

	FILE *streamt; // ストリームのポインタ設定
	char buffert[30];
	sprintf(buffert, "data/temp/1d%d.csv", step);
	streamt = fopen(buffert, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{

		fprintf(streamt, "%lf   ", temp[ix][0][0]);
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

	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			fprintf(stream, "%lf   ", phi[0][ix][iy][0]);
			fprintf(stream, "\n");
		}
	}
	fclose(stream); // ファイルをクローズ

	FILE *streaml; // ストリームのポインタ設定
	char bufferl[30];
	sprintf(bufferl, "data/conl/2d%d.csv", step);
	streaml = fopen(bufferl, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			fprintf(streaml, "%lf   ", conp[0][ix][iy][0]);
			fprintf(streaml, "\n");
		}
	}
	fclose(streaml); // ファイルをクローズ

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

	FILE *streamt; // ストリームのポインタ設定
	char buffert[30];
	sprintf(buffert, "data/temp/2d%d.csv", step);
	streamt = fopen(buffert, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			fprintf(streamt, "%lf   ", temp[ix][iy][0]);
			fprintf(streamt, "\n");
		}
	}
	fclose(streamt); // ファイルをクローズ
}

double calKv(double vi)
{
	double kv;
	double vd = 0.68;
	double vdb = 2.6;

	if (vi > vdb)
	{
		kv = 0.99;
	}
	else
	{
		kv = (kap1 * (1.0 - pow((vi / vdb), 2.0)) + vi / vd) / (1.0 - pow((vi / vdb), 2.0) + vi / vd);
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

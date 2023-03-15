// Adapted from Prof. Koyama's MPF code in his textbook
// phase field coupled with concentration field
// Author: Chuanqi Zhu
// Created on: 2022/11/11

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <omp.h>

#define NDX 256
#define NDY 256
#define NDZ 1
#define NTH 8

#define N 2

int ndx = NDX;
int ndy = NDY;
int ndz = NDZ;
int ndmx = NDX - 1;
int ndmy = NDY - 1;
int ndmz = NDZ - 1;
int nm = N - 1;
int rows = NDX / NTH;
double PI = 3.141592;
double RR = 8.3145;

double phi[N][NDX][NDY][NDZ], phi2[N][NDX][NDY][NDZ];
double conp[N][NDX][NDY][NDZ], conp2[N][NDX][NDY][NDZ];
double con[NDX][NDY][NDZ], temp[NDX][NDY][NDZ], coni[NDX][NDY][NDZ], vii[NDX][NDY][NDZ];
int phiIdx[N + 1][NDX][NDY][NDZ], phiNum[NDX][NDY][NDZ];
double aij[N][N], wij[N][N], mij[N][N], fij[N][N];
double anij[N][N], thij[N][N], vpij[N][N], etaij[N][N];

int nstep, pstep;
int ni, nj, ix, iy, iz;
double dtime, dx;
double M0, W0, A0;
double suma, sumb, sumc, c0;

double cl, temp0, Tg, Tr;
double Te, ce, ml1, ml2, kap1;
double c1e, c10e, c2e, c20e;
double Dl, Ds;

double gt0, delta, mobi;
double astre, astrem;
double rr;
double mophi, pho;
int mstep;

//----------------------------------------------------------------------------------------------

void datasave1d(int step), datasave2d(int step), datasave3d(int step);
double calCle(int pi, double tp), calCse(int pi, double tp);
double calDT(int pi, int pj, double cl, double tp);
double caltpc(int pi, double cl, double veli);

//******* メインプログラム ******************************************
int main()
{

	nstep = 64001;
	pstep = 8000;

	// thermodynamic and kinetic parameters
	Te = 850.0;
	ce = 0.122;
	kap1 = 0.131;
	ml1 = -680.0;

	mobi = 1.0e-3;
	Dl = 3.0e-9;
	Ds = 3.0e-12;

	dx = 1.25e-7;
	dtime = 0.1 * dx * dx / Dl;

	temp0 = 922.0;
	Tg = 0.0;
	Tr = 0.0;
	cl = 0.01;

	gt0 = 2.4e-7;
	delta = 6.0 * dx;
	mophi = 0.5;
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

	A0 = 8.0 * delta * gt0 / PI / PI;	 // 勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gt0 / delta;				 // ペナルティー項の係数[式(4.40)]
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
				if (ix + iy < NDX / 32)
				// if ((ix - NDX / 2) * (ix - NDX / 2) + (iy - NDY / 2) * (iy - NDY / 2) < NDX / 50 * NDX / 50)
				{
					phi[1][ix][iy][iz] = 1.0;
					conp[1][ix][iy][iz] = calCse(1, temp[ix][iy][iz]);
					phi[0][ix][iy][iz] = 0.0;
					conp[0][ix][iy][iz] = calCle(1, temp[ix][iy][iz]);
				}
				else
				{
					phi[1][ix][iy][iz] = 0.0;
					conp[1][ix][iy][iz] = calCse(1, temp[ix][iy][iz]);
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

	c0 = sumc / NDX / NDY / NDZ;

	std::cout << "####### Computation Start! ########" << std::endl;
	std::cout << "phi stability is: " << mobi * gt0 * dtime / dx / dx << std::endl;
	std::cout << "diffusion stability is: " << Dl * dtime / dx / dx << std::endl;
	std::cout << "dT stability 1: " << ((cl - ce) * ml1 + (Te - temp0)) * delta / PI / gt0 << std::endl;
	std::cout << "dT stability 2: " << ((cl - ce) * ml2 + (Te - temp0)) * delta / PI / gt0 << std::endl;
	std::cout << "###################################" << std::endl;

	FILE *streamc; // ストリームのポインタ設定
	char bufferc[30];
	sprintf(bufferc, "data/ave_con.csv");
	streamc = fopen(bufferc, "a"); // 書き込む先のファイルを追記方式でオープン
	fprintf(streamc, "################# Computation Start! ###############\n");
	fprintf(streamc, "##   phi stability is:                 %lf   ##\n", mobi * gt0 * dtime / dx / dx);
	fprintf(streamc, "##   diffusion stability is:           %lf   ##\n", Dl * dtime / dx / dx);
	fprintf(streamc, "##   dT stability 1:                   %lf   ##\n", ((cl - ce) * ml1 + (Te - temp0)) * delta / PI / gt0);
	fprintf(streamc, "##   dT stability 2:                   %lf   ##\n", ((cl - ce) * ml2 + (Te - temp0)) * delta / PI / gt0);
	fprintf(streamc, "####################################################\n");
	fclose(streamc);

#pragma omp parallel num_threads(NTH)
	{
		int start, end, offset, th_id;

		th_id = omp_get_thread_num();
		offset = th_id * rows;
		start = offset;
		end = offset + rows - 1;

		int istep = 0;
		int i, j, k, ip, im, jp, jm, kp, km;
		int ii, jj, kk;
		int n1, n2, n3, phinum;

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
		double phidxm, phidym, phidzm, phidxpm, phidypm, phidzpm, phiabsm2;
		double pddtt, sum1, dT;

		double flxci;
		double phidxi, phidyi, phidzi;
		double nxi, nyi, nzi;
		int di;
		int xdi, ydi, zdi;
		int ixdi, jydi, kzdi;
		int xdip, ydip, zdip;
		int ixdip, jydip, kzdip;

		double conpii_ip, conpii_jp, conpii_kp, conpii_im, conpii_jm, conpii_km;
		double phiii_ip, phiii_im;
		double fwii, feii;
		double piiddc, dpiidc;
		double Dii, miijj;
		double veli;

	start:;

		if (th_id == 0)
		{
			suma = 0.0;
			sumb = 0.0;
			sumc = 0.0;
			for (i = 0; i <= ndmx; i++)
			{
				for (j = 0; j <= ndmy; j++)
				{
					for (k = 0; k <= ndmz; k++)
					{
						suma += phi[1][i][j][k];
						sumb += phi[2][i][j][k];
						sumc += con[i][j][k];
						temp[i][j][k] -= Tr * dtime;
					}
				}
			}

			if (istep % pstep == 0)
			{
				datasave2d(istep);
				std::cout << "------------------------------------" << std::endl;
				std::cout << "step number:      " << istep << std::endl;
				std::cout << "average concentration is: " << c0 << std::endl;
				std::cout << "current concentration is: " << sumc / NDX / NDY / NDZ << std::endl;
				std::cout << "Phase Fractions: " << std::endl;
				std::cout << " Al          " << suma / NDX / NDY / NDZ * 100 << "%" << std::endl;
				std::cout << " Si          " << sumb / NDX / NDY / NDZ * 100 << "%" << std::endl;
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
						kp = ndmz;
					}
					if (k == 0)
					{
						km = 0;
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
					vii[i][j][k] = 0.0;
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
						kp = ndmz;
					}
					if (k == 0)
					{
						km = 0;
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

								phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
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

							dT = calDT(ii, jj, coni[i][j][k], temp[i][j][k]);

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
							pddtt += -2.0 * miijj / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * dT * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
						}
						phi2[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime;
						if (phi2[ii][i][j][k] >= 1.0)
						{
							phi2[ii][i][j][k] = 1.0;
						}
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
						kp = ndmz;
					}
					if (k == 0)
					{
						km = 0;
					}
					if (phi[0][i][j][k] > 0.0 && phi2[0][i][j][k] < 1.0)
					{
						// calculate the inteface normal vector
						if (phi2[0][i][j][k] == 0.0)
						{
							phidxi = (phi[0][ip][j][k] - phi[0][im][j][k]) / 2.0 / dx;
							phidyi = (phi[0][i][jp][k] - phi[0][i][jm][k]) / 2.0 / dx;
							phidzi = (phi[0][i][j][kp] - phi[0][i][j][km]) / 2.0 / dx;
						}
						else
						{
							phidxi = (phi2[0][ip][j][k] - phi2[0][im][j][k]) / 2.0 / dx;
							phidyi = (phi2[0][i][jp][k] - phi2[0][i][jm][k]) / 2.0 / dx;
							phidzi = (phi2[0][i][j][kp] - phi2[0][i][j][km]) / 2.0 / dx;
						}

						nxi = phidxi / sqrt(phidxi * phidxi + phidyi * phidyi + phidzi * phidzi);
						nyi = phidyi / sqrt(phidxi * phidxi + phidyi * phidyi + phidzi * phidzi);
						nzi = phidzi / sqrt(phidxi * phidxi + phidyi * phidyi + phidzi * phidzi);

						// update the local solid concentration &&  calculate the ejection flux
						flxci = 0.0;
						for (ii = 1; ii <= nm; ii++)
						{
							veli = (phi2[ii][i][j][k] - phi[ii][i][j][k]) / dtime / sqrt(phidxi * phidxi + phidyi * phidyi + phidzi * phidzi);

							vii[i][j][k] = (phi2[ii][i][j][k] - phi[ii][i][j][k]) / dtime / sqrt(phidxi * phidxi + phidyi * phidyi + phidzi * phidzi);

							if (phi2[ii][i][j][k] != 0.0)
							{
								conp[ii][i][j][k] = (conp[ii][i][j][k] * phi[ii][i][j][k] + (phi2[ii][i][j][k] - phi[ii][i][j][k]) * caltpc(ii, conp[0][i][j][k], veli)) / phi2[ii][i][j][k];
							}
							else
							{
								conp[ii][i][j][k] = calCse(ii, temp[i][j][k]);
							}
							flxci += (phi2[ii][i][j][k] - phi[ii][i][j][k]) * (conp[0][i][j][k] - caltpc(ii, conp[0][i][j][k], veli));
						}

						if (phi2[0][i][j][k] < mophi)
						{
							// Search the front liquid and cast out the ejected solute
							di = 0;
							do
							{
								di++;
								xdi = int(round(nxi * di));
								ydi = int(round(nyi * di));
								zdi = int(round(nzi * di));

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

								kzdi = k + zdi;
								if (kzdi > ndmz)
								{
									kzdi = kzdi - ndmz - 1;
								}
								if (kzdi < 0)
								{
									kzdi = ndmz + kzdi + 1;
								}

								if (phi2[0][ixdi][jydi][kzdi] >= mophi)
								{
									con[i][j][k] -= flxci;
									con[ixdi][jydi][kzdi] += flxci;
									break;
								}
							} while (phi2[0][ixdi][jydi][kzdi] < mophi && di <= (int(delta / dx)));
						}
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
						conp[0][i][j][k] = 0.0;
						for (ii = 1; ii <= nm; ii++)
						{
							conp[0][i][j][k] += calCle(ii, temp[i][j][k]) * phi[ii][i][j][k];
						}
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
						kp = ndmz;
					}
					if (k == 0)
					{
						km = 0;
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
						}
						else
						{
							piiddc = 0.0;
						}

						if (phi[ii][i][j][k] >= pho)
						{
							dpiidc = Dii * ((phi[ii][ip][j][k] - phi[ii][im][j][k]) * (conpii_ip - conpii_im) + (phi[ii][i][jp][k] - phi[ii][i][jm][k]) * (conpii_jp - conpii_jm) + (phi[ii][i][j][kp] - phi[ii][i][j][km]) * (conpii_kp - conpii_km)) / 4.0 / dx / dx;
						}
						else
						{
							piiddc = 0.0;
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
		for (i = start; i <= end; i++)
		{
			for (j = 0; j <= ndmy; j++)
			{
				for (k = 0; k <= ndmz; k++)
				{
					con[i][j][k] = 0.0;
					for (ii = 0; ii <= nm; ii++)
					{
						conp[ii][i][j][k] = conp2[ii][i][j][k];
						con[i][j][k] += conp[ii][i][j][k] * phi[ii][i][j][k];
					}
					if (con[i][j][k] < 0.0 && phi[0][i][j][k] == 0.0)
					{
						con[i][j][k] = calCse(1, temp[i][j][k]);
						conp[1][i][j][k] = calCse(1, temp[i][j][k]);
					}
					coni[i][j][k] = 0.0;
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
						kp = ndmz;
					}
					if (k == 0)
					{
						km = 0;
					}

					if (phi[0][i][j][k] < 1.0 && phi[0][i][j][k] > 0.0)
					{
						coni[i][j][k] = conp[0][i][j][k];

						if (phi[0][i][j][k] >= 0.5 &&
							(conp[0][i][j][k] < calCle(1, temp[i][j][k]) && phi[1][i][j][k] > 0.0))
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

								if (phi[0][ixdi][jydi][kzdi] < 0.5)
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

double calCle(int pi, double tp)
{
	double cc;
	if (pi == 1)
	{
		cc = (ce + (tp - Te) / ml1);
	}
	return cc;
}

double calCse(int pi, double tp)
{
	double cc;
	if (pi == 1)
	{
		cc = (ce + (tp - Te) / ml1) * kap1;
	}
	return cc;
}

double calDT(int pi, int pj, double cl, double tp)
{
	double dtt;

	if (pi == 1 && pj == 0)
	{
		dtt = ((cl - ce) * ml1 + Te - tp);
	}
	else if (pi == 0 && pj == 1)
	{
		dtt = -(((cl - ce) * ml1 + Te - tp));
	}
	else
	{
		dtt = 0.0;
	}
	return dtt;
}
double caltpc(int pi, double cl, double vi)
{
	double tpcc;
	double vd = 0.68;
	double vdb = 2.6;
	if (pi == 1)
	{
		// tpcc = kap1 * cl;
		if (vi > vdb)
		{
			vi = 1.0;
		}
		tpcc = (kap1 * (1.0 - pow((vi / vdb), 2.0)) + vi / vd) / (1.0 - pow((vi / vdb), 2.0) + vi / vd) * cl;
	}
	return tpcc;
}

void datasave1d(int step)
{
	FILE *stream; // ストリームのポインタ設定
	char buffer[30];
	sprintf(buffer, "data/phi/1d%d.csv", step);
	stream = fopen(buffer, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		fprintf(stream, "%lf   ", vii[ix][0][0]);
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
	fclose(streamc); // ファイルをクローズ

	FILE *streamcl; // ストリームのポインタ設定
	char buffercl[30];
	sprintf(buffercl, "data/conl/1d%d.csv", step);
	streamcl = fopen(buffercl, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		fprintf(streamcl, "%lf   ", conp[0][ix][0][0]);
		fprintf(streamcl, "\n");
	}
	fclose(streamcl); // ファイルをクローズ

	FILE *streamcs; // ストリームのポインタ設定
	char buffercs[30];
	sprintf(buffercs, "data/cons/1d%d.csv", step);
	streamcs = fopen(buffercs, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		fprintf(streamcs, "%lf   ", conp[1][ix][0][0]);
		fprintf(streamcs, "\n");
	}
	fclose(streamcs); // ファイルをクローズ
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

	FILE *streamcl; // ストリームのポインタ設定
	char buffercl[30];
	sprintf(buffercl, "data/conl/2d%d.csv", step);
	streamcl = fopen(buffercl, "a"); // 書き込む先のファイルを追記方式でオープン

	for (ix = 0; ix <= ndmx; ix++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			fprintf(streamcl, "%lf   ", conp[0][ix][iy][0]);
			fprintf(streamcl, "\n");
		}
	}
	fclose(streamcl); // ファイルをクローズ

	FILE *streamcc; // ストリームのポインタ設定
	char buffercc[30];
	sprintf(buffercc, "data/ave_con.csv");
	streamcc = fopen(buffercc, "a"); // 書き込む先のファイルを追記方式でオープン
	fprintf(streamcc, "step number:      %d\n", step);
	fprintf(streamcc, "average concentration is:          %lf\n", c0);
	fprintf(streamcc, "current concentration is:          %lf\n", sumc / NDX / NDY / NDZ);
	fprintf(streamcc, "Al          %lf\n", suma / NDX / NDY / NDZ * 100);
	fprintf(streamcc, "Si          %lf\n", sumb / NDX / NDY / NDZ * 100);
	fclose(streamcc); // ファイルをクローズ
}

void datasave3d(int step)
{
	FILE *stream;
	char buffer[30];
	sprintf(buffer, "data/con/3d%d.vtk", step);
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

	for (iz = 0; iz <= ndmz; iz++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			for (ix = 0; ix <= ndmx; ix++)
			{
				fprintf(stream, "%e\n", con[ix][iy][iz]);
			}
		}
	}
	fclose(stream);

	FILE *streaml;
	char bufferl[30];
	sprintf(bufferl, "data/conl/3d%d.vtk", step);
	streaml = fopen(bufferl, "a");

	fprintf(streaml, "# vtk DataFile Version 1.0\n");
	fprintf(streaml, "phi_%d.vtk\n", step);
	fprintf(streaml, "ASCII\n");
	fprintf(streaml, "DATASET STRUCTURED_POINTS\n");
	fprintf(streaml, "DIMENSIONS %d %d %d\n", NDX, NDY, NDZ);
	fprintf(streaml, "ORIGIN 0.0 0.0 0.0\n");
	fprintf(streaml, "ASPECT_RATIO 1.0 1.0 1.0\n");
	fprintf(streaml, "\n");
	fprintf(streaml, "POINT_DATA %d\n", NDX * NDY * NDZ);
	fprintf(streaml, "SCALARS scalars float\n");
	fprintf(streaml, "LOOKUP_TABLE default\n");

	for (iz = 0; iz <= ndmz; iz++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			for (ix = 0; ix <= ndmx; ix++)
			{
				fprintf(streaml, "%e\n", conp[0][ix][iy][iz]);
			}
		}
	}
	fclose(streaml);

	FILE *streami;
	char bufferi[30];
	sprintf(bufferi, "data/phi/3d%d.vtk", step);
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

	for (iz = 0; iz <= ndmz; iz++)
	{
		for (iy = 0; iy <= ndmy; iy++)
		{
			for (ix = 0; ix <= ndmx; ix++)
			{
				fprintf(streami, "%e\n", coni[ix][iy][iz]);
			}
		}
	}
	fclose(streami);

	FILE *streamc; // ストリームのポインタ設定
	char bufferc[30];
	sprintf(bufferc, "data/ave_con.csv");
	streamc = fopen(bufferc, "a"); // 書き込む先のファイルを追記方式でオープン
	fprintf(streamc, "------------------------------\n");
	fprintf(streamc, "step number:      %d\n", step);
	fprintf(streamc, "average concentration is:          %lf\n", c0);
	fprintf(streamc, "current concentration is:          %lf\n", sumc / NDX / NDY / NDZ);
	fprintf(streamc, "Al          %lf\n", suma / NDX / NDY / NDZ * 100);
	fprintf(streamc, "Si          %lf\n", sumb / NDX / NDY / NDZ * 100);
	fclose(streamc);
}

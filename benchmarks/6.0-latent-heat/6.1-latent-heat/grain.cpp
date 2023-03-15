// Adapted from Prof. Koyama's MPF code in his textbook
// Author: Chuanqi Zhu
// Created on: 2022/11/11

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NDX 351
#define NDY 351
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
double aij[N][N], wij[N][N], mij[N][N], sij[N][N];

int phinum;
int i, j, k, ip, im, jp, jm, kp, km;
int ii, jj, kk;
int n1, n2, n3;

int istep, nstep, pstep;
double dtime, L, dx;
double M0;				 //粒界の易動度
double W0;				 //ペナルティー項の係数
double A0;				 //勾配エネルギー係数
double S0;				 //粒界移動の駆動力
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

//----------- variables for latent heat-----------

double temp[NDX][NDY][NDZ], temp2[NDX][NDY][NDZ];
double Tm, cndct, speht, rlate, Tini, tddtt, Tr, Tg;
int intpos, intdis;
double dtp, dtt;

//------------------------------------------------

void datasave1d(int step), datasave2d(int step), datasave3d(int step);

//******* メインプログラム ******************************************
int main()
{

	//----------- variables for latent heat-----------
	cndct = 84.01;	   //熱伝導率[W/mK]
	speht = 5.42e+06;  //比熱[J/Km^3]
	rlate = 2.350e+09; //潜熱[J/m^3]
	Tm = 1728.0;	   //融点[K]
	Tini = 1511.2;	   //初期温度[K]
	//------------------------------------------------

	nstep = 100000001;
	pstep = 5000;
	dtime = 5.0;
	temp0 = 1000.0;
	L = 3510.0 * 3.0;
	vm0 = 7.0e-6;
	delta = 3.0;
	mobi = 2.0;

	dx = L / double(NDX) * 1.0e-9;		   //差分プロック１辺の長さ(m)
	gamma0 = 0.37 * vm0 / RR / temp0 / dx; //粒界エネルギ密度（0.5J/m^2）を無次元化
	A0 = 8.0 * delta * gamma0 / PI / PI;   //勾配エネルギー係数[式(4.40)]
	W0 = 4.0 * gamma0 / delta;			   //ペナルティー項の係数[式(4.40)]
	M0 = mobi * PI * PI / (8.0 * delta);   //粒界の易動度[式(4.40)]
	S0 = rlate / Tm / RR / temp0;		   //粒界移動の駆動力

	dtp = dx * dx / (5.0 * M0 * A0);	   //時間きざみ
	dtt = dx * dx / (5.0 * cndct / speht); //時間きざみ
	if (dtp > dtt)
	{
		dtime = dtt;
	}
	else
	{
		dtime = dtp;
	}
	printf("delt= %e \n", dtime);

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
				sij[ii][jj] = S0;
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
				if (i * i + j * j < NDX / 8 * NDX / 8)
				// // if ((i - NDX / 2) * (i - NDX / 2) + (j - NDY / 2) * (j - NDY / 2) + (k - NDZ / 2) * (k - NDZ / 2) < NDX / 8 * NDX / 8)
				{
					phi[1][i][j][k] = 1.0;
					phi[0][i][j][k] = 0.0;
				}
				// else
				// {
				// 	phi[1][i][j][k] = 0.0;
				// 	phi[0][i][j][k] = 1.0;
				// }
				// if (i < NDX / 8 && j < NDY / 8)
				// {
				// 	phi[1][i][j][k] = 1.0;
				// 	phi[0][i][j][k] = 0.0;
				// }
				// else if (i < NDX / 7 && j < NDY * 2 / 8 && j >= NDY / 8)
				// {
				// 	phi[1][i][j][k] = 1.0;
				// 	phi[0][i][j][k] = 0.0;
				// }
				// else if (i < NDX / 8 && j < NDY * 3 / 8 && j >= NDY * 2 / 8)
				// {
				// 	phi[1][i][j][k] = 1.0;
				// 	phi[0][i][j][k] = 0.0;
				// }
				// else if (i < NDX / 7 && j < NDY * 4 / 8 && j >= NDY * 3 / 8)
				// {
				// 	phi[1][i][j][k] = 1.0;
				// 	phi[0][i][j][k] = 0.0;
				// }
				// else if (i < NDX / 8 && j < NDY * 5 / 8 && j >= NDY * 4 / 8)
				// {
				// 	phi[1][i][j][k] = 1.0;
				// 	phi[0][i][j][k] = 0.0;
				// }
				// else if (i < NDX / 7 && j < NDY * 6 / 8 && j >= NDY * 5 / 8)
				// {
				// 	phi[1][i][j][k] = 1.0;
				// 	phi[0][i][j][k] = 0.0;
				// }
				// else if (i < NDX / 8 && j < NDY * 7 / 8 && j >= NDY * 6 / 8)
				// {
				// 	phi[1][i][j][k] = 1.0;
				// 	phi[0][i][j][k] = 0.0;
				// }
				// else if (i < NDX / 7 && j >= NDY * 7 / 8)
				// {
				// 	phi[1][i][j][k] = 1.0;
				// 	phi[0][i][j][k] = 0.0;
				// }
				else
				{
					phi[1][i][j][k] = 0.0;
					phi[0][i][j][k] = 1.0;
				}
				temp[i][j][k] = Tini; // + Tg * i;
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
							phidx = (phi[kk][ip][j][k] - phi[kk][im][j][k]) / 2.0 / dx;
							phidy = (phi[kk][i][jp][k] - phi[kk][i][jm][k]) / 2.0 / dx;
							phidz = (phi[kk][i][j][kp] - phi[kk][i][j][km]) / 2.0 / dx;

							phidxx = (phi[kk][ip][j][k] + phi[kk][im][j][k] - 2.0 * phi[kk][i][j][k]) / dx / dx; //フェーズフィールドの空間２階微分
							phidyy = (phi[kk][i][jp][k] + phi[kk][i][jm][k] - 2.0 * phi[kk][i][j][k]) / dx / dx;
							phidzz = (phi[kk][i][j][kp] + phi[kk][i][j][km] - 2.0 * phi[kk][i][j][k]) / dx / dx;

							termiikk = aij[ii][kk] * (phidxx + phidyy + phidzz);
							termjjkk = aij[jj][kk] * (phidxx + phidyy + phidzz);

							sum1 += 0.5 * (termiikk - termjjkk) + (wij[ii][kk] - wij[jj][kk]) * phi[kk][i][j][k];
						}
						pddtt += -2.0 * mij[ii][jj] / (double)phiNum[i][j][k] * (sum1 - 8.0 / PI * sij[ii][jj] * (Tm - temp[i][j][k]) * sqrt(phi[ii][i][j][k] * phi[jj][i][j][k]));
						//フェーズフィールドの発展方程式[式(4.31)]
					}
					phi2[ii][i][j][k] = phi[ii][i][j][k] + pddtt * dtime; //フェーズフィールドの時間発展（陽解法）

					tddtt = cndct * (temp[ip][j][k] + temp[im][j][k] + temp[i][jp][k] + temp[i][jm][k] - 4.0 * temp[i][j][k]) / dx / dx + 30.0 * phi[1][i][j][k] * (1.0 - phi[1][i][j][k]) * phi[1][i][j][k] * (1.0 - phi[1][i][j][k]) * rlate * pddtt / speht;

					temp2[i][j][k] = temp[i][j][k] + tddtt * dtime;

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
				temp[i][j][k] = temp2[i][j][k];
				temp[i][j][k] -= Tr * dtime;
			}
		}
	}

	// for (i = 0; i <= ndmx; i++)
	// {
	// 	for (j = 0; j <= ndmy; j++)
	// 	{
	// 		for (k = 0; k <= ndmz; k++)
	// 		{
	// 			if (phi[0][i][0][0] < 1.0 && phi[0][i + 1][0][0] == 1.0)
	// 			{
	// 				intpos = i;
	// 				break;
	// 			}
	// 		}
	// 	}
	// }

	// intdis = intpos - NDX / 2;
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
	// 					temp[i][j][k] = temp[i + intdis][j][k];
	// 				}
	// 				else
	// 				{
	// 					phi[0][i][j][k] = 1.0;
	// 					phi[1][i][j][k] = 0.0;
	// 					temp[i][j][k] = temp[i - 1][j][k] + Tg * 2.0;
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
		fprintf(stream, "%lf   ", phi[0][i][0][0]);
		fprintf(stream, "\n");
	}
	fclose(stream); //ファイルをクローズ

	FILE *streamt; //ストリームのポインタ設定
	char buffert[30];
	sprintf(buffert, "data/temp/1d%d.csv", step);
	streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		fprintf(streamt, "%lf   ", temp[i][0][0]);
		fprintf(streamt, "\n");
	}
	fclose(streamt); //ファイルをクローズ
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
			fprintf(stream, "%lf   ", phi[1][i][j][0]);
			fprintf(stream, "\n");
		}
	}
	fclose(stream); //ファイルをクローズ

	FILE *streamt; //ストリームのポインタ設定
	char buffert[30];
	sprintf(buffert, "data/temp/2d%d.csv", step);
	streamt = fopen(buffert, "a"); //書き込む先のファイルを追記方式でオープン

	for (i = 0; i <= ndmx; i++)
	{
		for (j = 0; j <= ndmy; j++)
		{
			fprintf(streamt, "%lf   ", temp[i][j][0]);
			fprintf(streamt, "\n");
		}
	}
	fclose(streamt); //ファイルをクローズ
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

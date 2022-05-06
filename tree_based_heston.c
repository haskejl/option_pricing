#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "./include/math_funcs.h"

//#define _DEBUG_

// Parameters
const double kappa = 5.f;
const double theta = 0.16;
const double omega = 0.9;
const double rho = 0.1;
const double r = 0.1;
const double T = 0.25;
const double V_0 = 0.0625;

double calc_z(const double y, const double v, const double z, const double v_pos,
							const double dt) {
	return z + (r-0.5*v)*dt + y*sqrt(v_pos*dt);
}

double calc_v(const double y, const double v, const double v_pos,
							const double dt) {
	return v + kappa*(theta-v_pos)*dt + y*omega*sqrt(v_pos*dt);
}

double bilin_interp(const double x, const double x0, const double x1, 
										const double y, const double y0, const double y1,
										const double z00, const double z10, 
										const double z01, const double z11) {
	#ifdef _DEBUG_
	if(x > x1 || x < x0) {
		printf("Invalid value of X: x=%f, x0=%f, x1=%f\n", x, x0, x1);
		exit(-1);
	}
	if(y > y1 || y < y0) {
 		printf("Invalid value of Y: y=%f, y0=%f, y1=%f\n", y, y0, y1);
		exit(-1);
	}
	#endif

	const double y_frac = (y-y0)/(y1-y0);
	#ifdef _DEBUG_
	if(y_frac > 1.f || y_frac < 0.f) {
		printf("Invalid y_frac");
		exit(-1);
	}
	#endif
	const double x_frac = (x-x0)/(x1-x0);
	#ifdef _DEBUG_
	if(x_frac > 1.f || x_frac < 0.f) {
		printf("Invalid x_frac");
		exit(-1);
	}
	#endif
	return (1-y_frac)*((1-x_frac)*z00 + x_frac*z10) + y_frac*((1-x_frac)*z01 + x_frac*z11);
}

double european_16(const double S_0, const double E, 
									 const int n, const int mv, const int mz) {
	const double dt = T/n;
	double max_z[n];
	double min_z[n];
	double max_v[n];
	double min_v[n];

	double v = V_0;
	double z = log(S_0);
	
	double v_pos_max = v;
	double v_pos_min = v;

	max_z[0] = min_z[0] = z;
	max_v[0] = min_v[0] = v;

	for(int t=1; t < n; t++) {
		v_pos_max = max(max_v[t-1], 0);
		v_pos_min = max(min_v[t-1], 0);

		max_z[t] = max(max(max(calc_z(-1.f, max_v[t-1], max_z[t-1], v_pos_max, dt),
		 								calc_z(1.f, max_v[t-1], max_z[t-1], v_pos_max, dt) ),
		 								calc_z(-1.f, min_v[t-1], max_z[t-1], v_pos_min, dt)),
		 								calc_z(1.f, min_v[t-1], max_z[t-1], v_pos_min, dt));

		min_z[t] = min(min(min(calc_z(-1.f, max_v[t-1], min_z[t-1], v_pos_max, dt),
		 								calc_z(1.f, max_v[t-1], min_z[t-1], v_pos_max, dt)),
		 								calc_z(-1.f, min_v[t-1], min_z[t-1], v_pos_min, dt)),
		 								calc_z(1.f, min_v[t-1], min_z[t-1], v_pos_min, dt));

		max_v[t] = max(calc_v(-1.f, max_v[t-1], v_pos_max, dt),
		 								calc_v(1.f, max_v[t-1], v_pos_max, dt));

		min_v[t] = min(calc_v(-1.f, min_v[t-1], v_pos_min, dt),
		 								calc_v(1.f, min_v[t-1], v_pos_min, dt));
	}

	double dv[n];
	double dz[n];

	for(int t=0; t < n; t++) {
		dz[t] = (max_z[t]-min_z[t])/mz;
		dv[t] = (max_v[t]-min_v[t])/mv;
	}

	double v_t[n][mv];
	double z_t[n][mz];
	for(int t=0; t<n; t++) {
		for(int i=0; i<mv-1; i++) {
			v_t[t][i] = min_v[t] + dv[t]*i;
		}
		for(int j=0; j<mz-1; j++) {
			z_t[t][j] = min_z[t] + dz[t]*j;
		}
		//avoid floating point error
		v_t[t][mv-1] = max_v[t];
		z_t[t][mz-1] = max_z[t];
	}

	double* St[n][mv];
	double* Stvz = (double*)malloc(n*mv*mz*sizeof(double));
	for(int i=0; i<n; i++) {
		for(int j=0; j<mv; j++)
			St[i][j] = Stvz + i*mv*mz + j*mz;
	}

	//double St[n][mv][mz];
	for(int i=0; i<mv; i++) {
		for(int j=0; j<mz; j++) {
			St[n-1][i][j] = fmax(E-exp(z_t[n-1][j]), 0.f);
			//St[n-1][i][j] = fmax(exp(z_t[n-1][j])-E, 0.f);
		}
	}

	for(int t=n-2; t >= 0; t--) {
		for(int i=0; i<mv; i++) {
			for(int j=0; j<mz; j++) {
				//project forward
				const double v_pos = fmax(v_t[t][i], 0.f);
				const double v0 = calc_v(-1.f, v_t[t][i], v_pos, dt);
				const double v1 = calc_v(1.f, v_t[t][i], v_pos, dt);
				const double z0 = calc_z(-1.f, v_t[t][i], z_t[t][j], v_pos, dt);
				const double z1 = calc_z(1.f, v_t[t][i], z_t[t][j], v_pos, dt);
				#ifdef _DEBUG_
				if(v1 > v_t[t+1][mv-1] || v1 < v_t[t+1][0]) {
					printf("Invalid value of v1\n");
					exit(-1);
				};
				if(v0 > v_t[t+1][mv-1] || v0 < v_t[t+1][0]) {
					printf("Invalid value of v0\n");
					exit(-1);
				};
				if(z1 > z_t[t+1][mz-1] || z1 < z_t[t+1][0]) {
					printf("Invalid value of z1\n");
					exit(-1);
				};
				if(z0 > z_t[t+1][mz-1] || z0 < z_t[t+1][0]) {
					printf("Invalid value of z0\n");
					exit(-1);
				};
				#endif

				//find indices
				const int v0_i = (int)fmin((v0-v_t[t+1][0])/dv[t+1], mv-2);
				const int v1_i = (int)fmin((v1-v_t[t+1][0])/dv[t+1], mv-2);
				const int z0_i = (int)fmin((z0-z_t[t+1][0])/dz[t+1], mz-2);
				const int z1_i = (int)fmin((z1-z_t[t+1][0])/dz[t+1], mz-2);
				
				//interpolate
				const double nn = bilin_interp(v0, v_t[t+1][v0_i], v_t[t+1][v0_i+1], 
											z0, z_t[t+1][z0_i], z_t[t+1][z0_i+1],
											St[t+1][v0_i][z0_i], St[t+1][v0_i+1][z0_i], 
											St[t+1][v0_i][z0_i+1], St[t+1][v0_i+1][z0_i+1]);
				const double np = bilin_interp(v0, v_t[t+1][v0_i], v_t[t+1][v0_i+1], 
											z1, z_t[t+1][z1_i], z_t[t+1][z1_i+1],
											St[t+1][v0_i][z1_i], St[t+1][v0_i+1][z1_i], 
											St[t+1][v0_i][z1_i+1], St[t+1][v0_i+1][z1_i+1]);
				const double pn = bilin_interp(v1, v_t[t+1][v1_i], v_t[t+1][v1_i+1], 
											z0, z_t[t+1][z0_i], z_t[t+1][z0_i+1],
											St[t+1][v1_i][z0_i], St[t+1][v1_i+1][z0_i], 
											St[t+1][v1_i][z0_i+1], St[t+1][v1_i+1][z0_i+1]);
				const double pp = bilin_interp(v1, v_t[t+1][v1_i], v_t[t+1][v1_i+1], 
											z1, z_t[t+1][z1_i], z_t[t+1][z1_i+1],
											St[t+1][v1_i][z1_i], St[t+1][v1_i+1][z1_i], 
											St[t+1][v1_i][z1_i+1], St[t+1][v1_i+1][z1_i+1]);
				St[t][i][j] = 0.25*((1.f+rho)*(nn+pp) + (1.f-rho)*(pn+np));
			  const double z_nn = z0*St[t+1][v0_i][z0_i]/nn;
			 	const double z_np = z1*St[t+1][v0_i][z1_i]/np;
			  const double z_pn = z0*St[t+1][v1_i][z0_i]/pn;
			  const double z_pp = z1*St[t+1][v1_i][z1_i]/pp;
				St[t][i][j] = exp(-r*dt)*max(St[t][i][j], max(E-exp(0.5*(z1+z0)), 0.f));
			}
		}
	}
	const float ret_val = St[0][0][0];
	free(Stvz);
	return ret_val;
}

int main() {
	#ifdef _DEBUG_
	printf("Debug Enabled\n");
	#endif
	const double S0[] = {8.f, 9.f, 10.f, 11.f, 12.f};
	const double E = 10.f;
	const int n[] = {25, 35, 50, 71};
	const int mz[] = {125, 250, 500, 1000};
	const int mv[] = {6, 12, 24, 48};

	double results[5][4];
	double times[5][4];
	for(int i=0; i<5; i++) {
		for(int j=0; j<4; j++) {
			clock_t start_t, end_t;
			start_t = clock();
			results[i][j] = european_16(S0[i], E, n[j], mv[j], mz[j]);
			end_t = clock();
			times[i][j] = (double)(end_t-start_t)/CLOCKS_PER_SEC;
		}
	}
	printf("PUT OPTION\n");
	printf("mz\tmv\tn\t");
	for(int i=0; i<5; i++) {
		printf("%6.1f\t", S0[i]);
	}
	printf("avg time\n");
	for(int i=0; i<4; i++) {
		double avg_time = 0;
		printf("%d\t%d\t%d\t", mz[i], mv[i], n[i]);
		for(int j=0; j<5; j++) {
			avg_time += times[j][i];
			printf("%6.4f\t", results[j][i]);
		}
		printf("%6.4f\n", avg_time/5.f);
	}
}

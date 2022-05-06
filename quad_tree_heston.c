#include <math.h>
#include <stdio.h>
#include <time.h>

#include "./include/math_funcs.h"

// IBM parameters
const double alpha = 5.f;
const double nu = 0.16;
const double psi = 0.9;
const double mu = 0.04588;

const int M = 300;
const double h = 1.f/252.f;
const int n = 1000;
const double T = 0.25;
const double r = 0.10;
const int N = 100;
const double E = 10.f;
double p = 0.135;
const double dt = T/(double)N;
const double X_0 = log(10.f);

//const int x_len = 1329;
const int x_len = 42;
double xvals[] = {10.f, 10.10476689, 10.00654793, 10.01964379,  9.92797276,
        9.9528549, 10.10214772, 10.09690938,  9.89392352, 10.06286014,
       10.12964903,  9.92535359,  9.82189628,  9.82713463,  9.79570456,
        9.81272918,  9.7917758,  9.82844421,  9.80749083,  9.99214248,
       10.09036145, 10.00392876, 10.02488214, 10.00654793, 10.113934,
        9.87558931,  9.69224725,  9.67522263,  9.86118387,  9.78653745,
        9.71712939,  9.77867994,  9.79439497,  9.92797276, 10.13357779,
       10.38501833, 10.3404924, 10.48192771, 10.66657936, 10.79360922,
       10.78837087, 10.71372446};

double calc_phi(const double x) {
	if(x > -1.f && x < 1.f)
		return 1-fabs(x);
	else
		return 0.f;
}

void calc_cdf(double* cdf, const double* X, const double xn, const int n) {
	double C = 0.f;
	for(int i=0; i<n; i++) {
		C += calc_phi(X[i]-xn);
	}
	cdf[0] = calc_phi(X[0]-xn)/C;
	for(int i=1; i<n; i++) {
		cdf[i] = calc_phi(X[i]-xn)/C + cdf[i-1];
	}
}

double calc_sigma(const double x) {
	return log(1.f/sqrt(x));
}

double payoff_func(const double S, const double E, const double r, const double T) {
	return max(E-S, 0)*exp(-r*T);
}

double gen_y_val(const double* cdf, const double* y, const int n) {
	double val = gen_unif_dist_rn(0.f, 1.f);
	int i=0;
	while(i < n) {
		if(cdf[i] == 0.f) i++;
		else break;
	}
	for(int i=0; i<n; i++) {
		if(cdf[i] >= val) {
			return y[i];
		}
	}
	return y[n-1];
}

void find_qs(const double d1, const double d2, const double mul_pt, double* q) {
	if(-d1 < d2) {
		double b = d1/mul_pt;
		double bsq = pow(b, 2);
		q[3] = p;
		q[0] = 0.5*(1.f+b+bsq)-p;
		q[1] = 3.f*p-bsq;
		q[2] = 0.5*(1.f-b+bsq)-3*p;
	}
	else {
		double b = d2/mul_pt;
		double bsq = pow(b, 2);
		q[0] = p;
		q[1] = 0.5*(1.f+b+bsq)-3*p;
		q[2] = 3.f*p-bsq;
		q[3] = 0.5*(1.f-b+bsq)-p;
	}
}

double gen_tree(const double* Y_bar, double** xlist, double** xlist2, 
								double** qlist, double** qlist2, const int list_size) {
	double* curr_x = *xlist2;
	double* curr_q = *qlist2;
	double* next_x = *xlist;
	double* next_q = *qlist;
	static int size_mul = 1;

	double sig = calc_sigma(Y_bar[0]);
	double add_pt = (r-pow(sig, 2.f)/2.f)*dt;
	double mul_pt = sig*sqrt(dt);
	double q[4];

	//Calculate the successors of the first
	int node_j = (int)ceil(X_0/mul_pt);

	double d1 = X_0 - node_j*mul_pt;
	double d2 = X_0 - (node_j-1)*mul_pt;

	find_qs(d1, d2, mul_pt, q);

	for(int i=0; i<4; i++) {
		next_x[i] = (node_j-2.f+i)*mul_pt+add_pt;
		next_q[i] = q[3-i];
	}
	
	int top = 4;
	int j_upp = 4;
	int j_downn = 0;

	for(int i=1; i<N; i++) {
		sig = calc_sigma(Y_bar[i]);
		add_pt = (r-pow(sig, 2.f)/2.f)*dt;
		mul_pt = sig*sqrt(dt);

		top = j_upp-j_downn; //j_upp is over the highest point by 1
		j_upp = ceil(next_x[top-1]/mul_pt)+2; //use next_x for now
		j_downn = ceil(next_x[0]/mul_pt)-2;

		// resize the array if it won't fit
		if(j_upp-j_downn+1 > list_size*size_mul) {
			printf("\nRESIZING...\n");
			while(j_upp-j_downn+1 > list_size*size_mul)
				size_mul++;
			*xlist = realloc(*xlist, sizeof(double)*list_size*size_mul);
			*qlist = realloc(*qlist, sizeof(double)*list_size*size_mul);
			*xlist2 = realloc(*xlist2, sizeof(double)*list_size*size_mul);
			*qlist2 = realloc(*qlist2, sizeof(double)*list_size*size_mul);
			if(*xlist == 0 || *qlist == 0 || *xlist2 == 0 || *qlist2 == 0) {
				printf("ERROR: memory reallocation failed");
				exit(-1);
			}
		}
		if(i % 2 == 0) {
			curr_x = *xlist2;
			curr_q = *qlist2;
			next_x = *xlist;
			next_q = *qlist;
		}
		else {
			curr_x = *xlist;
			curr_q = *qlist;
			next_x = *xlist2;
			next_q = *qlist2;
		} 
		// Clear out the old vals
		for(int j=0; j<j_upp-j_downn; j++) {
			next_q[j] = 0.f;
		}
		// find the new values for x
		for(int j=0; j<j_upp-j_downn; j++) {
			next_x[j] = (j_downn + j)*mul_pt+add_pt;
		}
		
		// find the new probabilities
		for(int j=0; j<top; j++) {
			if(curr_q[j] == 0.f) continue;
			node_j = ceil(curr_x[j]/mul_pt);

			d1 = curr_x[j] - node_j*mul_pt;
			d2 = curr_x[j] - (node_j-1)*mul_pt;

			find_qs(d1, d2, mul_pt, q);

			if(node_j <= j_downn+1) {
				printf("Invalid node value\n");
				exit(-1);
			}

			next_q[node_j+1-j_downn] += curr_q[j]*q[0]; 
			next_q[node_j-j_downn] += curr_q[j]*q[1];
			next_q[node_j-1-j_downn] += curr_q[j]*q[2];
			next_q[node_j-2-j_downn] += curr_q[j]*q[3];
		}
	}
	
	// Use the last calculated values
	curr_x = next_x;
	curr_q = next_q;
	double expected_val = 0;
	for(int j=0; j<j_upp-j_downn; j++) {
		expected_val += payoff_func(exp(curr_x[j]), E, r, T)*curr_q[j];
	}
	return expected_val;
}

int main() {
	clock_t start_prob, end_prob;
	start_prob = clock();
	time_t tim;
	srand((unsigned) time(&tim));

	for(int i=0; i<x_len; i++) {
		xvals[i] = log(xvals[i]);
	}
	
	double x_res[n];
	double y_res[n];
	double y_res2[n];
	double cdf[n];
	// Generate probability distribution
	printf("Calculating cdf...\n");

	// Step 0 - Mutation
	double X;
	double Y;
	for(int k=0; k<n; k++) {
		X = xvals[0];
		Y = log(1/0.25);
		for(int j=0; j<M; j++) {
			const double sig = calc_sigma(Y);
			Y = Y + h/M * alpha * (nu-Y) + sqrt(h/M*Y) * psi * gen_norm_dist_rn(0.f, 1.f) ;
			Y = min(fabs(Y), 1.f);
			X = X + h/M * (mu-pow(sig,2)/2) + sqrt(h/M) * sig * gen_norm_dist_rn(0.f, 1.f);
		}
		x_res[k] = X;
		y_res[k] = Y;
	}

	// Step 0 - Selection
	for(int k=0; k<n; k++) {
		calc_cdf(&cdf[0], &x_res[0], xvals[1], n);
	}
	double* y_curr;
	double* y_past;
	// Steps 1 to K-1
	for(int i=1; i<x_len-1; i++) {
		if(i % 2 == 0) {
			y_past = &y_res2[0];
			y_curr = &y_res[0];
		} else {
			y_past = &y_res[0];
			y_curr = &y_res2[0];
		}
		// Mutation
		for(int k=0; k<n; k++) {
			X = xvals[i];
			Y = gen_y_val(&cdf[0], y_past, n);
			for(int j=0; j<M; j++) {
				const double sig = calc_sigma(Y);
				Y = Y + h/M * alpha * (nu-Y) + sqrt(h/M*Y) * psi * gen_norm_dist_rn(0.f, 1.f) ;
				Y = min(fabs(Y), 1.f);
				X = X + h/M * (mu-pow(sig,2)/2) + sqrt(h/M) * sig * gen_norm_dist_rn(0.f, 1.f);
			}
			x_res[k] = X;
			y_curr[k] = Y;
		}

		// Selection
		for(int k=0; k<n; k++) {
			calc_cdf(&cdf[0], &x_res[0], xvals[i+1], n);
		}
	}

	if(x_len % 2 == 0) {
		y_curr = &y_res[0];
	}
	else {
		y_curr = &y_res2[0];
	}
	end_prob = clock();
	double avg_sigma = cdf[0];
	for(int i=1; i<n; i++) {
		avg_sigma += (cdf[i]-cdf[i-1])*calc_sigma(y_curr[i]);
	}
	printf("Average Sigma = %f\n", avg_sigma);
	printf("Generation time = %f\n", (double)(end_prob-start_prob)/CLOCKS_PER_SEC);
	
	// Generate the tree
	// allocate 10,000 doubles initially (80kb)
	clock_t start_mc, end_mc;
	start_mc = clock();
	const int list_size = 10000;
	
	double* xlist = malloc(sizeof(double)*list_size);
	double* qlist = malloc(sizeof(double)*list_size);
	double* xlist2 = malloc(sizeof(double)*list_size);
	double* qlist2 = malloc(sizeof(double)*list_size);

	if(xlist == 0 || qlist == 0 || xlist2 == 0 || qlist2 == 0) {
		printf("ERROR: memory allocation failed");
		exit(-1);
	}
	double expected_val[17];
	double Y_bar[N];
	int nruns = 100;
	for(int i=0; i<nruns; i++) {
		// Generate the values of Y_bar before making the tree
		for(int j=0; j<N; j++) {
			double rn = gen_unif_dist_rn(0.f, 1.f);
			for(int k=0; k<n; k++) {
				while(cdf[j] == 0)	{
					k++;
					if(k == n) {
						printf("Problem exists with generating probability distribution");
						exit(-1);
					}
				}
				if(rn <= cdf[k]) {
					Y_bar[j] = y_curr[k];
					break;
				}
			}
		}
		printf("\rRun %d/%d", i+1, nruns);
		for(int numer = 0; numer<17; numer++) {
			p = (16.f+numer)/192.f;
			expected_val[numer] += gen_tree(Y_bar, &xlist, &xlist2, &qlist, &qlist2, list_size);
		}
	}
	end_mc = clock();
	printf("\n");
	for(int i = 0; i<17; i++)
		printf("Expected Value: %f\n", expected_val[i]/(double)nruns);
	printf("MC runtime: %f\n", (double)(end_mc-start_mc)/CLOCKS_PER_SEC);

	free(xlist);
	free(qlist);
	free(xlist2);
	free(qlist2);
}

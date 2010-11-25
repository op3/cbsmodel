//    CBS - a program to fit parameters and calculate excitation energies 
//    and transition strengths within the "Confined Beta Soft" rotor model 
//    N. Pietralla and O. M. Gorbachenko, Phys. Rev. C 70, 011304(R) (2004) 
//
//    Copyright (C) 2008  Michael Reese
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

//#include "QuantizationCondition.hpp"
#include "CBSNucleus.hpp"

#include <iostream>
#include <cmath>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_coupling.h>

//#include <boost/math/special_functions/bessel.hpp>



namespace CBS
{


// quantization condition
double Q_nu_r_beta(double nu, double r_beta, double z_M);

// this is a version of Q compatible with all the GSL algorithms
double Q_f(double x, void *params);

// this is the first derivative of Q compatible with all the GSL algorithms
// it is computed numerically, because I don't know the analytical form of this
double Q_df(double x, void *params);

// this fuction provides both, Q and it's 1st derivative
// the GLS algorithms insist on having this because in some
// cases it's possible to optimize things by calculating both
// at the same time.
void Q_fdf(double x, void *params, double *f, double *df);


// This function will find the first N zeros of 
// Q_nu,r_beta(z_M) = J_nu(z_M) Y_nu(r_beta z_M) - J_nu(r_beta z_M) Y_nu(z_M)
// and store them in the vector zeros
bool find_zeros(std::vector<std::vector<double> > &all_zeros, double nu, double r_beta, unsigned s, unsigned L, unsigned max_iterations = 100);




CBSNucleus::CBSNucleus(int Z, int A, double r_beta, double Bbeta_max2, double chi, double beta_max, double e0, double epsil, double a, double bf, double bs)
  : Wu(false), Z_(Z), A_(A), r_b(r_beta), Bbm2(Bbeta_max2), Chi(chi), bM(beta_max), E_0(e0), eps(epsil), a_(a), beta_factor(bf), beta_shift(bs),
	  fix_r_b(true), fix_Bbm2(true), fix_Chi(true), fix_bM(true), fix_E_0(true), fix_eps(true), fix_beta_factor(true), fix_beta_shift(true)
{
}

bool &CBSNucleus::use_Wu()
{
	return Wu;
}

int &CBSNucleus::Z()
{
	return Z_;
}

int &CBSNucleus::A()
{
	return A_;
}

double &CBSNucleus::r_beta() 
{ 
	return r_b; 
}

double &CBSNucleus::Bbeta_max2() 
{
	return Bbm2; 
}

double &CBSNucleus::chi()
{
	return Chi;
}

double &CBSNucleus::beta_max()
{
	return bM;
}

double &CBSNucleus::E0()
{
	return E_0;
}

double &CBSNucleus::epsilon()
{
	return eps;
}

double &CBSNucleus::a()
{
	return a_;
}

double &CBSNucleus::beta_f()
{
	return beta_factor;
}

double &CBSNucleus::beta_s()
{
	return beta_shift;
}

bool &CBSNucleus::fix_r_beta() 
{ 
	return fix_r_b; 
}

bool &CBSNucleus::fix_Bbeta_max2() 
{
	return fix_Bbm2; 
}

bool &CBSNucleus::fix_chi()
{
	return fix_Chi;
}

bool &CBSNucleus::fix_beta_max()
{
	return fix_bM;
}

bool &CBSNucleus::fix_E0()
{
	return fix_E_0;
}

bool &CBSNucleus::fix_epsilon()
{
	return fix_eps;
}

bool &CBSNucleus::fix_a()
{
	return fix_a_;
}

bool &CBSNucleus::fix_beta_f()
{
	return fix_beta_factor;
}

bool &CBSNucleus::fix_beta_s()
{
	return fix_beta_shift;
}

double CBSNucleus::qE0()
{
	double factor = 1e3; // to the the result in units of 10^{-3}
	double X = 3.0 / 4.0 / M_PI * beta_max()*beta_max() * Z();
	return X*X * factor;
}

double CBSNucleus::qE2()
{
	double R0 = 1.22; // fm
	double fm2_barn = 1e-2; // factor to express the result in barn^2 instead of fm^4
	double X = 3.0 * Z() * R0*R0 * exp(2.0/3.0 * log(A())) * beta_max() / 4.0 / M_PI * fm2_barn;
	return X*X;
}

double CBSNucleus::mean_beta(unsigned L, unsigned s)
{
	std::vector< std::vector< double > > &bm(beta_mean[r_beta()]);

	if (bm.size() <= L)
		bm.resize(L+1);
	if (bm[L].size() <= s)
		bm[L].resize(s+1,0.0);
	
	if (bm[L][s] == 0.0)
	{	
		double chi_save = chi();
		chi() = 0;
		bm[L][s] = matrix_element(L,s,L,s,1);
		chi() = chi_save;
	}
	return bm[L][s] * beta_max();
}

double CBSNucleus::mean_beta2(unsigned L, unsigned s)
{
	std::vector< std::vector< double > > &bm(beta_mean2[r_beta()]);

	if (bm.size() <= L)
		bm.resize(L+1);
	if (bm[L].size() <= s)
		bm[L].resize(s+1,0.0);
	
	if (bm[L][s] == 0.0)
	{
		double chi_save = chi();
		chi() = 0;
		bm[L][s] = matrix_element(L,s,L,s,2);
		chi() = chi_save;
	}
	return bm[L][s] * beta_max();
}


double CBSNucleus::nu(double L)
{	

	return sqrt(L*(L+1.0)/3.0 + 9.0/4.0);
}

double CBSNucleus::z(unsigned L, unsigned s)
{
  
	std::vector< std::vector< double > > &zr(zeros[r_beta()]);

//	if (zr.size() <= L)
//		zr.resize(L+1);
	unsigned zrs = zr.size();
	if (zrs > 0) 
		--zrs;
		
	if (zr.size() <= L)
		zr.resize(L+1);
		
	for (unsigned l = zrs; l < zr.size(); ++l)
		find_zeros(zr, nu(l), r_beta(), 2, l);	
	if (zr[L].size() <= s)
		find_zeros(zr, nu(L), r_beta(), s, L);

	return zr[L][s];
}

// calculates level energy
double CBSNucleus::E(unsigned L, unsigned s)
{
	if (L > 400)
		return -1;
	double zero = z(L,s);
	double zero_gs = z(0,0); // ground state
	
	double E = E0() + (zero*zero - zero_gs*zero_gs)/ 2.0 / Bbeta_max2() + epsilon()*L*(L+1);
	if (s == 1)
	{
		E += beta_shift;
		E *= beta_factor;
	}	
	return E;
}

double CBSNucleus::Efull(unsigned L, unsigned s, unsigned n_gamma, int K)
{
	if (L > 400)
		return -1;
//	std::cout << L << "," << s << std::endl;
	double zero = z(L,s);
	double zero_gs = z(0,0); // ground state
	
	double X = 0.5 / (Bbeta_max2()/beta_max()/beta_max());
	double A = 0;
	if (n_gamma != 0)
		A = a() * 3.0 * X / sqrt(mean_beta2(L,s));
	double C = 0;
	if (K != 0)
		C = A / a() / sqrt(mean_beta2(L,s)) / 9.0;
	
	double E = E0() 
			+ (zero*zero - zero_gs*zero_gs)/ 2.0 / Bbeta_max2() 
			+ A * n_gamma
			+ C * K*K
			+ epsilon()*L*(L+1);
	return E;
}

// calculates level energy OR transition strength of a coordinate. 
// This is used by the fit algorithm only
double CBSNucleus::operator()(const Coordinate &c, const double *pbegin, const double *pend)
{
	int i = 0;
	if (!fix_r_beta()) 		r_beta() = std::min(std::max(pbegin[i++], 0.001), 0.999);
	if (!fix_Bbeta_max2()) 	Bbeta_max2() = pbegin[i++];
	if (!fix_beta_max()) 	beta_max() = std::max(pbegin[i++], 0.00001);
	if (!fix_chi()) 		chi() = pbegin[i++];
	if (!fix_E0()) 			E0() = pbegin[i++];
	if (!fix_epsilon()) 	epsilon() = pbegin[i++];
	if (!fix_beta_f()) 		beta_f() = pbegin[i++];
	if (!fix_beta_s()) 		beta_s() = pbegin[i++];

	switch(c.type)
	{
		case Coordinate::E:
			return E(c.L1, c.s1);
		break;
		case Coordinate::BE2:
			return BE2(c.L1, c.s1, c.L2, c.s2);
		break;
		case Coordinate::rho2E0:
			return rho2E0(c.L1, c.s1, c.L2, c.s2);
		break;
		default :
			return 0;
	}	
}

// integrand without any operator. This is used for wave function normalization
double normalization_integrand(double beta, void *params)
{
	double *p = (double*)params;
	
	double nu = p[0];
	double z = p[1];
	double gamma_Y = p[2];
	double beta_max = p[3];
	double c = p[4];

	double psi = c * (gsl_sf_bessel_Jnu(nu, z*beta/beta_max) + gamma_Y * gsl_sf_bessel_Ynu(nu, z*beta/beta_max));
	return beta/beta_max * psi*psi;
}

// returns the normalization factor for the wave function
double CBSNucleus::normalization(unsigned L, unsigned s)
{
	double n = nu(L);
	double z_Ls = z(L, s);
	double gamma_Y =  - gsl_sf_bessel_Jnu(n, r_beta()*z_Ls) / gsl_sf_bessel_Ynu(n, r_beta()*z_Ls);
	double beta_max = 1.0;
	double beta_min = r_beta()*beta_max;
	
	// divide the integration into k parts
	unsigned k = s + (unsigned)(1.0/sqrt(r_beta()));
	double one_div_c_Ls_sqr = 0;
	for (unsigned i = 0; i <= k; ++i)
	{
		// integration borders
		double beta_0 = beta_min + i*(beta_max-beta_min)/(k+1);
		double beta_1 = beta_min + (i+1)*(beta_max-beta_min)/(k+1);
		
		double params[] = {n, z_Ls, gamma_Y, beta_max, 1.0};
		gsl_function f = {normalization_integrand, (void*)params};
		double result, abserr;
		size_t neval;
		gsl_integration_qng(&f, beta_0, beta_1, 0, 1e-6, &result, &abserr, &neval);
		one_div_c_Ls_sqr += result;
	}	
	return sqrt(1.0/one_div_c_Ls_sqr);
}


// writes the wave function into the vector wf
void CBSNucleus::wavefunction(std::vector<double> &wf, unsigned L, unsigned s, unsigned samples)
{
	wf.resize(samples);

	double n = nu(L);
	double z_Ls = z(L, s);
	double gamma_Y =  - gsl_sf_bessel_Jnu(n, r_beta()*z_Ls) / gsl_sf_bessel_Ynu(n, r_beta()*z_Ls);
	double beta_max = 1.0;
	double beta_min = r_beta()*beta_max;
	
	double c_Ls = normalization(L,s);
	
	for (unsigned i = 0; i < samples; ++i)
	{
		double beta = beta_min + i*(beta_max-beta_min)/(samples-1);

		wf[i] = sqrt(beta) * c_Ls * ( gsl_sf_bessel_Jnu(n, z_Ls*beta/beta_max) 
						  + gamma_Y * gsl_sf_bessel_Ynu(n, z_Ls*beta/beta_max) );
	}
}

// transition matrix element integrand with operator: T = (b + chi*b^2)^k 
double matrix_element_integrand(double beta, void *params)
{
	double *p = (double*)params;
	double nu1 = p[0];
	double z1 = p[1];
	double gamma_Y1 = p[2];
	double c1 = p[3];
	
	double nu2 = p[4];
	double z2 = p[5];
	double gamma_Y2 = p[6];
	double c2 = p[7];
	
	double beta_max = p[8];
	
	double chi = p[9];
	
	double k = p[10];
	
	double v = beta/beta_max;
	
	double bra = sqrt(v)/v/v * c1 * (gsl_sf_bessel_Jnu(nu1, z1*v) + gamma_Y1 * gsl_sf_bessel_Ynu(nu1, z1*v));
	double ket = sqrt(v)/v/v * c2 * (gsl_sf_bessel_Jnu(nu2, z2*v) + gamma_Y2 * gsl_sf_bessel_Ynu(nu2, z2*v));
	
	double v2 = v*v;
	double v4 = v2*v2;
	
	return bra * v4 * pow(( v + chi*v2 ),k) * ket;
}

// matrix element of operator: T = (b + chi*b^2)^k
double CBSNucleus::matrix_element(unsigned L1, unsigned s1, unsigned L2, unsigned s2, int k)
{
	double beta_max = 1.0;
	double beta_min = r_beta()*beta_max;

	double n1 = nu(L1);
	double z1 = z(L1, s1);
	double gamma_Y1 =  - gsl_sf_bessel_Jnu(n1, r_beta()*z1) / gsl_sf_bessel_Ynu(n1, r_beta()*z1);
	double c1 = normalization(L1,s1);
	
	double n2 = nu(L2);
	double z2 = z(L2, s2);
	double gamma_Y2 =  - gsl_sf_bessel_Jnu(n2, r_beta()*z2) / gsl_sf_bessel_Ynu(n2, r_beta()*z2);
	double c2 = normalization(L2,s2);

	unsigned intervals = std::max(s1,s2) + (unsigned)(1.0/sqrt(r_beta())); 
	double ME = 0;
	for (unsigned i = 0; i <= intervals; ++i)
	{
		// integration borders
		double beta_0 = beta_min + i*(beta_max-beta_min)/(intervals+1);
		double beta_1 = beta_min + (i+1)*(beta_max-beta_min)/(intervals+1);
		
		double params[] = {n1, z1, gamma_Y1, c1,   n2, z2, gamma_Y2, c2,   beta_max, chi(), k};
		gsl_function f = {matrix_element_integrand, (void*)params};
		double result, abserr;
		size_t neval;
		gsl_integration_qng(&f, beta_0, beta_1, 0, 1e-6, &result, &abserr, &neval);
		ME += result;
	}	
	
	return ME;
}

double clebsch_gordan(double j1, double m1, double j2, double m2, double j3, double m3)
{
	int _j1 = static_cast<int>(2*j1+0.5);
	int _j2 = static_cast<int>(2*j2+0.5);
	int _j3 = static_cast<int>(2*j3+0.5);

	int _m1 = static_cast<int>(2*m1+0.5);
	int _m2 = static_cast<int>(2*m2+0.5);
	int _m3 = static_cast<int>(2*m3+0.5);
	
	return ((((_j1-_j2+_m3)/2)%2)?(-1):(1)) * sqrt(_j3 + 1)
			* gsl_sf_coupling_3j(_j1, _j2, _j3, _m1, _m2, _m3);
}

double CBSNucleus::BE2(unsigned L1, unsigned s1, unsigned L2, unsigned s2)
{
	int k = 1; // for E2 transition
	double ME = matrix_element(L1,s1,L2,s2, k);
	double cl = clebsch_gordan(L1, 0, 2, 0, L2, 0);
	double Wu_factor = use_Wu()?(1.0/(5.94e-6*pow(A(),4.0/3.0))):1;
	return ME*ME * cl*cl* qE2() * Wu_factor;
}

double CBSNucleus::rho2E0(unsigned L1, unsigned s1, unsigned L2, unsigned s2)
{
	if (L1 != L2)
		return 0;
	int k = 2; // for E0 transition
	double ME = matrix_element(L1,s1,L2,s2, k);
	double Wu_factor = use_Wu()?(1.0):1;
	return ME*ME * qE0() * Wu_factor;
}

double CBSNucleus::Qt(unsigned L1, unsigned s1, unsigned L2, unsigned s2)
{
	int k = 1; // for E2 transition
	double ME = matrix_element(L1,s1,L2,s2, k);
	return sqrt(ME*ME * qE2() * 16.0*M_PI/5.0);
}

double CBSNucleus::QtJ_Qt2(unsigned L)
{
	int mode = 2; // for E2 transition
	double MEJ = matrix_element(L,0,L-2,0, mode);
	double ME2 = matrix_element(2,0,0,0, 2);
	return MEJ/ME2;
}




// quantization condition
double Q_nu_r_beta(double nu, double r_beta, double z_M)
{
//	using namespace boost::math;
//	return cyl_bessel_j<long double, long double>(nu, z_M) * cyl_neumann<long double, long double>(nu, r_beta*z_M) -
//		   cyl_bessel_j<long double, long double>(nu, r_beta*z_M) * cyl_neumann<long double, long double>(nu, z_M);
	return gsl_sf_bessel_Jnu(nu, z_M) * gsl_sf_bessel_Ynu(nu, r_beta*z_M) -
		   gsl_sf_bessel_Jnu(nu, r_beta*z_M) * gsl_sf_bessel_Ynu(nu, z_M);
}

// this is a version of Q compatible with all the GSL algorithms
double Q_f(double x, void *params)
{
	double *p = (double*)params;
	double nu = p[0];
	double r_beta = p[1];
	return Q_nu_r_beta(nu, r_beta, x);
}

// this is the first derivative of Q compatible with all the GSL algorithms
// it is computed numerically, because I don't know the analytical form of this
double Q_df(double x, void *params)
{
	const double h = 1e-8;
	const gsl_function f = {Q_f, params};
	double result, abserr;

	gsl_deriv_forward(&f, x, h, &result, &abserr);
	return result;
}

// this fuction provides both, Q and it's 1st derivative
// the GLS algorithms insist on having this because in some
// cases it's possible to optimize things by calculating both
// at the same time.
void Q_fdf(double x, void *params, double *f, double *df)
{
	*f = Q_f(x, params);
	*df = Q_df(x, params);
}


// Zero finding is a tricky thing. But in the special case of the solutions of
// the CBS hamiltonian, the following procedure seems to be successful:
//
//   in the beginning low = 0, high = 0.01. While f(high) doesn't change sign: high -> high*1.1
//                             high2 = 0.01. While f(high+high2) doesn't change sign: high2 -> high2*1.1
//   1st zero: find z_1 by bisection of the interval [low,high]
//   2nd zero: find z_2 by bisection of the interval [high,high2] 
//   ith zero (for i > 2): find z_i by bisection of the interval 
//                             [(3*z_(i-1)-z_(i-2))/2, (5*z_(i-1)-3*z_(i-2))/2] 

// This fuction fills the vector &zeros with the first s zeros of the CBS model quantization condition.
// If the vector &zeros already contains n zeros, only the (s+1)-n remaining zeros will be calculated.
// This only(!) works if the values in &zeros really are zeros of the quantization condition.
// Returns true on success.
bool find_zeros(std::vector< std::vector<double> >  &all_zeros, double nu, double r_beta, unsigned s, unsigned L, unsigned max_iterations)
{
	std::vector<double> & zeros = all_zeros[L];
	double params[2] = {nu, r_beta};
	gsl_function f = {Q_f, params};

//	std::cerr << "nu = " << nu << std::endl;

//double dz = 5e-1;
//for (double z = 1e-1; z < 500; z+=dz)
//std::cout << z << " " << Q_f(z, &params) << " " << Q_df(z, &params) << "\n";

	// finding the zeros
	double low = 0.0001;
	double up = 0.0001; 	
	double up2 = 0.0001;
	
	if (zeros.size() < 2)
	{
		if (L == 0)
		{
			while (Q_f(up, &params)*Q_f(low, &params) > 0) // no change in sign
				up *= 1.01;
	
			while (Q_f(up+up2, &params)*Q_f(up, &params) > 0) // no change in sign
				up2 *= 1.01;
		}
		else
		{
			low = all_zeros[L-1][0];
			up = all_zeros[L-1][1];
			up2 = all_zeros[L-1][2]-up;
		}		
	}		
	
	for (unsigned n = zeros.size(); n <= s; ++n)
	{
		// implementing strategy at top of the page
		double lower, upper;
		switch(n)
		{
			case 0:
				lower = low;
				upper = up;
			break;
			case 1:
				lower = up;
				upper = up+up2;
			break;
			default:
				if (zeros.size() < 2)
					return false; // shouldn't happen
					
				lower = (3.0*zeros[n-1] - zeros[n-2])/2.0;
				upper = (5.0*zeros[n-1] - 3.0*zeros[n-2])/2.0;
			break;
		}
		
		gsl_root_fsolver* s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
		gsl_root_fsolver_set(s, &f, lower, upper);
	
		int status;
		unsigned iter = 0;
		double z0;
		do
		{
			iter++;
			gsl_root_fsolver_iterate(s);
			z0 = gsl_root_fsolver_root(s);
			lower = gsl_root_fsolver_x_lower(s);
			upper = gsl_root_fsolver_x_upper(s);
			status = gsl_root_test_interval(lower, upper, 0, 1e-15);
		}
		while (status == GSL_CONTINUE && iter < max_iterations);

		if (iter == max_iterations)
			std::cerr << "warning: reached maximum number of iterations" << std::endl;

		zeros.push_back(z0);
		
		gsl_root_fsolver_free(s);
		
	}
	return true;
}



} // namespace CBS



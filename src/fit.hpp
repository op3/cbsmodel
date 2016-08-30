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

#ifndef TFIT_HPP
#define TFIT_HPP

#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_deriv.h>


// C coordinate type
template <class C>
struct Dp	// data point
{
	Dp(const C &coordinate = C(), 
	   const double &value = 0, 
	   const double &sigma = 1)
		: c(coordinate), v(value), s(sigma) {}
	C c;		// coordinates
	double v;	// value
	double s;	// sigma
};

// F function type
// Di iterator on container of data points
template <class F, class Di>
struct Fs			// fuction structure
{
	Fs(const F &function,
	   const Di &dbegin,
	   const Di &dend,
	   const double &hstep)
		: f(function), db(dbegin), de(dend), h(hstep)
	{}
	
	F f;
	Di db, de;
	double h; 			// step size during derivation calculation
};

template <class F, class Di>
int fit_f_pp(const gsl_vector *x, void *params, gsl_vector *f)
{
	Fs<F,Di> *fs = static_cast<Fs<F,Di>* >(params);
	const double *pb = gsl_vector_const_ptr(x,0);
	const double *pe = gsl_vector_const_ptr(x,x->size-1); ++pe;


	unsigned int j = 0;
	for (Di i = fs->db; i != fs->de; ++i, ++j)
		gsl_vector_set(f, j, 
			(fs->f(i->c, pb, pe) - i->v) / i->s);

	return GSL_SUCCESS;
}

// helper struct to build a gsl-derivation-calculation compatible function
// F function type
// C coordinate type
template <class F, class C>
struct Ds			// drivation structure
{
	Ds(const F &function,
	   const C *coordinate,
	   const double *pbegin, 
	   const double *pend,
	   const unsigned int n)
		: f(function), c(coordinate), pb(pbegin), pe(pend), n(n)
	{}
	   
	F f;			// function
	const C *c;		// coordintate
	const double *pb, *pe;// parameters
	unsigned int n; // the n-th parameter will be modified
};

// gsl derivation of a function F with respect to a certain parameter
template <class F, class C>
double deriv_f(double x, void *p)
{
	Ds<F,C> *ds = static_cast<Ds<F,C>* >(p);
	double tmp = ds->pb[ds->n];		// save n-th parameter
	((double*)ds->pb)[ds->n] = x;	// overwrite n-th parameter with x
	double v = ds->f(*ds->c, ds->pb, ds->pe); // calculate t				fit_step();
	((double*)ds->pb)[ds->n] = tmp;	// put the n-th parameter back on its place
	return v;
}

template <class F, class C, class Di>
int fit_df_pp(const gsl_vector *x, void *params, gsl_matrix *J)
{
	Fs<F,Di> *fs = static_cast<Fs<F,Di>* >(params);
	const double *pb = gsl_vector_const_ptr(x,0);
	const double *pe = gsl_vector_const_ptr(x,x->size-1);	++pe;
	Ds<F,C> ds(fs->f, 0, pb, pe, 0);
	
	unsigned int m = 0;
	for (Di i = fs->db; i != fs->de; ++i, ++m)
	{
		ds.c = &(i->c);
		gsl_function gslf = {deriv_f<F,C>, &ds};
		
		for (unsigned int n = 0; n < x->size; ++n)
		{
			ds.n = n;
			double result, abserr;
			gsl_deriv_central(&gslf, ds.pb[n], fs->h, &result, &abserr);
			gsl_matrix_set(J, m, n, result/i->s);
		}
	}
	return GSL_SUCCESS;
}

template <class F, class C, class Di>
int fit_fdf_pp(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J)
{
	if (fit_f_pp<F,Di>(x,params,f) == GSL_SUCCESS &&
		fit_df_pp<F,C,Di>(x,params,J) == GSL_SUCCESS)
		return GSL_SUCCESS;
	return GSL_SUCCESS;
}

// F function type
// D container of data points
// P container of parameters
template <class F, class C, class Di>
class Fit
{
	public:
		Fit(F &function,
			const Di &dbegin,
			const Di &dend,
			double *pbegin,
			double *pend,
			bool v = false,
			double hstep = 1e-10)
			  : f(function), db(dbegin), de(dend), pb(pbegin), pe(pend), verbose(v),
			    dsize(0), psize(0), fs(function, db,de, hstep)
		{
			Di i = db;
			while(i++ != de) ++dsize;
			double *j = pb;
			while(j++ != pe) ++psize;
			
			gslf.f = fit_f_pp<F,Di>;
			gslf.df = fit_df_pp<F,C,Di>;
			gslf.fdf = fit_fdf_pp<F,C,Di>;
			gslf.n = dsize;
			gslf.p = psize;
			gslf.params = &fs;
			
			s = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, dsize, psize);
			param = gsl_vector_alloc(psize);
			for (unsigned int i = 0; i < psize; ++i)
				gsl_vector_set(param, i, pb[i]);
				
			gsl_multifit_fdfsolver_set(s, &gslf, param);

			iter = 0;
			result_covar = gsl_matrix_alloc(psize, psize);
			
			if (verbose)
			{
				for (unsigned int i = 0; i < psize; ++i)
					std::cerr << gsl_vector_get(param, i) << " ";
				std::cerr << "(start)" << std::endl;
			}
			
		}
		~Fit()
		{
			gsl_vector_free(param);
			gsl_multifit_fdfsolver_free(s);
			gsl_matrix_free(result_covar);
		}
		
		bool fit_continue(double epsilon)
		{
			return GSL_CONTINUE == gsl_multifit_test_delta(s->dx, s->x, 0.0, epsilon);
		}
		
		int fit_step()
		{
			status = gsl_multifit_fdfsolver_iterate(s);
			++iter;

			if (verbose)
			{
				for (unsigned int i = 0; i < psize; ++i)
					std::cerr << gsl_vector_get(s->x, i) << " ";
				std::cerr << "(" << iter << ")" << std::endl;
			}	
			
			return iter;
		}
		

		int fit(double epsilon = 1e-5, int limit = 100)
		{ 
			do 
			{
				if (limit >= 0 && iter >= limit)
					break;
				fit_step();
			}
			while(fit_continue(epsilon));
			calc_covar();
			calc_chi();
			result_params();
			return status;
		}
		
		void calc_covar()
		{
#if GSL_MAJOR_VERSION >= 2
			gsl_matrix *J = gsl_matrix_alloc(dsize, psize);
			gsl_multifit_fdfsolver_jac(s, J);
			gsl_multifit_covar(J, 0.0, result_covar);
			gsl_matrix_free(J);
#else
			gsl_multifit_covar(s->J, 0.0, result_covar);
#endif		}
		void calc_chi()
		{
			result_chi = gsl_blas_dnrm2(s->f);
		}
		void result_params()
		{
			for (unsigned int i = 0; i < psize; ++i)
				pb[i] = gsl_vector_get(s->x, i);
		}
		void result_errors(double *eb, double *ee)
		{
			for (unsigned int i = 0; i < psize; ++i)
				eb[i] = sqrt(gsl_matrix_get(result_covar,i,i));
		}
		
		double chi()
		{
			return result_chi;
		}
		gsl_matrix covar()
		{
			return result_covar;
		}

	private:
		F f;
		Di db, de;
		double *pb;
		double *pe;

		bool verbose;
		
		unsigned int dsize;
		unsigned int psize;
		
		Fs<F,Di> fs;
		gsl_multifit_function_fdf gslf;
		gsl_multifit_fdfsolver *s;
		
		gsl_vector *param;
		
		int status;
		int iter;
		
		double result_chi;
		gsl_matrix *result_covar;

};

#endif

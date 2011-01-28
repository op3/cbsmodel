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

#ifndef CBSNUCLEUS_HPP
#define CBSNUCLEUS_HPP

#include <vector>
#include <map>
#include <iostream>

namespace CBS
{


// the coordinate structure for the fitting procedure
struct Coordinate
{
	enum Type
	{
		E, 
		BE2,
		rho2E0
	};
	
	Type type;
	unsigned L1, s1, L2, s2;
	Coordinate(Type t, unsigned L1_, unsigned s1_, unsigned L2_, unsigned s2_)
		: type(t), L1(L1_), s1(s1_), L2(L2_), s2(s2_)
	{}
	
	Coordinate (std::istream &in)
		: type(E), L1(0), s1(0), L2(0), s2(0)
	{	
		std::string c;
		in >> c;
		if (c == "E")
		{
			in >> L1 >> s1;
			type = E;
		}
		else if (c == "BE2" || c == "T") // c == "T" for compatibility reasons
		{
			in >> L1 >> s1 >> L2 >> s2;
			type = BE2;
		}
		else if (c == "rho2E0")	
		{
			in >> L1 >> s1 >> L2 >> s2;
			type = rho2E0;
		}
	}
};


class CBSNucleus
{
	public:
		// constructs a CBSNucleus
		CBSNucleus(int Z, int A, double r_beta, double Bbeta_max2, double chi = 0, double qE2 = 1, double e0 = 0, double epsil = 0, double a = 0.52228, double bf = 1, double bs = 0);
		
		// use Weisskopf units as unit for transition strengths
		bool &use_Wu();
		
		// non fitable parameters
		int &Z();
		int &A();
		
		// getter for the current model parameters
		double &r_beta();
		double &Bbeta_max2();
		double &chi();
		double &beta_max();
		double &E0();
		double &epsilon();
		double &a();
		
		double &beta_f();
		double &beta_s();
		
		// to fix or unfix parameters for a fit		
		bool &fix_r_beta();
		bool &fix_Bbeta_max2();
		bool &fix_chi();
		bool &fix_beta_max();
		bool &fix_E0();
		bool &fix_epsilon();
		bool &fix_a();

		bool &fix_beta_f();
		bool &fix_beta_s();
		
		// derived parameters
		double qE2();
		double qE0();
		double mean_beta(unsigned L, unsigned s);
		double mean_beta2(unsigned L, unsigned s);
	
		// returns nu = sqrt(L*(L+1.0)/3.0 + 9.0/4.0)
		double nu(double L);
		
		// returns zero corresponding to angular momentum L in band s
		double z(unsigned L, unsigned s);

		// returns energy corresponding to angular momentum L in band s
		double E(unsigned L, unsigned s);
		
		// returns energy corresponding to all quantum numbers
		double Efull(unsigned L, unsigned s, unsigned n_gamma, int K);
		
		// this gives the fitting procedure access to the nuclear level scheme
		double operator()(const Coordinate &c, const double *pbegin, const double *pend);

		double normalization(unsigned L, unsigned s);

		// writes the normalized wave function into the vector wf
		void wavefunction(std::vector<double> &wf, unsigned L, unsigned s, unsigned samples = 100);
		
		double matrix_element(unsigned L1, unsigned s1, unsigned L2, unsigned s2, int mode);
		
		double transition_strength(unsigned L1, unsigned s1, unsigned L2, unsigned s2);
		
		// reduced transition strength
		double BE2(unsigned L1, unsigned s1, unsigned L2, unsigned s2);

		// reduced matrix element
		double reducedME2(unsigned L1, unsigned s1, unsigned L2, unsigned s2);
		
		// quadrupole moment
		double Qt(unsigned L1, unsigned s1, unsigned L2, unsigned s2);
		
		// quadrupole moment ratio Qt(L)/Qt(2)
		double QtJ_Qt2(unsigned L);
		
		// reduced E0 transition strength
		double rho2E0(unsigned L1, unsigned s1, unsigned L2, unsigned s2);
		
	private:
		bool Wu; // if true, transition strengths are clculated in Weisskopf units
		
		int Z_, A_;
	
		double r_b;
		double Bbm2;
		double Chi;
		double bM;
		double E_0;
		double eps;
		double a_;
		
		// new experimental parameters
		double beta_factor;
		double beta_shift;

		double b_m;

		
		bool fix_r_b;
		bool fix_Bbm2;
		bool fix_Chi;
		bool fix_bM;
		bool fix_E_0;
		bool fix_eps;
		bool fix_a_;
		
		bool fix_beta_factor;
		bool fix_beta_shift;
	
		std::map < double, 
			std::vector< std::vector<double> > > zeros; // zeros[r_beta][L][s] contains the s-th zero for
														// angualar momentum L for nucleus with r_beta.
														// Like this we need to calculate the zeros for one value
														// of r_beta only once even if someone changes r_beta and
														// puts it back to it's old value afterwards (as fitting
														// algorithms usually do).
														
		std::map < double, std::vector< std::vector<double> > > beta_mean;
		std::map < double, std::vector< std::vector<double> > > beta_mean2;												
};

} // namespace CBS

#endif

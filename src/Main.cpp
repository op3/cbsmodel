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

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>

#ifdef HAVE_LIBREADLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_errno.h>

#include "CBSNucleus.hpp"
#include "fit.hpp"



void help(std::ostream &hout)
{
	hout << R"(------------ CBS help  ----------
List of all parameters:
A       : parameter            number of Nucleons
Z       : parameter            number of Protons
rb      : parameter, fitable   stiffness in the interval ]0,1[ (unitless)
Bbm2    : parameter, fitable   B*bmax^2, energy scaling in units of hbar^2 / keV
chi     : parameter, fitable   next order in quadrupole operator (unitless)
bmax    : parameter, fitable   maximal deformation (unitless) 
E0      : parameter, fitable   energy of the ground state in units of keV
eps     : parameter, fitable   admixture of rigid rotor (unitless)
bf      : parameter, fitable   factor for energies of the beta band (unitless)
bs      : parameter, fitable   offset on the beta band energies in units of keV
L0      : parameter            offset for angular momentum
s0      : parameter            offset for beta excitation
maxfitsteps: parameter         maximum number of fit steps

Examples:

Type: 'Z' to see the current value of parameter Z
Type: 'chi 0.5' to set the current value of
      parameter chi to 0.5


IMPORTANT: differently as in the paper
N. Pietralla and O. M. Gorbachenko, Phys. Rev. C 70, 011304(R) (2004)
the quantum number 's' starts counting from zero!!


List of commands:
energy  L s : calculate energy of level (L,s)
Erange  L1 s1  L2 s2 : calculate energies of level (L,s) with 
                            L1 <= L <= L2   and   s1 <= s <= s2
Erange  L s n_gamma K : calculate energy with all quantum numbers
ME2  L1 s1  L2 s2 : caluculate matrix element <L1,s1|E2|L2,s2>
BE2  L1 s1  L2 s2 : calculate B(E2;(L1,s1)->(L2,s2))
rho2E0  L1 s1  L2 s2 : calculate rho^2(E0;(L1,s1)->(L2,s2))
fit  (filename|'begindata DATA 'enddata) {parameter}
                  : fits the given parameters to the
                    data in file 'filename' or the
                    data between the keywords 'begindata'
                    and 'enddata' 
store  filename : store all parameters in file 'filename'
load  fileanme : load all parameters from file 'filename'
Wu : use Weisskopf units for transition strengths
noWu : don't use Weisskopf units for transition strengths
       use e^2 b^2 for B(E2) strength
       and 10^{-3} for rho^2(E0) strength
simpleoutput : use simplified output scheme
fulloutput : use default output scheme
outprecision : precision of number output
wavedat L s filename samples : write wave function data
                               into file 'filename'
waveeps L s filename : create eps file of wavefunction
wavedisp L s : display wavefunction
help : ouptut this help   
info : output program information   
license : output license information
exit : quit the program   

Examples:

Type: 'energy 4 0' to calculate the level
      energy of the (L=4,s=0)-state in 
      units of keV
Type: 'BE2 4 0 2 0' to calculate the transition
       strength from the (L=4,s=0)-state to the
       (L=2,s=0)-state in units of e^2*b^2
Type: 'rho2E0 0 1 0 0' to calculate the transition
       strength from the (L=0,s=1)-state to the
       (L=0,s=0)-state in units of 10^{-3}
Type: 'fit 154Sm.dat rb Bbm2 bmax' to fit the
      parameters rb, Bbm2 and bmax to the data
      in file "154Sm.dat"
Type: 'fit begindata E 2 0 87.0 1.0 E 4 0 286.0 1.0 enddata rb Bbm2'
      to fit the parameters directly to the given datapoints
      without reading datapoints from a file
Type: 'store 154Sm.params' to write the current
      values of all parameters to the file
      "154Sm.params"
Type: 'wavedat 10 0 154Sm.wave 100' to write the
      wavefunction data of the (L=10,s=0)-state
      with 100 samples to the file "154Sm.wave"
Type: 'waveeps 10 0 154Sm.eps' to plot the wavefunction
      of the (L=10,s=0)-state to the file "154Sm.eps"
Type: 'wavedisp 10 1' to display the wavefunction
      of the data of the (L=10,s=0)-state

Example of an inputfile for a fit:
E 2 0       87.7 0.2  #for energy level (L=2,s=0)
E 4 0      286.5 0.3  #for energy level (L=4,s=0)
BE2 2 0 0 0 1.15 0.02 #for B(E2) transition strength
                      #from (L=2,s=0) to (L=0,s=0)


Every command can be given as command line argument, too.

To end the program type 'exit')" << '\n';
}

void info(std::ostream &iout)
{
    iout << R"(CBS - a program to fit parameters and calculate excitation energies
and transition strengths within the "Confined Beta Soft" rotor model 
N. Pietralla and O. M. Gorbachenko, Phys. Rev. C 70, 011304(R) (2004))" << '\n';
}

void license(std::ostream &lout)
{
    lout << R"(CBS - a program to fit parameters and calculate excitation energies
and transition strengths within the \"Confined Beta Soft\" rotor model 
N. Pietralla and O. M. Gorbachenko, Phys. Rev. C 70, 011304(R) (2004) 

Copyright (C) 2008  Michael Reese (email: reese@ikp.tu-darmstadt.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.)" << '\n';
}

// this text will be shown after calling the program
// without any command line arguments
void startup(char *argv0)
{
	std::cout << R"(
CBS Copyright (C) 2008 Michael Reese (email: reese@ikp.tu-darmstadt.de)
This program comes with ABSOLUTELY NO WARRANTY;
This is free software, and you are welcome to redistribute it
under certain conditions; type 'license' for details.

Type 'info' for more information or 'help' for instructions)" << "\n\n";
}

bool simple_output = false;
bool use_weisskopf_units = false;
int outprecision = 7;
int maxfitsteps = 50;
int L0 = 0;
int s0 = 0;

void wavefunction(CBS::CBSNucleus &ncl, unsigned L, unsigned s, std::string filename)
{
	bool disp = false;
	if (filename == "__disp")
		disp = true;
		
// calculate wavefunction
	std::vector<double> wf;
	std::ofstream out("__wavefunction.dat");
	ncl.wavefunction(wf, L, s, 1000);
	for (unsigned i = 0; i < wf.size(); ++ i)
		out << ncl.r_beta()*1.0+(1.0-ncl.r_beta()*1.0)*i/(1000.0-1) << " " << wf[i] << '\n';
		
// gnuplot eps output
	out.close();
	out.open("__wavefunction.gnuplot");
	out << "set grid\n";
	if (disp)
	{
		out << "set terminal png\n";
		out << "set xlabel \"b\"\n";
		out << "set output \"__disp.png\"\n";
		out << "plot["<<1.0*ncl.r_beta()<<":"<<1.0<<"] \"__wavefunction.dat\" using 1:($2**2) w l lt 3 lw 4 title \"b^4*|xi(b)|^2 (L = " << L << ", s = " << s << ")\",\\\n";
		out << "\"__wavefunction.dat\" using 1:2 w l lt 1 lw 4 title \"b^2 * xi(b) (L = " << L << ", s = " << s << ")\", 0 lt -1 lw 4 title \"\"\n";
	}	
	else
	{
		// check for file extension
		int last = filename.size()-1;
		if (last < 3 || filename[last-3] != '.' || filename[last-2] != 'e' || filename[last-1] != 'p' || filename[last] != 's')
			filename += ".eps"; 
	
		out << "set terminal postscript enhanced eps solid color" << '\n';
		out << "set xlabel \"{/Symbol b} / {/Symbol b}_M\" " << '\n';
		out << "set output \"" << filename << '\n';
		out << "plot["<<1.0*ncl.r_beta()<<":"<<1.0<<"] \"__wavefunction.dat\" using 1:($2**2) w l lt 3 lw 4 title \"{/Symbol b}^4|{/Symbol x}({/Symbol b})|^2 (L = " << L << ", s = " << s << ")\",\\" << '\n';
		out << "\"__wavefunction.dat\" using 1:2 w l lt 1 lw 4 title \"{/Symbol b}^2{/Symbol x}({/Symbol b}) (L = " << L << ", s = " << s << ")\", 0 lt 7 lw 4 title \"\"" << '\n';
	}	

//	out << "plot["<<1.0*ncl.r_beta()<<":"<<1.0<<"] \"wavefunction.dat\" using 1:($2*$3) w l lt 3 lw 4 title \"{/Symbol b}^4|{/Symbol x}({/Symbol b})|^2 (L = " << L << ", s = " << s << ")\",\\" << '\n';
//	out << "0 lt 7 lw 4 title \"\"" << '\n';
	out.close();
	system("gnuplot __wavefunction.gnuplot");
	if (disp)
		system("rm __wavefunction.gnuplot __wavefunction.dat && display __disp.png && rm __disp.png");
	else 
		system("rm __wavefunction.gnuplot __wavefunction.dat");
}

void fit(std::istream &in, CBS::CBSNucleus &ncl)
{
	std::vector<Dp<CBS::Coordinate> > data;
	

	for (;;)
	{
		// read lines of the data inputfile
		std::string line;
		getline(in, line);
		
		if (!in)
			break;
			
		std::istringstream lin(line);
//		std::cout << line << '\n';
		
		// check for commented lines
		std::string c;
		lin >> c;

		if (c[0] == '#')
			continue;

		CBS::Coordinate coord(CBS::Coordinate::E, 0,0, 0,0); // coordinate
		if (c == "E")
		{
			coord.type = CBS::Coordinate::E;
			lin >> coord.L1 >> coord.s1;
			coord.L1 += L0;
			coord.s1 += s0;
		}
		if (c == "BE2" || c == "T" )
		{
			coord.type = CBS::Coordinate::BE2;
			lin >> coord.L1 >> coord.s1 >> coord.L2 >> coord.s2;
			coord.L1 += L0;
			coord.s1 += s0;
			coord.L2 += L0;
			coord.s2 += s0;
		}
		if (c == "rho2E0")
		{
			coord.type = CBS::Coordinate::rho2E0;
			lin >> coord.L1 >> coord.s1 >> coord.L2 >> coord.s2;
			coord.L1 += L0;
			coord.s1 += s0;
			coord.L2 += L0;
			coord.s2 += s0;
		}
	
		double v, dv; // value and its error
		lin >> v >> dv;
	
		if (!lin)
			continue;
		
		data.push_back(Dp<CBS::Coordinate>(coord , v, dv));
	}

	// setup vector with nonfixed parameters only
	std::vector<double> params;
	if (!ncl.fix_r_beta()) params.push_back(ncl.r_beta());
	if (!ncl.fix_Bbeta_max2()) params.push_back(ncl.Bbeta_max2());
	if (!ncl.fix_beta_max()) params.push_back(ncl.beta_max());
	if (!ncl.fix_chi()) params.push_back(ncl.chi());
	if (!ncl.fix_E0()) params.push_back(ncl.E0());
	if (!ncl.fix_epsilon()) params.push_back(ncl.epsilon());
	if (!ncl.fix_beta_f()) params.push_back(ncl.beta_f());
	if (!ncl.fix_beta_s()) params.push_back(ncl.beta_s());

	bool abort = false;
	if (data.size()	< 1)
	{
		std::cout << "no data points found" << '\n';
		abort = true;
	}	
	if (params.size() < 1)
	{
		std::cout << "no fit parameters defined" << '\n';
		abort = true;
	}	
	if (data.size() < params.size())
	{
		std::cout << "too few data points in file" << '\n';
		abort = true;
	}	
	if (abort)
	{
		std::cout << '\n';
		std::cout << "expecting: \'fit filename {parameter}'" << '\n';
		std::cout << "where parameter has to be fitable. Type \'help\' for details." << '\n';
		return;
	}	
	bool verbose = true;
	if (simple_output)
		verbose = false;
		
	Fit<CBS::CBSNucleus, CBS::Coordinate, std::vector<Dp<CBS::Coordinate> >::iterator> 
		fit(ncl, data.begin(), data.end(), &(params[0]), &(params[params.size()]), verbose);
		
	int status = fit.fit(1e-6, maxfitsteps); 
	//std::cerr << status << '\n';
	
	std::vector<double> errors(params.size());
	
	fit.result_errors(&(errors[0]), &(errors[errors.size()]));
	
	int i = 0; // points to the current fit parameter
	if (simple_output)
	{
		std::cout << "chi_red^2: " << fit.chi()*fit.chi()/(data.size()-params.size()) << " " ;
		if (!ncl.fix_r_beta()) std::cout << "rbeta: " << (ncl.r_beta() = params[i]) << " " << errors[i] << " ", ++i;
		if (!ncl.fix_Bbeta_max2()) std::cout << "Bbm2: " << (ncl.Bbeta_max2() = params[i]) << " " << errors[i] << " ", ++i;
		if (!ncl.fix_beta_max()) std::cout << "bmax: " << (ncl.beta_max() = params[i]) << " " << errors[i] << " ", ++i;
		if (!ncl.fix_chi()) std::cout << "chi: " << (ncl.chi() = params[i]) << " " << errors[i] << " ", ++i;
		if (!ncl.fix_E0()) std::cout << "E0: " << (ncl.E0() = params[i]) << " " << errors[i] << " ", ++i;
		if (!ncl.fix_epsilon()) std::cout << "eps: " << (ncl.epsilon() = params[i]) << " " << errors[i] << " ", ++i;
		if (!ncl.fix_beta_f()) std::cout << "bf: " << (ncl.beta_f() = params[i]) << " " << errors[i] << " ", ++i;
		if (!ncl.fix_beta_s()) std::cout << "bs: " << (ncl.beta_s() = params[i]) << " " << errors[i] << " ", ++i;
	}
	else
	{
		if (status)
		{
			std::cout << '\n';
			std::cout << "Fit not successful" << '\n';
			return;
		}
		std::cout << '\n';
		std::cout << "Fit successful" << '\n';
		std::cout << '\n';
		if (ncl.fix_r_beta())	std::cout << "rb    : " << ncl.r_beta() << " (fixed) " << '\n';
		else					std::cout << "rb    : " << (ncl.r_beta() = params[i]) << " +- " << errors[i] << '\n', ++i;
			
		if (ncl.fix_Bbeta_max2())	std::cout << "Bbm2  : " << ncl.Bbeta_max2() << " (fixed) " << '\n';
		else						std::cout << "Bbm2  : " << (ncl.Bbeta_max2() = params[i]) << " +- " << errors[i] << '\n', ++i;
				
		if (ncl.fix_beta_max())	std::cout << "bmax  : " << ncl.beta_max() << " (fixed) " << '\n';
		else				std::cout << "bmax  : " << (ncl.beta_max() = params[i]) << " +- " << errors[i] << '\n', ++i;
			
		if (ncl.fix_chi())	std::cout << "chi   : " << ncl.chi() << " (fixed) " << '\n';
		else				std::cout << "chi   : " << (ncl.chi() = params[i]) << " +- " << errors[i] << '\n', ++i;
			
		if (ncl.fix_E0())	std::cout << "E0    : " << ncl.E0() << " (fixed) " << '\n';
		else				std::cout << "E0    : " << (ncl.E0() = params[i]) << " +- " << errors[i] << '\n', ++i;
		
		if (ncl.fix_epsilon())	std::cout << "eps   : " << ncl.epsilon() << " (fixed) " << '\n';
		else					std::cout << "eps   : " << (ncl.epsilon() = params[i]) << " +- " << errors[i] << '\n', ++i;
		if (ncl.fix_beta_f())	std::cout << "bf    : " << ncl.beta_f() << " (fixed) " << '\n';
		else					std::cout << "bf    : " << (ncl.beta_f() = params[i]) << " +- " << errors[i] << '\n', ++i;
		if (ncl.fix_beta_s())	std::cout << "bs    : " << ncl.beta_s() << " (fixed) " << '\n';
		else					std::cout << "bs    : " << (ncl.beta_s() = params[i]) << " +- " << errors[i] << '\n', ++i;
		
		std::cout << '\n';
		std::cout << "reduced chisquare : " << fit.chi()*fit.chi()/(data.size()-params.size()) << '\n';
		std::cout << '\n';
	
		std::cout << std::setw(24) << "" << "   " << std::setw(10) << " CBS fit" << "     " << std::setw(11) << "data points     " << std::setw(10) << "residues" << '\n';
		std::cout << std::setw(24) << "---------------------" << "---" << std::setw(10) << "----------" << "-----" << std::setw(11) << "---------------" << std::setw(10) << "-------------------" << '\n';
		for (unsigned i = 0; i < data.size(); ++i)
		{
			std::ostringstream coord, cbs_value, exp_value, residues;
			if (data[i].c.type == CBS::Coordinate::BE2)
				coord << "B(E2; " << data[i].c.L1 << " " << data[i].c.s1 << " -> " << data[i].c.L2 << " " << data[i].c.s2 << " )", 
				cbs_value << ncl.BE2(data[i].c.L1,data[i].c.s1, data[i].c.L2,data[i].c.s2),
				exp_value << data[i].v << " +- " << data[i].s,
				residues << (ncl.BE2(data[i].c.L1,data[i].c.s1, data[i].c.L2,data[i].c.s2) - data[i].v)/data[i].s;
//				std::cout << data[i].c.L1 << " " << data[i].c.s1 << " -> " << data[i].c.L2 << " " << data[i].c.s2 << " : " << ncl.BE2(data[i].c.L1,data[i].c.s1, data[i].c.L2,data[i].c.s2) << " <=> " << data[i].v << " +- " << data[i].s << '\n';
			else if (data[i].c.type == CBS::Coordinate::rho2E0)
				coord << "rho^2(E0; " << data[i].c.L1 << " " << data[i].c.s1 << " -> " << data[i].c.L2 << " " << data[i].c.s2 << " )", 
				cbs_value << ncl.rho2E0(data[i].c.L1,data[i].c.s1, data[i].c.L2,data[i].c.s2),
				exp_value << data[i].v << " +- " << data[i].s,
				residues << (ncl.rho2E0(data[i].c.L1,data[i].c.s1, data[i].c.L2,data[i].c.s2) - data[i].v)/data[i].s;
//				std::cout << data[i].c.L1 << " " << data[i].c.s1 << " -> " << data[i].c.L2 << " " << data[i].c.s2 << " : " << ncl.rho2E0(data[i].c.L1,data[i].c.s1, data[i].c.L2,data[i].c.s2) << " <=> " << data[i].v << " +- " << data[i].s << '\n';
			else
				coord << "energy( " << data[i].c.L1 << " " << data[i].c.s1 << " )",
				cbs_value << ncl.E(data[i].c.L1,data[i].c.s1),
				exp_value << data[i].v << " +- " << data[i].s,
				residues << (ncl.E(data[i].c.L1,data[i].c.s1) - data[i].v)/data[i].s;
//				std::cout << data[i].c.L1 << " " << data[i].c.s1 << " : " << ncl.E(data[i].c.L1,data[i].c.s1) << " <=> " << data[i].v << " +- " << data[i].s << '\n';

			std::cout << std::setw(24) << coord.str() << " : " 
					  << std::setw(10) << cbs_value.str() << "     " 
					  << std::setw(10) << exp_value.str() << "     "
					  << std::setw(10) << residues.str() << '\n';
		}	
		std::cout << '\n';
	}
}



enum TokenValue
{
	tok_RBETA = 0,
	tok_Z,
	tok_A,
	tok_BBMAX2,
	tok_CHI,
	tok_BMAX,
	tok_E0,
	tok_EPS,
	tok_a,
	tok_BETA_F,
	tok_BETA_S,
	tok_FIT,
	tok_ENERGY,
	tok_ERANGE,
	tok_EFULL,
	tok_EALL,
	tok_ME2,
	tok_BE2,
	tok_BE2RANGE,
	tok_RHO2E0,
	tok_STORE,
	tok_LOAD,
	tok_WAVEDAT,
	tok_WAVEEPS,
	tok_WAVEDISP,
	tok_WU,
	tok_NOWU,
	tok_SIMPLEOUT,
	tok_FULLOUT,
	tok_OUTPRECISION,
	tok_MAXFITSTEPS,
	tok_L0,
	tok_s0,
	
	tok_HELP,
	tok_INFO,
	tok_LICENSE,
	
	tok_EXIT,
	
	tok_LAST // has to be the last in the list
};

std::string TokenStrings[] = 
{
	"rb",
	"Z",
	"A",
	"Bbm2",
	"chi",
	"bmax",
	"E0",
	"eps",
	"a",
	"bf",
	"bs",
	"fit",
	"energy",
	"Erange",
	"Efull",
	"Eall",
	"ME2",
	"BE2",
	"BE2range",
	"rho2E0",
	"store",
	"load",
	"wavedat",
	"waveeps",
	"wavedisp",
	"Wu",
	"noWu",
	"simpleoutput",
	"fulloutput",
	"outprecision",
	"maxfitsteps",
	"L0",
	"s0",
	"help",
	"info",
	"license",
	"exit"
};

TokenValue get_token(std::istringstream &in)
{
	char c;
	std::string token;
	for (int i = 0;;++i)
	{
		if (i == 0)
			in >> c;
		else
			c = in.get();	
		
//		std::cout << "char: " << c << '\n';
		
//		if (!in)
//			return tok_LAST;
		if (!in || c == ' ' || c == '\n')
			break;	
			
		token.push_back(c);
//		std::cout << "progress: " << token << '\n';
	}
//	std::cout << "final: " << token << '\n';
	
	// find the matching token string		
	for (int i = 0; i < (int)tok_LAST; ++i)
		if (token == TokenStrings[i])
			return (TokenValue)i;
	return tok_LAST;		
}

std::string process_command(std::string &command, CBS::CBSNucleus &ncl, int &n_errors, bool &end)
{
	std::istringstream com_in(command); // one command line
	
	double value;
	int intvalue;
	
	std::ostringstream ret_str;
	
	switch(get_token(com_in))
	{
		case tok_RBETA:
			com_in >> value;
			if (com_in)
			{
				if (value < 0.01) value = 0.01;
				if (value > 0.99) value = 0.99;
				ncl.r_beta() = value;
			}	
			ret_str << ncl.r_beta();
		break;
		case tok_Z:
			com_in >> value;
			if (com_in)
			{
				if (value < 0) value = 0;
				ncl.Z() = static_cast<int>(value);
			}	
			ret_str << ncl.Z();
		break;
		case tok_A:
			com_in >> value;
			if (com_in)
			{
				if (value < 0) value = 0;
				ncl.A() = static_cast<int>(value);
			}	
			ret_str << ncl.A();
		break;
		case tok_BBMAX2:
			com_in >> value;
			if (com_in)
			{
				if (value < 1e-8) value = 1e-8;
				ncl.Bbeta_max2() = value;
			}	
			ret_str << ncl.Bbeta_max2();
		break;
		case tok_CHI:
			com_in >> value;
			if (com_in)
				ncl.chi() = value;
			ret_str << ncl.chi();
		break;
		case tok_BMAX:
			com_in >> value;
			if (com_in)
			{
				if (value < 0) value = 0;
				ncl.beta_max() = value;
			}	
			ret_str << ncl.beta_max();
		break;
		case tok_E0:
			com_in >> value;
			if (com_in)
				ncl.E0() = value;
			ret_str << ncl.E0();
		break;
		case tok_EPS:
			com_in >> value;
			if (com_in)
				ncl.epsilon() = value;
			ret_str << ncl.epsilon();
		break;
		case tok_a:
			com_in >> value;
			if (com_in)
				ncl.a() = value;
			ret_str << ncl.a();
		break;
		case tok_BETA_F:
			com_in >> value;
			if (com_in)
				ncl.beta_f() = value;
			ret_str << ncl.beta_f();
		break;
		case tok_BETA_S:
			com_in >> value;
			if (com_in)
				ncl.beta_s() = value;
			ret_str << ncl.beta_s();
		break;
		case tok_FIT:
		{
			// read name of file that contains data
			std::string filename;
			com_in >> filename;
			std::ifstream fin;
			std::istringstream sin;
			std::string buffer;
			bool read_from_file = true;
			bool read_error = false;
			if (filename == "begindata")
			{
				read_from_file = false;
				std::string str;
				do
				{
					com_in >> str;
					if (!com_in)
					{
						std::cout << "expecting \'enddata\' after \'begindata\'" << '\n';
						read_error = true;
						break;
					}
					if (str == "E" || str == "BE2" || str == "rho2E0")
						buffer += "\n" + str + " ";
					else
						buffer += str + " ";	
				}
				while(str != "enddata");
			}
			ncl.fix_r_beta() = true;
			ncl.fix_Bbeta_max2() = true;
			ncl.fix_beta_max() = true;
			ncl.fix_chi() = true;
			ncl.fix_E0() = true;
			ncl.fix_epsilon() = true;
			ncl.fix_beta_f() = true;
			ncl.fix_beta_s() = true;
			for (bool go_on = true; go_on;)
			{
				switch (get_token(com_in))
				{
					case tok_RBETA:
						ncl.fix_r_beta() = false;
					break;
					case tok_BBMAX2:
						ncl.fix_Bbeta_max2() = false;
					break;
					case tok_BMAX:
						ncl.fix_beta_max() = false;
					break;
					case tok_CHI:
						ncl.fix_chi() = false;
					break;
					case tok_E0:
						ncl.fix_E0() = false;
					break;
					case tok_EPS:
						ncl.fix_epsilon() = false;
					break;
					case tok_BETA_F:
						ncl.fix_beta_f() = false;
					break;
					case tok_BETA_S:
						ncl.fix_beta_s() = false;
					break;
					default:
						go_on = false;
					break;
				}
			}

			if (read_from_file)
			{
				fin.open(filename.c_str());
				if (!fin)
					std::cout << "cannot open file \"" << filename << "\"" << '\n';
				fit(fin, ncl);
			}
			else
			{
				sin.str(buffer);
				fit(sin, ncl);
			}
		}	
		break;
		case tok_ENERGY:
		{
			unsigned L, s;
			com_in >> L >> s;
			if (!com_in) 
				std::cout << " expecting: \'energy L s\'" << '\n';
			else	
				std::cout << ncl.E(L,s);
		}
		break;
		case tok_ERANGE:
		{
			unsigned L1, s1;
			unsigned L2, s2;
			com_in >> L1 >> s1 >> L2 >> s2;
			if (!com_in) 
				std::cout << " expecting: \'Erange L1 s1  L2 s2\'" << '\n';
			else	
			{
				for (unsigned L = L1; L <= L2; L += 2)
					for (unsigned s = s1; s <= s2; s += 1)
						std::cout << ncl.E(L,s) << '\n';
			}		
		}
		break;
		case tok_EFULL:
		{
			unsigned L, s, n_gamma;
			int K;
			com_in >> L >> s >> n_gamma >> K ;
			if (!com_in) 
				std::cout << " expecting: \'Efull L s n_gamma K\'" << '\n';
			else	
				std::cout << ncl.Efull(L,s,n_gamma,K);
		}
		break;
		case tok_EALL:
		{
			unsigned L_max = 0;
			unsigned E_max;
			com_in >> E_max;
			if (!com_in) 
				std::cout << " expecting: \'Eall E_max\'" << '\n';
			else	
			{
				for (unsigned ng = 0; ncl.Efull(0,0,ng,0) < E_max; ++ng)
				{
//					std::cerr << "n_gamma = " << ng << '\n';
					for (unsigned s = 0; ncl.Efull(0,s,ng,0) < E_max; ++s)
					{
//						std::cerr << "s = " << s << '\n';
						for (unsigned K = 0 + 2*(ng%2); K <= 2*ng; K += 4)
						{
//							std::cerr << "K = " << K << '\n';
							if (K == 0)
							{
								for (unsigned L = 0; L < 400 && ncl.Efull(L,s,ng,K) < E_max; L += 2)
								{
//									std::cerr << "L = " << L << '\n';
									double E = ncl.Efull(L,s,ng,K);
									if (E < E_max)
									{
										std::cout 
											<< E << " " 
											<< L << " " << s << " " << ng << " " << K << " " 
											<< '\n';
										if (L_max < L)
										{
											L_max = L;
											std::cerr << L_max << '\n';
										}	
									}		
								}		
							}
							else
							{
								for (unsigned L = K; L < 400 && ncl.Efull(L,s,ng,K) < E_max; L += 1)
								{
//									std::cerr << "L = " << L << '\n';
									double E = ncl.Efull(L,s,ng,K);
									if (E < E_max)
									{
										std::cout 
											<< E << " " 
											<< L << " " << s << " " << ng << " " << K << " " 
											<< '\n';
										if (L_max < L)
										{
											L_max = L;
											std::cerr << L_max << '\n';
										}	
									}		
								}		
							}
						}
					}	
				}	
//				std::cerr << "L_max = " << L_max << '\n';	
			}	
		}
		break;
		case tok_ME2:
		{
			unsigned L1, s1, L2, s2;
			com_in >> L1 >> s1 >> L2 >> s2;
			if (!com_in) 
				std::cout << " expecting: \'ME2 L1 s1  L2 s2\'" << '\n';
			else	
				std::cout << ncl.reducedME2(L1,s1, L2,s2);
		}
		break;
		case tok_BE2:
		{
			unsigned L1, s1, L2, s2;
			com_in >> L1 >> s1 >> L2 >> s2;
			if (!com_in) 
				std::cout << " expecting: \'BE2 L1 s1  L2 s2\'" << '\n';
			else	
				std::cout << ncl.BE2(L1,s1, L2,s2);
		}
		break;
		case tok_BE2RANGE:
		{
			unsigned L1, s1, L2, s2, dL, ds;
			com_in >> L1 >> s1 >> L2 >> s2 >> dL >> ds;
			if (!com_in) 
				std::cout << " expecting: \'BE2 L1 s1  L2 s2  dL ds\'" << '\n';
			else	
			{
				for (unsigned L = L1; L <= L2; L += 2)
					for (unsigned s = s1; s <= s2; s += 1)
						std::cout << ncl.BE2(L,s, L-dL,s-ds) << '\n';
			}			
		}
		break;
		case tok_RHO2E0:
		{
			unsigned L1, s1, L2, s2;
			com_in >> L1 >> s1 >> L2 >> s2;
			if (!com_in) 
				std::cout << " expecting: \'rho2E0 L1 s1  L2 s2\'" << '\n';
			else	
				std::cout << ncl.rho2E0(L1,s1, L2,s2);
		}
		break;
		case tok_STORE:
		{
			std::string filename;
			com_in >> filename;
			if (!com_in)
				std::cout << " expecting: \'store filename\'" << '\n';
			else
			{
				std::ofstream store(filename.c_str());
				store << ncl.A() << '\n' 
					  << ncl.Z() << '\n' 
					  << ncl.r_beta() << '\n'
					  << ncl.Bbeta_max2() << '\n'
					  << ncl.beta_max() << '\n'
					  << ncl.chi() << '\n'
					  << ncl.E0() << '\n'
					  << ncl.epsilon() << '\n';
			}
		}
		break;
		case tok_LOAD:
		{
			std::string filename;
			com_in >> filename;
			if (!com_in)
				std::cout << " expecting: \'load filename\'" << '\n';
			else
			{
				std::ifstream load(filename.c_str());
				if (!load)
					std::cout << "cannot open file \"" << filename << "\"" << '\n';
				else 
					load >> ncl.A()
						 >> ncl.Z()
						 >> ncl.r_beta() 
						 >> ncl.Bbeta_max2()
						 >> ncl.beta_max()
						 >> ncl.chi()
						 >> ncl.E0() 
						 >> ncl.epsilon();
			}
		}	
		break;
		case tok_WU:
			ncl.use_Wu() = true;
		break;
		case tok_NOWU:
			ncl.use_Wu() = false;
		break;
		case tok_WAVEDAT:
		{
			int L, s, samples;
			std::string filename;
			com_in >> L >> s >> filename >> samples;
			if (!com_in)
				std::cout << " expecting: \'wavedat L s filename samples\'" << '\n';
			else
			{
				std::vector<double> wfdat;
				ncl.wavefunction(wfdat, L, s, samples);
				
				std::ofstream wav(filename.c_str());
				if (!wav)
					std::cout << "cannot open file \"" << filename << "\"" << '\n';
				else 
					for (int i = 0; i < samples; ++i)
						wav << ncl.r_beta() + i*(1-ncl.r_beta())/(samples-1) << " " << wfdat[i] << '\n';
			}
		}	
		break;
		case tok_WAVEEPS:
		{
			int L, s;
			std::string filename;
			com_in >> L >> s >> filename;
			if (!com_in)
				std::cout << " expecting: \'wavedat L s filename\'" << '\n';
			else
				wavefunction(ncl, L, s, filename);
		}	
		break;
		case tok_WAVEDISP:
		{
			int L, s;
			std::string filename = "__disp";
			com_in >> L >> s;
			if (!com_in)
				std::cout << " expecting: \'wavedat L s \'" << '\n';
			else
				wavefunction(ncl, L, s, filename);
		}	
		break;
		case tok_SIMPLEOUT:
			simple_output = true;
		break;
		case tok_FULLOUT:
			simple_output = false;
		break;
		case tok_OUTPRECISION:
			com_in >> intvalue;
			if (com_in)
			{
				if (intvalue < 1) intvalue = 1;
				outprecision = intvalue;
			}
			ret_str << outprecision;
			std::cout.precision(outprecision);
		break;
		case tok_L0:
			com_in >> intvalue;
			if (com_in)
			  L0 = intvalue;
			ret_str << L0;
		break;
		case tok_s0:
			com_in >> intvalue;
			if (com_in)
			  s0 = intvalue;
			ret_str << s0;
		break;
		case tok_MAXFITSTEPS:
			com_in >> intvalue;
			if (com_in)
			  maxfitsteps = intvalue;
			if (maxfitsteps < 0)
			  maxfitsteps = 0;
			ret_str << maxfitsteps;
		break;
		case tok_HELP:
			help(ret_str);
		break;
		case tok_INFO:
			info(ret_str);
		break;
		case tok_LICENSE:
			license(ret_str);
		break;
		case tok_EXIT:
			end = true;
		break;
		case tok_LAST:
			++n_errors;
			return std::string ("unknown identifier");
		break;
	}

	return ret_str.str();
}


#ifdef HAVE_LIBREADLINE
// A static variable for holding the line
static char *line_read = (char*)0;

// Read a string, and return a pointer to it.
//   Returns NULL on EOF. 
char *rl_gets()
{
	// If the buffer has already been allocated,
	//  return the memory to the free pool. 
	if (line_read)
	{
		free (line_read);
		line_read = (char*)0;
	}

	// Get a line from the user. 
	line_read = readline (">>> ");

	// If the line has any text in it,
	//   save it on the history. 
	if (line_read && *line_read)
		add_history (line_read);

	return (line_read);
}
#endif

void add_command_to_history(std::string &cmd)
{
#ifdef HAVE_LIBREADLINE
	if (line_read)
	{
		free(line_read);
		line_read = (char*)0;
	}
	line_read = (char*)malloc(cmd.size()+1);
	for (unsigned i = 0; i < cmd.size(); ++i) line_read[i] = cmd[i];

	if (line_read && *line_read)
		add_history (line_read);
#endif
}

int main(int argc, char *argv[])
{
	// output license
	if (argc == 1)
		startup(argv[0]);

	gsl_set_error_handler_off();

	// initial parameters for the CBS nucleus
	const int Z = 70; 
	const int A = 168;
	const double r_beta_init = 0.1;
	const double Bbeta_max2_init = 0.01;
	const double chi_init = 0;
	const double bmax_init = 0.3;
	const double E0_init = 0;
	const double eps_init = 0;

	// our CBS nucleus
	CBS::CBSNucleus ncl(Z, A, r_beta_init, Bbeta_max2_init, chi_init, bmax_init,
						E0_init, eps_init);

	// number of errors will be the return value of the program
	int n_errors = 0;
	bool end = false;

	// process command line arguments
	for (int i = 1; i < argc; ++i)
	{
		std::string arg(argv[i]);
		if (arg == "-h" || arg == "--help" || arg == "help")
		{
			std::string command = "help";
			add_command_to_history(command);
			std::cout << process_command(command, ncl, n_errors, end) << '\n';	
		}
		
		if (arg == "simpleoutput" || arg == "fulloutput" || arg == "exit" ||
		    arg == "help" || arg == "license" || arg == "info" || arg == "Wu" || arg == "noWu")
		{
			std::string command = arg;
			add_command_to_history(command);
			std::cout << process_command(command, ncl, n_errors, end);
			if (arg != "exit" && !simple_output)
				std::cout << '\n';
		}
		else if (arg == "Z" || arg == "A" || 
			arg == "rb" || arg == "Bbm2" || arg == "chi" || 
			arg == "bmax" || arg == "E0" || arg == "eps" ||
			arg == "a" ||
			arg == "bf" || arg == "bs" ||
			arg == "outprecision" || arg == "maxfitsteps" || arg == "L0" || arg == "s0")
		{
			bool do_output = true;
			// commands with 1 or 0 arguments (numerical)
			std::string command = arg;
			if (i+1 < argc)
			{
				std::istringstream test_dbl(argv[i+1]); // test for numerical argument
				double x;
				test_dbl >> x;
				if (!!test_dbl)
				{
					command += " ";
					command += argv[++i];
					if (simple_output)
						do_output = false;
				}
			}
			add_command_to_history(command);
			std::string output = process_command(command, ncl, n_errors, end);
			if (do_output)
			{
				std::cout << output;
				std::cout << '\n';
			}
			
			if (end)
				return n_errors;
		}
		else if (arg == "load" || arg == "store" || arg == "Eall") 
		{
			// commands with 1 argument
			std::string command = arg;
			command += " ";
			if (i+1 < argc)
				command += argv[++i];
			add_command_to_history(command);
			std::cout << process_command(command, ncl, n_errors, end);
			if (!simple_output)
				std::cout << '\n';
			
			if (end)
				return n_errors;
		}
		else if (arg == "energy" || arg == "wavedisp")
		{
			// commands with 2 arguments
			std::string command = arg;
			if (i+2 < argc)
			{
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
			}
			add_command_to_history(command);
			std::cout << process_command(command, ncl, n_errors, end);
			std::cout << '\n';
		}
		else if (arg == "waveeps")
		{
			// commands with 3 arguments
			std::string command = arg;
			if (i+3 < argc)
			{
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
			}	
			add_command_to_history(command);
			std::cout << process_command(command, ncl, n_errors, end);
			std::cout << '\n';
		}
		else if (arg == "BE2" || arg == "rho2E0" || arg == "ME2" || arg == "wavedat" || 
				 arg == "Erange" || arg == "Efull")
		{
			// commands with 4 arguments
			std::string command = arg;
			if (i+4 < argc)
			{
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
			}
			add_command_to_history(command);
			std::cout << process_command(command, ncl, n_errors, end);
			std::cout << '\n';
		}
		else if (arg == "BE2range")
		{
			// commands with 6 arguments
			std::string command = arg;
			if (i+6 < argc)
			{
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
				command += " "; command += argv[++i];
			}	
			add_command_to_history(command);
			std::cout << process_command(command, ncl, n_errors, end);
			std::cout << '\n';
		}
		else if (arg == "fit")
		{
			// this command has a variable number of arguments 
			// and has to be treated specially 
			bool data_mode = false;
			int i0 = i+1;
			std::string command = arg;
			while (++i < argc)
			{
				std::string par_check(argv[i]);
				if (par_check == "rb" || par_check == "Bbm2" || par_check == "chi" || 
					par_check == "bmax" || par_check == "E0" || par_check == "eps" ||
					par_check == "bf" || par_check == "bs" ||
					i == i0 || data_mode) // for the filename
				{	
					if (i == i0 && data_mode == false) // for direct given datapoints
					{
						if (std::string(argv[i]) == "begindata")
							data_mode = true;
					}
					else if (data_mode)
					{
						if (std::string(argv[i]) == "enddata")
							data_mode = false;
					}
					command += " ";
					command += argv[i];
				}
				else
				{
					--i;
					break;
				}
			}
			
			add_command_to_history(command);
			std::cout << process_command(command, ncl, n_errors, end);
			std::cout << '\n';
		}
		else
		{
			std::cerr << "unknown command line option: " << argv[i] << '\n';
			
			++n_errors;
			return n_errors;
		}
	}
	
	// main loop for command input
	try
	{
		while (!end)
		{
#ifdef HAVE_LIBREADLINE
			// version with GNU readline
			char *line_rl = rl_gets();
			std::string command(line_rl);
#else
			// version without GNU readline
			std::string command;
			std::cout << ">>> ";
			std::getline(std::cin,command);
			if (!std::cin)
			{
				std::cout << '\n';
				break;
			}
#endif
			// check for comments
			char first_char;
			std::istringstream com_in(command);
			com_in >> first_char; 
			if (first_char == '#' || !com_in)
				continue;

			// remove a trailing comment
			size_t comment_pos = command.find('#');
			if (comment_pos != command.npos)
				command.erase(comment_pos);

			// process command and output
			std::cout << process_command(command, ncl, n_errors, end) << '\n';
		}
	}
	catch (...) // there is always a sdt::logic_error when terminating the
				// program using EOF (typing C-d). This catch avoids nasty 
				// error messages
	{
		std::cout << '\n';
	}
	
	return n_errors;
}


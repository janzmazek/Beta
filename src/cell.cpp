#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

#include "cell.h"

#define RANKSZ_MAX 2147483648.0		// constant for random number generation (not seed)

// Noise generating variables
long ksz_num[256];
int ksz_i;
long seme;
double fac, rsq, v1, v2;

// NOISE ROUTINES
// randl -> RANDOM LONG FROM 0 TO (num)
// randd -> RANDOM DOBULE FROM 0 TO 1
// gasdev -> RANDOM DOUBLE GAUSS DISTR. FROM 0 TO 1
//
void Cell::seed(long seed)
{
	ksz_num[0] = (long)fmod(16807.0*(double)seed, 2147483647.0);
	for (ksz_i = 1; ksz_i<256; ksz_i++)
	{
		ksz_num[ksz_i] = (long)fmod(16807.0 * (double)ksz_num[ksz_i - 1], 2147483647.0);
	}
}

long Cell::randl(long num)
{
	ksz_i = ++ksz_i & 255;
	ksz_num[ksz_i] = ksz_num[(ksz_i - 103) & 255] ^ ksz_num[(ksz_i - 250) & 255];
	ksz_i = ++ksz_i & 255;
	ksz_num[ksz_i] = ksz_num[(ksz_i - 30) & 255] ^ ksz_num[(ksz_i - 127) & 255];
	return(ksz_num[ksz_i] % num);
}

double Cell::randd(void)
{
	ksz_i = ++ksz_i & 255;
	ksz_num[ksz_i] = ksz_num[(ksz_i - 103) & 255] ^ ksz_num[(ksz_i - 250) & 255];
	ksz_i = ++ksz_i & 255;
	ksz_num[ksz_i] = ksz_num[(ksz_i - 30) & 255] ^ ksz_num[(ksz_i - 127) & 255];
	return((double)ksz_num[ksz_i] / RANKSZ_MAX);
}

double Cell::gasdev(void)
{
	static int iset = 0;
	static double gset;

	if (iset == 0)
	{
		do
		{
			v1 = (float)(2.0*randd() - 1.0);
			v2 = (float)(2.0*randd() - 1.0);
			rsq = v1*v1 + v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = (float)(sqrt(-2.0*log(rsq) / rsq));
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}

double Cell::variability(double mean) {
    return mean;
}


Cell::Cell(double glucose, double dt) {
    m_glucose = glucose;
    m_dt = dt;
    double init_vals[6] = {-63, 0, 0.25, 0.1, 0.05, 0};
    m_results.insert(m_results.end(), begin(init_vals), end(init_vals));
}

Cell::~Cell() {

}

void Cell::set_initial(vector<double> initial) {
    m_results = initial;
}

void Cell::update(void) {
    double V = m_results[0];
    double n = m_results[1];
    double s = m_results[2];
    double ci = m_results[3];
    double cer = m_results[4];
    double binser = m_results[5];
// ------------------------------ Parameters ------------------------------- //
    double gama = 1.25; // Hill coef. glc vs. GKatp
    double VCa = 25.0; // Rev. pot. Ca2+ channel
    double VK = -75.0; // Rev. pot. K+ channel
    double Vm = -20.0; // half-max pot. for minf curve
    double Tm = 12.0; // Voltage const.
    double Vn = -16.0; // half-max pot. for ninf curve
    double Tn = 5.6; // Voltage const.
    double Vs = -52.0; // half-max pot. for sinf curve
    double Ts = 10.0; // Voltage const.
    double taun = 20.0; // time const
    double taus = 20000.0; // time const
    double fc = 0.001; // fraction of free Ca2+ in intracellular space

    double cm = variability(5300.0); // Plasma membrane capacitance
    double gca = variability(1000.0); 
    double gkatp = variability(220.0); // 150
    double gk = variability(2700.0);
    double gs = variability(200.0);
    double gkatp50 = variability(10.0); // 10
    double alfa = variability(0.000005); // conversion from electrical into chemical gradient
    double kca = variability(0.6); // 0.2 / removal rate of Ca2+ from intracellular space

    double gcran = 1.5;
    double Vcran = 0.0;

    double C = 0.0;

// ------------------------------- Functions ------------------------------- //
    // Fluxes
    double Okatp, ICa, IKatp, IK, Is, Icran;

    // Activating functions
    double minf, ninf, sinf, rinf;

    // empirical relation between glc concentration and Gkatp conductance
    Okatp = 1.0/(1.0 + pow((m_glucose/gkatp50), gama));

    // simple electrophysiological model, Stamper J Theor Biol 475 2019; 
    // based on Nittala et al.,PLoS ONE, 2 (2007), p. e983
    minf = 1.0/(1.0 + (exp((Vm-V)/Tm)));
    ninf = 1.0/(1.0 + (exp((Vn-V)/Tn)));
    sinf = 1.0/(1.0 + (exp((Vs-V)/Ts)));

    ICa = gca*minf*(V - VCa);
    IKatp = gkatp*Okatp*(V - VK);
    IK = gk*n*(V - VK);
    Is = gs*s*(V - VK);

    // Icran activation function
    rinf = 1.0/(1.0+(exp((cer-10.0)/3.0)));
    // calcium release activated nonselective cation current,
    // Mears et al., J. Membrane Biol. 155, 47, 1997
    Icran = gcran*rinf*cer*(V-Vcran);
    Icran = 0.0;

// ------------------------ Differential equations ------------------------- //
    V += m_dt*( -(IK + ICa + Is + IKatp + Icran + C)/cm);
    n += m_dt*((ninf-n)/taun);
    s += m_dt*((sinf-s)/taus);
    // original, simplistic description of Ca2+ dynamics, 
    // Stamper J Theor Biol 475 2019
    ci += m_dt*(-fc*(alfa*ICa + kca*ci));

    m_results[0] = V;
    m_results[1] = n;
    m_results[2] = s;
    m_results[3] = ci;
    m_results[4] = cer;
    m_results[5] = binser;
}

void Cell::print_results(void) {
    double V = m_results[0];
    cout << "V: " << V << "\n";
}

vector<double> Cell::get_result(void) {
    return m_results;
}

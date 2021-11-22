#include <iostream>
#include <vector>
#include <math.h>
#include <random>
using namespace std;

#include "cell.h"

double Cell::variability(double mean, double std, mt19937 generator) {
    normal_distribution<double> distribution(mean, std);
    double number = distribution(generator);

    return number;
}

Cell::Cell(double glucose, double dt, double initial[6], mt19937 generator) {
    m_glucose = glucose;
    m_dt = dt;
    for (size_t i = 0; i < 6; ++i)
        m_results.push_back(initial[i]);
    m_cm = variability(5300, 1, generator);
    m_gca = variability(1000, 1, generator);
    m_gkatp = variability(220, 1, generator); // 150
    m_gk = variability(2700, 1, generator);
    m_gs = variability(200, 1, generator);
    m_gkatp50 = variability(10, 1, generator);
    m_alfa = variability(5e-6, 1e-7, generator);
    m_kca = variability(0.6, 1e-3, generator);
    cout << m_cm << " ";
}
Cell::~Cell() {}

void Cell::update(double Icoupling) {
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

    double cm = m_cm; // Plasma membrane capacitance
    double gca = m_gca; 
    double gkatp = m_gkatp; // 150
    double gk = m_gk;
    double gs = m_gs;
    double gkatp50 = m_gkatp50; // 10
    double alfa = m_alfa; // conversion from electrical into chemical gradient
    double kca = m_kca; // 0.2 / removal rate of Ca2+ from intracellular space

    double gcran = 1.5;
    double Vcran = 0.0;

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
    V += m_dt*( -(IK + ICa + Is + IKatp + Icran + Icoupling)/cm);
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
    for (int i=0; i<6; i++) {
        cout << m_results[i];
        
    }
    cout << "\n";
    // double V = m_results[0];
    // cout << "V: " << V << "\n";
}

vector<double> Cell::get_result(void) {
    return m_results;
}

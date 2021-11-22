#include <iostream>
#include <vector>
#include <math.h>
#include <random>
using namespace std;

#include "cell.h"

double Cell::variability(double mean, double std) {
    // random_device rd;
    // mt19937 e2(rd());
    // normal_distribution<> dist(mean, std);
    // return dist(e2);
    return mean;
}


Cell::Cell(double glucose, double dt, double initial[6]) {
    m_glucose = glucose;
    m_dt = dt;
    for (size_t i = 0; i < 6; ++i)
        m_results.push_back(initial[i]);
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

    double cm = variability(5300.0, 1); // Plasma membrane capacitance
    double gca = variability(1000.0, 1); 
    double gkatp = variability(220.0, 1); // 150
    double gk = variability(2700.0, 1);
    double gs = variability(200.0, 1);
    double gkatp50 = variability(10.0, 1); // 10
    double alfa = variability(0.000005, 1e-7); // conversion from electrical into chemical gradient
    double kca = variability(0.6, 1e-3); // 0.2 / removal rate of Ca2+ from intracellular space

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

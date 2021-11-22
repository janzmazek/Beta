#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include "cell.cpp"

int main(void) {
    int L = 5;
	// Integration parameters
    double dt = 0.1;
    double t_fin = 300000;
	double write_interval = 10;
    double init_vals[6] = {-63, 0, 0.25, 0.1, 0.05, 0};

    ofstream MyFile("filename.txt");

	vector<double> result;
    vector<Cell> cells;

    for (int i=0; i<L; i++) {
        Cell cell(8, dt, init_vals);
        cells.push_back(cell);
        }
    
    // Loop over time
    for (double t=0; t<t_fin; t+=dt) {
        // Loop over cells
        for (int i=0; i<L; i++) { 
            cells[i].update(0);
        }
        if (fmod(t, write_interval)<0.1) {
            MyFile << t << " ";
            for (int i = 0; i < L; i++)
            {
                result = cells[i].get_result();
                MyFile << result[0] << " ";
            }
            MyFile << "\n";
        }
    }

    MyFile.close();

	return 0;
}
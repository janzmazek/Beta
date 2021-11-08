#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "cell.cpp"

int main(void) {
	// Integration parameters
    double dt = 0.1;
    double t_fin = 300000;
	double write_interval = 10;

    ofstream MyFile("filename.txt");

	vector<double> result;

	Cell cell(8, dt);
    for (double t=0; t<t_fin; t+=dt) {
		cell.update();
        if (fmod(t, write_interval)<0.1) {
			result = cell.get_result();
			MyFile << t << " " << result[0] << "\n";
        }
    }

    MyFile.close();

	return 0;
}
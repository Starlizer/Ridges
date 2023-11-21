
#include "bj_flow.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>


double delta(int i, double eta_max, int Nl) 
{

	double PI = 3.14;
	double sigma2 = 2*.005;
	double eta = -eta_max + 2 * eta_max * i / (Nl-1);
	double norm = 0;
	for (int ti = 0; ti < Nl; ti++)
	{
		double teta = -eta_max + 2 * eta_max * ti / (Nl - 1);
		norm += exp(-(teta * teta) / (sigma2));
	}
	return 1. / norm * exp(-(eta * eta) / (sigma2));
}
void make_greens_functions(std::string Q_source, int source_num, double tstart, double tend, int Ntstart, double dt, int Nl, double eta_max)
{

	std::ofstream file11;
	std::ofstream file12;
	std::ofstream file13;
	std::ofstream file2;
	std::ofstream file3;
	file11.open("Green-" + Q_source + "-B.dat");
	file12.open("Green-" + Q_source + "-Q.dat");
	file13.open("Green-" + Q_source + "-S.dat");
	file2.open("Green-" + Q_source + "-E.dat");
	file3.open("Green-" + Q_source + "-Uz.dat");

	bj_flow system;
	for (int ti = 0; ti < Ntstart; ti++)
	{
		double t1 = tstart + (tend - tstart) * ti / (Ntstart - 1);
		system.t = t1;
		std::vector<std::vector<double>> IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++)
		{
			IC[source_num][i] = delta(i, eta_max, Nl) / t1;
		}
		system.initial_conditions(IC);
		int Nt = (tend - t1) / dt;
		for (int i = 0; i < (tend - t1) / dt; i++)
		{
			system.step();
		}
		std::vector<double> t_B = system.get_B();
		std::vector<double> t_Q = system.get_Q();
		std::vector<double> t_S = system.get_S();
		std::vector<double> t_e = system.get_e();
		std::vector<double> t_u = system.get_u();
		for (int j = 0; j < Nl; j++)
		{
			file11 << t_B[j] << "\t";
		}
		file11 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file12 << t_Q[j] << "\t";
		}
		file12 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file13 << t_S[j] << "\t";
		}
		file13 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file2 << t_e[j] << "\t";
		}
		file2 << std::endl;
		for (int j = 0; j < Nl; j++)
		{
			file3 << t_u[j] << "\t";
		}
		file3 << std::endl;

	}


	file11.close();
	file12.close();
	file13.close();
	file2.close();
	file3.close();




}
void make_greens_functions_transverse(std::string Q_source, int source_num, double tstart, double tend, int Ntstart, double dt, int Nl, double eta_max)
{


	std::ofstream file3;
	file3.open("Green-" + Q_source + "-Ux.dat");


	bj_flow system;
	for (int ti = 0; ti < Ntstart; ti++)
	{
		double t1 = tstart + (tend - tstart) * ti / (Ntstart - 1);
		system.t = t1;
		std::vector<std::vector<double>> IC(5, std::vector<double>(Nl, 0));
		for (int i = 0; i < Nl; i++)
		{
			IC[source_num][i] = delta(i, eta_max, Nl) / t1;
		}
		system.initial_conditions(IC);
		int Nt = (tend - t1) / dt;
		for (int i = 0; i < (tend - t1) / dt; i++)
		{
			system.step_transverse();
		}
		std::vector<double> t_u = system.get_u();
		for (int j = 0; j < Nl; j++)
		{
			file3 << t_u[j] << "\t";
		}
		file3 << std::endl;

	}



	file3.close();




}


int main(int argc, char* argv[])
{
	std::string name = "n0";
	int Nl = 501;
	double eta_max=5.;
	int Ntstart = 500;
	double tstart = 1;
	double tend = 11;
	double dt = 10. / 500 * 0.001;
	
	//make_greens_functions(argv[1], atoi(argv[2]),tstart, tend, Ntstart, dt, Nl, eta_max);

	make_greens_functions_transverse("Ux", 4, tstart, tend, Ntstart, dt, Nl, eta_max);

	return 1;
}

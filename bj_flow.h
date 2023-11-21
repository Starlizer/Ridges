#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
class bj_flow
{
	const double eta_max = 5;
	const int Nl = 501;
	const int pow_pol = 11;
	std::string name_rho = "0n0";
	const double dt = 0.001 * 2 * eta_max /(Nl-1);
	const double deta = 2 * eta_max / (Nl - 1);
	const double t_start = 1;
	double rho_startB = 0 * 0.16;
	double rho_startQ = .4*rho_startB;


	std::vector<double> pol_Mst;
	std::vector<double> pol_dMst;
	std::vector<double> pol_cs2;
	std::vector<double> pol_dBP;
	std::vector<double> pol_dQP;
	std::vector<double> pol_dSP;
	std::vector<double> pol_etas;
	std::vector<double> pol_D;

	std::vector<double> u;
	std::vector<double> e;
	std::vector<double> B;
	std::vector<double> Q;
	std::vector<double> S;
	std::vector<double> tu;
	std::vector<double> te;
	std::vector<double> tB;
	std::vector<double> tQ;
	std::vector<double> tS;
	double cs2(double);
	double dBP(double);
	double dQP(double);
	double dSP(double);
	double etas(double);
	double Mst(double);
	double dMst(double);
	double D(double);
	double w(double);
	double rhoB(double);
	double dt_rhoB(double);
	double rhoQ(double);
	double dt_rhoQ(double);

	void read_cs2();
	void read_dBP();
	void read_dQP();
	void read_dSP();
	void read_etas();
	void read_Mst();
	void read_dMst();
	void read_D();

public:


	void initial_conditions(std::vector<std::vector<double>>& IC);
	void step();	
	void step_transverse();
	double t;
	std::vector<double> get_B();
	std::vector<double> get_Q();
	std::vector<double> get_S();
	std::vector<double> get_u();
	std::vector<double> get_e();
	bj_flow()
	{
		tB = std::vector<double>(Nl, 0);
		tQ = std::vector<double>(Nl, 0);
		tS = std::vector<double>(Nl, 0);
		te = std::vector<double>(Nl, 0);
		tu = std::vector<double>(Nl, 0);
		B = std::vector<double>(Nl, 0);
		Q = std::vector<double>(Nl, 0);
		S = std::vector<double>(Nl, 0);
		e = std::vector<double>(Nl, 0);
		u = std::vector<double>(Nl, 0);
//       read_cs2();
//		read_dBP();
//		read_dQP();
//		read_dSP();
		read_etas();
		read_Mst();
		read_dMst();
//		read_D();
	}

};


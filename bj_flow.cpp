#include "bj_flow.h"
void bj_flow::initial_conditions(std::vector<std::vector<double>>& IC)
{
	for (int i = 0; i < Nl; i++)
	{
		B[i] = IC[0][i];
		Q[i] = IC[1][i];
		S[i] = IC[2][i];
		e[i] = IC[3][i];
		u[i] = IC[4][i];
	}
	B[0] = 0;
	Q[0] = 0;
	S[0] = 0;
	e[0] = 0;
	u[0] = 0;
	u[Nl - 1] = 0;
	e[Nl - 1] = 0; 
	B[Nl - 1] = 0;
	Q[Nl - 1] = 0;
	S[Nl - 1] = 0;
}

double bj_flow::rhoB(double t)
{
	return rho_startB * t_start / t;
}
double bj_flow::dt_rhoB(double t)
{
	return -rhoB(t) / t;
}
double bj_flow::rhoQ(double t)
{
	return rho_startQ * t_start / t;
}
double bj_flow::dt_rhoQ(double t)
{
	return -rhoQ(t) / t;
}



void bj_flow::step()
{

	double t2 = t * t;
	double ddt = 2 * deta;
	double deta2 = deta * deta;
	double dB_P = dBP(t);
	double dQ_P = dQP(t);
	double dS_P = dSP(t);
	double dP = cs2(t)*e[0] + dB_P * B[0] +dB_P * Q[0] + dS_P * S[0];
	double tMst = Mst(t);
	double tD = D(t);
	double trhoB = rhoB(t);
	double tdt_rhoB = dt_rhoB(t);
	double trhoQ = rhoQ(t);
	double tdt_rhoQ = dt_rhoQ(t);
	double tdMst = dMst(t);
	double tetas = etas(t);
	double tcs2 = cs2(t);
	double deta_u = (u[1] - u[Nl-1]) / (ddt);
	double deta_B = (B[1] - B[Nl - 1]) / (ddt);
	double deta_Q = (Q[1] - Q[Nl - 1]) / (ddt);
	double deta_S = (S[1] - S[Nl - 1]) / (ddt);
	double deta_e = (e[1] - e[Nl - 1]) / (ddt);
	double d2eta_B = (B[1] + B[Nl - 1] - 2 * B[0]) / deta2;
	double d2eta_Q = (Q[1] + Q[Nl - 1] - 2 * Q[0]) / deta2;
	double d2eta_S = (S[1] + S[Nl - 1] - 2 * S[0]) / deta2;
	double d2eta_u = (u[1] + u[Nl - 1] - 2 * u[0]) / deta2;
	tB[0] = -B[0] / t + tD * d2eta_B / t2 + (trhoB + tD * tdt_rhoB) * deta_u / t /tMst;
	tQ[0] = -Q[0] / t + tD * d2eta_Q / t2 + (trhoQ + tD * tdt_rhoQ) * deta_u / t/tMst;
	tS[0] = -S[0] / t + tD * d2eta_S / t2;
	te[0] = -(e[0] + dP) / t  - (1 - 4 * tetas / (3 * t)/tMst) * deta_u / t;
	tu[0] = -(2.  / t ) * u[0] - (tcs2 * deta_e + dB_P * deta_B + dQ_P * deta_Q + dS_P * deta_S) / t + 4./3. * tetas / t2 * d2eta_u / tMst;
	for (int i = 1; i < Nl-1; i++)
	{
		dP = tcs2 * e[i] + dB_P * B[i] + dB_P * Q[i] + dS_P * S[i];
		deta_u = (u[i + 1] - u[i - 1]) / (ddt);
		deta_B = (B[i + 1] - B[i - 1]) / (ddt);
		deta_Q = (Q[i + 1] - Q[i - 1]) / (ddt);
		deta_S = (S[i + 1] - S[i - 1]) / (ddt);
		deta_e = (e[i + 1] - e[i - 1]) / (ddt);
		d2eta_B = (B[i + 1] + B[i - 1] - 2. * B[i]) / deta2;
		d2eta_Q = (Q[i + 1] + Q[i - 1] - 2. * Q[i]) / deta2;
		d2eta_S = (S[i + 1] + S[i - 1] - 2. * S[i]) / deta2;
		d2eta_u = (u[i + 1] + u[i - 1] - 2. * u[i]) / deta2;
		tB[i] = -B[i] / t + tD * d2eta_B / t2 + (trhoB + tD * tdt_rhoB) * deta_u / t/tMst;
		tQ[i] = -Q[i] / t + tD * d2eta_Q / t2 + (trhoQ + tD * tdt_rhoQ) * deta_u / t/tMst;
		tS[i] = -S[i] / t + tD * d2eta_S / t2;
		te[i] = -(e[i]+dP)/t - (1-4.*tetas/(3.*t*tMst)) * deta_u/t;
		tu[i] = -(2. / t) *u[i] - (tcs2 * deta_e + dB_P * deta_B + dQ_P * deta_Q+ dS_P * deta_S) / t +4./3.*tetas/t2*d2eta_u/tMst;
	}
	dP = tcs2 * e[Nl - 1] + dB_P * B[Nl - 1] + dB_P * Q[Nl - 1] + dS_P * S[Nl - 1];
	deta_u = -(u[Nl - 2] - u[0]) / (ddt);
	deta_B = -(B[Nl - 2] - B[0]) / (ddt);
	deta_Q = -(Q[Nl - 2] - Q[0]) / (ddt);
	deta_S = -(S[Nl - 2] - S[0]) / (ddt);
	deta_e = -(e[Nl - 2] - e[0]) / (ddt);
	d2eta_B = (B[Nl - 2] + B[0] - 2 * B[Nl - 1]) / deta2;
	d2eta_Q = (Q[Nl - 2] + Q[0] - 2 * Q[Nl - 1]) / deta2;
	d2eta_S = (S[Nl - 2] + S[0] - 2 * S[Nl - 1]) / deta2;
	d2eta_u = (u[Nl - 2] + u[0] - 2 * u[Nl - 1]) / deta2;

	tB[Nl - 1] = -B[Nl - 1] / t + tD * d2eta_B / t2 + (trhoB + tD * tdt_rhoB) * deta_u / t/tMst;
	tQ[Nl - 1] = -Q[Nl - 1] / t + tD * d2eta_Q / t2 + (trhoQ + tD * tdt_rhoQ) * deta_u / t/tMst;
	tS[Nl - 1] = -S[Nl - 1] / t + tD * d2eta_S / t2;
	te[Nl - 1] = -(e[Nl - 1] + dP) / t - (1 - 4 * tetas / (3 * t * tMst)) * deta_u / t;
	tu[Nl - 1] = -(2./t) * u[Nl - 1]  - (tcs2 * deta_e + dB_P * deta_B + dQ_P * deta_Q + dS_P * deta_S) / t + 4./3. * tetas / t2 * d2eta_u / tMst;

	for (int i = 0; i < Nl; i++)
	{
		B[i] += dt * tB[i];
		Q[i] += dt * tQ[i];
		S[i] += dt * tS[i];
		e[i] += dt * te[i];
		u[i] += dt * tu[i];

	}


	t = t + dt;
}

void bj_flow::step_transverse()
{

	double t2 = t * t;
	double ddt = 2 * deta;
	double deta2 = deta * deta;
	
	double tMst = Mst(t);
	double tdMst = dMst(t);
	double tetas = etas(t);
	double deta_u = (u[1] - u[Nl - 1]) / (ddt);
	double d2eta_u = (u[1] + u[Nl - 1] - 2 * u[0]) / deta2;
	tu[0] = -(1. / t) * u[0]  + tetas / t2 * d2eta_u / tMst;
	for (int i = 1; i < Nl - 1; i++)
	{
		deta_u = (u[i + 1] - u[i - 1]) / (ddt);
        d2eta_u = (u[i+1] + u[i - 1] - 2 * u[i]) / deta2;
		tu[i] = -(1. / t) * u[i] + tetas / t2 * d2eta_u / tMst;
	}
	deta_u = -(u[Nl - 2] - u[0]) / (ddt);
	d2eta_u = (u[Nl - 2] + u[0] - 2 * u[Nl - 1]) / deta2;
	tu[Nl - 1] = -(1. / t) * u[Nl - 1]+ tetas / t2 * d2eta_u / tMst;

	for (int i = 0; i < Nl; i++)
	{
		u[i] += dt * tu[i];
	}


	t = t + dt;
}
std::vector<double> bj_flow::get_u()
{
	return u;
}

std::vector<double> bj_flow::get_e()
{
	return e;
}

std::vector<double> bj_flow::get_B()
{
	return B;
}
std::vector<double> bj_flow::get_Q()
{
	return Q;
}
std::vector<double> bj_flow::get_S()
{
	return S;
}

void bj_flow::read_Mst()
{
	std::ifstream Msts(name_rho + "-Mst.dat");
	for (int i = 0; i < pow_pol; i++)
	{
		double tmst;
		Msts >> tmst;
		pol_Mst.push_back(tmst);
	}
}
void bj_flow::read_dMst()
{
	std::ifstream Msts(name_rho + "-dMst.dat");
	for (int i = 0; i < pow_pol; i++)
	{
		double tmst;
		Msts >> tmst;
		pol_dMst.push_back(tmst);
	}
}
void bj_flow::read_etas()
{
	std::ifstream Msts(name_rho + "-etas.dat");
	for (int i = 0; i < pow_pol; i++)
	{
		double tmst;
		Msts >> tmst;
		pol_etas.push_back(tmst);
	}
}
void bj_flow::read_cs2()
{
	std::ifstream Msts(name_rho + "-dpde.dat");
	for (int i = 0; i < pow_pol; i++)
	{
		double tmst;
		Msts >> tmst;
		pol_cs2.push_back(tmst);
	}
}
void bj_flow::read_dBP()
{
	std::ifstream Msts(name_rho + "-dpdrhoB.dat");
	for (int i = 0; i < pow_pol; i++)
	{
		double tmst;
		Msts >> tmst;
		pol_dBP.push_back(tmst);
	}
}
void bj_flow::read_dQP()
{
	std::ifstream Msts(name_rho + "-dpdrhoQ.dat");
	for (int i = 0; i < pow_pol; i++)
	{
		double tmst;
		Msts >> tmst;
		pol_dQP.push_back(tmst);
	}
}
void bj_flow::read_dSP()
{
	std::ifstream Msts(name_rho + "-dpdrhoS.dat");
	for (int i = 0; i < pow_pol; i++)
	{
		double tmst;
		Msts >> tmst;
		pol_dSP.push_back(tmst);
	}
}
void bj_flow::read_D()
{
	std::ifstream Msts(name_rho + "-D.dat");
	for (int i = 0; i < pow_pol; i++)
	{
		double tmst;
		Msts >> tmst;
		pol_D.push_back(tmst);
	}
}

double bj_flow::cs2(double x)
{
	double tx = 1.;
	double ty = 0.;
	for (int i = 0; i < pow_pol; i++)
	{
		ty += pol_cs2[pow_pol - i - 1] * tx;
		tx = tx * x;
	}
	return ty;
}
double bj_flow::dBP(double x)
{
	double tx = 1.;
	double ty = 0.;
	for (int i = 0; i < pow_pol; i++)
	{
		ty += pol_dBP[pow_pol - i - 1] * tx;
		tx = tx * x;
	}
	return ty;
}

double bj_flow::dQP(double x)
{
	double tx = 1.;
	double ty = 0.;
	for (int i = 0; i < pow_pol; i++)
	{
		ty += pol_dQP[pow_pol - i - 1] * tx;
		tx = tx * x;
	}
	return ty;
}

double bj_flow::dSP(double x)
{
	double tx = 1.;
	double ty = 0.;
	for (int i = 0; i < pow_pol; i++)
	{
		ty += pol_dSP[pow_pol - i - 1] * tx;
		tx = tx * x;
	}
	return ty;
}
double bj_flow::etas(double x)
{
	double tx = 1;
	double ty = 0;
	for (int i = 0; i < pow_pol; i++)
	{
		ty += pol_etas[pow_pol - i - 1] * tx;
		tx = tx * x;
	}
	return ty;
}


double bj_flow::Mst(double x)
{
	double tx = 1;
	double ty = 0;
	for (int i = 0; i < pow_pol; i++)
	{
		ty += pol_Mst[pow_pol - i - 1] * tx;
		tx = tx * x;
	}
	return ty;
}

double bj_flow::dMst(double x)
{
	double tx = 1;
	double ty = 0;
	for (int i = 1; i < pow_pol; i++)
	{
		ty += i*pol_Mst[pow_pol - i - 1] * tx;
		tx = tx * x;
	}
	return ty;
}

double bj_flow::w(double x)
{
	return Mst(x) + 4. / 3. * etas(x);
}

double bj_flow::D(double x)
{
	double tx = 1;
	double ty = 0;
	for (int i = 0; i < pow_pol; i++)
	{
		ty += pol_D[pow_pol - i-1] * tx;
		tx = tx * x;
	}
	return ty;
}

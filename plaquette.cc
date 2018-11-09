#include "plaquette.h"

void MeasPlq(multi1d<LatticeColorMatrix>& u, Double& plaquette){

	int n_planes = Nd*(Nd-1)/2;
	LatticeColorMatrix plaq = zero;
	for(int mu=0; mu<Nd; ++mu){
		for(int nu=mu+1; nu<Nd; ++nu){
			LatticeColorMatrix tmp, tmp2, tmp3;
			tmp = shift(u[nu], FORWARD, mu);
			tmp2 = u[mu] * tmp;
			tmp = shift(u[mu], FORWARD, nu);
			tmp3 = u[nu] * tmp;

			plaq += tmp2 * adj(tmp3);
		}
	}

	Double normalize = Real(3) * Real(n_planes) * Layout::vol();
	plaquette = (Double(1)/normalize) * sum(real(trace(plaq)));
}

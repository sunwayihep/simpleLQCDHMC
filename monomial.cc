#include "monomial.h"

Double WilsonGaugeMonomial::S(const AbsFieldState<GaugeP, GaugeQ>& s) const {
	//Double S1 = zero;
	//Double S2 = zero;
	Double S = zero;
	const GaugeQ& u = s.getQ();

	//for(int mu=1; mu<Nd; ++mu){
	//	for(int nu=0; nu<mu; ++nu){
	//		S1 += Nc;
	//		S2 += sum(real(trace(u[mu]
	//						*shift(u[nu],FORWARD,mu)
	//						*adj(shift(u[mu],FORWARD,nu))
	//						*adj(u[nu]))));
	//	}
	//}
	////S *= Double(-beta)/Double(Nc);
	//S1 -= S2;
	//S1 *= Double(beta)/Double(Nc);
	//return S1;

	for(int mu=1; mu<Nd; ++mu){
		for(int nu=0; nu<mu; ++nu){
			S += sum(real(trace(u[mu]
							*shift(u[nu],FORWARD,mu)
							*adj(shift(u[mu],FORWARD,nu))
							*adj(u[nu]))));
		}
	}
	S *= Double(-beta)/Double(Nc);
	return S;
}

void WilsonGaugeMonomial:: dsdq(GaugeP& F, const GaugeQ& u) const {
	F.resize(Nd);
	LatticeColorMatrix tmp_0;
	F = zero;
	for(int mu=0; mu<Nd; mu++){
		for(int nu=mu+1; nu<Nd; nu++){
			tmp_0 = adj(shift(u[mu], FORWARD, nu))*adj(u[nu]);
			F[mu] += shift(u[nu], FORWARD, mu)*tmp_0;
			F[nu] += shift(tmp_0*u[mu], BACKWARD, mu);
			tmp_0 = adj(shift(u[nu], FORWARD, mu))*adj(u[mu]);
			F[mu] += shift(tmp_0*u[nu], BACKWARD, nu);
			F[nu] += shift(u[mu], FORWARD, nu)*tmp_0;
		}
		tmp_0 = Real(-beta)/(Real(2*Nc))*F[mu];
		F[mu] = u[mu] * tmp_0;
	}
}

void
TwoFlavorWilsonFermMonomial::refreshInternalFields(
		const AbsFieldState<GaugeP, GaugeQ>& s){
	const GaugeQ& u = s.getQ();
	UnprecWilsonLinOp M(u, Mass);

	LatticeDiracFermion eta;
	gaussian(eta);
	eta *= sqrt(0.5);
	// apply M^dagger: phi = M^+ * eta
	M(phi, eta, -1);
}

void
TwoFlavorWilsonFermMonomial:: getX(LatticeFermion& X, const GaugeQ& u) const {
	UnprecWilsonLinOp M(u, Mass);
	Real RsdCGOut;
	int n_count;

	InvCG_a(M, phi, X, RsdCG, MaxCG, n_count, RsdCGOut);
}

Double
TwoFlavorWilsonFermMonomial:: S(const AbsFieldState<GaugeP, GaugeQ>& s) const {
	const GaugeQ& u = s.getQ();
	LatticeFermion X = zero;
	getX(X, u);
	Double result = real(innerProduct(phi, X));
	return result;
}

void TwoFlavorWilsonFermMonomial:: dsdq(GaugeP& F, const GaugeQ& u) const {
	UnprecWilsonLinOp M(u, Mass);
	LatticeDiracFermion X, Y;

	getX(X, u);
	M(Y, X, 1);

	GaugeP F_tmp;
	M.deriv(F_tmp, X, Y, -1);
	M.deriv(F, Y, X, +1);
	for(int mu=0; mu<Nd; ++mu){
		F_tmp[mu] += F[mu];
		F_tmp[mu] *= Real(-1);
	}
	
	for(int mu=0; mu<F.size(); ++mu){
		F[mu] = u[mu]*F_tmp[mu];
	}
}

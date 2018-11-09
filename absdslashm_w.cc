#include "absdslashm_w.h"


void UnprecWilsonLinOp:: operator()(LatticeDiracFermion& result,
		const LatticeDiracFermion& source,
		int isign) const {
	Real mass_term = Real(Nd) + Mass;
	Real half = Real(0.5);

	result = mass_term * source;
	LatticeDiracFermion tmp;
	dslash(tmp, u, source, isign, 0);
	dslash(tmp, u, source, isign, 1);

	result -= half * tmp;
}

void UnprecWilsonLinOp:: deriv(multi1d<LatticeColorMatrix>& F,
		const LatticeDiracFermion& X,
		const LatticeDiracFermion& Y,
		int isign) const {
	F.resize(Nd);
	for(int mu=0; mu<Nd; mu++){ F[mu] = zero; }
	multi1d<LatticeColorMatrix> F_tmp(Nd);
	dslash_deriv(F, X, Y, isign, 0);
	//dslash_deriv(F_tmp, Y, X, isign, 1);
	dslash_deriv(F_tmp, X, Y, isign, 1);
	F += F_tmp;

	for(int mu=0; mu<Nd; mu++){
		F[mu] *= Real(-0.5);
	}
}

void InvMR(const LinearOperator<LatticeDiracFermion>& M,
		const LatticeDiracFermion& source,
		LatticeDiracFermion& target,
		const Real& MRovpar,
		const Real& RsdMR,
		int MaxMR,
		int isign,
		int& n_count,
		Double& resid) {

	int k = 0;
	Double d;
	DComplex c, a;
	LatticeDiracFermion Mr, r;
	QDPIO::cout<<"I'm in abs MR inverter"<<std::endl;
	M(Mr, target, isign);
	r = source - Mr;
	resid = sqrt(norm2(r));
	while((k<MaxMR) && (toBool(resid > RsdMR))) {
		++k;
		M(Mr, r, isign);
		c = innerProduct(Mr, r);
		d = norm2(Mr);
		a = c / d;
		a *= MRovpar;

		target += a * r;
		r -= a * Mr;
		resid = sqrt(norm2(r));
#ifdef QUIET
#if 0
		QDPIO::cout<<"abs # "<<k<<" inner residual = "<<resid<<std::endl;
#endif
#endif
	}
	n_count = k;
	QDPIO::cout<<"# of iteration = "<<n_count<<std::endl;
	QDPIO::cout<<"residual = "<<resid<<std::endl;
}


template<typename T>
void InvMR_a(const LinearOperator<T>& M,
		const T& source,
		T& target,
		const Real& MRovpar,
		const Real& RsdMR,
		int MaxMR,
		int isign,
		int& n_count,
		Double& resid) { 

	int k = 0;
	Double d;
	DComplex c, a;
	T Mr, r;
	QDPIO::cout<<"I'm in templated abs MR inverter"<<std::endl;
	M(Mr, target, isign);
	r = source - Mr;
	resid = sqrt(norm2(r));
	while((k<MaxMR) && (toBool(resid > RsdMR))) {
		++k;
		M(Mr, r, isign);
		c = innerProduct(Mr, r);
		d = norm2(Mr);
		a = c / d;
		a *= MRovpar;

		target += a * r;
		r -= a * Mr;
		resid = sqrt(norm2(r));
#ifdef QUIET
#if 0
		QDPIO::cout<<"abs # "<<k<<" inner residual = "<<resid<<std::endl;
#endif
#endif
	}
	n_count = k;
	QDPIO::cout<<"# of iteration = "<<n_count<<std::endl;
	QDPIO::cout<<"residual = "<<resid<<std::endl;
}

template void InvMR_a<LatticeDiracFermion>(const LinearOperator<LatticeDiracFermion>& M,
		const LatticeDiracFermion& source,
		LatticeDiracFermion& target,
		const Real& MRovpar,
		const Real& RsdMR,
		int MaxMR,
		int isign,
		int& n_count,
		Double& resid);

template<typename T>
void InvCG_a(const LinearOperator<T>& M,
		const T& source,
		T& target,
		const Real& RsdCG,
		int MaxCG,
		int& n_count,
		Double& resid){
	int k = 0;
	DComplex alpha, beta, a, b;
	T tmp_0, tmp_1;
	M(tmp_0, target, 1);
	M(tmp_1, tmp_0, -1);
	T r = source - tmp_1;
	T p = r;
	T Mp;
	resid = sqrt(norm2(r));
	QDPIO::cout<<"I'm in CG solver---"<<std::endl;

	while((k<MaxCG) && (toBool(resid > RsdCG))){
		++k;
		a = innerProduct(r, r);
		M(Mp, p, 1);
		b = innerProduct(Mp, Mp);
		alpha = a / b;
		target += alpha * p;
		M(tmp_0, p, 1);
		M(tmp_1, tmp_0, -1);
		r -= alpha*tmp_1;

		b = innerProduct(r, r);
		beta = b / a;
		p = r + beta*p;
		resid = sqrt(norm2(r));
#ifdef QUIET
#if 0
		QDPIO::cout<<"abs # "<<k<<" inner residual = "<<resid<<std::endl;
#endif
#endif

	}
	n_count = k;
	QDPIO::cout<<"# of iterations = "<<n_count<<", residual = "<<resid<<std::endl;
}

template void InvCG_a<LatticeDiracFermion>(const LinearOperator<LatticeDiracFermion>& M,
		const LatticeDiracFermion& source,
		LatticeDiracFermion& target,
		const Real& RsdCG,
		int MaxCG,
		int& n_count,
		Double& resid);

void dslash_deriv(multi1d<LatticeColorMatrix>& F,
		const LatticeDiracFermion& X,
		const LatticeDiracFermion& Y,
		int isign, int cb){
	F.resize(Nd);
	for(int mu=0; mu<Nd; ++mu){
		LatticeDiracFermion temp_ferm1;
		LatticeHalfFermion tmp_h;

		switch(isign){
			case 1:
				{
					switch(mu){
						case 0:
							tmp_h[rb[1-cb]] = spinProjectDir0Minus(Y);
							temp_ferm1[rb[1-cb]] = spinReconstructDir0Minus(tmp_h);
							break;
						case 1:
							tmp_h[rb[1-cb]] = spinProjectDir1Minus(Y);
							temp_ferm1[rb[1-cb]] = spinReconstructDir1Minus(tmp_h);
							break;
						case 2:
							tmp_h[rb[1-cb]] = spinProjectDir2Minus(Y);
							temp_ferm1[rb[1-cb]] = spinReconstructDir2Minus(tmp_h);
							break;
						case 3:
							tmp_h[rb[1-cb]] = spinProjectDir3Minus(Y);
							temp_ferm1[rb[1-cb]] = spinReconstructDir3Minus(tmp_h);
							break;
						default:
							break;
					};
				}break;

			case -1:
				{
					switch(mu) {
						case 0:
							tmp_h[rb[1-cb]] = spinProjectDir0Plus(Y);
							temp_ferm1[rb[1-cb]] = spinReconstructDir0Plus(tmp_h);
							break;
						case 1:
							tmp_h[rb[1-cb]] = spinProjectDir1Plus(Y);
							temp_ferm1[rb[1-cb]] = spinReconstructDir1Plus(tmp_h);
							break;
						case 2:
							tmp_h[rb[1-cb]] = spinProjectDir2Plus(Y);
							temp_ferm1[rb[1-cb]] = spinReconstructDir2Plus(tmp_h);
							break;
						case 3:
							tmp_h[rb[1-cb]] = spinProjectDir3Plus(Y);
							temp_ferm1[rb[1-cb]] = spinReconstructDir3Plus(tmp_h);
							break;
						default:
							break;
					};
				}break;
			default:
				QDP_error_exit("unknown case");
		};

		LatticeDiracFermion temp_ferm2 = shift(temp_ferm1, FORWARD, mu);

		F[mu][rb[cb]] = traceSpin(outerProduct(temp_ferm2, X));
		F[mu][rb[1-cb]] = zero;
	}
}

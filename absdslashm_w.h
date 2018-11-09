#ifndef _ABSDSLASHM_W_H_
#define _ABSDSLASHM_W_H_
#include "qdp.h"
#include "dslashm_w.h"

using namespace QDP;

template<typename T>
class LinearOperator {
	public:
		virtual ~LinearOperator() {}
		virtual void operator() (T& result, const T& source, int isign) const = 0;
		virtual Subset& subset() const = 0;
};

template<typename P, typename T>
class DiffLinearOperator: public LinearOperator<T> {
	public:
		virtual ~DiffLinearOperator(){}
		virtual void operator()(T& result, const T& source, int isign) const = 0;
		virtual const Subset& subset() const = 0;

		virtual void deriv(P& F, const T& X, const T& Y, int isign) const = 0;
};

template<typename T>
class UnprecLinearOperator: public LinearOperator<T> {
	public:
		virtual ~UnprecLinearOperator(){}
		virtual void operator()(T& result, const T& source, int isign) const = 0;
		Subset& subset() const {return all;}
};

class UnprecWilsonLinOp: public UnprecLinearOperator<LatticeDiracFermion> {
	public:
		~UnprecWilsonLinOp(){}
		// constructor
		UnprecWilsonLinOp(const multi1d<LatticeColorMatrix>& u_,
				const Real& Mass_): u(u_), Mass(Mass_) {}

		void operator()(LatticeDiracFermion& result, const LatticeDiracFermion& source, int isign) const;
		void deriv(multi1d<LatticeColorMatrix>& F,
				const LatticeDiracFermion& X,
				const LatticeDiracFermion& Y,
				int isign) const;
	private:
		multi1d<LatticeColorMatrix> u;
		Real Mass;
};

//MR Inverter
void InvMR(const LinearOperator<LatticeDiracFermion>& M,
		const LatticeDiracFermion& source,
		LatticeDiracFermion& target,
		const Real& MRovpar,
		const Real& RsdMR,
		int MaxMR,
		int isign,
		int& n_count,
		Double& resid);

//templated MR inverter
template<typename T>
void InvMR_a(const LinearOperator<T>& M,
		const T& source,
		T& target,
		const Real& MRovpar,
		const Real& RsdMR,
		int MaxMR,
		int isign,
		int& n_count,
		Double& resid);

//templated CG inverter
//M^+M phi = chi
template<typename T>
void InvCG_a(const LinearOperator<T>& M,
		const T& source,
		T& target,
		const Real& RsdCG,
		int MaxCG,
		int& n_count,
		Double& resid);

void dslash_deriv(multi1d<LatticeColorMatrix>& F,
		const LatticeDiracFermion& X,
		const LatticeDiracFermion& Y,
		int isign, int cb);

#endif

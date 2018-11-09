#ifndef _MONOMIAL_H_
#define _MONOMIAL_H_
#include "qdp.h"
#include "fieldstate.h"
#include "absdslashm_w.h"
#include "handle.h"
using namespace QDP;
using namespace Chroma;

template<typename P, typename Q>
class AbsMonomial {
	public:
		virtual ~AbsMonomial(){}
		//force computation
		virtual void dsdq(P& F, const Q& s) const = 0;

		// total action
		virtual Double S(const AbsFieldState<P,Q>& s) const = 0;

		// refresh pseudofermion field 
		virtual void refreshInternalFields(const AbsFieldState<P,Q>& field_state) = 0;
};

class WilsonGaugeMonomial: public AbsMonomial<GaugeP, GaugeQ> {
	public:
		~WilsonGaugeMonomial(){}
		WilsonGaugeMonomial(const Real& beta_):beta(beta_){}

		void dsdq(GaugeP& F, const GaugeQ& q) const;
		Double S(const AbsFieldState<GaugeP, GaugeQ>& s) const;
		//no need of refresh for gauge action, leave empty
		void refreshInternalFields(const AbsFieldState<GaugeP, GaugeQ>& s) {}

	private:
		Real beta;
};

class TwoFlavorWilsonFermMonomial: public AbsMonomial<GaugeP, GaugeQ> {
	public:
		~TwoFlavorWilsonFermMonomial(){}
		TwoFlavorWilsonFermMonomial(const Real& Mass_,
				const Real& RsdCG_,
				int MaxCG_):Mass(Mass_), RsdCG(RsdCG_), MaxCG(MaxCG_){}
		void dsdq(GaugeP& F, const GaugeQ& q) const;
		Double S(const AbsFieldState<GaugeP, GaugeQ> & s) const;

		void refreshInternalFields(const AbsFieldState<GaugeP, GaugeQ>& s);

	private:
		void getX(LatticeDiracFermion& X, const GaugeQ& u) const;
		Real Mass;
		Real RsdCG;
		int MaxCG;
		LatticeDiracFermion phi;
};

#endif

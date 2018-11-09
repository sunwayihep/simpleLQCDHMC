#ifndef _HAMILTONIAN_H_
#define _HAMILTONIAN_H_
#include "qdp.h"
#include "monomial.h"
using namespace QDP;

template<typename P, typename Q>
class AbsHamiltonian {
	public:
		virtual ~AbsHamiltonian(){}

		virtual int numMonomials() const = 0;
		virtual const AbsMonomial<P, Q>& getMonomial(int i) const = 0;
		virtual AbsMonomial<P, Q>& getMonomial(int i) = 0;

		// aggregate energy
		// kinetic energy
		virtual Double mesKE(const AbsFieldState<P, Q>& s) const {
			Double KE = norm2(s.getP());
			//return KE*Real(0.5);
			return KE;
			//LatticeColorMatrix tmp = zero;
			//for(int mu=0; mu<Nd; ++mu){
			//	tmp += s.getP()[mu] * s.getP()[mu];
			//}
			//return Double(sum(real(trace(tmp))));

			//multi1d<LatticeDouble> ke_per_site(Nd);
			//for(int mu=0; mu<Nd; ++mu){
			//	ke_per_site[mu] = -Double(4);
			//	ke_per_site[mu] += localNorm2(s.getP()[mu]);
			//}
			//Double KE = zero;
			//for(int mu=0; mu<Nd; ++mu){
			//	KE += sum(ke_per_site[mu]);
			//}
			//return KE;

		}

		//potential energy
		virtual Double mesPE(const AbsFieldState<P,Q>& s) const {
			Double PE;
			PE = getMonomial(0).S(s);
			for(int i=1; i<numMonomials(); ++i){
				PE += getMonomial(i).S(s);
			}
			return PE;
		}

		virtual void mesE(const AbsFieldState<P,Q>& s, Double& KE, Double& PE) const {
			KE = mesKE(s);
			PE = mesPE(s);
		}

		void dsdq(P& F, const Q& s) const {
			P F_tmp;
			getMonomial(0).dsdq(F, s);
			for(int i=1; i<numMonomials(); ++i) {
				(getMonomial(i)).dsdq(F_tmp, s);
				F += F_tmp;
			}
		}

		virtual void refreshInternalFields(const AbsFieldState<P,Q>& s) {
			getMonomial(0).refreshInternalFields(s);
			for(int i=1; i<numMonomials(); ++i) {
				getMonomial(i).refreshInternalFields(s);
			}
		}
};


class QCDHamiltonian: public AbsHamiltonian<GaugeP, GaugeQ> {
	public:
		~QCDHamiltonian(){}
		QCDHamiltonian(multi1d<Handle<AbsMonomial<GaugeP, GaugeQ> > >& m_) {
			monomials.resize(m_.size());
			for(int i=0; i<monomials.size(); ++i) {
				monomials[i] = (m_[i]);
			}
		}
		//using P=multi1d<LatticeColorMatrix>;
		//using Q=multi1d<LatticeColorMatrix>;
		//Double mesKE(const AbsFieldState<P, Q>& s) const {
		//	multi1d<LatticeDouble> ke_per_site(Nd);
		//	for(int mu=0; mu<Nd; ++mu){
		//		ke_per_site[mu] = -Double(4);
		//		ke_per_site[mu] += localNorm2(s.getP()[mu]);
		//	}
		//	Double KE = zero;
		//	for(int mu=0; mu<Nd; ++mu){
		//		KE += sum(ke_per_site[mu]);
		//	}
		//	return KE;
		//}

		int numMonomials() const {
			return monomials.size();
		}
		const AbsMonomial<GaugeP, GaugeQ>& getMonomial(int i) const {
			return *(monomials[i]);
		}
		AbsMonomial<GaugeP, GaugeQ>& getMonomial(int i) {
			return *(monomials[i]);
		}

	private:
		multi1d<Handle<AbsMonomial<GaugeP, GaugeQ> > > monomials;
};

























#endif

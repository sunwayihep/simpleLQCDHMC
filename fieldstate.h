#ifndef _FIELDSTATE_H_
#define _FIELDSTATE_H_
#include "qdp.h"
using namespace QDP;

template<typename P, typename Q>
class AbsFieldState {
	public:
		virtual ~AbsFieldState<P,Q>() {}
		virtual AbsFieldState<P,Q>* clone() const = 0;
		// read
		virtual const P& getP() const = 0;
		virtual const Q& getQ() const = 0;
		//write
		virtual P& getP() = 0;
		virtual Q& getQ() = 0;
};

typedef multi1d<LatticeColorMatrix> GaugeP;
typedef multi1d<LatticeColorMatrix> GaugeQ;

class GaugeFieldState: public AbsFieldState<GaugeP, GaugeQ> {
	public:
		GaugeFieldState(const GaugeP& p_, const GaugeQ& q_) {
			p.resize(Nd);
			q.resize(Nd);
			for(int mu=0; mu<Nd; ++mu) {
				p[mu] = p_[mu];
				q[mu] = q_[mu];
			}
		}
		GaugeFieldState(const GaugeFieldState& s) {
			p.resize(Nd);
			q.resize(Nd);
			for(int mu=0; mu<Nd; ++mu) {
				p[mu] = s.p[mu];
				q[mu] = s.q[mu];
			}
		}
		~GaugeFieldState(){};
		
		GaugeFieldState* clone() const {
			return new GaugeFieldState(*this);
		}

		const GaugeP& getP() const { return p; }
		const GaugeQ& getQ() const { return q; }

		GaugeP& getP() { return p; }
		GaugeQ& getQ() { return q; }

	private:
		GaugeP p;
		GaugeQ q;
};


#endif

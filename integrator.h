#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_
#include "qdp.h"
#include "chromabase.h"
#include "reunit.h"
#include "util/gauge/expmat.h"
#include "fieldstate.h"
#include "hamiltonian.h"

using namespace Chroma;
using namespace QDP;

template<typename P, typename Q>
class AbsIntegrator {
	public:
		virtual ~AbsIntegrator() {}

		// MD integration of length n*\delta_\tau
		virtual void operator()(AbsFieldState<P,Q>& s,
				const Real traj_length) const = 0;
};

template<typename P, typename Q>
class AbsLeapfrogIntegrator: public AbsIntegrator<P, Q> {
	public:
		virtual ~AbsLeapfrogIntegrator() {}
		virtual int getNumSteps() const = 0;
	
		virtual void operator()(AbsFieldState<P,Q>& s,
				const Real traj_length) const {
			int n_steps = getNumSteps();
			Real dt = traj_length / Real(n_steps);
			Real dtby2 = dt / Real(2);

			leapP(s, dtby2);
			leapQ(s, dt);
			for(int i=0; i<n_steps-1; i++){
				leapP(s, dt);
				leapQ(s, dt);
			}
			leapP(s, dtby2);
		}
	
	protected:
		virtual void leapP(AbsFieldState<P,Q>& s,
				const Real dt) const = 0;
		virtual void leapQ(AbsFieldState<P,Q>& s,
				const Real dt) const = 0;
		
};

// projection to traceless anti-hermitian
// forward declaration
void taproj(LatticeColorMatrix& a);

class QCDLeapfrog: public AbsLeapfrogIntegrator<GaugeP, GaugeQ> {
	public:
		~QCDLeapfrog(){}
		QCDLeapfrog(AbsHamiltonian<GaugeP, GaugeQ>& H_, int n_steps_):
			H(H_), n_steps(n_steps_) {}

		int getNumSteps() const { return n_steps; }

	protected:
		void leapP(AbsFieldState<GaugeP, GaugeQ>& s, Real dt) const {
			GaugeP F(Nd);
			H.dsdq(F, s.getQ());

			for(int mu=0; mu<Nd; ++mu) {
				taproj(F[mu]);
				// i*dot(p) = -F
				s.getP()[mu] += dt * F[mu];
			}
		}
		
		void leapQ(AbsFieldState<GaugeP, GaugeQ>& s, Real dt) const {
			LatticeColorMatrix tmp_1;
			LatticeColorMatrix tmp_2;

			for(int mu=0; mu<Nd; mu++) {
				tmp_1 = dt * (s.getP())[mu];
				expmat(tmp_1, EXP_EXACT);

				tmp_2 = tmp_1 * (s.getQ())[mu];
				(s.getQ())[mu] = tmp_2;

				reunit((s.getQ())[mu]);
			}
		}
	private:
		int n_steps;
		AbsHamiltonian<GaugeP, GaugeQ>& H;
};

#endif

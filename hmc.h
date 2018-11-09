#ifndef _HMC_H_
#define _HMC_H_
#include "qdp.h"
#include "chromabase.h"
#include "integrator.h"
using namespace Chroma;
using namespace QDP;

template<typename P, typename Q>
class AbsHMCTrj {
	public:
		virtual ~AbsHMCTrj() {}
		virtual void operator()(AbsFieldState<P,Q>& s, const bool WarmUpP){
			AbsIntegrator<P, Q>& MD = getMDIntegrator();
			AbsHamiltonian<P, Q>& H_MC = getMCHamiltonian();

			refreshP(s);
			H_MC.refreshInternalFields(s);

			Handle<AbsFieldState<P,Q> > s_old(s.clone());

			Double KE_old, PE_old, KE, PE;
			H_MC.mesE(s, KE_old, PE_old);
			MD(s, getMDTrajLength());
			H_MC.mesE(s, KE, PE);
			QDPIO::cout<<"old KE = "<<KE_old<<", old PE = "<<PE_old<<std::endl;
			QDPIO::cout<<"new KE = "<<KE<<", new PE = "<<PE<<std::endl;

			Double DeltaKE = KE - KE_old;
			Double DeltaPE = PE - PE_old;
			Double DeltaH = DeltaKE + DeltaPE;
			Double AccProb = where(DeltaH < 0.0, Double(1), exp(-DeltaH));
			QDPIO::cout<<"DeltaH = "<<DeltaH<<std::endl;
			QDPIO::cout<<"AccProb = "<<AccProb <<std::endl;
			if(!WarmUpP) {
				bool acceptTestResult = acceptReject(DeltaH);
				QDPIO::cout<<"AcceptP = "<<acceptTestResult<<std::endl;
				if(!acceptTestResult) {
					s.getQ() = s_old->getQ();
					s.getP() = s_old->getP();
				}
			}
		}
	protected:
		virtual AbsHamiltonian<P, Q>& getMCHamiltonian() = 0;
		virtual AbsIntegrator<P, Q>& getMDIntegrator() = 0;
		virtual Real getMDTrajLength() const = 0;

		virtual void refreshP(AbsFieldState<P, Q>& state) const = 0;
		virtual bool acceptReject(const Double& DeltaH) const = 0;
};

class QCDHMCTrj : public AbsHMCTrj<GaugeP, GaugeQ> {
	public:
		~QCDHMCTrj(){};
		QCDHMCTrj(Handle<AbsHamiltonian<GaugeP, GaugeQ> > H_,
				Handle<AbsIntegrator<GaugeP, GaugeQ> > integrator_,
				const Real& MD_traj_length_):
			H(H_), the_integrator(integrator_),
			MD_traj_length(MD_traj_length_) {}
	protected:
		AbsHamiltonian<GaugeP, GaugeQ>& getMCHamiltonian() { return *H; }
		AbsIntegrator<GaugeP, GaugeQ>& getMDIntegrator() {
			return *the_integrator;
		}
		Real getMDTrajLength() const {
			return MD_traj_length;
		}
		void refreshP(AbsFieldState<GaugeP, GaugeQ>& state) const {
			for(int mu=0; mu<Nd; ++mu){
				gaussian(state.getP()[mu]);
				state.getP()[mu] *= sqrt(Real(0.5));
				taproj(state.getP()[mu]);
			}
		}
		bool acceptReject(const Double& DeltaH) const {
			bool ret_val;
			if(toBool(DeltaH <= Double(0))) {
				ret_val = true;
			}else{
				Double AccProb = exp(-DeltaH);
				Double uni_dev;
				random(uni_dev);

				if(toBool(uni_dev <= AccProb)) { return true; }
				else{ ret_val = false; }
			}
			return ret_val;
		}
	private:
		Handle<AbsHamiltonian<GaugeP, GaugeQ> > H;
		Handle<AbsIntegrator<GaugeP, GaugeQ> > the_integrator;
		Real MD_traj_length;
};

#endif

#include "aggregate.h"

int main(int argc, char* argv[])
{
	QDP_initialize(&argc, &argv);
	multi1d<int> latt_size(Nd);
	latt_size[0] = 4;
	latt_size[1] = 4;
	latt_size[2] = 4;
	latt_size[3] = 8;

	Layout::setLattSize(latt_size);
	Layout::create();
	
	multi1d<LatticeColorMatrix> initial_q(Nd);
	multi1d<LatticeColorMatrix> initial_p(Nd);
	for(int mu=0; mu<Nd; mu++){
		gaussian(initial_q[mu]);
		reunit(initial_q[mu]);

		gaussian(initial_p[mu]);
		initial_p[mu] *= sqrt(Real(0.5));
		taproj(initial_p[mu]);
	}
	Double plaq;
	MeasPlq(initial_q, plaq);
	QDPIO::cout<<"original plaquette = "<<plaq<<std::endl;
	//create a field state(P, Q)
	GaugeFieldState s(initial_p, initial_q);

	//set up monomials, hamiltonian and integrator
	Real beta = Real(5.4);
	Real Mass = 0.02;
	int MaxCG = 100;
	Real RsdCG = Real(1.0e-8);
	int n_steps = 16;
	Real traj_length = 1;

	multi1d<Handle<AbsMonomial<GaugeP, GaugeQ> > > monomials(2);
	monomials[0] = new WilsonGaugeMonomial(beta);
	monomials[1] = new TwoFlavorWilsonFermMonomial(Mass, RsdCG, MaxCG);

	Handle<AbsHamiltonian<GaugeP, GaugeQ> > H(new QCDHamiltonian(monomials));
	Handle<AbsIntegrator<GaugeP, GaugeQ> > integrator(new QCDLeapfrog(*H, n_steps));

	//create HMC function object
	QCDHMCTrj hmc(H, integrator, traj_length);
	for(int i=0; i<100; i++){
		hmc(s, false);

		Double plaquette;
		MeasPlq(s.getQ(), plaquette);
		QDPIO::cout<<"i= "<<i<<" plaquette= "<<plaquette<<std::endl;
	}


	QDP_finalize();
	exit(0);
}

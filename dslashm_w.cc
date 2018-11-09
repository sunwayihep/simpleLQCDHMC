#include "dslashm_w.h"


void dslash(LatticeDiracFermion& out, const multi1d<LatticeColorMatrix>& u,
		const LatticeDiracFermion& phi, int isign, int cb){
	
	switch(isign){
		case 1:
		{
		//	LatticeHalfFermion tmp, tmp2;
		//	//Dir 0 Forward
		//	tmp[rb[cb]] = spinProjectDir0Minus(phi);
		//	tmp2[rb[cb]] = shift(tmp, FORWARD, 0);
		//	out[rb[cb]] = spinReconstructDir0Minus(u[0]*tmp2);

		//	//Dir 0 Backward
		//	tmp[rb[cb]] = adj(u[0])*spinProjectDir0Plus(phi);
		//	tmp2[rb[cb]] = shift(tmp, BACKWARD, 0);
		//	out[rb[cb]] += spinReconstructDir0Plus(tmp2);

		//	//Dir 1 Forward
		//	tmp[rb[cb]] = spinProjectDir1Minus(phi);
		//	tmp2[rb[cb]] = shift(tmp, FORWARD, 1);
		//	out[rb[cb]] += spinReconstructDir1Minus(u[1]*tmp2);

		//	//Dir 1 Backward
		//	tmp[rb[cb]] = adj(u[1])*spinProjectDir1Plus(phi);
		//	tmp2[rb[cb]] = shift(tmp,BACKWARD, 1);
		//	out[rb[cb]] += spinReconstructDir1Plus(tmp2);

		//	//Dir 2 Forward
		//	tmp[rb[cb]] = spinProjectDir2Minus(phi);
		//	tmp2[rb[cb]] = shift(tmp, FORWARD, 2);
		//	out[rb[cb]] += spinReconstructDir2Minus(u[2]*tmp2);

		//	//Dir 2 Backward
		//	tmp[rb[cb]] = adj(u[2])*spinProjectDir2Plus(phi);
		//	tmp2[rb[cb]] = shift(tmp, BACKWARD, 2);
		//	out[rb[cb]] += spinReconstructDir2Plus(tmp2);

		//	//Dir 3 Forward
		//	tmp[rb[cb]] = spinProjectDir3Minus(phi);
		//	tmp2[rb[cb]] = shift(tmp, FORWARD, 3);
		//	out[rb[cb]] += spinReconstructDir3Minus(u[3]*tmp2);

		//	//Dir 3 Backward
		//	tmp[rb[cb]] = adj(u[3])*spinProjectDir3Plus(phi);
		//	tmp2[rb[cb]] = shift(tmp, BACKWARD, 3);
		//	out[rb[cb]] +=spinReconstructDir3Plus(tmp2);
		
		
      			out[rb[cb]] = spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(phi), FORWARD, 0))
			+ spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(phi), BACKWARD, 0))
			+ spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(phi), FORWARD, 1))
			+ spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(phi), BACKWARD, 1))
			+ spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(phi), FORWARD, 2))
			+ spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(phi), BACKWARD, 2))
			+ spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(phi), FORWARD, 3))
			+ spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(phi), BACKWARD, 3));
		} break;

		case -1:
		{
 		        out[rb[cb]] = spinReconstructDir0Plus(u[0] * shift(spinProjectDir0Plus(phi), FORWARD, 0))
			+ spinReconstructDir0Minus(shift(adj(u[0]) * spinProjectDir0Minus(phi), BACKWARD, 0))
			+ spinReconstructDir1Plus(u[1] * shift(spinProjectDir1Plus(phi), FORWARD, 1))
			+ spinReconstructDir1Minus(shift(adj(u[1]) * spinProjectDir1Minus(phi), BACKWARD, 1))
			+ spinReconstructDir2Plus(u[2] * shift(spinProjectDir2Plus(phi), FORWARD, 2))
			+ spinReconstructDir2Minus(shift(adj(u[2]) * spinProjectDir2Minus(phi), BACKWARD, 2))
			+ spinReconstructDir3Plus(u[3] * shift(spinProjectDir3Plus(phi), FORWARD, 3))
			+ spinReconstructDir3Minus(shift(adj(u[3]) * spinProjectDir3Minus(phi), BACKWARD, 3));
		} break;

	}
}


void M_unprec_wils(LatticeDiracFermion& result,
			const LatticeDiracFermion& phi,
			const multi1d<LatticeColorMatrix>& u,
			int isign,
			const Real Mass){
			
			Real mass_term = Real(Nd) + Mass;
			Real half = Real(0.5);
			result = mass_term * phi;
			LatticeDiracFermion tmp;

			dslash(tmp, u, phi, isign, 0);
			dslash(tmp, u, phi, isign, 1);

			result -= half* tmp;
}
				
void InvUnprecWilsonMR(const LatticeDiracFermion& chi,
			LatticeDiracFermion& psi,
			const Real& MRovpar,
			const Real& RsdMR,
			int MaxMR,
			const multi1d<LatticeColorMatrix>& u,
			const Real& Mass,
			int isign,
			int& n_count,
			Double& resid){

			int k=0;
			Double d;
			DComplex c, a;
			LatticeDiracFermion Mr, r;
	
			QDPIO::cout<<"I'm in MR inverter"<<std::endl;
			M_unprec_wils(Mr, psi, u, isign, Mass);
			r = chi - Mr;
			resid = sqrt(norm2(r));

			QDPIO::cout<<"initial resid = "<<resid<<std::endl;
			while((k < MaxMR) && toBool(resid > RsdMR)){
				++k;
				M_unprec_wils(Mr, r, u, isign, Mass);
				c = innerProduct(Mr, r);
				d = norm2(Mr);
				a = c / d;
				a *= MRovpar;

				psi += a * r;
				r -= a * Mr;
				resid = sqrt(norm2(r));
				QDPIO::cout<<"# "<<k<<" inner residual = "<<resid<<std::endl;	
			}
			QDPIO::cout<<"# of iteration = "<<k<<std::endl;
			QDPIO::cout<<"residual = "<<resid<<std::endl;
} 

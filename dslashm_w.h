#ifndef _DSLASHM_H_
#define _DSLASHM_H_
#include "qdp.h"
using namespace QDP;

void dslash(LatticeDiracFermion& out, const multi1d<LatticeColorMatrix>& u,
		const LatticeDiracFermion& phi, int isign, int cb);	

// Unpreconditioned Wilson Opertor
void M_unprec_wils(LatticeDiracFermion& result,
			const LatticeDiracFermion& phi,
			const multi1d<LatticeColorMatrix>& u,
			int isign,
			const Real Mass);
			


// Minimal Residual Solver
void InvUnprecWilsonMR(const LatticeDiracFermion& chi,
			LatticeDiracFermion& psi,
			const Real& MRovpar,
			const Real& RsdMR,
			int MaxMR,
			const multi1d<LatticeColorMatrix>& u,
			const Real& Mass,
			int isign,
			int& n_count,
			Double& resid);
#endif

#include "integrator.h"
void taproj(LatticeColorMatrix& a) {
	LatticeColorMatrix aux_1 = a;
	a -= adj(aux_1);
	if(Nc>1) {
		LatticeReal tmp = imag(trace(a));
		tmp *= (Real(1)/Real(Nc));
		LatticeColorMatrix aux = cmplx(0, tmp);
		a -= aux;
	}
	a *= (Real(1)/Real(2));
}

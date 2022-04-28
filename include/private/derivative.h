#ifndef DERIVATIVE_H
#define DERIVATIVE_H


# include "osqp.h"
# include "types.h"

c_int adjoint_derivative(c_int m, c_int n, OSQPMatrix *P, OSQPVectorf *q, OSQPMatrix *A, OSQPVectorf *l, OSQPVectorf *u, OSQPVectorf *x, OSQPVectorf *y);


#endif /* ifndef DERIVATIVE_H */

#include "kkt.h"


/**
 * Form square symmetric KKT matrix of the form
 *
 * [P + scalar I,         A';
 * A             -1/scalar I]
 *
 * for the ADMM iterations (flag == 0). Otherwise, if flag == 1, it forms
 *
 * [P + scalar I,         A';
 * A               -scalar I]
 *
 * N.B. Only the upper triangular part is stuffed!
 *
 * @param  P         cost matrix (already just upper triangular part)
 * @param  A         linear constraint matrix
 * @param  scalar    ADMM step scalar (flag == 0), or other scalar (flag == 1)
 * @param  flag      boolean to for which kind of matrix is returned
 * @return           return status flag
 */
csc * form_KKT(const csc * P, const  csc * A, c_float scalar, c_int flag){
    c_int nKKT, nnzKKTmax; // Size, number of nonzeros and max number of nonzeros in KKT matrix
    csc *KKT_trip, *KKT;           // KKT matrix in triplet format and CSC format
    c_int ptr, i, j; // Counters for elements (i,j) and index pointer
    c_int z_P=0, z_KKT=0;   // Counter for total number of elements in P and in KKT

    // Get matrix dimensions
    nKKT = P->m + A->m;

    // Get maximum number of nonzero elements (only upper triangular part)
    nnzKKTmax = P->nzmax +           // Number of elements in P
                P->m +               // Number of elements in scalarI
                A->nzmax +           // Number of nonzeros in A
                A->m;                // Number of elements in -1/scalar I

    // Preallocate KKT matrix in triplet format
    KKT_trip = csc_spalloc(nKKT, nKKT, nnzKKTmax, 1, 1);

    if (!KKT_trip) return OSQP_NULL;  // Failed to preallocate matrix

    // Allocate Triplet matrices
    // P + scalar I
    for (j = 0; j < P->n; j++){ // cycle over columns
        // No elements in column j => add diagonal element scalar
        if (P->p[j] == P->p[j+1]){
            KKT_trip->i[z_KKT] = j;
            KKT_trip->p[z_KKT] = j;
            KKT_trip->x[z_KKT] = scalar;
            z_KKT++;
        }
        for (ptr = P->p[j]; ptr < P->p[j + 1]; ptr++) { // cycle over rows
            // Get current row
            i = P->i[ptr];

            // Add element of P
            KKT_trip->i[z_KKT] = i;
            KKT_trip->p[z_KKT] = j;
            KKT_trip->x[z_KKT] = P->x[z_P];
            if (i == j){ // P has a diagonal element, add scalar
                KKT_trip->x[z_KKT] += scalar;
            }
            z_P++;
            z_KKT++;

            // Add diagonal scalar in case
            if ((i < j) && // Diagonal element not reached
                (ptr + 1 == P->p[j+1])){ // last element of column j

                // Add diagonal element scalar
                KKT_trip->i[z_KKT] = j;
                KKT_trip->p[z_KKT] = j;
                KKT_trip->x[z_KKT] = scalar;
                z_KKT++;
            }
        }
    }


    // A' at top right
    for (j = 0; j < A->n; j++) {  // Cycle over columns of A
        for (ptr = A->p[j]; ptr < A->p[j + 1]; ptr++) {
            KKT_trip->p[z_KKT] = P->m + A->i[ptr];  // Assign column index from row index of A
            KKT_trip->i[z_KKT] = j; // Assign row index from column index of A
            KKT_trip->x[z_KKT] = A->x[ptr];  // Assign A value element
            z_KKT++;
        }
    }

    // flag = 0:  -1/scalar*I at bottom right
    // flag = 1:    -scalar*I at bottom right
    if (!flag) {
        for (j = 0; j < A->m; j++) {
            KKT_trip->i[z_KKT] = j + P->n;
            KKT_trip->p[z_KKT] = j + P->n;
            KKT_trip->x[z_KKT] = -1./scalar;
            z_KKT++;
        }
    } else {
        for (j = 0; j < A->m; j++) {
            KKT_trip->i[z_KKT] = j + P->n;
            KKT_trip->p[z_KKT] = j + P->n;
            KKT_trip->x[z_KKT] = -scalar;
            z_KKT++;
        }
    }


    // Allocate number of nonzeros
    KKT_trip->nz = z_KKT;

    // Convert triplet matrix to csc format
    KKT = triplet_to_csc(KKT_trip);

    // Clean matrix in triplet format and return result
    csc_spfree(KKT_trip);
    return KKT;

}

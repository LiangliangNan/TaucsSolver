#ifndef _TAUCS_UTILITIES_H_
#define _TAUCS_UTILITIES_H_

#include <vector>
#include <map>


struct taucs_ccs_matrix;

namespace TaucsUtil {

	// Assuming nothing about the result (the result is NOT stored symmetric).
	taucs_ccs_matrix* Mul2NonSymmetricMatrices(
		const taucs_ccs_matrix* matA,
		const taucs_ccs_matrix* matB);

	// For usage when it's known that the result is symmetric, like A^double * A.
	taucs_ccs_matrix* Mul2NonSymmMatSymmResult(
		const taucs_ccs_matrix* matA,
		const taucs_ccs_matrix* matB);

	// Computes the transpose of a matrix.
	taucs_ccs_matrix* MatrixTranspose(const taucs_ccs_matrix* mat);

	taucs_ccs_matrix* CreateTaucsMatrixFromColumns(
		const std::vector< std::map<int, double> >& cols, 
		int nRows,
		int flags);

	// Multiplies matA by x and stores the result in b. Assumes all memory has 
	// been allocated and the sizes match; assumes matA is not symmetric!!
	void MulNonSymmMatrixVector(
		const taucs_ccs_matrix* matA,
		const double* x,
		double* b);

	// Adds two vectors vecA and VecB and stores the result in vecResult.
	// Assumes all memory has been allocated and the sizes match!
	void Add2Vectors(
		int n, 
		const double* vecA, 
		const double* vecB, 
		double* vecResult);

	// Copy mat to a new matrix, memory will be allocated during the copy process.
	taucs_ccs_matrix* MatrixCopy(const taucs_ccs_matrix* mat);
};


#endif // _TAUCS_UTILITIES_H_
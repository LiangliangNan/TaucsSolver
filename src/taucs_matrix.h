#ifndef _TAUCS_MATRIX_H_
#define _TAUCS_MATRIX_H_

#include "sparse_matrix.h"


// The class TaucsMatrix is a C++ wrapper around TAUCS' matrix type taucs_ccs_matrix.
// This kind of matrix can be either symmetric or not. Symmetric matrices store only 
// the lower triangle.

struct taucs_ccs_matrix;

class TaucsMatrix : public SparseMatrix
{
public:
	/// Create a square matrix initialized with zeros.
	TaucsMatrix(int dim, bool is_symmetric = false);

	/// Create a rectangular matrix initialized with zeros.
	TaucsMatrix(int rows, int columns, bool is_symmetric = false);

	/// Delete this object and the wrapped TAUCS matrix.
	~TaucsMatrix();

	/// Construct and return the TAUCS matrix wrapped by this object.
	/// Note: the TAUCS matrix returned by this method is valid
	///       only until the next call to set_coef(), add_coef() or get_taucs_matrix().
	const taucs_ccs_matrix* get_taucs_matrix() const;

private:
	/// The actual TAUCS matrix wrapped by this object.
	// This is in fact a COPY of the columns array
	mutable taucs_ccs_matrix* m_matrix;

}; // TaucsMatrix


#endif /*_TAUCS_MATRIX_H_*/
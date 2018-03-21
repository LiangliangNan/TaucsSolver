#include "taucs_matrix.h"

// Taucs is a C library
extern "C" {
#include "taucs.h"
}


	/// Create a square matrix initialized with zeros.
	TaucsMatrix::TaucsMatrix(int dim, bool is_symmetric /* = false*/)	
		: SparseMatrix(dim, is_symmetric)
		, m_matrix(0)
	{
	}

	/// Create a rectangular matrix initialized with zeros.
	TaucsMatrix::TaucsMatrix(int rows, int columns, bool is_symmetric /* = false*/)		
		: SparseMatrix(rows, columns, is_symmetric)
		, m_matrix(0)
	{
	}


	/// Delete this object and the wrapped TAUCS matrix.
	TaucsMatrix::~TaucsMatrix()
	{
		// Delete the the wrapped TAUCS matrix
		if (m_matrix != NULL) {
			taucs_ccs_free(m_matrix);
			m_matrix = NULL;
		}
	}

	/// Construct and return the TAUCS matrix wrapped by this object.
	/// Note: the TAUCS matrix returned by this method is valid
	///       only until the next call to set_coef(), add_coef() or get_taucs_matrix().
	const taucs_ccs_matrix* TaucsMatrix::get_taucs_matrix() const
	{
		if (m_matrix != NULL) {
			taucs_ccs_free(m_matrix);
			m_matrix = NULL;
		}

		// Convert matrix's double type to the corresponding TAUCS constant
		int flags = TAUCS_DOUBLE;

		// We store only the lower triangle of symmetric matrices
		if (m_is_symmetric)
			flags |= TAUCS_TRIANGULAR | TAUCS_SYMMETRIC | TAUCS_LOWER;

		// Compute the number of non null elements in the matrix
		int nb_max_elements = 0;
		for (int col=0; col < m_column_dimension; col++)
			nb_max_elements += m_columns[col].dimension();

		// Create the TAUCS matrix wrapped by this object
		m_matrix = taucs_ccs_create(m_row_dimension, m_column_dimension, nb_max_elements, flags);

		// Fill m_matrix's colptr[], rowind[] and values[] arrays
		// Implementation note:
		// - rowind[] = array of non null elements of the matrix, ordered by columns
		// - values[] = array of row index of each element of rowind[]
		// - colptr[j] is the index of the first element of the column j (or where it
		//   should be if it doesn't exist) + the past-the-end index of the last column
		m_matrix->colptr[0] = 0;
		for (int col=0; col < m_column_dimension; col++)
		{
			// Number of non null elements of the column
			int nb_elements = m_columns[col].dimension();

			// Fast copy of column indices and values
			memcpy(&m_matrix->rowind[m_matrix->colptr[col]], &m_columns[col].m_indices[0], nb_elements*sizeof(int));
			double* taucs_values = (double*) m_matrix->values.v;
			memcpy(&taucs_values[m_matrix->colptr[col]], &m_columns[col].m_values[0],  nb_elements*sizeof(double));

			// Start of next column will be:
			m_matrix->colptr[col+1] = m_matrix->colptr[col] + nb_elements;
		}

		return m_matrix;
	}

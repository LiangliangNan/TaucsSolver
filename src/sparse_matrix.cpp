#include "sparse_matrix.h"
#include <cassert>



/// Create a square matrix initialized with zeros.
SparseMatrix::SparseMatrix(int dim, bool is_symmetric /* = false*/)
{
	assert(dim > 0);

	m_row_dimension     = dim;
	m_column_dimension  = dim;
	m_columns           = new Column[m_column_dimension];
	m_is_symmetric      = is_symmetric;
}

/// Create a rectangular matrix initialized with zeros.
SparseMatrix::SparseMatrix(int rows, int columns, bool is_symmetric /* = false*/)
{
	assert(rows > 0);
	assert(columns > 0);
	if (is_symmetric) {
		assert(rows == columns);
	}

	m_row_dimension     = rows;
	m_column_dimension  = columns;
	m_columns           = new Column[m_column_dimension];
	m_is_symmetric      = is_symmetric;
}

SparseMatrix::~SparseMatrix()
{
	// Delete the columns array
	delete[] m_columns;
	m_columns = NULL;
}

/// Read access to a matrix coefficient.
/// Preconditions:
/// - 0 <= i < row_dimension().
/// - 0 <= j < column_dimension().
double  SparseMatrix::get_coef(int i, int j) const
{
	// For symmetric matrices, we store only the lower triangle
	// => swap i and j if (i, j) belongs to the upper triangle
	if (m_is_symmetric && (j > i))
		std::swap(i, j);

	assert(i < m_row_dimension);
	assert(j < m_column_dimension);
	return m_columns[j].get_coef(i);
}

/// Write access to a matrix coefficient: a_ij <- val.
/// Optimization:
/// For symmetric matrices, SparseMatrix stores only the lower triangle
/// set_coef() does nothing if (i, j) belongs to the upper triangle.
/// Preconditions:
/// - 0 <= i < row_dimension().
/// - 0 <= j < column_dimension().
void SparseMatrix::set_coef(int i, int j, double val)
{
	if (m_is_symmetric && (j > i))
		return;

	assert(i < m_row_dimension);
	assert(j < m_column_dimension);
	m_columns[j].set_coef(i, val);
}

/// Write access to a matrix coefficient: a_ij <- a_ij + val.
/// Optimization:
/// For symmetric matrices, SparseMatrix stores only the lower triangle
/// add_coef() does nothing if (i, j) belongs to the upper triangle.
/// Preconditions:
/// - 0 <= i < row_dimension().
/// - 0 <= j < column_dimension().
void SparseMatrix::add_coef(int i, int j, double val)
{
	if (m_is_symmetric && (j > i))
		return;

	assert(i < m_row_dimension);
	assert(j < m_column_dimension);
	m_columns[j].add_coef(i, val);
}



//////////////////////////////////////////////////////////////////////////



// column{index} <- column{index} + val
void Column::add_coef(int index, double val)
{
	// Search for element in vectors
	std::vector<int>::iterator		index_it;
	std::vector<double>::iterator   value_it;
	for (index_it = m_indices.begin(), value_it = m_values.begin();
		index_it != m_indices.end();
		++index_it, ++value_it)
	{
		if(*index_it == index) {
			*value_it += val;       // +=
			return;
		}
	}

	// Element doesn't exist yet if we reach this point
	m_indices.push_back(index);
	m_values.push_back(val);
}

// column{index} <- val
void Column::set_coef(int index, double val)
{
	// Search for element in vectors
	std::vector<int>::iterator		index_it;
	std::vector<double>::iterator   value_it;
	for (index_it = m_indices.begin(), value_it = m_values.begin();
		index_it != m_indices.end();
		++index_it, ++value_it)
	{
		if(*index_it == index) {
			*value_it = val;        // =
			return;
		}
	}

	// Element doesn't exist yet if we reach this point
	m_indices.push_back(index);
	m_values.push_back(val);
}

// return column{index} (0 by default)
double Column::get_coef(int index) const
{
	// Search for element in vectors
	std::vector<int>::const_iterator	index_it;
	std::vector<double>::const_iterator value_it;
	for (index_it = m_indices.begin(), value_it = m_values.begin();
		index_it != m_indices.end();
		++index_it, ++value_it)
	{
		if(*index_it == index)
			return *value_it;       // return value
	}

	// Element doesn't exist yet if we reach this point
	return 0;
}

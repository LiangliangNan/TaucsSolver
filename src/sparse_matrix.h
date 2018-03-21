#ifndef _SPARST_MATRIX_H_
#define _SPARST_MATRIX_H_


// The class SparseMatrix can be either symmetric or not. 
// Symmetric matrices store only the lower triangle.
// This codes is copied and modified a little from CGAL/Taucs_matrix.h

#include <vector>



/*
* A column of a SparseMatrix. The column is compressed, and stored in the form of
* (a vector of values) + (a vector of indices).
* NOTE: client code should use this class
*/
class Column
{
public:
	// Return the number of elements in the column
	int dimension() const    { return static_cast<int>(m_values.size()); }

	// column{index} <- column{index} + val
	void add_coef(int index, double val);

	// column{index} <- val
	void   set_coef(int index, double val);
	double get_coef(int index) const;

public:
	// (Vector of values) + (vector of indices) (linked)
	std::vector<double>	m_values;
	std::vector<int>	m_indices;
}; // class Column





class SparseMatrix
{
public:
	/// Create a square matrix initialized with zeros.
	SparseMatrix(int dim, bool is_symmetric = false);
	/// Create a rectangular matrix initialized with zeros.
	SparseMatrix(int rows, int columns, bool is_symmetric = false);
	~SparseMatrix();

	/// Return the matrix number of rows
	int row_dimension() const    { return m_row_dimension; }
	/// Return the matrix number of columns
	int column_dimension() const { return m_column_dimension; }

	/// Read access to a matrix coefficient.
	/// Preconditions:
	/// - 0 <= i < row_dimension().
	/// - 0 <= j < column_dimension().
	double  get_coef(int i, int j) const;

	/// Write access to a matrix coefficient: a_ij <- val.
	/// Optimization:
	/// For symmetric matrices, SparseMatrix stores only the lower triangle
	/// set_coef() does nothing if (i, j) belongs to the upper triangle.
	/// Preconditions:
	/// - 0 <= i < row_dimension().
	/// - 0 <= j < column_dimension().
	void set_coef(int i, int j, double val);

	/// Write access to a matrix coefficient: a_ij <- a_ij + val.
	/// Optimization:
	/// For symmetric matrices, SparseMatrix stores only the lower triangle
	/// add_coef() does nothing if (i, j) belongs to the upper triangle.
	/// Preconditions:
	/// - 0 <= i < row_dimension().
	/// - 0 <= j < column_dimension().
	void add_coef(int i, int j, double val);

private:
	/// SparseMatrix cannot be copied (yet)
	SparseMatrix(const SparseMatrix& rhs);
	SparseMatrix& operator=(const SparseMatrix& rhs);

protected:
	// Matrix dimensions
	int     m_row_dimension;
	int		m_column_dimension;

	// Columns array
	Column* m_columns;

	// Symmetric/hermitian?
	bool    m_is_symmetric;

}; // SparseMatrix



#endif /*_SPARST_MATRIX_H_*/
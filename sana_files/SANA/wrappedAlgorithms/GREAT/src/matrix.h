#ifndef MATRIX_H_
#define MATRIX_H_
/*
 *   Copyright (c) 2007 John Weaver
 *
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
 */
template<class T> class Matrix {
public:
	Matrix();
	Matrix(int rows, int columns);
	Matrix(const Matrix<T> &other);
	Matrix<T> & operator=(const Matrix<T> &other);
	~Matrix();
	// all operations except product modify the matrix in-place.
	void resize(int rows, int columns);
	void clear(void);
	T& operator ()(int x, int y);
	T get(int const x, int const y) const;
	void set(int const x, int const y, T const entry);
	int minsize(void) const {
		return ((m_rows < m_columns) ? m_rows : m_columns);
	}
	int columns(void) const {
		return m_columns;
	}
	int rows(void) const {
		return m_rows;
	}
private:
	T **m_matrix;
	int m_rows;
	int m_columns;
};

#ifndef USE_EXPORT_KEYWORD
#include "matrix.cpp"
//#define export /*export*/
#endif

#endif /*MATRIX_H_*/

// Part of OptProject: A basis project for RO problems
// By Sohaib Afifi <me [at] sohaibafifi.com>
// Copyright (c) 2012-2013 University of Technology of Compiegne, France

// Origin
//! \file	SimpleMatrix.hpp
//! \brief	Generic Matrix Class for with fast access and reading/writing capabilities
//! \author	Houssame Yahiaoui
//! \version	1.7
//! \date	April, 11, 2011

#ifndef	__SIMPLEMATRIX_H_
#define	__SIMPLEMATRIX_H_

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

//! \def NO_CSV_READING
/**
 *	When defined, this macro disables reading from .csv files, and permits simpler compilation instructions.
 *	If not defined, the compilation of csvReader.y and csvReader.l are mandatory.
 */
#define	NO_CSV_READING

//! \class SimpleMatrix
/**
 * \brief A generic matrix class.
 *
 * This generic class possess advanced reading and writing capabilities, to import to/export from csv, mtl,... files. It is also able to
 * dynamically change its dimensions.
 * The matrix uses a set of mathematical operations: additions, subtraction, multiplication, scalar multiplication,
 * (bi)stochastic testing, ...
 *
 * \todo Read matrix data from 'mtl' file: Matlab format file for sparse matrix. Only contains cells!=0, with their coordinates.
 */
namespace Utility{
template<class _T>
    class SimpleMatrix
    {
        private:
            _T *_data;
            unsigned _rows, _columns;

            void
            quickSortRow(unsigned row, unsigned beginCol, unsigned endCol);

        public:
            /**
             * \brief Constructor, which creates a vector with the desired size.
             * \param size width of the wanted vector. If omitted, a 1x1 matrix is created
             */
            SimpleMatrix(const unsigned size = 1);

            /**
             * \brief Constructor, which creates a vector with the desired size.
             * \param rows height of the wanted matrix
             * \param cols width of the wanted matrix
             */
            SimpleMatrix(const unsigned rows, const unsigned cols);

            /**
             * \brief Constructor, which creates a vector with the desired size.
             * \param rows height of the wanted matrix
             * \param cols width of the wanted matrix
             * \param initValue initial value for each cell of the new matrix
             */
                    SimpleMatrix(const unsigned rows, const unsigned cols, const _T initValue);

            /**
             * \brief Constructor, which creates a vector with the desired size.
             * \param rows height of the wanted matrix
             * \param cols width of the wanted matrix
             * \param initValuesArray an array of at least rows*cols values of type _T, which will used as initialization values.
             */
                    SimpleMatrix(const unsigned rows, const unsigned cols, const _T * initValuesArray);

            /**
             * \brief Destructor.
             */
            ~SimpleMatrix();

            /**
             * \brief Copy Constructor
             * \param other An object of class SimpleMatrix
             */
            SimpleMatrix(const SimpleMatrix<_T> &other);

            /**
             * \brief Vertical size of the matrix
             * \return Number of rows in Matrix
             */
            inline unsigned
            getRows() const
            {
                return _rows;
            }

            /**
             * \brief Horizontal size of the matrix
             * \return Number of columns in Matrix
             */
            inline unsigned
            getColumns() const
            {
                return _columns;
            }

            /**
             * \brief Global size of the matrix
             * \return Number of columns in Matrix * Number of Rows
             */
            inline unsigned
            getSize() const
            {
                return _rows * _columns;
            }

            /**
             * \brief Reveals if the matrix contains a unique vector of data.
             * \return True, if the matrix contains only one vector of data
             */
            inline bool
            isVector() const
            {
                return (_rows == 1);
            }

            /**
             * \brief Assignment operator
             * \param other An object of class SimpleMatrix
             * \return reference to the actual SimpleMatrix instance
             */
            SimpleMatrix&
            operator=(const SimpleMatrix<_T> &other);

            /**
             * \brief Comparison operator
             * \param other An object of class SimpleMatrix
             * \return true, if both matrix are of the same size, and have same content. false, else.
             */
            bool
            operator==(const SimpleMatrix<_T> &other) const;

            /**
             * \brief Comparison operator
             * \param other An object of class SimpleMatrix
             * \return false, if both matrix are of the same size, and have same content. true, else.
             */
            bool
            operator!=(const SimpleMatrix<_T> &other) const;

            /**
             * \brief Addition Compound Assignment operator
             * \param other An object of class SimpleMatrix, to be added to the actual object
             * \return reference to the actual SimpleMatrix instance, newly modified by the addition.
             */
            SimpleMatrix&
            operator+=(const SimpleMatrix<_T> &other);

            /**
             * \brief Subtraction Compound Assignment operator
             * \param other An object of class SimpleMatrix, to be removed from the actual object
             * \return reference to the actual SimpleMatrix instance, newly modified by the substraction.
             */
            SimpleMatrix&
            operator-=(const SimpleMatrix<_T> &other);

            /**
             * \brief Multiplication by a scalar Compound Assignment operator
             * \param scalar An object of generic type, which will be used for a multiplication of each cell of the matrix.
             * \return reference to the actual SimpleMatrix instance, newly modified by the multiplication.
             */
            SimpleMatrix&
            operator*=(const _T &scalar);

            /**
             * \brief Access (reference) to matrix cell
             * \param row vertical position of wanted cell
             * \param col horizontal position of wanted cell
             * \return reference to  matrix cell
             */
            inline _T&
            operator()(const unsigned row, const unsigned col);

            /**
             * \brief Access to matrix cell
             * \param row vertical position of wanted cell
             * \param col horizontal position of wanted cell
             * \return data copy of matrix cell
             */
            inline _T
            operator()(const unsigned row, const unsigned col) const;

            /**
             * \brief Access (reference) to vector cell
             * \param index position of wanted cell
             * \return reference to wanted vector cell
             */
            inline _T&
            operator()(const unsigned index);

            /**
             * \brief Access to vector cell
             * \param index position of wanted cell
             * \return copy of wanted vector cell
             */
            inline _T
            operator()(const unsigned index) const;

            /**
             * \brief Writes matrix data to textual file
             * \param file output file name
             * \param separator cell separator used  when writing matrix data
             */
            void
            write(const char * file, const char separator = ';') const;

            /**
             * \brief Writes matrix data to output stream
             * \param ost output stream
             * \param separator cell separator used  when writing matrix data
             */
            void
            write(std::ostream &ost, const char separator = ';') const;

            /**
             * \brief Writes sparse matrix data to text file. Only non-zero cells are written.
             * \param file output file name
             * \param commentChar used character to write comments in output stream
             */
            void
            writeSparse(const char * file, const char commentChar = '%') const;

            /**
             * \brief Writes sparse matrix data to output stream. Only non-zero cells are written.
             * \param ost output stream
             * \param commentChar used character to write comments in output stream
             */
            void
            writeSparse(std::ostream &ost, const char commentChar = '%') const;

            /**
             * \brief Adds several, 0 filled, rows to the 'top' of the matrix.
             * \param nb number of rows to be added.
             */
            void
            addRowsTop(const unsigned nb);

            /**
             * \brief Adds several, 0 filled, rows to the 'bottom' of the matrix.
             * \param nb number of rows to be added.
             */
            void
            addRowsBottom(const unsigned nb);

            /**
             * \brief Adds several, 0 filled, columns to the 'right' of the matrix.
             * \param nb number of columns to be added.
             */
            void
            addColumnsRight(const unsigned nb);

            /**
             * \brief Adds several, 0 filled, columns to the 'left' of the matrix.
             * \param nb number of columns to be added.
             */
            void
            addColumnsLeft(const unsigned nb);

            /**
             * \brief Removes several rows from the 'top' of the matrix.
             * \param nb number of rows to be removed.
             */
            void
            removeRowsTop(const unsigned nb);

            /**
             * \brief Removes several rows from the 'bottom' of the matrix.
             * \param nb number of rows to be removed.
             */
            void
            removeRowsBottom(const unsigned nb);

            /**
             * \brief Removes several columns from the 'right' of the matrix.
             * \param nb number of columns to be removed.
             */
            void
            removeColumnsRight(const unsigned nb);

            /**
             * \brief Removes several columns from the 'left' of the matrix.
             * \param nb number of columns to be removed.
             */
            void
            removeColumnsLeft(const unsigned nb);

            /**
             * \brief Sorts the given row, from greater to smaller value (_T must have an ordering operator < or >)
             * \param sortedRow matrix row's index, to be sorted.
             */
            void
            sortRow(const unsigned sortedRow = 0);

            /**
             * \brief Sorts the matrix, row by row, from greater to smaller value (_T must have an ordering operator < or >)
             */
            void
            sortByRow();

            /**
             * \brief	Allows access to a single row of the matrix
             * \param	rn index of the requested row.
             * \return	a SimpleMatrix formed of a single row, containing a copy of the requested matrix row.
             */
            const SimpleMatrix
            row(const unsigned rn) const;

            /**
             * \brief	Allows access to a single column of the matrix
             * \param	cn index of the requested row.
             * \return	a SimpleMatrix formed of a single column, containing a copy of the requested matrix column.
             */
            const SimpleMatrix
            column(const unsigned cn) const;

            /**
             * \brief	Produces the sum of all rows
             * \return	a SimpleMatrix vector, containing the sum of all rows.
             */
            const SimpleMatrix
            rowsSum() const;

            /**
             * \brief	Produces the sum of all columns
             * \return	a SimpleMatrix vector, containing the sum of all columns.
             */
            const SimpleMatrix
            columnsSum() const;

            /**
             * \brief	Builds a transposition of this matrix: each item at position (i,j), is placed at position (j,i).
             * \return	a SimpleMatrix containing the transposition of this actual matrix.
             */
            const SimpleMatrix
            transposition() const;

    };

/**
 *	\brief Writes matrix data to an output stream, separating cells with blank (tabulation).
 *	\param ost output stream to write matrix data
 *	\param m matrix to be wrote
 *	\return reference to output stream
 */
template<class _T>
    std::ostream &
    operator <<(std::ostream &ost, const SimpleMatrix<_T> &m);

/**
 * \brief Addition operator
 * \param one An object of class SimpleMatrix, first operand of the addition
 * \param other An object of class SimpleMatrix, second operand of the addition
 * \return a new SimpleMatrix instance, containing addition of the two matrix
 */
template<class _T>
    const SimpleMatrix<_T>
    operator+(const SimpleMatrix<_T> &one, const SimpleMatrix<_T> &other);

/**
 * \brief Subtraction operator
 * \param one An object of class SimpleMatrix, first operand of the subtraction
 * \param other An object of class SimpleMatrix, second operand of the subtraction
 * \return a new SimpleMatrix instance, containing subtraction of the two matrix
 */
template<class _T>
    const SimpleMatrix<_T>
    operator-(const SimpleMatrix<_T> &one, const SimpleMatrix<_T> &other);

//===========================================================================//
//=================== SubTypes and specific functions =======================//
//===========================================================================//

/**
 *	\brief A matrix containing double values.
 */
typedef SimpleMatrix<double> DoubleMatrix;

/**
 * \brief Reveals if the matrix is stochastic (sum of each row is 1)
 * \param M tested matrix
 * \param epsilon error margin for comparison between probabilities sum and 1
 * \return True, if the matrix is stochastic, Else, false.
 */
bool
isStochastic(DoubleMatrix & M, double epsilon = 1e-15);

/**
 * \brief Reveals if the matrix is bi-stochastic (sum of each row and each column is 1)
 * \param M tested Matrix
 * \param epsilon error margin for comparison between probabilities sum and 1
 * \return True, if the matrix is bi-stochastic. Else, false.
 */
bool
isBiStochastic(DoubleMatrix & M, double epsilon = 1e-15);

/**
 * \brief modifies a matrix (by normalization of rows) to create a stochastic matrix.
 * \param M input matrix, to be normalized
 */
void
normalize(DoubleMatrix & M);



//===========================================================================//
//========================= Implementations =================================//
//===========================================================================//
template<class _T>
    SimpleMatrix<_T>::SimpleMatrix(unsigned size) :
        _rows(1), _columns(size)
    {
        if(size == 0)
            throw std::runtime_error(
                    "Empty matrix (0*X or X*0) creation is impossible.");
        _data = new _T[size];
        memset(_data, 0, size * sizeof(_T));
    }

template<class _T>
    SimpleMatrix<_T>::SimpleMatrix(unsigned rows, unsigned columns) :
        _rows(rows), _columns(columns)
    {
        if(rows == 0 || columns == 0)
            throw std::runtime_error(
                    "Empty matrix (0*X or X*0) creation is impossible.");
        _data = new _T[rows * columns];
        memset(_data, 0, rows * columns * sizeof(_T));
    }

template<class _T>
    SimpleMatrix<_T>::SimpleMatrix(unsigned rows, unsigned columns, const _T initValue) :
        _rows(rows), _columns(columns)
    {
        if(rows == 0 || columns == 0)
            throw std::runtime_error(
                    "Empty matrix (0*X or X*0) creation is impossible.");
        _data = new _T[rows * columns];

        //put ALL cells to initial value
        for (unsigned i = 0; i < _rows; i++)
            for (unsigned j = 0; j < _columns; j++)
                _data[i * _columns + j] = initValue;
    }

template<class _T>
    SimpleMatrix<_T>::SimpleMatrix(unsigned rows, unsigned columns, const _T * initValuesArray) :
        _rows(rows), _columns(columns)
    {
        if(rows == 0 || columns == 0)
            throw std::runtime_error(
                    "Empty matrix (0*X or X*0) creation is impossible.");
        _data = new _T[rows * columns];

        //put ALL cells to initial value from initValuesArray
        memcpy(_data, initValuesArray, _rows * _columns * sizeof(_T));
    }

template<class _T>
    SimpleMatrix<_T>::~SimpleMatrix()
    {
        if(_data != NULL)
            delete[] _data;
    }

template<class _T>
    SimpleMatrix<_T>::SimpleMatrix(const SimpleMatrix<_T> &other) :
        _data(NULL), _rows(other._rows), _columns(other._columns)
    {
        (*this) = other;
    }

template<class _T>
    SimpleMatrix<_T>&
    SimpleMatrix<_T>::operator=(const SimpleMatrix<_T> &other)
    {
        if(&other != this)
        {
            // lose the actual data
            if(_data != NULL)
            {
                delete[] _data;
                _data = NULL;
            }

            // allocate a new one
            _rows = other._rows;
            _columns = other._columns;
            _data = new _T[_rows * _columns];
            memcpy(_data, other._data, _rows * _columns * sizeof(_T));
        }
        return *this;
    }

template<class _T>
    bool
    SimpleMatrix<_T>::operator==(const SimpleMatrix<_T> &other) const
    {
        if(&other != this)
        {
            if((getRows() != other.getRows()) || (getColumns()
                    != other.getColumns()))
                return false;
            for (unsigned i = 0; i < _rows; i++)
                for (unsigned j = 0; j < _columns; j++)
                    if((*this)(i, j) != other(i, j))
                        return false;
            return true;
        }
        else
            return true;
    }

template<class _T>
    bool
    SimpleMatrix<_T>::operator!=(const SimpleMatrix<_T> &other) const
    {
        return !(*this == other);
    }
template<class _T>
    SimpleMatrix<_T>&
    SimpleMatrix<_T>::operator+=(const SimpleMatrix<_T> &other)
    {
        if(getRows() == other.getRows() && getColumns() == other.getColumns())
        {
            for (unsigned i = 0; i < _rows; i++)
                for (unsigned j = 0; j < _columns; j++)
                    (*this)(i, j) += other(i, j);
            return *this;
        }
        else
            throw std::runtime_error("Adding Matrix of different sizes");
    }

template<class _T>
    SimpleMatrix<_T>&
    SimpleMatrix<_T>::operator-=(const SimpleMatrix<_T> &other)
    {
        if(getRows() == other.getRows() && getColumns() == other.getColumns())
        {
            for (unsigned i = 0; i < _rows; i++)
                for (unsigned j = 0; j < _columns; j++)
                    (*this)(i, j) -= other(i, j);
            return *this;
        }
        else
            throw std::runtime_error("Subtracting Matrix of different sizes");
    }

template<class _T>
    SimpleMatrix<_T>&
    SimpleMatrix<_T>::operator*=(const _T &scalar)
    {
        for (unsigned i = 0; i < _rows; i++)
            for (unsigned j = 0; j < _columns; j++)
                (*this)(i, j) *= scalar;
        return *this;
    }

template<class _T>
    const SimpleMatrix<_T>
    operator+(const SimpleMatrix<_T> &one, const SimpleMatrix<_T> &other)
    {
        SimpleMatrix<_T> resultingMatrix(one); // copy one of the matrix ...
        resultingMatrix += other; // ... add other to it, ...
        return (resultingMatrix); // ... and return it!
    }

template<class _T>
    const SimpleMatrix<_T>
    operator-(const SimpleMatrix<_T> &one, const SimpleMatrix<_T> &other)
    {
        SimpleMatrix<_T> resultingMatrix(one); // copy one of the matrix ...
        resultingMatrix -= other; // ... Subtract other from it, ...
        return (resultingMatrix); // ... and return it!
    }

template<class _T>
    inline _T&
    SimpleMatrix<_T>::operator ()(const unsigned row, const unsigned column)
    {
        if(row >= _rows || column >= _columns)
            throw std::runtime_error("Matrix access is out of bound");
        return _data[row * _columns + column];
    }

template<class _T>
    inline _T
    SimpleMatrix<_T>::operator()(const unsigned index) const
    {
        if(!this->isVector())
            std::cerr << "Warning: Accessing Matrix as Vector" << std::endl;
        return _data[index];
    }

template<class _T>
    inline _T&
    SimpleMatrix<_T>::operator ()(const unsigned index)
    {
        if(!this->isVector())
            std::cerr << "Warning: Accessing Matrix as Vector" << std::endl;
        return _data[index];
    }

template<class _T>
    inline _T
    SimpleMatrix<_T>::operator()(const unsigned row, const unsigned column) const
    {
        if(row >= _rows || column >= _columns)
            throw std::runtime_error("Matrix access is out of bound");
        return _data[row * _columns + column];
    }

template<class _T>
    void
    SimpleMatrix<_T>::write(const char * file, const char separator) const
    {
        std::fstream f;

        f.open(file, std::fstream::out | std::fstream::trunc);
        if(f.is_open())
        {
            write(f, separator);
            f.close();
        }
        else
            throw "Unable to open matrix write file";
    }

template<class _T>
    void
    SimpleMatrix<_T>::write(std::ostream &ost, const char separator) const
    {
        for (unsigned i = 0; i < _rows; i++)
        {
            for (unsigned j = 0; j < _columns; j++)
                ost << (*this)(i, j) << separator;
            ost << std::endl;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::writeSparse(const char * file, const char commentChar) const
    {
        std::fstream f;

        f.open(file, std::fstream::out | std::fstream::trunc);
        if(f.is_open())
        {
            writeSparse(f, commentChar);
            f.close();
        }
        else
            throw "Unable to open sparse matrix write file";
    }

template<class _T>
    void
    SimpleMatrix<_T>::writeSparse(std::ostream &ost, const char commentChar) const
    {
        // Write header
        /**
         * \todo Write a real header for every sparse matrix data file.
         */

        // write matrix size
        ost << getRows() << '\t' << getColumns() << std::endl;
        // write actual data. only non-zero data are to be written
        for (unsigned i = 0; i < _rows; i++)
        {
            for (unsigned j = 0; j < _columns; j++)
                if((*this)(i, j) != 0)
                    ost << i << '\t' << j << '\t' << (*this)(i, j) << std::endl;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::addRowsBottom(const unsigned nb)
    {
        if((nb > 0) && ((_rows + nb) * _columns > 0))
        {
            _T *tmp_data = new _T[(_rows + nb) * _columns * sizeof(_T)];
            memset(tmp_data, 0, (_rows + nb) * _columns * sizeof(_T));

            // Cells to be copied are on the top part of the matrix. A global copy should work (faster)
            memcpy(tmp_data, _data, _rows * _columns * sizeof(_T));

            _rows += nb;
            if(_data != NULL)
                delete[] _data;
            _data = tmp_data;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::addRowsTop(const unsigned nb)
    {
        if((nb > 0) && ((_rows + nb) * _columns > 0))
        {
            _T *tmp_data = new _T[(_rows + nb) * _columns * sizeof(_T)];
            memset(tmp_data, 0, (_rows + nb) * _columns * sizeof(_T));

            // Cells to be copied are on the top part of the matrix. A global copy should work (faster)
            memcpy(tmp_data + (nb * _columns), _data, _rows * _columns
                    * sizeof(_T));
            _rows += nb;
            if(_data != NULL)
                delete[] _data;
            _data = tmp_data;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::addColumnsRight(const unsigned nb)
    {
        if((nb > 0) && (_rows * (nb + _columns) > 0))
        {
            _T *tmp_data = new _T[_rows * (_columns + nb) * sizeof(_T)];
            memset(tmp_data, 0, _rows * (_columns + nb) * sizeof(_T));

            for (unsigned i = 0; i < _rows; i++)
                // Cells to be copied are on the same row, a global copy should speed things
                memcpy(tmp_data + i * (_columns + nb), _data + i * _columns,
                        _columns * sizeof(_T));

            _columns += nb;
            if(_data != NULL)
                delete[] _data;
            _data = tmp_data;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::addColumnsLeft(const unsigned nb)
    {
        if((nb > 0) && (_rows * (nb + _columns) > 0))
        {
            _T *tmp_data = new _T[_rows * (_columns + nb) * sizeof(_T)];
            memset(tmp_data, 0, _rows * (_columns + nb) * sizeof(_T));

            for (unsigned i = 0; i < _rows; i++)
                // Cells to be copied are on the same row, a global copy should speed things
                memcpy(tmp_data + i * (_columns + nb) + nb, _data + i
                        * _columns, _columns * sizeof(_T));
            _columns += nb;
            if(_data != NULL)
                delete[] _data;
            _data = tmp_data;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::removeRowsBottom(const unsigned nb)
    {
        if((nb > 0) && ((_rows - nb) * _columns > 0))
        {
            _T *tmp_data = new _T[(_rows - nb) * _columns * sizeof(_T)];

            // Cells to be copied are on the top part of the matrix. A global copy should work (faster)
            memcpy(tmp_data, _data, (_rows - nb) * _columns * sizeof(_T));

            _rows -= nb;
            delete[] _data;
            _data = tmp_data;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::removeRowsTop(const unsigned nb)
    {
        if((nb > 0) && ((_rows - nb) * _columns > 0))
        {
            _T *tmp_data = new _T[(_rows - nb) * _columns * sizeof(_T)];

            // Cells to be copied are on the top part of the matrix. A global copy should work (faster)
            memcpy(tmp_data, _data + (nb * _columns), (_rows - nb) * _columns
                    * sizeof(_T));

            _rows -= nb;
            delete[] _data;
            _data = tmp_data;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::removeColumnsRight(const unsigned nb)
    {
        if((nb > 0) && (_rows * (nb - _columns) > 0))
        {
            _T *tmp_data = new _T[_rows * (_columns - nb) * sizeof(_T)];

            for (unsigned i = 0; i < _rows; i++)
                // Cells to be copied are on the same row, a global copy should speed things
                memcpy(tmp_data + i * (_columns - nb), _data + i * _columns,
                        (_columns - nb) * sizeof(_T));

            _columns -= nb;
            delete[] _data;
            _data = tmp_data;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::removeColumnsLeft(const unsigned nb)
    {
        if((nb > 0) && (_rows * (nb - _columns) > 0))
        {
            _T *tmp_data = new _T[_rows * (_columns - nb) * sizeof(_T)];

            for (unsigned i = 0; i < _rows; i++)
                // Cells to be copied are on the same row, a global copy should speed things
                memcpy(tmp_data + i * (_columns - nb), _data + i * _columns
                        + nb, (_columns - nb) * sizeof(_T));
            _columns -= nb;
            delete[] _data;
            _data = tmp_data;
        }
    }

template<class _T>
    void
    SimpleMatrix<_T>::sortRow(const unsigned sortedRow)
    {
        //  Sort row, using quick sort. elements to be sorted are all the cells in the row.
        quickSortRow(sortedRow, 0, getColumns() - 1);
    }

template<class _T>
    void
    SimpleMatrix<_T>::sortByRow()
    {
        for (unsigned i = 0; i < getRows(); i++)
            sortRow(i);
    }

template<class _T>
    void
    SimpleMatrix<_T>::quickSortRow(unsigned row, unsigned beginCol, unsigned endCol)
    {
        unsigned pivot = beginCol, left = beginCol, right = endCol;
        _T tmp;

        while (left < right)
        {
            if((*this)(row, left) < (*this)(row, right))
            {
                // switch
                tmp = (*this)(row, left);
                (*this)(row, left) = (*this)(row, right);
                (*this)(row, right) = tmp;

                pivot = (pivot == left) ? right : left;
            }
            if(pivot == left)
                right--;
            else
                left++;
        }
        // recursive calls
        if(beginCol + 1 < left)
            quickSortRow(row, beginCol, left - 1);
        if(endCol > right + 1)
            quickSortRow(row, right + 1, endCol);
    }

template<class _T>
    std::ostream &
    operator <<(std::ostream &ost, const SimpleMatrix<_T> &m)
    {
        m.write(ost, '\t');
        return ost;
    }

template<class _T>
    const SimpleMatrix<_T>
    SimpleMatrix<_T>::row(const unsigned rn) const
    {
        SimpleMatrix<_T> resultingMatrix(1, getColumns()); // make a matrix of size 1*Columns
        for (unsigned i = 0; i < getColumns(); i++)
            resultingMatrix(0, i) = (*this)(rn, i);
        return (resultingMatrix);
    }

template<class _T>
    const SimpleMatrix<_T>
    SimpleMatrix<_T>::column(const unsigned cn) const
    {
        SimpleMatrix<_T> resultingMatrix(getRows(), 1); // make a matrix of size rows*1
        for (unsigned i = 0; i < getRows(); i++)
            resultingMatrix(i, 0) = (*this)(i, cn);
        return (resultingMatrix);
    }

template<class _T>
    const SimpleMatrix<_T>
    SimpleMatrix<_T>::rowsSum() const
    {
        SimpleMatrix<_T> resultingMatrix = row(0); // make a matrix of size 1*Columns, and put first row there (for optimisation)
        for (unsigned i = 1; i < getRows(); i++)
        {
            for (unsigned j = 0; j < getColumns(); j++)
                resultingMatrix(0, j) += (*this)(i, j);
        }
        return (resultingMatrix);
    }

template<class _T>
    const SimpleMatrix<_T>
    SimpleMatrix<_T>::columnsSum() const
    {
        SimpleMatrix<_T> resultingMatrix(1, getRows(), 0); // make a matrix of size 1*rows, filled with 0s
        for (unsigned i = 0; i < getRows(); i++)
        {
            for (unsigned j = 0; j < getColumns(); j++)
                resultingMatrix(0, i) += (*this)(i, j);
        }
        return (resultingMatrix);
    }

template<class _T>
    const SimpleMatrix<_T>
    SimpleMatrix<_T>::transposition() const
    {
        SimpleMatrix<_T> resultingMatrix(getColumns(), getRows()); // make a matrix of size columns*rows.
        for (unsigned i = 0; i < getRows(); i++)
        {
            for (unsigned j = 0; j < getColumns(); j++)
                resultingMatrix(i, j) = (*this)(j, i);
        }
        return (resultingMatrix);
    }

} // end namespace
#endif // __SIMPLEMATRIX_H_

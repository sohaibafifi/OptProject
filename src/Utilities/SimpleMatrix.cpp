//! \file	SimpleMatrix.cpp
//! \brief	Generic Matrix Class for with fast access and reading/writing capabilities
//! \author	Houssame Yahiaoui
//! \version	0.7
//! \date	September, 16, 2010

#include "SimpleMatrix.hpp"
#include <iostream>

using namespace Utility;
bool	isStochastic( DoubleMatrix & M, double epsilon)
{
    double sum;

    for (unsigned i=0; i<M.getRows(); i++)
    {
        sum=0;
        for (unsigned j=0; j<M.getColumns(); j++)
        {
            sum += M(i,j);
            if ( sum>1+epsilon )
                return false;
        }
        if ( (sum<=1-epsilon) || (sum>1+epsilon) )
            return false;
    }
    return true;
}

bool	isBiStochastic( DoubleMatrix & M, double epsilon)
{
    if ( isStochastic(M) )
    {
        double sum=0;

        for (unsigned j=0; j<M.getColumns(); j++)
        {
            sum=0;
            for (unsigned i=0; i<M.getRows(); i++)
            {
                sum += M(i,j);
                if ( sum>1+epsilon )
                    return false;
            }
            if ( (sum<=1-epsilon) || (sum>1+epsilon) )
                return false;
        }
        return true;
    }
    else
        return false;
}

void normalize( DoubleMatrix & M )
{
    double sum;

    for (unsigned i=0; i<M.getRows(); i++)
    {
        sum=0;
        for (unsigned j=0; j<M.getColumns(); j++)
            sum += M(i,j);
        for (unsigned j=0; j<M.getColumns(); j++)
            M(i,j)/=sum;
    }
}

namespace Utility{

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
}

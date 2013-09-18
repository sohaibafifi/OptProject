// Part of OptProject: A basis project for RO problems
// By Sohaib Afifi <me [at] sohaibafifi.com>
// Copyright (c) 2012-2013 University of Technology of Compiegne, France

// Origin
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


// Part of OptProject: A basis project for RO problems
// By Sohaib Afifi <me [at] sohaibafifi.com>
// Copyright (c) 2012-2013 University of Technology of Compiegne, France


#pragma once
#include "../Models/ISolution.hpp"
namespace optimization{
    class IOperator{
        public :
            static virtual bool Apply(ISolution & solution);
    };
}// end namespace

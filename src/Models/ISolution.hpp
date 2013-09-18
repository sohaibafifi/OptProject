// Part of OptProject: A basis project for RO problems
// By Sohaib Afifi <me [at] sohaibafifi.com>
// Copyright (c) 2012-2013 University of Technology of Compiegne, France


#pragma once
#include "IProblem.hpp"
namespace optimization {
    class ISolution{
        public :
            virtual double getCost() = 0;
            IProblem * getProblem(){return _problem;}
        private :
            IProblem *_problem;
    };
}// end namespace

// Part of OptProject: A basis project for RO problems
// By Sohaib Afifi <me [at] sohaibafifi.com>
// Copyright (c) 2012-2013 University of Technology of Compiegne, France


#pragma once
#include "../Models/IProblem.hpp"
#include "../Models/ISolution.hpp"
namespace optimization{
    template <class Solution>
    class ISolver {
        public :
            virtual void solve() = 0;
            IProblem * getProblem(){return _problem;}
            Solution getSolution(){return _solution;}
        private :
            IProblem *_problem;
            Solution _solution;
    };
}// end namespace

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

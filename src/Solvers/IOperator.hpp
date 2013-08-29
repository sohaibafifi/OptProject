#pragma once
#include "../Models/ISolution.hpp"
namespace optimization{
    class IOperator{
        public :
            static virtual bool Apply(ISolution & solution);
    };
}// end namespace

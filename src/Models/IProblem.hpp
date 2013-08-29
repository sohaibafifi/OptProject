#pragma once
#include <string>
namespace optimization{
    class IProblem{
        public :
            virtual void readInstance(std::string filepath) = 0;
            virtual void preprocess() = 0;
            virtual std::string getName() = 0;
    };
}// end namespace

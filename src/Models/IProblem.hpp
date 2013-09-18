// Part of OptProject: A basis project for RO problems
// By Sohaib Afifi <me [at] sohaibafifi.com>
// Copyright (c) 2012-2013 University of Technology of Compiegne, France


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

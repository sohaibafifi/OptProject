#include "Algorithms.hpp"
namespace Utility {
    string Algorithms::itos(int i)
    {
        return static_cast<ostringstream*>( &(ostringstream() << i) )->str();;
    }
    string Algorithms::dtos(double i, unsigned precision)
    {
        return static_cast<ostringstream*>( &(ostringstream() << fixed << setprecision(precision)<< i) )->str();;
    }
    void Algorithms::debug(string message) {
#ifdef DEBUG
        cout << "::" << message.c_str() << endl;
#endif
    }

    string Algorithms::stripWhiteSpace(const string & s, const string & ws) {
        string str = s;
        size_t found;
        found = str.find_last_not_of(ws);
        if (found != string::npos) str.erase(found+1);
        else str.clear();
        return str;
    }

    void Algorithms::splitString(string line, vector<string> & split, char sep) {
        split.clear();
        size_t pos;
        do {
            pos = line.find_first_of(sep);
            if (pos != 0)
                split.push_back(line.substr(0, pos));
            line = line.substr(pos + 1);
        } while (pos < line.npos);
    }

    double Algorithms::getVectorLength(vector<double> v) {
        double length = 0.0;
        unsigned int i, n = v.size();
        for(i=0 ; i<n ; i++) length += v[i]*v[i];
        return sqrt(length);
    }

    double Algorithms::getAverage(vector<double> v) {
        double avg = 0.0;
        unsigned int i, n = v.size();
        for(i=0 ; i<n ; i++) avg += v[i];
        return avg/n;
    }

    double Algorithms::getStandardDeviation(vector<double> v) {
        unsigned int i, n = v.size();
        double avg = getAverage(v), sdv = 0.0;
        for(i=0 ; i<n ; i++) sdv += pow(v[i] - avg, 2);
        return sqrt(sdv);
    }

    double Algorithms::getRelativeDeviation(vector<double> v) {
        unsigned int i, n = v.size();
        double avg = getAverage(v),sdv = 0.0;
        for(i=0 ; i<n ; i++) sdv += pow(v[i] - avg, 2);
        return sqrt(sdv)/avg;
    }

    void Algorithms::generateRandomIndices(vector<unsigned int> & indices, const unsigned int & n) {
        indices.clear();
        if (n==0) return;
        vector<unsigned int> rpick(n);
        unsigned int i, ri;
        for (i=0; i<n; i++) rpick[i] = i;
        while (rpick.size()>0) {
            ri = __::irand() % rpick.size();
            indices.push_back(rpick[ri]);
            rpick[ri] = rpick[rpick.size()-1];
            rpick.pop_back();
        }
    }
    MTRand_int32 __::irand;
    MTRand __::drand;
    unsigned long Algorithms::seed = 0;
    unsigned long Algorithms::initializeRandomSeed(unsigned long _seed) {
        seed = _seed == 0 ? time(NULL):_seed;
        srand(seed);
        __::irand.seed(seed);
        __::drand.seed(seed);
        return seed;
    }

    void Algorithms::generateRandomSeed(vector<unsigned int> & seed, unsigned int n) {
        if (n>=RAND_MAX) return;
        set<unsigned int> check;
        unsigned int current;
        seed.clear();
        while (seed.size()<n) {
            current = __::irand();
            if (check.find(current)==check.end()) {
                check.insert(current);
                seed.push_back(current);
            }
        }
    }

    string Algorithms::directory(string path) {
#if defined(WIN32) || defined(_WIN32)
        size_t pos = path.find_last_of("\\");
#else
        size_t pos = path.find_last_of("/");
#endif
        return path.substr(0,pos);
    }

    bool Algorithms::fileExists(string path) {
        std::ifstream my_file(path.c_str());
        return (my_file.good());
    }

    string Algorithms::filename(string path) {
        std::string filename;
#if defined(WIN32) || defined(_WIN32)
        size_t pos = path.find_last_of("\\");
#else
        size_t pos = path.find_last_of("/");
#endif
        if (pos != std::string::npos)
            filename.assign(path.begin() + pos + 1, path.end());
        else
            filename = path;
        return filename.substr(0, filename.find_last_of('.'));
    }
} // end namespace

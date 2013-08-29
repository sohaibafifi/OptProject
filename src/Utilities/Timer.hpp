#pragma once
namespace Utility {


#if defined(_WIN32)
#include <Windows.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>	/* POSIX flags */
#include <time.h>	/* clock_gettime(), time() */
#include <sys/time.h>	/* gethrtime(), gettimeofday() */
#include <sys/times.h>
#include <sys/resource.h>
#if defined(__MACH__) && defined(__APPLE__)
#include <mach/mach.h>
#include <mach/mach_time.h>
#endif

#else
#error "Unable to define getRealTime( ) for an unknown OS."
#endif
    /**
     * @brief The Timer class represents a timer influenced by the cplex IloTimer
     * It works like a stop watch. The timer reports the CPU time.
     */
    class Timer{
        public :
            enum ClockType{ cpu, real};
            /**
             * @brief This constructor creates a timer.
             */
            Timer(){reset(); setClockType(cpu);}

            /**
             * @brief This member function returns the accumulated time, in seconds, since one of these conditions:
             * - the first call of the member function start after construction of the invoking timer;
             * - the most recent call to the member function restart;
             * - a call to reset.
             */
            double getTime() {
                if ( _running ) return (_total + _gettime() - _start);
                else            return (_total);
            }
            /**
             * @brief This member function stops the invoking timer so that it no longer accumulates time.
             */
            double stop(){
                _total = _total + _gettime() - _start;
                _running = false;
                return _total;
            }
            /**
             * @brief This member function makes the invoking timer resume accumulating time. It returns the current clock time.
             */
            double start(){
                _start = _gettime();
                _running = true;
                return _start;
            }
            /**
             * @brief This member function returns the accumulated time, resets the invoking timer to 0.0, and starts the timer again. In other words, the member function restart is equivalent to the member function reset followed by start.
             */
            double restart(){
                double time = getTime();
                reset();
                start();
                return (time);
            }
            /**
             * @brief Decides how computation times are measured (CPU or wall clock)
             */
            void setClockType(ClockType type){_clockType = type;}
            ClockType getClockType( ){return _clockType;}
            /**
             * @brief This member function sets the elapsed time of the invoking timer to 0.0. It also stops the clock.
             */
            void reset(){
                _total   =  0.0;
                _start   = -0.0;
                _running = false;
            }
        private :
            double _gettime(){if(getClockType() == cpu) return getCPUtime();else return getRealTime();}
            double getCPUtime();
            double getRealTime();
            bool _running;
            double _start, _total;
            ClockType _clockType;

    };
} // end namespace

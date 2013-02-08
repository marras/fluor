#ifndef _CORRELATOR_H_
#define _CORRELATOR_H_

#include "common.h"
#include "MersenneTwister.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <vector>
#include "log.h"

#include <assert.h>

// extern MTRand rng;				//random number generator, instantiated in main.cpp

/**********************************
 * Photon Arrival Times structure *
 **********************************/
class PAT {
 public:
  PAT () {b=k=0;}
  ~PAT () {}

  int b;	//the index of the bin in which photon counts have been observed
  int k;	//number of photon count in this bin

  void Zero () {b=k=0;}
};

/*************************
 * Correlator base class *
 *************************/
class Correlator {
  protected:
    int steps_in_a_bin;		//binning resolution for the calculation of the correlation function (how many frames fit in one bin?)
    int B;			//total number of bins
    long long K;		//total number of photon counts
    double bin_size;		//bin size in time units [ps]
    double avg_intensity;	//average number of photons per bin

  public:
    Correlator (const double step_size, const int num_steps, const int _steps_per_bin) {
	if (num_steps%_steps_per_bin!=0) LOG ("*The number of steps (%d) is not an integer multiple of bin size (%d).\n",num_steps,_steps_per_bin);

	this->steps_in_a_bin = _steps_per_bin;
	B = num_steps / steps_in_a_bin;		//save the total number of bins
	bin_size = steps_in_a_bin * step_size;
	K=0;
     }

    virtual ~Correlator () {}

   /***********************************************************
    * Returns the value of the correlation function.          *
    * Arguments: double time - correlation delay time in [ps] *
    ***********************************************************/
    virtual double Calculate_correlation (double time)=0;
};

/********************************************************
 * Correlator implementation using the direct algorithm *
 ********************************************************/
class CorrelatorSimple : public Correlator {
 private:
    int * photons;				//the array of photon counts in each bin

 public:
    CorrelatorSimple (int * distr, const double step_size, const int num_steps, const int _steps_per_bin);
    CorrelatorSimple (const char* fluor_file, const double step_size, const int num_steps, const int _steps_per_bin);

    ~CorrelatorSimple ();

    double Calculate_correlation (double time);
};

class CrossCorrelator : public Correlator {
 private:
    int * photons1, * photons2;				//the array of photon counts in each bin in each trajectory
    long long K2;
    double avg_intensity2;
    
 public:
//    CrossCorrelator (int * distr, const double step_size, const int num_steps, const int _steps_per_bin);
    CrossCorrelator (const char* fluor_file1, const char* fluor_file2, const double step_size, const int num_steps, const int _steps_per_bin);
    ~CrossCorrelator ();

    double Calculate_correlation (double time);
};

#endif //_CORRELATOR_H_

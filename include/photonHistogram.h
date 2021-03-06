#ifndef _PHOTON_HISTOGRAM_H_
#define _PHOTON_HISTOGRAM_H_

#include <fstream>
#include <vector>

#include <assert.h>

/******************************
 * Histogram generating class *
 ******************************/
class Histogram {
 private:
  std::vector<int> traj;
  int res;
  
 public:
   Histogram (int resolution, const char * input_file);
   ~Histogram () {}

   void HistogramInterdetectionTimes (int nbins, const char * output_file);
   void HistogramPhotonCounts (int num_steps, const char * output_file);
};

#endif //_PHOTON_HISTOGRAM_H_

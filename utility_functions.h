#ifndef UTILITY_FUNCTIONS_H
#define UTILITY_FUNCTIONS_H

#include <vector>
#include "TVector3.h"

//A large number of these are simple functions needed to create the distributions
//such as beta decay or a poisson distribution.
//The other half of the code deals with parameterisations that Diego created
//to express the time it takes photons to propagate through the detector.

  namespace utility{

    int poisson(double mean, double draw, double eng);
    double SpectrumFunction(double *x, double *par);
    double fsn(double *x, double *par);
    double Rn_function(double *x, double *par);
    double Scintillation_function(double *t, double *par);
    std::vector<double> GetVUVTime(double distance, int number_photons);
    std::vector<double> GetVisibleTimeOnlyCathode(double t0, int number_photons);
    std::vector<double> GetVisibleTimeFullConfig1(double t0, double tmean, double distance, int number_photons);
    std::vector<double> GetVisibleTimeFullConfig2(double t0, double tmean, double distance, int number_photons);
    double TimingParamReflected(TVector3 ScintPoint, TVector3 OpDetPoint );
    double finter_d(double *x, double *par);
    double LandauPlusExpoFinal(double *x, double *par);
    double finter_r(double *x, double *par);
    double LandauPlusLandauFinal(double *x, double *par);

  }

#endif

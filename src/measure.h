#include "spectrum.h"
#include "resonances.h"

#ifndef MEASURE
#define MEASURE

//class adapter
class Measure {
    const Spectrum& spectr;
    const Carrier& carrier;
public:
    // capture spectr and carrier with resonances for further calculations
    Measure(const Spectrum& spectr, const Carrier& carrier): spectr(spectr), carrier(carrier) {};

    double operator()(double cs) const;
    std::pair<double, double> minmax_interval() const;
    LinearGrid calculate(const LinearGrid& cs_points) const;
};

#endif

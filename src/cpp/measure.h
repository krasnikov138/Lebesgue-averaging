#include "spectrum.h"
#include "resonances.h"

#ifndef MEASURE
#define MEASURE

//class adapter
class Measure {
    const Spectrum& spectr;
    const Carrier& carrier;

    double minS, maxS;
public:
    // capture spectr and carrier with resonances for further calculations
    Measure(const Spectrum& spectr, const Carrier& carrier);

    double operator()(double cs) const;
    LinearGrid calculate(const LinearGrid& cs_points) const;

    double min_cs() const {
    	return minS;
    }
    double max_cs() const {
    	return maxS;
    }
    std::pair<double, double> minmax_interval() const {
    	return {minS, maxS};
    }

    LinearGrid optimal_grid_exp(size_t nodes, double char_size) const;
    LinearGrid optimal_grid_rational(size_t nodes, double char_size) const;

    // pair <energies with same cs, derivatives dEi/dm>
    std::pair<std::vector<double>, std::vector<double>> const_section_params(double cs) const;
};

#endif

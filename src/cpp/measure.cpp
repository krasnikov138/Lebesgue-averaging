#include <algorithm>
#include "measure.h"

double Measure::operator()(double cs) const {
    double mes = 0;
    double e_l, e_r, cs_l, cs_r;

    for (int i = 0; i < carrier.size(); i++) {
        e_l = spectr.get_e(carrier[i].l);
        e_r = spectr.get_e(carrier[i].r);
        cs_l = spectr.get_cs(carrier[i].l);
        cs_r = spectr.get_cs(carrier[i].r);

        if (cs_l < cs && cs_r < cs)
            mes += (e_r - e_l);
        else if (cs_l < cs && cs_r >= cs)
            mes += (e_r - e_l) * (cs - cs_l) / (cs_r - cs_l);
        else if (cs_l >= cs && cs_r < cs)
            mes += (e_r - e_l) * (cs - cs_r) / (cs_l - cs_r);
    }
    return mes;
}

std::pair<double, double> Measure::minmax_interval() const {
    double cs_min = 1.0e9;
    double cs_max = 0.0;
    double cs_l, cs_r;
    for (int i = 0; i < carrier.size(); i++) {
        cs_l = spectr.get_cs(carrier[i].l);
        cs_r = spectr.get_cs(carrier[i].r);
        cs_min = std::min(cs_min, std::min(cs_l, cs_r));
        cs_max = std::max(cs_max, std::max(cs_l, cs_r));
    }
    return std::make_pair(cs_min, cs_max);
}

LinearGrid Measure::calculate(const LinearGrid& cs_points) const {
    LinearGrid res(cs_points.size());
    std::transform(cs_points.begin(), cs_points.end(), res.begin(), 
                   [this](double cs){return (*this)(cs);});

    return res;
}

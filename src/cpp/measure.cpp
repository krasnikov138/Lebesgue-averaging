#include <algorithm>
#include <numeric>
#include "measure.h"

Measure::Measure(const Spectrum& spectr, const Carrier& carrier): 
    spectr(spectr), carrier(carrier) {
    minS = 1.0e10;
    maxS = 0.0;
    double cs_l, cs_r;
    for (int i = 0; i < carrier.size(); i++) {
        cs_l = spectr.get_cs(carrier[i].l);
        cs_r = spectr.get_cs(carrier[i].r);
        minS = std::min(minS, std::min(cs_l, cs_r));
        maxS = std::max(maxS, std::max(cs_l, cs_r));
    }
};

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

LinearGrid Measure::calculate(const LinearGrid& cs_points) const {
    LinearGrid res(cs_points.size());
    std::transform(cs_points.begin(), cs_points.end(), res.begin(), 
                   [this](double cs){return (*this)(cs);});

    return res;
}

LinearGrid Measure::optimal_grid_exp(size_t nodes, double char_length) const {
    LinearGrid grid(nodes + 1);

    for (int i = 0; i <= nodes; ++i) {
        double trig_arg = M_PI * (2 * i + 1) / 4 / (nodes + 1);

        grid[i] = -std::log(
            std::exp(-minS * char_length) * std::pow(std::sin(trig_arg), 2) +
            std::exp(-maxS * char_length) * std::pow(std::cos(trig_arg), 2)
        ) / char_length;
    }
    std::reverse(grid.begin(), grid.end());
    return grid;
}

LinearGrid Measure::optimal_grid_rational(size_t nodes, double char_length) const {
    LinearGrid grid(nodes + 1);

    for (int i = 0; i <= nodes; ++i) {
        double trig_arg = M_PI * (2 * i + 1) / 4 / (nodes + 1);
        double cos_part = std::pow(std::cos(trig_arg), 2);
        double sin_part = std::pow(std::sin(trig_arg), 2);

        grid[i] = (maxS * cos_part + minS * sin_part + minS * maxS * char_length) /
            (1 + (maxS * sin_part + minS * cos_part) * char_length);
    }
    std::reverse(grid.begin(), grid.end());
    return grid;
}

std::pair<std::vector<double>, std::vector<double>>
Measure::const_section_params(double cs) const {
    std::vector<double> energies, derivatives;

    for (int i = 0; i < carrier.size(); i++) {
        double e_l = spectr.get_e(carrier[i].l);
        double e_r = spectr.get_e(carrier[i].r);
        double cs_l = spectr.get_cs(carrier[i].l);
        double cs_r = spectr.get_cs(carrier[i].r);

        if (!(cs_l <= cs && cs < cs_r) && !(cs_r <= cs && cs < cs_l))
            continue;
        double coef = (e_r - e_l) / (cs_r - cs_l);
        double energy = e_l + (cs - cs_l) * coef;

        energies.push_back(energy);
        derivatives.push_back(std::fabs(coef));
    }

    double normalizer = std::accumulate(derivatives.begin(), derivatives.end(), 0.0);
    for (auto& item: derivatives)
        item /= normalizer;

    return {energies, derivatives};
}

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "spectrum.h"

#ifndef RESONANCES
#define RESONANCES

struct Segment {
    double l, r;
    Segment(double l = 0, double r = 0): l(l), r(r) {};
};

using Carrier = std::vector<Segment>;

void set_union_grid(std::vector<Spectrum> &spectra);
int find_dominant_mat(const std::vector<Spectrum>& spectra, std::vector<bool>& mask, double point);
bool is_oscillated(const Spectrum& spectrum, int left_idx, double right);

std::vector<std::vector<Carrier>> build_carriers(std::vector<Spectrum>& spectra, double start_energy = 1e-5, double K = 0.1);

#endif

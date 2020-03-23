#include <tuple>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "spectrum.h"

struct Segment {
    double l, r;
    Segment(double l = 0, double r = 0): l(l), r(r) {};
};

typedef std::vector<Segment> Carrier;

int find_dominant_mat(std::vector<Spectrum> &spectra, std::vector<bool> &mask, int index, double e = -1.0) {
    int index_max = -1;
    double cs_max = 0.0;
    for (int i = 0; i < spectra.size(); i++) {
        if (!mask[i])
            continue;
        double cs;
        if (e > 0 )
            cs = spectra[i].linear(e, index);
        else
            cs = spectra[i].cs[index];
        if (cs_max < cs) {
            index_max = i;
            cs_max = cs;
        }
    }
    return index_max;
}

void set_union_grid(std::vector<Spectrum> &spectra) {
    std::vector<double> grid;
    size_t max_size = 0;
    for (const auto& spec: spectra)
        max_size = (max_size < spec.size)? spec.size : max_size;
    grid.reserve(max_size);

    std::vector<int> offsets(spectra.size(), 0);
    while (true) {
        double e = -1.0;
        for (int i = 0; i < spectra.size(); i++) {
            if (offsets[i] < spectra[i].size && (e < 0 || spectra[i].e[offsets[i]] < e))
                e = spectra[i].e[offsets[i]];
        }
        if (e > 0) {
            grid.push_back(e);
            /* skip points equal to e */
            for (int i = 0; i < spectra.size(); i++) {
                if (offsets[i] < spectra[i].size && is_close(e, spectra[i].e[offsets[i]]))
                    offsets[i]++;
            }
        } else
            break;
    }

    /* count cross sections on union grid */
    for (int i = 0; i < spectra.size(); i++) {
        int offset = 0;
        std::vector<double> cs;
        cs.reserve(grid.size());

        /* before main part */
        while (grid[offset] < spectra[i].e[0]) {
            cs.push_back(0.0);
            offset++;
        }

        /* main part */
        for (int j = 0; j < spectra[i].size - 1; j++) {
            while (grid[offset] < spectra[i].e[j + 1])
                cs.push_back(spectra[i].linear(grid[offset++], j));
        }
        /* after main part */
        for (; offset < grid.size(); offset++)
            cs.push_back(0.0);
        spectra[i].cs = cs;
        spectra[i].e = grid;
        spectra[i].size = grid.size();
    }
}

std::vector<std::vector<Carrier>> build_carriers(std::vector<Spectrum> &spectra, 
    double start_energy = 1e-5, double K = 0.1) {

    std::vector<std::vector<Carrier>> carriers(spectra.size());

    int mat, prev_mat;
    double prev_point;

    int index = 0;
    std::vector<double>& grid = spectra[0].e;
    while (grid[index] < start_energy)
        index++;
    index--;

    double left = start_energy;
    while (index < grid.size() ) {
        assert(grid[index] < left);
        double right = (1 + K / 2) / (1 - K / 2) * left;

        std::vector<Carrier> stage(spectra.size());
        std::vector<bool> mask(spectra.size(), false);
        // count volatility
        for (int i = 0; i < spectra.size(); i++) {
            int j = index;
            double cs_mean = 0, cs_sq_mean = 0;
            while (j < grid.size() && grid[j] < right) {
                cs_mean += spectra[i].cs[j];
                cs_sq_mean += pow(spectra[i].cs[j], 2);
                j++;
            }
            cs_mean /= (j - index);
            cs_sq_mean /= (j - index);
            if (cs_sq_mean / cs_mean - 1 > 0.2)
                mask[i] = true;
        }
        if (std::all_of(mask.begin(), mask.end(), [](bool el){return !el;})) 
            mask = std::vector<bool>(spectra.size(), true);
        prev_mat = find_dominant_mat(spectra, mask, index, left);
        prev_point = spectra[prev_mat].get_point(left, index++) * (1 - 1e-5);
        while (index < grid.size() && grid[index] < right) {
            mat = find_dominant_mat(spectra, mask, index);
            if (mat == -1) {
                index = grid.size();
                break;
            }
            if (mat == prev_mat)
                stage[mat].push_back(Segment(prev_point, index));
            else {
                double ypl = spectra[prev_mat].get_cs(prev_point);
                double ypr = spectra[prev_mat].cs[index];
                double yl = spectra[mat].get_cs(prev_point);
                double yr = spectra[mat].cs[index];

                double intersect = ((yl - ypl) * index + (ypr - yr) * prev_point) / (yl - ypl + ypr - yr);
                assert(intersect >= prev_point);
                assert(intersect <= index);

                stage[prev_mat].push_back(Segment(prev_point, intersect));
                stage[mat].push_back(Segment(intersect, index));
            }
            prev_mat = mat;
            prev_point = index;
            index++;
        }
        // write result in carriers;
        for (int i = 0; i < spectra.size(); i++) {
            if (!stage[i].empty())
                carriers[i].push_back(std::move(stage[i]));
        }
        // if not last point shift index before right for next stage
        if (index < grid.size()) {
            left = right;
            index--;
        }
    }
    return carriers;
}

double measure(double s, Spectrum& spectrum, Carrier& carrier) {
    double mes = 0;
    double e_l, e_r, cs_l, cs_r;

    for (int i = 0; i < carrier.size(); i++) {
        e_l = spectrum.get_energy(carrier[i].l);
        e_r = spectrum.get_energy(carrier[i].r);
        cs_l = spectrum.get_cs(carrier[i].l);
        cs_r = spectrum.get_cs(carrier[i].r);

        if (cs_l < s && cs_r < s)
            mes += (e_r - e_l);
        else if (cs_l < s && cs_r >= s)
            mes += (e_r - e_l) * (s - cs_l) / (cs_r - cs_l);
        else if (cs_l >= s && cs_r < s)
            mes += (e_r - e_l) * (s - cs_r) / (cs_l - cs_r);
    }
    return mes;
}

std::tuple<double, double> cross_section_interval(Spectrum& spectrum, Carrier& carrier) {
    double cs_min = 1.0e9;
    double cs_max = 0.0;
    double cs_l, cs_r;
    for (int i = 0; i < carrier.size(); i++) {
        cs_l = spectrum.get_cs(carrier[i].l);
        cs_r = spectrum.get_cs(carrier[i].r);
        cs_min = std::min(cs_min, std::min(cs_l, cs_r));
        cs_max = std::max(cs_max, std::max(cs_l, cs_r));
    }
    return std::make_tuple(cs_min, cs_max);
}

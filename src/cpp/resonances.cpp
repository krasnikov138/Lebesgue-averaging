#include "resonances.h"

void set_union_grid(std::vector<Spectrum> &spectra) {
    LinearGrid grid;
    size_t max_size = std::max_element(spectra.begin(), spectra.end(), 
        [](const auto& a, const auto& b){return a.size < b.size;})->size;
    grid.reserve(max_size);

    std::vector<size_t> offsets(spectra.size(), 0);
    double e;
    do {
        e = -1.0;
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
    } while (e > 0);

    /* count cross sections on union grid */
    for (auto& spectr: spectra) {
        spectr.set_grid(grid);
    }
}

int find_dominant_mat(const std::vector<Spectrum>& spectra, std::vector<bool>& mask, double point) {
    int index_max = -1;
    double cs_max = 0.0;
    for (int i = 0; i < spectra.size(); i++) {
        if (!mask[i])
            continue;
        double cs = spectra[i].get_cs(point);
        if (cs_max < cs) {
            index_max = i;
            cs_max = cs;
        }
    }
    return index_max;
}

bool is_oscillated(const Spectrum& spectrum, int left_idx, double right) {
    int count = 0;
    double cs_mean = 0, cs_sq_mean = 0;
    for (int j = left_idx; j < spectrum.size && spectrum.e[j] < right; j++) {
        count++;
        cs_mean += spectrum.cs[j];
        cs_sq_mean += pow(spectrum.cs[j], 2);
    }
    cs_mean /= count;
    cs_sq_mean /= count;

    // condition on high volatility
    if (cs_sq_mean > 1.5 * pow(cs_mean, 2))
        return true;
    return false;
}

std::vector<std::vector<Carrier>> 
build_carriers(std::vector<Spectrum>& spectra, double start_energy, double K) {
    std::vector<std::vector<Carrier>> carriers(spectra.size());

    int mat, prev_mat;
    double prev_point;

    int index = 0;
    LinearGrid& grid = spectra[0].e;
    while (grid[index] < start_energy)
        index++;
    index--;

    double left = start_energy;
    while (index < grid.size()) {
        assert(grid[index] < left);
        double right = (1 + K / 2) / (1 - K / 2) * left;

        std::vector<Carrier> stage(spectra.size());
        std::vector<bool> mask(spectra.size(), false);
        std::transform(spectra.begin(), spectra.end(), mask.begin(), 
                       [index, right](auto& spectrum) {return is_oscillated(spectrum, index, right); });
        if (std::all_of(mask.begin(), mask.end(), [](bool el){return !el;})) 
            std::fill(mask.begin(), mask.end(), true);
        prev_mat = find_dominant_mat(spectra, mask, grid.get_point(left, index));
        prev_point = spectra[prev_mat].get_point(left, index++) - 1e-5;
        while (index < grid.size() && grid[index] < right) {
            mat = find_dominant_mat(spectra, mask, index);
            if (mat == -1) {
                index = grid.size();
                break;
            }
            if (mat == prev_mat)
                stage[mat].push_back(Segment(prev_point, index));
            else {
                double ipoint = intersect(spectra[prev_mat], spectra[mat], prev_point, index);
                assert(ipoint >= prev_point);
                assert(ipoint <= index);

                stage[prev_mat].push_back(Segment(prev_point, ipoint));
                stage[mat].push_back(Segment(ipoint, index));
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

#include "spectrum.h"

/* some common functions for approximative comparison */
bool is_close(double a, double b, double rel_tol) {
    if (fabs(a - b) < rel_tol * b) 
        return true;
    return false;
}

bool is_less(double a, double b, double rel_tol) {
    if (a - b < rel_tol * b)
        return true;
    return false;
}

double intersect(const Spectrum& lspec, const Spectrum& rspec, double lpoint, double rpoint) {
    double ypl = lspec.get_cs(lpoint);
    double ypr = lspec.get_cs(rpoint);
    double yl = rspec.get_cs(lpoint);
    double yr = rspec.get_cs(rpoint);

    return ((yl - ypl) * rpoint + (ypr - yr) * lpoint) / (yl - ypl + ypr - yr);
}

double LinearGrid::get_value(double point) const {
    int i = static_cast<int>(point);
    if (i < 0 || i + 1 >= size())
        return 0.0;
    double ratio = point - i;
    return ratio * (*this)[i + 1] + (1 - ratio) * (*this)[i];
}

double LinearGrid::get_point(double value, std::optional<int> hint) const {
    double point;
    int i;
    if (!hint.has_value()) {
        // find point using binary search
        auto bound = std::lower_bound(begin(), end(), value);
        if (bound == begin() || bound == end())
            i = 0;
        else 
            i = std::distance(bound, begin()) - 1;
    } else 
        i = hint.value();
    point = i;
    if (i < 0) return point;
    if (i + 1 < size())
        point += (value - (*this)[i]) / ((*this)[i + 1] - (*this)[i]);
    return point;
}

LinearGrid linspace(double left, double right, size_t points) {
    assert(points > 2);

    LinearGrid grid(points);
    std::generate(grid.begin(), grid.end(), 
        [left, right, points, n = 0]() mutable {
            return left + (right - left) / (points - 1) * (n++);
        });
    return grid;
}

Spectrum::Spectrum(const std::vector<double>& xs, const std::vector<double>& ys):
    concentration(1.0) {
    if (xs.size() != ys.size())
        throw std::runtime_error("Not equal size of input vectors.");
    e = xs; cs = ys; size = xs.size();
}

void Spectrum::load(const std::string& file_name, int mt) {
    EndfFile endf(file_name);
    load(std::get<CrossSectionData>(endf.get_section(3, mt)));
}

void Spectrum::load(CrossSectionData&& endf_data) {
    e = std::move(endf_data.cs.xs);
    cs = std::move(endf_data.cs.ys);

    for (auto& el: cs)
        el *= concentration;
    size = e.size();
}

void Spectrum::save(const std::string& file_name) const {
    std::ofstream output(file_name);
    for (int i = 0; i < size; i++)
        output << e[i] << ' ' << cs[i] << '\n';
}

void Spectrum::set_grid(const LinearGrid& grid) {
    LinearGrid new_cs(grid.size(), 0.0);

    size_t pos = 0;
    for (size_t i = 0; i < size - 1; i++) {
        while (grid[pos] >= e[i] && grid[pos] < e[i + 1]) {
            new_cs[pos] = linear(grid[pos], i);
            ++pos;
        }
    }

    e = grid;
    cs = std::move(new_cs);
    size = grid.size();
}

double Spectrum::resonance_begin(double threshold, int use_points) const {
    double der_mean = 0;
    size_t i = 1;
    while (i < size) {
        double der = fabs((log10(cs[i]) - log10(cs[i - 1])) / (log10(e[i]) - log10(e[i - 1])));
        if (i < use_points) 
            der_mean += der / use_points;
        else if (der > threshold * der_mean)
            break;
        i++;
    }
    return e[i];
}

Spectrum reactions_sum(const std::vector<Spectrum>& reactions) {
    Spectrum result;
    result.e = reactions[0].e;
    result.size = reactions[0].size;
    result.cs.resize(result.size);
    for (auto& reaction: reactions) {
        for (size_t i = 0; i < reaction.size; ++i)
            result.cs[i] += reaction.cs[i];
    }
    return result;
}

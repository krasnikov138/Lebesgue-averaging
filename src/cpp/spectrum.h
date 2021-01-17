#include <cmath>
#include <cassert>
#include <fstream>
#include <vector>
#include <optional>
#include <iterator>
#include <algorithm>

#include "endf.h"

#ifndef SPECTRUM
#define SPECTRUM

bool is_close(double a, double b, double rel_tol=1e-8);
bool is_less(double a, double b, double rel_tol=1e-8);

struct LinearGrid: public std::vector<double> {
    // array which support linearization
    LinearGrid() {};
    LinearGrid(size_t size, double value = 0.0): std::vector<double>(size, value) {};

    double get_value(double point) const;
    double get_point(double value, std::optional<int> hint = {}) const;

    template<typename Iter>
    std::vector<double> get_points(Iter begin, Iter end) const {
        size_t processed = 0;
        size_t dist = std::distance(begin, end);
        std::vector<double> result(dist, -1);
        for (int i = 1; i < this->size(); ++i) {
            if (processed >= dist) break;
            while (processed < dist && (*this)[i] > *(begin + processed)) {
                double val = *(begin + processed);
                result[processed] = get_point(val, i - 1);
                ++processed;
            }
        }
        return result;
    }

    LinearGrid(const std::vector<double>& obj): std::vector<double>(obj) {};
    LinearGrid(std::vector<double>&& obj): std::vector<double>(std::move(obj)) {};
    LinearGrid& operator=(const std::vector<double>& obj) {
        std::vector<double>::operator=(obj);
        return *this;
    };
    LinearGrid& operator=(std::vector<double>&& obj) {
        std::vector<double>::operator=(std::move(obj));
        return *this;
    }
};

LinearGrid linspace(double left, double right, size_t points);

struct Spectrum {
public:
    std::string name;
    size_t size;
    double concentration;
    LinearGrid e, cs;

    Spectrum(double concentration = 1.0): size(0), concentration(concentration) {};
    Spectrum(const std::vector<double>& xs, const std::vector<double>& ys);

    void load(const std::string& file_name, int mt = 1);
    void load(CrossSectionData&& endf_data);
    void save(const std::string& file_name) const;

    // access values via ratio point or energy point
    inline double get_e(double point) const {return e.get_value(point);};
    inline double get_cs(double point) const {return cs.get_value(point);};
    inline double get_point(double energy, std::optional<int> hint = {}) const {
        return e.get_point(energy, hint);
    };
    template<typename Iter>
    std::vector<double> get_points(Iter begin, Iter end) const {
        return e.get_points(begin, end);
    }

    inline double linear(double energy, std::optional<int> hint = {}) const {
        return cs.get_value(e.get_point(energy, hint));
    };

    void set_grid(const LinearGrid& grid);
    double resonance_begin(double threshold = 5.0, int use_points = 10) const;
};

double intersect(const Spectrum& lspec, const Spectrum& rspec, double lpoint, double rpoint);
// this function consider that all Spectrum in array have the same grid
Spectrum reactions_sum(const std::vector<Spectrum>& reactions);

#endif

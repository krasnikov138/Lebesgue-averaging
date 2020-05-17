#include <cmath>
#include <cassert>
#include <fstream>
#include <vector>
#include <optional>
#include <algorithm>

#ifndef SPECTRUM
#define SPECTRUM

bool is_close(double a, double b, double rel_tol=1e-8);
bool is_less(double a, double b, double rel_tol=1e-8);

struct LinearGrid: public std::vector<double> {
    // array which support linearization
    LinearGrid() {};
    LinearGrid(size_t size, double value = 0.0): std::vector<double>(size, value) {};

    double get_value(double point) const;
    double get_point(double value, std::optional<unsigned int> hint = {}) const;

    LinearGrid& operator=(const std::vector<double>& obj);
};

LinearGrid linspace(double left, double right, size_t points);

struct Spectrum {
public:
    std::string name;
    size_t size;
    double concentration;
    LinearGrid e, cs;

    Spectrum(double concentration = 1.0): 
        name("Unnamed"), size(0), concentration(concentration) {};
    Spectrum(std::string name, double concentration = 1.0): 
        name(name), size(0), concentration(concentration) {};

    void load(std::string file_name);
    void save(std::string file_name) const;

    // access values via ratio point or energy point
    inline double get_e(double point) const {return e.get_value(point);};
    inline double get_cs(double point) const {return cs.get_value(point);};
    inline double get_point(double energy, std::optional<unsigned int> hint = {}) const {
        return e.get_point(energy, hint);
    };
    inline double linear(double energy, std::optional<unsigned int> hint = {}) const {
        return cs.get_value(e.get_point(energy, hint));
    };

    void set_grid(const LinearGrid& grid);
    double resonance_begin(double threshold = 5.0, int use_points = 10) const;
};

double intersect(const Spectrum& lspec, const Spectrum& rspec, double lpoint, double rpoint);

#endif

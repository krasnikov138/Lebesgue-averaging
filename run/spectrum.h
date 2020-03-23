#include <vector>
#include <fstream>
#include <cmath>

/* some common funtions for approximative comparison */
bool is_close(double a, double b, double rel_tol=1e-8) {
    if (fabs(a - b) < rel_tol * b) 
        return true;
    return false;
}

bool is_less(double a, double b, double rel_tol=1e-8) {
    if (a - b < rel_tol * b)
        return true;
    return false;
}

struct Spectrum {
public:
    std::string name;
    size_t size;
    std::vector<double> e, cs;

    Spectrum(): name("Unnamed"), size(0) {};
    Spectrum(std::string name): name(name), size(0) {};

    void load(std::string file_name);
    void save(std::string file_name);

    double get_energy(double point);
    double get_point(double energy, int index);
    double get_cs(double point);
    double get_cs(int i, double ratio);
    double linear(double energy, int i);
    double resonance_begin(double threshold, int use_points);
};

void Spectrum::load(std::string file_name) {
    std::ifstream input(file_name);
    double e_val, cs_val;
    while (input) {
        input >> e_val >> cs_val;
        e.push_back(e_val);
        cs.push_back(cs_val);
    }
    size = e.size();
}

void Spectrum::save(std::string file_name) {
    std::ofstream output(file_name);
    for (int i = 0; i < size; i++)
        output << e[i] << ' ' << cs[i] << '\n';
}

double Spectrum::get_energy(double point) {
    int i = static_cast<int>(point);
    if (i + 1 >= size)
        return 0.0;
    double ratio = point - i;
    return ratio * e[i + 1] + (1 - ratio) * e[i];
}

double Spectrum::get_point(double energy, int index) {
    double point = index;
    if (index + 1 < size)
        point += (energy - e[index]) / (e[index + 1] - e[index]);
    return point;
}

double Spectrum::get_cs(int i, double ratio) {
    if (i + 1 >= size)
        return 0.0;
    return ratio * cs[i + 1] + (1 - ratio) * cs[i];
}

double Spectrum::get_cs(double point) {
    int i = static_cast<int>(point);
    double ratio = point - i;
    return get_cs(i, ratio);
}

double Spectrum::linear(double energy, int i) {
    if (i + 1 >= size)
        return 0.0;
    double ratio = (energy - e[i]) / (e[i + 1] - e[i]);
    return get_cs(i, ratio);
}

double Spectrum::resonance_begin(double threshold = 5.0, int use_points = 10) {
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

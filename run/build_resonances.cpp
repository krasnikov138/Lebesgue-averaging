#include <iostream>
#include <iomanip>
#include <vector>
#include <tuple>
#include <fstream>
#include <cmath>
#include <cassert>
#include <filesystem>
#include <algorithm>

/* use convinient lib for parsing command line */
#include <args.hxx>

namespace fs = std::filesystem;

struct Segment {
    double l, r;
    Segment(double l = 0, double r = 0): l(l), r(r) {};
};

typedef std::vector<Segment> Carrier;

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

int main(int argc, char** argv) {
    /* read command lines arguments */
    args::ArgumentParser parser("Splitting the spectrs on carriers of resonances.", 
                                "This goes after the options.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<double> k_arg(parser, "K", "Set K value - relative width of carriers.", {'k', "width"});
    args::ValueFlag<double> start_energy_arg(parser, "start-energy", 
        "Set the left border of first carrier", {'s', "start"});
    args::ValueFlag<std::string> path_arg(parser, "path", 
        "Path for searching input materials / writing output values. Default value is ./", {'p', "path"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    double K = 0.1, start_energy = 1e-5;
    if (k_arg)
        K = args::get(k_arg);
    if (start_energy_arg)
        start_energy = args::get(start_energy_arg);
    fs::path path(".");
    if (path_arg)
        path = fs::path(args::get(path_arg));

    /* load all files from material directory */
    std::vector<Spectrum> spectra;
    fs::path materials_dir = path / "materials";
    for (auto& entry : fs::directory_iterator(materials_dir)) {
        std::string name(entry.path().filename());

        Spectrum sp(name);
        sp.load(entry.path());
        spectra.push_back(sp);
    }

    set_union_grid(spectra);
    std::vector<std::vector<Carrier>> res = build_carriers(spectra, start_energy, K);

    fs::path carriers_dir = path / "carriers";
    // write carriers in files 
    fs::remove_all(carriers_dir);
    for (int i = 0; i < res.size(); i++) {
        fs::create_directories(carriers_dir / spectra[i].name);
        for (int j = 0; j < res[i].size(); j++) {
            std::string fname = std::to_string(j) + ".csv";
            std::ofstream f(carriers_dir / spectra[i].name / fname);
            f << std::setprecision(10);
            for (int k = 0; k < res[i][j].size(); k++) {
                f << spectra[i].get_energy(res[i][j][k].l) << ',';
                f << spectra[i].get_cs(res[i][j][k].l) << ',';
                f << spectra[i].get_energy(res[i][j][k].r) << ',';
                f << spectra[i].get_cs(res[i][j][k].r) << '\n';
            }
        }
    }

    /* calculate measure variable for some cases */
    /*
    const int N = 100;
    fs::path measures_dir = path / "measures";
    fs::remove_all(measures_dir);
    for (int i = 0; i < res.size(); i++) {
        fs::create_directories(measures_dir / spectra[i].name);
        for (int j = 0; j < res[i].size(); j++) {
            std::string fname = std::to_string(j) + ".csv";
            std::ofstream f(measures_dir / spectra[i].name / fname);
            f << std::setprecision(10);
            auto [cs_min, cs_max] = cross_section_interval(spectra[i], res[i][j]);
            for (int k = 0; k < N + 1; k++) {
                double s = cs_min + (cs_max - cs_min) * k / N;
                f << s << ',' << measure(s, spectra[i], res[i][j]) << '\n';
            }
        }
    } 
    */

    return 0;
}

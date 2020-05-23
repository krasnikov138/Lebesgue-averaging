#include <iostream>
#include <iomanip>
#include <map>
#include <filesystem>

#include "resonances.h"
#include "measure.h"
#include "config_parser.h"
#include "endf.h"

/* use convinient lib for parsing command line */
#include <args.hxx>

namespace fs = std::filesystem;

std::vector<Spectrum> load_spectra(const std::string& config_name, const fs::path& path) {
    ConfigParser config_parser;
    auto info = config_parser.load(config_name);

    // load all files from material directory
    std::vector<Spectrum> spectra;

    fs::path std_dir = path / "materials";
    for (auto& [name, values]: info) {
        double concentration;
        std::string cross_section_file;
        if (auto iter = values.find("concentration"); iter != values.end()) {
            try {
                concentration = std::stod(iter->second);
            } catch (std::invalid_argument& exp) {
                throw std::runtime_error("Can't convert '" + iter->second + "' to double in config file.");
            }
        } else 
            throw std::runtime_error("For material '" + name + "' concentration is required.");

        if (auto iter = values.find("cross_section_file"); iter != values.end()) 
            cross_section_file = iter->second;
        else 
            cross_section_file = std_dir / name;

        if (fs::exists(cross_section_file)) {
            Spectrum sp(name, concentration);
            sp.load(cross_section_file);
            spectra.push_back(sp);
        } else
            throw std::runtime_error("File '" + cross_section_file + "' doesn't exist.");
    }

    return spectra;
}

void load_distribution(const std::string& fname) {
    EndfFile endf(fname);

    std::unique_ptr<EndfData> data = endf.get_section(6, 5);
    if (data) {
        auto& distr = *(reinterpret_cast<EnergyAngleData*>(data.get()));
        std::cout << distr;
    }
    exit(0);
}

int main(int argc, char** argv) {
    /* read command lines arguments */
    args::ArgumentParser argparser("Splitting the spectrs on carriers of resonances.");
    args::HelpFlag help(argparser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<double> k_arg(argparser, "K", "Set K value - relative width of carriers.", {'k', "width"});
    args::ValueFlag<double> start_energy_arg(argparser, "start-energy", 
        "Set the left border of first carrier.", {'s', "start"});
    args::ValueFlag<std::string> path_arg(argparser, "path", 
        "Path for searching input materials / writing output values. "
        "Default value is current directory.", {'p', "path"});
    args::ValueFlag<std::string> config_arg(argparser, "config", 
        "Set name of config file. This file includes parameters of materials "
        "like concentration and some others. Default file is 'materials.conf'.", {"config"});
    args::ValueFlag<std::string> distr_fname_arg(argparser, "distribution", 
        "File with energy-angle distribution.", {"distr-fname"});

    try {
        argparser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << argparser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << argparser;
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
    std::string config_name("materials.conf");
    if (config_arg)
        config_name = args::get(config_arg);

    /*
    if (distr_fname_arg)
        load_distribution(args::get(distr_fname_arg));
    */

    auto spectra = load_spectra(config_name, path);
    set_union_grid(spectra);
    std::vector<std::vector<Carrier>> res = build_carriers(spectra, start_energy, K);

    fs::path carriers_dir = path / "carriers";
    // write carriers in files 
    fs::remove_all(carriers_dir);
    for (int i = 0; i < res.size(); i++) {
        // print num of carriers for given material
        std::cout << spectra[i].name << ": " << res[i].size() << '\n';

        fs::create_directories(carriers_dir / spectra[i].name);
        for (int j = 0; j < res[i].size(); j++) {
            std::string fname = std::to_string(j) + ".csv";
            std::ofstream f(carriers_dir / spectra[i].name / fname);
            f << std::setprecision(10);
            for (int k = 0; k < res[i][j].size(); k++) {
                f << spectra[i].get_e(res[i][j][k].l) << ',';
                f << spectra[i].get_cs(res[i][j][k].l) << ',';
                f << spectra[i].get_e(res[i][j][k].r) << ',';
                f << spectra[i].get_cs(res[i][j][k].r) << '\n';
            }
        }
    }

    /* calculate measure variable for some cases */
    const int N = 100;
    fs::path measures_dir = path / "measures";
    fs::remove_all(measures_dir);
    for (int i = 0; i < res.size(); i++) {
        fs::create_directories(measures_dir / spectra[i].name);
        for (int j = 0; j < res[i].size(); j++) {
            std::string fname = std::to_string(j) + ".csv";
            std::ofstream f(measures_dir / spectra[i].name / fname);
            f << std::setprecision(10);
            Measure lebesgue_variable(spectra[i], res[i][j]);
            auto [cs_min, cs_max] = lebesgue_variable.minmax_interval();
            for (int k = 0; k < N + 1; k++) {
                double cs = cs_min + (cs_max - cs_min) * k / N;
                f << cs << ',' << lebesgue_variable(cs) << '\n';
            }
        }
    }

    return 0;
}

#include <iostream>
#include <iomanip>
#include <filesystem>
#include "resonances.cpp"

/* use convinient lib for parsing command line */
#include <args.hxx>

namespace fs = std::filesystem;

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

    return 0;
}

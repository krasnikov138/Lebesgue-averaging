#include <iostream>
#include <iomanip>
#include <map>
#include <iterator>
#include <numeric>
#include <filesystem>

#include "resonances.h"
#include "measure.h"
#include "yaml_config_parser.h"
#include "endf.h"
#include "distribution2d.h"

/* use convinient lib for parsing command line */
#include <args.hxx>

namespace fs = std::filesystem;

struct GroupIndex {
    size_t size;
    std::vector<size_t> offsets;

    GroupIndex() {};
    GroupIndex(const std::vector<size_t>& group_sizes) {
        from_groups(group_sizes);
    }

    void from_groups(const std::vector<size_t>& group_sizes) {
        offsets.resize(group_sizes.size() + 1, 0);
        for (size_t i = 1; i < offsets.size(); ++i)
            offsets[i] = offsets[i - 1] + group_sizes[i - 1];
        size = offsets.back();
    }

    GroupIndex operator*(size_t mul) const {
        GroupIndex result;
        result.size = size * mul;
        result.offsets = offsets;
        for (auto& off: result.offsets)
            off *= mul;
        return result;
    }

    GroupIndex& operator*=(size_t mul) {
        size *= mul;
        for (auto& off: offsets)
            off *= mul;
        return *this;
    }

    inline size_t operator()(size_t i, size_t j) const {
        return offsets[i] + j;
    }

    inline size_t groups() const {
        return offsets.size() - 1;
    }
};

class LebesgueGrid {
    // optimal grid parameters
    size_t nodes;
    double char_length;
    LinearGrid (Measure::*optimize_grid)(size_t, double) const;
    // indexing tools
    GroupIndex main_idx;
    GroupIndex mat_energy_idx;
    GroupIndex mes_energy_idx;

    std::vector<double> optimal_cs, optimal_mes;
    std::vector<double> const_section_energies;
    std::vector<double> const_section_weights;

    std::vector<double> sorted_const_energies;
    std::vector<size_t> sorted_energies_index;

    void sort_energies() {
        sorted_energies_index.resize(const_section_energies.size());
        std::iota(sorted_energies_index.begin(), sorted_energies_index.end(), 0);
        std::sort(sorted_energies_index.begin(), sorted_energies_index.end(),
            [this](size_t a, size_t b){
                return const_section_energies[a] < const_section_energies[b];
            }
        );
        sorted_const_energies.resize(const_section_energies.size());
        for (size_t i = 0; i < sorted_const_energies.size(); ++i)
            sorted_const_energies[i] = const_section_energies[sorted_energies_index[i]];
    }

    std::vector<double> get_ratio_points(const LinearGrid& grid) const {
        std::vector<double> result(const_section_energies.size());
        std::vector<double> ratio_points = grid.get_points(
            sorted_const_energies.cbegin(), sorted_const_energies.cend());
        for (size_t i = 0; i < sorted_energies_index.size(); ++i)
            result[sorted_energies_index[i]] = ratio_points[i];
        return result;
    }
public:
    LebesgueGrid(size_t nodes, double char_length, char grid_type):
        nodes(nodes), char_length(char_length), optimize_grid(nullptr) {
        if (grid_type != 0 && grid_type != 1)
            throw std::runtime_error("Unknown grid type.");
        optimize_grid = (grid_type == 0) ? &Measure::optimal_grid_exp : &Measure::optimal_grid_rational;
    }

    void from_measures(const std::vector<Measure>& measures, const GroupIndex& idx) {
        main_idx = idx * (nodes + 1);
        size_t size = measures.size() * (nodes + 1);

        for (int i = 0; i < idx.groups(); ++i)
            std::cout << idx(i, 0) << ' ' << idx(i + 1, 0) << std::endl;

        optimal_cs.resize(size);
        optimal_mes.resize(size);

        size_t offset = 0;
        std::vector<size_t> mes_energy_groups;
        for (const auto& measure: measures) {
            LinearGrid grid = (measure.*optimize_grid)(nodes, char_length);
            LinearGrid lebvar = measure.calculate(grid);
            std::copy(grid.cbegin(), grid.cend(), optimal_cs.begin() + offset);
            std::copy(lebvar.cbegin(), lebvar.cend(), optimal_mes.begin() + offset);
            for (auto cs: grid) {
                auto [energies, weights] = measure.const_section_params(cs);
                std::copy(energies.begin(), energies.end(), std::back_inserter(const_section_energies));
                std::copy(weights.begin(), weights.end(), std::back_inserter(const_section_weights));
                mes_energy_groups.push_back(energies.size());
            }
            offset += (nodes + 1);
        }
        sort_energies();

        std::vector<size_t> mat_energy_groups;
        for (size_t i = 0; i < idx.groups(); ++i) {
            size_t mat_group_size = 0;
            for (size_t j = main_idx(i, 0); j < main_idx(i + 1, 0); ++j)
                mat_group_size += mes_energy_groups[j];
            mat_energy_groups.push_back(mat_group_size);
        }
        mat_energy_idx.from_groups(mat_energy_groups);
        mes_energy_idx.from_groups(mes_energy_groups);
    }

    std::vector<double> calculate_neutron_escape(const Spectrum& cs) const {
        std::vector<double> result(mes_energy_idx.groups(), 0);
        std::vector<double> points = get_ratio_points(cs.e);
        for (size_t i = 0; i < mes_energy_idx.groups(); ++i) {
            for (size_t j = mes_energy_idx(i, 0); j <= mes_energy_idx(i + 1, 0); ++j)
                result[i] += cs.get_cs(points[j]) * const_section_weights[j];
        }
        return result;
    }

    std::vector<double> calculate_transition_matrix(
        const std::vector<Spectrum>& cs, const std::vector<ProductSubsection>& energy_angle, int legandre_dim = 5) {
        std::vector<Distribution2D> indicatrix(energy_angle.size());
        for (size_t i = 0; i < energy_angle.size(); ++i)
            indicatrix[i] = Distribution2D(energy_angle[i].distr, legandre_dim);

        std::vector<LinearGrid> yields(energy_angle.size());
        for (size_t i = 0; i < energy_angle.size(); ++i)
            yields[i] = LinearGrid(energy_angle[i].energy_yield.ys);

        // calculate all interpolation points
        std::vector<std::vector<double>> cs_points(cs.size());
        for (size_t i = 0; i < cs.size(); ++i)
            cs_points[i] = get_ratio_points(cs[i].e);

        std::vector<std::vector<double>> primary_points(indicatrix.size());
        for (size_t i = 0; i < indicatrix.size(); ++i)
            primary_points[i] = get_ratio_points(indicatrix[i].primary_energies);

        std::vector<std::vector<std::vector<double>>> secondary_points;
        for (size_t i = 0; i < indicatrix.size(); ++i) {
            std::vector<std::vector<double>> temp;
            for (size_t j = 0; j < indicatrix[i].primary_energies.size(); ++j)
                temp.push_back(get_ratio_points(indicatrix[i].secondary_energies[j]));
            secondary_points.push_back(std::move(temp));
        }

        // transitions from groups g -> h
        size_t matrix_size = main_idx.size * main_idx.size * legandre_dim;
        std::vector<double> transition_matrix(matrix_size);
        for (size_t mat = 0; mat < main_idx.groups(); ++mat) {
            for (size_t g = main_idx(mat, 0); g < main_idx(mat + 1, 0); ++g) {
                for (size_t h = 0; h < mes_energy_idx.groups(); ++h) {
                    std::vector<double> transition(legandre_dim);
                    for (size_t i = mes_energy_idx(g, 0); i <= mes_energy_idx(g + 1, 0); ++i) {
                        for (size_t j = mes_energy_idx(h, 0); j <= mes_energy_idx(h + 1, 0); ++j) {
                            double weight = const_section_weights[i] * const_section_weights[j];
                            double cs_value = cs[mat].get_cs(cs_points[mat][i]);
                            double yield = yields[mat].get_value(primary_points[mat][i]);
                            double scalar_coef = weight * cs_value * yield;

                            int primary_pos = static_cast<int>(primary_points[mat][i]);
                            double secondary_left = -1;
                            if (primary_pos >= 0 & primary_pos <= secondary_points[mat].size())
                                secondary_left = secondary_points[mat][primary_pos][j];
                            double secondary_right = -1;
                            if (primary_pos >= 0 & primary_pos <= secondary_points[mat].size() - 1)
                                secondary_right = secondary_points[mat][primary_pos + 1][j];

                            std::vector<double> legandre = indicatrix[mat].interpolate(
                                primary_points[mat][i], secondary_left, secondary_right);
                            transition = sum_vectors_with_coefs(
                                std::move(transition), 1.0, std::move(legandre), scalar_coef);
                        }
                    }
                    size_t offset = (g * main_idx.size + h) * legandre_dim;
                    std::copy(transition.begin(), transition.end(), transition_matrix.begin() + offset);
                }
            }
        }

        return transition_matrix;
    }
};

class NeutronPlatform {
    double K;
    double start_energy;
    std::map<std::string, MaterialInfo> config;

    std::vector<std::string> material_names;

    Spectrum total_cs;
    std::vector<Spectrum> total_cs_components;
    std::vector<Spectrum> cross_sections;
    std::vector<ProductSubsection> energy_angle;

    std::vector<std::vector<Carrier>> carriers;

    GroupIndex idx;
    std::vector<Measure> measures;

    void print_warning(const std::string& name, const std::string& msg) {
        std::cout << "Warning " << name << ": " << msg << std::endl;
    }

    void raise_error(const std::string& name, const std::string& msg) {
        std::string error = "Error " + name + ": " + msg;
        throw std::runtime_error(error);
    }

    Spectrum load_cs(EndfFile& endf, double concentration, int mt, const std::string& name) {
        Spectrum sp(concentration);
        auto data = std::get<CrossSectionData>(endf.get_section(3, mt));

        if (!data.is_linear())
            print_warning(name, "cross section (mt = " + std::to_string(mt) + ") interpolation isn't linear");
        sp.load(std::move(data));
        return sp;
    }

    ProductSubsection load_energy_angle(EndfFile& endf, const MaterialInfo& info, const std::string& name) {
        auto data = std::get<EnergyAngleData>(endf.get_section(6, 5));

        if (info.energy_angle_index >= 0 && info.energy_angle_zap >= 0)
            print_warning(name, "in energy_angle section both index and ZAP are provided - index is used");
        if (info.energy_angle_index < 0 && info.energy_angle_zap < 0)
            raise_error(name, "in energy_angle section neither index nor ZAP are provided");
        if (info.energy_angle_index >= 0) {
            if (info.energy_angle_index >= data.subsections.size())
                raise_error(name, "subsection with given index in energy_angle data isn't found");
            return data.subsections[info.energy_angle_index];
        } else {
            for (auto& sub: data.subsections)
                if (sub.ZAP == info.energy_angle_zap)
                    return sub;
            raise_error(name, "subsection with given ZAP in energy_angle data isn't found");
        }
    }

    void add_material(const std::string& name, const MaterialInfo& info) {
        material_names.push_back(name);

        EndfFile tcs_endf(info.total_section_file);
        total_cs_components.push_back(load_cs(tcs_endf, info.average_concentration, 1, name));
        EndfFile ea_endf(info.energy_angle_file);
        energy_angle.push_back(load_energy_angle(ea_endf, info, name));
        
        EndfFile cs_endf(info.cross_section_file);
        std::vector<Spectrum> reactions;
        if (!info.reactions_mt.size())
            raise_error(name, "reactions section is empty");
        for (auto mt: info.reactions_mt)
            reactions.push_back(load_cs(cs_endf, info.average_concentration, mt, name));
        set_union_grid(reactions);
        cross_sections.push_back(reactions_sum(reactions));
    }

    template<typename Func>
    void general_save(const std::string& path, Func&& saver, bool stat = false) {
        fs::path dirname = fs::path(path);
        fs::remove_all(dirname);
        for (size_t i = 0; i < carriers.size(); i++) {
            if (stat)
                std::cout << material_names[i] << ": " << carriers[i].size() << '\n';
            fs::create_directories(dirname / material_names[i]);
            for (int j = 0; j < carriers[i].size(); j++) {
                std::string fname = std::to_string(j) + ".csv";
                fs::path full_name = dirname / material_names[i] / fname;
                saver(full_name, i, j);
            }
        }
    }
public:
    NeutronPlatform(const std::string& config_fname) {
        config = load_yaml_config(config_fname);
        for (auto& [name, info]: config) {
            std::cout << name << std::endl;
            add_material(name, info);
        }
        set_union_grid(total_cs_components);
        total_cs = reactions_sum(total_cs_components);
    }

    void create_carriers(double start_energy, double K) {
        carriers = build_carriers(total_cs_components, start_energy, K);

        // create indexation for measures / transition matrix
        std::vector<size_t> groups(carriers.size());
        std::transform(carriers.cbegin(), carriers.cend(), groups.begin(), 
            [](const auto& car) {return car.size();});
        idx.from_groups(groups);

        measures.reserve(idx.size);
        for (size_t i = 0; i < carriers.size(); ++i)
            for (size_t j = 0; j < carriers[i].size(); ++j)
                measures.push_back(Measure(total_cs_components[i], carriers[i][j]));
    }

    auto calculate_coefs(size_t nodes, double char_length, char grid_type, int legandre_dim) {
        LebesgueGrid grid(nodes, char_length, grid_type);
        grid.from_measures(measures, idx);
        auto neutron_escape = grid.calculate_neutron_escape(total_cs);
        auto transition_matrix = grid.calculate_transition_matrix(
            total_cs_components, energy_angle, legandre_dim);

        return std::make_pair(std::move(neutron_escape), std::move(transition_matrix));
    }

    static void save_carrier(const std::string& fname, const Spectrum& spectr, const Carrier& carrier) {
        std::ofstream output(fname);
        output << std::setprecision(10);
        for (size_t i = 0; i < carrier.size(); i++) {
            output << spectr.get_e(carrier[i].l) << ',';
            output << spectr.get_cs(carrier[i].l) << ',';
            output << spectr.get_e(carrier[i].r) << ',';
            output << spectr.get_cs(carrier[i].r) << '\n';
        }
    }

    void save_carriers(const std::string& path) {
        general_save(
            fs::path(path) / "carriers", 
            [this](const std::string& fname, size_t i, size_t j){
                NeutronPlatform::save_carrier(fname, total_cs_components[i], carriers[i][j]);
            }, true
        );
    }

    static void save_measure(const std::string& fname, const Measure& mes, const LinearGrid& grid) {
        LinearGrid lebesgue_variable = mes.calculate(grid);
        std::ofstream output(fname);
        output << std::setprecision(10);
        for (size_t i = 0; i < grid.size(); i++)
            output << grid[i] << ',' << lebesgue_variable[i] << '\n';
    }

    void save_measures(const std::string& path, size_t nodes = 10) {
        general_save(
            fs::path(path) / "measures", 
            [nodes, this](const std::string& fname, size_t i, size_t j){
                LinearGrid optimal_grid = measures[idx(i, j)].optimal_grid_exp(nodes, 1.0);
                NeutronPlatform::save_measure(fname, measures[idx(i, j)], optimal_grid);
            }, false
        );
    }
};

void save_vector(const std::string& fname, const std::vector<double>& vec) {
    std::ofstream output(fname, std::ios::binary);
    output.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(double));
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
        "like concentration and some others. Default file is 'materials.yaml'.", {"config"});
    args::ValueFlag<int> nodes_arg(argparser, "nodes",
        "Number of discrete nodes in grid for measure variable in each carrier. Default: 5.", {"nodes"});
    args::ValueFlag<double> char_length_arg(argparser, "char-length",
        "Approximate characteristic length. Default value: 1.0.", {"char-length"});
    args::ValueFlag<int> legandre_arg(argparser, "legandre-dim",
        "Maximum degree of Legandre polinom to describe "
        "energy-angle distribution. Default: 5.", {"legandre-dim"});

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

    int nodes = 5, legandre_dim = 5;
    double K = 0.1, start_energy = 1.0, char_length = 1.0;
    if (k_arg)
        K = args::get(k_arg);
    if (start_energy_arg)
        start_energy = args::get(start_energy_arg);
    if (nodes_arg)
        nodes = args::get(nodes_arg);
    if (legandre_arg)
        legandre_dim = args::get(legandre_arg);
    if (char_length_arg)
        char_length = args::get(char_length_arg);
    fs::path path(".");
    if (path_arg)
        path = fs::path(args::get(path_arg));
    std::string config_name("materials.yaml");
    if (config_arg)
        config_name = args::get(config_arg);

    NeutronPlatform platform(config_name);
    platform.create_carriers(start_energy, K);
    auto [escape, transitions] = platform.calculate_coefs(nodes, char_length, 0, legandre_dim);
    save_vector(path / "neutron_escape.dat", escape);
    save_vector(path / "transition_matrix.dat", transitions);
    platform.save_carriers(path);
    platform.save_measures(path, nodes);

    return 0;
}

#include <iostream>
#include <fstream>
#include <cassert>
#include <memory>
#include <cctype>
#include <cmath>
#include <type_traits>
#include <algorithm>
#include <variant>

#ifndef ENDF
#define ENDF

class EndfData {
    std::string info;

public:
    EndfData(const std::string& data_info): info(data_info) {};

    std::string get_info() const {return info;}
    virtual ~EndfData() = default;
};

struct InterpolationTable {
    std::vector<size_t> boundaries;
    std::vector<unsigned char> interpolation_types;

    std::vector<double> xs, ys;

    bool is_linear() {
        return (interpolation_types.size() == 1 && interpolation_types[0] == 2);
    };
};

struct CrossSectionData: public EndfData {
    // metainfo
    float za, awr;
    // mass difference
    double qm, qi;
    // break up flag
    unsigned char lr;

    InterpolationTable cs;

    CrossSectionData(): EndfData("cross sections") {};
    inline bool is_linear() {
        return cs.is_linear();
    }
};

struct ContinuumDistribution {
    unsigned char angular_repr, lep;

    std::vector<double> primary_energies;
    std::vector<size_t> discrete_energies;

    // interpolation parameters
    std::vector<size_t> boundaries;
    std::vector<unsigned char> interpolation_types;

    std::vector<std::vector<double>> secondary_energies;
    std::vector<std::vector<std::vector<double>>> coefs;
};

struct ProductSubsection {
    // information about reaction product
    float ZAP, atomic_weight_ratio;
    unsigned char product_modifier_flag;

    // distribution law (continuum only supported = 1)
    int law;

    InterpolationTable energy_yield;
    ContinuumDistribution distr;
};

struct EnergyAngleData : public EndfData {
    // metainfo
    float ZA, atomic_weight_ratio;
    unsigned char yield_multiplicity_flag, reference_frame;

    std::vector<ProductSubsection> subsections;
    EnergyAngleData(): EndfData("energy-angle distr") {};
};

//std::ostream& operator<<(std::ostream& out, const CrossSectionData& obj);
std::ostream& operator<<(std::ostream& out, const ProductSubsection& obj);
std::ostream& operator<<(std::ostream& out, const EnergyAngleData& obj);

class EndfFile {
    bool eof;
    bool is_open;
    std::ifstream pfile;

    std::string line;
    unsigned int line_number;
    unsigned int cur_mf, cur_mt;

    std::vector<std::pair<uint16_t, uint16_t>> table_of_content;

    bool next_line(bool allow_stops = true);
    template<typename T = double>
    T get_number(char i);

    class Sentinel {}; Sentinel skip; // to perform skipping values
    template<typename T, typename... Ts>
    void _get_numbers(char pos, T& first, Ts&... args);
    template<typename... Ts>
    void get_numbers(Ts&... args);

    template<typename FirstType, typename SecondType>
    void load_relationship(size_t size, std::vector<FirstType>& seq1, std::vector<SecondType>& seq2);

    void read(InterpolationTable& data);
    void read(CrossSectionData& data);
    void read(ContinuumDistribution& data);
    void read(ProductSubsection& section);
    void read(EnergyAngleData& distr);
public:
    EndfFile() noexcept: is_open(false), eof(false), line_number(0) {};
    EndfFile(const std::string& fname): is_open(false), eof(false), line_number(0) {open(fname);}

    void open(const std::string& fname);
    void close();
    void reset();

    std::vector<std::pair<uint16_t, uint16_t>>& get_table_of_content();
    std::variant<CrossSectionData, EnergyAngleData> get_section(unsigned int mf, unsigned int mt);
};

#endif

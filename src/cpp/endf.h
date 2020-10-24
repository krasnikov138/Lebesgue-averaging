#include <iostream>
#include <fstream>
#include <cassert>
#include <memory>
#include <cctype>
#include <cmath>
#include <type_traits>
#include <algorithm>

#ifndef ENDF
#define ENDF

class EndfData {
    std::string info;

public:
    EndfData(const std::string& data_info): info(data_info) {};

    std::string get_info() const {return info;}
    virtual ~EndfData() = default;
};

struct CrossSectionData: public EndfData {
    // metainfo
    float za, awr;
    // mass difference
    double qm, qi;
    // break up flag
    unsigned char lr;

    // interpolation info
    unsigned int np;
    unsigned char int_type;

    std::vector<double> energies;
    std::vector<double> cross_sections;

    CrossSectionData(): EndfData("cross sections") {};
};

struct SecondaryEnergyTab {
    unsigned char na;
    unsigned int nd, nep;

    std::vector<double> secondary_energies;
    std::vector<std::vector<double>> coefs;
};

struct ContinuumSubsection {
    // information about reaction product
    float zap, awp;
    unsigned char lip;
    // interpolation parameters
    int law;
    unsigned char lang, lep, int_type_primary, int_type_secondary;
    unsigned int np, ne;

    std::vector<double> yields;
    std::vector<double> primary_energies;

    std::vector<SecondaryEnergyTab> distribution;
};

struct EnergyAngleData : public EndfData {
    // metainfo
    float za, awr;
    unsigned char jp, lct;

    std::vector<ContinuumSubsection> subsections;
    EnergyAngleData(): EndfData("energy-angle distr") {};
};

std::ostream& operator<<(std::ostream& out, const CrossSectionData& obj);
std::ostream& operator<<(std::ostream& out, const ContinuumSubsection& obj);
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
    double get_number(char i);

    class Sentinel {}; Sentinel skip; // to perform skipping values
    template<typename T, typename... Ts>
    void _get_numbers(char pos, T& first, Ts&... args);
    template<typename... Ts>
    void get_numbers(Ts&... args);

    void read(CrossSectionData& data);
    void read(SecondaryEnergyTab& data);
    void read(ContinuumSubsection& section);
    void read(EnergyAngleData& distr);
public:
    EndfFile() noexcept: is_open(false), eof(false), line_number(0) {};
    EndfFile(const std::string& fname): is_open(false), eof(false), line_number(0) {open(fname);}

    void open(const std::string& fname);
    void close();
    void reset();

    std::vector<std::pair<uint16_t, uint16_t>>& get_table_of_content();
    std::unique_ptr<EndfData> get_section(unsigned int mf, unsigned int mt);
};

#endif

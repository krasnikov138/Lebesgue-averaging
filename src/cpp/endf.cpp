#include "endf.h"

void EndfFile::open(const std::string& fname) {
    pfile.open(fname);
    if (pfile.is_open())
        is_open = true;
    else
        throw std::runtime_error("Can't open endf file.");
    table_of_content.clear();
}

void EndfFile::close() {
    if (is_open) {
        pfile.close();
        is_open = false;
    }
}

void EndfFile::reset() {
    pfile.clear();
    pfile.seekg(0);
    eof = false;
}

bool EndfFile::next_line(bool allow_stops) {
    line_number++;
    eof = !static_cast<bool>(std::getline(pfile, line));
    if (!allow_stops && eof)
        throw std::runtime_error("Endf file suddenly end. Stop parsing.");
    if (!eof) {
        cur_mf = std::stoi(line.substr(70, 2));
        cur_mt = std::stoi(line.substr(72, 3));
    }
    return !eof;
}

template<typename T = double>
T EndfFile::get_number(char i) {
    if (eof || i > 5 || i < 0)
        throw std::runtime_error("Endf parsing error, line: " + std::to_string(line_number));
    std::string snum = line.substr(i * 11, 11);
    // check empty string
    if (std::all_of(snum.begin(), snum.end(), [](char ch) {return std::isspace(ch);}))
        throw std::runtime_error("Endf parsing error, line: " + std::to_string(line_number) + 
                                 ", cell number: " + std::to_string((int) i));

    double value;
    auto pos = snum.find_last_of("-+");
    if (pos == std::string::npos || pos == 0)
        value = std::stod(snum);
    else {
        value = std::stod(snum.substr(0, pos));
        int order = std::stoi(snum.substr(pos));

        value *= pow(10.0, order);
    }
    return value;
}

template<typename T, typename... Ts>
void EndfFile::_get_numbers(char pos, T& first, Ts&... args) {
    if constexpr (!std::is_same_v<Sentinel, T>)
        first = get_number<T>(pos);
    if constexpr(sizeof...(args) > 0)
        _get_numbers(++pos, args...);
}

template<typename... Ts>
void EndfFile::get_numbers(Ts&... args) {
    _get_numbers(0, args...);
}

template<typename FirstType, typename SecondType>
void EndfFile::load_relationship(size_t size, std::vector<FirstType>& seq1, std::vector<SecondType>& seq2) {
    seq1.reserve(size);
    seq2.reserve(size);

    size_t count = 0;
    while (count < 2 * size) {
        if (count % 6 == 0) next_line(false);
        seq1.push_back(get_number<FirstType>(count++ % 6));
        seq2.push_back(get_number<SecondType>(count++ % 6));
    }
}

void EndfFile::read(InterpolationTable& data) {
    size_t ranges, num_points;
    get_numbers(skip, skip, skip, skip, ranges, num_points);
    load_relationship(ranges, data.boundaries, data.interpolation_types);
    load_relationship(num_points, data.xs, data.ys);
}

void EndfFile::read(CrossSectionData& data) {
    get_numbers(data.za, data.awr);
    next_line(false);
    get_numbers(data.qm, data.qi, skip, data.lr);
    read(data.cs);
}

void EndfFile::read(EnergyAngleData& data) {
    get_numbers(data.ZA, data.atomic_weight_ratio, data.yield_multiplicity_flag, data.reference_frame);

    auto num_of_sections = static_cast<size_t>(get_number(4));
    for (size_t i = 0; i < num_of_sections; i++) {
        ProductSubsection section;
        next_line(false);
        read(section);
        data.subsections.push_back(std::move(section));
    }
}

void EndfFile::read(ProductSubsection& section) {
    get_numbers(section.ZAP, section.atomic_weight_ratio, section.product_modifier_flag, section.law);
    // by far only continuum subsection is supported
    assert(section.law == 1);
    read(section.energy_yield);

    next_line(false);
    read(section.distr);
}

void EndfFile::read(ContinuumDistribution& data) {
    size_t ranges, ne;
    get_numbers(skip, skip, data.angular_repr, data.lep, ranges, ne);
    load_relationship(ranges, data.boundaries, data.interpolation_types);

    data.primary_energies.resize(ne);
    data.discrete_energies.resize(ne);
    data.secondary_energies.resize(ne);
    data.coefs.resize(ne);

    for (size_t i = 0; i < ne; ++i) {
        next_line(false);
        size_t na, nw, nep;
        get_numbers(skip, data.primary_energies[i], data.discrete_energies[i], na, nw, nep);
        assert(nw == nep * (na + 2));

        size_t pos = 0;
        while (pos < nw) {
            if (pos % 6 == 0) next_line(false);
            data.secondary_energies[i].push_back(get_number(pos++ % 6));

            std::vector<double> legendre_coefs;
            legendre_coefs.reserve(na + 1);
            for (unsigned int j = 0; j < na + 1; ++j) {
                if (pos % 6 == 0) next_line(false);
                legendre_coefs.push_back(get_number(pos++ % 6));
            }
            data.coefs[i].push_back(std::move(legendre_coefs));
        }
    }
}

std::vector<std::pair<uint16_t, uint16_t>>& EndfFile::get_table_of_content() {
    if (!table_of_content.size()) {
        std::pair<uint16_t, uint16_t> prev;
        do {
            next_line();
            if (cur_mf == 0 || cur_mt == 0)
                continue;
            if (cur_mf != prev.first || cur_mt != prev.second) {
                prev = std::make_pair(cur_mf, cur_mt);
                table_of_content.push_back(prev);
            }
        } while (!eof);
    }
    reset();
    return table_of_content;
}

std::variant<CrossSectionData, EnergyAngleData> EndfFile::get_section(unsigned int mf, unsigned int mt) {
    if (!is_open)
        throw std::runtime_error("File isn't opened.");
    if (mf != 3 && mf != 6)
        throw std::runtime_error("Sections 3 and 6 are only supported.");
    do {
        next_line();
    } while (!eof && ((cur_mf < mf) || ((cur_mf == mf) && (cur_mt < mt))));
    if (eof || cur_mf != mf || cur_mt != mt) {
        reset();
        throw std::runtime_error("Can't find given section in file.");
    }

    std::variant<CrossSectionData, EnergyAngleData> result;
    if (mf == 3) {
        auto data = CrossSectionData();
        read(data);
        result = std::move(data);
    } else if (mf == 6) {
        auto data = EnergyAngleData();
        read(data);
        result = std::move(data);
    }
    reset();
    return result;
}

std::ostream& operator<<(std::ostream& out, const ProductSubsection& obj) {
    out << "\tProduct charge/mass parameters (ZAP) = " << (int) obj.ZAP << std::endl;
    out << "\tAtomic weight ratio (AWP) = " << obj.atomic_weight_ratio << std::endl;
    out << "\tProduct modifier flag (LIP) = " << (int) obj.product_modifier_flag << std::endl;
    out << "\tDistribution law (LAW) = " << (int) obj.law << std::endl;
    out << "\tAngular represenation (LANG) = " << (int) obj.distr.angular_repr << std::endl;

    return out;
}

std::ostream& operator<<(std::ostream& out, const EnergyAngleData& obj) {
    out << "Charge and mass parameters (ZA) = " << (int) obj.ZA << std::endl;
    out << "Atomic weight ratio (AWR) = " << obj.atomic_weight_ratio << std::endl;
    out << "Information about promt fission neutrons and photons (JP) = " << (int) obj.yield_multiplicity_flag << std::endl;
    out << "Reference system for energy and angle (LCT) = " << (int) obj.reference_frame << std::endl;
    out << "Subsections (total number is " << obj.subsections.size() << "):" << std::endl;

    for (const auto& subsection: obj.subsections) {
        out << subsection;
        std::cout << std::string(60, '-') << std::endl;
    }
    return out;
}

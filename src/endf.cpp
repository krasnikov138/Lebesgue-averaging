#include "endf.h"

void EndfFile::open(const std::string& fname) {
    pfile.open(fname);
    if (pfile.is_open())
        is_open = true;
    else
        throw std::runtime_error("Can't open endf file.");
}

void EndfFile::reset() {
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

double EndfFile::get_number(char i) {
    if (eof | i > 5 | i < 0)
        throw std::runtime_error("Endf parsing error, line: " + std::to_string(line_number));
    std::string snum = line.substr(i * 11, 11);
    // check empty string
    if (std::all_of(snum.begin(), snum.end(), [](char ch) {return std::isspace(ch);}))
        throw std::runtime_error("Endf parsing error, line: " + std::to_string(line_number) + 
                                 ", cell number: " + std::to_string((int) i));

    double value;
    auto pos = snum.find_last_of("-+");
    if (pos == std::string::npos | pos == 0)
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
        first = static_cast<T>(get_number(pos));
    if constexpr(sizeof...(args) > 0)
        _get_numbers(++pos, args...);
}

template<typename... Ts>
void EndfFile::get_numbers(Ts&... args) {
    _get_numbers(0, args...);
}

void EndfFile::read(CrossSectionData& data) {
    size_t ranges;

    get_numbers(data.za, data.awr);
    next_line(false);
    get_numbers(data.qm, data.qi, skip, data.lr, ranges, data.np);
    // number of interpolation ranges (only 1 is supported)
    assert(ranges == 1);
    next_line(false);
    get_numbers(skip, data.int_type);

    size_t count = 0;
    data.energies.reserve(data.np);
    data.cross_sections.reserve(data.np);
    for (size_t i = 0; i < data.np; i++) {
        if (count % 6 == 0) next_line(false);
        data.energies.push_back(get_number(count++ % 6));
        data.cross_sections.push_back(get_number(count++ % 6));
    }
}

void EndfFile::read(EnergyAngleData& distr) {
    get_numbers(distr.za, distr.awr, distr.jp, distr.lct);

    auto num_of_sections = static_cast<size_t>(get_number(4));
    for (size_t i = 0; i < num_of_sections; i++) {
        ContinuumSubsection section;
        next_line(false);
        read(section);
        distr.subsections.push_back(std::move(section));
    }
}

void EndfFile::read(ContinuumSubsection& section) {
    size_t ranges;
    get_numbers(section.zap, section.awp, section.lip, section.law, ranges, section.np);
    // by far only continuum subsection is supported
    assert(section.law == 1);
    // number of interpolation ranges (only 1 is supported)
    assert(ranges == 1);

    next_line(false);
    get_numbers(skip, section.int_type_primary);
    for (unsigned int i = 0; i < section.np; i++) {
        if (2 * i % 6 == 0)
            next_line(false);
        section.primary_energies.push_back(get_number(2 * i % 6));
        section.yields.push_back(get_number((2 * i + 1) % 6));
    }
    next_line(false);
    get_numbers(skip, skip, section.lang, section.lep, ranges, section.ne);
    assert(ranges == 1);
    next_line(false);
    get_numbers(skip, section.int_type_primary);

    for (unsigned int i = 0; i < section.ne; i++) {
        SecondaryEnergyTab data;
        next_line(false);
        read(data);
        section.distribution.push_back(std::move(data));
    }
}

void EndfFile::read(SecondaryEnergyTab& data){
    get_numbers(skip, skip, data.nd, data.na, skip, data.nep);

    unsigned int pos = 0;
    data.secondary_energies.reserve(data.nep);
    data.coefs.reserve(data.nep);
    for (unsigned int i = 0; i < data.nep; i++) {
        if (pos % 6 == 0) next_line(false);
        data.secondary_energies.push_back(get_number(pos++ % 6));

        std::vector<double> legendre_coefs;
        legendre_coefs.reserve(data.na + 1);
        for (unsigned int j = 0; j < data.na + 1; j++) {
            if (pos % 6 == 0) next_line(false);
            legendre_coefs.push_back(get_number(pos++ % 6));
        }
        data.coefs.push_back(std::move(legendre_coefs));
    }
}

std::unique_ptr<EndfData> EndfFile::get_section(unsigned int mf, unsigned int mt) {
    if (!is_open)
        throw std::runtime_error("File isn't opened.");
    if (mf != 3 && mf != 6)
        throw std::runtime_error("Sections 3 and 6 are only supported.");
    do {
        next_line();
    } while (!eof && (cur_mf == 0 || cur_mf < mf));
    if (eof || cur_mf > mf) {
        std::cout << "Can't find given MF in file." << std::endl;
        return {};
    }

    if (mf == 3) {
        auto data = new CrossSectionData();
        read(*data);
        return std::unique_ptr<EndfData>(data);
    } else if (mf == 6) {
        auto data = new EnergyAngleData();
        read(*data);
        return std::unique_ptr<EndfData>(data);
    }
}

std::ostream& operator<<(std::ostream& out, const CrossSectionData& obj) {
    out << "Charge and mass parameters: ZA=" << (int) obj.za << ", AWR=" << obj.awr << std::endl;
    out << "Mass differences: QM=" << obj.qm << ", QI=" << obj.qi << std::endl;
    out << "Interpolation type: " << (int) obj.int_type << std::endl;

    return out;
}

std::ostream& operator<<(std::ostream& out, const ContinuumSubsection& obj) {
    out << "\tProduct of reaction info: ZAP=" << (int) obj.zap << 
           ", AWP=" << obj.awp << ", LIP=" << (int) obj.lip << std::endl;
    out << "\tWay of angular representation: LANG=" << (int) obj.lang << std::endl;
    out << "\tInterpolation parameters: " << (int) obj.int_type_primary << ' ' << 
           (int) obj.int_type_secondary << ' ' << (int) obj.lep << std::endl;
    out << "\tPoints number: NP=" << obj.np << ", NE=" << obj.ne << std::endl;

    return out;
}

std::ostream& operator<<(std::ostream& out, const EnergyAngleData& obj) {
    out << "Charge and mass parameters: ZA=" << (int) obj.za << ", AWR=" << obj.awr << std::endl;
    out << "Information about promt fission neutrons and photons: JP=" << (int) obj.jp << std::endl;
    out << "Reference system for energy and angle: LCT=" << (int) obj.lct << std::endl;
    out << "Subsections: " << std::endl;

    for (const auto& subsection: obj.subsections) {
        out << subsection;
        std::cout << std::string(65, '-') << std::endl;
    }
    return out;
}

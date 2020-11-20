#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <sstream>
#include "endf.h"

namespace py = pybind11;
using namespace py::literals;

template<typename Sequence>
py::array_t<typename Sequence::value_type> move_to_np(Sequence&& seq) {
	Sequence* seq_ptr = new Sequence(std::move(seq));
	auto capsule = py::capsule(seq_ptr, [](void *p) {delete reinterpret_cast<Sequence*>(p);});
    return py::array(seq_ptr->size(), seq_ptr->data(), capsule);
}

template<template<typename> class Seq, typename T>
py::array_t<T> copy2d_to_np(const Seq<Seq<T>>& obj) {
    if (!obj.size())
        return {};

    size_t N = obj.size();
    size_t M = obj[0].size();

    auto copy = py::array_t<T>(N * M);
    T *ptr = static_cast<T*>(copy.request().ptr);
    for (int i = 0; i < N; ++i)
        std::copy(obj[i].cbegin(), obj[i].cend(), ptr + i * M);
    copy.resize({(int) N, (int) M});

    return copy;
}

std::string int_schema(unsigned char int_type) {
	const static std::vector<std::string> prefixes = {
		"", "corresponding_points_", "unit_base_"
	};
	const static std::vector<std::string> schemas = {
		"histogram", "linlin", "linlog", "loglin", "loglog", "chargedparticles"
	};
	return prefixes[int_type / 10] + schemas[int_type % 10 - 1];
}

struct PyInterpolationTable {
    py::array_t<size_t> boundaries;
    py::array_t<unsigned char> interpolation_types;

    py::array_t<double> xs, ys;

    PyInterpolationTable(InterpolationTable&& obj) {
        boundaries = move_to_np(std::move(obj.boundaries));
        interpolation_types = move_to_np(std::move(obj.interpolation_types));
        xs = move_to_np(std::move(obj.xs));
        ys = move_to_np(std::move(obj.ys));
    }
};

py::dict cross_section_to_py(CrossSectionData&& obj) {
    py::dict data;

    data["za"] = obj.za;
    data["awr"] = obj.awr;
    data["lr"] = obj.lr;
    data["qm"] = obj.qm;
    data["qi"] = obj.qi;
    data["cs"] = PyInterpolationTable(std::move(obj.cs));

    return data;
}

py::dict continuum_distribution_to_py(ContinuumDistribution&& obj) {
    py::dict data;

    data["angular_repr"] = obj.angular_repr;
    data["lep"] = obj.lep;

    data["primary_energies"] = move_to_np(std::move(obj.primary_energies));
    data["discrete_energies"] = move_to_np(std::move(obj.discrete_energies));
    data["boundaries"] = move_to_np(std::move(obj.boundaries));
    data["interpolation_types"] = move_to_np(std::move(obj.interpolation_types));

    py::list secondary_energies, coefs;
    for (auto& el: obj.secondary_energies)
        secondary_energies.append(move_to_np(std::move(el)));

    for (auto& el: obj.coefs)
        coefs.append(copy2d_to_np(el));

    data["secondary_energies"] = secondary_energies;
    data["coefs"] = coefs;

    return data;
}

py::dict product_subsection_to_py(ProductSubsection&& obj) {
    py::dict data;

    data["ZAP"] = obj.ZAP;
    data["atomic_weight_ratio"] = obj.atomic_weight_ratio;
    data["product_modifier_flag"] = obj.product_modifier_flag;
    data["law"] = obj.law;

    data["energy_yield"] = PyInterpolationTable(std::move(obj.energy_yield));
    data["distr"] = continuum_distribution_to_py(std::move(obj.distr));

    return data;
}

py::dict energy_angle_to_py(EnergyAngleData&& obj) {
    py::dict data;

    data["ZA"] = obj.ZA;
    data["atomic_weight_ratio"] = obj.atomic_weight_ratio;
    data["yield_multiplicity_flag"] = obj.yield_multiplicity_flag;
    data["reference_frame"] = obj.reference_frame;

    py::list subsections;
    for (auto& subsection: obj.subsections)
        subsections.append(product_subsection_to_py(std::move(subsection)));
    data["subsections"] = subsections;

    return data;
}

class PyEndfFile {
    EndfFile endf;

    py::array_t<uint16_t> content;
public:
    PyEndfFile() noexcept {};
    PyEndfFile(const std::string& fname): endf(fname) {}

    void open(const std::string& fname) {
    	endf.open(fname);
    }

    void close() {
    	endf.close();
    }

    py::array_t<uint16_t>& table_of_content() {
    	if (content.size())
    		return content;

    	const auto& endf_content = endf.get_table_of_content();
    	content = py::array_t<uint16_t>(2 * endf_content.size());
    	uint16_t *ptr = static_cast<uint16_t*>(content.request().ptr);
    	for (int i = 0; i < endf_content.size(); ++i) {
    		ptr[2 * i] = endf_content[i].first;
    		ptr[2 * i + 1] = endf_content[i].second;
    	}
    	content.resize({(int) endf_content.size(), 2});

    	return content;
    }

    py::dict get_cross_section(unsigned int mt) {
    	auto data = std::get<CrossSectionData>(endf.get_section(3, mt));
    	return cross_section_to_py(std::move(data));
    }

    py::dict get_energy_angle() {
        auto data = std::get<EnergyAngleData>(endf.get_section(6, 5));
        return energy_angle_to_py(std::move(data));
    }
};

PYBIND11_MODULE(bindings, m) {
    py::class_<PyInterpolationTable>(m, "PyInterpolationTable")
        .def_readonly("boundaries", &PyInterpolationTable::boundaries)
        .def_readonly("interpolation_types", &PyInterpolationTable::interpolation_types)
        .def_readonly("xs", &PyInterpolationTable::xs)
        .def_readonly("ys", &PyInterpolationTable::ys);

    py::class_<PyEndfFile>(m, "PyEndfFile")
    	.def(py::init<>())
        .def(py::init<const std::string&>())
        .def("open", &PyEndfFile::open)
        .def("close", &PyEndfFile::close)
        .def_property("table_of_content", &PyEndfFile::table_of_content, nullptr)
        .def("get_cross_section", &PyEndfFile::get_cross_section)
        .def("get_energy_angle", &PyEndfFile::get_energy_angle);
}


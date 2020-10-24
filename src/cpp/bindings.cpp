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

std::string int_schema(unsigned char int_type) {
	const static std::vector<std::string> prefixes = {
		"", "corresponding_points_", "unit_base_"
	};
	const static std::vector<std::string> schemas = {
		"histogram", "linlin", "linlog", "loglin", "loglog", "chargedparticles"
	};
	return prefixes[int_type / 10] + schemas[int_type % 10 - 1];
}

struct PyCrossSectionData {
    float za, awr;
    double qm, qi;
    unsigned char lr;

    unsigned int np;
    unsigned char int_type;

    py::array_t<double> energies;
    py::array_t<double> cross_sections;

	PyCrossSectionData(CrossSectionData&& cdata) {
		za = cdata.za;
		awr = cdata.awr;
		lr = cdata.lr;
		qm = cdata.qm;
		qi = cdata.qi;
		int_type = cdata.int_type;

		energies = move_to_np(std::move(cdata.energies));
		cross_sections = move_to_np(std::move(cdata.cross_sections));
	}

	std::string interpolation_type() const {
		return int_schema(int_type);
	}

	std::string to_string() const {
		std::stringstream buf;
		buf << "ENDF cross section data object\n";
		buf << "Material charge and mass parameters ZA=" << za << ", mass=" << awr << '\n';
    	buf << "Mass differences QM=" << qm << ", QI=" << qi << '\n';
    	buf << "Breakup flag LR=" << (int)lr << '\n';
    	buf << "Interpolation type: " << interpolation_type();
    	return buf.str();
	}
};

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

    PyCrossSectionData get_cross_section(unsigned int mt) {
    	std::unique_ptr<EndfData> data = endf.get_section(3, mt);
    	auto cs_data_ptr = reinterpret_cast<CrossSectionData*>(data.release());
    	PyCrossSectionData result(std::move(*cs_data_ptr));
    	return result;
    }
    // std::unique_ptr<EndfData> get_section(unsigned int mf, unsigned int mt);
};

PYBIND11_MODULE(bindings, m) {
	py::class_<PyCrossSectionData>(m, "PyCrossSectionData")
		.def_readonly("za", &PyCrossSectionData::za)
		.def_readonly("awr", &PyCrossSectionData::awr)
		.def_readonly("lr", &PyCrossSectionData::lr)
		.def_readonly("qm", &PyCrossSectionData::qm)
		.def_readonly("qi", &PyCrossSectionData::qi)
		.def_readonly("energies", &PyCrossSectionData::energies)
		.def_readonly("cross_sections", &PyCrossSectionData::cross_sections)
		.def_property("interpolation_type", &PyCrossSectionData::interpolation_type, nullptr)
		.def("__repr__", &PyCrossSectionData::to_string);

    py::class_<PyEndfFile>(m, "PyEndfFile")
    	.def(py::init<>())
        .def(py::init<const std::string&>())
        .def("open", &PyEndfFile::open)
        .def("close", &PyEndfFile::close)
        .def_property("table_of_content", &PyEndfFile::table_of_content, nullptr)
        .def("get_cross_section", &PyEndfFile::get_cross_section);
}


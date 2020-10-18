#include <iostream>
#include <map>
#include "yaml-cpp/yaml.h"

struct MaterialInfo {
    double concentration;
    std::string cross_section_file;

    MaterialInfo(double conc = 1.0, const std::string& cross_section_file = ""):
    	concentration(conc), cross_section_file(cross_section_file) {}
    
    MaterialInfo(const YAML::Node& node) {
    	concentration = node["concentration"].as<double>();
    	cross_section_file = node["cross_section_file"].as<std::string>();
    }
};

std::map<std::string, MaterialInfo> load_yaml_config(const std::string& fname) {
	YAML::Node doc = YAML::LoadFile(fname);

    std::string data_path;
    if (auto dp_node = doc["data_path"])
	    data_path = dp_node.as<std::string>();

	if (auto material_node = doc["materials"])
	    doc = material_node;
	else
		throw std::runtime_error("File doesn't contain 'materials' section.");

    std::map<std::string, MaterialInfo> res;
    for (auto it = doc.begin(); it != doc.end(); ++it) {
	    std::string name = it->first.as<std::string>();
	    MaterialInfo info(it->second);

	    info.cross_section_file = data_path + '/' + info.cross_section_file;
	    res[name] = std::move(info);
	}

    return res;
}

int main() {
	const std::string test_file_name = "materials.yaml";
	std::map<std::string, MaterialInfo> config = load_yaml_config(test_file_name);

	for (const auto& [name, info]: config) {
		std::cout << name << " " << info.concentration << " " << info.cross_section_file << '\n';
	}
	return 0;
}

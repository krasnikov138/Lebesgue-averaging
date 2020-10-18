#include <iostream>
#include <map>
#include "yaml-cpp/yaml.h"

#ifndef YAML_CONFIG_PARSER
#define YAML_CONFIG_PARSER

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

std::map<std::string, MaterialInfo> load_yaml_config(const std::string& fname);

#endif

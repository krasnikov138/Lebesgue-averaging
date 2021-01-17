#include <iostream>
#include <map>
#include "yaml-cpp/yaml.h"

#ifndef YAML_CONFIG_PARSER
#define YAML_CONFIG_PARSER

struct MaterialInfo {
    double average_concentration;
    int energy_angle_zap;
    int energy_angle_index;

    std::string total_section_file;
   	std::string energy_angle_file;
   	std::string cross_section_file;

    std::vector<int> reactions_mt;

    MaterialInfo(
    	double average_concentration = 1.0,
    	int energy_angle_zap = -1, 
    	int energy_angle_index = -1,
    	const std::string& total_section_file = {}, 
    	const std::string& energy_angle_file = {}, 
    	const std::string& cross_section_file = {}, 
    	const std::vector<int>& reactions_mt = {}
    ):
		average_concentration(average_concentration),
		energy_angle_zap(energy_angle_zap),
		energy_angle_index(energy_angle_index),
		total_section_file(total_section_file),
		energy_angle_file(energy_angle_file),
		cross_section_file(cross_section_file),
		reactions_mt(reactions_mt) {}

	MaterialInfo(const YAML::Node& node);
};

std::map<std::string, MaterialInfo> load_yaml_config(const std::string& fname);

#endif

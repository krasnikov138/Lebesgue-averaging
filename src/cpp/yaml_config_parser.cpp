#include "yaml_config_parser.h"

MaterialInfo::MaterialInfo(const YAML::Node& node):
	average_concentration(1.0),
    energy_angle_zap(-1), 
    energy_angle_index(-1) {
	average_concentration = node["average_concentration"].as<double>();
	total_section_file = node["total_cross_section"]["fname"].as<std::string>();
	const auto& energy_angle = node["energy_angle"];
	if (auto index_node = energy_angle["index"])
		energy_angle_index = index_node.as<int>();
	if (auto zap_node = energy_angle["ZAP"])
		energy_angle_zap = zap_node.as<int>();
	energy_angle_file = energy_angle["fname"].as<std::string>();
	cross_section_file = node["reactions"]["fname"].as<std::string>();
	reactions_mt = node["reactions"]["mt"].as<std::vector<int>>();
}

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

	    if (data_path.length()) {
	    	info.total_section_file = data_path + '/' + info.total_section_file;
	    	info.energy_angle_file = data_path + '/' + info.energy_angle_file;
	    	info.cross_section_file = data_path + '/' + info.cross_section_file;
	    }
	    res[name] = std::move(info);
	}

    return res;
}

#include "yaml_config_parser.h"

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

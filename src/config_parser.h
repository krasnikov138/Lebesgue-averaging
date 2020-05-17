#include <iostream>
#include <fstream>
#include <cctype>
#include <string>
#include <string_view>
#include <exception>
#include <map>

std::string_view strip(const std::string_view& str);
std::string_view strip(const std::string& str);
bool is_empty(const std::string& str);
bool startswith(const std::string& str, const std::string& prefix, unsigned int N = 1);

class ConfigParser {
public:
	using Block = std::map<std::string, std::string>;

	class ParseError: public std::exception {
		std::string message;
	public:
		ParseError(const std::string& error): message(error) {};
		ParseError(size_t line_number, const std::string& error) {
			message = "Error on line " + std::to_string(line_number) + ": " + error;
		};
		const char* what() const noexcept {
			return message.data();
		}
	};
private:
	std::string indent;
	// parsing state descriptors
	std::ifstream conf_file;
	bool eof_flag;
	std::string line;
	size_t line_number;

	bool next_line();
	std::string get_block_name(unsigned int level);
	Block read_block(unsigned int level);
public:
	ConfigParser(const std::string& indent = "\t"): indent(indent) {};

	std::map<std::string, Block> load(const std::string& fname);
};

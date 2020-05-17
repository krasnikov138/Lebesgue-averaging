#include "config_parser.h"

std::string_view strip(const std::string_view& str) {
    const char spaces[] = "\n\t\f\r\v ";
    std::string_view res(str);
    auto pos = std::min(str.find_first_not_of(spaces), str.size());
    res.remove_prefix(pos);
    pos = str.find_last_not_of(spaces);
    if (pos != std::string::npos)
        res.remove_suffix(str.length() - pos - 1);

    return res;
}

std::string_view strip(const std::string& str) {
    return strip(std::string_view(str));
}

bool is_empty(const std::string& str) {
    return strip(str).length() == 0 ? true : false;
}

bool startswith(const std::string& str, const std::string& prefix, unsigned int N) {
    for (unsigned int i = 0; i < N; i++)
        if (str.compare(prefix.length() * i, prefix.length(), prefix))
            return false;
    return true;
}

bool ConfigParser::next_line() {
    while (std::getline(conf_file, line)) {
        ++line_number;
        if (!is_empty(line) && strip(line).front() != '#')
            return true;
    }
    eof_flag = true;
    return false;
}

std::string ConfigParser::get_block_name(unsigned int level) {
    std::string name;

    if (!startswith(line, indent, level) || std::isspace(line[indent.length() * level]))
        throw ParseError(line_number, "wrong indentation");
    name = strip(line);
    if (name.back() != ':')
        throw ParseError(line_number, "no ':' in the end of block name");
    name = name.substr(0, name.length() - 1);

    return name;
}

ConfigParser::Block ConfigParser::read_block(unsigned int level) {
    Block block;

    while (next_line()) {
        if (!startswith(line, indent, level)) 
            break;
        if (std::isspace(line[indent.length() * level]))
            throw ParseError(line_number, "wrong indentation");

        auto pos = line.find(':');
        if (pos != std::string::npos) {
            std::string key(strip(line.substr(0, pos)));
            std::string value(strip(line.substr(pos + 1)));

            if (!key.length())
                throw ParseError(line_number, "empty key");
            if (!value.length())
                throw ParseError(line_number, "empty value");

            auto [iter, success] = block.emplace(key, value);
            if (!success)
                throw ParseError(line_number, "key already exists");
        } else
            throw ParseError(line_number, "empty value");
    }
    return block;
}

std::map<std::string, ConfigParser::Block> ConfigParser::load(const std::string& fname) {
    eof_flag = false;
    line_number = 0;

    std::map<std::string, Block> values;

    conf_file.open(fname);
    if (!conf_file.is_open())
        throw std::runtime_error("Can't open config file: '" + fname + "'");

    next_line();
    while (!eof_flag) {
        std::string name = get_block_name(0);
        Block block = read_block(1);

        if (name.length()) {
            if (block.size() == 0)
                throw ParseError("Block '" + name + "' has empty body");
            auto [iter, success] = values.emplace(name, std::move(block));
            if (!success)
                throw ParseError("Block '" + name + "' is duplicated");
        }
    }
    return values;
}

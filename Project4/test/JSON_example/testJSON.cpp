// test/JSON_example/testJSON.cpp
// run with make run JSON=test/JSON_example/testJSON.json

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "nlohmann/json.hpp"

using json = nlohmann::json;

struct Config {
    struct Positions { double x, y, z; } positions;
    struct Output { std::string directory, file_name; } output;
};

// json -> Config::Positions
void from_json(const json& j, Config::Positions& p) {
    // j is a JSON object representing positions
    // p is the Config::Positions struct
    j.at("x").get_to(p.x); 
    j.at("y").get_to(p.y);
    j.at("z").get_to(p.z);
}
// json -> Config::Output
void from_json(const json& j, Config::Output& o) {
    j.at("directory").get_to(o.directory);
    j.at("file_name").get_to(o.file_name);
}

// json -> Config
void from_json(const json& j, Config& c) {
    j.at("positions").get_to(c.positions);
    j.at("output").get_to(c.output);
}

int main(int argc, char** argv) { // argc and argv to get JSON file path (argc is number of arguments, argv is array of arguments)
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config.json>\n";
        return 1;
    }
    const std::string cfg_path = argv[1];

    std::ifstream f(cfg_path);
    if (!f) {
        std::cerr << "Could not open JSON file: " << cfg_path << "\n";
        return 1;
    }

    json j; f >> j; // parse JSON file into json object
    Config cfg = j.get<Config>(); // convert json object to Config struct

    std::filesystem::create_directories(cfg.output.directory);
    const std::string outpath = cfg.output.directory + cfg.output.file_name; 

    std::ofstream out(outpath);
    if (!out) { // check if file opened successfully
        std::cerr << "Could not open output file: " << outpath << "\n";
        return 2;
    }

    out << "# positions (x y z)\n" // write header
        << cfg.positions.x << " "
        << cfg.positions.y << " "
        << cfg.positions.z << "\n";
    out.close();

    std::cout << "Wrote positions to: " << outpath << "\n";
    return 0;
}

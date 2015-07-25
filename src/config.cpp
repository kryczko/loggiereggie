#include "yaml-cpp/yaml.h"
#include "storage.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include <string>

void parse_yaml_file(Storage& storage) {
    YAML::Node config = YAML::LoadFile("input.yaml");
    storage.set_filename(config["filename"].as<std::string>());
    storage.set_alpha(config["alpha"].as<double>());
    
    
    std::cout << "Read configuration file...\n";
}
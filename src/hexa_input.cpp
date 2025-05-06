#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "hexa.h"
#include "pml.h"
void parseMaterial(const std::string& line, std::vector<Material>& materials) {
    std::istringstream iss(line);
    Material mat;
    iss >> mat.type >> mat.vp >> mat.vs >> mat.rho;
    materials.push_back(mat);
}

Input readInputFile(const std::string& filePath) {
    Input input;
    std::ifstream file(filePath);
    std::string line;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        exit(EXIT_FAILURE);
    }

    while (std::getline(file, line)) {
        // Remove comments and trim whitespace
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        if (line.empty()) continue;

        // Parse key-value pairs
        if (line.find("topo") == 0) {
            input.topo = line.substr(line.find('=') + 1);
            input.topo.erase(0, input.topo.find_first_not_of(" \t"));
        } else if (line.find("interfaceNumber") == 0) {
            input.interfaceNumber = std::stoi(line.substr(line.find('=') + 1));
        } else if (line.find("inter") == 0) {
            input.inter = line.substr(line.find('=') + 1);
            input.inter.erase(0, input.inter.find_first_not_of(" \t"));
        } else if (line.find("ref") == 0) {
            input.ref = std::stoi(line.substr(line.find('=') + 1));
        } else if (line.find("nmat") == 0) {
            input.nmat = std::stoi(line.substr(line.find('=') + 1));
        } else if (line.find("S") == 0 || line.find("F") == 0) {
            parseMaterial(line, input.materials);
        } else if (line.find("PML") == 0) {
            input.PML = std::stoi(line.substr(line.find('=') + 1)) == 1;
        } else if (line.find("pmlx") == 0) {
            input.pmlx = std::stod(line.substr(line.find('=') + 1));
        } else if (line.find("nlayersx") == 0) {
            input.nlayersx = std::stoi(line.substr(line.find('=') + 1));
        } else if (line.find("pmly") == 0) {
            input.pmly = std::stod(line.substr(line.find('=') + 1));
        } else if (line.find("nlayersy") == 0) {
            input.nlayersy = std::stoi(line.substr(line.find('=') + 1));
        } else if (line.find("pmlz") == 0) {
            input.pmlz = std::stod(line.substr(line.find('=') + 1));
        } else if (line.find("nlayersz") == 0) {
            input.nlayersz = std::stoi(line.substr(line.find('=') + 1));
        } else if (line.find("A") == 0) {
            input.A = std::stoi(line.substr(line.find('=') + 1));
        } else if (line.find("npow") == 0) {
            input.npow = std::stoi(line.substr(line.find('=') + 1));
        } else if (line.find("meshOpt") == 0) {
            input.meshOpt = std::stoi(line.substr(line.find('=') + 1)) == 1;
        }
    }

    file.close();
    return input;
}

int inpreader(hexa_tree_t *mesh) {
    std::string filePath = "./HexMesh.input";
    Input input = readInputFile(filePath);
    mesh->input = input;
    // Output the parsed data for verification
    std::cout << "Topo: " << input.topo << std::endl;
    std::cout << "Interface Number: " << input.interfaceNumber << std::endl;
    std::cout << "Inter: " << input.inter << std::endl;
    std::cout << "Refinement Level: " << input.ref << std::endl;
    std::cout << "Number of Materials: " << input.nmat << std::endl;
    for (const auto& mat : input.materials) {
        std::cout << "Material: " << mat.type << " " << mat.vp << " " << mat.vs << " " << mat.rho << std::endl;
    }
    std::cout << "PML: " << (input.PML ? "Enabled" : "Disabled") << std::endl;
    std::cout << "PML X: " << input.pmlx << ", Layers X: " << input.nlayersx << std::endl;
    std::cout << "PML Y: " << input.pmly << ", Layers Y: " << input.nlayersy << std::endl;
    std::cout << "PML Z: " << input.pmlz << ", Layers Z: " << input.nlayersz << std::endl;
    std::cout << "Mesh Optimization: " << (input.meshOpt ? "Enabled" : "Disabled") << std::endl;

    return 0;
}
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "hexa.h"
#include "pml.h"
#include <regex>

void parseMaterial(const std::string &line, std::vector<Material> &materials)
{
    std::istringstream iss(line);
    Material mat;
    iss >> mat.type >> mat.vp >> mat.vs >> mat.rho;
    materials.push_back(mat);
}

static inline void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
        [](unsigned char ch){ return !std::isspace(ch); }));
}
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
        [](unsigned char ch){ return !std::isspace(ch); }).base(), s.end());
}
static inline void trim(std::string &s) { ltrim(s); rtrim(s); }

void parse_zcuts(const std::string &line, std::vector<double> &out) {
    out.clear();
    std::string s = line;

    size_t cpos = s.find('#');
    if (cpos != std::string::npos) s.erase(cpos);
    cpos = s.find("//");
    if (cpos != std::string::npos) s.erase(cpos);

    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ';')) {
        trim(token);
        if (token.empty()) continue;
        try {
            size_t idx = 0;
            double v = std::stod(token, &idx);
            std::string rest = token.substr(idx);
            trim(rest);
            if (!rest.empty()) {
                std::cerr << "Warning: trailing characters after number in zcuts token '" << token << "'\n";
                continue;
            }
            out.push_back(v);
        } catch (const std::exception &e) {
            std::cerr << "Warning: failed to parse zcuts token '" << token << "': " << e.what() << "\n";
        }
    }
}

Input readInputFile(const std::string &filePath)
{
    Input input;
    std::ifstream file(filePath);
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filePath << std::endl;
        exit(EXIT_FAILURE);
    }

    while (std::getline(file, line))
    {
        // Remove comments and trim whitespace
        size_t commentPos = line.find('#');
        if (commentPos != std::string::npos)
        {
            line = line.substr(0, commentPos);
        }
        // trim() uses std::isspace, so it also strips '\r' / '\n' — required
        // for HexMesh.input files saved with Windows (CRLF) line endings, where
        // a trailing '\r' would otherwise become part of the GTS file path and
        // make SurfaceRead fail (NULL surface → crash in GetMeshFromSurface).
        trim(line);

        if (line.empty())
            continue;

        if (line.find("topo") == 0)
        {
            input.topo = line.substr(line.find('=') + 1);
            trim(input.topo);
        }
        else if (line.find("interfaceNumber") == 0)
        {
            input.interfaceNumber = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("inter") == 0)
        {
            input.inter = line.substr(line.find('=') + 1);
            trim(input.inter);
        }
        else if (line.find("ref") == 0)
        {
            input.ref = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("movingNodes") == 0)
        {
            input.movingNodes = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("CgalUse") == 0)
        {
            input.CgalUse = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("nmat") == 0)
        {
            input.nmat = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("S") == 0 || line.find("F") == 0)
        {
            parseMaterial(line, input.materials);
        }
        else if (line.find("PML") == 0)
        {
            input.PML = std::stoi(line.substr(line.find('=') + 1)) == 1;
        }
        else if (line.find("pmlx") == 0)
        {
            input.pmlx = std::stod(line.substr(line.find('=') + 1));
        }
        else if (line.find("nlayersx") == 0)
        {
            input.nlayersx = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("pmly") == 0)
        {
            input.pmly = std::stod(line.substr(line.find('=') + 1));
        }
        else if (line.find("nlayersy") == 0)
        {
            input.nlayersy = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("pmlz") == 0)
        {
            input.pmlz = std::stod(line.substr(line.find('=') + 1));
        }
        else if (line.find("nlayersz") == 0)
        {
            input.nlayersz = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("A") == 0)
        {
            input.A = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("npow") == 0)
        {
            input.npow = std::stoi(line.substr(line.find('=') + 1));
        }
        else if (line.find("meshOpt") == 0)
        {
            input.meshOpt = std::stoi(line.substr(line.find('=') + 1)) == 1;
        }
        else if (line.find("zcuts") == 0)
        {
            parse_zcuts(line.substr(line.find('=') + 1), input.zcut);
        }
        else if (line.find("z") == 0) {
            input.z = std::stoi(line.substr(line.find('=') + 1));
        }
    }
    file.close();    
    return input;
}

int inpreader(hexa_tree_t *mesh)
{
    std::string filePath = "./HexMesh.input";
    Input input = readInputFile(filePath);
    mesh->input = input;

    std::cout << "Topo: " << input.topo << std::endl;
    std::cout << "Interface Number: " << input.interfaceNumber << std::endl;
    std::cout << "Inter: " << input.inter << std::endl;
    std::cout << "Refinement Level: " << input.ref << std::endl;
    std::cout << "Z-Depth: " << input.z << std::endl;
    std::cout << "Z-Cuts: ";
    for (const auto &zcut : input.zcut)
    {
        std::cout << zcut << " ";
    }
    std::cout << std::endl;
    std::cout << "Number of Materials: " << input.nmat << std::endl;
    for (const auto &mat : input.materials)
    {
        std::cout << "Material: " << mat.type << " " << mat.vp << " " << mat.vs << " " << mat.rho << std::endl;
    }
    
    std::cout << "Moving Nodes: " << input.movingNodes << std::endl;

    std::cout << "PML: " << (input.PML ? "Enabled" : "Disabled") << std::endl;
    std::cout << "PML X: " << input.pmlx << ", Layers X: " << input.nlayersx << std::endl;
    std::cout << "PML Y: " << input.pmly << ", Layers Y: " << input.nlayersy << std::endl;
    std::cout << "PML Z: " << input.pmlz << ", Layers Z: " << input.nlayersz << std::endl;
    std::cout << "Mesh Optimization: " << (input.meshOpt ? "Enabled" : "Disabled") << std::endl;

    // Verify material count matches declaration
    if (input.nmat != input.materials.size()) {
        std::cerr << "Warning: Declared number of materials (nmat=" << input.nmat 
                  << ") differs from actual number of materials defined (" 
                  << input.materials.size() << ")\n";
    }

    // use CGAL 
    std::cout << "Using CGAL exact kernel: " << input.CgalUse << std::endl;

    return 0;
}

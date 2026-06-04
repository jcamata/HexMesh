#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <algorithm>
#include <filesystem>
#include <iomanip>

#include <netcdf.h>
#include <curl/curl.h>

// --- Configuration Constants ---
const double EARTH_RADIUS_M = 6371000.0;
const double PI = 3.14159265358979323846;

// --- Data Structures ---

struct Point3D {
    double x, y, z; // Using double for metric precision
};

struct Triangle {
    int v1, v2, v3;
};

// --- CURL Callback ---

size_t write_data(void* ptr, size_t size, size_t nmemb, FILE* stream) {
    return fwrite(ptr, size, nmemb, stream);
}

// --- 1. Download Logic with Existence Check ---

bool download_file_if_missing(const std::string& url, const std::string& local_path) {
    if (std::filesystem::exists(local_path)) {
        std::cout << "File '" << local_path << "' already exists. Skipping download." << 
std::endl;
        return true;
    }

    std::cout << "File not found. Downloading from " << url << "..." << std::endl;
    
    CURL* curl = curl_easy_init();
    if (!curl) return false;

    FILE* fp = fopen(local_path.c_str(), "wb");
    if (!fp) {
        std::cerr << "Error opening file for writing: " << local_path << std::endl;
        return false;
    }

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_data);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, fp);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L); // Follow redirects
    curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 0L);     // Show progress

    CURLcode res = curl_easy_perform(curl);
    
    fclose(fp);
    curl_easy_cleanup(curl);

    if (res != CURLE_OK) {
        std::cerr << "Download failed: " << curl_easy_strerror(res) << std::endl;
        // Remove partially downloaded file if failed
        std::filesystem::remove(local_path);
        return false;
    }
    
    std::cout << "Download complete." << std::endl;
    return true;
}

// --- 2. NetCDF Extraction ---

bool extract_region(const std::string& nc_path, 
                    double lat_min, double lat_max, 
                    double lon_min, double lon_max,
                    std::vector<double>& elevations,
                    int& width, int& height,
                    double& pixel_res_meters) 
{
    int ncid, varid;
    
    if (nc_open(nc_path.c_str(), NC_NOWRITE, &ncid) != NC_NOERR) {
        std::cerr << "Error opening NetCDF file." << std::endl;
        return false;
    }

    int dimid_lat, dimid_lon;
    size_t dim_len_lat, dim_len_lon;

    // GEBCO standard names
    if (nc_inq_dimid(ncid, "lat", &dimid_lat) != NC_NOERR) nc_inq_dimid(ncid, "y", 
&dimid_lat);
    if (nc_inq_dimid(ncid, "lon", &dimid_lon) != NC_NOERR) nc_inq_dimid(ncid, "x", 
&dimid_lon);
    
    nc_inq_dimlen(ncid, dimid_lat, &dim_len_lat);
    nc_inq_dimlen(ncid, dimid_lon, &dim_len_lon);

    // GEBCO is typically 15 arc-seconds (~0.004166 deg)
    double lat_res_deg = 180.0 / dim_len_lat;
    double lon_res_deg = 360.0 / dim_len_lon;

    // Calculate indices (Assuming GEBCO global grid)
    // GEBCO usually starts at North (90) going South (-90)
    size_t start_lat = static_cast<size_t>((90.0 - lat_max) / lat_res_deg);
    size_t start_lon = static_cast<size_t>((lon_min + 180.0) / lon_res_deg);

    size_t count_lat = static_cast<size_t>((lat_max - lat_min) / lat_res_deg) + 1;
    size_t count_lon = static_cast<size_t>((lon_max - lon_min) / lon_res_deg) + 1;

    width = count_lon;
    height = count_lat;

    // Calculate approximate pixel resolution in meters at the center of the region
    double center_lat = (lat_min + lat_max) / 2.0;
    pixel_res_meters = lat_res_deg * (PI / 180.0) * EARTH_RADIUS_M; // meters per degree 
lat
    // Note: Longitudinal resolution is narrower, but for simple mesh gen we use a square 
approximation or average
    double lon_res_meters = lon_res_deg * (PI / 180.0) * EARTH_RADIUS_M * 
std::cos(center_lat * PI / 180.0);
    
    std::cout << "Estimated Pixel Resolution: " << pixel_res_meters << " x " << 
lon_res_meters << " meters." << std::endl;

    elevations.resize(width * height);

    // Get Variable ID
    if (nc_inq_varid(ncid, "elevation", &varid) != NC_NOERR) {
        if (nc_inq_varid(ncid, "z", &varid) != NC_NOERR) {
             std::cerr << "Cannot find elevation variable." << std::endl;
             nc_close(ncid);
             return false;
        }
    }

    size_t start[2] = {start_lat, start_lon};
    size_t count[2] = {count_lat, count_lon};
    ptrdiff_t stride[2] = {1, 1};

    short* raw_data = new short[width * height];
    int status = nc_get_vars(ncid, varid, start, count, stride, raw_data);
    
    if (status != NC_NOERR) {
        std::cerr << "Error reading data: " << nc_strerror(status) << std::endl;
        delete[] raw_data;
        nc_close(ncid);
        return false;
    }

    for (int i = 0; i < width * height; ++i) {
        elevations[i] = static_cast<double>(raw_data[i]);
    }

    delete[] raw_data;
    nc_close(ncid);
    return true;
}

// --- 3. Filter Logic ---

// Simple Moving Average filter to remove features smaller than 'h' meters.
void apply_size_filter(std::vector<double>& data, int width, int height, double 
min_feature_size_m, double pixel_res_m) {
    if (min_feature_size_m <= 0) return;

    // Calculate kernel radius in pixels
    int kernel_radius = static_cast<int>(std::ceil(min_feature_size_m / pixel_res_m / 
2.0));
    
    if (kernel_radius < 1) return; // Feature size smaller than pixel size

    std::cout << "Applying filter to remove features < " << min_feature_size_m << "m 
(Kernel radius: " << kernel_radius << " px)" << std::endl;

    std::vector<double> filtered_data = data; // Copy

    for (int r = 0; r < height; ++r) {
        for (int c = 0; c < width; ++c) {
            double sum = 0;
            int count = 0;

            // Iterate over kernel window
            for (int kr = -kernel_radius; kr <= kernel_radius; ++kr) {
                for (int kc = -kernel_radius; kc <= kernel_radius; ++kc) {
                    int nr = r + kr;
                    int nc = c + kc;

                    // Check bounds
                    if (nr >= 0 && nr < height && nc >= 0 && nc < width) {
                        sum += data[nr * width + nc];
                        count++;
                    }
                }
            }
            filtered_data[r * width + c] = sum / count;
        }
    }
    data = std::move(filtered_data);
}

// --- 4. Projection Logic (Degrees to Meters) ---

// Projects lat/lon to a local metric grid centered at center_lat, center_lon
Point3D project_to_meters(double lat, double lon, double elev, double center_lat, double 
center_lon) {
    // X: Easting
    // Approximation: 1 degree lon = cos(lat) * 111320 meters
    double x = (lon - center_lon) * (PI / 180.0) * EARTH_RADIUS_M * std::cos(center_lat * 
PI / 180.0);
    
    // Y: Northing
    // Approximation: 1 degree lat = 111320 meters (approx)
    double y = (lat - center_lat) * (PI / 180.0) * EARTH_RADIUS_M;
    
    return {x, y, elev};
}

// --- Mesh Generation ---

void generate_mesh(const std::vector<double>& elevations, 
                   int width, int height,
                   double lat_min, double lat_max, 
                   double lon_min, double lon_max,
                   std::vector<Point3D>& vertices,
                   std::vector<Triangle>& triangles) 
{
    double d_lon = (lon_max - lon_min) / (width - 1);
    double d_lat = (lat_max - lat_min) / (height - 1);
    double center_lat = (lat_min + lat_max) / 2.0;
    double center_lon = (lon_min + lon_max) / 2.0;

    vertices.resize(width * height);
    for (int r = 0; r < height; ++r) {
        for (int c = 0; c < width; ++c) {
            int idx = r * width + c;
            
            double lat = lat_max - r * d_lat;
            double lon = lon_min + c * d_lon;
            
            // Convert to meters
            vertices[idx] = project_to_meters(lat, lon, elevations[idx], center_lat, 
center_lon);
        }
    }

    triangles.reserve(2 * (width - 1) * (height - 1));
    for (int r = 0; r < height - 1; ++r) {
        for (int c = 0; c < width - 1; ++c) {
            int tl = r * width + c;
            int tr = tl + 1;
            int bl = (r + 1) * width + c;
            int br = bl + 1;

            // Triangle 1
            triangles.push_back({tl, bl, tr});
            // Triangle 2
            triangles.push_back({tr, bl, br});
        }
    }
}

// --- Writers ---

void write_stl(const std::string& filename, const std::vector<Point3D>& vertices, const 
std::vector<Triangle>& triangles) {
    std::cout << "Writing STL file: " << filename << std::endl;
    std::ofstream f(filename, std::ios::binary);
    
    char header[80] = "GEBCO Metric STL";
    f.write(header, 80);
    
    uint32_t num_tri = triangles.size();
    f.write(reinterpret_cast<char*>(&num_tri), 4);

    for (const auto& t : triangles) {
        const Point3D& v0 = vertices[t.v1];
        const Point3D& v1 = vertices[t.v2];
        const Point3D& v2 = vertices[t.v3];

        // Calculate Normal
        double ax = v1.x - v0.x; double ay = v1.y - v0.y; double az = v1.z - v0.z;
        double bx = v2.x - v0.x; double by = v2.y - v0.y; double bz = v2.z - v0.z;
        
        float nx = ay * bz - az * by;
        float ny = az * bx - ax * bz;
        float nz = ax * by - ay * bx;
        
        f.write(reinterpret_cast<char*>(&nx), 4);
        f.write(reinterpret_cast<char*>(&ny), 4);
        f.write(reinterpret_cast<char*>(&nz), 4);

        // Write Vertices
        float p[3];
        p[0]=v0.x; p[1]=v0.y; p[2]=v0.z; f.write(reinterpret_cast<char*>(p), 12);
        p[0]=v1.x; p[1]=v1.y; p[2]=v1.z; f.write(reinterpret_cast<char*>(p), 12);
        p[0]=v2.x; p[1]=v2.y; p[2]=v2.z; f.write(reinterpret_cast<char*>(p), 12);

        uint16_t attr = 0;
        f.write(reinterpret_cast<char*>(&attr), 2);
    }
    f.close();
}

void write_gts(const std::string& filename, const std::vector<Point3D>& vertices, const 
std::vector<Triangle>& triangles) {
    std::cout << "Writing GTS file: " << filename << std::endl;
    std::ofstream f(filename);
    
    // For simplicity, using non-shared edges (triangle soup).
    // Proper GTS requires checking unique edges, but that requires complex bookkeeping.
    
    uint64_t num_edges = triangles.size() * 3;
    
    // Header
    f << vertices.size() << " " << num_edges << " " << triangles.size() << "\n";

    // Vertices
    f << std::fixed << std::setprecision(4);
    for(const auto& v : vertices) {
        f << v.x << " " << v.y << " " << v.z << "\n";
    }

    // Edges (1-based indexing)
    // Triangle 0 uses edges 1, 2, 3. Triangle 1 uses edges 4, 5, 6...
    for(size_t i = 0; i < triangles.size(); ++i) {
        const auto& t = triangles[i];
        // Edge 1: v1-v2
        f << (t.v1 + 1) << " " << (t.v2 + 1) << "\n";
        // Edge 2: v2-v3
        f << (t.v2 + 1) << " " << (t.v3 + 1) << "\n";
        // Edge 3: v3-v1
        f << (t.v3 + 1) << " " << (t.v1 + 1) << "\n";
    }

    // Triangles (1-based edge indexing)
    for(size_t i = 0; i < triangles.size(); ++i) {
        uint64_t base = i * 3 + 1; // Starting edge index (1-based)
        f << base << " " << (base + 1) << " " << (base + 2) << "\n";
    }
    f.close();
}

// --- MAIN ---

int main() {
    // --- USER PARAMETERS ---
    std::string nc_file = "GEBCO_2023.nc"; // Change this to your actual filename
    
    // Region of Interest (Spain/Portugal example)
    double lat_min = 35.0, lat_max = 45.0;
    double lon_min = -10.0, lon_max = 0.0;

    // Minimum Feature Size Filter (in meters)
    // Features smaller than this (e.g. small bumps) will be smoothed out.
    // Set to 0 to disable.
    double h_min_feature = 5000.0; // 5 km smoothing

    // URL (Placeholder - replace with actual GEBCO download link if you want 
auto-download)
    // Note: GEBCO requires registration usually. This URL is illustrative.
    std::string url = "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2023/zip/"; 

    // 1. Check and Download
    // Note: You should manually download the file first or place it in the binary dir.
    // This code checks existence.
    if (!download_file_if_missing(url, nc_file)) {
        std::cerr << "Failed to locate or download GEBCO file." << std::endl;
        return 1;
    }

    // 2. Extract Data
    std::vector<double> elevations;
    int width, height;
    double pixel_res_m = 0;

    std::cout << "Extracting region..." << std::endl;
    if (!extract_region(nc_file, lat_min, lat_max, lon_min, lon_max, elevations, width, 
height, pixel_res_m)) {
        return 1;
    }

    // 3. Apply Filter
    apply_size_filter(elevations, width, height, h_min_feature, pixel_res_m);

    // 4. Generate Mesh (Metric)
    std::vector<Point3D> vertices;
    std::vector<Triangle> triangles;
    
    std::cout << "Generating metric mesh..." << std::endl;
    generate_mesh(elevations, width, height, lat_min, lat_max, lon_min, lon_max, vertices, 
triangles);

    // 5. Write Outputs
    write_stl("output_metric.stl", vertices, triangles);
    write_gts("output_metric.gts", vertices, triangles);

    std::cout << "Process completed successfully." << std::endl;
    return 0;
}

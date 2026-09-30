#include "climate_erosion.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace fs=std::filesystem;
std::vector<double> read(const fs::path& path,size_t count) {
    if (fs::file_size(path)!=count*sizeof(double)) throw std::runtime_error("Wrong field size: "+path.string());
    std::vector<double> values(count);std::ifstream stream(path,std::ios::binary);
    stream.read(reinterpret_cast<char*>(values.data()),static_cast<std::streamsize>(count*sizeof(double)));
    if (!stream) throw std::runtime_error("Cannot read "+path.string());return values;
}
void write(const fs::path& path,const std::vector<double>& values) {
    std::ofstream stream(path,std::ios::binary);stream.write(reinterpret_cast<const char*>(values.data()),
        static_cast<std::streamsize>(values.size()*sizeof(double)));
    if (!stream) throw std::runtime_error("Cannot write "+path.string());
}
int main(int argc,char** argv) {
    try {
        if (argc!=9) throw std::runtime_error("Usage: tectonic_climate_erosion INPUT OUTPUT WIDTH HEIGHT MYR RIVER WAVE SEDIMENT_KM");
        const fs::path input=argv[1],output=argv[2];
        const auto w=static_cast<uint32_t>(std::stoul(argv[3])),h=static_cast<uint32_t>(std::stoul(argv[4]));
        const size_t count=static_cast<size_t>(w)*h;
        climateerosion::Parameters p;p.duration_myr=std::stod(argv[5]);p.river_coefficient=std::stod(argv[6]);
        p.wave_m_myr=std::stod(argv[7]);p.sediment_travel_km=std::stod(argv[8]);
        climateerosion::Forcing f{read(input/"rain.f64",count),read(input/"temperature.f64",count),
            read(input/"wind.f64",count),read(input/"snow.f64",count),read(input/"ice.f64",count)};
        const auto result=climateerosion::evolve(w,h,read(input/"elevation.f64",count),f,p);
        if (fs::exists(output)) throw std::runtime_error("Refusing to overwrite erosion output");
        fs::create_directories(output);
        write(output/"elevation.f64",result.elevation_m);write(output/"erosion.f64",result.erosion_m);
        write(output/"deposition.f64",result.deposition_m);
        std::ofstream budget(output/"budget.json");budget<<std::setprecision(17)
            <<"{\"eroded_m3\":"<<result.eroded_m3<<",\"deposited_m3\":"<<result.deposited_m3
            <<",\"relative_volume_residual\":"<<result.relative_volume_residual
            <<",\"duration_myr\":"<<p.duration_myr<<",\"step_myr\":"<<p.step_myr
            <<",\"ice_m_myr\":"<<p.ice_m_myr<<",\"wind_m_myr\":"<<p.wind_m_myr
            <<",\"hillslope_m2_myr\":"<<p.hillslope_m2_myr<<"}\n";
        if (std::abs(result.relative_volume_residual)>1e-9) throw std::runtime_error("Surface erosion volume budget failed");
        return 0;
    } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}

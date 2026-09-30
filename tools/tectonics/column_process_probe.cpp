// Isolated operator experiments, not additional numbered planetary runs.
#include "tectonic_columns.hpp"
#include "tectonic_thermal.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc,char** argv) {
    using namespace tectonics::columns;
    if(argc!=2 || std::filesystem::exists(argv[1])) {
        std::cerr<<"Usage: tectonic_column_probe NEW_OUTPUT.json\n";return 1;
    }
    std::ofstream out(argv[1]);if(!out)return 1;
    constexpr uint32_t w=64,h=8;
    constexpr double pi=3.141592653589793,cell_m=156250;
    out<<std::setprecision(15)<<"{\"schema\":\"tectonic-column-probe/v1\",\"seed\":1,\"cell_m\":"<<cell_m<<",\"cases\":[";
    bool first_case=true;
    for(int kind=0;kind<8;++kind) {
        std::vector<Cell> initial(w*h);
        const bool land=kind>=5;
        const double submarine=kind==3?.1:kind==4?.02:0;
        const double river=kind==6?5:kind==7?15:0;
        for(size_t i=0;i<initial.size();++i) {
            const double wave=std::sin(2*pi*(i%w)/8);
            initial[i]=land?seed(1000+7000*std::pow(.5+.5*std::cos(2*pi*(i%w)/w),4),true,100):seed(-4000,false,40);
            if(kind==0 || kind==2 || kind==3 || kind==4) {
                initial[i].ocean+=3000*wave;
                initial[i].oceanic_age_volume=initial[i].ocean*40;
            }
            if(kind==1)initial[i].oceanic_age_volume=initial[i].ocean*(40+30*wave);
        }
        Model model;model.initialize(w,h,cell_m,initial,1,1);
        tectonics::margins::Parameters p;p.enabled=land;p.river_erosion=river;
        const double thermal=kind<=1?1:0;
        if(!first_case)out<<',';first_case=false;
        const char* names[]={"cooling_fixed_crust_relief","cooling_age_relief","ocean_no_erosion",
            "ocean_diffusion_0.1","ocean_diffusion_0.02","gravity_only","gravity_rivers_5","gravity_rivers_15"};
        out<<"{\"name\":\""<<names[kind]<<"\",\"thermal_strength\":"<<thermal
            <<",\"submarine_ratio\":"<<submarine<<",\"river_erosion\":"<<river<<",\"states\":[";
        bool first_state=true;
        for(int t=0;t<=600;++t) {
            if(t%100==0) {
                double mean=0,second=0,high=-1e99,low=1e99;
                for(const auto& c:model.cells()) { const double z=elevation_m(c,thermal);
                    mean+=z/(w*h);second+=z*z/(w*h);high=std::max(high,z);low=std::min(low,z); }
                if(!first_state)out<<',';first_state=false;
                out<<"{\"time_myr\":"<<t<<",\"mean_m\":"<<mean<<",\"std_m\":"<<std::sqrt(std::max(0.0,second-mean*mean))
                    <<",\"peak_m\":"<<high<<",\"minimum_m\":"<<low<<",\"mass_residual\":"<<model.budget().relative_mass_residual
                    <<",\"eroded_kg\":"<<model.budget().eroded_kg<<",\"deposited_kg\":"<<model.budget().deposited_kg
                    <<",\"root_redistribution_kg\":"<<model.budget().redistributed_continent_kg<<'}';
            }
            if(t<600)model.advance(std::vector<Displacement>(w*h),1,submarine>0?1:0,thermal,submarine,nullptr,p);
        }
        out<<"]}";
    }
    out<<"]}\n";std::cout<<"Wrote isolated cooling, diffusion, gravity and river probes\n";
}

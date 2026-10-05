#pragma once

#include "simulations/climate/climate_physics.hpp"
#include <algorithm>
#include <cmath>

namespace climatehydrology::detail
{
// Diagnostic parcel only: this does not modify stored water or heat. Dry
// adiabatic lifting conserves potential temperature. Above the LCL, integrate
// dT/dln(p) for a pseudoadiabat (Rd*T + Lv*rs) /
// (cp + Lv^2*rs*epsilon/(Rd*T^2)); see MetPy's documented moist_lapse equation.
// https://unidata.github.io/MetPy/latest/api/generated/metpy.calc.moist_lapse.html
// Use constant latent heat and the model's water/ice saturation curve.
// The compact Bolton LCL and four midpoint steps approximate a bulk layer,
// not a resolved sounding, entraining plume, or a CAPE calculation.
inline float liftedParcelTemperatureC(float temperatureC, float relativeHumidity,
    float initialPressureHpa, float targetPressureHpa)
{
    constexpr double rd=287.05, cp=1004.0, lv=2.5e6, epsilon=0.622;
    const double p0=std::max(10.0f,initialPressureHpa);
    const double target=std::clamp(static_cast<double>(targetPressureHpa),1.0,p0);
    if (target==p0) return temperatureC;
    const double initial=std::max(150.0,static_cast<double>(temperatureC)+273.15);
    const double rh=std::clamp(static_cast<double>(relativeHumidity),1.0e-6,1.0);
    const double gamma=std::log(rh)+17.67*temperatureC/(temperatureC+243.5);
    const double dewpoint=std::min(initial,243.5*gamma/(17.67-gamma)+273.15);
    const double lcl=1.0/(1.0/(dewpoint-56.0)+std::log(initial/dewpoint)/800.0)+56.0;
    const double pLcl=p0*std::pow(lcl/initial,cp/rd);
    if (target>=pLcl)
        return static_cast<float>(initial*std::pow(target/p0,rd/cp)-273.15);
    const double step=std::log(target/pLcl)/4.0;
    auto tendency=[&](double t,double logp) {
        const double pressure=std::exp(logp);
        const double vapour=std::min(pressure*.95,static_cast<double>(
            climatephysics::saturationVapourPressureHpa(static_cast<float>(t-273.15))));
        const double rs=epsilon*vapour/(pressure-vapour);
        return (rd*t+lv*rs)/(cp+lv*lv*rs*epsilon/(rd*t*t));
    };
    double t=lcl,logp=std::log(pLcl);
    for (int i=0;i<4;++i)
    {
        const double mid=t+.5*step*tendency(t,logp);
        t+=step*tendency(mid,logp+.5*step);
        logp+=step;
    }
    return static_cast<float>(t-273.15);
}
}

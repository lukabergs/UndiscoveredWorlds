"""Stable semantic paths shared by reference generation, migration and its index."""

from pathlib import PurePosixPath

SEASONS = {"jan": "1", "apr": "2", "jul": "3", "oct": "4"}

# Additional products retain a provider branch when two datasets measure the
# same quantity. OISST is the primary SST comparison; ERA5 wind is the baseline.
ADDITIONAL = {
    "oisst_sst": "sea_temp/{s}/s",
    "era5_sst": "sea_temp/{s}/s/era5",
    "glorys_thetao": "sea_temp/{s}/s/glorys",
    "era5_t2m": "air_temp/{s}/s/era5",
    "era5_d2m": "air_temp/dewpoint/{s}/s/era5",
    "era5_tcc": "cloud/fraction/{s}/column/era5",
    "era5_lcc": "cloud/fraction/{s}/low/era5",
    "era5_hcc": "cloud/fraction/{s}/high/era5",
    "era5_blh": "air/boundary_height/{s}/b/era5",
    "era5_sp": "pressure/absolute/{s}/s/era5",
    "era5_msnswrf": "energy/net_shortwave/{s}/s/era5",
    "era5_msnlwrf": "energy/net_longwave/{s}/s/era5",
    "era5_mtnswrf": "energy/net_shortwave/{s}/toa/era5",
    "era5_mtnlwrf": "energy/net_longwave/{s}/toa/era5",
    "era5_msshf": "energy/sensible/{s}/s/era5",
    "era5_mslhf": "energy/latent/{s}/s/era5",
    "era5_mer": "water/evaporation/{s}/s/era5",
    "era5_vimd": "moisture/convergence/{s}/column/era5_native",
    "glorys_uo": "ocean/current_east/{s}/s",
    "glorys_vo": "ocean/current_north/{s}/s",
    "glorys_mlotst": "ocean/mixed_layer_depth/{s}/s",
    "era5land_ro": "water/runoff/{s}/total/era5land",
    "era5land_sro": "water/runoff/{s}/s/era5land",
    "era5land_ssro": "water/runoff/{s}/subsurface/era5land",
    "era5land_e": "water/evaporation/{s}/s/era5land",
    "era5land_sd": "snow/water_equivalent/{s}/s",
    "era5land_snow_depth": "snow/depth/{s}/s",
    "era5land_snowc": "snow/cover/{s}/s",
    "ceres_cloud_fraction": "cloud/fraction/{s}/column/ceres",
    "ceres_cloud_pressure": "cloud/pressure/{s}/effective/ceres",
    "ceres_cloud_temp": "cloud/temperature/{s}/effective/ceres",
    "ceres_cloud_tau": "cloud/optical_depth/{s}/daytime/ceres",
    "ceres_solar": "energy/solar/{s}/toa/ceres",
}
for _i in range(1, 5):
    ADDITIONAL[f"era5land_swvl{_i}"] = f"soil/water/{{s}}/layer{_i}"
for _level in ("toa", "sfc"):
    for _quantity in ("sw_all", "lw_all", "net_all", "sw_down_all", "sw_up_all",
                      "lw_down_all", "lw_up_all", "net_sw_all", "net_lw_all",
                      "net_tot_all", "cre_sw", "cre_lw", "cre_net", "cre_net_sw",
                      "cre_net_lw", "cre_net_tot"):
        ADDITIONAL[f"ceres_{_level}_{_quantity}"] = (
            f"energy/{_quantity}/{{s}}/{'s' if _level == 'sfc' else 'toa'}/ceres")


def product_directory(name):
    annual = {"koppen": "koppen", "temperature": "air_temp/annual/s",
              "precipitation": "rain/annual/s"}
    if name in annual:
        return PurePosixPath("climate", annual[name])
    month, quantity = name.split("_", 1)
    season = SEASONS[month]
    if quantity in ADDITIONAL:
        directory = ADDITIONAL[quantity].format(s=season)
    elif "_wind_" in quantity:
        layer, style = quantity.split("_wind_", 1)
        layer = {"surface": "s", "upper": "u", "850hpa": "850"}[layer]
        directory = f"wind/{style}/{season}/{layer}"
    elif quantity.endswith("_divergence"):
        layer = {"surface": "s", "upper": "u", "850hpa": "850"}[quantity.removesuffix("_divergence")]
        directory = f"wind/divergence/{season}/{layer}"
    else:
        directory = {
            "column_water": "moisture/water/{s}/column",
            "pressure_anomaly": "pressure/anomaly/{s}/slp",
            "ascent": "wind/ascent/{s}/u",
            "era5_precipitation": "rain/{s}/s/era5",
            "column_moisture_flux": "moisture/flux/{s}/column",
            "column_moisture_flux_east": "moisture/flux_east/{s}/column",
            "column_moisture_flux_north": "moisture/flux_north/{s}/column",
            "moisture_flux_convergence": "moisture/convergence/{s}/column",
            "ocean_current_speed": "ocean/current/{s}/s",
        }[quantity].format(s=season)
    return PurePosixPath("climate", directory)


def product_path(name, width, kind="maps", category="climate"):
    suffix = {"maps": ".png", "fields": ".tif", "csv": ".csv"}[kind]
    if category != "climate":
        return PurePosixPath(category, kind, f"earth_{name}_{width}x{width//2}{suffix}")
    return PurePosixPath(category, kind, product_directory(name), f"{width}{suffix}")

"""Compare end-to-end Hay-Davies POA accuracy and speed against pvlib."""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pandas as pd
import pvlib
import solars

from common import THREADS, average_duration, make_inputs, report_error_and_speed
from compare_clearsky_with_pvlib import extra_radiation

PANEL_TILT = 30.0
PANEL_AZIMUTH = 180.0
ALBEDO = 0.25


def calculate_with_solars(inputs: tuple[np.ndarray, ...]) -> dict[str, np.ndarray]:
    latitude, longitude, time, elevation, pressure, temperature = inputs
    solar_position = solars.calculate_solar_position(
        latitude, longitude, time, elevation, pressure, temperature, THREADS
    )
    clearsky = solars.calculate_clearsky(
        latitude, longitude, time, solar_position.zenith, elevation, pressure, THREADS
    )
    aoi = solars.calculate_aoi(
        solar_position.zenith,
        solar_position.azimuth,
        PANEL_TILT,
        PANEL_AZIMUTH,
        THREADS,
    )
    result = solars.calculate_poa(
        time,
        solar_position.zenith,
        aoi.aoi,
        PANEL_TILT,
        clearsky.dni,
        clearsky.ghi,
        clearsky.dhi,
        ALBEDO,
        THREADS,
    )
    return {
        "global": getattr(result, "global"),
        "direct": result.direct,
        "diffuse": result.diffuse,
        "sky_diffuse": result.sky_diffuse,
        "ground_diffuse": result.ground_diffuse,
    }


def calculate_with_pvlib(inputs: tuple[np.ndarray, ...]) -> dict[str, np.ndarray]:
    latitude, longitude, time, elevation, pressure, temperature = inputs
    times = pd.DatetimeIndex(time).tz_localize("UTC")
    shape = (time.size, latitude.size, longitude.size)
    values = {
        name: np.empty(shape, dtype=np.float64)
        for name in ("global", "direct", "diffuse", "sky_diffuse", "ground_diffuse")
    }
    dni_extra = np.array([extra_radiation(value.dayofyear) for value in times])
    for latitude_index, latitude_value in enumerate(latitude):
        for longitude_index, longitude_value in enumerate(longitude):
            solar_position = pvlib.solarposition.get_solarposition(
                times,
                latitude_value,
                longitude_value,
                altitude=elevation[latitude_index, longitude_index],
                pressure=cast(Any, pressure[:, latitude_index, longitude_index] * 100.0),
                temperature=cast(Any, temperature[:, latitude_index, longitude_index]),
                method="nrel_numpy",
            )
            turbidity = pvlib.clearsky.lookup_linke_turbidity(
                times, latitude_value, longitude_value
            )
            airmass = pvlib.atmosphere.get_relative_airmass(
                solar_position["apparent_zenith"], model="kastenyoung1989"
            ) * pressure[:, latitude_index, longitude_index] / 1013.25
            clearsky = pvlib.clearsky.ineichen(
                solar_position["apparent_zenith"],
                airmass,
                turbidity,
                altitude=elevation[latitude_index, longitude_index],
                dni_extra=cast(Any, dni_extra),
            )
            poa = pvlib.irradiance.get_total_irradiance(
                PANEL_TILT,
                PANEL_AZIMUTH,
                solar_position["apparent_zenith"],
                solar_position["azimuth"],
                clearsky["dni"],
                clearsky["ghi"],
                clearsky["dhi"],
                dni_extra=dni_extra,
                albedo=ALBEDO,
                model="haydavies",
            )
            values["global"][:, latitude_index, longitude_index] = poa["poa_global"]
            values["direct"][:, latitude_index, longitude_index] = poa["poa_direct"]
            values["diffuse"][:, latitude_index, longitude_index] = poa["poa_diffuse"]
            values["sky_diffuse"][:, latitude_index, longitude_index] = poa["poa_sky_diffuse"]
            values["ground_diffuse"][:, latitude_index, longitude_index] = poa["poa_ground_diffuse"]
    return values


def main() -> None:
    inputs = make_inputs()
    solars_values = calculate_with_solars(inputs)
    pvlib_values = calculate_with_pvlib(inputs)
    report_error_and_speed(
        "Plane of array",
        solars_values,
        pvlib_values,
        average_duration(lambda: calculate_with_solars(inputs)),
        average_duration(lambda: calculate_with_pvlib(inputs)),
    )


if __name__ == "__main__":
    main()

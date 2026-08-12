"""Compare end-to-end Ineichen clear-sky accuracy and speed against pvlib."""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pandas as pd
import pvlib
import solars

from common import THREADS, average_duration, make_inputs, report_error_and_speed


def extra_radiation(day_of_year: int) -> float:
    """To mimic rust implementation"""
    x = 2.0 * np.pi * (day_of_year - 1) / 365.0
    return 1366.1 * (
        1.00011
        + 0.034221 * np.cos(x)
        + 0.00128 * np.sin(x)
        - 0.000719 * np.cos(2.0 * x)
        + 0.000077 * np.sin(2.0 * x)
    )


def calculate_with_solars(inputs: tuple[np.ndarray, ...]) -> dict[str, np.ndarray]:
    latitude, longitude, time, elevation, pressure, temperature = inputs
    solar_position = solars.calculate_solar_position(
        latitude, longitude, time, elevation, pressure, temperature, THREADS
    )
    result = solars.calculate_clearsky(
        latitude, longitude, time, solar_position.zenith, elevation, pressure, THREADS
    )
    return {"ghi": result.ghi, "dni": result.dni, "dhi": result.dhi}


def calculate_with_pvlib(inputs: tuple[np.ndarray, ...]) -> dict[str, np.ndarray]:
    latitude, longitude, time, elevation, pressure, temperature = inputs
    times = pd.DatetimeIndex(time).tz_localize("UTC")
    shape = (time.size, latitude.size, longitude.size)
    ghi = np.empty(shape, dtype=np.float64)
    dni = np.empty(shape, dtype=np.float64)
    dhi = np.empty(shape, dtype=np.float64)
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
            ghi[:, latitude_index, longitude_index] = clearsky["ghi"]
            dni[:, latitude_index, longitude_index] = clearsky["dni"]
            dhi[:, latitude_index, longitude_index] = clearsky["dhi"]
    return {"ghi": ghi, "dni": dni, "dhi": dhi}


def main() -> None:
    inputs = make_inputs()
    solars_values = calculate_with_solars(inputs)
    pvlib_values = calculate_with_pvlib(inputs)
    report_error_and_speed(
        "Clear sky",
        solars_values,
        pvlib_values,
        average_duration(lambda: calculate_with_solars(inputs)),
        average_duration(lambda: calculate_with_pvlib(inputs)),
    )


if __name__ == "__main__":
    main()

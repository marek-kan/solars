"""Compare solar-position accuracy and speed against pvlib."""

from __future__ import annotations

from typing import Any, cast

import numpy as np
import pandas as pd
import pvlib
import solars

from common import THREADS, average_duration, make_inputs, report_error_and_speed


def calculate_with_solars(inputs: tuple[np.ndarray, ...]) -> dict[str, np.ndarray]:
    latitude, longitude, time, elevation, pressure, temperature = inputs
    result = solars.calculate_solar_position(
        latitude, longitude, time, elevation, pressure, temperature, THREADS
    )
    return {"zenith": result.zenith, "azimuth": result.azimuth}


def calculate_with_pvlib(inputs: tuple[np.ndarray, ...]) -> dict[str, np.ndarray]:
    latitude, longitude, time, elevation, pressure, temperature = inputs
    times = pd.DatetimeIndex(time).tz_localize("UTC")
    shape = (time.size, latitude.size, longitude.size)
    zenith = np.empty(shape, dtype=np.float64)
    azimuth = np.empty(shape, dtype=np.float64)
    for latitude_index, latitude_value in enumerate(latitude):
        for longitude_index, longitude_value in enumerate(longitude):
            result = pvlib.solarposition.get_solarposition(
                times,
                latitude_value,
                longitude_value,
                altitude=elevation[latitude_index, longitude_index],
                pressure=cast(Any, pressure[:, latitude_index, longitude_index] * 100.0),
                temperature=cast(Any, temperature[:, latitude_index, longitude_index]),
                method="nrel_numpy",
            )
            zenith[:, latitude_index, longitude_index] = result["apparent_zenith"]
            azimuth[:, latitude_index, longitude_index] = result["azimuth"]
    return {"zenith": zenith, "azimuth": azimuth}


def main() -> None:
    inputs = make_inputs()
    solars_values = calculate_with_solars(inputs)
    pvlib_values = calculate_with_pvlib(inputs)
    report_error_and_speed(
        "Solar position",
        solars_values,
        pvlib_values,
        average_duration(lambda: calculate_with_solars(inputs)),
        average_duration(lambda: calculate_with_pvlib(inputs)),
    )


if __name__ == "__main__":
    main()

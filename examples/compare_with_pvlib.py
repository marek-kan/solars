"""Compare solars against pvlib over a small Europe grid.

Install the optional dependencies and build the editable extension first:

    python -m pip install -e '.[comparison]'
    maturin develop

pvlib accepts one location per call, whereas solars evaluates the entire grid
in one call. Pressure is converted from solars' millibars to pvlib's Pascals.
"""

from __future__ import annotations

from time import perf_counter
from typing import Any, cast

import numpy as np
import pandas as pd
import pvlib
import solars

THREADS = 4
REPEATS = 5000
MAX_ACCEPTABLE_DIFFERENCE_DEGREES = 0.01


def make_inputs() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Return a 2x2 grid and daylight inputs for 12 hours over central Europe."""
    latitude = np.array([45.4642, 48.8566], dtype=np.float64)
    longitude = np.array([9.1900, 2.3522], dtype=np.float64)
    time = np.datetime64("2024-06-21T05:00:00", "ns") + np.arange(48) * np.timedelta64(15, "m")
    elevation = np.array([[120.0, 180.0], [35.0, 55.0]], dtype=np.float64)

    time_variation = np.linspace(0.0, 1.0, time.size, dtype=np.float64)[:, None, None]
    pressure = 1013.25 - 8.0 * time_variation + np.array([[0.0, 2.0], [4.0, 6.0]])
    temperature = 18.0 + 6.0 * time_variation + np.array([[0.0, 0.5], [1.0, 1.5]])

    return latitude, longitude, time, elevation, pressure, temperature


def calculate_with_solars(
    latitude: np.ndarray,
    longitude: np.ndarray,
    time: np.ndarray,
    elevation: np.ndarray,
    pressure: np.ndarray,
    temperature: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    result = solars.calculate_solar_position(
        latitude,
        longitude,
        time,
        elevation,
        pressure,
        temperature,
        THREADS,
    )
    return result.zenith, result.azimuth


def calculate_with_pvlib(
    latitude: np.ndarray,
    longitude: np.ndarray,
    time: np.ndarray,
    elevation: np.ndarray,
    pressure: np.ndarray,
    temperature: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate each grid cell because pvlib's API takes scalar locations."""
    times = pd.DatetimeIndex(time).tz_localize("UTC")
    shape = (time.size, latitude.size, longitude.size)
    zenith = np.empty(shape, dtype=np.float64)
    azimuth = np.empty(shape, dtype=np.float64)

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
            zenith[:, latitude_index, longitude_index] = solar_position["apparent_zenith"]
            azimuth[:, latitude_index, longitude_index] = solar_position["azimuth"]

    return zenith, azimuth


def average_duration(function, *args) -> float:
    """Measure repeated calls after the first invocation has warmed local state."""
    function(*args)
    started_at = perf_counter()
    for _ in range(REPEATS):
        function(*args)
    return (perf_counter() - started_at) / REPEATS


def main() -> None:
    inputs = make_inputs()
    solars_zenith, solars_azimuth = calculate_with_solars(*inputs)
    pvlib_zenith, pvlib_azimuth = calculate_with_pvlib(*inputs)

    zenith_difference = np.abs(solars_zenith - pvlib_zenith)
    azimuth_difference = np.abs(solars_azimuth - pvlib_azimuth)
    max_difference = max(zenith_difference.max(), azimuth_difference.max())

    solars_duration = average_duration(calculate_with_solars, *inputs)
    pvlib_duration = average_duration(calculate_with_pvlib, *inputs)

    print(f"Grid shape: {solars_zenith.shape} (time, lat, lon)")
    print(
        "Zenith difference (degrees): "
        f"mean={zenith_difference.mean():.6f}, max={zenith_difference.max():.6f}"
    )
    print(
        "Azimuth difference (degrees): "
        f"mean={azimuth_difference.mean():.6f}, max={azimuth_difference.max():.6f}"
    )
    print(f"solars: {solars_duration * 1_000:.3f} ms per run")
    print(f"pvlib:  {pvlib_duration * 1_000:.3f} ms per run")
    print(f"Speedup: {pvlib_duration / solars_duration:.2f}x")

    if max_difference > MAX_ACCEPTABLE_DIFFERENCE_DEGREES:
        raise AssertionError(
            "solars and pvlib differ by more than "
            f"{MAX_ACCEPTABLE_DIFFERENCE_DEGREES} degrees (max={max_difference:.6f})"
        )


if __name__ == "__main__":
    main()
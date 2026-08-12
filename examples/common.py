"""Shared inputs and timing helpers for pvlib comparison scripts."""

from __future__ import annotations

from time import perf_counter
from typing import Callable

import numpy as np

THREADS = 4
REPEATS = 20
TIME_COUNT = 100
GRID_SIZE = 10


def make_inputs() -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Return a (TIME_COUNT, GRID_SIZE, GRID_SIZE) central-Europe grid."""
    latitude = np.linspace(45.0, 55.0, GRID_SIZE, dtype=np.float64)
    longitude = np.linspace(2.0, 16.0, GRID_SIZE, dtype=np.float64)

    time = np.datetime64("2024-06-21T04:00:00", "ns") + np.arange(TIME_COUNT) * np.timedelta64(10, "m")

    elevation = 50.0 + np.add.outer(
        np.linspace(0.0, 900.0, GRID_SIZE), np.linspace(0.0, 200.0, GRID_SIZE)
    )

    time_variation = np.linspace(0.0, 1.0, TIME_COUNT, dtype=np.float64)[:, None, None]

    pressure = 1013.25 - 18.0 * time_variation - elevation[None, :, :] / 100.0

    temperature = 16.0 + 8.0 * time_variation - elevation[None, :, :] / 200.0

    return latitude, longitude, time, elevation, pressure, temperature


def average_duration(function: Callable[[], object]) -> float:
    """Measure repeated runs after a warm-up invocation."""
    function()
    started_at = perf_counter()
    for _ in range(REPEATS):
        function()
    return (perf_counter() - started_at) / REPEATS


def report_error_and_speed(
    name: str,
    solars_values: dict[str, np.ndarray],
    pvlib_values: dict[str, np.ndarray],
    solars_duration: float,
    pvlib_duration: float,
) -> None:
    """Print component-wise absolute errors and average execution durations."""
    print(f"{name}: shape={next(iter(solars_values.values())).shape} (time, lat, lon)")

    for component, solars_value in solars_values.items():
        difference = np.abs(solars_value - pvlib_values[component])
        print(f"{component}: mean_error={difference.mean():.6f}, max_error={difference.max():.6f}")

    print(f"solars: {solars_duration * 1_000:.3f} ms per run")
    print(f"pvlib:  {pvlib_duration * 1_000:.3f} ms per run")
    print(f"speedup: {pvlib_duration / solars_duration:.2f}x")

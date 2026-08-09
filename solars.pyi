from typing import Union

import numpy as np
from numpy.typing import NDArray

Float64Array = NDArray[np.float64]
Datetime64Array = NDArray[np.datetime64]
SpatialInput = Union[float, Float64Array]


class SolarPositionResult:
    @property
    def zenith(self) -> Float64Array: ...

    @property
    def azimuth(self) -> Float64Array: ...


class AoiResult:
    @property
    def aoi(self) -> Float64Array: ...


def calculate_solar_position(
    latitude: Float64Array,
    longitude: Float64Array,
    time: Datetime64Array,
    elevation: SpatialInput,
    pressure: Float64Array,
    temperature: Float64Array,
    num_threads: int,
) -> SolarPositionResult: ...


def calculate_aoi(
    zenith: Float64Array,
    azimuth: Float64Array,
    panel_tilt: SpatialInput,
    panel_azimuth: SpatialInput,
    num_threads: int,
) -> AoiResult: ...
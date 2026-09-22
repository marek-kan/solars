from typing import Union

import numpy as np
from numpy.typing import NDArray

__all__ = (
    "SolarPositionResult",
    "AoiResult",
    "ClearSkyResult",
    "PoaResult",
    "calculate_solar_position",
    "calculate_aoi",
    "calculate_clearsky",
    "calculate_poa",
)

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


class ClearSkyResult:
    @property
    def ghi(self) -> Float64Array: ...

    @property
    def dni(self) -> Float64Array: ...

    @property
    def dhi(self) -> Float64Array: ...


class PoaResult:
    @property
    def direct(self) -> Float64Array: ...

    @property
    def diffuse(self) -> Float64Array: ...

    @property
    def sky_diffuse(self) -> Float64Array: ...

    @property
    def ground_diffuse(self) -> Float64Array: ...


def calculate_solar_position(
    latitude: Float64Array,
    longitude: Float64Array,
    time: Datetime64Array,
    elevation: SpatialInput,
    pressure: Float64Array | None = None,
    temperature: Float64Array | None = None,
    num_threads: int = 1,
) -> SolarPositionResult: 
    """
    Compute solar zenith and azimuth arrays shaped `(time, lat, lon)`.

    `latitude`, `longitude` are one-dimensional degree axes;
    `time` is a one-dimensional `numpy.datetime64[ns]` array. 
    `elevation` is scalar or `(lat, lon)`.
    `pressure` and `temperature` are optional `(time)` or `(time, lat, lon)`
    arrays. Omit either to skip refraction correction.
    """
    ...


def calculate_aoi(
    zenith: Float64Array,
    azimuth: Float64Array,
    panel_tilt: SpatialInput,
    panel_azimuth: SpatialInput,
    num_threads: int,
    apply_optical_loss: bool = False,
    refractive_index: float = 1.526,
    extinction_coefficient: float = 4.0,
    thickness: float = 0.002,
) -> AoiResult:
    """
    Compute angle of incidence from precomputed solar-position arrays.
    
    `zenith` and azimuth use `(time, lat, lon)`. 
    `panel_tilt` and `panel_azimuth` are scalar or `(lat, lon)` degree arrays and are constant over time.
    Set `apply_optical_loss=True` to return physical IAM instead of geometric AOI.
    """
    ...


def calculate_clearsky(
    latitude: Float64Array,
    longitude: Float64Array,
    time: Datetime64Array,
    zenith: Float64Array,
    elevation: SpatialInput,
    pressure: Float64Array | None = None,
    num_threads: int = 1,
) -> ClearSkyResult: 
    """
    Compute Ineichen/Perez clear-sky irradiance arrays shaped `(time, lat, lon)`.
    
    `time` is a one-dimensional `numpy.datetime64[ns]` array. 
    `zenith` uses `(time, lat, lon)`.
    `elevation` is scalar or `(lat, lon)`.
    `pressure` is optional, in millibars, and may be `(time)` or `(time, lat, lon)`.
    """
    ...


def calculate_poa(
    time: Datetime64Array,
    zenith: Float64Array,
    aoi: Float64Array,
    panel_tilt: SpatialInput,
    dni: Float64Array,
    ghi: Float64Array,
    dhi: Float64Array,
    albedo: Float64Array,
    num_threads: int = 1,
) -> PoaResult: 
    """
    Compute Hay-Davies plane-of-array irradiance arrays shaped `(time, lat, lon)`.
    
    Zenith, AOI, DNI, GHI, and DHI use `(time, lat, lon)`. 
    Panel tilt is scalar or `(lat, lon)`.
    Albedo is `(time)` or `(time, lat, lon)`. For default behavior you can use 0.25
    `time` derives extraterrestrial DNI.
    """
    ...
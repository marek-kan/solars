# TODO

1) refractive/Optical loss in aoi, thermal loss in python (xarray)
2) clearsky GHI
3) POA

# How to compile

`RUSTFLAGS="-C target-cpu=x86-64-v3" maturin develop --release`
It can improve code generation, but it would make a distributed wheel fail on older CPUs. For a local deployment, it can be a benchmark configuration:

# Refactor

Understand and refactor `ObserverLatitudeGeometry` and `grid.rs` 

# Clearsky

Use ineichen model, TL from `LinkeTurbidities.h5`, AM calculates as `am = 1 / (cos(zenith_rad) + 0.50572 * ((06.07995 - zenith_angle_deg)**-1.6364))`
for I_0 use [book](https://www.osti.gov/servlets/purl/1039404) equations 12, 13 with x being `(2pi(n-1))/365`

# POA

POA = E_b + E_g + E_d

Component Calculations:

1. POA Beam ( E_b )
Calculated by projecting the Direct Normal Irradiance (DNI) onto the module plane using the angle of incidence (AOI):
E_b = DNI×cos(AOI)

2. POA Ground Reflected ( E_g )
Estimated using Global Horizontal Irradiance (GHI), ground albedo ( ρ ), and the module tilt angle ( β ):
E_g = GHI × ρ × (1−cos(β))/2

3. POA Sky Diffuse ( E_d  )
Computed by transposing the Diffuse Horizontal Irradiance (DHI) from the horizontal plane to the tilted module surface. Common industry models for this transposition include the Perez model (bankable standard) and the Hay model. pvlib uses isotropic model.

 - isotropic: Assumes diffuse radiation is uniformly distributed across the sky dome.  Simple but often less accurate.
 - haydavies: The Hay & Davies (1980) model, which accounts for circumsolar brightening. 
 - perez: The Perez et al. (1987/1990) model, considered the industry standard for bankable simulations. It accounts for circumsolar brightening and horizon brightening. 
 - perez-driesse: A continuous version of the Perez model that removes discontinuities found in the standard implementation. 
 - reindl: The Reindl et al. (1990) model, an enhancement of Hay-Davies. 
 - klucher: The Klucher (1979) model, which modifies the isotropic model for cloudy conditions.
 
 There is no single "best" model for all scenarios; the optimal choice depends on the specific application, data resolution, and orientation.

### **Perez** (Industry Standard for Bankability)
The **Perez** model is widely considered the **benchmark for bankable energy yield assessments** and is the default in major commercial software like PVsyst. It offers the highest accuracy for standard south-facing arrays by accounting for circumsolar and horizon brightening. However, it relies on a discrete look-up table for sky clearness, which introduces **mathematical discontinuities** (small jumps) in the output, particularly noticeable in sub-hourly simulations or reverse transposition tasks.

### **Perez-Driesse** (Best for High-Resolution & Optimization)
The **Perez-Driesse** model is a modern reformulation that replaces the discrete look-up table with continuous quadratic splines. It maintains the **same high accuracy** as the standard Perez model but eliminates discontinuities. This makes it the **superior choice** for:
*   **Sub-hourly simulations** (e.g., minute-resolution data) where smooth transitions are critical.
*   **Reverse transposition** (deriving horizontal irradiance from tilted measurements), where the discrete nature of the standard Perez model often causes solution failures.
*   **Optimization algorithms** that require continuous derivatives.

### **Hay & Davies** (Best for Stability & Low Bias)
The **Hay & Davies** model is often found to have the **lowest mean bias error (MBE)** and is statistically the "least risky" model, particularly for **non-south orientations** (East/West) and vertical surfaces. While it may slightly under-predict total irradiance compared to Perez in some clear-sky conditions (by ~1-2%), it is more robust when input data (specifically Diffuse Horizontal Irradiance) has higher uncertainty. It is an excellent choice for conservative estimates or when data quality is questionable.

### **Reindl**
The **Reindl** model performs very similarly to Hay & Davies, often trading places with it as the most accurate model for East/West orientations. It is a robust alternative but generally does not offer a distinct advantage over Hay & Davies for general utility-scale modeling.

### Summary Recommendation
*   **For Bankable Yield Reports:** Use **Perez** (standard industry practice) or **Perez-Driesse** (technically superior, gaining acceptance).
*   **For Research, Sub-hourly Data, or Optimization:** Use **Perez-Driesse**.
*   **For East/West Facades or Conservative Estimates:** Use **Hay & Davies**.

When using forecast data like **ECMWF IFS**, **Hay & Davies** is often a **superior practical choice** to Perez, despite Perez being the theoretical benchmark for measured data.

### Why Hay & Davies is Preferred for Forecasts
1.  **Lower Bias with Derived Inputs**: ECMWF IFS provides GHI and often DNI, but Diffuse Horizontal Irradiance (DHI) is frequently *derived* (calculated via a separation model) rather than directly forecasted with high precision. Studies show that when DHI is uncertain or derived, the complex **Perez** model tends to accumulate errors and exhibits a **positive bias** (overestimating energy by ~1-2% annually). **Hay & Davies** is more robust to this uncertainty and typically yields **lower Mean Bias Error (MBE)** in these conditions.
2.  **Stability**: Forecast data inherently contains noise and smoothing artifacts. The simpler geometry of Hay & Davies makes it less sensitive to small fluctuations in the input variables compared to the complex sky-clearness indexing of Perez.
3.  **Industry Practice for NWP**: Many operational forecasting chains using Numerical Weather Prediction (NWP) data default to Hay & Davies (or Reindl) specifically to minimize systematic overestimation risks, which are critical for grid balancing and energy trading.

### Recommendation
*   **Opt for Hay & Davies** if your priority is **minimizing bias** and ensuring conservative, stable predictions, especially if you are deriving DHI from GHI/DNI.
*   **Consider Perez** only if you have high-confidence, direct forecasts for all three components (GHI, DNI, DHI) and require the specific horizon-brightening physics for a specific site calibration, though the marginal gain is often negligible compared to the inherent uncertainty of the weather model itself.

---

```python
import numpy as np

def calculate_poa_hay_davies(
    ghi: np.ndarray,
    dhi: np.ndarray,
    dni: np.ndarray,
    solar_zenith_deg: np.ndarray,
    aoi_deg: np.ndarray,
    surface_tilt_deg: float,
    dni_extra: np.ndarray = 1361.0,
    albedo: float = 0.2,
    zenith_threshold_deg: float = 85.0
) -> dict:
    """
    Calculates Plane of Array (POA) irradiance using the Hay & Davies anisotropic model.

    Parameters:
    -----------
    ghi : np.ndarray
        Global Horizontal Irradiance (W/m²)
    dhi : np.ndarray
        Diffuse Horizontal Irradiance (W/m²)
    dni : np.ndarray
        Direct Normal Irradiance (W/m²)
    solar_zenith_deg : np.ndarray
        Solar zenith angle in degrees (0° = directly overhead)
    aoi_deg : np.ndarray
        Angle of incidence between solar vector and module surface normal in degrees
    surface_tilt_deg : float or np.ndarray
        Array tilt angle from horizontal in degrees (0° = flat)
    dni_extra : float or np.ndarray, optional
        Extraterrestrial normal irradiance (W/m²), defaults to solar constant 1361.0
    albedo : float, optional
        Ground reflectance factor (0.2 is standard grass/soil)
    zenith_threshold_deg : float, optional
        Zenith angle beyond which beam & circumsolar are forced to 0 (default 85.0°)

    Returns:
    --------
    dict containing arrays for:
        - 'poa_global': Total incident POA irradiance (W/m²)
        - 'poa_direct': Direct beam component (W/m²)
        - 'poa_diffuse': Sky diffuse component (W/m²)
        - 'poa_ground': Ground-reflected component (W/m²)
    """
    # Convert angles from degrees to radians
    zenith_rad = np.radians(solar_zenith_deg)
    aoi_rad = np.radians(aoi_deg)
    tilt_rad = np.radians(surface_tilt_deg)

    # Calculate trigonometric terms
    cos_zenith = np.cos(zenith_rad)
    cos_aoi = np.cos(aoi_rad)

    # Mask daytime / valid geometry conditions
    is_daylight = (solar_zenith_deg < zenith_threshold_deg) & (cos_zenith > 0)
    is_sun_facing = (cos_aoi > 0)

    # 1. Direct Beam Component
    # Beam is 0 if sun is below horizon threshold OR behind the module face
    poa_direct = np.where(is_daylight & is_sun_facing, dni * cos_aoi, 0.0)

    # 2. Geometric Projection Ratio (Rb = cos(AOI) / cos(Zenith))
    # Safeguard against division-by-zero near horizon
    rb = np.zeros_like(ghi)
    np.divide(cos_aoi, cos_zenith, out=rb, where=is_daylight)
    rb = np.maximum(rb, 0.0)  # Cannot be negative

    # 3. Anisotropy Index (Ai = DNI / DNI_extra)
    # Scaled atmospheric transmittance indicator constrained between [0, 1]
    anisotropy_index = np.clip(dni / dni_extra, 0.0, 1.0)

    # 4. Sky Diffuse Component (Hay & Davies equation)
    # Term A: Circumsolar diffuse component = DHI * Ai * Rb
    # Term B: Isotropic diffuse component   = DHI * (1 - Ai) * ((1 + cos(tilt)) / 2)
    view_factor_sky = (1.0 + np.cos(tilt_rad)) / 2.0
    
    circumsolar = dhi * anisotropy_index * rb
    isotropic_background = dhi * (1.0 - anisotropy_index) * view_factor_sky

    poa_diffuse = circumsolar + isotropic_background

    # 5. Ground-Reflected Component
    view_factor_ground = (1.0 - np.cos(tilt_rad)) / 2.0
    poa_ground = ghi * albedo * view_factor_ground

    # Total POA Irradiance
    poa_global = poa_direct + poa_diffuse + poa_ground

    return {
        "poa_global": poa_global,
        "poa_direct": poa_direct,
        "poa_diffuse": poa_diffuse,
        "poa_ground": poa_ground
    }
```
 Mathematical FormulationTotal 
 POA Irradiance ($E_{POA}$) is the sum of three distinct components: direct beam, sky diffuse, and ground-reflected radiation:
 
 $$E_{POA} = E_{beam} + E_{sky,diffuse} + E_{ground,reflected}$$
 
 ## 1. Direct Beam Component ($E_{beam}$)
 
 $$E_{beam} = DNI \cdot \cos(\theta)$$
 
 Where $DNI$ is Direct Normal Irradiance and $\theta$ is the Angle of Incidence (AOI) between the solar vector and panel normal.
 
 ## 2. Sky Diffuse Component ($E_{sky,diffuse}$) — Hay & Davies Model
 
 $$E_{sky,diffuse} = DHI \left[ A_i R_b + (1 - A_i) \left( \frac{1 + \cos\beta}{2} \right) \right]$$

Where:
$DHI$: Diffuse Horizontal Irradiance  
$A_i$: Anisotropy Index, defined as $A_i = \frac{DNI}{E_{extra}}$ (ratio of $DNI$ to extraterrestrial normal irradiance 
$E_{extra}$). $A_i$ represents atmospheric transmittance for beam radiation.
$R_b$: Direct beam geometric projection ratio, 
$R_b = \frac{\cos\theta}{\cos\theta_z}$ (where $\theta_z$ is the solar zenith angle).  
$\beta$: Array tilt angle from horizontal.  How the equation behaves:Clear Sky Conditions ($A_i \to 1$): 
Atmospheric transmittance is high. $E_{sky,diffuse}$ becomes dominated by the circumsolar term $DHI \cdot R_b$, matching beam geometry.Overcast Conditions ($A_i \to 0$): $DNI$ drops to near zero. The circumsolar portion vanishes, reducing the equation to the classic isotropic sky model $DHI \cdot \frac{1 + \cos\beta}{2}$.  
## 3. Ground-Reflected Component ($E_{ground,reflected}$)

$$E_{ground,reflected} = GHI \cdot \rho \cdot \left( \frac{1 - \cos\beta}{2} \right)$$

Where $GHI$ is Global Horizontal Irradiance and $\rho$ is the ground albedo.  

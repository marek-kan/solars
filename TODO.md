# How to compile

`RUSTFLAGS="-C target-cpu=x86-64-v3" maturin develop --release`
It can improve code generation, but it would make a distributed wheel fail on older CPUs. For a local deployment, it can be a benchmark configuration:

# Refactor

Understand and refactor `ObserverLatitudeGeometry` and `grid.rs` 

# Clearsky

Use ineichen model, TL from `LinkeTurbidities.h5`, AM calculates as `am = 1 / (cos(zenith_rad) + 0.50572 * ((06.07995 - zenith_angle_deg)**-1.6364))`
use [book](https://www.osti.gov/servlets/purl/1039404) equations 12, 13 with x being `(2pi(n-1))/365`

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

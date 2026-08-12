# Examples

The comparison scripts benchmark `solars` against `pvlib` on 100 timestamps and a 10x10 latitude/longitude grid. Each script prints mean and maximum absolute error for every returned component, the average run duration, and the resulting speedup.

## Setup

Install the extension and comparison dependencies in the project virtual environment:

```bash
source .venv/bin/activate
python -m pip install -e '.[comparison]'
maturin develop --release
```

For local benchmarking on CPUs that support it, compile with:

```bash
RUSTFLAGS="-C target-cpu=x86-64-v3" maturin develop --release
```

This flag improves code generation for a compatible local CPU, but the resulting wheel is not portable to older x86-64 processors.

## Scripts

Run scripts from this directory so their shared helper module is importable:

```bash
cd examples
python compare_solar_position_with_pvlib.py
python compare_clearsky_with_pvlib.py
python compare_poa_with_pvlib.py
```

`compare_solar_position_with_pvlib.py` compares apparent zenith and azimuth.

`compare_clearsky_with_pvlib.py` compares end-to-end Ineichen GHI, DNI, and DHI. Both paths use Kasten-Young relative airmass, the project extraterrestrial-radiation formula, and pvlib-compatible Linke-turbidity lookup.

`compare_poa_with_pvlib.py` compares the end-to-end solar-position, Ineichen clear-sky, AOI, and Hay-Davies POA pipeline, including direct, sky-diffuse, ground-diffuse, diffuse, and global POA irradiance.

## Results

The following results were captured with CPython 3.12, `pvlib 0.15.2`, `THREADS=4`, and a release build using `RUSTFLAGS="-C target-cpu=x86-64-v3"`:

| Calculation | solars | pvlib | Speedup |
| --- | ---: | ---: | ---: |
| Solar position | 0.836 ms | 130.815 ms | 156.41x |
| Clear sky | 1.408 ms | 443.275 ms | 314.86x |
| Plane of array | 2.060 ms | 615.923 ms | 298.99x |

The reported speedup is `pvlib duration / solars duration`; values above `1x` mean `solars` is faster for that workload. Exact timing varies with CPU, compilation flags, Python environment, and the configured number of threads.

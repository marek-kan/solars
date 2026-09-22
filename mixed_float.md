# Basic design

```rust

use std::ops::{Add, Sub, Mul, Div, Neg};

macro_rules! define_and_impl_float_trait {
    (
        pub(crate) trait $trait_name:ident {
            // Associated conversion functions taking a single argument
            $( fn $ctor:ident ( $c_arg:ident : $c_arg_ty:ty ) -> Self; )*
            ---
            // Instance methods (take self)
            $( fn $method:ident ( $( $arg:ident : $arg_ty:ty ),* ) -> $ret:ty; )*
        }
    ) => {
        // 1. Trait definition
        pub(crate) trait $trait_name:
            Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self>
            + Neg<Output = Self>
            + Sized + Copy + Clone + PartialEq + PartialOrd
        {
            $( fn $ctor($c_arg: $c_arg_ty) -> Self; )*
            $( fn $method(self, $( $arg: $arg_ty ),* ) -> $ret; )*
        }

        // 2. Implementation for f32
        impl $trait_name for f32 {
            $(
                #[inline(always)]
                fn $ctor($c_arg: $c_arg_ty) -> Self {
                    $c_arg as f32
                }
            )*
            $(
                #[inline(always)]
                fn $method(self, $( $arg: $arg_ty ),* ) -> $ret {
                    self.$method($( $arg ),*)
                }
            )*
        }

        // 3. Implementation for f64
        impl $trait_name for f64 {
            $(
                #[inline(always)]
                fn $ctor($c_arg: $c_arg_ty) -> Self {
                    $c_arg as f64
                }
            )*
            $(
                #[inline(always)]
                fn $method(self, $( $arg: $arg_ty ),* ) -> $ret {
                    self.$method($( $arg ),*)
                }
            )*
        }
    };
}

define_and_impl_float_trait! {
    pub(crate) trait MixedFloat {
        fn from_f32(val: f32) -> Self;
        fn from_f64(val: f64) -> Self;
        ---
        fn sqrt() -> Self;
        fn exp() -> Self;
        fn max(other: Self) -> Self;
        fn min(other: Self) -> Self;
        fn clamp(min: Self, max: Self) -> Self;
        fn powi(n: i32) -> Self;
        fn powf(n: Self) -> Self;
        fn sin() -> Self;
        fn cos() -> Self;
        fn acos() -> Self;
        fn to_degrees() -> Self;
        fn to_radians() -> Self;
        fn abs() -> Self;
    }
}

```

# Shortcomings

Searched for text `ndarray::`, 3 results

Searched for text `PyReadonlyArray`, 28 results

Not quite yet. While `MixedFloat` is a solid foundation, there are several missing traits, constants, and Python binding considerations before you can seamlessly pass `f32` arrays end-to-end.

Here is the checklist of what is missing:

---

### 1. Missing Trait Bounds & Methods on `MixedFloat`

1. **`ndarray::ScalarOperand`**:
   To perform array arithmetic like `array * 2.0` or `&array1 + &array2` with `ndarray`, `T` must implement `ndarray::ScalarOperand`.
2. **`numpy::Element`**:
   To return or accept NumPy arrays (`PyArray3<T>`, `PyReadonlyArrayDyn<T>`), `T` must implement `numpy::Element`.
3. **`num_traits::Float` or Mathematical Constants**:
   Constants like $\pi$, $\infty$, and integer casts:
   - `std::f64::consts::PI` $\rightarrow$ Needs `T::PI` or `T::from_f64(std::f64::consts::PI)`.
   - `f64::INFINITY` (used in `relative_airmass_kasten`) $\rightarrow$ Needs `T::INFINITY` or `fn infinity() -> Self`.
   - `atan2`, `asin`, `tan` (used in solar position and angle calculations) will need to be added to the trait macro if they aren't there.
4. **`Send + Sync + 'static`**:
   Required by Rayon (`par_iter_mut`, `par_chunks_mut`).
5. **`FromPrimitive` or generic integer casting**:
   `(day_of_year - 1) as f64` $\rightarrow$ Needs `T::from_i64(...)` or `T::from_f64(...)`.

---

### 2. Physical / Periodic Table Precision Considerations

- Tables in `periodic_tables` (like `LinkeTurbidities.h5`, `L0_TABLE`, `nutation`, earth coefficients) use double precision (`f64`).
- It is usually best to keep internal high-precision astronomical calculations in `f64`, or cast table outputs via `T::from_f64(...)` when interacting with generic pipeline arrays.

---

### 3. Recommended `MixedFloat` Trait Definition

```rust
use std::ops::{Add, Sub, Mul, Div, Neg};
use ndarray::ScalarOperand;
use numpy::Element;

macro_rules! define_and_impl_float_trait {
    (
        pub(crate) trait $trait_name:ident {
            const PI: Self;
            const INFINITY: Self;

            $( fn $ctor:ident ( $c_arg:ident : $c_arg_ty:ty ) -> Self; )*
            ---
            $( fn $method:ident ( $( $arg:ident : $arg_ty:ty ),* ) -> $ret:ty; )*
        }
    ) => {
        pub(crate) trait $trait_name:
            Add<Output = Self> + Sub<Output = Self> + Mul<Output = Self> + Div<Output = Self>
            + Neg<Output = Self>
            + Sized + Copy + Clone + PartialEq + PartialOrd + Send + Sync + 'static
            + ScalarOperand + Element
        {
            const PI: Self;
            const INFINITY: Self;

            $( fn $ctor($c_arg: $c_arg_ty) -> Self; )*
            $( fn $method(self, $( $arg: $arg_ty ),* ) -> $ret; )*
        }

        impl $trait_name for f32 {
            const PI: Self = std::f32::consts::PI;
            const INFINITY: Self = f32::INFINITY;

            $(
                #[inline(always)]
                fn $ctor($c_arg: $c_arg_ty) -> Self {
                    $c_arg as f32
                }
            )*
            $(
                #[inline(always)]
                fn $method(self, $( $arg: $arg_ty ),* ) -> $ret {
                    self.$method($( $arg ),*)
                }
            )*
        }

        impl $trait_name for f64 {
            const PI: Self = std::f64::consts::PI;
            const INFINITY: Self = f64::INFINITY;

            $(
                #[inline(always)]
                fn $ctor($c_arg: $c_arg_ty) -> Self {
                    $c_arg as f64
                }
            )*
            $(
                #[inline(always)]
                fn $method(self, $( $arg: $arg_ty ),* ) -> $ret {
                    self.$method($( $arg ),*)
                }
            )*
        }
    };
}
```

---

### 4. How Python Dispatches `f32` vs `f64` in PyO3

PyO3 functions cannot be generic directly in `#[pyfunction]`. In [src/lib.rs](src/lib.rs), you accept dynamic types or inspect the NumPy dtype (e.g., using `&Bound<'py, PyAny>`), then dispatch:

```rust
if let Ok(arr_f32) = array.extract::<PyReadonlyArrayDyn<f32>>() {
    // Run generic Rust function with T = f32
} else if let Ok(arr_f64) = array.extract::<PyReadonlyArrayDyn<f64>>() {
    // Run generic Rust function with T = f64
}
```
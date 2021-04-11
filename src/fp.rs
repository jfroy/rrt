//! Provides floating point utilities.

/// MACHINE_EPSILON is the numerical analysis definition of machine epsilon,
/// where it is interchangeble with unit roundoff. In C++ and Rust, the
/// built-in EPSILON constant is defined as the difference between 1 and
/// the next larger floating point number, i.e. twice the unit roundoff
/// or the traditional academic definition of machine epsilon. This constant
/// matches pbrt and most academic papers and literature.
///
/// See https://en.wikipedia.org/wiki/Machine_epsilon and
/// http://www.pbr-book.org/3ed-2018/Shapes/Managing_Rounding_Error.html#x1-ArithmeticOperations.
pub const MACHINE_EPSILON: f32 = f32::EPSILON * 0.5_f32;

/// gamma_eb computes a tight bound for products of (1 +/- machine epislon)
/// error terms. See http://www.pbr-book.org/3ed-2018/Shapes/Managing_Rounding_Error.html#x1-ErrorPropagation.
pub const fn gamma_eb(n: i32) -> f32 {
    let n = n as f32;
    (n * MACHINE_EPSILON) / (1_f32 - n * MACHINE_EPSILON)
}

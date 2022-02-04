use crate::impl_Interpolate;

use cgmath::{Quaternion, Vector1, Vector2, Vector3, Vector4, MetricSpace};
use crate::interpolate::InterpolateDerivative;
use cgmath::num_traits::Pow;

impl_Interpolate!(f32, Vector1<f32>, std::f32::consts::PI);
impl_Interpolate!(f32, Vector2<f32>, std::f32::consts::PI);
impl_Interpolate!(f32, Vector3<f32>, std::f32::consts::PI);
impl_Interpolate!(f32, Vector4<f32>, std::f32::consts::PI);
impl_Interpolate!(f32, Quaternion<f32>, std::f32::consts::PI);

impl_Interpolate!(f64, Vector1<f64>, std::f64::consts::PI);
impl_Interpolate!(f64, Vector2<f64>, std::f64::consts::PI);
impl_Interpolate!(f64, Vector3<f64>, std::f64::consts::PI);
impl_Interpolate!(f64, Vector4<f64>, std::f64::consts::PI);
impl_Interpolate!(f64, Quaternion<f64>, std::f64::consts::PI);

impl InterpolateDerivative<f64> for Vector2<f64> {
    fn cubic_hermite(t: f64, x: (f64, Self), a: (f64, Self), b: (f64, Self), y: (f64, Self)) -> Option<Self> {
        let (t0, p0) = x;
        let (_, p1) = a;
        let (_, p2) = b;
        let (_, p3) = y;

        let t1 = p0.distance(p1).pow(0.5) + t0;
        let t2 = p1.distance(p2).pow(0.5) + t1;
        let t3 = p2.distance(p3).pow(0.5) + t2;

        let a1 = (t1 - t) / (t1 - t0) * p0 + (t - t0) / (t1 - t0) * p1;
        let a2 = (t2 - t) / (t2 - t1) * p1 + (t - t1) / (t2 - t1) * p2;
        let a3 = (t3 - t) / (t3 - t2) * p2 + (t - t2) / (t3 - t2) * p3;

        let b1 = (t2 - t) / (t2 - t0) * a1 + (t - t0) / (t2 - t0) * a2;
        let b2 = (t3 - t) / (t3 - t1) * a2 + (t - t1) / (t3 - t1) * a3;

        let a1_der = 1f64 / (t1 - t0) * (p1 - p0);
        let a2_der = 1f64 / (t2 - t1) * (p2 - p1);
        let a3_der = 1f64 / (t3 - t2) * (p3 - p2);

        let b1_der_t = 1f64 / (t2 - t0) * (a2 - a1) + (t2 - t) / (t2 - t0) * a1_der + (t - t0) / (t2 - t0) * a2_der;
        let b2_der_t = 1f64 / (t3 - t1) * (a3 - a2) + (t3 - t) / (t3 - t1) * a2_der + (t - t1) / (t3 - t1) * a3_der;

        Some(
            1f64 / (t2 - t1) * (b2 - b1) + (t2 - t) / (t2 - t1) * b1_der_t + (t - t1) / (t2 - t1) * b2_der_t
        )
    }

    fn quadratic_bezier_derivative(t: f64, a: Self, u: Self, b: Self) -> Self {
        2. * (1. - t) * (u - a) + 2. * t * (b - u)
    }

    // a -> P0
    // u -> P1
    // v -> P2
    // b -> P3
    fn cubic_bezier_derivative(t: f64, a: Self, u: Self, v: Self, b: Self) -> Self {
        let one_t = 1. - t;
        let one_t2 = one_t * one_t;
        let t2 = t * t;

        3. * one_t2 * (u - a) + 6. * one_t * t * (v - u) + 3. * t2 * (b - v)
    }

    fn cubic_bezier_mirrored_derivative(t: f64, a: Self, u: Self, v: Self, b: Self) -> Self {
        Self::cubic_bezier_derivative(t, a, u, b + b - v, b)
    }
}
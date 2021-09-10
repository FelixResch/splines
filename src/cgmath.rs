use crate::impl_Interpolate;

use cgmath::{Quaternion, Vector1, Vector2, Vector3, Vector4};
use crate::interpolate::InterpolateDerivative;

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
        let (t1, p1) = a;
        let (t2, p2) = b;
        let (t3, p3) = y;

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
}
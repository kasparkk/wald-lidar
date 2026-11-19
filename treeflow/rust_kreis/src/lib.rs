//! Schnelle Teile von treeflow in Rust.
//!
//! `fit_circle` ist dieselbe Rechnung wie `treeflow.detect.fit_circle` in Python:
//! RANSAC-Kreis durch je 3 zufaellige Punkte, dann kleinste Quadrate auf den
//! Punkten, die zum besten Kreis passen. Unterschied: die Punkte werden vorher
//! auf ihren Mittelpunkt verschoben. Das aendert das Ergebnis mathematisch
//! nicht, vermeidet aber Rundungsfehler bei grossen UTM-Koordinaten (~400 000 m).

use numpy::PyReadonlyArray2;
use pyo3::prelude::*;

/// Kleiner, reproduzierbarer Zufallsgenerator (xorshift64*), damit gleiche
/// Eingaben immer gleiche Kreise liefern.
struct Zufall(u64);

impl Zufall {
    fn neu(seed: u64) -> Self {
        Zufall(seed ^ 0x9E37_79B9_7F4A_7C15 | 1)
    }
    fn naechste(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    fn index(&mut self, n: usize) -> usize {
        (self.naechste() % n as u64) as usize
    }
}

/// Kreis durch drei Punkte: (Mittelpunkt x, y, Radius) oder None, wenn sie auf einer Linie liegen.
fn kreis_durch_drei(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> Option<(f64, f64, f64)> {
    let d = 2.0 * (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1));
    if d.abs() < 1e-9 {
        return None;
    }
    let aa = a.0 * a.0 + a.1 * a.1;
    let bb = b.0 * b.0 + b.1 * b.1;
    let cc = c.0 * c.0 + c.1 * c.1;
    let ux = (aa * (b.1 - c.1) + bb * (c.1 - a.1) + cc * (a.1 - b.1)) / d;
    let uy = (aa * (c.0 - b.0) + bb * (a.0 - c.0) + cc * (b.0 - a.0)) / d;
    Some((ux, uy, (a.0 - ux).hypot(a.1 - uy)))
}

/// Kleinste Quadrate (algebraischer Kreis-Fit) ueber die Normalgleichungen 3x3.
fn ausgleich(q: &[(f64, f64)]) -> Option<(f64, f64, f64)> {
    // Unbekannte: cx, cy, k mit x^2 + y^2 = 2 cx x + 2 cy y + k
    let mut m = [[0.0f64; 3]; 3];
    let mut v = [0.0f64; 3];
    for &(x, y) in q {
        let zeile = [2.0 * x, 2.0 * y, 1.0];
        let rechts = x * x + y * y;
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] += zeile[i] * zeile[j];
            }
            v[i] += zeile[i] * rechts;
        }
    }
    let loesung = loese3(m, v)?;
    let (cx, cy, k) = (loesung[0], loesung[1], loesung[2]);
    let r = (k + cx * cx + cy * cy).max(0.0).sqrt();
    Some((cx, cy, r))
}

/// 3x3-Gleichungssystem mit Gauss und Spaltenpivot.
fn loese3(mut m: [[f64; 3]; 3], mut v: [f64; 3]) -> Option<[f64; 3]> {
    for s in 0..3 {
        let p = (s..3).max_by(|&i, &j| m[i][s].abs().total_cmp(&m[j][s].abs()))?;
        if m[p][s].abs() < 1e-12 {
            return None;
        }
        m.swap(s, p);
        v.swap(s, p);
        for z in s + 1..3 {
            let f = m[z][s] / m[s][s];
            for k in s..3 {
                m[z][k] -= f * m[s][k];
            }
            v[z] -= f * v[s];
        }
    }
    let mut x = [0.0f64; 3];
    for s in (0..3).rev() {
        let mut summe = v[s];
        for k in s + 1..3 {
            summe -= m[s][k] * x[k];
        }
        x[s] = summe / m[s][s];
    }
    Some(x)
}

/// RANSAC-Kreis-Fit. Gibt (cx, cy, r, rmse) zurueck oder None, wie die Python-Fassung.
#[pyfunction]
#[pyo3(signature = (p, iters=200, tol=0.02, seed=0))]
fn fit_circle(
    py: Python<'_>,
    p: PyReadonlyArray2<'_, f64>,
    iters: usize,
    tol: f64,
    seed: u64,
) -> PyResult<Option<(f64, f64, f64, f64)>> {
    let arr = p.as_array();
    if arr.ncols() < 2 {
        return Err(pyo3::exceptions::PyValueError::new_err("p braucht zwei Spalten (x, y)"));
    }
    let n = arr.nrows();
    if n < 6 {
        return Ok(None);
    }
    // Auf den Mittelpunkt verschieben (bessere Rechengenauigkeit).
    let (mx, my) = (
        arr.column(0).iter().sum::<f64>() / n as f64,
        arr.column(1).iter().sum::<f64>() / n as f64,
    );
    let pts: Vec<(f64, f64)> = arr.rows().into_iter().map(|r| (r[0] - mx, r[1] - my)).collect();

    let ergebnis = py.allow_threads(|| {
        let mut zufall = Zufall::neu(seed);
        let mut beste: Option<(usize, (f64, f64, f64))> = None;
        for _ in 0..iters {
            let i = zufall.index(n);
            let mut j = zufall.index(n);
            while j == i {
                j = zufall.index(n);
            }
            let mut k = zufall.index(n);
            while k == i || k == j {
                k = zufall.index(n);
            }
            let Some((ux, uy, r)) = kreis_durch_drei(pts[i], pts[j], pts[k]) else { continue };
            if r > 1.0 {
                continue;
            }
            let treffer = pts.iter().filter(|&&(x, y)| ((x - ux).hypot(y - uy) - r).abs() < tol).count();
            if beste.map_or(true, |(b, _)| treffer > b) {
                beste = Some((treffer, (ux, uy, r)));
            }
        }
        let (anzahl, (ux, uy, r)) = beste?;
        if anzahl < 6 {
            return None;
        }
        let q: Vec<(f64, f64)> = pts
            .iter()
            .copied()
            .filter(|&(x, y)| ((x - ux).hypot(y - uy) - r).abs() < tol)
            .collect();
        let (cx, cy, r) = ausgleich(&q)?;
        let rmse = (q.iter().map(|&(x, y)| ((x - cx).hypot(y - cy) - r).powi(2)).sum::<f64>() / q.len() as f64).sqrt();
        Some((cx + mx, cy + my, r, rmse))
    });
    Ok(ergebnis)
}

#[pymodule]
fn treeflow_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fit_circle, m)?)?;
    Ok(())
}

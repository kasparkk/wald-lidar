# treeflow_rust (optional)

The stem circle fit (`treeflow.detect.fit_circle`) in Rust. Same method as the Python version
(RANSAC through 3 random points, then least squares on the inliers), about 8-9x faster.
treeflow uses it automatically when installed and falls back to Python otherwise.

Install (needs Rust, https://rustup.rs):

```
pip install maturin
pip install ./rust_kreis          # from inside treeflow/
```

Force the Python version (e.g. to compare): set `TREEFLOW_OHNE_RUST=1`.

Check that both give the same trees: `python tests/compare_rust.py SCAN.laz`

Measured on the synthetic ground scan (`tests/make_synthetic.py`, 727,000 points, 53 stems):

| | Python | Rust |
| --- | --- | --- |
| Circle fit, 53 stems | 290 ms | 34 ms |
| Whole run (no viewer) | 1.18 s | 0.76 s |
| trees.csv | identical | identical |

Difference to the Python version: points are shifted to their centre before fitting. That changes
nothing mathematically but avoids rounding error with large UTM coordinates. The random choice of
points uses its own generator, so on noisy real data single stems can come out minimally different;
`compare_rust.py` reports how much.

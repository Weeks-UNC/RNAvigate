ShapeMapper2 profiles
=====================

ShapeMapper2 produces a useful, recognizable plot with 3 panels: Normalized
SHAPE reactivities, Mutation rates with measurement errors for modified,
untreated, and optional denatured samples, and read depths and effective read
depths for all samples.

Unlike other plot types, ShapeMapper2 profiles have only 1 creation
method, and very limited options. Below is that method, with all optional
arguments and their default values.

```python
plot1 = sample.plot_sm(
    plots=["profile", "rates", "depth"]
)
```

`sample` here is a hypothetical rnavigate.Sample object containing `"shapemap"`
data.

---

`plots`

* A list containing any of these values:
  * `"profile"`: display the Normalized SHAPE reactivities panel
  * `"rates"`: display raw mutation rates and measurement errors for all
    samples.
  * `"depth"`: display read depth and effective read depth for all
    samples.

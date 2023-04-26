Receiver operator characteristic curves
=======================================

Receiver operator characteristic (ROC) curves illustrate the ability of
per-nucleotide data to predict base-paired vs. single stranded nucleotides.
RNAvigate's ROC plot draws the ROC curve for all values and for each
nucleotide. RNAvigate calculates the area under the ROC curve (AUC) and
displays it in the legend. Higher values, approaching 1, are better, while 0.5
is not better than random guessing.

In short, a threshold value splits the per-nucleotide data. The fraction of
nucleotides above this threshold that are single-stranded is the True Positive
Rate (TPR). The fraction of base-paired nucleotides that fall above the
threshold is the False Positive Rate (FPR). The ROC curve connects the TPR and
FPR of every possible threshold value. Perfect predictors, at some point,
acheive a TPR of 1 and FPR of 0. The area under this curve would be 1.

There are two ways to quickly make ROC plots:

```python
plot1 = sample.plot_roc()
plot2 = rnavigate.plot_roc(samples=[sample])
```

`sample` here is a hypothetical rnavigate.Sample object containing data. As
written, these two lines of code do exactly the same thing: create a
five panel (all, A, U, C, G) plot using the data from `sample`. With the second
method, you can also plot multiple datasets on the same axes by passing
multiple samples to the samples argument. e.g.:

```python
plot3 = rnavigate.plot_arcs(samples=[sample, another_sample])
```

Below are all of the optional arguments that work with each of the methods
above, along with their default values. `plot4` below would produce exactly the
same result as `plot1` and `plot2`.

```python
plot4 = sample.plot_roc(
    ct="ct",
    profile="profile",
    labels=None,
    region="all",
    plot_kwargs={"figsize": None},
)
```

`ct` and `profile` accept a key of `sample.data`. These are typically the
argument names given to the `rnavigate.Sample()` method when data are loaded.
To see these keys, run: `print(sample.data.keys())`.

---

`ct`

* A sample.data key that points to a secondary structure, e.g.: `"ct"`,
  `"compct"`, `"ss"`, etc.
* Base-pairing status is taken from this value.

---

`profile`

* A sample.data key that points to per-nucleotide data.
* The default is "profile" which uses the first valid value in this list:
  ["shapemap", "dmsmap", "dancemap", "rnpmap"]
* These data are mapped to `ct` and used to calculat TPR and FPR for the ROC
  curve.

---

`labels`

* A list of strings, one for each sample.
* Defaults to using the sample name by retrieving `sample.sample`.

---

`plot_kwargs`

* A dictionary of keyword argument pairs passed to `rnavigate.ROC()`.
* These values are automatically determined by the plotting function if not
  provided.
* `"figsize"` specifies the total size of the matplotlib figure in inches.

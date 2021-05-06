
Secondary Structure options test script
=======================================
This test script is to test implemented features with 2D secondary structures.

You can also use it to view the different options and how they look.
The function calls are likely to stay the same, but the default plots
may look a little different in the future.

Currently broken:
* coloring by 3d distance.
* plotting pairs:
  * draws one line from i to j.
  * should draw 3 lines:
    * i to j+2
    * i+1 to j+1
    * i+2 to j.

Notebook set-up
---------------


```python
# This sets plots to display in-line by default
%matplotlib inline

# Import module, for high-level functions, no additional modules are needed
import plotmapper as MaP

# Creates an HTML button that hides/shows code cells
# Useful for lab notebook reports and research updates
# NOTE: this does not display well on GitHub.
MaP.create_code_button()
```


<script>
                 code_show=true;
                 function code_toggle() {
                 if (code_show) {$('div.input').hide();}
                 else {$('div.input').show();}
                 code_show = !code_show
                 }
                 $( document ).ready(code_toggle);
                 </script>
                 <form action="javascript:code_toggle()">
                 <input type="submit" value="Hide/show raw code.">
                 </form>


Initializing MaP sample
-----------------------
If you have consistently named files, (which you should), you can use a function to create a dictionary of keyword arguments (kwargs). Then, "unpack" the dictionary using the double asterisk.


```python
path = 'data/'
def kwargs(sample):
    kwargs = {}
    kwargs["profile"] = path+sample+"_rnasep_profile.txt"
    kwargs["ct"] = path+"RNaseP.ct"
    kwargs["ss"] = path+"RC_CRYSTAL_STRUCTURE.xrna"
    kwargs["rings"] = path+sample+"-rnasep.corrs"
    kwargs["pairs"] = path+sample+"-rnasep-pairmap.txt"
    kwargs["log"] = path+sample+"_shapemapper_log.txt"
    kwargs["dance_reactivities"] = path+sample+"_rnasep-reactivities.txt"
    kwargs["deletions"] = path+"example-rnasep-deletions.txt"
    kwargs["fasta"] = path+"RNaseP-noSC.fasta"
    kwargs["pdb"] = path+"3dhs_Correct.pdb"
    kwargs["pdb_name"] = "3dhs"
    return kwargs

example = MaP.Sample(sample="example2", **kwargs("example2"))
```


```python
example.make_ss(attribute="rings", metric="Statistic", Statistic=20, cdAbove=20)
```


![svg](ss_test_files/ss_test_5_0.svg)



```python
example.make_ss(attribute="deletions", Percentile=0.95)
```


![svg](ss_test_files/ss_test_6_0.svg)



```python
example.make_ss(attribute="rings", metric="Zij", Zij=10)
```


![svg](ss_test_files/ss_test_7_0.svg)



```python
example.make_ss(attribute="rings", cdBelow=30)
```


![svg](ss_test_files/ss_test_8_0.svg)



```python
example.make_ss(attribute="pairs", all_pairs=True)
```


![svg](ss_test_files/ss_test_9_0.svg)



```python
example.make_ss(attribute="pairs")
```


![svg](ss_test_files/ss_test_10_0.svg)



```python
example.make_ss(attribute="pairs", metric="distance")
```


![svg](ss_test_files/ss_test_11_0.svg)



```python
example.make_ss(attribute="deletions", metric="distance", Percentile=0.98)
```


![svg](ss_test_files/ss_test_12_0.svg)



```python

```

Secondary Structure options test script
=======================================
This test script is to test implemented features with 2D secondary structures.

You can also use it to view the different options and how they look.
The function calls are likely to stay the same, but the default plots
may look a little different in the future.

Currently broken:
* nothing that I'm aware of

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
    kwargs["sample"] = sample
    kwargs["profile"] = path+sample+"_rnasep_profile.txt"
    kwargs["ct"] = path+"RNaseP.ct"
    kwargs["ss"] = path+"RC_CRYSTAL_STRUCTURE.xrna"
    kwargs["rings"] = path+sample+"-rnasep.corrs"
    kwargs["pairs"] = path+sample+"-rnasep-pairmap.txt"
    kwargs["log"] = path+sample+"_shapemapper_log.txt"
    kwargs["dance_prefix"] = path+sample+"_rnasep"
    kwargs["deletions"] = path+"example-rnasep-deletions.txt"
    kwargs["fasta"] = path+"RNaseP-noSC.fasta"
    kwargs["pdb"] = path+"3dhs_Correct.pdb"
    kwargs["pdb_name"] = "3dhs"
    return kwargs

example = MaP.Sample(**kwargs("example2"))
```


```python
example.make_ss(ij_data="pairs")
```




    <AxesSubplot:title={'center':'example2'}>




    
![svg](ss_test_files/ss_test_5_1.svg)
    



```python
example.make_ss(ij_data="rings", metric="Statistic", Statistic=20, cdAbove=20)
```




    <AxesSubplot:title={'center':'example2'}>




    
![svg](ss_test_files/ss_test_6_1.svg)
    



```python
example.make_ss(ij_data="deletions", Percentile=0.95, colorby="sequence", ss_only=True)
```




    <AxesSubplot:title={'center':'example2'}>




    
![svg](ss_test_files/ss_test_7_1.svg)
    



```python
example.make_ss(ij_data="rings", metric="Zij", Zij=10, colorby='position')
```




    <AxesSubplot:title={'center':'example2'}>




    
![svg](ss_test_files/ss_test_8_1.svg)
    



```python
example.make_ss(ij_data="rings", cdBelow=30, profAbove=0.8)
```




    <AxesSubplot:title={'center':'example2'}>




    
![svg](ss_test_files/ss_test_9_1.svg)
    



```python
example.make_ss(ij_data="pairs", all_pairs=True)
```




    <AxesSubplot:title={'center':'example2'}>




    
![svg](ss_test_files/ss_test_10_1.svg)
    



```python
example.make_ss(ij_data="deletions", metric="Distance", Percentile=0.98)
```




    <AxesSubplot:title={'center':'example2'}>




    
![svg](ss_test_files/ss_test_11_1.svg)
    



```python

```

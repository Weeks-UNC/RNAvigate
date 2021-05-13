Longleaf/OpenOnDemand instructions
----------------------------------
From Longleaf
=============
These code blocks should be run from the bash command line.
```bash
module load anaconda/2019.10
module list
```
Make sure that there no other python environments: `module rm python` if necessary.
```bash
vi JNBTools/env.yaml
```
Change the path on the last line of the file, only change the 'psirving' directory to match your ONYEN. Save and quit.
```bash
conda env create -f JNBTools/env.yaml -n mapexplors
source activate mapexplors
python -m ipykernel install --user --name=mapexplors
```
Once you do this, exit longleaf. This will restore your previous environment when you log back in.

From OpenOnDemand
=================
1. Go to https://ondemand.rc.unc.edu/ and start a new jupyter notebook.
2. Navigate to JNBTools/JNB-examples and open 3d-test.ipynb.
3. In the upper right corner, click on "Python 3" and change it to "mapexplors" from the menu. This is the environment we created and installed for ipynb from longleaf.
4. To test that this works, restart the kernel and run all cells.
5. On the JupyterLab home page, there should be an option to start a new notebook using the "mapexplors" environment.

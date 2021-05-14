Longleaf/OpenOnDemand instructions
==================================
From Longleaf
-------------
These code blocks should be run from the bash command line in Longleaf.

---

## 1. From the home directory, clone the git repo for JNBTools

```bash
cd $HOME
git clone https://github.com/Weeks-UNC/JNBTools.git
```
It will prompt you for a GitHub username and password. For now, you must be a
member of the Weeks-UNC organization. Email psirving@email.unc.edu for access.

---

## 2. Add JNBTools to your PYTHONPATH

Add this to your .bash_profile, or .bash_{ONYEN} if using longleaf-dotfiles.
```bash
export PYTHONPATH="$PYTHONPATH:$HOME/JNBTools"
```

---

## 3. Load Anaconda and  create an environment for Jupyter Notebooks.
```bash
module load anaconda/2019.10
module list
```
Make sure that there no other python environments. Remove them if necessary.
```bash
module rm python pymol pyrosetta
```
Create a new environment from the yaml file, activate, and add to Jupyter.
```bash
conda env create -f JNBTools/env.yaml -n mapexplors
source activate mapexplors
python -m ipykernel install --user --name=mapexplors
```
Once you do this, exit longleaf. Your previous environment and default packages
will be restored when you log back in.

---

From OpenOnDemand
-----------------
OpenOnDemand is UNC's platform for using interactive programs within Longleaf.

---

1. Go to https://ondemand.rc.unc.edu/, log in, and start a jupyter notebook.
   1. Click on interactive apps, under servers, click on Jupyter Notebook.
   2. Enter the number of hours you will be working, 1 CPU, other fields blank.
   3. DONT FORGET TO SAVE BEFORE YOUR NOTEBOOK SHUTS DOWN.
2. Using the file navigation tree on the left, navigate to JNBTools/JNB-examples
   and open 3d-test.ipynb.
3. In the upper right corner, click on "Python 3" and change it to "mapexplors"
   from the menu.
4. To test that this works, restart the kernel and run all cells.
5. On the JupyterLab home page, there should be an option to start a new python
   notebook using the "mapexplors" environment.

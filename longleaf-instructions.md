Longleaf/OpenOnDemand instructions
==================================
From Longleaf
-------------
These code blocks should be run from the bash command line in Longleaf.

---

## 1. From the home directory, clone the git repo for StarMapper

```bash
cd $HOME
git clone https://github.com/Weeks-UNC/StarMapper.git
```
It will prompt you for a GitHub username and password. For now, you must be a
member of the Weeks-UNC organization. Email psirving@email.unc.edu for access.

---

## 2. Add StarMapper to your PYTHONPATH

Add this to your .bash_profile, or .bash_{ONYEN} if using longleaf-dotfiles.
```bash
export PYTHONPATH="$PYTHONPATH:$HOME/StarMapper/"
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
conda env create -f StarMapper/env.yaml -n starmapper
source activate starmapper
python -m ipykernel install --user --name=starmapper
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
2. Using the file navigation tree on the left, navigate to 
   StarMapper/JNB-examples and open starmapper-example.ipynb.
3. In the upper right corner, click on "Python 3" and change it to "starmapper"
   from the menu.
4. To test that this works, restart the kernel and run all cells.
5. On the JupyterLab home page, there should be an option to start a new python
   notebook using the "starmapper" environment.

To Manage environments
----------------------
```bash
## Create the virtual environment
conda create -n environment_name

## Activate the virtual environment
conda activate environment_name

## Make sure that ipykernel is installed
pip install --user ipykernel

## Add the new virtual environment to Jupyter
python -m ipykernel install --user --name=environment_name

## To list existing Jupyter virtual environments
jupyter kernelspec list

## To list existing conda environments
conda env list

## To remove conda environment
conda env remove -n environment_name

## To remove the environment from Jupyter
jupyter kernelspec uninstall environment_name
```

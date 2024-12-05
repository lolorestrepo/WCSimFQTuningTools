# Create a conda environment

Gonzalo installed everything on the pbs

```jsx
1) create a .condarc and add the installation path for your environment (for example)
pkgs_dirs:
 - $HOME/.conda/pkgs
envs_dirs:
  - $HOME/.conda/envs/
2)  setup common conda: source $THRONG_DIR/miniconda3/etc/profile.d/conda.sh
3) create environment conda create --name py311 python=3.11
4) activate it to check it is fine conda activate py311
5) remove it if you want conda env remove -n py311 (edited)
```

At step 3 I preferred to create my own conda environment based on hk-pilot requirements. So changed step one to point to where to find the packages…

Cloned hk-pilot in /sps/t2k/mguigue/hk_software

```jsx
conda env create -f env.yaml
conda activate hk_software
which python 
python install.py
source setup.sh

hkp --system
-> Linux_x86_64-gcc_8-python_3.10.13 # doesn't work with python 3.11 (root installation)
```
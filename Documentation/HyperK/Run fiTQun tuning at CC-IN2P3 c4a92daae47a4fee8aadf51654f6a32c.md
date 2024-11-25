# Run fiTQun tuning at CC-IN2P3

```bash
function setup_fitqun(){
eval $(ssh-agent)
ssh-add -k ~/.ssh/id_rsa_github
module add Compilers/gcc/8.3.1 git cmake 
export PATH=/pbs/throng/hyperk/Software/Julia/julia-1.10.0/bin:$PATH
source $THRONG_DIR/miniconda3/etc/profile.d/conda.sh
conda activate hk_software
cd /sps/t2k/mguigue/hk_software/hk-pilot
source setup.sh
source ../ROOT/install-Linux_x86_64-gcc_8-python_3.10.13/setup.sh
cd -
};
setup_fitqun
```

[Open questions](Open questions d3840be03ed048da864fb84f092e85cf.md)

[Create a conda environment](Create a conda environment 24c6a1df110747f08ced11a438ad2f01.md)

[Installing WCSim and tuning tools](Installing WCSim and tuning tools 9a90f3654bc04e2a85ad9a1ad3a780b3.md)

[Run tuning](Run tuning 11c34b4a53114aa99e4381c351c68a5c.md)

[Selecting the HyperK geometry in WCSim](Selecting the HyperK geometry in WCSim f2a5dad21ab1460b9c5146fa751eb438.md)
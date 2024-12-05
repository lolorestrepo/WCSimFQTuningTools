# Installing WCSim and tuning tools

On Mac:

```bash
cd ~/Work/HK/fiTQun
git clone git@github.com:gondiaz/WCSimFQTuningTools.git
# setup deployment in CLion
# use main branch
git clone git@github.com:gondiaz/WCSimFQTuner.git
# setup deployment in CLion
# for the compilation use compilation branch
```

On CC

```bash
hkp install -r WCSim

# install WCSimFQTuner
cd /sps/t2k/mguigue/hk_software/WCSimFQTuner
# installed via CLion and the new install command
```

- **I had to edit the CMakeLists a bit to adapt things**
    
    Changed to WCSim 1.12.6 (updated hyperk version)
    
    Added the install_WCSimFQTuner target
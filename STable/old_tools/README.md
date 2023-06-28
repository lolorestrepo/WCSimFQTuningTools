Instructions for procesing a single file:

1) **runScatTableLooper.C**: creates and fills scatering table histograms
```
root -b
.L runScatTableLooper.C++
runScatTableLooper("/path/to/file", 1, 1)
```

2) **mergetables.C**
3) makeScatTable.C: Divides 6D scattering table by the 4D direct isotropic table
```
root -b
.L makeScatTable.C++
makeScatTable()
```

4) **convScatTable.C**: Converts to float scattering table
(same as in 3 with makeScatTable --> convScatTable)
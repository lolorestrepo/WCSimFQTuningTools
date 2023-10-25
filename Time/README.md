

Make sure **fiTQun** is installed and you have produced the following tunning files with the proper **fiTQun** naming standard:
- **Cherenkov profile** for each particle you wish to produce the timePDF.
- The **Charge PDF** for the PMTs of the detector.
- The **Scattering Table** for the desired geometry (used to compute the **Time PDF** parameters for indirect light).
- A **Time PDF** which does not need to be the one for the detector and particle you want to produce. This is just needed because fiTQun loads a timePDF when initialized. Nevertheless, the naming convention applies, so you have to rename it to the fiTQun standard for your desired production.

The files must be placed in the **fiTQun** **const/** directory. The `gpdgTime.root` and `gpdgRange.root` files must also be present.

About the fiTQun naming standad for the tuning files:
- **Cherenkov Profiles**: CProf_{PDG MC particle code}_fit_WCSim.root
- **Charge PDF**: cPDFpar_{`fiTQun.WCSimPMTType`}.root
- **Scattering Table**: fiTQun_scattablesF_{`fiTQun.WCSimConfig`}.root
- **Time PDF**: 11_tpdfpar_{`fiTQun.WCSimConfig`}_{`fiTQun.WCSimPMTType`}.root

where the `fiTQun.{parameter}` are set in the `parameters_file` (see below) for the particular detector.


### Create 2D histogram
The first step for the **timePDF** production is to create the 2D histogram of residual time ($t^{res}$) versus true or predicted charge ($\mu^{dir}$) from direct light. This is done with the `create_2D_histogram.cc` script.

To simplify the process, a bash script *compile_create_2D_histogram.sh* is provided so you only need to modify the parameters therein. In particular, set the correct paths to your **WCSim**, **fiTQun**, and **ROOT** installations. To run the compilled executable *create_2D_histogram*, first define FITQUN_ROOT as the path to your **fiTQun** directory containing the **const/** directory with the tunning files.

    export FITQUN_ROOT=/path/to/fiTQun/sourcedir/
    ./create_2D_histogram input_file particle_name parameters_file water_refractive_index vertical_axis max-statistics binning_parameters

**TODO**: Explain the parameters
 - `statistics-limit` uses at least this many no-PC events in the computation
 - `water_refractive_index` should be the same as the fiTQun.WaterRefractiveIndex used in fiTQun

### The fits

To produce the fits, run the `fit_timepdf.py` script as:

`fit_timepdf.py [-v] [--npars_gauss n1] [--npars_pars n2]`,

where `npars_gauss` is the number of parameters (degree + 1) of the polynomia fitting the mean and sigma of the fitted gaussians; and `npars_pars` is the number of parameters used to fit the previous parameters.

$$ (\mu, \sigma) (p, \mu^{dir}) = \sum_{i=0}^{n_{gauss}} a_i(p)~\mu^{dir}$$

$$ a_i(p) = \sum_{j=0}^{n_{pars}} b_{ij}~p^j$$


### The indirect time pdf
The parameters of the indirect time pdf are computed in the `fit_timepdf.py`. Nevertheless, since the values are hardcoded in **fiTQun** you should set them manually, threfore my recommended way to check the parameters is through the `histogram.ipynb` Jupyter-Notebook.
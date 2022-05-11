# Bayesian regression model for layer-counted proxy records
 
This Repository contains all the data/code to reproduce the results for 

*Myrvoll-Nilsen, E., Riechers, K., Rypdal, M. & Boers, N. (2022). Comprehensive uncertainty estimation of the timing of Greenland warmings of the Greenland Ice core records. Climate of the Past, x(x),x-x. doi:xxxxx*

To use this code, install all files in the same directory. Use function `main` located in `main.R` as a starting point. All other functions are loaded from there.

The code should run a Bayesian regression analysis on the GICC05 time scale (*NGRIP_d18O_and_dust_5cm.xls*) from [iceandclimate.nbi.ku.dk](https://www.iceandclimate.nbi.ku.dk/data/), but should be able to analyse other layer-counted proxy records as well.

**Warning:** This code produces (by default) 30,000 simulated chronologies with individual lengths exceeding 18,000. This will require a significant amount of memory. If needed, decrease the number of samples down or free up memory before use.

## Attribution
This code is associated and written for the paper *Myrvoll-Nilsen et al. 2022* mentioned above. Feel free to use the code, but please cite the accompanying paper.

## License
The code in this repository is made available under the terms of the MIT License. For details, see LICENSE file.
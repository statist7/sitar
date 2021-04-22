## Test environments
* local MacOS install, R 4.0.5
* windows-x86_64-devel (r-devel)
* ubuntu-gcc-release (r-release)
* debian-gcc-release (r-release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), debian-gcc-release (r-release)
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Tim Cole <tim.cole@ucl.ac.uk>'
  
  Found the following (possibly) invalid URLs:
    URL: https://academic.oup.com/ajcn/article/53/4/839/4715058
      From: man/deren.Rd
      Status: 403
      Message: Forbidden
      From: man/deren.Rd
      Message: Forbidden
    URL: https://doi.org/10.1136/bmj.320.7244.1240
  
  Found the following (possibly) invalid DOIs:
    DOI: 10.1093/ije/dyq115
      From: DESCRIPTION
      Status: Forbidden
      Message: 403
      Status: 403

  I've checked the URLs and DOI and am mystified why they are flagged
  
> On ubuntu-gcc-release (r-release)
* checking for future file timestamps ... NOTE
  unable to verify current time
  

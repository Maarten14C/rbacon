## Test environments
* Fedora and Mac, R 4.4.1
* rhub and win-builder (devel, old_release and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

Using "R CMD check --as-cran rbacon_3.3.0.tar.gz" on OSX (14.6.1) (R version 4.4.1), a Note appears about HTML validation problems when the .Rd files are made into .html (e.g, "Error: <main> is not recognized!"). The .Rd files look OK though so this is probably just a problem in the html rendering. 

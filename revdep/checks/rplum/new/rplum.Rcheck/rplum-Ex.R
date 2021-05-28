pkgname <- "rplum"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rplum')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Plum")
### * Plum

flush(stderr()); flush(stdout())

### Name: Plum
### Title: Main 210Pb age-depth modelling function
### Aliases: Plum

### ** Examples




cleanEx()
nameEx("Plum.agedepth")
### * Plum.agedepth

flush(stderr()); flush(stdout())

### Name: Plum.agedepth
### Title: Plot an age-depth model
### Aliases: Plum.agedepth

### ** Examples

## Don't show: 
  Plum(run=FALSE, ask=FALSE, coredir=tempfile())
  Plum.agedepth(yr.res=50, d.res=50, d.by=10)
 
## End(Don't show)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')

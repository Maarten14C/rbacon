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
nameEx("agedepth")
### * agedepth

flush(stderr()); flush(stdout())

### Name: agedepth
### Title: Plot an age-depth model
### Aliases: agedepth

### ** Examples

  Plum(run=FALSE, ask=FALSE, coredir=tempfile())
  agedepth(age.res=50, d.res=50, d.by=1)



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

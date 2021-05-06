## Test environments
* local Fedora install, R 4.0.4
* local Mac OS, R 4.0.5
* rhub and win-builder (devel, old_release and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a resubmission, correcting 1 note and commenting on 1 false error:

* Corrected NOTE: the revdep folder has been removed

* Installation ERROR: this package installs without problems on devtools' win_oldrelease, win_release and win_devel.

The installation error reported on
https://www.r-project.org/nosvn/R.check/r-devel-windows-ix86+x86_64/rbacon-00check.html
must be a server error. No information is provided in the details page.
Similar packages relying on coda, such as spBayes, report the same error:
https://www.r-project.org/nosvn/R.check/r-devel-windows-ix86+x86_64/spBayes-00check.html


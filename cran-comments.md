## Test environments
* local Fedora install, R 4.1.3
* rhub and win-builder (devel, old_release and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

* This version solves a Warning under r-devel-linux-x86_64-debian-gcc reported by Kurt Hornik: `delete x;` has now been replaced by `delete[] x;`

* On Windows Server 2022, R-devel, 64 bit there is a note "checking for detritus in the temp directory ... NOTE
  Found the following files/directories: 'lastMiKTeXException'". This is likely a bug in [MikTex](https://github.com/r-hub/rhub/issues/503).

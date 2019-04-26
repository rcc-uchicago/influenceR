## Resubmission
I have:
* Added the author (year) <doi> references in the description. 
* Updated the CRAN comments to be more reflective of the changes.
* Added `revdep/` to `.gitignore`.
* Registered dynamic routines in `nfluenceR/src/packagename_init.c`.
* There is a note about possible spelling errors:
Possibly mis-spelled words in DESCRIPTION:
  Borgatti (13:32)
  Fujimoto (14:39)
  Valente (14:27)
These are author names and are spelled correcly.
 * I received a question from Uwe Ligges that:
  ```Flavor: r-devel-linux-x86_64-debian-gcc
Check: examples, Result: NOTE
   Examples with CPU time > 2.5 times elapsed time
            user system elapsed ratio
   bridging  1.1  0.004   0.376 2.936
```
"Are you using more than 2 cores?"  
 We use OpenMP parallelization which by default will use available cores.   
  
## Test environments
* local OS X install, R 3.5.3
* Scientific Linux 7.4, R 3.5.1
* win-builder

## R CMD check results
There were no ERRORs or WARNINGs.  

## Downstream dependencies
There are currently no downstream dependencies on this package.

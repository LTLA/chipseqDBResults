# chipseqDB results

This repository contains compiled Markdown files for the [_chipseqDB_](http://bioconductor.org/packages/devel/workflows/html/chipseqDB.html) workflow.
Each file represents the compiled version of one of the vignettes.
The aim is to allow the maintainers to easily diagnose changes to the workflow results upon updates to dependent packages.
To check the results on a fresh system:

1. Install relevant packages with `BiocManager::install("chipseqDB", dependencies=TRUE)`.
2. Run `make all`. 
This will take a bit over an hour.
3. Run `git diff` will then highlight the changes to the compiled results.
Images are not diff'd but most code chunks will report numbers for checking.

Please report any discrepancies as issues at the [`chipseqDB`](https://github.com/LTLA/chipseqDB) repository, rather than this repository.

SV Properties
=============

Computes various SV statistics from an input VCF file.

Installation
------------

`git clone --recursive https://github.com/tobiasrausch/svprops.git`

`cd svprops/`

`make all`

Usage
-----

`./src/svprops [sv.vcf|sv.bcf|sv.vcf.gz] > sv.tab`

`./src/sampleprops [sv.vcf|sv.bcf|sv.vcf.gz] > sample.tab`

Plotting
--------

`Rscript R/svprops.R sv.tab`

`Rscript R/svprops.R sample.tab`



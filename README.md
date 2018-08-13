# somaticwrapper
Detect somatic variants from tumor and normal exome data

SomaticWrapper pipeline is a fully automated and modular software package
designed for detection of somatic variants from tumor and normal exome data. 
It was developed from GenomeVIP. Multiple standard
variant callings are included in the pipeline such as varscan, strelka and
pindel. 

## Installation

See [SomaticWrapper.CPTAC3.b1](https://github.com/ding-lab/SomaticWrapper.CPTAC3.b1) for details
about installation and usage of SomaticWrapper

## Implementation

### Strelka
![Somatic Wrapper Strelka Details](docs/SomaticWrapper.CWL.Strelka.png)
### Varscan
![Somatic Wrapper Varscan Details](docs/SomaticWrapper.CWL.Varscan.png)
### Pindel
![Somatic Wrapper Pindel Details](docs/SomaticWrapper.CWL.Pindel.png)
### Merging
![Somatic Wrapper Overview](docs/SomaticWrapper.CWL.Merge.png)

## Branches

`docker` branch has work on version of SomaticWrapper which runs in dockerized container at 
MGI or DC2 (uses SomaticWrapper.Workflow for help)

`cwl` branch makes changes to make SomaticWrapper operate in CWL environment. Specific changes:
  * all arguments are passed on command line, rather than configuraiton file
  * Output directory is passed as an argument explicitly, so that directry structure is not
    dependent on run name
  * inputs and outputs are more explicitly defined
  * All steps need to have input data passed as an argument
    * in cases where a tool writes its output to the same directory as input data (`pindel_filter` and `varscan` do this)
      in order to control where output goes, we create a link to the input data in the output directory
  * One of the steps in `parse_pindel`, the `grep ChrID` step, was moved to `pindel_run`
  * All references of `genomevip_label` were removed, simplifying internal data file structure

### `run_vep`

The script `run_vep` has been replaced by a more CWL-friendly version, `annotate_vep`.  The latter
takes one VCF file as input and writes an annotated VCF (or VEP) file. As such, this "step" may be used
to process any number of files by placing it in their workflow

## Authors

* Song Cao
* Matthew Wyczalkowski

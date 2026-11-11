# ForHot: the demonstration dataset

Sixteen public soil metatranscriptomes from the ForHot geothermal warming
experiment in Iceland (Söllinger *et al.*, 2022, *Sci Adv* 8:eabm3230):
grassland soil at ambient temperature (AT) and warmed by about 6 °C (ET), at
two sites warmed for 8 years (MTW) and for more than 50 years (LTW), with four
paired gradients per site. Unstranded total-RNA libraries (NEBNext Ultra).

| File | What it is |
|---|---|
| `runs_forhot.tsv` | The 20 SRA runs of BioProject PRJNA663238 and the sample each belongs to, with ENA download links and MD5 sums. Three samples were sequenced in several runs: MTW_AT_G1 (3), MTW_AT_G3 (2), MTW_ET_G3 (2). |
| `download_forhot.sh` | Downloads the runs over HTTPS, checks every MD5 and joins the runs of each sample into `fastq/<sample>_1/_2.fastq.gz`. |
| `samplesheet_forhot.csv` | The 16 samples: treatment (`group`), `site`, `gradient` (the AT/ET pair) and library type. |
| `forhot_soil_chemistry.xlsx` | Soil variables from `001_metadata.txt` (DataverseNO doi:10.18710/HZBZQE, CC0). Sheet `soil` is read by the pipeline; `dictionary` gives each column's source, unit and use; `source` gives the file and its MD5. |
| `config_forhot.yaml` | The settings of the paper's run (only keys that differ from `config/config.yaml`). |

**Sample names.** The published files call the sites GO (= LTW) and GN
(= MTW), and number the MTW gradients 1, 2, 4 and 5; here they are G1–G4. The
mapping was checked by correlating each sample's KO profile with the
published per-sample counts (`validation/forhot_checks.R`). The source file's
"GN-4_A" is read as GN-4-A.

**Running it.**

```bash
bash examples/forhot/download_forhot.sh /path/to/forhot      # about 232 GB
# set fastq_dir in examples/forhot/config_forhot.yaml to /path/to/forhot/fastq
ulimit -n 65536
snakemake --cores 40 --keep-going --configfile config/config.yaml examples/forhot/config_forhot.yaml
```

The comparison with the published analysis and the checks behind it are in
`validation/` (see its README).

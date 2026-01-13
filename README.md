Tool to detect if genes/regions of interest within mice and human genomes are within evolutionary breakpoints.

Note, this is a test version and not yet submitted to Bioconductor.
You need to set up the environment manually until that happens.
You can get this using Conda/Mamba (make a new environment or add to existing one). <br>
After installing conda/mamba, run:
```
conda install -c bioconda bioconductor-rtracklayer bioconductor-genomicranges
conda install -c conda-forge r-data.table r-ggplot2 r-ggnewscale r-igraph r-stringr
```

To download from github, get the remotes package in R:
```
install.packages("remotes")
```

Then download the ConservationAnalysisTool. <br>
When prompted to update packages, skip or say NO. <br>
Use the authentication token given. Feel free to email me at htetm@uchicago.edu.: <br>
```
library(remotes)
install_github("Minty-Htet/ConservationAnalysisTool", auth_token = '<Token>')
```

The tools needs: <br>
species : "Mmusculus" or "Hsapiens" <br>
genome: "mm10" or "mm39" or "hg38" <br>
genes: a vector of gene symbols such as c("OASL", "P2RX7") for Hsapiens and c("Oasl1", "P2rx7") for Mmusculus <br>
peaks: a bedfile read in as a data frame. Make sure the columns are chromosome, start, end, symbol. In this case, symbol is just whatever id you wish to associate with the peak. <br>
chromosomes: a vector of chromosome of interest such as c("chr1", "chr2") or "All" <br>

Use the ConservationAnalysisTool:
```
library(ConservationAnalysisTool)
run_CAT(
species = "Hsapiens",
genome = "hg38",
genes = c("OASL", "P2RX7"),
peaks = bed_file,
chromosome = c("chr12", "chr13")
)
```


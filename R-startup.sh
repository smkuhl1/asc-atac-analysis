alias ll='ls -l'

apt-get update
apt-get install sudo -y
apt-get install cmake -y
apt-get install libxt-dev -y
apt-get install libcairo2-dev -y
conda install -c conda-forge r-curl -y
conda install -c conda-forge libuv -y
conda install -c conda-forge r-httpuv -y
conda install -c conda-forge r-ragg -y
conda install -c conda-forge r-shiny -y
conda install -c conda-forge r-pkgdown -y
conda install -c conda-forge r-miniui -y
conda install -c conda-forge r-devtools -y

conda install -c conda-forge r-dplyr r-ggplot2 r-tidyr r-pROC r-data.table r-glmnet r-Seurat=4.3.0 r-IRkernel r-ggpubr r-cairo r-proc r-rstatix r-metr r-hdf5r r-units r-sf r-terra r-tidyverse r-data.table r-devtools r-viridis r-signac conda install r-harmony  -y
conda install  -y -c bioconda bioconductor-txdb.hsapiens.ucsc.hg38.refgene bioconductor-org.hs.eg.db bioconductor-nebulosa bioconductor-mofa2 bioconductor-scran bioconductor-dittoseq bioconductor-complexheatmap
conda install -y -c bioconda bioconductor-mofa2 bioconductor-scran bioconductor-dittoseq bioconductor-complexheatmap bioconductor-deseq2 bioconductor-edger bioconductor-genomicranges bioconductor-qvalue
conda install -y -c bioconda bioconductor-genomicfeatures bioconductor-plyranges bioconductor-ggbio bioconductor-organismdbi bioconductor-biovizbase bioconductor-ensembldb
conda install -c bioconda bioconductor-chromvar -y
conda install -c bioconda bioconductor-gsva -y
conda install r-Seurat=4.3.0 


Rscript -e "install.packages('remotes', repos = 'http://cran.us.r-project.org')"
Rscript -e "remotes::install_version('Seurat', version = '4.3.0', repos = 'http://cran.us.r-project.org')"

Rscript -e 'install.packages("BiocManager")'



Rscript -e 'BiocManager::install(c("RaggedExperiment", "ComplexHeatmap",
                       "BSgenome.Hsapiens.UCSC.hg38",
                       "qvalue",
                       "DESeq2","edgeR",
                       "ggrepel","ggbio",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "org.Hs.eg.db","OrganismDbi",
                       "Signac",
                       "TxDb.Hsapiens.UCSC.hg38.refGene", "Nebulosa"
                       ))'


Rscript -e 'ip <- installed.packages()
if(!"BarMixer" %in% rownames(ip)) {
    devtools::install_github(
        "alleninstitute/BarMixer",
        upgrade = "never"
    )
}'

Rscript -e 'devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())'


install.packages("H5weaver", type = "source", repos = NULL)

required_packages <- c('ggplot2', 'dplyr', 'viridis', 'harmony', 'Matrix', 
                       'grid', 'reshape2', 'ggpubr', 'scales', 'stringr', 
                       'purrr', 'tidyr', 'tidyverse', 'janitor', 
                       'RColorBrewer', 'tibble', 'gridExtra', 'data.table')

missing_packages <- required_packages[!(required_packages %in% rownames(ip))]

# Install missing packages
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

ArchR::installExtraPackages()

quit()


R 'install.packages("lubridate")'
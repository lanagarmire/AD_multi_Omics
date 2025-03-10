# Integrated Metabolomics and Transcriptomics Analysis Identifies Molecular Subtypes within the Early and Late Mild Cognitive Impairment Stages of Alzheimer's Disease

## Description

This is the github repository for the project **Integrated Metabolomics and Transcriptomics Analysis Identifies Molecular Subtypes within the Early and Late Mild Cognitive Impairment Stages of Alzheimer's Disease** by Bowei Li, Shashank Yadav, Shu Zhou, Yuheng Du, Leyuan Qian, Zachary Karas, Yueyang Zhang and Lana Garmire et al.. It contains code and data for generating Figure 2,3,4 and label transfering in the paper. 

## Getting Started

### Dependencies
* Linux Working Environment
* [Python 3](https://www.python.org/downloads/)
* [R](https://www.R-project.org)
* [Anaconda3](https://www.anaconda.com/)
* [Jupyter](https://jupyter.org)
* Python libraries:
  * [Numpy](https://numpy.org/)
  * [pandas](https://pandas.pydata.org/docs/index.html)
  * [Scipy](https://scipy.org/)
  * [Scikit-learn](http://scikit-learn.org/)
  * [matplotlib](https://matplotlib.org/)
  * [snfpy](https://github.com/rmarkello/snfpy)
  * [seaborn](https://seaborn.pydata.org/)
  * [elpigraph](https://github.com/j-bac/elpigraph-python)
  * [trimap](https://pypi.org/project/trimap/)
  * [UMAP](https://github.com/lmcinnes/umap)
  * [pymde](https://pymde.org/)
  * [tensorflow](https://www.tensorflow.org/)
  * [Unioncom](https://github.com/caokai1073/UnionCom/)

* R libraries:
  * [limma](http://bioconductor.org/packages/release/bioc/html/limma.html)
  * [Biobase](https://bioconductor.org/packages/release/bioc/html/Biobase.html)
  * [dplyr](https://dplyr.tidyverse.org/)
  * [ggplot2](https://ggplot2.tidyverse.org/)
  * [magrittr](https://magrittr.tidyverse.org/)
  * [ggrepel](https://github.com/slowkow/ggrepel)
  * [SNFtool](https://github.com/maxconway/SNFtool)
  * [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
  * [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)
  * [org.Hs.eg.db](http://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)
  * [ggnewscale](https://cran.r-project.org/web/packages/ggnewscale/index.html)
  * [scales](https://scales.r-lib.org/)
  * [ADNIMERGE](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://adni.loni.usc.edu/wp-content/uploads/2012/08/instruction-ADNIMERGE-packages.pdf)


### Installing

Installing the R kernel on the jupyter
```R
install.packages('IRkernel')
IRkernel::installspec()  # to register the kernel in the current R installation
```
Use the [Bioconductor](https://www.bioconductor.org/install/) to install R packages.
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(*PackageName* = "*Version*")
```

For [ADNIMERGE](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://adni.loni.usc.edu/wp-content/uploads/2012/08/instruction-ADNIMERGE-packages.pdf) please send a request to [ADNI](https://adni.loni.usc.edu/) for a data usage permission



Use the package manager [pip](https://pip.pypa.io/en/stable/) to install python packages.
```bash
pip install Numpy
```

## Repository directories & files

The directories are as follows:
+ [`figs`](figs) contains subfigures of Figure 2,3,4 showing in the paper.
+ [`extdata`](extdata) contains the external data that required for preprocessing and further figure generation
+ [`clintra`](clintra) contains the python script of [ClinTrajan](https://github.com/auranic/ClinTrajan) for [`4_clintrajan.ipynb`](4_clintrajan.ipynb)

The other files are as follows.
+ [`2_limma_EMCI_META.ipynb`](2_limma_EMCI_META.ipynb) contains the differential analysis based on the metabolomics on the EMCI subtypes.
+ [`2_limma_EMCI_MRNA.ipynb`](2_limma_EMCI_MRNA.ipynb) contains the differential analysis based on the MRNA on the EMCI subtypes.
+ [`2_limma_LMCI_META.ipynb`](2_limma_LMCI_META.ipynb) contains the differential analysis based on the metabolomics on the LMCI subtypes.
+ [`2_limma_EMCI_MRNA.ipynb`](2_limma_LMCI_MRNA.ipynb) contains the differential analysis based on the MRNA on the LMCI subtypes.
+ [`2_SNF_MCI.ipynb`](2_SNF_MCI.ipynb) contains the Similarity Network Fusion (SNF) clustering on the MCI data to get subtypes.
+ [`2_snf_python_EMCI.ipynb`](2_SNF_EMCI.ipynb) contains the script to determine the best number of SNF clusters on EMCI data.
+ [`2_snf_python_LMCI.ipynb`](2_SNF_LMCI.ipynb) contains the script to determine the best number of SNF clusters on LMCI data.
+ [`3_GSEA.ipynb`](3_GSEA.ipynb)contains the gene set enrichment analysis based on differential analysis results
+ [`4_clintrajan.ipynb`](4_clintrajan.ipynb)contains the trajecotry analysis of different subgroups based on their clinical variables.
+ [`4_prepare_time_to_diagnosis.ipynb`](4_prepare_time_to_diagnosis.ipynb) is for preparing the time to diagnoisis analysis plot
+ [`4_time_to_diagnosis.ipynb`](4_time_to_diagnosis.ipynb)contains the Time-to-diagnosis analysis for the percentage of dementia in every six-month check up, among the EMCI and LMCI subtypes.
+ [`5_Preprocess_EFIGA.R`](5_Preprocess_EFIGA.R) contains EFIGA data preprocessing for later label transfering.
+ [`5_Unioncom_label_transfer.py`](5_Unioncom_label_transfer.py) contains label transfering using ADNI labels and EFIGA metabolites data using Unioncom for MCI subtypes.
+ [`6_Unioncom_refined_label_transfer.py`](5_Unioncom_refined__label_transfer.py) contains label transfering using ADNI labels and EFIGA metabolites data using Unioncom for refined MCI subtypes.

### Local execution
+ for `.ipynb` files: Using jupyter lab to execute the codes


## Data Availability
The external data used in this study (located in the `extdata` folder) excluded the EFIGA dataset. Access to the EFIGA data is available upon reasonable request. EFIGA Data Request can be made at https://www.neurology.columbia.edu/research/research-centers-and-programs/alzheimers-disease-research-center-adrc/investigators/investigator-resources


## Authors

### Contributors names and contact info

+ [@Lana Garmire](https://github.com/lanagarmire)
+ [@Shashank Yadav](https://github.com/xinformatics)
+ [@Shu Zhou](https://github.com/Sukumaru)
+ [@Bowei Li](https://github.com/GHBLA-NI)
+ [@Yuheng Du](https://github.com/yhdu36)

### Current Maintainer
* Bowei Li - https://github.com/GHBLA-NI

## License

This project is licensed under the `GNU General Public License v3.0` License - see the LICENSE.md file for details

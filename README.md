# Integrated Metabolomics and Transcriptomics Analysis Identifies Molecular Subtypes within the Early and Late Mild Cognitive Impairment Stages of Alzheimer's Disease

## Description

This is the github repository for the project **Integrated Metabolomics and Transcriptomics Analysis Identifies Molecular Subtypes within the Early and Late Mild Cognitive Impairment Stages of Alzheimer's Disease** by Shashank Yadav, Shu Zhou, Leyuan Qian, Zachary Karas, Yueyang Zhang and Lana Garmire et al.. It contains code and data for generating Figure 2,3,4 in the paper. 

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
  * [coxnnet](http://garmiregroup.org/cox-nnet/docs/)
  * [Theano](https://github.com/Theano/Theano)
  * [tqdm](https://github.com/tqdm/tqdm)
  * [pickle](https://docs.python.org/3/library/pickle.html)
  * [seaborn](https://seaborn.pydata.org/)
  * [holoviews](https://holoviews.org/)
  * [plotly](https://plotly.com/)
  * [lifelines](https://lifelines.readthedocs.io/en/latest/)
  * [UnionCom](https://github.com/caokai1073/UnionCom)
* R libraries:
  * [NMF](https://cran.r-project.org/web/packages/NMF/index.html)
  * [circlize](https://github.com/jokergoo/circlize)
  * [BBmisc](https://cran.rstudio.com/web/packages/BBmisc/index.html)


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

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install python packages.
```bash
pip install Numpy
```

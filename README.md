# Antibody affinity birth through somatic hypermutation

**VARIANCE** (Variation Analysis and Rational Investigation in Antibody sequeNCe Evolution) pipeline is used to analyze Sanger sequencing datasets of **Antibody affinity birth through somatic hypermutation** publication.

## Introduction:
This pipeline is divided to 8 sections. At the beginning of each section there is a comment which indicates which figures of the publication are generated based on that section.

VARIANCE ideally should be used on Linux or MacOS platforms. Running this pipeline on other operating systems is not guaranteed.

The pipeline script is located in "./script" folder. Aligned Sanger sequencing input files for these analyses are uploaded in "./data/input" folder. After a successful run of the pipeline, the analysis results will be saved in "./data/output" folder.

## Installation:
To perfectly reproduce the results of this publication, install Miniconda package manager and environment management system using its documentation which can be found [here](https://docs.conda.io/projects/miniconda/en/latest/)

Next, create a Conda environment using the following command:

`conda create -n VARIANCE python=3.11.3 pandas=2.0 numpy=1.24 matplotlib-base=3.7.1`

Next, activate the Conda environment using the following command:

`conda activate VARIANCE`

Next, install logomaker python package dependency to generate sequence logos using the following command:

`pip install logomaker`

Download VARIANCE Github repository using the following command:

`git clone https://github.com/Shahab-Sa/VARIANCE-affinity-birth`

## Usage

`cd VARIANCE-affinity-birth/scripts`

`python run.py`

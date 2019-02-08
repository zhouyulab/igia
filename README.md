# Integrative Gene Isoform Assembler (IGIA)

## Contents

- [Overview](#overview)
- [Repo contents](#repo-contents)
- [System requirements](#system-requirements)
- [Installation guide](#installation-guide)
- [Demo](#demo)
- [Instructions for use](#instuctions-for-use)
- [License](./LICENSE.txt)
- [Issues](https://github.com/zhouyulab/igia/issues)

# Overview
Currently there are multiple high-throughput sequencing techniques for transcriptome profiling. The next generation sequencing (NGS) based RNA-seq which generates millions of short reads, is often used for gene expression profiling, but it doesn't have the capability to identify accurate full-length transcripts, not mentioning potential amplification biases introduced during library construction. Pacbio sequencing offers long reads, with average read lengths over 10 kb but is hindered by lower throughput, higher error rate (11%-15%) and larger cost.  We devised a computational pipeline named Integrative Gene Isoform Assembler (IGIA) to reconstruct accurate gene structures from improved Pacbio long reads with ssRNA-seq correction, and TSS/TES boundary information.

# Repo contents

- [igia](./igia): `Python` package code.
- [docs](./docs): `IGIA` package documentation.
- [tests](./tests): `Python` unit tests written using the `unittest` package.

# System requirements

## Hardware requirements

The IGIA package can run on a standard computer or server cluster. For single-process mode, we recommend a computer with more than 32 GB RAM. For the MPI mode, we recommend preparing 16 GB of memory for each core.

## Software Requirements

### OS requirements

The package has been tested on the following *Linux* operating systems.

- Linux: Ubuntu 16.04
- Linux: Red Hat 4.8.3
- MacOS: 

# Installation guide

Before setting up the `IGIA` package, users should have `python` version 3.5.2 or higher, and several packages installed from [PyPi](https://pypi.org/). Here, we recommend to use [Conda](https://conda.io/projects/conda/en/latest/) to install and use IGIA.

## Download IGIA

  ```bash
  git clone https://github.com/zhouyulab/igia.git path/to/igia
  ```
replacing the 'path/to/' with the path to your local copy of the repo.

## Prepare the virtual environment

IGIA is implemented in Python, and depends on several packages. With the installation and activation of virtual environment (with Conda) as shown below, you can ensure that the tools run properly.

### Step 1: Download Miniconda3

  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  ```

### Step 2: Create a new Python3 environment

  ```bash
  conda create -n igia python=3.5
  ```

## Install IGIA

### Software dependencies

- [numpy](https://www.numpy.org/) (v1.11.3)
- [scipy](https://www.scipy.org/) (v0.18.1)
- [pysam](https://pysam.readthedocs.io/) (v0.10)
- [networkx](https://networkx.github.io/) (v1.11)
- [deepTools](https://deeptools.readthedocs.io/) (v2.5.1)
- [bx-python](https://pypi.org/project/bx-python/) (v0.7.3)
- [pybedtools](https://daler.github.io/pybedtools/) (v0.7.10)
- [pyBigWig](https://pypi.org/project/pyBigWig/) (v0.3.4)
- [mpi4py](https://pypi.org/project/mpi4py/) (v2.0.0) # required only for MPI mode

You can use following way to install the dependencies first, then install IGIA package. The typical install time on a "normal" desktop computer is several minutes, depending on the network speed. 

  ```bash
  cd path/to/igia
  source activate igia
  pip install -r requirements.txt
  python setup.py install
  ```

# Demo

If you have successfully installed IGIA, you can use the following command to run IGIA on test data.

  ```bash
  cd /path/to/igia/tests
  bash ./demo.sh
  ```

The expected run time for demo on a "normal" desktop computer is about 4 minutes, and the results from IGIA will be generated in `/path/to/igia/tests/igia_res`.

The expected output include several `iso*.bed12` files, a set of assembled transcripts in BED12 format, and 4 `*.bed6` files for different genomic elements identified.

For IGIA assembled transcripts, isoF and isoA are the most reliable annotations. For details, please refer IGIA manuscript.

# Instructions for use

To run IGIA with single-threaded mode, you can execute:

  ```bash
  source activate igia
  igia --tgs tgs1.bam --tss tss.csv --tes tes.csv --ngs ngs1.bam ngs2.bam -o igia_res
  ```

  See `/path/to/igia/tests/example.sh` for usage.

To run IGIA with MPI mode in a cluster, you must first ensure that Openmpi/Mpich is installed and already configured in the cluster. Then you can execute:

  ```bash
  source activate igia
  mpirun -genv I_MPI_DEVICE ssm -n 8 igiampi \
    --tgs tgs1.bam --tss tss.csv --tes tes.csv --ngs ngs1.bam ngs2.bam -o igia_res
  ```



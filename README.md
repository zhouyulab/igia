# IGIA

Integrative Gene Isoform Assembler (IGIA), a computational pipeline to reconstruct accurate gene structures from improved Pacbio long reads with ssRNA-seq correction, and TSS/TES boundary information. 


## Download IGIA

  ```bash
  git clone https://github.com/zhouyulab/igia.git path/to/igia
  ```

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

  ```bash
  cd path/to/igia
  source activate igia
  pip install -r requirements.txt
  python setup.py install
  ```

With the setup of the environment, the following packages are installed automatically.

- [pysam](https://pysam.readthedocs.io/) (v0.10)
- [numpy](https://www.numpy.org/) (v1.11.3)
- [scipy](https://www.scipy.org/) (v0.18.1)
- [networkx](https://networkx.github.io/) (v1.11)
- [deepTools](https://deeptools.readthedocs.io/) (v2.5.1)
- [bx-python](https://pypi.org/project/bx-python/) (v0.7.3)
- [pybedtools](https://daler.github.io/pybedtools/) (v0.7.10)
- [pyBigWig](https://pypi.org/project/pyBigWig/) (v0.3.4)
- [mpi4py](https://pypi.org/project/mpi4py/) (v2.0.0)


## Execute IGIA

To run IGIA with single-threaded mode, you can execute:
  ```bash
  source activate igia
  igia --tgs tgs1.bam --tss tss.csv --tes tes.csv --ngs ngs1.bam ngs2.bam -o igia_res
  ```

To run IGIA with MPI mode in a cluster, you must first ensure that Openmpi/Mpich is installed and already configured in the cluster. Then you can execute:
  ```bash
  source activate igia
  mpirun -genv I_MPI_DEVICE ssm -n 8 igiampi \
    --tgs tgs1.bam --tss tss.csv --tes tes.csv --ngs ngs1.bam ngs2.bam -o igia_res
  ```

## Execute IGIA with test data in the package

If you have successfully installed IGIA, you can use the following command to run IGIA on test data.
  ```bash
  cd /path/to/igia/test
  bash ./example.sh
  ```

This demo run on example data will execute for a few minutes, and the results from IGIA will be generated in /path/to/igia/test/igia_res.

The results contain iso*.bed12 files, a set of assembled transcripts in BED12 format, and 4 *.bed6 files for different genomic elements identified.

For IGIA assembled transcripts, isoF and isoA are the most reliable annotations. For details, please refer IGIA manuscript.



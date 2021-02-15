![logo](QBRC.jpg)
# BepiTBR
Improved B cell epitope prediction using T cell-based prediction
## Introduction
The ability to predict B cell epitopes from antigen sequences is critical for biomedical research and many clinical applications. However, despite substantial efforts over the past 20 years, the performance of even the best B cell epitope prediction software is still modest. Based on the idea of T-B reciprocity, BepiTBR is a B cell epitope prediction model that demonstrates improved performance by incorporating CD4+ T cell epitope prediction. 

Researchers interested in more information about BepiTBR and other bioinformatics tools can visit Dr. Tao Wang's [lab website](https://qbrc.swmed.edu/labs/wanglab/index.php). 
## Getting started
### System requirements
BepiTBR requires a linux x86-64 operating system with basic utilities (tested on RHEL 6, kernel 3.10.0-693 and Ubuntu 18.04, 20.04).
### Installation
BepiTBR is written in python, raku, and R and can be downloaded from [github](https://github.com/zzhu33/BepiTBR/blob/main/BepiTBR.zip). Note that some dependecies need to be manually installed.
### Dependencies
Raku v6.d or later<br/>
python 3.6.4+<br/>
R 3.6+, glmnet package<br/>
conda 4.4.10+<br/>
java 1.6+<br/>
PERL 5.0+<br/>
gcc 5.4.0+<br/>
g++ 4.8.5+

#### BepiPred 1.0
download site: [http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?bepipred](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?bepipred)<br/>
follow the instructions in `bepipred-1.0.readme` to complete installation.
#### BepiPred 2.0
Follow the instructions in `bp2_env_install_instructions.txt` under `install` to install the bp2 conda environment.
#### LBEEP 1.0
download site: [https://github.com/brsaran/LBEEP](https://github.com/brsaran/LBEEP)<br/>
follow the instructions in `README.md` to install LBEEP 1.0.
#### NetMHCIIpan 3.2 
download site: [http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan)<br/>
follow the instructions in `netMHCIIpan-3.2.readme` to complete installation.
#### MixMHC2pred 
download site: [https://github.com/GfellerLab/MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred)<br/>
follow the instructions in `README.md` to install.
## Tutorial
This tutorial will guide the user in running BepiTBR in several different modes. 
### Epitope mode
Analyzes specific epitopes within proteins<br/>
Requires epitope sequences and their corresponding full protein sequences<br/>
Example command using provided example files:<br/>
(replace paths with appropriate paths according to user's particular installation)<br/>
```
raku BepiTBR.raku \
--motif0_file=examples/test_data_BepiTBR/Ind-positive.txt \
--full0_file=examples/test_data_BepiTBR/peptide_with_full_length.txt \
--bepipred2=bp2/bin/activate \
--bepipred1=/home/exampleUser/bp1/bepipred-1.0/bepipred \
--LBEEP=/home/exampleUser/LBEEP/ \
--MixMHC2pred=/home/exampleUser/MixMHC2pred/MixMHC2pred_unix \
--netMHCIIpan=/home/exampleUser/netMHCIIpan-3.2/netMHCIIpan \
--dir=example/test_output_BepiTBR \
--thread=20
```
`--motif0_file`: candidate epitopes file<br/>
`--full0_file`: full protein sequences file<br/>
`--bepipred2`: the conda environment activation file for BepiPred 2.0 (`<bp2_env_directory>/bin/activate`)<br/>
`--bepipred1`: path to the BepiPred 1.0 executable<br/>
`--LBEEP`:  path to the LBEEP directory<br/>
`--MixMHC2pred`: path to the MixMHC2pred executable<br/>
`--netMHCIIpan`: path to the NetMHCIIpan 3.2 executable<br/>
`--dir`: output directory<br/>
`--thread`: number of CPU thread to use; system dependent<br/>

### Full protein mode
Identify potential B cell epitopes from a protein sequence<br/>
Requires the protein sequence to be input directly<br/>
Example command:<br/>
(replace paths with appropriate paths according to user's particular installation)<br/>
```
raku BepiTBR_full.raku \
--full0=$(cat examples/test_data_BepiTBR_full/example_full.txt) \
--length=15 \
--bepipred2=bp2/bin/activate \
--bepipred1=/home/exampleUser/bp1/bepipred-1.0/bepipred \
--LBEEP=/home/exampleUser/LBEEP/ \
--MixMHC2pred=/home/exampleUser/MixMHC2pred/MixMHC2pred_unix \
--netMHCIIpan=NA \
--dir=example/test_output_BepiTBR \
--thread=20
```
`--full0`: protein sequence<br/>
`--length`: length of the B cell epitopes to scan in a moving window. Recommended length: 15<br/>
`--bepipred2`: the conda environment activation file for BepiPred 2.0 (`<bp2_env_directory>/bin/activate`)<br/>
`--bepipred1`: path to the BepiPred 1.0 executable<br/>
`--LBEEP`:  path to the LBEEP directory<br/>
`--MixMHC2pred`: path to the MixMHC2pred executable<br/>
`--netMHCIIpan`: recommended to set to NA to greatly improve run speed<br/>
`--dir`: output directory<br/>
`--thread`: number of CPU thread to use; system dependent<br/>

### Fasta mode
Identify potential B cell epitopes from a file containing multiple protein sequences<br/>
Requires protein sequence file in fasta format<br/>
Example command:<br/>
(replace paths with appropriate paths according to user's particular installation)<br/>
```
raku BepiTBR_full.raku \
--fasta0=examples/test_data_BepiTBR_fasta/peptide_with_full_length.txt \
--length=15 \
--bepipred2=bp2/bin/activate \
--bepipred1=/home/exampleUser/bp1/bepipred-1.0/bepipred \
--LBEEP=/home/exampleUser/LBEEP/ \
--MixMHC2pred=/home/exampleUser/MixMHC2pred/MixMHC2pred_unix \
--netMHCIIpan=NA \
--dir=example/test_output_BepiTBR \
--thread=20
```
`--fasta0`: protein sequence fasta<br/>
`--length`: length of the B cell epitopes to scan in a moving window. Recommended length: 15<br/>
`--bepipred2`: the conda environment activation file for BepiPred 2.0 (`<bp2_env_directory>/bin/activate`)<br/>
`--bepipred1`: path to the BepiPred 1.0 executable<br/>
`--LBEEP`:  path to the LBEEP directory<br/>
`--MixMHC2pred`: path to the MixMHC2pred executable<br/>
`--netMHCIIpan`: recommended to set to NA to greatly improve run speed<br/>
`--dir`: output directory; results for each antigen sequence will be placed in a sub-directory <br/>
`--thread`: number of CPU thread to use; system dependent<br/>


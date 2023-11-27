![logo](QBRC.jpg)
# BepiTBR
Leveraging T-B reciprocity to enhance B cell epitope prediction
## Introduction
The ability to predict B cell epitopes from antigen sequences is critical for biomedical research and many clinical applications. However, despite substantial efforts over the past 20 years, the performance of even the best B cell epitope prediction software is still modest. Based on the idea of T-B reciprocity, BepiTBR is a B cell epitope prediction model that demonstrates improved performance by incorporating prediction of nearby CD4+ T cell epitopes close to the B cell epitopes. 

Researchers interested in more information about BepiTBR and other bioinformatics tools can visit Dr. Tao Wang's [lab website](https://qbrc.swmed.edu/labs/wanglab/index.php). 
## Online version
A free [online version](https://dbai.biohpc.swmed.edu/bepitbr/) of BepiTBR is provided as a part of the Database for Actionable Immunology ([dbAI](https://dbai.biohpc.swmed.edu/)) to facilitate its use by researchers. However, running a local version of BepiTBR is recommended for users with large numbers of samples.
## Installation guide
### System requirements
BepiTBR requires a linux x86-64 operating system with basic utilities (tested on RHEL 6, kernel 3.10.0-693 and Ubuntu 18.04, 20.04).
### Installation
BepiTBR is written in python, raku, and R. It and can be downloaded from [github](https://github.com/zzhu33/BepiTBR/releases). Simply extract and place the BepiTBR directory in a desired location. Note that some dependecies need to be manually installed. Although actual compile/install times are low, obtaining and installing the dependencies may take approximately 1 hour. ~5 GB of free disk space is also required during installation. 
### Dependencies
Raku v6.d or later<br/>
python 3.6.4+<br/>
R 3.6+, glmnet package<br/>
conda 4.4.10+<br/>
java 1.6+<br/>
PERL 5.0+<br/>
gcc 5.4.0+<br/>
g++ 4.8.5+<br/>
tar

#### BepiPred 1.0
[download link](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=bepipred&version=1.0&packageversion=1.0c&platform=Linux)<br/>
follow the instructions in `bepipred-1.0.readme` to complete installation.
#### BepiPred 2.0
follow the instructions in `bp2_env_install_instructions.txt` under the `install` directory to install the bp2 conda environment.
#### LBEEP 1.0
download site: [https://github.com/brsaran/LBEEP](https://github.com/brsaran/LBEEP)<br/>
follow the instructions in `README.md` to install LBEEP 1.0.
#### NetMHCIIpan 3.2 (optional)
download site: [http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?netMHCIIpan)<br/>
follow the instructions in `netMHCIIpan-3.2.readme` to complete the installation.
#### MixMHC2pred 
download site: [https://github.com/GfellerLab/MixMHC2pred](https://github.com/GfellerLab/MixMHC2pred)<br/>
follow the instructions in `README.md` to install.
## Tutorial
This tutorial will guide the user in running BepiTBR in several different modes. The generated results can be compared to the provided results in `examples/example_output`
### Epitope mode
Analyzes specific epitopes within proteins<br/>
Requires epitope sequences and their corresponding full antigen protein sequences<br/>
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
--dir=/home/exampleUser/BepiTBR/examples/test_output_BepiTBR \
--thread=4
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
#### Output
The main output is the `predictions.txt` file in the output directory. Each line corresponds to an epitope and consists of tab-seperated data with the format:
epitope name, base bepipred1.0,	base bepipred2.0,	base LBEEP, enhanced bepipred1.0,	enhanced bepipred2.0,	enhanced LBEEP,	ensemble,	epitope sequence. For example:
```
base_bepipred1.0	base_bepipred2.0	base_LBEEP	enhanced_bepipred1.0	enhanced_bepipred2.0	enhanced_LBEEP	ensemble	epitope
Negative_IEDB_ID_946463	1.927	0.6338	0.67	-0.5321	-0.3189	-0.4096	0.8883	KKLIPNPNKIRKPPK
Negative_IEDB_ID_947078	1.013	0.6201	0.51	-0.5706	-0.4672	-0.6792	0.5832	VTRLRYRSVREVWQS
...
```
Here, `base` indicates original prediction scores using existing B cell epitope prediction software, `enhanced` indicates predictions made by considering both the base models and T cell epitope predictions, for each of the three base models, and `ensemble` is the final BepiTBR prediction score that is the aggregate of the three enhanced models.<br/>
The B cell and T cell (MHC II) predictions used to calculate the final scores are removed by defaul in order to save space; they can be kept using `--keep=True`.<br/>
The example code should take ~4 minutes to complete. Performance can be increased by increasing `--thread` if hardware resources allows. Peak memory usage is approximately 3 GB/thread in this mode.
### Full protein mode
Identify potential B cell epitopes from a full antigen protein sequence<br/>
Requires the antigen protein sequence to be input directly<br/>
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
--dir=/home/exampleUser/BepiTBR/example/test_output_BepiTBR_full \
--thread=8
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
#### Output
The main output is the `predictions.txt` file in the output directory. Each line corresponds to the start location of a potential B cell epitope on the protein and consists of tab-seperated data with the format:
epitope name, base bepipred1.0,	base bepipred2.0,	base LBEEP, enhanced bepipred1.0,	enhanced bepipred2.0,	enhanced LBEEP,	ensemble,	epitope sequence. For example:
```
base_bepipred1.0	base_bepipred2.0	base_LBEEP	enhanced_bepipred1.0	enhanced_bepipred2.0	enhanced_LBEEP	ensemble	epitope
job_pos=0	1.855	0.6331	0.47	-0.4154	-0.5009	-0.8526	0.5792	MTENSTSAPAAKPKR
job_pos=1	1.855	0.6418	0.48	-0.4242	-0.4819	-0.8472	0.5914	TENSTSAPAAKPKRA
job_pos=2	1.855	0.6418	0.48	-0.433	-0.4867	-0.8582	0.5743	ENSTSAPAAKPKRAK
...
```
same as epitope mode, `base` indicates original prediction scores using existing B cell epitope prediction software, `enhanced` indicates predictions made by considering both the base models and T cell epitope predictions, for each of the three base models, and `ensemble` is the final BepiTBR prediction score that is the aggregate of the three enhanced models. <br/>
Same as epitope mode, the B cell and T cell (MHC II) predictions used to calculate the final scores are compressed in order to save space; they are kept as `Bepi.tar.gz` and `Tepi.tar.gz`, respectively.<br/>
The example code should take 1-2 minutes to complete. Unlike epitope mode, memory use does not scale with `--thread`, thus using the maximum number of CPU cores available is recommended.

### Fasta mode
Identify potential B cell epitopes from a .fasta file containing multiple antigen protein sequences<br/>
Requires protein sequence file in fasta format<br/>
Example command:<br/>
(replace paths with appropriate paths according to user's particular installation)<br/>
```
raku BepiTBR_fasta.raku \
--fasta0=examples/test_data_BepiTBR_fasta/peptide_with_full_length.txt \
--length=15 \
--bepipred2=bp2/bin/activate \
--bepipred1=/home/exampleUser/bp1/bepipred-1.0/bepipred \
--LBEEP=/home/exampleUser/LBEEP/ \
--MixMHC2pred=/home/exampleUser/MixMHC2pred/MixMHC2pred_unix \
--netMHCIIpan=NA \
--dir=/home/exampleUser/BepiTBR/examples/test_output_BepiTBR_fasta \
--thread=4
```
`--fasta0`: protein sequence file in fasta format<br/>
`--length`: length of the B cell epitopes to scan in a moving window. Recommended length: 15<br/>
`--bepipred2`: the conda environment activation file for BepiPred 2.0 (`<bp2_env_directory>/bin/activate`)<br/>
`--bepipred1`: path to the BepiPred 1.0 executable<br/>
`--LBEEP`:  path to the LBEEP directory<br/>
`--MixMHC2pred`: path to the MixMHC2pred executable<br/>
`--netMHCIIpan`: recommended to set to NA to greatly improve run speed<br/>
`--dir`: output directory; results for each antigen sequence will be placed in a sub-directory <br/>
`--thread`: number of CPU thread to use; system dependent<br/>
#### Output
For each protein in the fasta file labelled with a header line, a sub-directory for that protein will be created in the output directory. For example, for the input in `examples/test_data_BepiTBR_fasta`, two sub-directories, `Positive_IEDB_ID_114414` and `Positive_IEDB_ID_123439` will be created in the output directory. In each sub-directory, the main output is the `predictions.txt` file. Same as for full protein mode, each line corresponds to the start location of a potential epitope on the protein and consists of tab-seperated data with the format:
epitope name, base bepipred1.0,	base bepipred2.0,	base LBEEP, enhanced bepipred1.0,	enhanced bepipred2.0,	enhanced LBEEP,	ensemble,	epitope sequence. For `Positive_IEDB_ID_114414`, the example output is as follows:
```
base_bepipred1.0	base_bepipred2.0	base_LBEEP	enhanced_bepipred1.0	enhanced_bepipred2.0	enhanced_LBEEP	ensemble	epitope
job_pos=0	0.656	0.3888	0.43	-0.8591	-1.2612	-1.2655	-0.6383	MASSSSVLLVVALFA
job_pos=1	0.656	0.3888	0.32	-0.844	-1.2423	-1.583	-0.7871	ASSSSVLLVVALFAV
job_pos=2	0.288	0.3888	0.29	-0.9434	-1.235	-1.6714	-0.9076	SSSSVLLVVALFAVF
...
```
and for `Positive_IEDB_ID_123439`:
```
base_bepipred1.0	base_bepipred2.0	base_LBEEP	enhanced_bepipred1.0	enhanced_bepipred2.0	enhanced_LBEEP	ensemble	epitope
job_pos=0	-1.393	0.3664	0.28	-1.891	-1.1241	-1.5599	-1.4901	MRLLQCVLLCVSLSL
job_pos=1	-0.917	0.3923	0.3	-1.6216	-1.0555	-1.4938	-1.1855	RLLQCVLLCVSLSLV
job_pos=2	-0.719	0.4382	0.32	-1.5136	-0.9351	-1.4293	-0.9637	LLQCVLLCVSLSLVL
...
```
same as the other two modes, `base` indicates original prediction scores using existing B cell epitope prediction software, `enhanced` indicates predictions made by considering both the base models and T cell epitope predictions, for each of the three base models, and `ensemble` is the final BepiTBR prediction score that is the aggregate of the three enhanced models.<br/>
Also the same as the other two modes, the B cell and T cell (MHC II) predictions used to calculate the final scores are compressed in order to save space; they are kept as `Bepi.tar.gz` and `Tepi.tar.gz`, respectively.<br/>
A combined output is also provided in the main output directory.<br/>
The example code should take ~7 minutes to complete. Using the maximum number of CPU cores available / 4 is recommended.

## Data 
### Training/validation data
The training/validation data used for BepiTBR were downloaded from [IEDB](https://www.iedb.org/) and then filtered according to a number of criteria. It is provided on the BepiTBR github [page](https://github.com/zzhu33/BepiTBR/blob/main/data.zip).
### SARS-CoV-2 epitope data
The epitope predictions for select SARS-CoV-2 sequences are available on the BepiTBR github page ([reference sequence predictions](https://github.com/zzhu33/BepiTBR/blob/main/SARS-CoV-2_B-epitope_predictions.txt), [B epitopes](https://github.com/zzhu33/BepiTBR/blob/main/unique_variants_bepiTBR_epitopes.tar.xz), [T epitopes](https://github.com/zzhu33/BepiTBR/blob/main/unique_variants_T_epitopes.tar.xz)). For the reference sequence predictions file, "label" refers to the protein/ORF name, "epitope_start_position" the protein position of the first amino acid in potential epitopes, and the other column being are the same as those in bepiTBR outputs. The suggested cutoff for epitopes are 0 for general epitopes, and 0.5 for highly likely epitopes (referring to the "ensemble" column). The original viral sequences were downladed from [NGDC](https://bigd.big.ac.cn/ncov) and [GISAID](https://www.gisaid.org/). The original sequences are not provided here as per GISAID's terms of use. Due to large file sizes, these two files were uploaded using [git lfs](https://git-lfs.github.com/). To download the files, use `git clone` or `git pull` once the git lfs extension is installed.
## Citing BepiTBR
If you use BepiTBR in your research, please cite our [paper](https://www.cell.com/iscience/fulltext/S2589-0042(22)00034-7).<br/>
Zhu J, Gouru A, Wu F, Berzofsky JA. Xie Y, and Wang T. (2022). BepiTBR: T-B reciprocity enhances B cell epitope prediction. iScience. 25. 103764. 10.1016/j.isci.2022.103764. 

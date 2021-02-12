![logo](QBRC.jpg)
# BepiTBR
Improved B cell epitope prediction using T cell-based prediction
## Introduction
The ability to predict B cell epitopes from antigen sequences is critical for biomedical research and many clinical applications. However, despite substantial efforts over the past 20 years, the performance of even the best B cell epitope prediction software is still modest. Based on the idea of T-B reciprocity, BepiTBR is a proof of concept B cell epitope prediction model that demonstrates improved performance by incorporating CD4+ T cell epitope prediction. 

Researchers interested in more information about BepiTBR and other bioinformatics tools can visit Dr. Tao Wang's [lab website](https://qbrc.swmed.edu/labs/wanglab/index.php). 
## Getting started
### System requirements
BepiTBR requires a linux x86-64 operating system with basic utilities (tested on RHEL 6, kernel 3.10.0-693 and Ubuntu 18.04, 20.04).
### Installation
BepiTBR is written in python and raku and can be downloaded from [github](https://github.com/zzhu33/BepiTBR/blob/main/BepiTBR.zip). Note that dependecies need to be manually installed.
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
follow the instructions in `bepipred-1.0.readme` to complete installation
#### BepiPred 2.0
Follow the instructions in `bp2_env_install_instructions.txt` under `install` to install the bp2 conda environment
#### LBEEP 1.0
download site: [https://github.com/brsaran/LBEEP](https://github.com/brsaran/LBEEP)<br/>
follow the instructions in `README.md` to install LBEEP 1.0


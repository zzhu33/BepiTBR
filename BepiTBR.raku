#!/usr/bin/env raku

# Prediction of B cell epitopes using BepiTBR given candidate B cell epitopes
# predicted B cell epitopes have ensemble scores of >0.5 (stringent cutoff) or >0 (lenient cutoff)
# note: check carefully that the full antigen sequences do not have "*".

# dependency: raku/perl6 (>=2020.10), R (>=3.6, glmnet package available), python (>=3.6.4)
# BepiPred 1.0, BepiPred 2.0 (under a python environment, check the ./install folder), LBEEP,
# MixMHC2pred, netMHCIIpan (>=3.2, optional), tar

# motif0_file: candidate B cell epitopes in the format of a fasta file. The job id (no " " allowed) is on the ">" line.
#             It's best that the B cell epitopes are shorter than 15aa. The epitope is to be searched within the antigen sequence to
#             locate the epitope. Sometimes the same epitope appears more than once in the antigen sequence. To
#             handle such cases, a number can be optionally added after the job id, separated by "_pos=". this is
#             the 0-based index of the epitope in the antigen sequence. More on this below
# full0_file: full protein sequences in the format of a fasta file. The same job id, as in motif0_file, is on the ">" line.
#             If the motif file has the index for some records, this file also needs to have the index for the same records
# bepipred2: conda env activate script for bepipred2.
# bepipred1: executable of bepipred1.
# LBEEP: folder to LBEEP. IMPORTANT!!! This is not the executable, but the folder
# MixMHC2pred: executable to MixMHC2pred
#              (/project/shared/xiao_wang/software/MixMHC2pred/MixMHC2pred_unix)
# netMHCIIpan: executable to netMHCIIpan. this is optional. default is NA, to disable netMHCIIpan analysis.
#            (NA or /project/DPDS/Xiao_lab/shared/bcell_epitope_prediction/netMHCIIpan-3.2/netMHCIIpan)
# dir: work folder to store all jobs' work folders. This folder will be cleaned before BepiTBR starts to work
# thread: how many threads to start for the T/B cell epitope prediction jobs (BepiPred 2.0 is very slow)
# keep: keep intermediate folders of Tepi and Bepi. "true" or "false". Default is false. When set to true, will compress these two folders
#
#raku BepiTBR/BepiTBR.raku \
#--motif0_file=BepiTBR/examples/test_data_BepiTBR/Ind-positive.txt \
#--full0_file=BepiTBR/examples/test_data_BepiTBR/peptide_with_full_length.txt \
#--bepipred2=/home/conda_envs/bp2/bin/activate \
#--bepipred1=/home/bcell_epitope_prediction/bp1/bepipred-1.0/bepipred \
#--LBEEP=/home/bcell_epitope_prediction/LBEEP/ \
#--MixMHC2pred=/home/bcell_epitope_prediction/MixMHC2pred/MixMHC2pred_unix \
#--netMHCIIpan=/home/bcell_epitope_prediction/netMHCIIpan-3.2/netMHCIIpan \
#--dir=BepiTBR/examples/test_output_BepiTBR \
#--thread=20 \
#--keep=false

# example motif and full files without the index
# motif file:
# >job1
# >VDEAVM
# full file:
# >job1
# >VVDEAVMWYNL

# example motif and full files with the index
# motif file:
# >job1_pos=1
# >VDEAVM
# full file:
# >job1_pos=1
# >VVDEAVMWYNL

sub MAIN(Str :$motif0_file,Str :$full0_file,Str :$bepipred2,Str :$bepipred1,Str :$LBEEP,
    Str :$MixMHC2pred,Str :$netMHCIIpan,Str :$dir,Int :$thread,Str :$keep="false")
{
    my ($path,$dir0);
    my $length=15;
    my $expand=210;

    # find executing folders/files
    $path=$?FILE;
    $path~~s/BepiTBR\.raku//;
    $path~="/script";

    $dir0=abs_path($dir,True);
    shell("rm -f -r "~$dir0~"/Tepi");mkdir($dir0~"/Tepi");
    shell("rm -f -r "~$dir0~"/Bepi");mkdir($dir0~"/Bepi");
    shell("rm -f -r "~$dir0~"/binned_csv1.5");mkdir($dir0~"/binned_csv1.5");
    shell("rm -f -r "~$dir0~"/binned_csv2");mkdir($dir0~"/binned_csv2");

    # call T epitopes
    shell("raku "~$path~"/wrapper_T_cell_epitopes.raku "~
        "--script="~$path~"/T_cell_epitopes.raku "~
        "--motif0_file="~abs_path($motif0_file)~
        " --full0_file="~abs_path($full0_file)~
        " --expand="~$expand~" --length="~$length~
        " --MixMHC2pred="~abs_path($MixMHC2pred)~
        " --netMHCIIpan="~abs_path($netMHCIIpan)~
        " --dir="~$dir0~"/Tepi --thread="~$thread);

    # call B epitopes
    shell("raku "~$path~"/wrapper_B_cell_epitopes.raku "~
        "--script="~$path~"/B_cell_epitopes.raku "~
        "--motif0_file="~abs_path($motif0_file)~
        " --full0_file="~abs_path($full0_file)~
        " --bepipred2="~abs_path($bepipred2)~
        " --bepipred1="~abs_path($bepipred1)~
        " --LBEEP="~abs_path($LBEEP)~
        " --dir="~$dir0~"/Bepi --thread="~$thread~
        " --mode=motif");

    # organize T cell epitopes, cutoff=1.5 and 2
    shell("python "~$path~"/processed2csv_production.py "~
        "-s "~abs_path($motif0_file)~" -o "~$dir0~"/binned_csv1.5 "~
        "-t "~$dir0~"/Tepi -b "~$dir0~"/Bepi --c 1.5");
    shell("python "~$path~"/processed2csv_production.py "~
        "-s "~abs_path($motif0_file)~" -o "~$dir0~"/binned_csv2 "~
        "-t "~$dir0~"/Tepi -b "~$dir0~"/Bepi --c 2");

    # generate predictions
    shell("Rscript "~$path~"/prediction.R "~$path~" "~$dir0);

    # clean up
    if ($keep eq "true")
    {
      1;
    }else
    {
      shell("cd "~$dir0~";rm -f -r Bepi");
      shell("cd "~$dir0~";rm -f -r Tepi");
    }
}

sub abs_path(Str $path0,Bool $create=False)
{
    my $path=$path0;
    if ($path~~s/^\~//) {$path=%*ENV{"HOME"}~$path;}
    if ((!$path.IO.e) && $create) {mkdir $path;}
    $path;
}

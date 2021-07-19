#!/usr/bin/env raku

# Prediction of B cell epitopes using BepiTBR for one full antigen sequence
# dependency: same as BepiTBR.raku
# note: if "*" is detected in the middle of the protein sequence, the protein
# sequence will be truncated before "*"

# full0: full protein sequence. better to quote the sequence as there could be empty spaces by mistake
# length: length of the B cell epitopes to scan in a moving window manner. Recommend length=15
# bepipred2: same as BepiTBR.raku
# bepipred1: same as BepiTBR.raku
# LBEEP: same as BepiTBR.raku
# MixMHC2pred: same as BepiTBR.raku
# netMHCIIpan: same as BepiTBR.raku. Strongly recommend to set to NA, as it's slow
# dir: same as BepiTBR.raku
# thread: same as BepiTBR.raku
# keep: same as BepiTBR.raku

#raku /project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/BepiTBR_full.raku \
#--full0="mtenstsapaakpkrakaskkstdhpkysdmivaaiqaeknragSSRQSIQKYIKSHYKvgenadsqiklsikrlvttgvlkqtkgvgag\
#sfrlaksdepkksvafkktkkeikkvatpkkaskpkkaaskaptkkpkatpvkkakkklaatpkkakkpktvkakpvkaskpkkakpvkpkakssakragkkk" \
#--length=15 \
#--bepipred2=/project/shared/xiao_wang/projects/Bcell_epitope/code/conda_envs/bp2/bin/activate \
#--bepipred1=/project/DPDS/Xiao_lab/shared/bcell_epitope_prediction/bp1/bepipred-1.0/bepipred \
#--LBEEP=/project/shared/xiao_wang/software/LBEEP/ \
#--MixMHC2pred=/project/shared/xiao_wang/software/MixMHC2pred/MixMHC2pred_unix \
#--netMHCIIpan=NA \
#--dir=/project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/example/test_output_BepiTBR_full \
#--thread=20 \
#--keep=false

sub MAIN(Str :$full0,Int :$length,Str :$bepipred2,Str :$bepipred1,Str :$LBEEP,
    Str :$MixMHC2pred,Str :$netMHCIIpan,Str :$dir,Int :$thread,Str :$keep="false")
{
    my ($path,$i,@tmp,$dir0,$full);
    my $tepi_length=15;

    # find executing folders/files
    $path=$?FILE;
    $path~~s/BepiTBR_full\.raku//;
    $path~="/script";

    $dir0=abs_path($dir,True);
    shell("rm -f -r "~$dir0~"/input");mkdir($dir0~"/input");
    shell("rm -f -r "~$dir0~"/Tepi");mkdir($dir0~"/Tepi");
    shell("rm -f -r "~$dir0~"/Bepi");mkdir($dir0~"/Bepi");
    shell("rm -f -r "~$dir0~"/binned_csv1.5");mkdir($dir0~"/binned_csv1.5");
    shell("rm -f -r "~$dir0~"/binned_csv2");mkdir($dir0~"/binned_csv2");

    # handle early truncation
    $full=$full0;
    $full~~s/\*.*//;
    if ( $full.chars-$length < 0 ) {die "Warning: Abort as the antigen sequence is too short!";}

    # call T epitopes
    shell("raku "~$path~"/T_cell_epitopes.raku"~
        " --motif0="~substr($full,1,min(10,$full.chars-1))~" --full0="~
        $full~" --expand="~($full.chars+100)~" --length="~$tepi_length~
        " --MixMHC2pred="~abs_path($MixMHC2pred)~" --netMHCIIpan="~
        abs_path($netMHCIIpan)~" --dir="~$dir0~"/input/Tepi --id=job");

    for (0..$full.chars-$length) -> $i
    {
        shell("cp -r "~$dir0~"/input/Tepi "~$dir0~"/Tepi/job_pos="~$i);
        shell("tail -n 3 "~$dir0~"/Tepi/job_pos="~$i~"/information.txt >"~
                $dir0~"/Tepi/job_pos="~$i~"/tmp.txt");
        spurt $dir0~"/Tepi/job_pos="~$i~"/tmp1.txt", "Bepi_pos\t"~$i~"\n";
        shell("cat "~$dir0~"/Tepi/job_pos="~$i~"/tmp1.txt "~
                $dir0~"/Tepi/job_pos="~$i~"/tmp.txt >"~
                $dir0~"/Tepi/job_pos="~$i~"/information.txt");
        shell("rm -f "~$dir0~"/Tepi/job_pos="~$i~"/tmp*.txt");
    }

    # call B epitopes
    @tmp=(0..$full.chars-$length).map({">job_pos="~$_~"\n"~$full~"\n"});
    spurt $dir0~"/input/full.txt", [~] @tmp;
    @tmp=(0..$full.chars-$length).map({">job_pos="~$_~"\n"~substr($full,$_,$length)~"\n"});
    spurt $dir0~"/input/motif.txt", [~] @tmp;

    shell("raku "~$path~"/wrapper_B_cell_epitopes.raku "~
        "--script="~$path~"/B_cell_epitopes.raku "~
        "--motif0_file="~$dir0~"/input/motif.txt --full0_file="~$dir0~"/input/full.txt "~
        " --bepipred2="~abs_path($bepipred2)~" --bepipred1="~abs_path($bepipred1)~
        " --LBEEP="~abs_path($LBEEP)~" --dir="~$dir0~"/Bepi --thread="~$thread~
        " --mode=full");

    # organize T cell epitopes, cutoff=1.5 and 2
    shell("python "~$path~"/processed2csv_production.py "~
        "-s "~$dir0~"/input/motif.txt -o "~$dir0~"/binned_csv1.5 "~
        "-t "~$dir0~"/Tepi -b "~$dir0~"/Bepi --c 1.5");
    shell("python "~$path~"/processed2csv_production.py "~
        "-s "~$dir0~"/input/motif.txt -o "~$dir0~"/binned_csv2 "~
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

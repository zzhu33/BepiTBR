#!/usr/bin/env raku

# Prediction of B cell epitopes using BepiTBR for a fasta file of >=1 full antigen sequences
# dependency: same as BepiTBR_full.raku
# note: if "*" is detected in the middle of the protein sequence, the protein
# sequence will be truncated before "*"

# fasta0: fasta file.
# length: same as BepiTBR_full.raku
# bepipred2: same as BepiTBR.raku
# bepipred1: same as BepiTBR.raku
# LBEEP: same as BepiTBR.raku
# MixMHC2pred: same as BepiTBR.raku
# netMHCIIpan: same as BepiTBR_full.raku.
# dir: output dir. for each antigen sequence, a sub-dir will be built
# thread: same as BepiTBR.raku
# keep: same as BepiTBR.raku

#raku BepiTBR_fasta.raku \
#--fasta0=examples/test_data_BepiTBR_fasta/peptide_with_full_length.txt \
#--length=15 \
#--bepipred2=/home/conda_envs/bp2/bin/activate \
#--bepipred1=/home/bcell_epitope_prediction/bp1/bepipred-1.0/bepipred \
#--LBEEP=/home/bcell_epitope_prediction/LBEEP/ \
#--MixMHC2pred=/home/bcell_epitope_prediction/MixMHC2pred/MixMHC2pred_unix \
#--netMHCIIpan=NA \
#--dir=examples/test_output_BepiTBR_fasta \
#--thread=8 \
#--keep=false

sub MAIN(Str :$fasta0,Int :$length,Str :$bepipred2,Str :$bepipred1,Str :$LBEEP,
    Str :$MixMHC2pred,Str :$netMHCIIpan,Str :$dir,Int :$thread,Str :$keep="false")
{
    my ($fasta,$path,%full0,$a,$dir0,$id,$combinePath);
    my $tepi_length=15;

    # find executing folders/files
    $path=$?FILE;
    $path~~s/fasta/full/;

    # input/output
    $dir0=abs_path($dir,True);
    $fasta=abs_path($fasta0);

    for $fasta.IO.lines -> $a
    {
        if ($a.starts-with(">"))
        {
            $id=$a.substr(1).chomp;
            %full0{$id}="";
        }else
        {
            %full0{$id}~=$a.chomp;
        }
    }

    # execute BepiTBR.raku
    for %full0.keys -> $id
    {
        say "Executing "~$id;
        shell("raku "~$path~"  --full0="~%full0{$id}~" --length="~$length~
            " --bepipred2="~$bepipred2~" --bepipred1="~$bepipred1~
            " --LBEEP="~$LBEEP~" --MixMHC2pred="~$MixMHC2pred~
            " --netMHCIIpan="~$netMHCIIpan~" --dir="~$dir0~"/"~$id~
            " --thread="~$thread~" --keep="~$keep);
    }
	# combine results
	$combinePath = $?FILE;
	$combinePath.=subst(/BepiTBR_fasta.raku/, "combine_fasta_results.py", :nth(*));
	shell("python "~$combinePath~" -i "~$dir0~" "~$fasta);
    # final cleanup
    shell("tar -zcvf "~$dir0~"/../output.tar.gz -C "~$dir0~" .");
    shell("mv "~$dir0~"/../output.tar.gz "~$dir0);
    for %full0.keys -> $id {shell("rm -f -r "~$dir0~"/"~$id);}
}

sub abs_path(Str $path0,Bool $create=False)
{
    my $path=$path0;
    if ($path~~s/^\~//) {$path=%*ENV{"HOME"}~$path;}
    if ((!$path.IO.e) && $create) {mkdir $path;}
    $path;
}

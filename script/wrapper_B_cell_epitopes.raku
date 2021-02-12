# script: path to the B_cell_epitopes.raku script
#         (/project/shared/xiao_wang/projects/Bcell_epitope/code/T_cell_epitopes.raku)
# motif0_file: candidate B cell epitopes in the format of a fasta file. follows the format of the iBCE data
#             Importantly, the sequence of the same epitope must be on the same line
# full0_file: full protein sequences in the format of a fasta file. follows iBCE.
#             Importantly, the sequence of the same protein must be on the same line
# bepipred2: conda env activate script for bepipred2 (created by James Zhu).
# bepipred1: executable of bepipred1.
# LBEEP: folder to LBEEP. IMPORTANT!!! This is not the executable, but the folder
# dir: work folder to store all jobs' work folders.
# thread: #threads to start
# mode: motif (the base mode), or full (the full mode)

#raku /project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/script/wrapper_B_cell_epitopes.raku \
#--script=/project/shared/xiao_wang/projects/Bcell_epitope/code/B_cell_epitopes.raku \
#--motif0_file=/project/shared/xiao_wang/projects/Bcell_epitope/data/downloaded_raw/iBCEData/Ind-positive.txt \
#--full0_file=/project/shared/xiao_wang/projects/Bcell_epitope/data/downloaded_raw/iBCEData/peptide_with_full_length.txt \
#--bepipred2=/project/shared/xiao_wang/projects/Bcell_epitope/code/conda_envs/bp2/bin/activate \
#--bepipred1=/project/bioinformatics/Xiao_lab/shared/bcell_epitope_prediction/bp1/bepipred-1.0/bepipred \
#--LBEEP=/project/shared/xiao_wang/software/LBEEP/ \
#--dir=/home2/twang6/iproject/test \
#--thread=20 \
#--mode=motif

#!/usr/bin/env raku

sub MAIN(Str :$script,Str :$motif0_file,Str :$full0_file,Str :$bepipred2,
    Str :$bepipred1,Str :$LBEEP,Str :$dir,Int :$thread,Str :$mode)
{
    my (%motif0,%full0,$a,$b,@promises,$i,@keys,$n_jobs,$id);

    # read the B cell epitope and full protein seq files
    for $motif0_file.IO.lines -> $a
    {
        if ($a.starts-with(">"))
        {
            $id=$a.substr(1).chomp;
            %motif0{$id}="";
        }else
        {
            %motif0{$id}~=$a.chomp;
        }
    }

    for $full0_file.IO.lines -> $a
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

    # special handling for full mode
    if ($mode leg "full"==Same)
    {
        spurt $dir~"/../input/full.fasta",">job\n"~%full0.values()[0];
        shell($bepipred1~" "~$dir~"/../input/full.fasta | grep -v \"#\" | awk \'\{print \$6\}\' >"~
            $dir~"/../input/bepipred1.0.out.tmp");
        shell("source "~$bepipred2~";BepiPred-2.0 "~$dir~
            "/../input/full.fasta | grep -v \"#\" | awk \'\{print \$8\}\' >"~
            $dir~"/../input/bepipred2.0.out.tmp;source deactivate");
    }

    # run B_cell_epitopes.raku
    @keys=%motif0.keys();
    $n_jobs=@keys.elems div $thread;
    if ($n_jobs * $thread <@keys.elems) {$n_jobs++;}

    for (1..$thread) -> $i
    {
        @promises.push(start
        {
            for @keys[($i-1)*$n_jobs..($i*$n_jobs-1,@keys.elems-1).min] -> $a
            {
                ("Working on "~$a).say;
                shell("raku "~$script~" --motif0="~%motif0{$a}~" --full0="~
                    %full0{$a}~" --bepipred2="~$bepipred2~" --bepipred1="~$bepipred1~
                    " --LBEEP="~$LBEEP~" --dir="~$dir~"/"~$a~" --id="~$a~" --mode="~$mode);
            }
        })
    }

    await(|@promises);
}


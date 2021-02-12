# script: path to the T_cell_epitopes.raku script
#         (/project/shared/xiao_wang/projects/Bcell_epitope/code/T_cell_epitopes.raku)
# motif0_file: candidate B cell epitopes in the format of a fasta file. follows the format of the iBCE data
#             Importantly, the sequence of the same epitope must be on the same line
# full0_file: full protein sequences in the format of a fasta file. follows iBCE.
#             Importantly, the sequence of the same protein must be on the same line
# expand: search for MHC class II epitopes between (-expand, expand) of the candidate B cell epitope
# length: length of MHC class II epitopes (15)
# MixMHC2pred: executable to MixMHC2pred
#              (/project/shared/xiao_wang/software/MixMHC2pred/MixMHC2pred_unix)
# netMHCIIpan: executable to netMHCIIpan
#            (/project/bioinformatics/Xiao_lab/shared/bcell_epitope_prediction/netMHCIIpan-3.2/netMHCIIpan)
# dir: work folder to store all jobs' work folders.
# thread: #threads to start

#raku /project/shared/xiao_wang/projects/Bcell_epitope/code/wrapper_T_cell_epitopes.raku \
#--script=/project/shared/xiao_wang/projects/Bcell_epitope/code/T_cell_epitopes.raku \
#--motif0_file=/project/shared/xiao_wang/projects/Bcell_epitope/data/downloaded_raw/iBCEData/Ind-positive.txt \
#--full0_file=/project/shared/xiao_wang/projects/Bcell_epitope/data/downloaded_raw/iBCEData/peptide_with_full_length.txt \
#--expand=50 --length=15 \
#--MixMHC2pred=/project/shared/xiao_wang/software/MixMHC2pred/MixMHC2pred_unix \
#--netMHCIIpan=/project/bioinformatics/Xiao_lab/shared/bcell_epitope_prediction/netMHCIIpan-3.2/netMHCIIpan \
#--dir=/home2/twang6/iproject/test \
#--thread=20

#!/usr/bin/env raku

sub MAIN(Str :$script,Str :$motif0_file,Str :$full0_file,Int :$expand,Int :$length,Str :$MixMHC2pred,
        Str :$netMHCIIpan,Str :$dir,Int :$thread)
{
    my ($dir0,%motif0,%full0,$a,$b,@promises,$i,@keys,$n_jobs,$id);

    # prepare
    $dir0=$dir;
    unless ($dir0.IO.e) {mkdir $dir0;}

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

    # run T_cell_epitopes.raku
    @keys=%motif0.keys();
    $n_jobs=@keys.elems div $thread;
    if ($n_jobs * $thread <@keys.elems) {$n_jobs++;}

    for (1..$thread) -> $i
    {
        @promises.push(start
        {
            for @keys[($i-1)*$n_jobs..($i*$n_jobs-1,@keys.elems-1).min] -> $a
            {
                ("\nWorking on "~$a).say;
                ("Motif=\""~%motif0{$a}~"\"").say;
                ("Full=\""~%full0{$a}~"\"").say;
                shell("raku "~$script~" --motif0="~%motif0{$a}~" --full0="~
                    %full0{$a}~" --expand="~$expand~" --length="~$length~
                    " --MixMHC2pred="~$MixMHC2pred~" --netMHCIIpan="~
                    $netMHCIIpan~" --dir="~$dir0~"/"~$a~" --id="~$a);
            }
        })
    }

    await(|@promises);
}


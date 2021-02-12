# motif0: candidate B cell epitope.
#         for training/testing mode: follow T_cell_epitopes.raku. will carry out computation for this epitope
#         for production mode: specify as "PRODUCTION". this will carry out computation for all potential epitopes
# full0: full protein sequence. follow T_cell_epitopes.raku
# bepipred2: conda env activate script for bepipred2 (created by James Zhu).
# bepipred1: executable of bepipred1.
# LBEEP: folder to LBEEP. IMPORTANT!!! This is not the executable, but the folder
# dir: work folder. Unique to each job.
# id: job id
# mode: motif (the base mode), or full (the full mode)

#raku /project/shared/xiao_wang/projects/Bcell_epitope/code/BepiTBR/script/B_cell_epitopes.raku \
#--motif0="SSRQSIQKYIKSHYK" \
#--full0="mtenstsapaakpkrakaskkstdhpkysdmivaaiqaeknr\
#agSSRQSIQKYIKSHYKvgenadsqiklsikrlvttgvlkqtkgvgag\
#sfrlaksdepkksvafkktkkeikkvatpkkaskpkkaaskaptkkpk\
#atpvkkakkklaatpkkakkpktvkakpvkaskpkkakpvkpkakssa\
#kragkkk" \
#--bepipred2=/project/shared/xiao_wang/projects/Bcell_epitope/code/conda_envs/bp2/bin/activate \
#--bepipred1=/project/DPDS/Xiao_lab/shared/bcell_epitope_prediction/bp1/bepipred-1.0/bepipred \
#--LBEEP=/project/shared/xiao_wang/software/LBEEP/ \
#--dir=/home2/twang6/iproject/test \
#--id=job1 \
#--mode=base

#!/usr/bin/env raku

sub MAIN(Str :$motif0,Str :$full0,Str :$bepipred2,Str :$bepipred1,Str :$LBEEP,Str :$dir,Str :$id,
    Str :$mode)
{
    my ($motif,$full,$pos,$output,@output);

    # clean up the full information
    $full=clean_sequence($full0,"full sequence",$motif0);
    if (!$dir.IO.e) {mkdir $dir;}
    spurt $dir~"/full.fasta",">"~$id~"\n"~$full~"\n";

    # clean up the motif and pos information
    $motif=clean_sequence($motif0,"B cell epitope",$motif0);
    if ( $id ~~ / _pos\= (\d+) $ / )
    {
        $pos=$0;
        if (substr($full,$pos,$motif.chars) leg $motif != Same)
            {die "Error in "~$id~": the retrieved motif by pos doesn't match the specified motif!\n";}
    }else
    {
        $pos=$full.index($motif); # if there are >=2 matches, only the first one will be reported
        unless ($pos.defined) {die "Error in "~$id~": cannot find the B cell epitope!\n";}
    }
    spurt $dir~"/motif.fasta",">"~$id~"\n"~$motif~"\n";

    #LBEEP
    shell("cd "~$LBEEP~";./LBEEP -i "~$dir~"/motif.fasta -m pep -M C -o "~$dir~"/LBEEP.out.tmp -t 0.001");
    shell("tail -n 1 "~$dir~"/LBEEP.out.tmp | cut -d , -f 3 > "~$dir~"/LBEEP.out");
    shell("rm -f "~$dir~"/LBEEP.out.tmp"); # the output is the confidence score of LBEEP

    # bepipred1.0
    if ($mode leg "full"==Same)
    {
        shell("cp "~$dir~"/../../input/bepipred1.0.out.tmp "~$dir~"/bepipred1.0.out.tmp")
    }else
    {
        if (! $bepipred1.IO.e) {die "BepiPred1.0 not found!";}
        shell($bepipred1~" "~$dir~"/full.fasta | grep -v \"#\" | awk \'\{print \$6\}\' >"~
            $dir~"/bepipred1.0.out.tmp");
    }
    $output=slurp ($dir~"/bepipred1.0.out.tmp");
    @output=$output.split("\n")[$pos..($pos+$motif.chars-1)].map({+$_});
    $output=@output.max;
    #$output=([+] @output)/@output.elems;
    spurt $dir~"/bepipred1.0.out",$output;
    shell("rm -f "~$dir~"/bepipred1.0.out.tmp"); # the output is the confidence score of LBEEP

    # bepipred2.0
    if ($mode leg "full"==Same)
    {
        shell("cp "~$dir~"/../../input/bepipred2.0.out.tmp "~$dir~"/bepipred2.0.out.tmp")
    }else
    {
        if (! $bepipred2.IO.e) {die "BepiPred2.0 not found!";}
        shell("source "~$bepipred2~";BepiPred-2.0 "~$dir~
            "/full.fasta | grep -v \"#\" | awk \'\{print \$8\}\' >"~
            $dir~"/bepipred2.0.out.tmp;source deactivate");
    }
    $output=slurp ($dir~"/bepipred2.0.out.tmp");
    @output=$output.split("\n")[$pos..($pos+$motif.chars-1)].map({+$_});
    $output=@output.max;
    #$output=([+] @output)/@output.elems;
    spurt $dir~"/bepipred2.0.out",$output;
    shell("rm -f "~$dir~"/bepipred2.0.out.tmp"); # the output is the confidence score of LBEEP

    # cleanup
    #shell("rm -f "~$dir~"/full.fasta");
    #shell("rm -f "~$dir~"/motif.fasta");
}

sub clean_sequence(Str $protein0,Str $type,Str $job)
{
    my $protein=$protein0.uc.chomp;
    #my $temp=$protein;
    #$temp~~s/[A|R|N|D|C|E|Q|G|H|I|L|K|M|F|P|S|T|W|Y|V|X]+//;
    #if ($temp.chars>0) {die $job~": the "~$type~" is not a protein sequence!"}
    $protein;
}

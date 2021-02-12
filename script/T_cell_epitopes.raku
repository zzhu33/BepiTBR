# motif0: candidate B cell epitope. better to quote the sequence as there could be empty spaces by mistake
# full0: full protein sequence. better to quote the sequence as there could be empty spaces by mistake
# expand: search for MHC class II epitopes between (-expand, expand) of the candidate B cell epitope
# length: length of MHC class II epitopes (15)
# MixMHC2pred: executable to MixMHC2pred
#              (/project/shared/xiao_wang/software/MixMHC2pred/MixMHC2pred_unix)
# netMHCIIpan: executable to netMHCIIpan
#            (/project/bioinformatics/Xiao_lab/shared/bcell_epitope_prediction/netMHCIIpan-3.2/netMHCIIpan)
# dir: work folder. Unique to each job.
# id: job id

#raku /project/shared/xiao_wang/projects/Bcell_epitope/code/T_cell_epitopes.raku \
#--motif0="SSRQSIQKYIKSHYK" \
#--full0="mtenstsapaakpkrakaskkstdhpkysdmivaaiqaeknr\
#agSSRQSIQKYIKSHYKvgenadsqiklsikrlvttgvlkqtkgvgag\
#sfrlaksdepkksvafkktkkeikkvatpkkaskpkkaaskaptkkpk\
#atpvkkakkklaatpkkakkpktvkakpvkaskpkkakpvkpkakssa\
#kragkkk" --expand=50 --length=15 \
#--MixMHC2pred=/project/shared/xiao_wang/software/MixMHC2pred/MixMHC2pred_unix \
#--netMHCIIpan=/project/bioinformatics/Xiao_lab/shared/bcell_epitope_prediction/netMHCIIpan-3.2/netMHCIIpan \
#--dir=/home2/twang6/iproject/test \
#--id=job1

#!/usr/bin/env perl6

sub MAIN(Str :$motif0,Str :$full0,Int :$expand,Int :$length,Str :$MixMHC2pred,Str :$netMHCIIpan,Str :$dir,
    Str :$id)
{
    my ($motif,$full,$pos,$i,$epitope,$start,$end,$dir0,$temp);
    my ($line,@items,$output,$a,@b,$allele,@promises,$file);

    # hidden parameters
    my $rank_cutoff=2;

    # prepare
    $dir0=$dir;
    unless ($dir0.IO.e) {mkdir $dir0;}
    $motif=clean_sequence($motif0,"B cell epitope",$motif0);
    $full=clean_sequence($full0,"full sequence",$motif0);

    # find the position of the B cell epitope for MixMHC2pred
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

    $start=0 max ($pos-$expand);
    $end=($full.chars-$length) min ($pos+$motif.chars+$expand-$length);
    spurt $dir0~"/information.txt", "Bepi_pos\t"~$pos~"\nlength_full\t"~$full.chars~
        "\nstart\t"~$start~"\nend\t"~$end~"\n";
    if ($start>$end) {die "No candidate T cell epitopes\n";}

    # find the T cell epitopes for netMHCIIpan
    $epitope=substr($full,$start..$end+$length-1);
    spurt $dir0~"/netMHCIIpan.txt", ">epitope\n"~$epitope;

    # find the T cell epitopes for MixMHC2pred
    $epitope="";
    for $start..$end -> $i
      {$epitope~=(substr($full,$i..($i+$length-1))~"\n");}
    spurt $dir0~"/MixMHC2pred.txt", $epitope;

    # run MixMHC2pred
    shell($MixMHC2pred~" -i "~$dir0~"/MixMHC2pred.txt -o "~$dir0~
        "/MixMHC2pred_out.txt -a "~alleles("MixMHC2pred"));
    shell("rm -f "~$dir0~"/MixMHC2pred.txt");

    $output=""; # cleanup for MixMHC2pred
    $i=$start-1;
    for ($dir0~"/MixMHC2pred_out.txt").IO.lines -> $line
    {
        if ($line~~m/^\#/) {next;} # header
        @items=$line.split("\t")[(0,(6,*+2...*))];
        $a=@items[0];
        @b=@items[1].flat;
        if ($a eq "Peptide")
        {
            @b=@b.map({$_~~s/.*Rank_//;$_});
            $a~=" Pos"; # attach position of the T cell epitope
        }else
        {
            $i++;
            if (@b[1] eq "NA") {next;} # caused by amino acid of X
            if (@b.map({+$_}).min>$rank_cutoff) {next;}
            $a~=(" "~$i);
        }

        $output~=($a~" "~@b.join(" ")~"\n");
    };
    spurt $dir0~"/tmp.txt", $output;
    shell("mv "~$dir0~"/tmp.txt "~$dir0~"/MixMHC2pred_out.txt");

    # run netMHCIIpan
    if ($netMHCIIpan leg "NA" != Same)
    {
        for alleles("netMHCIIpan") -> $allele
        {
            @promises.push(start
            {
                shell($netMHCIIpan~" -tdir "~$dir0~"/"~$allele~".tmp -inptype 0 -length "~
                   $length~" -a "~$allele~" -f "~$dir0~"/netMHCIIpan.txt | grep epitope > "~
                    $dir0~"/netMHCIIpan_"~$allele~".txt");
            })
        }

        await(|@promises);
        @items=(); # binding affinities
        $i=0;

        for alleles("netMHCIIpan") -> $allele
        {
            $file=$dir0~"/netMHCIIpan_"~$allele~".txt";
            @items[$i]=Array.new; # read output
            for ($file).IO.lines
            {
                push @items[$i],+.trim.split(/\s+/)[9];
            }
            if ($i==0) {@b=$file.IO.lines.hyper.map({.split(/\s+/)[3]});} # epitope names
            shell("rm -f "~$file);
            $i++;
        }

        @items=[Z] @items; # print everything
        spurt $dir0~"/netMHCIIpan_out.txt","Epitope Pos "~alleles("MixMHC2pred")~"\n";
        for (0..@b.elems-1) -> $i
        {
            if (@items[$i].min>$rank_cutoff) {next;}
            spurt $dir0~"/netMHCIIpan_out.txt",@b[$i]~" "~($i+$start)~" "~
                    @items[$i].join(" ")~"\n",:append;
        }

        shell("rm -f "~$dir0~"/netMHCIIpan.txt");
    }else
    {
        # just give a blank txt file
        spurt $dir0~"/netMHCIIpan_out.txt","Epitope Pos "~alleles("MixMHC2pred")~"\n";
        shell("rm -f "~$dir0~"/netMHCIIpan.txt");
    }
}

sub clean_sequence(Str $protein0,Str $type,Str $job)
{
    my $protein=$protein0.uc.chomp;
    #my $temp=$protein;
    #$temp~~s/[A|R|N|D|C|E|Q|G|H|I|L|K|M|F|P|S|T|W|Y|V|X]+//;
    #if ($temp.chars>0) {die $job~": the "~$type~" is not a protein sequence!"}
    $protein;
}

sub alleles(Str $type)
{
    my $alleles="DRB1_01_01 DRB1_01_02 DRB1_01_03 DRB1_03_01 DRB1_04_01 DRB1_04_04 DRB1_04_05 DRB1_04_08 "~
    "DRB1_07_01 DRB1_08_01 DRB1_10_01 DRB1_11_01 DRB1_11_04 DRB1_12_01 DRB1_13_01 DRB1_13_03 "~
    "DRB1_15_01 DRB1_16_01 DRB3_01_01 DRB3_02_02 DRB4_01_01 DRB4_01_03 DRB5_01_01 DRB5_02_02 "~
    "DPA1_01_03__DPB1_02_01 DPA1_01_03__DPB1_03_01 DPA1_01_03__DPB1_04_01 DPA1_01_03__DPB1_06_01 "~
    "DPA1_01_03__DPB1_20_01 DQA1_02_01__DQB1_02_02 "~
    "DQA1_03_01__DQB1_03_01 DQA1_03_03__DQB1_03_01 DQA1_05_05__DQB1_03_01";
    my $allele;

    if ($type eq "MixMHC2pred")
    {
        return $alleles;
    }else # netMHCIIpan
    {
        return $alleles.split(" ").map:
        {
            $allele=$_;
            if ($allele~~s/__/-/) # DQ*/DP*
            {
                $allele~~s:g/_//;
                $allele="HLA-"~$allele;
            }else # DRB*
            {
                $allele~~m/(DRB\d)_(\d+)_(\d+)/;
                $allele=$0~"_"~$1~$2;
            }
            $allele;
        }
    }
}

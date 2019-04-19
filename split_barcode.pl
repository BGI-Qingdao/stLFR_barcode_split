#!/usr/bin/perl

use strict;
if(@ARGV != 4)
{
    print "#-------------------------------------------------------------------------------\n";
    print "Usage    : perl split_barcode.pl barcode.list r1 r2 output_prefix \n";
    print "Example  : perl split_barcode.pl barcode.list r1.fq.gz r2.fq.gz split_reads \n";
    print "#-------------------------------------------------------------------------------\n";
    exit(1);
}

my @line;

#
# Step 1.
#       0. make barcode fault-tolerant hash
#           that map single barcode to a num ID ;
#

print "step 1 : load barcodes ... \n";
$|=1;
my %barcode_hash;
open IN,"$ARGV[0]" or die "cann't not open barcode.list";
my $index = 0;
while(<IN>)
{
    $index ++;
    my @line = split;
    my @barcode = split(//,$line[0]);
    my $barcode_ID = $line[1];
    for(my $num = 0; $num <= 9; $num++)
    {
        my @barcode_mis = @barcode;
        $barcode_mis[$num] = "A";
        my $barcode_mis = join("",@barcode_mis);
        $barcode_hash{$barcode_mis} = $barcode_ID;
        @barcode_mis = @barcode;
        $barcode_mis[$num] = "G";
        my $barcode_mis = join("",@barcode_mis);
        $barcode_hash{$barcode_mis} = $barcode_ID;
        @barcode_mis = @barcode;
        $barcode_mis[$num] = "C";
        my $barcode_mis = join("",@barcode_mis);
        $barcode_hash{$barcode_mis} = $barcode_ID;
        @barcode_mis = @barcode;
        $barcode_mis[$num] = "T";
        my $barcode_mis = join("",@barcode_mis);
        $barcode_hash{$barcode_mis} = $barcode_ID;
    }
}
close IN;

my $line_num = 0;
#
# Step 2.
#       0. detect barcodes type from r2
print "step 2 : detect barcodes type from r2 ... \n";
$|=1;
my $n1 = 10;
my $n2 = 6;
my $n3 = 10;
my $n4 = 0;
my $n5 = 10;
my $valid_read_len = 100  ;
my $r2_len = 0 ;
open IN2,"gzip -dc $ARGV[2] |" or die "cannot open file";
while(<IN2>)
{
    chomp;
    @line = split;
    $line_num ++;
    if( $line_num == 2 )
    {
        my $len = length($line[0]);
        $r2_len = $len ;
        if( $len == 154 )
        {
            $n4 = 18 ;
        }
        elsif ( $len == 142 )
        {
            $n4 = 6 ;
        }
        elsif ( $len == 126 )
        {
            $n4 = 0 ;
            $n5 = 0;
        }
        else
        {
            print "Unknow read len of r2 : $len . Please double check the $ARGV[2] file.";
            exit(1);
        }
        last;
    }
}
close IN2 ; 
my $barcode_types = $index * $index *$index;
my $barcode_each = $index;
print "Barcode_types = $barcode_each * $barcode_each * $barcode_each = $barcode_types\n";
$|=1;


#   Step 3
#       2. make barcode_string --> barcode_num map.
#       3. make barcode_num --> barcode_string map.
#       4. calculate barcode_num frequence .
#       5. print stats.

my %barcode_str_2_num_hash;
my %barcode_num_2_str_hash;
my %barcode_freq_hash;

$barcode_str_2_num_hash{"0_0_0"} = 0;
$barcode_num_2_str_hash{0} = "0_0_0";

my $reads_num = 0;
my $split_reads_num = 0 ;
my $split_barcode_num = 0;
my $progress = 0 ;

print "step 3 : parse barcodes .... \n";
$|=1;

open IN_r1,"gzip -dc $ARGV[1] |" or die "cannot open $ARGV[1] for read \n";
open OUT_r1, "| gzip > $ARGV[3].1.fq.gz" or die "Can't open $ARGV[3].1.fq.gz for write";
open IN_r2,"gzip -dc $ARGV[2] |" or die "cannot open $ARGV[2] for read \n";
open OUT_r2, "| gzip > $ARGV[3].2.fq.gz" or die "Can't open $ARGV[3].2.fq.gz for write \n";
$line_num = 0;
while(<IN_r2>)
{
    chomp;
    my $R1_head=$_
    my $R2_seq=<IN_r2>
    my $R2_3=<IN_r2>
    my $R2_qua=<IN_r2>

    my $R1_head=<IN_r1>
    my $R1_seq=<IN_r1>
    my $R1_3=<IN_r1>
    my $R1_qua=<IN_r1>

    my @heads1  = split(/\//,$R1_head);
    my $id1 = $heads1[0];
    my $flag1=$heads1[1];
    my @heads2  = split(/\//,$R2_head);
    my $id2 = $heads2[0];
    my $flag2=$heads2[1];

    if( $id1 ne $id2 )
    {
        die "head not match error !!! $id1 != $id2";
    }
    if( $flag1 != 1 || $flag2 != 2 )
    {
        die "head not match error !!! $R1_head != $R2_head";
    }

    $reads_num ++;
    $line_num = $line_num + 4 ;
    # print process ...
    if($line_num % 4000000 == 1)
    {
        print "parse barcodes processed $progress (M) reads ...\n";
        $|=1;
        $progress ++ ;
    }

    # check barcodes 
    my $b1 = substr($R2_seq, $valid_read_len, $n1);
    my $b2 = substr($R2_seq, $valid_read_len+$n1+$n2, $n3);
    my $b3 = substr($R2_seq, $valid_read_len+$n1+$n2+$n3+$n4, $n5);
    my $R2_true_seq=substr($R2_seq,0,$valid_read_len);
    my $R2_true_qua=substr($R2_qua,0,$valid_read_len);
    my $barcode_str="0_0_0";
    if((exists $barcode_hash{$b1}) && (exists $barcode_hash{$b2}) && ($n5 != 0 && (exists $barcode_hash{$b3})) )
    {
        my $str = $barcode_hash{$b1}."_".$barcode_hash{$b2};
        if( $n5 != 0 )
        {
            $str = $str."_".$barcode_hash{$b3};
        }
        if(!(exists $barcode_str_2_num_hash{$str})) {
            $split_barcode_num ++;
            $barcode_str_2_num_hash{$str} = $split_barcode_num;
            $barcode_num_2_str_hash{$split_barcode_num} = $str;
            $barcode_freq_hash{$str} = 0;
        }
        $split_reads_num ++;
        $barcode_freq_hash{$str} ++;
        $barcode_str=$str;
    }
    print OUT_r1 $id1."\#$barcode_str\/1\t$barcode_str_2_num_hash{$barcode_str}\t1\n";
    print OUT_r1 "$R1_seq\n"
    print OUT_r1 "$R1_3\n"
    print OUT_r1 "$R1_qua\n"

    print OUT_r2 $id2."\#$barcode_str\/2\t$barcode_str_2_num_hash{$barcode_str}\t1\n";
    print OUT_r2 "$R2_true_seq\n"
    print OUT_r2 "$R2_3\n"
    print OUT_r2 "$R2_true_qua\n"
}
close IN_r1 ;
close OUT_r1 ;
close IN_r2 ;
close OUT_r2 ;

# print stat
my $r1 = 0 ;
my $r2 = 0 ;
$r1 = 100 *  $split_barcode_num/$barcode_types;
$r2 = 100 *  $split_reads_num/$reads_num;
open OUT, ">split_stat_read.log" or die "ERROR : Can't open split_stat_read.log for write \n";
print OUT "Barcode_types = $barcode_each * $barcode_each * $barcode_each = $barcode_types\n";
print OUT "Real_Barcode_types = $split_barcode_num ($r1 %)\n";
print OUT "Reads_pair_num  = $reads_num \n";
print OUT "Reads_pair_num(after split) = $split_reads_num ($r2 %)\n";
print OUT "read2_read_len : $r2_len \n";
print OUT "barcode type : $valid_read_len-$n1-$n2-$n3-$n4-$n5 \n";
close OUT;
#print barcode details 
open  OUT1, ">barcode_freq.txt" or die "ERROR : Can't open barcode_freq.txt for write !!! \n";
print OUT1 "bacode_str\tbarcode_count\tbarcode_num\n";
for(my $i=1;$i<=$split_barcode_num;$i++){
    my $str = $barcode_num_2_str_hash{$i} ;
    print OUT1 "$str\t$barcode_freq_hash{$str}\t$i\n";
}
close OUT1;

print "all done!\n";
$|=1;

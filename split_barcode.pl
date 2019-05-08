#!/usr/bin/perl -w
use strict;

if(@ARGV != 3){
    print "#-------------------------------------------------------------------------------\n";
    print "Usage    : perl split_barcode.pl barcode.list r1 r2 \n";
    print "Example  : perl split_barcode.pl barcode.list r1.fq.gz r2.fq.gz \n";
    print "#-------------------------------------------------------------------------------\n";
    exit(1);
}

#
# Step 1.
#       0. make barcode fault-tolerant hash
#           that map single barcode to a num ID ;
#

print "Step 1: Loading barcodes ... \n";
$|=1;
my %all_barcode;
open IN,"$ARGV[0]" or die "cann't not open barcode.list";
my $index = 0;
my $chek;
while(<IN>){
    chomp;
    $index++;
    my @line=split;
    foreach my $k(0..9){
        my $tmp_barcode=$line[0];
        substr($tmp_barcode,$k,1)="A";
        $all_barcode{$tmp_barcode}=$line[1];
        $tmp_barcode=$line[0];
        substr($tmp_barcode,$k,1)="T";
        $all_barcode{$tmp_barcode}=$line[1];
        $tmp_barcode=$line[0];
        substr($tmp_barcode,$k,1)="G";
        $all_barcode{$tmp_barcode}=$line[1];
        $tmp_barcode=$line[0];
        substr($tmp_barcode,$k,1)="C";
        $all_barcode{$tmp_barcode}=$line[1];
    }
    last;
}
close IN;

#
# Step 2.
#       0. detect barcodes type from r2
print "Step 2 : detect barcodes type from r2 ... \n";
if($ARGV[1]=~/\.gz/){
    open IN1,"gzip -dc $ARGV[1] |" or die "cannot open r1 file";
}else{
    open IN1,"$ARGV[1]" or die "cannot open r1 file";
}
open OUT1,"| gzip > $ARGV[1].split.1.fq.gz" or die "Can't open $ARGV[3].1.fq.gz for write";

if($ARGV[2]=~/\.gz/){
    open IN2,"gzip -dc $ARGV[2] |" or die "cannot open r2 file";
}else{
    open IN2,"$ARGV[2]" or die "cannot open r2 file";
}
open OUT2,"| gzip > $ARGV[2].split.1.fq.gz" or die "Can't open $ARGV[4].1.fq.gz for write";

my ($n1,$n2,$n3,$n4,$n5)=(10,6,10,0,10);
my $valid_read_len = 100;
my $r2_length = 0;
my $total_barcode_types = $index * $index *$index;
print "Barcode_types: $index X $index X $index = $total_barcode_types\n";
my $flag=0;
my $line_num = 0;
my $total_reads_num = 0;
my $valid_barcode_num = 0;
my $progress = 0 ;

print "Step 3 : parse barcodes .... \n";

my %new_barcode;
my %new_barcode_freq;
$new_barcode{"0_0_0"}=0;
$new_barcode_freq{"0_0_0"}=0;
my ($r1_1,$r1_2,$r1_3,$r1_4);
my ($r2_1,$r2_2,$r2_3,$r2_4);
while(<IN2>){
    chomp($r1_1=<IN1>);
    chomp($r1_2=<IN1>);
    chomp($r1_3=<IN1>);
    chomp($r1_4=<IN1>);
    chomp($r2_1=$_);
    chomp($r2_2=<IN2>);
    chomp($r2_3=<IN2>);
    chomp($r2_4=<IN2>);
#    chomp(my ($r2_1,$r2_2,$r2_3,$r2_4)=($_,<IN2>,<IN2>,<IN2>));
#    chomp(my ($r1_1,$r1_2,$r1_3,$r1_4)=(<IN1>,<IN1>,<IN1>,<IN1>));
    if($flag==0){
        $r2_length=length$r2_2;
        if( $r2_length == 154 ){
            $n4 = 18 ;
        }elsif ( $r2_length == 142 ){
            $n4 = 6 ;
        }elsif ( $r2_length == 126 ){
            $n4 = 0 ;
            $n5 = 0 ;
        }else{
            print "Unknow read length of r2:$r2_length.\nPlease check the $ARGV[2] file!\n";
            exit(1);
        }
        $flag=1;
    }

    my ($r1_id,$r2_id);
    $r1_id=$1 if($r1_1=~/^(.+)\/1/);
    $r2_id=$1 if($r2_1=~/^(.+)\/2/);
    die "Error: r1 r2 not match" if($r1_id ne $r2_id);
    $total_reads_num ++;
    $line_num += 4;
# print process ...
    if($line_num % 4000000 == 1){
        print "parse barcodes processed $progress (M) reads ...\n";
        $|=1;
        $progress ++ ;
    }

# check barcodes 
    my $b1 = substr($r2_2, $valid_read_len, $n1);
    my $b2 = substr($r2_2, $valid_read_len+$n1+$n2, $n3);
    my $b3 = substr($r2_2, $valid_read_len+$n1+$n2+$n3+$n4, $n5) if ($n5 !=0);
    my $r2_true_seq=substr($r2_2,0,$valid_read_len);
    my $r2_true_qua=substr($r2_4,0,$valid_read_len);
    my $str;
    if(exists $all_barcode{$b1} && exists $all_barcode{$b2} && exists $all_barcode{$b3}){
        $str = $all_barcode{$b1}."_".$all_barcode{$b2}."_".$all_barcode{$b3};
        $valid_barcode_num++;
    }else{
        $str = "0_0_0";
    }
    if(exists $new_barcode{$str}){
        $new_barcode_freq{$str}++;
    }else{
        $new_barcode{$str}=$flag;
        $new_barcode_freq{$str}=1;
        $flag++;
    }
    print OUT1 "$r1_id#$str\/1\t$new_barcode{$str}\t1\n";
    print OUT1 "$r1_2\n";
    print OUT1 "$r1_3\n";
    print OUT1 "$r1_4\n";
    print OUT2 "$r2_id#$str\/2\t$new_barcode{$str}\t1\n";
    print OUT2 "$r2_true_seq\n";
    print OUT2 "$r2_3\n";
    print OUT2 "$r2_true_qua\n";
}
close IN1;
close IN2;
close OUT1;
close OUT2;

# print stat
my $stat1 = 0 ;
my $stat2 = 0 ;
$stat1 = 100 * $valid_barcode_num/$total_barcode_types;
$stat2 = 100 * $valid_barcode_num/$total_reads_num;
open LOG, ">split_read_stat.log" or die "Can't open split_stat_read.log for write\n";
print LOG "Total_Barcode_num = $index X $index X $index = $total_barcode_types\n";
print LOG "Valid_Barcode_num = $valid_barcode_num($stat1% of all barcode types)\n";
print LOG "Total_reads_pair_num = $total_reads_num \n";
print LOG "Valid_reads_pair_num = $valid_barcode_num ($stat2 % of all)\n";
print LOG "read2_read_len : $r2_length \n";
print LOG "True_Barcode_type : $valid_read_len-$n1-$n2-$n3-$n4-$n5 \n";
close LOG;
#print barcode details 
open FREQ,">barcode_freq.txt" or die "Can't open barcode_freq.txt!\n";
print FREQ "Barcode_seq\tBarcode_count\tBarcode_num\n";
$flag=1;
foreach my $key(sort keys%new_barcode){
    print FREQ "$key\t$new_barcode_freq{$key}\t$flag\n";
    $flag++;
}
close FREQ;
print "all jobs done!\n";
$|=1;

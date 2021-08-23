#!/usr/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;
my ($inFile1,$inFile2,$field1,$field2,$out1,$tab,$suppress,%hash);
GetOptions(
            'i1|in1=s'          => \$inFile1,
	    'i2|in2=s'          => \$inFile2,
	    'f1|field1=i'       => \$field1,
            'f2|field2=i'       => \$field2,
            'o1|out1'           => \$out1,
			't|tab'             => \$tab,
			'v|supp'            => \$suppress,
			'h|help'            => sub{usage()}
	  ) || usage();
if(! defined $field1){
    $field1=1;
}
if(! defined $field2){
    $field2=1;
}
my ($IN1,$IN2);
if(defined $inFile1){
    open $IN1, "$inFile1" or die "Can't open file $inFile1:$!"; 
    }else{
    $IN1=\*STDIN;
}
if(defined $inFile2){
    open $IN2, "$inFile2" or die "Can't open file $inFile2:$!"; 
    }else{
    $IN2=\*STDIN;
}
while(<$IN1>){
    chomp;
    next if /^#/;
	my @split;
    if(defined $tab){
		@split=split /\t/;
	}else{
		@split=split /\s+/;
	}
    if(exists $hash{$split[$field1-1]}){
        $hash{$split[$field1-1]}="$hash{$split[$field1-1]}".";;;"."$_";
    }else{
        $hash{$split[$field1-1]}=$_;
    }
}
while(<$IN2>){
    chomp;
    next if /^#/;
    my @split;
	if(defined $tab){
		@split=split /\t/;
	}else{
		@split=split /\s+/;
	}
	if(defined $suppress){
		if(exists $hash{$split[$field2-1]}){
			$hash{$split[$field2-1]}="remove";
		}
	}else{
		if(exists $hash{$split[$field2-1]}){
			if(defined $out1){
				if($hash{$split[$field2-1]}=~/;;;/){
					my @fields=split /;;;/,$hash{$split[$field2-1]};
					for(my $i=0;$i<=$#fields;$i++){
						say "$fields[$i]";
					}
				}else{
					say "$hash{$split[$field2-1]}";
				}
			}else{
				if($hash{$split[$field2-1]}=~/;;;/){
					my @fields=split /;;;/,$hash{$split[$field2-1]};
					for(my $i=0;$i<=$#fields;$i++){
						say "$fields[$i]\t$_";
					}
				}else{
					say "$hash{$split[$field2-1]}\t$_";
				}
			}
		}
	}
}
if(defined $suppress){
	foreach my $key(keys %hash){
		if($hash{$key} ne "remove"){
			if($hash{$key}=~/;;;/){
				my @fields=split /;;;/,$hash{$key};
				for(my $i=0;$i<=$#fields;$i++){
					say "$fields[$i]";
				}
			}else{
				say "$hash{$key}";
			}
		}
	}
}
sub usage{
print <<HELP;
Usage: perl $0 -i1 file1 -i2 file2 -f1 2 -f2 3  >result
Decription: For each pair of input lines with identical join fields(the field can be used as hash key), write a line to standard output. The default join field is the first column.If you don't set file1 or file2, then it will read from STDIN.
Revision 1: Change separator ';' to ';;;' in case there is ; in the input file;
Revision 2: Add -v parameter. 2018/8/28.
Output:file1 file2
Options:
    -i1 FILE           The first file to join(seperated by whitespace or tab).
    -i2 FILE           The second file to join(seperated by whitespace or tab).
    -f1 INT            Join on this FIELD of file 1[default 1].
    -f2 INT            Join on this FIELD of file 2[default 1].
    -t	LOGIC		   Input files are tab-delimited.
	-v  LOGIC          Suppress the joined lines of file1.
	-o1 LOGIC          Only output file1.
    -h --help          Print this help information.
HELP
    exit(-1);
}

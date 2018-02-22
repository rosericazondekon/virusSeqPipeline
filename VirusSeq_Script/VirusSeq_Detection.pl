#! /usr/bin/perl -w

#my $sample_name_file=$ARGV[0]; ##filename_virusLog from MosaikSort
#my $min_readCount=$ARGV[1]; #the empirical cutoff for virus existence=1000
#my $output_filename=$ARGV[2]; #define by user

VirusSeq_Detection();
process_virus_profiling();

sub VirusSeq_Detection {

 
    $file_list[1]=$ARGV[0];
    $file_number=1;
    @virus_data=();
    @virus_data_row_number=();


    for (my $i=1; $i<=$file_number; $i++) {
	my $virus_indicator=0;
	my $read_indicator=0;
	$virus_data_row_number[$i]=0;
	my $in_data_file=$file_list[$i];
	open (INPUT, "<$in_data_file")
	    || die "Can't open $in_data_file $!";

	while (<INPUT>) {
	    $text=$_;
	    $text=~s/\n//g;
 
	    if ($text=~/Processing/) {
		last;
	    }
	    if ($virus_indicator==0) { $text=~s/ //g; }

	    if ($text=~/alignmentcountreferencesequence/) {
		$virus_indicator=1;
	    }
	    if (($virus_indicator==1) && (length($text)>5) && (!($text=~/chr/)) &&
		(!($text=~/\-\-/)) && (!($text=~/alignmentcountreferencesequence/))) {

		$virus_data_row_number[$i]++;
		my (@temp_data)=split(" ", $text);
		$virus_data[$i][$virus_data_row_number[$i]][1]=$temp_data[1];
		$virus_data[$i][$virus_data_row_number[$i]][2]=$temp_data[0];
	    }
	    
	}
	close(INPUT);
	$virus_data_row_number[$i]--;
    }

    @aggregate_data=();
    $aggregate_data_row_number=0;
    $aggregate_data_column_number=$file_number;
    %virus_index=();
    ##get the overall number of virus
    for (my $i=1; $i<=$file_number; $i++) {
	for (my $j=1; $j<=$virus_data_row_number[$i]; $j++) {
	    if (!(exists($virus_index{$virus_data[$i][$j][1]}))) {
		$aggregate_data_row_number++;
		$virus_index{$virus_data[$i][$j][1]}=$aggregate_data_row_number;
		$aggregate_data[$aggregate_data_row_number][0]=$virus_data[$i][$j][1];
	    }
	}
    }
    #initialize the aggregate data
    $aggregate_data[0][0]="VirusName";
    for (my $j=1; $j<=$file_number; $j++) {
	$aggregate_data[0][$j]=$file_list[$j];
    }
    for (my $i=1; $i<=$aggregate_data_row_number; $i++) {
	for (my $j=1; $j<=$file_number; $j++) {
	    $aggregate_data[$i][$j]=0;
	}
    }
    #assign the value into aggregate data
    for (my $i=1; $i<=$file_number; $i++) {
	for (my $j=1; $j<=$virus_data_row_number[$i]; $j++) {
	    if (exists($virus_index{$virus_data[$i][$j][1]})) {
		my $temp_index=$virus_index{$virus_data[$i][$j][1]};
		$aggregate_data[$temp_index][$i]=$virus_data[$i][$j][2];
	    }
	}
    }

    ##detect common number
    ##the read number has to be larger than 1000 to be counted
    $min_readCount=$ARGV[1];
    $aggregate_data[0][$file_number+1]="SampleNo_withVirus";
    for (my $i=1; $i<=$aggregate_data_row_number; $i++) {
	my $temp_comm=0;
	for (my $j=1; $j<=$aggregate_data_column_number; $j++) {
	    if ($aggregate_data[$i][$j]>$min_readCount) {
		$temp_comm++;
	    }
	    $aggregate_data[$i][$file_number+1]=$temp_comm;
	}
    }


    #OUTPUT THE DATA
    #my $out_data_file=$ARGV[2];
    #open (OUTPUT, ">$out_data_file")
	#|| die "Can't open $out_data_file $!";

    #for (my $i=0; $i<=0; $i++) {
	#for (my $j=0; $j<=$aggregate_data_column_number; $j++) {
	 #   print OUTPUT "$aggregate_data[$i][$j]\t";
	#}
	#print OUTPUT "$aggregate_data[$i][$aggregate_data_column_number+1]\n";
    #}

   #OUTPUT THE DATA
    for (my $i=1; $i<=$aggregate_data_row_number; $i++) {
	if ((!($aggregate_data[$i][0]=~/J02482/)) &&
	    (!($aggregate_data[$i][0]=~/HepatitisC/)) &&
	    (!($aggregate_data[$i][0]=~/Entero/)) &&
	    (!($aggregate_data[$i][0]=~/Coliphage/)) &&
	    (!($aggregate_data[$i][0]=~/Genomeofphage/)) &&
	    (!($aggregate_data[$i][0]=~/Citrus/)) &&
	    (!($aggregate_data[$i][0]=~/Cercopithecine/)) &&
	    (!($aggregate_data[$i][0]=~/AY653733/)) &&
	    (!($aggregate_data[$i][0]=~/Choristoneura/)) &&
	    (!($aggregate_data[$i][0]=~/Vaccinia/)) &&
	    (!($aggregate_data[$i][0]=~/Acanthamoeba/)) &&
	    (!($aggregate_data[$i][0]=~/Cafeteria/)) &&
	    (!($aggregate_data[$i][0]=~/Phage/)) &&
	    (!($aggregate_data[$i][0]=~/Spring/)) &&
	    (!($aggregate_data[$i][0]=~/Lacto/)) &&
	    (!($aggregate_data[$i][0]=~/Bacteriophage/))) {
	    my $virus_indicator=0;
	    for (my $j=1; $j<=$aggregate_data_column_number; $j++) {
		if ($aggregate_data[$i][$j]>=$min_readCount) {
		    $virus_indicator=1;
		}
	    }
	    if ($virus_indicator==1) {
		#print OUTPUT "$aggregate_data[$i][0]\t";
		for (my $j=1; $j<=$aggregate_data_column_number; $j++) {
		    if ($aggregate_data[$i][$j]<$min_readCount) {
			$aggregate_data[$i][$j]=0;
		    }
		    #printf OUTPUT "%0.0f\t", $aggregate_data[$i][$j];
		    #print OUTPUT "$aggregate_data[$i][$j]\t";
		}
		#print OUTPUT "$aggregate_data[$i][$aggregate_data_column_number+1]\n";
	    }
	}    
    }
    #close(OUTPUT);
}

sub process_virus_profiling {
    #$text=~s/|\n//g;
    @sampleID_virus=();
    for (my $j=1; $j<=$aggregate_data_column_number; $j++) {
	$sampleID_virus[$j][1]=$aggregate_data[0][$j];
	$sampleID_virus[$j][2]="NA";
	$sampleID_virus[$j][3]=0;
	for (my $i=1; $i<=$aggregate_data_row_number; $i++) {
	    if ((!($aggregate_data[$i][0]=~/J02482/)) &&
		(!($aggregate_data[$i][0]=~/HepatitisC/)) &&
		(!($aggregate_data[$i][0]=~/Entero/)) &&
		(!($aggregate_data[$i][0]=~/Coliphage/)) &&
		(!($aggregate_data[$i][0]=~/Genomeofphage/)) &&
		(!($aggregate_data[$i][0]=~/Citrus/)) &&
		(!($aggregate_data[$i][0]=~/Cercopithecine/)) &&
		(!($aggregate_data[$i][0]=~/AY653733/)) &&
		(!($aggregate_data[$i][0]=~/Choristoneura/)) &&
		(!($aggregate_data[$i][0]=~/Vaccinia/)) &&
		(!($aggregate_data[$i][0]=~/Acanthamoeba/)) &&
		(!($aggregate_data[$i][0]=~/Cafeteria/)) &&
		(!($aggregate_data[$i][0]=~/Phage/)) &&
		(!($aggregate_data[$i][0]=~/Spring/)) &&
		(!($aggregate_data[$i][0]=~/Lacto/)) &&
		(!($aggregate_data[$i][0]=~/Bacteriophage/)) &&
		($aggregate_data[$i][$j]>=$min_readCount) &&
		($aggregate_data[$i][$j]>$sampleID_virus[$j][3])) {
		
		$sampleID_virus[$j][2]=$aggregate_data[$i][0];
		$sampleID_virus[$j][3]=$aggregate_data[$i][$j];
	    }		
	}
    }       
    
    #output HPV virus data
    
    my $out_data_file=$ARGV[2];
    open (OUTPUT, ">$out_data_file")
	|| die "Can't open $out_data_file $!";

    print OUTPUT "SampleID\tVirusType\tMappedReads\n";
    for (my $i=1; $i<=$aggregate_data_column_number; $i++) {
	if ($sampleID_virus[$i][3]>=$min_readCount) {
	    print OUTPUT "$sampleID_virus[$i][1]\t$sampleID_virus[$i][2]\t$sampleID_virus[$i][3]\n";
	}
    }
    close(OUTPUT);
}

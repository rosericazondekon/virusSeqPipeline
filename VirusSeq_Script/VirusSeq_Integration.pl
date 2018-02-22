#! /usr/bin/perl

####! /usr/bin/perl -w
# my $in_data_file=$ARGV[0]; crossChip file
# my $in_data_file=$ARGV[1]; #refSeq annotation file=hg19_refGene_RIS
# my $genome_version=$ARGV[2]; #version of genome=hg19 or hg18
# my $fragment_length=$ARGV[3]; #fragment length
# my $fragment_stdev=$ARGV[4]; #fragment stdev
# my $PEread_length=$ARGV[5]; #=50; #this one is for PE read length
# my $final_output_filename=$ARGV[6]; #Virus integration output filename (named by user)

#hgGenome refGene location
#start to deal with UCSC refGene data
get_hgGenome_refGene();
regenerate_hgGenome_refGene_from_fiveprime_to_threeprime();

sort_hgGenome_refGene_by_geneName();
#remove the redundant gene symbols
remove_redundant_genesymbols_byExonCount_GeneSize();

get_chromosome_genomic_location_hgGenome_refGene();
generate_hgGenome_refGene_Pointer_by_chrom();

###generate the gene relative location();
generate_cDNA_order();

#start to deal with cross data
get_reduced_cross_data();

#add exon info to each gene (put them together)
######add_exon_location_to_gene_for_each_read();


##re-sort the data by anchor0_anchor1
filter_reduced_cross_data();
resort_reduced_crossdata_by_reduced_anchor1();
in_silico_clustering_generator();
################output_reduced_cross_data();

#Count the read number for each gene (cluster)
count_reads_for_each_gene();

#claculate the in-silico fragment length
get_Spanner_stats_data();
in_silico_fragment_by_anchor0_anchor1();
filter_reduced_cross_data_by_insilicoFrag();

#start blat search
replace_blat_parsing();
#output_fasta_for_blat_from_reduced_annotate_cross_data();
#call_blat_alignment();
#get_blat_data();

##for reduced cross data
build_gene_exon_pair_reduced_cross_data();
#remove again
remove_exon_from_gene_reduced_cross_data();

add_exon_location_for_each_read();

#add exon back to gene
add_exon_backto_gene_reduced_cross_data();

#generate the consensus sequence
#output_reduced_annotate_cross_data_after_blat();

#aggregtae the reduced data
aggregate_blat_crossdata_by_anchor0_anchor1_BlatHit();

#output the result
#aggregate blat cross data
build_gene_exon_pair_aggregate_blat_cross_data();
remove_exon_from_gene_aggregate_blat_cross_data();
assess_orientation_both_reads_refGenes();

#add in-silico sequence for predicted fusion
extract_sequence_for_predicted_fusions();

#output the final table
add_exon_backto_gene_aggregate_blat_cross_data();
output_aggregate_reduced_annotate_cross_data_after_blat();
######################generate_predicted_fusion_reference_fasta();

#output the fasta file for bed file generation by Mosaik
#output_fasta_for_MosaikBed_from_aggregate_reduced_annotate_cross_data_after_blat();

sub get_reduced_cross_data {
    my $in_data_file=$ARGV[0];
    open (INPUT, "<$in_data_file")
	|| die "Can't open $in_data_file $!";
    
    
    @reduced_cross_data=();
    
    $reduced_cross_row_number=0; #genes number
    #input data from files
    
    while (<INPUT>) {
	$text=$_;
	$text=~s/|\n//g;
	
	(@line)=split("\t", $text);
	
	$reduced_cross_column_number=0;
	for $read_data(@line) {
	    # $read_data=~s/ //g;
	    $reduced_cross_data[$reduced_cross_row_number][$reduced_cross_column_number]=$read_data;
	    $reduced_cross_column_number++;
	}
	
	$reduced_cross_row_number++; 
	
   } #while(INPUT)
    
    close(INPUT);
    
    $reduced_cross_row_number--;
    $reduced_cross_column_number--;
    
    #recover the original column number from Spanner
    $reduced_cross_column_number=$reduced_cross_column_number-2;
    
    #output all data
    
    #for (my $i=0; $i<=$reduced_cross_row_number; $i++) {
    # for (my $j=0; $j<=$reduced_cross_column_number; $j++) {
    #   print "$reduced_cross_data[$i][$j]	";
    #}
    # print "\n";
    #}
    
}

sub count_reads_for_each_gene {
    %gene_reads_count=();
    
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	if (!(exists($gene_reads_count{$reduced_cross_data[$i][$reduced_cross_column_number+1]}))) {
	    $gene_reads_count{$reduced_cross_data[$i][$reduced_cross_column_number+1]}=1;
	}
	else {
	    $gene_reads_count{$reduced_cross_data[$i][$reduced_cross_column_number+1]}++;
	}
	
	if (!(exists($gene_reads_count{$reduced_cross_data[$i][$reduced_cross_column_number+2]}))) {
	    $gene_reads_count{$reduced_cross_data[$i][$reduced_cross_column_number+2]}=1;
	}
	else {
	    $gene_reads_count{$reduced_cross_data[$i][$reduced_cross_column_number+2]}++;
	}
   }
}

sub get_cDNA_coverage {
    my $in_data_file=$ARGV[1];
    open (INPUT, "<$in_data_file")
	 || die "Can't open $in_data_file $!";
    
    
    @cDNA_coverage=();
    
    $cDNA_coverage_row_number=0; #genes number
    $line_indicator=0;
    $first_line=1;
    
    #input data from files
    
    while (<INPUT>) {
	$text=$_;
	$text=~s/|\n//g;
	
	(@line)=split("\t", $text);
	
	$cDNA_coverage_column_number=0;
	for $read_data(@line) {
	    # $read_data=~s/ //g;
	    $cDNA_coverage[$cDNA_coverage_row_number][$cDNA_coverage_column_number]=$read_data;
	    $cDNA_coverage_column_number++;
	}
	
	# detect if there is any missing data point in a specific row
	if ($first_line==1) {
	    $cDNA_coverage_column_number_previous=$cDNA_coverage_column_number;
	}
	else {
	    if ($cDNA_coverage_column_number_previous!=$cDNA_coverage_column_number) {
		$line_indicator=1;
	   }
	    $cDNA_coverage_column_number_previous=$cDNA_coverage_column_number;
	}
	
	if ($line_indicator==1) {
	    
	    #print "there is missing data point in line= $cDNA_coverage_row_number\n";
	    $line_indicator=0;
	}
	# end of detection of missing data point
	
	#if (!($cDNA_coverage[$cDNA_coverage_row_number][3]=~/\_/)) { 
	$cDNA_coverage_row_number++; 
	#}
	
	#$cDNA_coverage_row_number++;
	$first_line=0;
	
   } #while(INPUT)
    
    close(INPUT);
    
    $cDNA_coverage_row_number--;
    $cDNA_coverage_column_number--;
    
    #output all data
    
    #%cDNA_depth=();
    #for (my $i=1; $i<=$cDNA_coverage_row_number; $i++) {
    #$cDNA_depth{$cDNA_coverage[$i][7]}=$cDNA_coverage[$i][$cDNA_coverage_column_number];
    #}
    #print "Gene_Number=$cDNA_coverage_row_number\n";
    #for (my $i=0; $i<=$cDNA_coverage_row_number; $i++) {
    # for (my $j=0; $j<=$cDNA_coverage_column_number; $j++) {
    #  print "$cDNA_coverage[$i][$j]\t";
    #}
    #print "\n";
    #}
    
}

sub replace_blat_parsing {

    #initialize blat hit
    $reduced_cross_data[0][$reduced_cross_column_number+3]="Read0_BlatHit";
    $reduced_cross_data[0][$reduced_cross_column_number+4]="Read1_BlatHit";
    $reduced_cross_data[0][$reduced_cross_column_number+5]="BlatCall";
    $reduced_cross_data[0][$reduced_cross_column_number+6]="Read0_Loc";
    $reduced_cross_data[0][$reduced_cross_column_number+7]="Read1_Loc";
    
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	$reduced_cross_data[$i][$reduced_cross_column_number+3]=0;
	$reduced_cross_data[$i][$reduced_cross_column_number+4]=0;
	$reduced_cross_data[$i][$reduced_cross_column_number+5]="Yes";
	$reduced_cross_data[$i][$reduced_cross_column_number+6]="NA";	
	$reduced_cross_data[$i][$reduced_cross_column_number+7]="NA";	
    }

}

sub output_fasta_for_blat_from_reduced_annotate_cross_data {
    #output the data
    #my $temp_fa_file="temp_blat.fa";
    my $temp_fa_file=$ARGV[0];
    #$temp_fa_file=~s/Chip/Blat/;
    $temp_fa_file=~s/\.txt/\.fa/;
    open (OUTPUT, ">$temp_fa_file")
	|| die "Can't open $temp_fa_file $!";
    
    #initialize blat hit
    $reduced_cross_data[0][$reduced_cross_column_number+3]="Read0_BlatHit";
    $reduced_cross_data[0][$reduced_cross_column_number+4]="Read1_BlatHit";
    $reduced_cross_data[0][$reduced_cross_column_number+5]="BlatCall";
    
    $reduced_cross_data[0][$reduced_cross_column_number+6]="Read0_Loc";
    $reduced_cross_data[0][$reduced_cross_column_number+7]="Read1_Loc";
    
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	$reduced_cross_data[$i][$reduced_cross_column_number+3]=0;
	$reduced_cross_data[$i][$reduced_cross_column_number+4]=0;
	$reduced_cross_data[$i][$reduced_cross_column_number+5]="No";	
	$reduced_cross_data[$i][$reduced_cross_column_number+6]="NA"; #used for duplicate	
	$reduced_cross_data[$i][$reduced_cross_column_number+7]="NA"; #used for duplicate	
   }
    
    %readsID_index_gene0=();
    %readsID_index_gene1=();
    
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	if (($reduced_cross_data[$i][3]=~/\>/) &&
	    ($reduced_cross_data[$i][4]=~/A|C|G|T|a|c|g|t/)) {
	    my $temp_ID=$reduced_cross_data[$i][3];
	    $temp_ID=~s/\>//;
	    $readsID_index_gene0{$temp_ID}=$i;
	    print OUTPUT "$reduced_cross_data[$i][3]\n";
	    print OUTPUT "$reduced_cross_data[$i][4]\n";
	}

	if (($reduced_cross_data[$i][8]=~/\>/) && 
	    ($reduced_cross_data[$i][9]=~/A|C|G|T|a|c|g|t/)) {
	    my $temp_ID=$reduced_cross_data[$i][8];
	    $temp_ID=~s/\>//;
	    $readsID_index_gene1{$temp_ID}=$i;
	    print OUTPUT "$reduced_cross_data[$i][8]\n";
	    print OUTPUT "$reduced_cross_data[$i][9]\n";
	}

   }
    
    close(OUTPUT);
}

sub call_blat_alignment {
    my $temp_fa_file=$ARGV[0];
    my $minScore=30; #=30; #this one is for Blat search sensitivity
    #$temp_fa_file=~s/Chip/Blat/;
    $temp_fa_file=~s/\.txt/\.fa/;
    my $temp_psl_file=$temp_fa_file;
    $temp_psl_file=~s/\.fa/\.psl/;
    #system("blat ../../RefSeq/hgGenome/hgGenome.2bit temp_blat.fa temp_blat.psl");
    my $genome_version=$ARGV[2];
    my $data_path_file="\/RIS\/home\/xsu1\/Mosaik\/hg19_Mosaik\/".$genome_version."\_2bit"."\/".$genome_version."OV\.2bit";
    #system("blat /home2/xsu/Mosaik/RefSeq/hgGenome/hgGenome.2bit $temp_fa_file $temp_psl_file");
    system("/RIS/home/xsu1/bin/blat -minScore=$minScore -minIdentity=90 $data_path_file $temp_fa_file $temp_psl_file");    
}

sub get_blat_data {
    #my $blat_file="temp_blat.psl";
    my $temp_fa_file=$ARGV[0];
    #$temp_fa_file=~s/Chip/Blat/;
    $temp_fa_file=~s/\.txt/\.fa/;
    my $blat_file=$temp_fa_file;
    $blat_file=~s/\.fa/\.psl/;
    #my $blat_file="temp_blat.psl";
    open (INPUT, "<$blat_file")
	|| die "Can't open $blat_file $!";

    my $start_line=0;
    my $refGene_RNA_Genome=$ARGV[1]; #=hg19_refGene_RIS.txt or hg19_TruecDNA_refGene_RIS.txt
    while (<INPUT>) {
	$text=$_;
	#$text=~s/|\n//g;
	
	 $start_line++;

	 if ($start_line>=6) {
	     (@line)=split("\t", $text);
	     my $temp_data=();
	     my $blat_column_number=0;
	     for $read_data(@line) {
		 #$read_data=~s/\"//g;
		 $temp_data[$blat_column_number]=$read_data;
		 $blat_column_number++;
	    }

	     #detect the blat hit number
	     ##[9]=readID, [13]=chrnumber [15]=t_start [16]=t_end
	     #if ((exists($readsID_index_gene0{$temp_data[9]})) && ($temp_data[7]<10)) {
	     #transcriptome mapping
	     if ($refGene_RNA_Genome=~/cDNA/) {
		 if (exists($readsID_index_gene0{$temp_data[9]})) {
		     my $gene0_index=$readsID_index_gene0{$temp_data[9]};
		     $reduced_cross_data[$gene0_index][$reduced_cross_column_number+3]++;
		 }
		 #if ((exists($readsID_index_gene1{$temp_data[9]})) && ($temp_data[7]<10)) {
		 if (exists($readsID_index_gene1{$temp_data[9]})) {
		     my $gene1_index=$readsID_index_gene1{$temp_data[9]};
		     $reduced_cross_data[$gene1_index][$reduced_cross_column_number+4]++;
		 }
	     }
	     #genome mapping
	     else {
		 if (exists($readsID_index_gene0{$temp_data[9]})) {
		     my $gene0_index=$readsID_index_gene0{$temp_data[9]};
		     $reduced_cross_data[$gene0_index][$reduced_cross_column_number+3]++;
		     #annotate the mapping location with geneID
		     my (@temp_partner0)=split("\/\/", $reduced_cross_data[$gene0_index][$reduced_cross_column_number+1]);
		     my (@temp_partner1)=split("\/\/", $reduced_cross_data[$gene0_index][$reduced_cross_column_number+2]);
		     if ($temp_partner0[0] ne $temp_partner1[0]) {
			 if ((!($temp_data[13]=~/ran/)) && (!($temp_data[13]=~/hap/)) && 
			     (!($temp_data[13]=~/Un/))) {
			     #(!($temp_data[13]=~/M/)) && (!($temp_data[13]=~/Un/))) {
			     $temp_data[13]=~s/X/23/;
			     $temp_data[13]=~s/Y/24/;
			     $temp_data[13]=~s/M/25/;
			     $temp_data[13]=~s/chr//;
			     my $geneID=annotate_blat_data_with_refGene($temp_data[13], $temp_data[15]);
			     #print "readID=$temp_data[9]  outside_geneID=$geneID\n";
			     #test if the geneID from first readID is in second geneID
			     my (@temp_fusion_partner)=split("\/\/", $reduced_cross_data[$gene0_index][$reduced_cross_column_number+2]);
			     if ($geneID eq $temp_fusion_partner[0]) {
				 #print "geneID=$geneID\n";
				 $reduced_cross_data[$gene0_index][$reduced_cross_column_number+6]="Duplicate";
			     }
			 }
		     }
		 }
		 #if ((exists($readsID_index_gene1{$temp_data[9]})) && ($temp_data[7]<10)) {
		 if (exists($readsID_index_gene1{$temp_data[9]})) {
		     my $gene1_index=$readsID_index_gene1{$temp_data[9]};
		     $reduced_cross_data[$gene1_index][$reduced_cross_column_number+4]++;
		     #annotate the mapping location with geneID
		     my (@temp_partner0)=split("\/\/", $reduced_cross_data[$gene1_index][$reduced_cross_column_number+1]);
		     my (@temp_partner1)=split("\/\/", $reduced_cross_data[$gene1_index][$reduced_cross_column_number+2]);
		     if ($temp_partner0[0] ne $temp_partner1[0]) {
			 if ((!($temp_data[13]=~/ran/)) && (!($temp_data[13]=~/hap/)) &&
			     (!($temp_data[13]=~/Un/))) {
			     #(!($temp_data[13]=~/M/)) && (!($temp_data[13]=~/Un/))) {
			     $temp_data[13]=~s/X/23/;
			     $temp_data[13]=~s/Y/24/;
			     $temp_data[13]=~s/M/25/;
			     $temp_data[13]=~s/chr//;
			     my $geneID=annotate_blat_data_with_refGene($temp_data[13], $temp_data[15]);
			     #print "readID=$temp_data[9]  outside_geneID=$geneID\n";
			     #test if the geneID from second readID is in first geneID
			     my (@temp_fusion_partner)=split("\/\/", $reduced_cross_data[$gene1_index][$reduced_cross_column_number+1]);
			     if ($geneID eq $temp_fusion_partner[0]) {
				 $reduced_cross_data[$gene1_index][$reduced_cross_column_number+7]="Duplicate";
			     }
			 }
		     }
		 }
	     } #for genome ampping
	}

    } #while(INPUT)

     close(INPUT);

     #finalize the consensus blat hit
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 if (($reduced_cross_data[$i][$reduced_cross_column_number+3]==1) &&
	     ($reduced_cross_data[$i][$reduced_cross_column_number+4]==1)) {
	     $reduced_cross_data[$i][$reduced_cross_column_number+5]="Yes";
	}
    }

     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 if (($reduced_cross_data[$i][$reduced_cross_column_number+6]=~/Duplicate/) ||
	     ($reduced_cross_data[$i][$reduced_cross_column_number+7]=~/Duplicate/)) {
	     $reduced_cross_data[$i][$reduced_cross_column_number+5]="Duplicate";
	}
    }
}


sub annotate_blat_data_with_refGene {

    my ($chromosome, $start_location)=@_;
    my $hgGenome_start_pointer=$hgGenome_refGene_pointer[$chromosome][1];
     my $hgGenome_end_pointer=$hgGenome_refGene_pointer[$chromosome][2];

     my $final_down_index=detect_start_location_hgGenome($start_location,
						     $hgGenome_start_pointer, $hgGenome_end_pointer);

     my $temp_geneID="Unkn".$chromosome;

     if (($start_location>=$hgGenome_refGene[$final_down_index][5]) && 
	 ($start_location<=$hgGenome_refGene[$final_down_index][6])) {
	 $temp_geneID=$hgGenome_refGene[$final_down_index][7];
    }

     #move up the pointer by 5
     if ($temp_geneID=~/Unkn/) {
	 my $temp_pointer=$final_down_index-1;
	 while ($temp_pointer>=($final_down_index-5)) {
	     if ($temp_pointer<=1) {
		 last;
		 $temp_pointer=2;
	    }
	     if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		 ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		 $temp_geneID=$hgGenome_refGene[$temp_pointer][7];
		 last;
	    }
	     #move the pointer
	     $temp_pointer--;
	}  
    }

     #move down the pointer by 5
     if ($temp_geneID=~/Unkn/) {
	 my $temp_pointer=$final_down_index+1;
	 while ($temp_pointer<=($final_down_index+5)) {
	     if ($temp_pointer>=$hgGenome_refGene_row_number) {
		 last;
		 $temp_pointer--;
	    }

	     if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		 ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		 $temp_geneID=$hgGenome_refGene[$temp_pointer][7];
		 last;
	    }
	     #move the pointer
	     $temp_pointer++;
	}  
    }

     #assign a nearby gene
     if ($temp_geneID=~/Unkn/) {
	 if ($final_down_index<1) {
	     $final_down_index=1;
	}
	 if ($final_down_index>$hgGenome_refGene_row_number) {
	     $final_down_index=$hgGenome_refGene_row_number;
	}	
	 $temp_geneID=$hgGenome_refGene[$final_down_index][7];
    }

     return($temp_geneID);
}


sub detect_start_location_hgGenome {
    my ($start_location, $hgGenome_start_pointer, $hgGenome_end_pointer)=@_;

     #detect the start location of hgGenome
     my $moving_start_pointer=$hgGenome_start_pointer;
     my $moving_end_pointer=$hgGenome_end_pointer;	
     my $moving_mid_pointer=$moving_start_pointer+
	 int(($moving_end_pointer-$moving_start_pointer+1)/2);

     while (($moving_mid_pointer>$moving_start_pointer) && ($moving_mid_pointer<$moving_end_pointer)) {	    
	 if ($start_location<=$hgGenome_refGene[$moving_mid_pointer][5]) {
	     #print "moving_start_pointer=$moving_start_pointer\n";
	     #print "moving_mid_pointer=$moving_mid_pointer\n";
	     #print "moving_end_pointer=$moving_end_pointer\n";
	     #reset index
	     $moving_end_pointer=$moving_mid_pointer;
	     #reset the mid-index
	     $moving_mid_pointer=$moving_start_pointer+
		 int(($moving_end_pointer-$moving_start_pointer+1)/2);
	}
	 elsif ($start_location>$hgGenome_refGene[$moving_mid_pointer][5]) {
	     #reset index
	     $moving_start_pointer=$moving_mid_pointer;	  
	     #reset the mid-index
	     $moving_mid_pointer=$moving_start_pointer+
		 int(($moving_end_pointer-$moving_start_pointer+1)/2);
	}
    }

     return($moving_mid_pointer);

}	


sub get_blat_data_old {
    #my $blat_file="temp_blat.psl";
    my $temp_fa_file=$ARGV[0];
     #$temp_fa_file=~s/Chip/Blat/;
     $temp_fa_file=~s/\.txt/\.fa/;
     my $blat_file=$temp_fa_file;
     $blat_file=~s/\.fa/\.psl/;
     #my $blat_file="temp_blat.psl";
     open (INPUT, "<$blat_file")
	 || die "Can't open $blat_file $!";

     my $start_line=0;

     while (<INPUT>) {
	 $text=$_;
	 #$text=~s/|\n//g;

	 $start_line++;

	 if ($start_line>=6) {
	     (@line)=split("\t", $text);
	     my $temp_data=();
	     my $blat_column_number=0;
	     for $read_data(@line) {
		 #$read_data=~s/\"//g;
		 $temp_data[$blat_column_number]=$read_data;
		 $blat_column_number++;
	    }

	     #detect the blat hit number
	     #if ((exists($readsID_index_gene0{$temp_data[9]})) && ($temp_data[7]<10)) {
	     if (exists($readsID_index_gene0{$temp_data[9]})) {
		 my $gene0_index=$readsID_index_gene0{$temp_data[9]};
		 $reduced_cross_data[$gene0_index][$reduced_cross_column_number+3]++;
	    }
	     #if ((exists($readsID_index_gene1{$temp_data[9]})) && ($temp_data[7]<10)) {
	     if (exists($readsID_index_gene1{$temp_data[9]})) {
		 my $gene1_index=$readsID_index_gene1{$temp_data[9]};
		 $reduced_cross_data[$gene1_index][$reduced_cross_column_number+4]++;
	    }
	}

    } #while(INPUT)

     close(INPUT);

     #finalize the consensus blat hit
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 if (($reduced_cross_data[$i][$reduced_cross_column_number+3]==1) &&
	     ($reduced_cross_data[$i][$reduced_cross_column_number+4]==1)) {
	     $reduced_cross_data[$i][$reduced_cross_column_number+5]="Yes";
	}
    }

}

sub add_exon_location_for_each_read {
    #index hgGenome by gene name
    index_hgGenome_by_genename();
     my $aligned_read_length=60; #read_length=60bp

     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 my $gene0_name=$reduced_cross_data[$i][$reduced_cross_column_number+1];
	 if (exists($genename_to_hgGenome_index{$gene0_name})) {
	     my $gene0_index=$genename_to_hgGenome_index{$gene0_name};
	     my $gene0_location=$reduced_cross_data[$i][1];
	     my $first_location=detect_location_site_within_gene_byLocation($gene0_index, $gene0_location);
	     my $second_location=detect_location_site_within_gene_byLocation($gene0_index, $gene0_location+$aligned_read_length);
	     my $third_location=detect_location_site_within_gene_byLocation($gene0_index, $gene0_location-$aligned_read_length);
	     if ($first_location=~/exon/) {
		 $reduced_cross_data[$i][$reduced_cross_column_number+6]=$first_location;
	    }
	     elsif ($second_location=~/exon/) {
		 $reduced_cross_data[$i][$reduced_cross_column_number+6]=$second_location;
	    }
	     else {
		 $reduced_cross_data[$i][$reduced_cross_column_number+6]=$third_location;
	    }		

	}

	 my $gene1_name=$reduced_cross_data[$i][$reduced_cross_column_number+2];
	 if (exists($genename_to_hgGenome_index{$gene1_name})) {
	     my $gene1_index=$genename_to_hgGenome_index{$gene1_name};
	     my $gene1_location=$reduced_cross_data[$i][6];
	     my $first_location=detect_location_site_within_gene_byLocation($gene1_index, $gene1_location);
	     my $second_location=detect_location_site_within_gene_byLocation($gene1_index, $gene1_location+$aligned_read_length);
	     my $third_location=detect_location_site_within_gene_byLocation($gene1_index, $gene1_location-$aligned_read_length);
	     if ($first_location=~/exon/) {
		 $reduced_cross_data[$i][$reduced_cross_column_number+7]=$first_location;
	    }
	     elsif ($second_location=~/exon/) {
		 $reduced_cross_data[$i][$reduced_cross_column_number+7]=$second_location;
	    }
	     else {
		 $reduced_cross_data[$i][$reduced_cross_column_number+7]=$third_location;
	    }	
	}
    } 

}

sub index_hgGenome_by_genename {
    %genename_to_hgGenome_index=();
    for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
	$genename_to_hgGenome_index{$hgGenome_refGene[$i][7]}=$i;
   }
}

sub detect_location_site_within_gene_byLocation {
    my ($temp_index, $del_data_location)=@_;
    
    my $exonCount=$hgGenome_refGene[$temp_index][1];
     #start to deal exon location
     my $exonlocation=$hgGenome_refGene[$temp_index][2];
     $exonlocation=~s/\|\|//;
     my (@temp_exonlocation)=split(",", $exonlocation);

     #put the exon location into an array
     my @true_exonlocation=();
     for (my $i=1; $i<=$exonCount; $i++) {
	 $true_exonlocation[$i][1]=$temp_exonlocation[$i-1];
	 $true_exonlocation[$i][2]=$temp_exonlocation[$i-1+$exonCount];
    }

     #detect the location of del_data in gene
     my $temp_del_data_location=$del_data_location;
     my $location="NA";
     for (my $i=1; $i<=$exonCount; $i++) {
	 #print "genename=$hgGenome_refGene[$temp_index][7]\n";

	 if (($temp_del_data_location>=$true_exonlocation[$i][1]) &&
	     ($temp_del_data_location<=$true_exonlocation[$i][2])) {
	     if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		 $location="exon".$i;
	    }
	     else { #reverse
		 my $exon_order=$exonCount-$i+1;
		 $location="exon".$exon_order;
	    }
	     return($location);
	}
	 elsif (($i<$exonCount) && ($temp_del_data_location>$true_exonlocation[$i][2]) &&
		($temp_del_data_location<$true_exonlocation[$i+1][1])) {
	     if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		 $location="intron".$i;
	    }
	     else { #reverse
		 my $exon_order=$exonCount-$i;
		 $location="intron".$exon_order;
	    }		
	     return($location);
	}
    }  

     #no detection
     #del_data site is located at upstream of gene
     if ($temp_del_data_location<$hgGenome_refGene[$temp_index][5]) {
	 #location
	 if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	     $location="5prime";
	}
	 else {		    
	     $location="3prime";
	}
	 #$location="5prime";
	 return($location);
    }
     #del_data site is located at upstream of gene
     if ($temp_del_data_location>$hgGenome_refGene[$temp_index][6]) {
	 #location
	 if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	     $location="3prime";
	}
	 else {		    
	     $location="5prime";
	}
	 #$location="3prime";
	 return($location);
    }

     #no detection
     return("NA");	    	    
}


sub output_reduced_annotate_cross_data_after_blat {
    #output the data

    #output the data
    my $temp_blat_file=$ARGV[0];
     $temp_blat_file=~s/Chip/Blat/;

     open (OUTPUT, ">$temp_blat_file")
	 || die "Can't open $temp_blat_file $!";

     for (my $i=0; $i<=$reduced_cross_row_number; $i++) {
	 for (my $j=0; $j<=($reduced_cross_column_number+7); $j++) {
	     print OUTPUT "$reduced_cross_data[$i][$j]\t";
	}
	 print OUTPUT "$reduced_cross_data[$i][$reduced_cross_column_number+8]\n";
    }

     close(OUTPUT);
}


sub aggregate_blat_crossdata_by_anchor0_anchor1_BlatHit_Old {

    @aggregate_blat_cross_data=();
     $aggregate_blat_cross_row_number=0;
     $aggregate_blat_cross_column_number=$reduced_cross_column_number;

     for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
	 $aggregate_blat_cross_data[0][$j]=$reduced_cross_data[0][$j];
    }
     $aggregate_blat_cross_data[0][$reduced_cross_column_number+3]="Mate_Number";

     my $temp_start_pointer=1;
     my $temp_end_pointer=1;
     my $gene_ID_0=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
     my $gene_ID_1=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2];

     my $moving_pointer=0;
     while ($moving_pointer<$reduced_cross_row_number) {
	 while (($moving_pointer<$reduced_cross_row_number) && 
		($gene_ID_0 eq $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+1]) &&
		($gene_ID_1 eq $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+2])) {
	     $moving_pointer++;
	}

	 #start to assign genename the ordinal number
	 $temp_end_pointer=$moving_pointer;


	 #calculate the MADs
	 #for gene_0
	 @temp_row_data_forMAD=();
	 $temp_row_number_forMAD=0;
	 for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	     if ($reduced_cross_data[$h][$reduced_cross_column_number+5]=~/Yes/) {
		 $temp_row_number_forMAD++;
		 $temp_row_data_forMAD[$temp_row_number_forMAD-1]=$reduced_cross_data[$h][1];
	    }
	}	    
	 my $gene_0_MAD=blat_cross_data_MAD($temp_row_number_forMAD);

	 #for gene_1
	 @temp_row_data_forMAD=();
	 $temp_row_number_forMAD=0;
	 for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	     if ($reduced_cross_data[$h][$reduced_cross_column_number+5]=~/Yes/) {
		 $temp_row_number_forMAD++;
		 $temp_row_data_forMAD[$temp_row_number_forMAD-1]=$reduced_cross_data[$h][6];
	    }
	}	    
	 my $gene_1_MAD=blat_cross_data_MAD($temp_row_number_forMAD);


	 my $temp_BlatHit_yes=0;
	 my $temp_BlatHit_no=0;
	 for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	     if ($reduced_cross_data[$h][$reduced_cross_column_number+5]=~/Yes/) {
		 $temp_BlatHit_yes++;
	    }
	     else {
		 $temp_BlatHit_no++;
	    }
	}

	 if (($temp_BlatHit_yes>$temp_BlatHit_no) &&
	     ($gene_0_MAD>0) && ($gene_1_MAD>0)) {
	     $aggregate_blat_cross_row_number++;

	     for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		 $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$j]=
		     $reduced_cross_data[$temp_start_pointer][$j];
	    }
	     $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+3]=$temp_BlatHit_yes;
	}

	 #re-start the pointer;
	 $moving_pointer++;
	 $temp_start_pointer=$moving_pointer;
	 $temp_end_pointer=$moving_pointer;
	 $gene_ID_0=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
	 $gene_ID_1=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2];

    }

}

sub in_silico_fragment_by_anchor0_anchor1 {

    $reduced_cross_data[0][$reduced_cross_column_number+8]="InsilicoFrag";
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 $reduced_cross_data[$i][$reduced_cross_column_number+8]="Outof";	
    }
     my $MSK_stats_3stdev_distance=$MSK_stats_means+3*$MSK_stats_stdev; #for in-silico fragment length

     my $temp_start_pointer=1;
     my $temp_end_pointer=1;
     my $gene_ID_0=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
     my $gene_ID_1=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2];
    my (@strand_0)=split("\/\/", $reduced_cross_data[$temp_start_pointer][2]);
    my (@strand_1)=split("\/\/", $reduced_cross_data[$temp_start_pointer][7]);

     my $moving_pointer=0;
     while ($moving_pointer<$reduced_cross_row_number) {
	 my (@temp_strand_0)=split("\/\/", $reduced_cross_data[$moving_pointer+1][2]);
	 my (@temp_strand_1)=split("\/\/", $reduced_cross_data[$moving_pointer+1][7]);			  
	 while (($moving_pointer<$reduced_cross_row_number) && 
		($temp_strand_0[0] eq $strand_0[0]) && ($temp_strand_1[0] eq $strand_1[0]) && 
		($gene_ID_0 eq $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+1]) &&
		($gene_ID_1 eq $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+2])) {
	     $moving_pointer++;
	     if ($moving_pointer<$reduced_cross_row_number) {
		 @temp_strand_0=split("\/\/", $reduced_cross_data[$moving_pointer+1][2]);
		 @temp_strand_1=split("\/\/", $reduced_cross_data[$moving_pointer+1][7]);	
	     }		
	 }

	 #start to assign genename the ordinal number
	 $temp_end_pointer=$moving_pointer;


	 #put cluster of reads into a temp array
	 my %readID_to_index=();
	 my @forward_data=();
	 my $forward_row_number=0;
	 my @reverse_data=();
	 my $reverse_row_number=0;

	 for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	     if (($reduced_cross_data[$h][3]=~/\>/) &&
		 ($reduced_cross_data[$h][4]=~/A|C|G|T|a|c|g|t/) &&
		 ($reduced_cross_data[$h][8]=~/\>/) && 
		 ($reduced_cross_data[$h][9]=~/A|C|G|T|a|c|g|t/)) {
		 $readID_to_index{$reduced_cross_data[$h][3]}=$h;
		 $readID_to_index{$reduced_cross_data[$h][8]}=$h;

		 #forward read pairs
		 if ($reduced_cross_data[$h][2]=~/\+/) {
		     $forward_row_number++;
		     for (my $j=0; $j<=$reduced_cross_column_number; $j++) {
			 $forward_data[$forward_row_number][$j]=$reduced_cross_data[$h][$j];
		    }
		}
		 #reverse
		 else {
		     $reverse_row_number++;
		     for (my $j=0; $j<=$reduced_cross_column_number; $j++) {
			 $reverse_data[$reverse_row_number][$j]=$reduced_cross_data[$h][$j];
		    }
		}
	    } #real read data
	}

	 #generate the in-silico fragment length for each cluster
	 if ($forward_row_number>=1) {
	     @read_cluster_data=@forward_data;
	     $cluster_row_number=$forward_row_number;
	     breakpoint_and_distance_generator();
	     for (my $i=1; $i<=$forward_row_number; $i++) {
		 if (exists($readID_to_index{$forward_data[$i][3]})) {
		     my $index_anchor0=$readID_to_index{$forward_data[$i][3]};
		     #$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+5]=
			 #$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+5]."\-".
			 #$read_cluster_data[$i][$reduced_cross_column_number+5];
		     $reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=$read_cluster_data[$i][$reduced_cross_column_number+5];
		     if ($reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]<=$MSK_stats_3stdev_distance) {
			 $reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=
			     "InsilicoFrag\-".$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8];
		    }
		     else {
			 $reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=
			     "Outof\-".$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8];
		    }
		     #co-ordinate
		     $reduced_cross_data[$index_anchor0][2]=$read_cluster_data[$i][2];
		     $reduced_cross_data[$index_anchor0][7]=$read_cluster_data[$i][7];

		}
	    }

	     #re-generate the boundary
	     my @left_anchor0=();
	     my @right_anchor0=();
	     my @left_anchor1=();
	     my @right_anchor1=();
	     my $temp_anchor_number=0;
	     for (my $i=1; $i<=$forward_row_number; $i++) {
		 if (exists($readID_to_index{$forward_data[$i][3]})) {
		     my $index_anchor0=$readID_to_index{$forward_data[$i][3]};
		     my (@temp_info_anchor0)=split("\/\/", $read_cluster_data[$i][2]);
		     my (@temp_info_anchor1)=split("\/\/", $read_cluster_data[$i][7]);		    
		     if ($reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=~/InsilicoFrag/) {
			 $temp_anchor_number++;
			 $left_anchor0[$temp_anchor_number]=$temp_info_anchor0[1];
			 $right_anchor0[$temp_anchor_number]=$temp_info_anchor0[2];

			 $left_anchor1[$temp_anchor_number]=$temp_info_anchor1[1];
			 $right_anchor1[$temp_anchor_number]=$temp_info_anchor1[2];
		    }
		}
	    }

	     if ($temp_anchor_number>=2) {
		 $leftest_anchor0=$left_anchor0[1];
		 $rightest_anchor0=$right_anchor0[1];
		 $leftest_anchor1=$left_anchor1[1];
		 $rightest_anchor1=$right_anchor1[1];
		 for (my $i=2; $i<=$temp_anchor_number; $i++) {
		     if ($leftest_anchor0>$left_anchor0[$i]) { $leftest_anchor0=$left_anchor0[$i];}
		     if ($rightest_anchor0<$right_anchor0[$i]) { $rightest_anchor0=$right_anchor0[$i];}
		     if ($leftest_anchor1>$left_anchor1[$i]) { $leftest_anchor1=$left_anchor1[$i];}
		     if ($rightest_anchor1<$right_anchor1[$i]) { $rightest_anchor1=$right_anchor1[$i];}
		}

		 #re-assign
		 for (my $i=1; $i<=$forward_row_number; $i++) {
		     if (exists($readID_to_index{$forward_data[$i][3]})) {
			 my $index_anchor0=$readID_to_index{$forward_data[$i][3]};
			 if ($reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=~/InsilicoFrag/) {
			     #co-ordinate
			     my (@temp_info_anchor0)=split("\/\/",$reduced_cross_data[$index_anchor0][2]); 
			     $reduced_cross_data[$index_anchor0][2]=$temp_info_anchor0[0]."\/\/".
				 $leftest_anchor0."\/\/".$rightest_anchor0."\/\/".length($read_cluster_data[$i][4]);

			      my (@temp_info_anchor1)=split("\/\/",$reduced_cross_data[$index_anchor0][7]); 
			     $reduced_cross_data[$index_anchor0][7]=$temp_info_anchor1[0]."\/\/".
				 $leftest_anchor1."\/\/".$rightest_anchor1."\/\/".length($read_cluster_data[$i][9]);

			}
		    }
		}
	    } #end of re-boundary

	}

	 if ($reverse_row_number>=1) {
	     @read_cluster_data=@reverse_data;
	     $cluster_row_number=$reverse_row_number;
	     breakpoint_and_distance_generator();
	     for (my $i=1; $i<=$reverse_row_number; $i++) {
		 if (exists($readID_to_index{$reverse_data[$i][3]})) {
		     my $index_anchor0=$readID_to_index{$reverse_data[$i][3]};
		     #$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+5]=
			 #$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+5]."\-".
			 #$read_cluster_data[$i][$reduced_cross_column_number+5];
		     $reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=$read_cluster_data[$i][$reduced_cross_column_number+5];
		     if ($reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]<=$MSK_stats_3stdev_distance) {
			 $reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=
			     "InsilicoFrag\-".$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8];
		    }
		     else {
			 $reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=
			     "Outof\-".$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8];
		    }

		     #co-ordinate
		     $reduced_cross_data[$index_anchor0][2]=$read_cluster_data[$i][2];
		     $reduced_cross_data[$index_anchor0][7]=$read_cluster_data[$i][7];		    
		}
	    }
	     #re-generate the boundary
	     my @left_anchor0=();
	     my @right_anchor0=();
	     my @left_anchor1=();
	     my @right_anchor1=();
	     my $temp_anchor_number=0;
	     for (my $i=1; $i<=$reverse_row_number; $i++) {
		 if (exists($readID_to_index{$reverse_data[$i][3]})) {
		     my $index_anchor0=$readID_to_index{$reverse_data[$i][3]};
		     my (@temp_info_anchor0)=split("\/\/", $read_cluster_data[$i][2]);
		     my (@temp_info_anchor1)=split("\/\/", $read_cluster_data[$i][7]);		    
		     if ($reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=~/InsilicoFrag/) {
			 $temp_anchor_number++;
			 $left_anchor0[$temp_anchor_number]=$temp_info_anchor0[1];
			 $right_anchor0[$temp_anchor_number]=$temp_info_anchor0[2];

			 $left_anchor1[$temp_anchor_number]=$temp_info_anchor1[1];
			 $right_anchor1[$temp_anchor_number]=$temp_info_anchor1[2];
		    }
		}
	    }

	     if ($temp_anchor_number>=2) {
		 $leftest_anchor0=$left_anchor0[1];
		 $rightest_anchor0=$right_anchor0[1];
		 $leftest_anchor1=$left_anchor1[1];
		 $rightest_anchor1=$right_anchor1[1];
		 for (my $i=2; $i<=$temp_anchor_number; $i++) {
		     if ($leftest_anchor0>$left_anchor0[$i]) { $leftest_anchor0=$left_anchor0[$i];}
		     if ($rightest_anchor0<$right_anchor0[$i]) { $rightest_anchor0=$right_anchor0[$i];}
		     if ($leftest_anchor1>$left_anchor1[$i]) { $leftest_anchor1=$left_anchor1[$i];}
		     if ($rightest_anchor1<$right_anchor1[$i]) { $rightest_anchor1=$right_anchor1[$i];}
		}

		 #re-assign
		 for (my $i=1; $i<=$reverse_row_number; $i++) {
		     if (exists($readID_to_index{$reverse_data[$i][3]})) {
			 my $index_anchor0=$readID_to_index{$reverse_data[$i][3]};
			 if ($reduced_cross_data[$index_anchor0][$reduced_cross_column_number+8]=~/InsilicoFrag/) {
			     #co-ordinate
			     my (@temp_info_anchor0)=split("\/\/",$reduced_cross_data[$index_anchor0][2]);  
			     $reduced_cross_data[$index_anchor0][2]=$temp_info_anchor0[0]."\/\/".
				 $leftest_anchor0."\/\/".$rightest_anchor0."\/\/".length($read_cluster_data[$i][4]);
			     my (@temp_info_anchor1)=split("\/\/",$reduced_cross_data[$index_anchor0][7]);  
			     $reduced_cross_data[$index_anchor0][7]=$temp_info_anchor1[0]."\/\/".
				 $leftest_anchor1."\/\/".$rightest_anchor1."\/\/".length($read_cluster_data[$i][9]);

			}
		    }
		}
	    } #end of re-boundary

	}	

	 #re-start the pointer;
	 $moving_pointer++;
	 $temp_start_pointer=$moving_pointer;
	 $temp_end_pointer=$moving_pointer;
	if ($temp_start_pointer<=$reduced_cross_row_number) {
	    $gene_ID_0=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
	    $gene_ID_1=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2];
	    @strand_0=split("\/\/", $reduced_cross_data[$temp_start_pointer][2]);
	    @strand_1=split("\/\/", $reduced_cross_data[$temp_start_pointer][7]);
	}

    }

}

sub breakpoint_and_distance_generator {
    # my (@read_cluster_data, $cluster_row_number)=@_;
    
    #anchor0 & anchor1
    for (my $i=1; $i<=$cluster_row_number; $i++) {
	 #right end co-ordinate of read for anchor0
	 $read_cluster_data[$i][$reduced_cross_column_number+1]=$read_cluster_data[$i][1]+length($read_cluster_data[$i][4]);
	 #right end co-ordinate of read for anchor1
	 $read_cluster_data[$i][$reduced_cross_column_number+2]=$read_cluster_data[$i][6]+length($read_cluster_data[$i][9]);

	 #initialize
	 $read_cluster_data[$i][$reduced_cross_column_number+3]=0; #in-silico distance to breakpoint for anchor0
	 $read_cluster_data[$i][$reduced_cross_column_number+4]=0; #in-silico distance to breakpoint for anchor1
	 $read_cluster_data[$i][$reduced_cross_column_number+5]=0; #in-silico fragment length
    }

     #anchor0
     #forward
     @sort_cluster_data_index=();
     if ($read_cluster_data[1][2]=~/\+/) {
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $sort_cluster_data_index[$i][1]=$read_cluster_data[$i][$reduced_cross_column_number+1];
	}
	 quicksort_read_cluster_data(1, $cluster_row_number);
	 #rightest one for virtual breakpoint
	 my $breakpoint_for_anchor0=$sort_cluster_data_index[$cluster_row_number][1];

	 #generate the leftest boundary
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $sort_cluster_data_index[$i][1]=$read_cluster_data[$i][1];
	}
	 quicksort_read_cluster_data(1, $cluster_row_number);
	 #rightest one for virtual breakpoint
	 my $leftest_for_anchor0=$sort_cluster_data_index[1][1];	
	 my $rightest_for_anchor0=$breakpoint_for_anchor0;

	 #longest_read_length
	 my $longest_read_length=0;
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     if (length($read_cluster_data[$i][4])>$longest_read_length) {
		 $longest_read_length=length($read_cluster_data[$i][4]);
	    }
	}

	 #distance between virtual breakpoint and left end of read for anchor0
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $read_cluster_data[$i][$reduced_cross_column_number+3]=$breakpoint_for_anchor0-$read_cluster_data[$i][1];
	     if ($read_cluster_data[$i][$reduced_cross_column_number+3]<0) {
		 $read_cluster_data[$i][$reduced_cross_column_number+3]=0;
	    }

	     #add leftest and rightest into the co-ordinate for anchor0
	     my $temp_coordinate=$read_cluster_data[$i][2]."\/\/".$leftest_for_anchor0."\/\/".$rightest_for_anchor0.
		 "\/\/".$longest_read_length;
	     $read_cluster_data[$i][2]=$temp_coordinate;
	}
    }
     #reverse
     else {
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $sort_cluster_data_index[$i][1]=$read_cluster_data[$i][1];
	}
	 quicksort_read_cluster_data(1, $cluster_row_number);
	 #leftest one for virtual breakpoint
	 my $breakpoint_for_anchor0=$sort_cluster_data_index[1][1];

	 #generate the rightest boundary
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $sort_cluster_data_index[$i][1]=$read_cluster_data[$i][$reduced_cross_column_number+1];
	}
	 quicksort_read_cluster_data(1, $cluster_row_number);
	 #rightest one for virtual breakpoint
	 my $rightest_for_anchor0=$sort_cluster_data_index[$cluster_row_number][1];
	 my $leftest_for_anchor0=$breakpoint_for_anchor0;

	 #longest_read_length
	 my $longest_read_length=0;
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     if (length($read_cluster_data[$i][4])>$longest_read_length) {
		 $longest_read_length=length($read_cluster_data[$i][4]);
	    }
	}

	 #distance between virtual breakpoint and left end of read for anchor0
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $read_cluster_data[$i][$reduced_cross_column_number+3]=$read_cluster_data[$i][$reduced_cross_column_number+1]-$breakpoint_for_anchor0;
	     if ($read_cluster_data[$i][$reduced_cross_column_number+3]<0) {
		 $read_cluster_data[$i][$reduced_cross_column_number+3]=0;
	    }

	     #add leftest and rightest into the co-ordinate for anchor0
	     my $temp_coordinate=$read_cluster_data[$i][2]."\/\/".$leftest_for_anchor0."\/\/".$rightest_for_anchor0.
		 "\/\/".$longest_read_length;
	     $read_cluster_data[$i][2]=$temp_coordinate;
	}
    }

     #anchor1
     #forward
     if ($read_cluster_data[1][7]=~/\+/) {
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $sort_cluster_data_index[$i][1]=$read_cluster_data[$i][$reduced_cross_column_number+2];
	}	
	 quicksort_read_cluster_data(1, $cluster_row_number);
	 #rightest one for virtual breakpoint
	 my $breakpoint_for_anchor1=$sort_cluster_data_index[$cluster_row_number][1];

	 #generate the leftest boundary
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $sort_cluster_data_index[$i][1]=$read_cluster_data[$i][6];
	}
	 quicksort_read_cluster_data(1, $cluster_row_number);
	 my $leftest_for_anchor1=$sort_cluster_data_index[1][1];
	 my $rightest_for_anchor1=$breakpoint_for_anchor1;

	 #longest_read_length
	 my $longest_read_length=0;
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     if (length($read_cluster_data[$i][9])>$longest_read_length) {
		 $longest_read_length=length($read_cluster_data[$i][9]);
	    }
	}

	 #distance between virtual breakpoint and left end of read for anchor1
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $read_cluster_data[$i][$reduced_cross_column_number+4]=$breakpoint_for_anchor1-$read_cluster_data[$i][6];
	     if ($read_cluster_data[$i][$reduced_cross_column_number+4]<0) {
		 $read_cluster_data[$i][$reduced_cross_column_number+4]=0;
	    }

	     #add leftest and rightest into the co-ordinate for anchor1
	     my $temp_coordinate=$read_cluster_data[$i][7]."\/\/".$leftest_for_anchor1."\/\/".$rightest_for_anchor1.
		 "\/\/".$longest_read_length;
	     $read_cluster_data[$i][7]=$temp_coordinate;
	}
    }
     #reverse
     else {
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $sort_cluster_data_index[$i][1]=$read_cluster_data[$i][6];
	}
	 quicksort_read_cluster_data(1, $cluster_row_number);
	 #leftest one for virtual breakpoint
	 my $breakpoint_for_anchor1=$sort_cluster_data_index[1][1];

	 #generate the rightest boundary
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $sort_cluster_data_index[$i][1]=$read_cluster_data[$i][$reduced_cross_column_number+2];
	}	
	 quicksort_read_cluster_data(1, $cluster_row_number);
	 #rightest one for virtual breakpoint
	 my $rightest_for_anchor1=$sort_cluster_data_index[$cluster_row_number][1];
	 my $leftest_for_anchor1=$breakpoint_for_anchor1;

	 #longest_read_length
	 my $longest_read_length=0;
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     if (length($read_cluster_data[$i][9])>$longest_read_length) {
		 $longest_read_length=length($read_cluster_data[$i][9]);
	    }
	}

	 #distance between virtual breakpoint and left end of read for anchor1
	 for (my $i=1; $i<=$cluster_row_number; $i++) {
	     $read_cluster_data[$i][$reduced_cross_column_number+4]=$read_cluster_data[$i][$reduced_cross_column_number+2]-$breakpoint_for_anchor1;
	     if ($read_cluster_data[$i][$reduced_cross_column_number+4]<0) {
		 $read_cluster_data[$i][$reduced_cross_column_number+4]=0;
	    }

	     #add leftest and rightest into the co-ordinate for anchor1
	     my $temp_coordinate=$read_cluster_data[$i][7]."\/\/".$leftest_for_anchor1."\/\/".$rightest_for_anchor1.
		 "\/\/".$longest_read_length;
	     $read_cluster_data[$i][7]=$temp_coordinate;
	}
    }

     #calculate the in-silico fragment length
     for (my $i=1; $i<=$cluster_row_number; $i++) {
	 $read_cluster_data[$i][$reduced_cross_column_number+5]=$read_cluster_data[$i][$reduced_cross_column_number+3]+
	     $read_cluster_data[$i][$reduced_cross_column_number+4];
    }

     #debugging
     #my $index_anchor0=$readID_to_index{$read_cluster_data[1][3]};
     #my $gene_pair=$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+1]."\-".
	 #$reduced_cross_data[$index_anchor0][$reduced_cross_column_number+2];

     #if ($gene_pair=~/BCR\-JAK2/) {
	 #for (my $i=1; $i<=$cluster_row_number; $i++) {
	  #   for (my $j=0; $j<=($reduced_cross_column_number+5); $j++) {
	 #	print "$read_cluster_data[$i][$j]\t";
	 #   }
	 #    print "\n";
	 #}
     #}

     #return(@read_cluster_data);

}

sub quicksort_read_cluster_data {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1);

     if ($right>=$left) {
	 $pivot=$sort_cluster_data_index[$right][1];
	 # print "pivot= $pivot\n";
	 $pointer_left=$left-1;
	 $pointer_right=$right;
	 do
	 {
	     do
	     {
		 $pointer_left=$pointer_left+1;
	    } until ($sort_cluster_data_index[$pointer_left][1]>=$pivot);

	     do
	     { 
		 $pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_cluster_data_index[$pointer_right][1]<=$pivot));

	     $temp_1=$sort_cluster_data_index[$pointer_left][1];
	     $sort_cluster_data_index[$pointer_left][1]=$sort_cluster_data_index[$pointer_right][1];
	     $sort_cluster_data_index[$pointer_right][1]=$temp_1;

	     #$temp_2=$sort_cluster_data_index[$pointer_left][2];
	     #$sort_cluster_data_index[$pointer_left][2]=$sort_cluster_data_index[$pointer_right][2];
	     #$sort_cluster_data_index[$pointer_right][2]=$temp_2;

	} until ($pointer_right<=$pointer_left);

	 $sort_cluster_data_index[$pointer_right][1]=$sort_cluster_data_index[$pointer_left][1];
	 $sort_cluster_data_index[$pointer_left][1]=$sort_cluster_data_index[$right][1];
	 $sort_cluster_data_index[$right][1]=$temp_1;

	 #$sort_cluster_data_index[$pointer_right][2]=$sort_cluster_data_index[$pointer_left][2];
	 #$sort_cluster_data_index[$pointer_left][2]=$sort_cluster_data_index[$right][2];
	 #$sort_cluster_data_index[$right][2]=$temp_2;


	 quicksort_read_cluster_data($left, $pointer_left-1);
	 quicksort_read_cluster_data($pointer_left+1, $right);
    }
}

sub get_Spanner_stats_data {
    #my $in_data_file="MSK.stats";
    
    #open (INPUT, "<$in_data_file")
    #|| die "Can't open $in_data_file $!";

    $MSK_stats_means=$ARGV[3];
    $MSK_stats_stdev=$ARGV[4];    

     #$line_indicator=0;
     #input data from files

  #   while (<INPUT>) {
 #	$text=$_;
 #	$text=~s/|\n//g;
 #	(@line)=split(" ", $text);

 #	$line_indicator++;

	 #readin the data
 #	if ($line_indicator==3) {
 #	    $MSK_stats_means=$line[1];
 #	    $MSK_stats_stdev=$line[2];
 #	    last;
 #	}	
 #   } #while(INPUT)
 #    
 #    close(INPUT);
     ##print "$MSK_stats_means $MSK_stats_stdev\n";
}


sub filter_reduced_cross_data_by_insilicoFrag {

    my $temp_data_row_number=0;
     
     #initialize blat hit
     $reduced_cross_data[0][$reduced_cross_column_number+3]="Read0_BlatHit";
     $reduced_cross_data[0][$reduced_cross_column_number+4]="Read1_BlatHit";
     $reduced_cross_data[0][$reduced_cross_column_number+5]="BlatCall";

     $reduced_cross_data[0][$reduced_cross_column_number+6]="Read0_Loc";
     $reduced_cross_data[0][$reduced_cross_column_number+7]="Read1_Loc";

     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 $reduced_cross_data[$i][$reduced_cross_column_number+3]=0;
	 $reduced_cross_data[$i][$reduced_cross_column_number+4]=0;
	 $reduced_cross_data[$i][$reduced_cross_column_number+5]="No";	
	 $reduced_cross_data[$i][$reduced_cross_column_number+6]="NA";	
	 $reduced_cross_data[$i][$reduced_cross_column_number+7]="NA";	
    }

     #header
     #filter the data by in-silico fragment
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 if ($reduced_cross_data[$i][$reduced_cross_column_number+8]=~/InsilicoFrag/) {
	     $temp_data_row_number++;
	     for (my $j=0; $j<=($reduced_cross_column_number+8); $j++) {
		 $reduced_cross_data[$temp_data_row_number][$j]=$reduced_cross_data[$i][$j];
	    }
	}
    }

     #re-assign the array row number
     $reduced_cross_row_number=$temp_data_row_number;

}

sub aggregate_blat_crossdata_by_anchor0_anchor1_BlatHit {

    @aggregate_blat_cross_data=();
    $aggregate_blat_cross_row_number=0;
    $aggregate_blat_cross_column_number=$reduced_cross_column_number;
    my $PEread_length=$ARGV[5]; #read_length=50 or 75  
    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
	$aggregate_blat_cross_data[0][$j]=$reduced_cross_data[0][$j];
    }
    $aggregate_blat_cross_data[0][$reduced_cross_column_number+3]="Mate_Number";
    $aggregate_blat_cross_data[0][$reduced_cross_column_number+4]="Read0_count";
    $aggregate_blat_cross_data[0][$reduced_cross_column_number+5]="Read1_count";
    
    $blast_sampleID="T";
    $aggregate_blat_cross_data[0][$reduced_cross_column_number+6]="Junction_Number";
    $aggregate_blat_cross_data[0][$reduced_cross_column_number+7]="Blast_count";
    #initialize
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	$aggregate_blat_cross_data[$i][$reduced_cross_column_number+6]=0;
	$aggregate_blat_cross_data[$i][$reduced_cross_column_number+7]=0;
    }

    $aggregate_blat_cross_data[0][$reduced_cross_column_number+9]="Read0_Loc";
    $aggregate_blat_cross_data[0][$reduced_cross_column_number+10]="Read1_Loc";

    my $temp_start_pointer=1;
    my $temp_end_pointer=1;
    my $gene_ID_0=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
    my $gene_ID_1=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2];
    my (@strand_0)=split("\/\/", $reduced_cross_data[$temp_start_pointer][2]);
    my (@strand_1)=split("\/\/", $reduced_cross_data[$temp_start_pointer][7]);
 
     my $moving_pointer=0;
     while ($moving_pointer<$reduced_cross_row_number) {
	my (@temp_strand_0)=split("\/\/", $reduced_cross_data[$moving_pointer+1][2]); 
	my (@temp_strand_1)=split("\/\/", $reduced_cross_data[$moving_pointer+1][7]); 
	while (($moving_pointer<$reduced_cross_row_number) && 
	       ($temp_strand_0[0] eq $strand_0[0]) && ($temp_strand_1[0] eq $strand_1[0]) && 
	       ($gene_ID_0 eq $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+1]) &&
	       ($gene_ID_1 eq $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+2])) {
	    $moving_pointer++;
	    if ($moving_pointer<$reduced_cross_row_number) {
		@temp_strand_0=split("\/\/", $reduced_cross_data[$moving_pointer+1][2]); 
		@temp_strand_1=split("\/\/", $reduced_cross_data[$moving_pointer+1][7]); 
	    }
	}
	
	 #start to assign genename the ordinal number
	 $temp_end_pointer=$moving_pointer;

	 #calculate the MADs
	 #for gene_0
 #	@temp_row_data_forMAD=();
 #	$temp_row_number_forMAD=0;
 #	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
 #	    if ($reduced_cross_data[$h][$reduced_cross_column_number+5]=~/Yes/) {
 #		$temp_row_number_forMAD++;
 #		$temp_row_data_forMAD[$temp_row_number_forMAD-1]=$reduced_cross_data[$h][1];
 #	   }
 #	}	    
 #	my $gene_0_MAD=blat_cross_data_MAD($temp_row_number_forMAD);

	 #for gene_1
 #	@temp_row_data_forMAD=();
 #	$temp_row_number_forMAD=0;
 #	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
 #	    if ($reduced_cross_data[$h][$reduced_cross_column_number+5]=~/Yes/) {
 #		$temp_row_number_forMAD++;
 #		$temp_row_data_forMAD[$temp_row_number_forMAD-1]=$reduced_cross_data[$h][6];
 #	   }
 #	}	    
 #	my $gene_1_MAD=blat_cross_data_MAD($temp_row_number_forMAD);


	#put cluster of reads into a temp array
	#my %readID_to_index=();
	my @forward_data=();
	my $forward_row_number=0;
	my $forward_junctionSpan=0;
	my $forward_blast=0;
	my @reverse_data=();
	my $reverse_row_number=0;
	my $reverse_junctionSpan=0;
	my $reverse_blast=0;
	
	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	    if (($reduced_cross_data[$h][3]=~/\>/) &&
		($reduced_cross_data[$h][4]=~/A|C|G|T|a|c|g|t/) &&
		($reduced_cross_data[$h][8]=~/\>/) && 
		($reduced_cross_data[$h][9]=~/A|C|G|T|a|c|g|t/)) {
		#$readID_to_index{$reduced_cross_data[$h][3]}=$h;
		#$readID_to_index{$reduced_cross_data[$h][8]}=$h;
		
		#forward read pairs
		if ($reduced_cross_data[$h][2]=~/\+/) {
		    $forward_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+7); $j++) {
			$forward_data[$forward_row_number][$j]=$reduced_cross_data[$h][$j];
		    }
		}
		#reverse
		else {
		    $reverse_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+7); $j++) {
			$reverse_data[$reverse_row_number][$j]=$reduced_cross_data[$h][$j];
		    }
		}
	    } #real read data
	}
	
	#forward
	#generate the Blatcall for both yes and no
	my $temp_BlatHit_forward_yes=0;
	my $temp_BlatHit_forward_no=0;
	my $temp_BlatHit_forward_dup=0;
	for (my $h=1; $h<=$forward_row_number; $h++) {
	    #count the junction spanning reads
	    my (@tempSeq1)=split("", $forward_data[$h][4]);
	    my (@tempSeq2)=split("", $forward_data[$h][9]);
	    if ((($#tempSeq1+1)<=($PEread_length*0.9)) || (($#tempSeq2+1)<=($PEread_length*0.9))) { 
		$forward_junctionSpan++;
	    }
	    if ($forward_data[$h][$reduced_cross_column_number+5]=~/Yes/) {
		$temp_BlatHit_forward_yes++;
		#distinguish normal from blast 
		if ($forward_data[$h][3]=~/$blast_sampleID/) {
		    $forward_blast++;
		}
		#else {
		#    $forward_normal++;
		#}
	    }
	    elsif ($forward_data[$h][$reduced_cross_column_number+5]=~/Duplicate/) {
		$temp_BlatHit_forward_dup++;
	    }
	    else {
		$temp_BlatHit_forward_no++;
	    }
	    
	    
	}
	
	#extract the reads in both most left and right for both mates
	my $forward_seq="";
	if ($forward_row_number>=1) {
	    #for mate1
	    for (my $i=1; $i<=$forward_row_number; $i++) {
		$sort_fusioncluster_data_index[$i][1]=$forward_data[$i][1];
		$sort_fusioncluster_data_index[$i][2]=$i;
	    }
	    quicksort_read_fusioncluster_data(1, $forward_row_number);
	    my $next_to_left_index=2;
	    if ($next_to_left_index>$forward_row_number) { $next_to_left_index=$forward_row_number;}
	    my $left_seq_mate_1=$forward_data[$sort_fusioncluster_data_index[$next_to_left_index][2]][4];
	    my $right_seq_mate_1=$forward_data[$sort_fusioncluster_data_index[$forward_row_number][2]][4];
	    
	    #for mate2
	    for (my $i=1; $i<=$forward_row_number; $i++) {
		$sort_fusioncluster_data_index[$i][1]=$forward_data[$i][6];
		$sort_fusioncluster_data_index[$i][2]=$i;
	    }
	    quicksort_read_fusioncluster_data(1, $forward_row_number);
	    ##if ($forward_data[$forward_row_number][7]=~/\+/)
	    $next_to_left_index=2;
	    if ($next_to_left_index>$forward_row_number) { $next_to_left_index=$forward_row_number;}
	    my $left_seq_mate_2=$forward_data[$sort_fusioncluster_data_index[$next_to_left_index][2]][9];
	    my $right_seq_mate_2=$forward_data[$sort_fusioncluster_data_index[$forward_row_number][2]][9];
	    
	    if ($forward_data[$forward_row_number][7]=~/\-/) {
		$left_seq_mate_2=$forward_data[$sort_fusioncluster_data_index[1][2]][9];
		my $next_to_right_index=$forward_row_number-1;
		if ($next_to_right_index<1) { $next_to_right_index=1;}
		$right_seq_mate_2=$forward_data[$sort_fusioncluster_data_index[$next_to_right_index][2]][9];
	    }
	    $forward_seq=$left_seq_mate_1."\/\/".$right_seq_mate_1."\/\/".
		$left_seq_mate_2."\/\/".$right_seq_mate_2;
	}
	
	
	if (($forward_row_number>=1) && ($temp_BlatHit_forward_yes>($temp_BlatHit_forward_no/2)) && 
	    ($temp_BlatHit_forward_dup<=1)) {
	    #($gene_0_MAD>0) && ($gene_1_MAD>0)) {
	    $aggregate_blat_cross_row_number++;
	    
	    
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$j]=
		    $forward_data[1][$j];
	    }
	    #read location
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+9]=
		$forward_data[1][$reduced_cross_column_number+6];
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+10]=
		$forward_data[1][$reduced_cross_column_number+7];
	    
	    #sequence
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+11]=$forward_seq;	
	    
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+3]=$temp_BlatHit_forward_yes;
	    
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+4]=
		$gene_reads_count{$forward_data[1][$reduced_cross_column_number+1]};
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+5]=
		$gene_reads_count{$forward_data[1][$reduced_cross_column_number+2]};	
	    
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+6]=$forward_junctionSpan;
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+7]=$forward_blast;
	    
	}
	
	#reverse
	#generate the Blatcall for both yes and no
	my $temp_BlatHit_reverse_yes=0;
	my $temp_BlatHit_reverse_no=0;
	my $temp_BlatHit_reverse_dup=0;
	for (my $h=1; $h<=$reverse_row_number; $h++) {
	    #count the junction spanning reads
	    my (@tempSeq1)=split("", $reverse_data[$h][4]);
	    my (@tempSeq2)=split("", $reverse_data[$h][9]);
	    if ((($#tempSeq1+1)<=($PEread_length*0.9)) || (($#tempSeq2+1)<=($PEread_length*0.9))) { 
		$reverse_junctionSpan++;
	    }

	    if ($reverse_data[$h][$reduced_cross_column_number+5]=~/Yes/) {
		$temp_BlatHit_reverse_yes++;
		
		#distinguish normal from blast 
		if ($reverse_data[$h][3]=~/$blast_sampleID/) {
		    $reverse_blast++;
		}
		#else {
		#    $reverse_normal++;
		#}
	    }
	    elsif ($reverse_data[$h][$reduced_cross_column_number+5]=~/Duplicate/) {
		$temp_BlatHit_reverse_dup++;
	    }
	    else {
		$temp_BlatHit_reverse_no++;
	    }
	}

	#extract the reads in both most left and right for both mates
	my $reverse_seq="";
	if ($reverse_row_number>=1) {
	    
	    #for mate1
	    for (my $i=1; $i<=$reverse_row_number; $i++) {
		$sort_fusioncluster_data_index[$i][1]=$reverse_data[$i][1];
		$sort_fusioncluster_data_index[$i][2]=$i;
	    }
	    quicksort_read_fusioncluster_data(1, $reverse_row_number);
	    my $left_seq_mate_1=$reverse_data[$sort_fusioncluster_data_index[1][2]][4];
	    my $next_to_right_index=$reverse_row_number-1;
	    if ($next_to_right_index<1) { $next_to_right_index=1;}
	    my $right_seq_mate_1=$reverse_data[$sort_fusioncluster_data_index[$next_to_right_index][2]][4];
	    
	    #for mate2
	    for (my $i=1; $i<=$reverse_row_number; $i++) {
		$sort_fusioncluster_data_index[$i][1]=$reverse_data[$i][6];
		$sort_fusioncluster_data_index[$i][2]=$i;
	    }
	    quicksort_read_fusioncluster_data(1, $reverse_row_number);
	    ##if ($reverse_data[$reverse_row_number][7]=~/\+/)
	    my $next_to_left_index=2;
	    if ($next_to_left_index>$reverse_row_number) { $next_to_left_index=$reverse_row_number;}
	    my $left_seq_mate_2=$reverse_data[$sort_fusioncluster_data_index[$next_to_left_index][2]][9];
	    my $right_seq_mate_2=$reverse_data[$sort_fusioncluster_data_index[$reverse_row_number][2]][9];
	    
	    if ($reverse_data[$reverse_row_number][7]=~/\-/) {
		$left_seq_mate_2=$reverse_data[$sort_fusioncluster_data_index[1][2]][9];
		my $next_to_right_index=$reverse_row_number-1;
		if ($next_to_right_index<1) { $next_to_right_index=1;}
		$right_seq_mate_2=$reverse_data[$sort_fusioncluster_data_index[$next_to_right_index][2]][9];
	    }
	    
	    $reverse_seq=$left_seq_mate_1."\/\/".$right_seq_mate_1."\/\/".
		$left_seq_mate_2."\/\/".$right_seq_mate_2;
	}
	
	
	if (($reverse_row_number>=1) &&  ($temp_BlatHit_reverse_yes>($temp_BlatHit_reverse_no/2)) && 
	    ($temp_BlatHit_reverse_dup<=1)) {
	    #($gene_0_MAD>0) && ($gene_1_MAD>0)) {
	    $aggregate_blat_cross_row_number++;
	    
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$j]=
		    $reverse_data[1][$j];
	    }
	    #read location
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+9]=
		$reverse_data[1][$reduced_cross_column_number+6];
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+10]=
		$reverse_data[1][$reduced_cross_column_number+7];	
	    
	    #sequence
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+11]=$reverse_seq;	
	    
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+3]=$temp_BlatHit_reverse_yes;

	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+4]=
		$gene_reads_count{$reverse_data[1][$reduced_cross_column_number+1]};
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+5]=
		$gene_reads_count{$reverse_data[1][$reduced_cross_column_number+2]};	
	    
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+6]=$reverse_junctionSpan;
	    $aggregate_blat_cross_data[$aggregate_blat_cross_row_number][$reduced_cross_column_number+7]=$reverse_blast;
	    
	}
	
	#re-start the pointer;
	$moving_pointer++;
	$temp_start_pointer=$moving_pointer;
	$temp_end_pointer=$moving_pointer;
	if ($temp_start_pointer<=$reduced_cross_row_number) {
	    $gene_ID_0=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
	    $gene_ID_1=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2];
	    @strand_0=split("\/\/", $reduced_cross_data[$temp_start_pointer][2]);
	    @strand_1=split("\/\/", $reduced_cross_data[$temp_start_pointer][7]);
	}
    }
    
}

sub quicksort_read_fusioncluster_data {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
     my ($pivot, $temp_1);

     if ($right>=$left) {
	 $pivot=$sort_fusioncluster_data_index[$right][1];
	 # print "pivot= $pivot\n";
	 $pointer_left=$left-1;
	 $pointer_right=$right;
	 do
	 {
	     do
	     {
		 $pointer_left=$pointer_left+1;
	    } until ($sort_fusioncluster_data_index[$pointer_left][1]>=$pivot);

	     do
	     { 
		 $pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_fusioncluster_data_index[$pointer_right][1]<=$pivot));

	     $temp_1=$sort_fusioncluster_data_index[$pointer_left][1];
	     $sort_fusioncluster_data_index[$pointer_left][1]=$sort_fusioncluster_data_index[$pointer_right][1];
	     $sort_fusioncluster_data_index[$pointer_right][1]=$temp_1;

	     $temp_2=$sort_fusioncluster_data_index[$pointer_left][2];
	     $sort_fusioncluster_data_index[$pointer_left][2]=$sort_fusioncluster_data_index[$pointer_right][2];
	     $sort_fusioncluster_data_index[$pointer_right][2]=$temp_2;

	} until ($pointer_right<=$pointer_left);

	 $sort_fusioncluster_data_index[$pointer_right][1]=$sort_fusioncluster_data_index[$pointer_left][1];
	 $sort_fusioncluster_data_index[$pointer_left][1]=$sort_fusioncluster_data_index[$right][1];
	 $sort_fusioncluster_data_index[$right][1]=$temp_1;

	 $sort_fusioncluster_data_index[$pointer_right][2]=$sort_fusioncluster_data_index[$pointer_left][2];
	 $sort_fusioncluster_data_index[$pointer_left][2]=$sort_fusioncluster_data_index[$right][2];
	 $sort_fusioncluster_data_index[$right][2]=$temp_2;


	 quicksort_read_fusioncluster_data($left, $pointer_left-1);
	 quicksort_read_fusioncluster_data($pointer_left+1, $right);
    }
}

sub generate_predicted_fusion_reference_fasta {
    
    #output the data
    my $temp_FusionRef_file=$ARGV[0];
    $temp_FusionRef_file=~s/Chip/FusionRef/;
     $temp_FusionRef_file=~s/\.txt/\.fa/;
     open (OUTPUT, ">$temp_FusionRef_file")
	 || die "Can't open $temp_FusionRef_file $!";

     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 if (!($aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=~/NA/)) {
	     my $temp_read_ID="\>chr\_".$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12];
	     print OUTPUT "$temp_read_ID\n";

	     #first half
	     my @predicted_sequence=split("\/\/", $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]);
	     my @temp_read1=split("\/\/", $aggregate_blat_cross_data[$i][2]);
	     #print "temp_readlength_1=$temp_read1[1] $temp_read1[2] $temp_read1[3]\n";
	     my $sequence_length1=$temp_read1[3]-int($temp_read1[3]*0.1);
	     if ($sequence_length1>length($predicted_sequence[0])) {
		 $sequence_length1=length($predicted_sequence[0]);
	    }

	     my @temp_predicted_first_half=split("", $predicted_sequence[0]);
	     my $start_position=$#temp_predicted_first_half-$sequence_length1+1;
	     my $predicted_first_half=join("", @temp_predicted_first_half[$start_position..$#temp_predicted_first_half]);

	     print OUTPUT "$predicted_first_half";
	     print OUTPUT "N";

	     #second half
	     my @temp_read2=split("\/\/", $aggregate_blat_cross_data[$i][7]);
	     my $sequence_length2=$temp_read2[3]-int($temp_read2[3]*0.1);
	     $sequence_length2--;
	     if ($sequence_length2>length($predicted_sequence[1])) {
		 $sequence_length2=length($predicted_sequence[1]);
	    }

	     my @temp_predicted_second_half=split("", $predicted_sequence[1]);
	     my $predicted_second_half=join("", @temp_predicted_second_half[0..$sequence_length2]);

	     print OUTPUT "$predicted_second_half\n";


	}
    }

    close(OUTPUT);
}


sub blat_cross_data_MAD {

     my ($temp_cross_row_number)=@_;

     if ($temp_cross_row_number>=2) {
	 my @sort_row_data=sort {$a<=>$b} @temp_row_data_forMAD;

	 my $m1=$sort_row_data[int($temp_cross_row_number/2)];
	 my $m2=$sort_row_data[int($temp_cross_row_number/2)-1];
	 my $median=($m1+$m2)/2;

	 my @MAD=();
	 for (my $i=0; $i<$temp_cross_row_number; $i++) {	  
	     $MAD[$i]=abs($temp_row_data_forMAD[$i]-$median);
	}

	 my @sort_MAD=sort {$a<=>$b} @MAD;

	 my $a1=$sort_MAD[int($temp_cross_row_number/2)];
	 my $a2=$sort_MAD[int($temp_cross_row_number/2)-1];

	 my $absolute_deviation=($a1+$a2)/2;

	 return($absolute_deviation);
    }
     else {
	 return(0);
    }

} #end ofsubroutine

sub assess_orientation_both_reads_refGenes {

     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+8]="Orientation";

     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+8]="NA";
	 if ((exists($cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]})) &&
	     (exists($cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]}))) {
	     my $ori_index1=$cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]};
	     my $ori_index2=$cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]};
	     my $ori_gene1=$hgGenome_refGene[$ori_index1][4]; #GENE STRAND
	     my $ori_gene2=$hgGenome_refGene[$ori_index2][4]; #gene strand
	     if ((($aggregate_blat_cross_data[$i][2]=~/\+/) && ($aggregate_blat_cross_data[$i][7]=~/\-/)) ||
		 (($aggregate_blat_cross_data[$i][2]=~/\-/) && ($aggregate_blat_cross_data[$i][7]=~/\+/))) {
		 if ((($ori_gene1=~/\+/) && ($ori_gene2=~/\+/)) ||
		     (($ori_gene1=~/\-/) && ($ori_gene2=~/\-/))) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+8]="Good";
		}
		 else {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+8]="Bad";
		}
	    } #if ((($aggregate_blat_cross_data[$i][2]=~/F/)
	     elsif ((($aggregate_blat_cross_data[$i][2]=~/\+/) && ($aggregate_blat_cross_data[$i][7]=~/\+/)) ||
		    (($aggregate_blat_cross_data[$i][2]=~/\-/) && ($aggregate_blat_cross_data[$i][7]=~/\-/))) {
		 if ((($ori_gene1=~/\+/) && ($ori_gene2=~/\-/)) ||
		     (($ori_gene1=~/\-/) && ($ori_gene2=~/\+/))) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+8]="Good";
		}
		 else {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+8]="Bad";
		}
	    } #if ((($aggregate_blat_cross_data[$i][2]=~/F/)

	} #if ((exists($cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]}))
    } #for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++)

}

sub output_aggregate_reduced_annotate_cross_data_after_blat {
     #output the data

     #output the data
     my $temp_blataggregate_file=$ARGV[6]; #Virus integration output filename (named by user)
     #$temp_blataggregate_file=~s/Chip/BlatAggregate/;

     #add the true or false remark
     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+13]="Status";
     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 if (($aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]=~/prime/) ||
	     ($aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]=~/prime/) ||
	     ($aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+3]<=2) || 
	     ($aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+8]=~/Bad/)) {
	     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+13]="False";
	}
	 else {
	     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+13]="True";
	}	    
    }

     open (OUTPUT, ">$temp_blataggregate_file")
	 || die "Can't open $temp_blataggregate_file $!";

     $aggregate_blat_cross_data[0][0]="Virus_ChrID";
     $aggregate_blat_cross_data[0][1]="Virus_Integration_Location";

     $aggregate_blat_cross_data[0][5]="Gene_ChrID";
     $aggregate_blat_cross_data[0][6]="Gene_Integration_Location";

     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+1]="Viral_Transcript";
     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+2]="Host_Gene";
     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+3]="Discordant_Read_Pairs";
     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+6]="JunctionSpanning_Reads";
     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+11]="PCR_primer_template";

     #header
     for (my $i=0; $i<=0; $i++) {
	 print OUTPUT "$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]\t";
	 print OUTPUT "$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]\t";

	 print OUTPUT "$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+3]\t";
	 print OUTPUT "$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+6]\t";

	 print OUTPUT "$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]\t";	 
	 
	 print OUTPUT "$aggregate_blat_cross_data[$i][0]\t";
	 print OUTPUT "$aggregate_blat_cross_data[$i][1]\t";
	 print OUTPUT "$aggregate_blat_cross_data[$i][5]\t";
	 print OUTPUT "$aggregate_blat_cross_data[$i][6]\n";
	 	 
     }
     #content
     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 if ($aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+3]>=4) {
	     my (@tempVirus)=split("\/\/", $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]);
	     print OUTPUT "$tempVirus[0]\t";
	     my (@preGene)=split("\/\/", $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]);
	     my (@preGeneExon)=split("\-", $preGene[1]);
	     my ($tempGene)=$preGene[0]."\/\/".$preGeneExon[0];
	     print OUTPUT "$tempGene\t";
	     
	     print OUTPUT "$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+3]\t";
	     print OUTPUT "$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+6]\t";	 
	     
	     print OUTPUT "$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]\t";
	     print OUTPUT "$aggregate_blat_cross_data[$i][0]\t";
	     print OUTPUT "$aggregate_blat_cross_data[$i][1]\t";
	     print OUTPUT "$aggregate_blat_cross_data[$i][5]\t";
	     print OUTPUT "$aggregate_blat_cross_data[$i][6]\n";
	 }
     }
     
     close(OUTPUT);
}

 #start to deal with UCSC refGene data
sub get_hgGenome_refGene {
     my $in_data_file=$ARGV[1];
     #my $in_data_file="hgGenome_cDNA_refGene_RIS.txt";
     #my $in_data_file="hgGenome_miRNA_RIS.txt";
     open (INPUT, "<$in_data_file")
	 || die "Can't open $in_data_file $!";


   @hgGenome_refGene=();

   $hgGenome_refGene_row_number=0; #genes number
   $line_indicator=0;
   $first_line=1;

   #input data from files

   while (<INPUT>) {
     $text=$_;
     $text=~s/|\n//g;

     (@line)=split("\t", $text);

     $hgGenome_refGene_column_number=0;
     for $read_data(@line) {
       # $read_data=~s/ //g;
       $hgGenome_refGene[$hgGenome_refGene_row_number][$hgGenome_refGene_column_number]=$read_data;
       $hgGenome_refGene_column_number++;
    }

     # detect if there is any missing data point in a specific row
     if ($first_line==1) {
       $hgGenome_refGene_column_number_previous=$hgGenome_refGene_column_number;
    }
     else {
       if ($hgGenome_refGene_column_number_previous!=$hgGenome_refGene_column_number) {
	 $line_indicator=1;
      }
       $hgGenome_refGene_column_number_previous=$hgGenome_refGene_column_number;
    }

     if ($line_indicator==1) {

       #print "there is missing data point in line= $hgGenome_refGene_row_number\n";
       $line_indicator=0;
    }
       # end of detection of missing data point

     if (!($hgGenome_refGene[$hgGenome_refGene_row_number][3]=~/\_/)) { 
	 $hgGenome_refGene_row_number++; 
    }

     #$hgGenome_refGene_row_number++;
     $first_line=0;

  } #while(INPUT)

   close(INPUT);

   $hgGenome_refGene_row_number--;
   $hgGenome_refGene_column_number--;

   #output all data

   #for (my $i=0; $i<=$hgGenome_refGene_row_number; $i++) {
    # for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
    #   print "$hgGenome_refGene[$i][$j]	";
    #}
    # print "\n";
   #}

}

sub regenerate_hgGenome_refGene_from_fiveprime_to_threeprime {

   for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
     $hgGenome_refGene[$i][3]=~s/chr//;
     if ($hgGenome_refGene[$i][3] eq "X") { $hgGenome_refGene[$i][3]=23;}
     elsif ($hgGenome_refGene[$i][3] eq "Y") { $hgGenome_refGene[$i][3]=24;}
     elsif ($hgGenome_refGene[$i][3] eq "M") { $hgGenome_refGene[$i][3]=25;}

     #convert reverse strand into normal strand

     #if ($hgGenome_refGene[$i][3] eq "-") {
     #  my $five_start=$chrom_length[$hgGenome_refGene[$i][2]][1]-$hgGenome_refGene[$i][5];
     #  my $five_end=$chrom_length[$hgGenome_refGene[$i][2]][1]-$hgGenome_refGene[$i][4];

     #  $hgGenome_refGene[$i][4]=$five_start;
     #  $hgGenome_refGene[$i][5]=$five_end;
     #}

  }

   #output the data
   #for (my $i=0; $i<=$hgGenome_refGene_row_number; $i++) {
   #  for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
   #    print "$hgGenome_refGene[$i][$j]	";
   # }
   #  print "\n";
   #}

}



sub get_chromosome_genomic_location_hgGenome_refGene {

   @hgGenome_refGene_probe_number_in_chrom=();
   @hgGenome_refGene_chrom_genome=();

   for (my $i=0; $i<=25; $i++) {
     $hgGenome_refGene_probe_number_in_chrom[$i]=0;
  }

   for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
     $hgGenome_refGene_probe_number_in_chrom[$hgGenome_refGene[$i][3]]++;
     $hgGenome_refGene_chrom_genome[$hgGenome_refGene[$i][3]][$hgGenome_refGene_probe_number_in_chrom[$hgGenome_refGene[$i][3]]][1]=$hgGenome_refGene[$i][5];
     $hgGenome_refGene_chrom_genome[$hgGenome_refGene[$i][3]][$hgGenome_refGene_probe_number_in_chrom[$hgGenome_refGene[$i][3]]][2]=$i;
  }

   #sort out each BAC in terms of genomic location chromosome by chromosome
   for (my $i=1; $i<=25; $i++) {
     quicksort_lookup_table_hgGenome_refGene($i, 1, $hgGenome_refGene_probe_number_in_chrom[$i]);
  }

}

sub quicksort_lookup_table_hgGenome_refGene {
   my ($local_chrom_number, $left, $right)=@_;
   my ($pointer_left, $pointer_right);
   my ($pivot, $temp_1, $temp_2);

   if ($right>=$left)
     {
       $pivot=$hgGenome_refGene_chrom_genome[$local_chrom_number][$right][1];
      # print "pivot= $pivot\n";
       $pointer_left=$left-1;
       $pointer_right=$right;
       do
	 {
	   do
	     {
	       $pointer_left=$pointer_left+1;
	    } until ($hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][1]>=$pivot);

	   do
	     { 
	       $pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_right][1]<=$pivot));

	   $temp_1=$hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][1];
	   $hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][1]=$hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_right][1];
	   $hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_right][1]=$temp_1;

	   $temp_2=$hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][2];
	   $hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][2]=$hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_right][2];
	   $hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_right][2]=$temp_2;

	} until ($pointer_right<=$pointer_left);

       $hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_right][1]=$hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][1];
       $hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][1]=$hgGenome_refGene_chrom_genome[$local_chrom_number][$right][1];
       $hgGenome_refGene_chrom_genome[$local_chrom_number][$right][1]=$temp_1;

       $hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_right][2]=$hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][2];
       $hgGenome_refGene_chrom_genome[$local_chrom_number][$pointer_left][2]=$hgGenome_refGene_chrom_genome[$local_chrom_number][$right][2];
       $hgGenome_refGene_chrom_genome[$local_chrom_number][$right][2]=$temp_2;


       quicksort_lookup_table_hgGenome_refGene($local_chrom_number, $left, $pointer_left-1);
       quicksort_lookup_table_hgGenome_refGene($local_chrom_number, $pointer_left+1, $right);
    }
}

sub generate_hgGenome_refGene_Pointer_by_chrom {

   my @temp_hgGenome_refGene=();

   #first header line
   for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
     $temp_hgGenome_refGene[0][$j]=$hgGenome_refGene[0][$j];
  #   print "$hgGenome_refGene[0][$j]	";
  }
  # print "\n";    

   #output data content
   my $row_index=0;
   for (my $i=1; $i<=25; $i++) {
     for (my $j=1; $j<=$hgGenome_refGene_probe_number_in_chrom[$i]; $j++) {
       $row_index++;
       for (my $h=0; $h<=$hgGenome_refGene_column_number; $h++) {
	 $temp_hgGenome_refGene[$row_index][$h]=$hgGenome_refGene[$hgGenome_refGene_chrom_genome[$i][$j][2]][$h];

 #	print "$hgGenome_refGene[$hgGenome_refGene_chrom_genome[$i][$j][2]][$h]	";
      }
 #      print "\n";
    }
  }

   undef(@hgGenome_refGene);
   @hgGenome_refGene=@temp_hgGenome_refGene;
   undef(@temp_hgGenome_refGene);

   #output all data

   #for (my $i=0; $i<=$hgGenome_refGene_row_number; $i++) {
    #   for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
	 #  print "$hgGenome_refGene[$i][$j]	";
       #}
       #print "\n";
   #}

   #remove the redundant gene symbols
   #remove_redundant_genesymbols_byExonCount();

   @hgGenome_refGene_pointer=();
   my $start_row_index=1;
   my $end_row_index=1;
   my $chromosome=$hgGenome_refGene[$start_row_index][3]; #start chromosome

   my $moving_pointer=1;
   while ($moving_pointer<=$hgGenome_refGene_row_number) {
       while (($moving_pointer<$hgGenome_refGene_row_number) &&
	      ($hgGenome_refGene[$moving_pointer+1][3]==$chromosome)) {
	   $moving_pointer++;
      }

       #start to assign blat the ordinal number
       $end_row_index=$moving_pointer;

       #save the chromosome index
       $hgGenome_refGene_pointer[$chromosome][1]=$start_row_index;
       $hgGenome_refGene_pointer[$chromosome][2]=$end_row_index;

       #re-start the pointer;
       $moving_pointer++;
       $start_row_index=$moving_pointer;
       $end_row_index=$moving_pointer;

       if ($start_row_index<=$hgGenome_refGene_row_number) {
	   $chromosome=$hgGenome_refGene[$start_row_index][3]; #reset the chromosome
      }
  } #while ($moving_pointer<=$hgGenome_refGene_row_number)

   #for (my $i=1; $i<=25; $i++) {
    #   print "$hgGenome_refGene_pointer[$i][1]	$hgGenome_refGene_pointer[$i][2]\n";
   #}

}

sub remove_redundant_genesymbols_byExonCount_GeneSize {

     my @temp_exon_hgGenome_refGene=();
     my @temp_size_hgGenome_refGene=();

     #first header line
     for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
	 $temp_exon_hgGenome_refGene[0][$j]=$hgGenome_refGene[0][$j];
	 #   print "$hgGenome_refGene[0][$j]	";
    }
     # print "\n";    

     my %genesymbol_indicator=();

     #output data content
     my $row_index=0;
     #my $exonCount=0;
     #my $geneSize=0;
     for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
	 if (!exists($genesymbol_indicator{uc($hgGenome_refGene[$i][7])})) {	    
	     $row_index++;
	     for (my $h=0; $h<=$hgGenome_refGene_column_number; $h++) {
		 $temp_exon_hgGenome_refGene[$row_index][$h]=$hgGenome_refGene[$i][$h];
		 $temp_size_hgGenome_refGene[$row_index][$h]=$hgGenome_refGene[$i][$h];		
	    }
	     $genesymbol_indicator{uc($hgGenome_refGene[$i][7])}=$i;
	     #$exonCount=$hgGenome_refGene[$i][1];
	     #$geneSize=$hgGenome_refGene[$i][6]-$hgGenome_refGene[$i][5]; #gene size
	}

	 else {
	     #filtering by gene genomic size
	     if (($temp_size_hgGenome_refGene[$row_index][6]-$temp_size_hgGenome_refGene[$row_index][5])<
		 ($hgGenome_refGene[$i][6]-$hgGenome_refGene[$i][5])) {
	       for (my $h=0; $h<=$hgGenome_refGene_column_number; $h++) {
		     $temp_size_hgGenome_refGene[$row_index][$h]=$hgGenome_refGene[$i][$h];
		}
		 #$genesymbol_indicator{uc($hgGenome_refGene[$i][7])}=$i;
	    }

	     #exonCount and geneSize
	     if ($temp_exon_hgGenome_refGene[$row_index][1]<$hgGenome_refGene[$i][1]) {
		 for (my $h=0; $h<=$hgGenome_refGene_column_number; $h++) {
		     $temp_exon_hgGenome_refGene[$row_index][$h]=$hgGenome_refGene[$i][$h];
		}
		 #$genesymbol_indicator{uc($hgGenome_refGene[$i][7])}=$i;
	    }
	}		
    } #for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++)

     #replace with much larger size refSeg gene (105% larger)
     for (my $i=1; $i<=$row_index; $i++) {

	 #start to deal cDNA size
	 #for exon-based
	 my $exonlocation=$temp_exon_hgGenome_refGene[$i][2];
	 my $exonCount=$temp_exon_hgGenome_refGene[$i][1];
	 $exonlocation=~s/\|\|//;
	 my (@temp_exonlocation)=split(",", $exonlocation);
	 my $exon_based_genesize=0;
	 for (my $h=1; $h<=$exonCount; $h++) {
	     $exon_based_genesize=$exon_based_genesize+$temp_exonlocation[$h-1+$exonCount]-
		 $temp_exonlocation[$h-1]+1;
	}

	 #for size-based
	 $exonlocation=$temp_size_hgGenome_refGene[$i][2];
	 $exonCount=$temp_size_hgGenome_refGene[$i][1];
	 $exonlocation=~s/\|\|//;
	 @temp_exonlocation=split(",", $exonlocation);
	 my $size_based_genesize=0;
	 for (my $h=1; $h<=$exonCount; $h++) {
	     $size_based_genesize=$size_based_genesize+$temp_exonlocation[$h-1+$exonCount]-
		 $temp_exonlocation[$h-1]+1;
	}

	 #my $exon_based_genesize=$temp_exon_hgGenome_refGene[$i][6]-
	 #    $temp_exon_hgGenome_refGene[$i][5];
	 #my $size_based_genesize=$temp_size_hgGenome_refGene[$i][6]-
	 #    $temp_size_hgGenome_refGene[$i][5];

	 my $genomic_exon_based_genesize=$temp_exon_hgGenome_refGene[$i][6]-
	     $temp_exon_hgGenome_refGene[$i][5];
	 my $genomic_size_based_genesize=$temp_size_hgGenome_refGene[$i][6]-
	     $temp_size_hgGenome_refGene[$i][5];

	 if ($size_based_genesize>(1.05*$exon_based_genesize)) {
	     for (my $h=0; $h<=$hgGenome_refGene_column_number; $h++) {
		 $temp_exon_hgGenome_refGene[$i][$h]=$temp_size_hgGenome_refGene[$i][$h];
	    }	    
	}
	 #same exonCount and larger size
	 elsif (($temp_exon_hgGenome_refGene[$i][1]==$temp_size_hgGenome_refGene[$i][1]) &&
		($genomic_size_based_genesize>$genomic_exon_based_genesize)) {
		#($size_based_genesize>$exon_based_genesize)) {
	     for (my $h=0; $h<=$hgGenome_refGene_column_number; $h++) {
		 $temp_exon_hgGenome_refGene[$i][$h]=$temp_size_hgGenome_refGene[$i][$h];
	    }
	}
    }

     undef(@hgGenome_refGene);
     @hgGenome_refGene=@temp_exon_hgGenome_refGene;
     undef(@temp_exon_hgGenome_refGene);
     undef(@temp_size_hgGenome_refGene);
     #re-set the number of genes
     $hgGenome_refGene_row_number=$row_index;


     %NM_index=();
     for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
	 $NM_index{$hgGenome_refGene[$i][0]}=$hgGenome_refGene[$i][3];
    }

     ##output all data
     #for (my $i=0; $i<=$hgGenome_refGene_row_number; $i++) {
	 #for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
	  #   print "$hgGenome_refGene[$i][$j]	";
	 #}
	 #print "\n";
     #}

}

sub sort_hgGenome_refGene_by_geneName {


     #output the data
     #for (my $i=0; $i<=$hgGenome_refGene_row_number; $i++) {
     #  for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
     #    print "$hgGenome_refGene[$i][$j]	";
     # }
     #  print "\n";
     #}


     #sort the data by anchor1 geneName
     @sort_hgGenome_index=();
     for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
	 $sort_hgGenome_index[$i][1]=$hgGenome_refGene[$i][7];
	 $sort_hgGenome_index[$i][2]=$i;
    }


     quicksort_hgGenome_index_by_geneName(1, $hgGenome_refGene_row_number);

     #re-order the data by geneName
     my @temp_sort_hgGenome_refGene=();
     for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
	 my $row_index=$sort_hgGenome_index[$i][2];
	 for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
	     $temp_sort_hgGenome_refGene[$i][$j]=$hgGenome_refGene[$row_index][$j];
	}
    }

     #header
     for (my $j=0; $j<=$hgGenome_refGene_column_number; $j++) {
	 $temp_sort_hgGenome_refGene[0][$j]=$hgGenome_refGene[0][$j];
    }

     undef(@hgGenome_refGene);
     @hgGenome_refGene=();
     @hgGenome_refGene=@temp_sort_hgGenome_refGene;  
     undef(@temp_sort_hgGenome_refGene);

}

sub quicksort_hgGenome_index_by_geneName {
     my ($left, $right)=@_;
     my ($pointer_left, $pointer_right);
     my ($pivot, $temp_1, $temp_2);

     if ($right>=$left) {
	 $pivot=$sort_hgGenome_index[$right][1];
	 # print "pivot= $pivot\n";
	 $pointer_left=$left-1;
	 $pointer_right=$right;
	 do
	 {
	     do
	     {
		 $pointer_left=$pointer_left+1;
	    } until ($sort_hgGenome_index[$pointer_left][1] ge $pivot);

	     do
	     { 
		 $pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_hgGenome_index[$pointer_right][1] le $pivot));

	     $temp_1=$sort_hgGenome_index[$pointer_left][1];
	     $sort_hgGenome_index[$pointer_left][1]=$sort_hgGenome_index[$pointer_right][1];
	     $sort_hgGenome_index[$pointer_right][1]=$temp_1;

	     $temp_2=$sort_hgGenome_index[$pointer_left][2];
	     $sort_hgGenome_index[$pointer_left][2]=$sort_hgGenome_index[$pointer_right][2];
	     $sort_hgGenome_index[$pointer_right][2]=$temp_2;

	} until ($pointer_right<=$pointer_left);

	 $sort_hgGenome_index[$pointer_right][1]=$sort_hgGenome_index[$pointer_left][1];
	 $sort_hgGenome_index[$pointer_left][1]=$sort_hgGenome_index[$right][1];
	 $sort_hgGenome_index[$right][1]=$temp_1;

	 $sort_hgGenome_index[$pointer_right][2]=$sort_hgGenome_index[$pointer_left][2];
	 $sort_hgGenome_index[$pointer_left][2]=$sort_hgGenome_index[$right][2];
	 $sort_hgGenome_index[$right][2]=$temp_2;


	 quicksort_hgGenome_index_by_geneName($left, $pointer_left-1);
	 quicksort_hgGenome_index_by_geneName($pointer_left+1, $right);
    }
}

sub generate_cDNA_order {
    #add the depth
     #%cDNA_depth=();
     %cDNA_order=();
     for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
	 #$cDNA_depth{$hgGenome_refGene[$i][7]}=$hgGenome_refGene[$i][$hgGenome_refGene_column_number];
	 $cDNA_order{$hgGenome_refGene[$i][7]}=$i;
    }
}

sub output_fasta_for_MosaikBed_from_aggregate_reduced_annotate_cross_data_after_blat {    
     #unlink the file
     my $temp_fa_file=$ARGV[0];
     #$temp_fa_file=~s/Chip/Blat/;
     $temp_fa_file=~s/\.txt/\.fa/;
     my $temp_psl_file=$temp_fa_file;
     $temp_psl_file=~s/\.fa/\.psl/;

     unlink($temp_fa_file);
     unlink($temp_psl_file);

     #output the data
     my $temp_file=$ARGV[0];
     $temp_file=~s/Chip/BlatBed/;
     $temp_file=~s/\.txt/\.fa/;
     open (OUTPUT, ">$temp_file")
	 || die "Can't open $temp_file $!";

     #generate the candidate index
     my %candidate_index=();
     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 #comment Within uses gene relative distance, but Cross doesn't
	 #if (abs($aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+8])>=2) {
	 my $gene_pair=$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1].
	     "\_".$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2];
	 $gene_pair=~s/ //g;
	 $candidate_index{$gene_pair}=1;
	 #}
    }

     #start to output the fasta data
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 my $gene_pair=$reduced_cross_data[$i][$reduced_cross_column_number+1].
	     "\_".$reduced_cross_data[$i][$reduced_cross_column_number+2];
	 $gene_pair=~s/ //g;
	 if (exists($candidate_index{$gene_pair})) {
	     #if ($reduced_cross_data[$i][$reduced_cross_column_number+5]=~/Yes/) {
	     if (($reduced_cross_data[$i][3]=~/\>/) &&
		 ($reduced_cross_data[$i][4]=~/A|C|G|T|a|c|g|t/)) {
		 print OUTPUT "$reduced_cross_data[$i][3]\n";
		 print OUTPUT "$reduced_cross_data[$i][4]\n";
	    }

	     if (($reduced_cross_data[$i][8]=~/\>/) && 
		 ($reduced_cross_data[$i][9]=~/A|C|G|T|a|c|g|t/)) {
		 print OUTPUT "$reduced_cross_data[$i][8]\n";
		 print OUTPUT "$reduced_cross_data[$i][9]\n";
	    }
	}
    }    
     close(OUTPUT);
}

sub extract_sequence_for_predicted_fusions {
     #my %chr_indicator=();
     #upload the reference sequence
     #upload_sequence_data();
     my $extract_seq_indicator="No";

     #for within, it has to be [+12]
     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+11]="insilico-sequence";

     $aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+12]="Partner_Order";
     #$aggregate_blat_cross_data[0][$aggregate_blat_cross_column_number+14]="NA";

     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 my $chr_1=$aggregate_blat_cross_data[$i][0];
	 my $chr_2=$aggregate_blat_cross_data[$i][5];
	 if ($chr_1=~/X/) { $chr_1=23;}
	 if ($chr_1=~/Y/) { $chr_1=24;}
	 if ($chr_1=~/M/) { $chr_1=25;}

	 if ($chr_2=~/X/) { $chr_2=23;}
	 if ($chr_2=~/Y/) { $chr_2=24;}
	 if ($chr_2=~/M/) { $chr_2=25;}

	 #$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]="NA";
	 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]="NA";
	 #$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+14]="NA";

	 if ((exists($cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]})) &&
	     (exists($cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]}))) {
	     my $ori_index1=$cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]};
	     my $ori_index2=$cDNA_order{$aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]};
	     my $ori_gene1=$hgGenome_refGene[$ori_index1][4]; #GENE STRAND
	     my $ori_gene2=$hgGenome_refGene[$ori_index2][4]; #gene strand

	     my @sequence1_loc=split("\/\/", $aggregate_blat_cross_data[$i][2]);
	     my @sequence2_loc=split("\/\/", $aggregate_blat_cross_data[$i][7]);

	     my (@read_seq)=split("\/\/", $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]);
	     my $left_seq_mate_1=$read_seq[0];
	     my $right_seq_mate_1=$read_seq[1];
	     my $left_seq_mate_2=$read_seq[2];
	     my $right_seq_mate_2=$read_seq[3];

	     #F vs R
	     if (($aggregate_blat_cross_data[$i][2]=~/\+/) && ($aggregate_blat_cross_data[$i][7]=~/\-/)) {
		 #predicted seq
		 my $short_seq=$right_seq_mate_1."\|\|".generate_complementary_sequence($left_seq_mate_2);
		 my $long_seq=$left_seq_mate_1."\|\|".generate_complementary_sequence($right_seq_mate_2);
		 my $final_seq="short_seq\=$short_seq\/\/\/long_seq=$long_seq";
		 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]=$final_seq;

		 #partner order
		 if (($ori_gene1=~/\+/) && ($ori_gene2=~/\+/)) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]."\_".
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2];		    
		}	    
		 elsif  (($ori_gene1=~/\-/) && ($ori_gene2=~/\-/)) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]."\_".
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1];
		}

		 if ($extract_seq_indicator=~/Yes/) {
		     #in-silico sequence
		     my $before_start_loc=$sequence1_loc[1];
		     my $before_end_loc=$sequence1_loc[2];

		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=$before_start_loc;
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+14]=$before_end_loc;

		     my @temp_before_seq=();
		     upload_sequence_data($chr_1);
		     my $before_length=$before_end_loc-$before_start_loc;
		     for (my $j=0; $j<=$before_length; $j++) {
			 $temp_before_seq[$j]=$chr_data[$chr_1][$j+$before_start_loc];
			 $temp_before_seq[$j]=$chr_data[$j+$before_start_loc];
		    }
		     #print "\n";

		     my $before_sequence=join("", @temp_before_seq[0..$before_length]);


		     my $after_start_loc=$sequence2_loc[1];
		     my $after_end_loc=$sequence2_loc[2];

		     my @temp_after_seq=();
		     if ($chr_2 ne $chr_1) { upload_sequence_data($chr_2);}
		     my $after_length=$after_end_loc-$after_start_loc;
		     for (my $j=0; $j<=$after_length; $j++) {
			 #$temp_after_seq[$j]=$chr_data[$chr_2][$j+$after_start_loc];
			 $temp_after_seq[$j]=$chr_data[$j+$after_start_loc];
		    }
		     my $after_sequence=join("", @temp_after_seq[0..$after_length]);
		     #my $after_sequence=join("", @chr_data[$chr_2][$after_start_loc..$after_end_loc]);

		     #in-silico sequence
		     my $temp_sequence=$before_sequence."\/\/".$after_sequence;
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]=$temp_sequence;
		}
	    }

	     #R vs F and swap R with F
	     elsif (($aggregate_blat_cross_data[$i][2]=~/\-/) && ($aggregate_blat_cross_data[$i][7]=~/\+/)) {
		 #predicted seq
		 my $short_seq=$right_seq_mate_2."\|\|".generate_complementary_sequence($left_seq_mate_1);
		 my $long_seq=$left_seq_mate_2."\|\|".generate_complementary_sequence($right_seq_mate_1);
		 my $final_seq="short_seq\=$short_seq\/\/\/long_seq=$long_seq";
		 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]=$final_seq;

		 #partner order
		 if (($ori_gene1=~/\+/) && ($ori_gene2=~/\+/)) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]."\_".
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1];
		}	    
		 elsif  (($ori_gene1=~/\-/) && ($ori_gene2=~/\-/)) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]."\_".
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2];
		}

		 #in-silico sequence
		 if ($extract_seq_indicator=~/Yes/) {
		     my $before_start_loc=$sequence2_loc[1];
		     my $before_end_loc=$sequence2_loc[2];
		     my @temp_before_seq=();
		     upload_sequence_data($chr_2);
		     my $before_length=$before_end_loc-$before_start_loc;
		     for (my $j=0; $j<=$before_length; $j++) {
			 #$temp_before_seq[$j]=$chr_data[$chr_2][$j+$before_start_loc];
			 $temp_before_seq[$j]=$chr_data[$j+$before_start_loc];
		    }
		     my $before_sequence=join("", @temp_before_seq[0..$before_length]);
		     #my $before_sequence=join("", @chr_data[$chr_2][$before_start_loc..$before_end_loc]);

		     my $after_start_loc=$sequence1_loc[1];
		     my $after_end_loc=$sequence1_loc[2];
		     my @temp_after_seq=();
		     if ($chr_1 ne $chr_2) { upload_sequence_data($chr_1);}
		     my $after_length=$after_end_loc-$after_start_loc;
		     for (my $j=0; $j<=$after_length; $j++) {
			 #$temp_after_seq[$j]=$chr_data[$chr_1][$j+$after_start_loc];
			 $temp_after_seq[$j]=$chr_data[$j+$after_start_loc];
		    }
		     my $after_sequence=join("", @temp_after_seq[0..$after_length]);
		     #my $after_sequence=join("", @chr_adat[$chr_1][$after_start_loc..$after_end_loc]);		
		     #in-silico sequence
		     my $temp_sequence=$before_sequence."\/\/".$after_sequence;
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]=$temp_sequence;
		}
	    }

	     #F vs F and reverse second F
	     elsif  (($aggregate_blat_cross_data[$i][2]=~/\+/) && ($aggregate_blat_cross_data[$i][7]=~/\+/)) {
		 #predicted seq ????
		 my $short_seq=$right_seq_mate_1."\|\|".generate_complementary_sequence($right_seq_mate_2);
		 my $long_seq=$left_seq_mate_1."\|\|".generate_complementary_sequence($left_seq_mate_2);
		 my $final_seq="short_seq\=$short_seq\/\/\/long_seq=$long_seq";
		 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]=$final_seq;

		 #partner order
		 if (($ori_gene1=~/\+/) && ($ori_gene2=~/\-/)) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]."\_".
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2];
		}	    
		 elsif  (($ori_gene1=~/\-/) && ($ori_gene2=~/\+/)) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]."\_".
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1];
		}

		 #in-silico sequence
		 if ($extract_seq_indicator=~/Yes/) {
		     my $before_start_loc=$sequence1_loc[1];
		     my $before_end_loc=$sequence1_loc[2];
		     my @temp_before_seq=();
		     upload_sequence_data($chr_1);
		     my $before_length=$before_end_loc-$before_start_loc;
		     for (my $j=0; $j<=$before_length; $j++) {
			 #$temp_before_seq[$j]=$chr_data[$chr_1][$j+$before_start_loc];
			 $temp_before_seq[$j]=$chr_data[$j+$before_start_loc];
		    }
		     my $before_sequence=join("", @temp_before_seq[0..$before_length]);

		     #my $before_sequence=join("", @chr_data[$chr_1][$before_start_loc..$before_end_loc]);

		     my $after_start_loc=$sequence2_loc[1];
		     my $after_end_loc=$sequence2_loc[2];
		     my @temp_after_seq=();
		     if ($chr_2 ne $chr_1) { upload_sequence_data($chr_2);}
		     my $after_length=$after_end_loc-$after_start_loc;
		     for (my $j=0; $j<=$after_length; $j++) {
			 #$temp_after_seq[$j]=$chr_data[$chr_2][$j+$after_start_loc];
			 $temp_after_seq[$j]=$chr_data[$j+$after_start_loc];
		    }
		     my $after_sequence=join("", @temp_after_seq[0..$after_length]);
		     #my $after_sequence=join("", @chr_data[$chr_2][$after_start_loc..$after_end_loc]);
		     #reverse the sequence
		     $after_sequence=generate_complementary_sequence($after_sequence);

		     #in-silico sequence
		     my $temp_sequence=$before_sequence."\/\/".$after_sequence;
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]=$temp_sequence;
		}
	    }

	     #R vs R and reverse first R
	     elsif  (($aggregate_blat_cross_data[$i][2]=~/\-/) && ($aggregate_blat_cross_data[$i][7]=~/\-/)) {
		 #predicted seq (complementary is different from alignment orientation
		 ##complementary both, and then reverse teh first
		 $left_seq_mate_1=generate_complementary_sequence($left_seq_mate_1);
		 $right_seq_mate_1=generate_complementary_sequence($right_seq_mate_1);
		 $left_seq_mate_2=generate_complementary_sequence($left_seq_mate_2);
		 $left_seq_mate_2=generate_complementary_sequence($left_seq_mate_2);

		 my $short_seq=generate_complementary_sequence($left_seq_mate_1)."\|\|".$left_seq_mate_2;
		 my $long_seq=generate_complementary_sequence($right_seq_mate_1)."\|\|".$right_seq_mate_2;
		 my $final_seq="short_seq\=$short_seq\/\/\/long_seq\=$long_seq";
		 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]=$final_seq;

		 #partner order
		 if (($ori_gene1=~/\+/) && ($ori_gene2=~/\-/)) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]."\_".
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1];
		}	    
		 elsif  (($ori_gene1=~/\-/) && ($ori_gene2=~/\+/)) {
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+12]=
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]."\_".
			 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2];
		}

		 #in-silico sequence
		 if ($extract_seq_indicator=~/Yes/) {
		     my $before_start_loc=$sequence1_loc[1];
		     my $before_end_loc=$sequence1_loc[2];
		     my @temp_before_seq=();
		     upload_sequence_data($chr_1);
		     my $before_length=$before_end_loc-$before_start_loc;
		     for (my $j=0; $j<=$before_length; $j++) {
			 #$temp_before_seq[$j]=$chr_data[$chr_1][$j+$before_start_loc];
			 $temp_before_seq[$j]=$chr_data[$j+$before_start_loc];
		    }
		     my $before_sequence=join("", @temp_before_seq[0..$before_length]);
		     #my $before_sequence=join("", @chr_data[$chr_1][$before_start_loc..$before_end_loc]);

		     #reverse the sequence
		     $before_sequence=generate_complementary_sequence($before_sequence);

		     my $after_start_loc=$sequence2_loc[1];
		     my $after_end_loc=$sequence2_loc[2];
		     my @temp_after_seq=();
		     if ($chr_2 ne $chr_1) { upload_sequence_data($chr_2);}
		     my $after_length=$after_end_loc-$after_start_loc;
		     for (my $j=0; $j<=$after_length; $j++) {
			 #$temp_after_seq[$j]=$chr_data[$chr_2][$j+$after_start_loc];
			 $temp_after_seq[$j]=$chr_data[$j+$after_start_loc];
		    }
		     my $after_sequence=join("", @temp_after_seq[0..$after_length]);
		     #my $after_sequence=join("", @chr_data[$chr_2][$after_start_loc..$after_end_loc]);
		     #in-silico sequence
		     my $temp_sequence=$before_sequence."\/\/".$after_sequence;
		     $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+11]=$temp_sequence;
		}
	    }

	}

    }

}

sub upload_sequence_data {
     my ($chr_number)=@_;
     @chr_data=();
     undef(@chr_data);
     @chr_data=();

     #for (my $chr=1; $chr<=25; $chr++) {

     #my ($chr_number)=$chr;
     my $temp_hg19_mrna="/RIS/home/xsu1/Mosaik/hg19_Mosaik/hg19_mrna/RefSeq_Processor/Fasta_generator/";
     my $in_data_file=$temp_hg19_mrna."chr".$chr_number."\.fa"; #it's not \.fasta";

     #my $temp_hg19_hgGenomerandom="/RIS/home/xsu1/Mosaik/hg19_Mosaik/hg19_hgGenomerandom/";
     #my $in_data_file=$temp_hg19_hgGenomerandom."chr".$chr_number."\.fa"; #it's not \.fasta";

     if ($in_data_file=~/chr23/) { $in_data_file=~s/chr23/chrX/;}
     elsif ($in_data_file=~/chr24/) { $in_data_file=~s/chr24/chrY/;}
     elsif ($in_data_file=~/chr25/) { $in_data_file=~s/chr25/chrM/;}
     
     open (INPUT, "<$in_data_file")
	 || die "Can't open $in_data_file $!";

     $chr_data_row_number=0;
     while (<INPUT>) {
	 $text=$_;
	 $text=~s/\n//g;
	 #$line_indicator++;
	 if (!($text=~/chr/)) {
	     my (@temp_letter)=split("", $text);
	     for (my $i=0; $i<=$#temp_letter; $i++) {
		 $chr_data_row_number++;
		 $chr_data[$chr_data_row_number]=$temp_letter[$i];
	    }
	}
	 #if ($line_indicator>5) { last;}
    }
     close(INPUT);
     #}
} 

sub generate_reverse_sequence {
     my ($temp_sequence)=@_;

     if ($temp_sequence=~/A|C|G|T|a|c|g|t/) {
	 my %basepair_conversion=();
	 $basepair_conversion{"A"}="T";
	 $basepair_conversion{"T"}="A";
	 $basepair_conversion{"G"}="C";
	 $basepair_conversion{"C"}="G";
	 $basepair_conversion{"N"}="N";
	 $basepair_conversion{"X"}="X";

	 my @base_pair=split("", $temp_sequence);
	 my $sequence_length=$#base_pair;

	 #for (my $i=0; $i<=$sequence_length; $i++) {
	 #    $base_pair[$i]=$basepair_conversion{uc($base_pair[$i])};
	 #}

	 my @converted_base=();
	 for (my $i=0; $i<=$sequence_length; $i++) {
	     $converted_base[$i]=$base_pair[$sequence_length-$i];
	}

	 my $true_sequence=join("", @converted_base[0..$sequence_length]);

	 return($true_sequence);
    }
     else {
	 return($temp_sequence);
    }

}

sub generate_complementary_sequence {
     my ($temp_sequence)=@_;

     if ($temp_sequence=~/A|C|G|T|a|c|g|t/) {
	 my %basepair_conversion=();
	 $basepair_conversion{"A"}="T";
	 $basepair_conversion{"T"}="A";
	 $basepair_conversion{"G"}="C";
	 $basepair_conversion{"C"}="G";
	 $basepair_conversion{"N"}="N";
	 $basepair_conversion{"X"}="X";

	 my @base_pair=split("", $temp_sequence);
	 my $sequence_length=$#base_pair;

	 for (my $i=0; $i<=$sequence_length; $i++) {
	     $base_pair[$i]=$basepair_conversion{uc($base_pair[$i])};
	}

	 my @converted_base=();
	 for (my $i=0; $i<=$sequence_length; $i++) {
	     $converted_base[$i]=$base_pair[$sequence_length-$i];
	}

	 my $true_sequence=join("", @converted_base[0..$sequence_length]);

	 return($true_sequence);
    }
     else {
	 return($temp_sequence);
    }

}

 ##for reduced cross data

sub add_exon_location_to_gene_for_each_read {
     #index hgGenome by gene name
     index_hgGenome_by_genename();
     my $aligned_read_length=60; #read_length=60bp

     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 my $gene0_name=$reduced_cross_data[$i][$reduced_cross_column_number+1];
	 if (exists($genename_to_hgGenome_index{$gene0_name})) {
	     my $gene0_index=$genename_to_hgGenome_index{$gene0_name};
	     my $gene0_location=$reduced_cross_data[$i][1];
	     my $first_location=detect_location_site_within_gene_byLocation($gene0_index, $gene0_location);
	     my $second_location=detect_location_site_within_gene_byLocation($gene0_index, $gene0_location+$aligned_read_length);
	     my $third_location=detect_location_site_within_gene_byLocation($gene0_index, $gene0_location-$aligned_read_length);
	     if ($first_location=~/exon/) {
		 $reduced_cross_data[$i][$reduced_cross_column_number+1]=
		     $reduced_cross_data[$i][$reduced_cross_column_number+1]."\/\/".$first_location;
	    }
	     elsif ($second_location=~/exon/) {
		 $reduced_cross_data[$i][$reduced_cross_column_number+1]=
		     $reduced_cross_data[$i][$reduced_cross_column_number+1]."\/\/".$second_location;
	    }
	     else {
		 $reduced_cross_data[$i][$reduced_cross_column_number+1]=
		     $reduced_cross_data[$i][$reduced_cross_column_number+1]."\/\/".$third_location;
	    }		

	}

	 my $gene1_name=$reduced_cross_data[$i][$reduced_cross_column_number+2];
	 if (exists($genename_to_hgGenome_index{$gene1_name})) {
	     my $gene1_index=$genename_to_hgGenome_index{$gene1_name};
	     my $gene1_location=$reduced_cross_data[$i][6];
	     my $first_location=detect_location_site_within_gene_byLocation($gene1_index, $gene1_location);
	     my $second_location=detect_location_site_within_gene_byLocation($gene1_index, $gene1_location+$aligned_read_length);
	     my $third_location=detect_location_site_within_gene_byLocation($gene1_index, $gene1_location-$aligned_read_length);
	     if ($first_location=~/exon/) {
		 $reduced_cross_data[$i][$reduced_cross_column_number+2]=
		     $reduced_cross_data[$i][$reduced_cross_column_number+2]."\/\/".$first_location;
	    }
	     elsif ($second_location=~/exon/) {
		 $reduced_cross_data[$i][$reduced_cross_column_number+2]=
		     $reduced_cross_data[$i][$reduced_cross_column_number+2]."\/\/".$second_location;
	    }
	     else {
		 $reduced_cross_data[$i][$reduced_cross_column_number+2]=
		     $reduced_cross_data[$i][$reduced_cross_column_number+2]."\/\/".$third_location;
	    }	
	}
    } 

}

sub build_gene_exon_pair_reduced_cross_data {
     @gene_exon_reduced_cross_data=();
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 my (@temp_name0)=split("\/\/", $reduced_cross_data[$i][$reduced_cross_column_number+1]);
	 $gene_exon_reduced_cross_data[$i][1]=$temp_name0[0]; #gene name
	 $gene_exon_reduced_cross_data[$i][3]=$temp_name0[1]; #exon
	 my (@temp_name1)=split("\/\/", $reduced_cross_data[$i][$reduced_cross_column_number+2]);
	 $gene_exon_reduced_cross_data[$i][2]=$temp_name1[0]; #gene name
	 $gene_exon_reduced_cross_data[$i][4]=$temp_name1[1]; #exon
    }    
}

sub remove_exon_from_gene_reduced_cross_data {
     #remove the exon location
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 $reduced_cross_data[$i][$reduced_cross_column_number+1]=$gene_exon_reduced_cross_data[$i][1];
	 $reduced_cross_data[$i][$reduced_cross_column_number+2]=$gene_exon_reduced_cross_data[$i][2];
    }
}

sub add_exon_backto_gene_reduced_cross_data {
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 $reduced_cross_data[$i][$reduced_cross_column_number+1]=$gene_exon_reduced_cross_data[$i][1].
	     "\/\/".$gene_exon_reduced_cross_data[$i][3];
	 $reduced_cross_data[$i][$reduced_cross_column_number+2]=$gene_exon_reduced_cross_data[$i][2].
	     "\/\/".$gene_exon_reduced_cross_data[$i][4];
    }
}

 ##aggregate blat cross data
sub build_gene_exon_pair_aggregate_blat_cross_data {
     @gene_exon_aggregate_blat_cross_data=();
     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 my (@temp_name0)=split("\/\/", $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]);
	 $gene_exon_aggregate_blat_cross_data[$i][1]=$temp_name0[0]; #gene name
	 $gene_exon_aggregate_blat_cross_data[$i][3]=$temp_name0[1]; #exon
	 my (@temp_name1)=split("\/\/", $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]);
	 $gene_exon_aggregate_blat_cross_data[$i][2]=$temp_name1[0]; #gene name
	 $gene_exon_aggregate_blat_cross_data[$i][4]=$temp_name1[1]; #exon
    }    
}

sub remove_exon_from_gene_aggregate_blat_cross_data {
     #remove the exon location
     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]=$gene_exon_aggregate_blat_cross_data[$i][1];
	 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]=$gene_exon_aggregate_blat_cross_data[$i][2];
    }
}

sub add_exon_backto_gene_aggregate_blat_cross_data {
     for (my $i=1; $i<=$aggregate_blat_cross_row_number; $i++) {
	 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+1]=$gene_exon_aggregate_blat_cross_data[$i][1].
	     "\/\/".$gene_exon_aggregate_blat_cross_data[$i][3];
	 $aggregate_blat_cross_data[$i][$aggregate_blat_cross_column_number+2]=$gene_exon_aggregate_blat_cross_data[$i][2]
	     ."\/\/".$gene_exon_aggregate_blat_cross_data[$i][4];
    }
}

sub resort_reduced_crossdata_by_reduced_anchor1 {

     my $temp_start_pointer=1;
     my $temp_end_pointer=1;
     my (@gene_ID)=split("\/\/", $reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1]);

     my $moving_pointer=0;
     while ($moving_pointer<$reduced_cross_row_number) {
	 while ($moving_pointer<$reduced_cross_row_number) {
	     my (@temp_gene_ID)=split("\/\/", $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+1]);
	     if ($gene_ID[0] eq $temp_gene_ID[0]) {
		 $moving_pointer++;
	    }
	     else {
		 last;
	    }
	}

	 #start to assign genename the ordinal number
	 $temp_end_pointer=$moving_pointer;


	 #re-sort the data by reduced_anchor1 genename
	 my @temp_reduced_anchorGene_data=();
	 my $temp_reduced_anchorGene_number=0;
	 for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	     $temp_reduced_anchorGene_number++;
	     for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		 $temp_reduced_anchorGene_data[$temp_reduced_anchorGene_number][$j]=$reduced_cross_data[$h][$j];
	    }
	}

	 @sort_reduced_anchorGene_index=();
	 for (my $h=1; $h<=$temp_reduced_anchorGene_number; $h++) {
	     $sort_reduced_anchorGene_index[$h][1]=$temp_reduced_anchorGene_data[$h][$reduced_cross_column_number+2];
	     $sort_reduced_anchorGene_index[$h][2]=$h;
	}
	 quicksort_by_reduced_anchor_geneName(1, $temp_reduced_anchorGene_number);

	 #put re-sorted data back into reduced_cross-reduced_anchordata
	 for (my $h=1; $h<=$temp_reduced_anchorGene_number; $h++) {
	     my $row_index=$sort_reduced_anchorGene_index[$h][2];
	     my $true_row_number=$temp_start_pointer+$h-1;
	     for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		 $reduced_cross_data[$true_row_number][$j]=$temp_reduced_anchorGene_data[$row_index][$j];
	    }
	}

	 #re-start the pointer;
	 $moving_pointer++;
	 $temp_start_pointer=$moving_pointer;
	 $temp_end_pointer=$moving_pointer;
	 if ($temp_start_pointer<=$reduced_cross_row_number) {
	     @gene_ID=split("\/\/", $reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1]);
	}	
    }

}

sub quicksort_by_reduced_anchor_geneName {
     my ($left, $right)=@_;
     my ($pointer_left, $pointer_right);
     my ($pivot, $temp_1, $temp_2);

     if ($right>=$left) {
	 $pivot=$sort_reduced_anchorGene_index[$right][1];
	 # print "pivot= $pivot\n";
	 $pointer_left=$left-1;
	 $pointer_right=$right;
	 do
	 {
	     do
	     {
		 $pointer_left=$pointer_left+1;
	    } until ($sort_reduced_anchorGene_index[$pointer_left][1] ge $pivot);

	     do
	     { 
		 $pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_reduced_anchorGene_index[$pointer_right][1] le $pivot));

	     $temp_1=$sort_reduced_anchorGene_index[$pointer_left][1];
	     $sort_reduced_anchorGene_index[$pointer_left][1]=$sort_reduced_anchorGene_index[$pointer_right][1];
	     $sort_reduced_anchorGene_index[$pointer_right][1]=$temp_1;

	     $temp_2=$sort_reduced_anchorGene_index[$pointer_left][2];
	     $sort_reduced_anchorGene_index[$pointer_left][2]=$sort_reduced_anchorGene_index[$pointer_right][2];
	     $sort_reduced_anchorGene_index[$pointer_right][2]=$temp_2;

	} until ($pointer_right<=$pointer_left);

	 $sort_reduced_anchorGene_index[$pointer_right][1]=$sort_reduced_anchorGene_index[$pointer_left][1];
	 $sort_reduced_anchorGene_index[$pointer_left][1]=$sort_reduced_anchorGene_index[$right][1];
	 $sort_reduced_anchorGene_index[$right][1]=$temp_1;

	 $sort_reduced_anchorGene_index[$pointer_right][2]=$sort_reduced_anchorGene_index[$pointer_left][2];
	 $sort_reduced_anchorGene_index[$pointer_left][2]=$sort_reduced_anchorGene_index[$right][2];
	 $sort_reduced_anchorGene_index[$right][2]=$temp_2;


	 quicksort_by_reduced_anchor_geneName($left, $pointer_left-1);
	 quicksort_by_reduced_anchor_geneName($pointer_left+1, $right);
    }
}

sub output_reduced_cross_data {
     #output the data
     $reduced_cross_data[0][$reduced_cross_column_number+1]="Gene0";
     $reduced_cross_data[0][$reduced_cross_column_number+2]="Gene1";

     for (my $i=0; $i<=$reduced_cross_row_number; $i++) {
	 for (my $j=0; $j<=($reduced_cross_column_number+1); $j++) {
	     print "$reduced_cross_data[$i][$j]\t";
	}
	 print "$reduced_cross_data[$i][$reduced_cross_column_number+2]\n";
    }

}

sub filter_reduced_cross_data {
     #output the data
     my @temp_reduced_cross_data=();
     my $temp_reduced_cross_row_number=0;

     #header
     for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
	 $temp_reduced_cross_data[0][$j]=$reduced_cross_data[0][$j];
    }

     #content
     for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	 if ((($reduced_cross_data[$i][2]=~/\+/) || ($reduced_cross_data[$i][2]=~/\-/)) &&
	     (($reduced_cross_data[$i][7]=~/\+/) || ($reduced_cross_data[$i][7]=~/\-/))) {
	     $temp_reduced_cross_row_number++;
	     for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		 $temp_reduced_cross_data[$temp_reduced_cross_row_number][$j]=$reduced_cross_data[$i][$j];
	    }
	}
    }

     undef(@reduced_cross_data);
     @reduced_cross_data=@temp_reduced_cross_data;
     $reduced_cross_row_number=$temp_reduced_cross_row_number;
     undef(@temp_reduced_cross_data);
}

sub in_silico_clustering_generator {
    #clustering by the first mate
#    master_clustering_generator("first");

    #clustering by the second mate after the first mate    
#    master_clustering_generator("second");

    master_clustering_generator();
}

sub master_clustering_generator {
    
    ##my $MSK_stats_3stdev_distance=$MSK_stats_means+3*$MSK_stats_stdev; #for in-silico fragment length
    my ($mate_indicator)=@_;
    my $temp_start_pointer=1;
    my $temp_end_pointer=1;
    my @gene_ID_0=();
#    if ($mate_indicator=~/first/) {
	@gene_ID_0=split("\/\/", $reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1]);
#    }
#    elsif ($mate_indicator=~/second/) {
#	$gene_ID_0[0]=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
#    }
    
    my (@gene_ID_1)=split("\/\/", $reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2]);
    
    my $moving_pointer=0;
    while ($moving_pointer<$reduced_cross_row_number) {
	while ($moving_pointer<$reduced_cross_row_number) {
#	    if ($mate_indicator=~/first/) {
		@temp_gene_ID_0=split("\/\/", $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+1]);
#	    }
#	    elsif ($mate_indicator=~/second/) {
#		$temp_gene_ID_0[0]=$reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+1];
#	    }
	    my (@temp_gene_ID_1)=split("\/\/", $reduced_cross_data[$moving_pointer+1][$reduced_cross_column_number+2]);
	    if (($gene_ID_0[0] eq $temp_gene_ID_0[0]) && ($gene_ID_1[0] eq $temp_gene_ID_1[0])) {
		$moving_pointer++;
	    }
	    else {
		last;
	    }
	}
	
	#start to assign genename the ordinal number
	$temp_end_pointer=$moving_pointer;
	
	#put cluster of reads into a temp array
	my %readID_to_index=();
	my @f_forward_data=();
	my $f_forward_row_number=0;
	my @f_reverse_data=();
	my $f_reverse_row_number=0;
	
	my @r_forward_data=();
	my $r_forward_row_number=0;
	my @r_reverse_data=();
	my $r_reverse_row_number=0;
	
	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	    if (($reduced_cross_data[$h][3]=~/\>/) &&
		($reduced_cross_data[$h][4]=~/A|C|G|T|a|c|g|t/) &&
		($reduced_cross_data[$h][8]=~/\>/) && 
		($reduced_cross_data[$h][9]=~/A|C|G|T|a|c|g|t/)) {
		$readID_to_index{$reduced_cross_data[$h][3]}=$h;
		$readID_to_index{$reduced_cross_data[$h][8]}=$h;
		
		#forward && forward read pairs
		if (($reduced_cross_data[$h][2]=~/\+/) && ($reduced_cross_data[$h][7]=~/\+/)) {
		    $f_forward_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$f_forward_data[$f_forward_row_number][$j]=$reduced_cross_data[$h][$j];
		    }
		}
		#forward && reverse read pairs
		elsif (($reduced_cross_data[$h][2]=~/\+/) && ($reduced_cross_data[$h][7]=~/\-/)) {
		    $f_reverse_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$f_reverse_data[$f_reverse_row_number][$j]=$reduced_cross_data[$h][$j];
		    }
		}
		#reverse && forward read pairs		
		elsif (($reduced_cross_data[$h][2]=~/\-/) && ($reduced_cross_data[$h][7]=~/\+/)) {
		    $r_forward_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$r_forward_data[$r_forward_row_number][$j]=$reduced_cross_data[$h][$j];
		    }
		}
		#reverse && reverse read pairs
		elsif (($reduced_cross_data[$h][2]=~/\-/) && ($reduced_cross_data[$h][7]=~/\-/)) {
		    $r_reverse_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$r_reverse_data[$r_reverse_row_number][$j]=$reduced_cross_data[$h][$j];
		    }
		}
	    } #real read data
	}
	
	#generate the in-silico clustering for each combinatorial
	#forward vs forward
	if ($f_forward_row_number>=1) {
	    @read_cluster_data=@f_forward_data;
	    $cluster_row_number=$f_forward_row_number;
	    #cluster_generator_for_each_combination($mate_indicator);
	    cluster_generator_for_each_combination();
	    @f_forward_data=@read_cluster_data;
	}
	#forward vs reverse
	if ($f_reverse_row_number>=1) {
	    @read_cluster_data=@f_reverse_data;
	    $cluster_row_number=$f_reverse_row_number;
	    #cluster_generator_for_each_combination($mate_indicator);
	    cluster_generator_for_each_combination();
	    @f_reverse_data=@read_cluster_data;
	}
        #reverse vs forward
	if ($r_forward_row_number>=1) {
	    @read_cluster_data=@r_forward_data;
	    $cluster_row_number=$r_forward_row_number;
	    #cluster_generator_for_each_combination($mate_indicator);
	    cluster_generator_for_each_combination();
	    @r_forward_data=@read_cluster_data;
	}
	#reverse vs reverse
	if ($r_reverse_row_number>=1) {
	    @read_cluster_data=@r_reverse_data;
	    $cluster_row_number=$r_reverse_row_number;
	    #cluster_generator_for_each_combination($mate_indicator);
	    cluster_generator_for_each_combination();
	    @r_reverse_data=@read_cluster_data;
	}
	###end of clustering

	#put clustered data back into reduced_cross_data by f-f f-r r-f r-r
	my $temp_clustered_index=$temp_start_pointer-1;
	#forward vs reverse
	for (my $h=1; $h<=$f_reverse_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$reduced_cross_data[$temp_clustered_index][$j]=$f_reverse_data[$h][$j];
	    }
	}
	#reverse vs forward
	for (my $h=1; $h<=$r_forward_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$reduced_cross_data[$temp_clustered_index][$j]=$r_forward_data[$h][$j];
	    }
	}
	#forward vs forward
	for (my $h=1; $h<=$f_forward_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$reduced_cross_data[$temp_clustered_index][$j]=$f_forward_data[$h][$j];
	    }
	}
	#reverse vs reverse
	for (my $h=1; $h<=$r_reverse_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$reduced_cross_data[$temp_clustered_index][$j]=$r_reverse_data[$h][$j];
	    }
	}
	###end of clustering

	#re-start the pointer;
	$moving_pointer++;
	$temp_start_pointer=$moving_pointer;
	$temp_end_pointer=$moving_pointer;
	if ($temp_start_pointer<=$reduced_cross_row_number) {
#	    if ($mate_indicator=~/first/) {
		@gene_ID_0=split("\/\/", $reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1]);
#	    }
#	    elsif ($mate_indicator=~/second/) {
#		$gene_ID_0[0]=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
#	    }
	    @gene_ID_1=split("\/\/", $reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2]);
	}
    }
}

sub cluster_generator_for_each_combination {
    #my (@read_cluster_data, $cluster_row_number)=@_;
#    my ($mate_indicator)=@_;
    my $temp_clustered_index=0;  #it's global in this procedure
    
    #detect the number of sampleIDs
    my @temp_sampleID=();
    my %temp_sampleID_indicator=();
    my $temp_sampleID_number=0;
    for (my $i=1; $i<=$cluster_row_number; $i++) {
	my (@temp)=split("\-", $read_cluster_data[$i][3]);
	if (!(exists($temp_sampleID_indicator{$temp[0]}))) {
	    $temp_sampleID_number++;
	    $temp_sampleID[$temp_sampleID_number]=$temp[0];
	    $temp_sampleID_indicator{$temp[0]}=$temp_sampleID_number;
	}
    }

    my @temp_read_cluster_data=();
    for (my $i=1; $i<=$cluster_row_number; $i++) {
	for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
	    $temp_read_cluster_data[$i][$j]=$read_cluster_data[$i][$j];
	}
    }

    #sample-based clustering
    for (my $k=1; $k<=$temp_sampleID_number; $k++) {
	my @exon_exon_data=();
	my $exon_exon_row_number=0;
	my @exon_intron_data=();
	my $exon_intron_row_number=0;
	
	my @intron_exon_data=();
	my $intron_exon_row_number=0;
	my @intron_intron_data=();
	my $intron_intron_row_number=0;
	
	my @else_else_data=();
	my $else_else_row_number=0;
	
	for (my $i=1; $i<=$cluster_row_number; $i++) {
	    #sample based clustering
	    if ($temp_read_cluster_data[$i][3]=~/$temp_sampleID[$k]/) {
		#exon vs exon
		if (($temp_read_cluster_data[$i][$reduced_cross_column_number+1]=~/exon/) &&
		    ($temp_read_cluster_data[$i][$reduced_cross_column_number+2]=~/exon/)) {
		    $exon_exon_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$exon_exon_data[$exon_exon_row_number][$j]=$temp_read_cluster_data[$i][$j];
		    }
		}
		#exon vs intron
		elsif (($temp_read_cluster_data[$i][$reduced_cross_column_number+1]=~/exon/) &&
		       ($temp_read_cluster_data[$i][$reduced_cross_column_number+2]=~/intron/)) {
		    $exon_intron_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$exon_intron_data[$exon_intron_row_number][$j]=$temp_read_cluster_data[$i][$j];
		    }
		}
		#intron vs exon
		elsif (($temp_read_cluster_data[$i][$reduced_cross_column_number+1]=~/intron/) &&
		       ($temp_read_cluster_data[$i][$reduced_cross_column_number+2]=~/exon/)) {
		    $intron_exon_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$intron_exon_data[$intron_exon_row_number][$j]=$temp_read_cluster_data[$i][$j];
		    }
		}
		#intron vs intron
		elsif (($temp_read_cluster_data[$i][$reduced_cross_column_number+1]=~/intron/) &&
		       ($temp_read_cluster_data[$i][$reduced_cross_column_number+2]=~/intron/)) {
		    $intron_intron_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$intron_intron_data[$intron_intron_row_number][$j]=$temp_read_cluster_data[$i][$j];
		    }
		}
		#else vs else
		else {
		    $else_else_row_number++;
		    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
			$else_else_data[$else_else_row_number][$j]=$temp_read_cluster_data[$i][$j];
		    }
		}
		
	    }
	} #sample-based clustering
	
	#exon vs exon
	if ($exon_exon_row_number>=1) {
	    #by anchor0
	    @unit_forCluster_data=@exon_exon_data;
	    $unit_forCluster_row_number=$exon_exon_row_number;
#	    if ($mate_indicator=~/first/) {
		sort_unitforClusterdata_by_anchor0_readLoc();
		unit_to_unit_clustering(1, 1);
#	    }
#	    elsif ($mate_indicator=~/second/) {	
		#by anchor1
		sort_unitforClusterdata_by_anchor1_readLoc();
		unit_to_unit_clustering(2, 6);		
		#sort unit forCluster data by anchor0 and anchor1
#	    sort_unitforClusterdata_by_anchor0_readLoc();
#	    resort_unitforClusterdata_by_anchor1_readLoc();
#	    }
	    
	    @exon_exon_data=@unit_forCluster_data;
	    undef(@unit_forCluster_data);
	    
	}
	#exon vs intron
	if ($exon_intron_row_number>=1) {
	    @unit_forCluster_data=@exon_intron_data;
	    $unit_forCluster_row_number=$exon_intron_row_number;
#	    if ($mate_indicator=~/first/) {
		sort_unitforClusterdata_by_anchor0_readLoc();
		unit_to_unit_clustering(1, 1);
#	    }
#	    elsif ($mate_indicator=~/second/) {
		sort_unitforClusterdata_by_anchor1_readLoc();
		unit_to_unit_clustering(2, 6);
		#sort unit forCluster data by anchor0 and anchor1
#	    sort_unitforClusterdata_by_anchor0_readLoc();
#	    resort_unitforClusterdata_by_anchor1_readLoc();
#	    }
	    
	    @exon_intron_data=@unit_forCluster_data;
	    undef(@unit_forCluster_data);
	}
	#intron vs exon
	if ($intron_exon_row_number>=1) {
	    @unit_forCluster_data=@intron_exon_data;
	    $unit_forCluster_row_number=$intron_exon_row_number;
#	    if ($mate_indicator=~/first/) {
		sort_unitforClusterdata_by_anchor0_readLoc();
		unit_to_unit_clustering(1, 1);
#	    }
#	    elsif ($mate_indicator=~/second/) {
		sort_unitforClusterdata_by_anchor1_readLoc();
		unit_to_unit_clustering(2, 6);		
		#sort unit forCluster data by anchor0 and anchor1
#	    sort_unitforClusterdata_by_anchor0_readLoc();
#	    resort_unitforClusterdata_by_anchor1_readLoc();
#	    }
	    
	    @intron_exon_data=@unit_forCluster_data;
	    undef(@unit_forCluster_data);
	    
	}
	#intron vs intron
	if ($intron_intron_row_number>=1) {
	    @unit_forCluster_data=@intron_intron_data;
	    $unit_forCluster_row_number=$intron_intron_row_number;
#	    if ($mate_indicator=~/first/) {
		sort_unitforClusterdata_by_anchor0_readLoc();
		unit_to_unit_clustering(1, 1);
#	    }
#	    elsif ($mate_indicator=~/second/) {
		sort_unitforClusterdata_by_anchor1_readLoc();
		unit_to_unit_clustering(2, 6);
		#sort unit forCluster data by anchor0 and anchor1
#	    sort_unitforClusterdata_by_anchor0_readLoc();
#	    resort_unitforClusterdata_by_anchor1_readLoc();
#	    }
	    @intron_intron_data=@unit_forCluster_data;
	    undef(@unit_forCluster_data);
	}
	
	#put clustered data back into reduced_cross_data by e-e e-i i-e i-i else-else
	#my $temp_clustered_index=0;   
	for (my $h=1; $h<=$exon_exon_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		if ($j<=$reduced_cross_column_number) {
		    $read_cluster_data[$temp_clustered_index][$j]=$exon_exon_data[$h][$j];
		}
		elsif ($j==$reduced_cross_column_number+1) {
		    $read_cluster_data[$temp_clustered_index][$j]=$exon_exon_data[$h][$j];
		}
		elsif ($j==$reduced_cross_column_number+2) {
		    $read_cluster_data[$temp_clustered_index][$j]=$exon_exon_data[$h][$j];
		}
	    }
	}
	#exon vs intron
	for (my $h=1; $h<=$exon_intron_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		if ($j<=$reduced_cross_column_number) {
		    $read_cluster_data[$temp_clustered_index][$j]=$exon_intron_data[$h][$j];
		}
		elsif ($j==$reduced_cross_column_number+1) {
		    $read_cluster_data[$temp_clustered_index][$j]=$exon_intron_data[$h][$j];
		}
		elsif ($j==$reduced_cross_column_number+2) {
		    $read_cluster_data[$temp_clustered_index][$j]=$exon_intron_data[$h][$j];
		}
	    }
	}
	#intron vs exon
	for (my $h=1; $h<=$intron_exon_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		if ($j<=$reduced_cross_column_number) {
		    $read_cluster_data[$temp_clustered_index][$j]=$intron_exon_data[$h][$j];
		}
		elsif ($j==$reduced_cross_column_number+1) {
		    $read_cluster_data[$temp_clustered_index][$j]=$intron_exon_data[$h][$j];
		}
		elsif ($j==$reduced_cross_column_number+2) {
		    $read_cluster_data[$temp_clustered_index][$j]=$intron_exon_data[$h][$j];
		}
	    }
	}

	#intron vs intron
	for (my $h=1; $h<=$intron_intron_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		if ($j<=$reduced_cross_column_number) {
		    $read_cluster_data[$temp_clustered_index][$j]=$intron_intron_data[$h][$j];
		}
		elsif ($j==$reduced_cross_column_number+1) {
		    $read_cluster_data[$temp_clustered_index][$j]=$intron_intron_data[$h][$j];
		}
		elsif ($j==$reduced_cross_column_number+2) {
		    $read_cluster_data[$temp_clustered_index][$j]=$intron_intron_data[$h][$j];
		}
	    }
	}
	
	#else vs else
	for (my $h=1; $h<=$else_else_row_number; $h++) {
	    $temp_clustered_index++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		if ($j<=$reduced_cross_column_number) {
		    $read_cluster_data[$temp_clustered_index][$j]=$else_else_data[$h][$j];
		}
		else {
		    $read_cluster_data[$temp_clustered_index][$j]=$else_else_data[$h][$j];
		}
	    }
	}

    } #for sample-based
    #end of clustering process
    
    #sort the first partner
    #sort the first partner after rename the cluster ID by cluster location
    sort_read_cluster_data_by_anchor0_geneExon_withoutClusterID();
    resort_read_cluster_data_by_anchor_readLoc_after_withoutClusterID("first"); 
    rename_read_cluster_data("first");

    #sort the second partner  
    resort_read_cluster_data_by_anchor1_geneExon_after_anchor0_withoutClusterID();
    resort_read_cluster_data_by_anchor_readLoc_after_withoutClusterID("second"); 
    rename_read_cluster_data("second");
}

sub rename_read_cluster_data {
    my ($mate_indicator)=@_;
    my @exon_unit=();
    my @gene_name=();

    #$MSK_stats_means=$ARGV[3]; #fragment length
    #$MSK_stats_stdev=$ARGV[4]; #fragment stdev
    #for in-silico fragment length
    #the half read length=50
    #my $MSK_stats_3stdev_minorhalfread_distance=$MSK_stats_means+3*$MSK_stats_stdev-50;
    #my $MSK_stats_3stdev_minorhalfread_distance=$MSK_stats_means+3*$MSK_stats_stdev;

    for (my $i=1; $i<=$cluster_row_number; $i++) {
	$exon_unit[$i][3]=1; #cluster ID for first cluster in exon-based unit
	if ($mate_indicator=~/first/) {
	    (@tempID)=split("\/\/", $read_cluster_data[$i][$reduced_cross_column_number+1]);
	}
	elsif ($mate_indicator=~/second/) {
	    (@tempID)=split("\/\/", $read_cluster_data[$i][$reduced_cross_column_number+2]);
	}
	
	$gene_name[$i]=$tempID[0];

	if ($tempID[1]=~/\-/) {
	    my (@tempexon)=split("\-", $tempID[1]);
	    $exon_unit[$i][1]=$tempexon[0];
	    $exon_unit[$i][2]=$tempexon[1];
	}
	else {
	    $exon_unit[$i][1]=$tempID[1];
	    $exon_unit[$i][2]=0;
	}	
    }
    
    #rename the cluster ID
    #initialize clusterID with zero
    if ($exon_unit[1][2]==0) {
	$exon_unit[1][3]=0;
    }
    #same exon same cluster ID with distance less than 100bp   
    for (my $i=2; $i<=$cluster_row_number; $i++) {
	#if ($exon_unit[$i][2]==0) {
	#$exon_unit[$i][3]=0;
	#}
	#same exon or intron and not cluster ID with zero
	#elsif ($exon_unit[$i][2]>0) {
	if ($exon_unit[$i][1] eq $exon_unit[$i-1][1]) {
	    if ($exon_unit[$i][2]==$exon_unit[$i-1][2]) {
		$exon_unit[$i][3]=$exon_unit[$i-1][3];
	    }
	    elsif ($exon_unit[$i][2]!=$exon_unit[$i-1][2]) {
		#estimate the distance between two cluster within same exon
		my ($tempDist)=0;
		if ($mate_indicator=~/first/) {
		    $tempDist=abs($read_cluster_data[$i][1]-$read_cluster_data[$i-1][1]);
		}
		elsif ($mate_indicator=~/second/) {
		    $tempDist=abs($read_cluster_data[$i][6]-$read_cluster_data[$i-1][6]);
		}
		
		#distance larger than 100bp means it's different cluster
		#if ($tempDist>$MSK_stats_3stdev_minorhalfread_distance) {
		if ($tempDist>100) {
		    $exon_unit[$i][3]=$exon_unit[$i-1][3]+1;
		}
		#elsif ($tempDist<=$MSK_stats_3stdev_minorhalfread_distance) {
		elsif ($tempDist<=100) {
		    if ($exon_unit[$i-1][3]>0) {
			$exon_unit[$i][3]=$exon_unit[$i-1][3];
		    }
		    #otherwise, do nothing
		}
	    }	  
	    #}
	    #different exon or intro
	    #do nothing since the first cluster already has ID=1
	} #eslif ($exon_unit[$i][2]>0)
    }

    #regenerate the cluster ID for all clusters
    for (my $i=1; $i<=$cluster_row_number; $i++) {
	my $clusterID=$exon_unit[$i][1]."\-".$exon_unit[$i][3];
	#if (($mate_indicator=~/first/) && ($clusterID=~/exon/)) {
	if ($mate_indicator=~/first/) {
	    $read_cluster_data[$i][$reduced_cross_column_number+1]=$gene_name[$i]."\/\/".$clusterID;
	}
	#elsif (($mate_indicator=~/second/) && ($clusterID=~/exon/)) {
	elsif ($mate_indicator=~/second/) {
	    $read_cluster_data[$i][$reduced_cross_column_number+2]=$gene_name[$i]."\/\/".$clusterID;
	}
    }
}

sub sort_read_cluster_data_by_anchor0_geneExon_withoutClusterID {
    my ($mate_indicator)=@_;
    #sort the data by anchor0 location
    @sort_unit_index=();
    for (my $i=1; $i<=$cluster_row_number; $i++) {
	#if ($mate_indicator=~/first/) {
	#$sort_unit_index[$i][1]=$read_cluster_data[$i][1];
	if ($read_cluster_data[$i][$reduced_cross_column_number+1]=~/\-/) {
	    my (@tempID)=split("\-", $read_cluster_data[$i][$reduced_cross_column_number+1]);
	    $sort_unit_index[$i][1]=$tempID[0];
	}
	else {
	    $sort_unit_index[$i][1]=$read_cluster_data[$i][$reduced_cross_column_number+1];
	}
	$sort_unit_index[$i][2]=$i;
    }

    quicksort_read_cluster_data_by_anchor0_geneExon_withoutClusterID(1, $cluster_row_number);
    
    #re-order the data by geneName
    my @temp_sort_unit_data=();
    for (my $i=1; $i<=$cluster_row_number; $i++) {
	my $row_index=$sort_unit_index[$i][2];
	for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
	    $temp_sort_unit_data[$i][$j]=$read_cluster_data[$row_index][$j];
	}
    }

    #put back the data
    for (my $i=1; $i<=$cluster_row_number; $i++) {
	for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
	    $read_cluster_data[$i][$j]=$temp_sort_unit_data[$i][$j];
	}
    }
    undef(@temp_sort_unit_data);
    
}

sub quicksort_read_cluster_data_by_anchor0_geneExon_withoutClusterID {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1, $temp_2);

    if ($right>=$left) {
	$pivot=$sort_unit_index[$right][1];
	# print "pivot= $pivot\n";
	$pointer_left=$left-1;
	$pointer_right=$right;
	do
	{
	    do
	    {
		$pointer_left=$pointer_left+1;
	    } until ($sort_unit_index[$pointer_left][1] ge $pivot);
	    
	    do
	    { 
		$pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_unit_index[$pointer_right][1] le $pivot));
	    
	    $temp_1=$sort_unit_index[$pointer_left][1];
	    $sort_unit_index[$pointer_left][1]=$sort_unit_index[$pointer_right][1];
	    $sort_unit_index[$pointer_right][1]=$temp_1;
	    
	    $temp_2=$sort_unit_index[$pointer_left][2];
	    $sort_unit_index[$pointer_left][2]=$sort_unit_index[$pointer_right][2];
	    $sort_unit_index[$pointer_right][2]=$temp_2;
	    
	} until ($pointer_right<=$pointer_left);
	
	$sort_unit_index[$pointer_right][1]=$sort_unit_index[$pointer_left][1];
	$sort_unit_index[$pointer_left][1]=$sort_unit_index[$right][1];
	$sort_unit_index[$right][1]=$temp_1;
	
	$sort_unit_index[$pointer_right][2]=$sort_unit_index[$pointer_left][2];
	$sort_unit_index[$pointer_left][2]=$sort_unit_index[$right][2];
	$sort_unit_index[$right][2]=$temp_2;
	
	
	quicksort_read_cluster_data_by_anchor0_geneExon_withoutClusterID($left, $pointer_left-1);
	quicksort_read_cluster_data_by_anchor0_geneExon_withoutClusterID($pointer_left+1, $right);
    }
}

sub resort_read_cluster_data_by_anchor_readLoc_after_withoutClusterID {
    my ($mate_indicator)=@_;
    my $read_index=1;
    if ($mate_indicator=~/second/) {
	#$mate_index=2; #always use cluster from mate1 (never use cluster from mate2
	$read_index=6;
    }
    my $temp_start_pointer=1;
    my $temp_end_pointer=1;
    if ($mate_indicator=~/first/) {
	if ($read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1]=~/\-/) {
	    my (@tempID)=split("\-", $read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1]);
	    $gene_ID=$tempID[0];
	}
	else {
	    $gene_ID=$read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1];
	}
    }
    elsif ($mate_indicator=~/second/) {
	$gene_ID=$read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1];
    }	
    
    my $moving_pointer=0;
    while ($moving_pointer<$cluster_row_number) {
	while ($moving_pointer<$cluster_row_number) {
	    if ($mate_indicator=~/first/) {
		if ($read_cluster_data[$moving_pointer+1][$reduced_cross_column_number+1]=~/\-/) {
		    my (@tempID)=split("\-", $read_cluster_data[$moving_pointer+1][$reduced_cross_column_number+1]);
		    $temp_gene_ID=$tempID[0];
		}
		else {
		    $temp_gene_ID=$read_cluster_data[$moving_pointer+1][$reduced_cross_column_number+1];
		}
	    }
	    elsif ($mate_indicator=~/second/) {
		$temp_gene_ID=$read_cluster_data[$moving_pointer+1][$reduced_cross_column_number+1];
	    }
	    
	    if ($gene_ID eq $temp_gene_ID) {
		$moving_pointer++;
	    }
	    else {
		last;
	    }
	}
	
	#start to assign genename the ordinal number
	$temp_end_pointer=$moving_pointer;
	
	#re-sort the data by reduced_anchor1 genename
	my @temp_read_cluster_data=();
	my $temp_read_cluster_number=0;
	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	    $temp_read_cluster_number++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$temp_read_cluster_data[$temp_read_cluster_number][$j]=$read_cluster_data[$h][$j];
	    }
	}
	
	@sort_read_cluster_index=();
	for (my $h=1; $h<=$temp_read_cluster_number; $h++) {
	    $sort_read_cluster_index[$h][1]=$temp_read_cluster_data[$h][$read_index];
	    $sort_read_cluster_index[$h][2]=$h;
	}
	quicksort_read_cluster_data_by_anchor_readLoc_after_withoutClusterID(1, $temp_read_cluster_number);
	
	#put re-sorted data back into reduced_cross-reduced_anchordata
	for (my $h=1; $h<=$temp_read_cluster_number; $h++) {
	    my $row_index=$sort_read_cluster_index[$h][2];
	    my $true_row_number=$temp_start_pointer+$h-1;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$read_cluster_data[$true_row_number][$j]=$temp_read_cluster_data[$row_index][$j];
	    }
	}
	
	#re-start the pointer;
	$moving_pointer++;
	$temp_start_pointer=$moving_pointer;
	$temp_end_pointer=$moving_pointer;
	if ($temp_start_pointer<=$cluster_row_number) {
	    if ($mate_indicator=~/first/) {
		if ($read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1]=~/\-/) {
		    my (@tempID)=split("\-", $read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1]);
		    $gene_ID=$tempID[0];
		}
		else {
		    $gene_ID=$read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1];
		}	
	    }
	    elsif ($mate_indicator=~/second/) {
		$gene_ID=$read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1];
	    }
	}
    }
    
}

sub quicksort_read_cluster_data_by_anchor_readLoc_after_withoutClusterID {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1, $temp_2);
    
    if ($right>=$left) {
	$pivot=$sort_read_cluster_index[$right][1];
	# print "pivot= $pivot\n";
	$pointer_left=$left-1;
	$pointer_right=$right;
	do
	{
	    do
	    {
		$pointer_left=$pointer_left+1;
	    } until ($sort_read_cluster_index[$pointer_left][1]>=$pivot);
	    
	    do
	    { 
		$pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_read_cluster_index[$pointer_right][1]<=$pivot));
	    
	    $temp_1=$sort_read_cluster_index[$pointer_left][1];
	    $sort_read_cluster_index[$pointer_left][1]=$sort_read_cluster_index[$pointer_right][1];
	    $sort_read_cluster_index[$pointer_right][1]=$temp_1;
	    
	    $temp_2=$sort_read_cluster_index[$pointer_left][2];
	    $sort_read_cluster_index[$pointer_left][2]=$sort_read_cluster_index[$pointer_right][2];
	    $sort_read_cluster_index[$pointer_right][2]=$temp_2;
	    
	} until ($pointer_right<=$pointer_left);
	
	$sort_read_cluster_index[$pointer_right][1]=$sort_read_cluster_index[$pointer_left][1];
	$sort_read_cluster_index[$pointer_left][1]=$sort_read_cluster_index[$right][1];
	$sort_read_cluster_index[$right][1]=$temp_1;
	
	$sort_read_cluster_index[$pointer_right][2]=$sort_read_cluster_index[$pointer_left][2];
	$sort_read_cluster_index[$pointer_left][2]=$sort_read_cluster_index[$right][2];
	$sort_read_cluster_index[$right][2]=$temp_2;
		
	quicksort_read_cluster_data_by_anchor_readLoc_after_withoutClusterID($left, $pointer_left-1);
	quicksort_read_cluster_data_by_anchor_readLoc_after_withoutClusterID($pointer_left+1, $right);
    }
}


sub resort_read_cluster_data_by_anchor1_geneExon_after_anchor0_withoutClusterID {
    
    my $temp_start_pointer=1;
    my $temp_end_pointer=1;

    my ($gene_ID)=$read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1];
    my $moving_pointer=0;
    while ($moving_pointer<$cluster_row_number) {
	while ($moving_pointer<$cluster_row_number) {
	    my ($temp_gene_ID)=$read_cluster_data[$moving_pointer+1][$reduced_cross_column_number+1];
	    if ($gene_ID eq $temp_gene_ID) {
		$moving_pointer++;
	    }
	    else {
		last;
	    }
	}
	
	#start to assign genename the ordinal number
	$temp_end_pointer=$moving_pointer;
	
	#re-sort the data by reduced_anchor1 genename
	my @temp_read_cluster_data=();
	my $temp_read_cluster_number=0;
	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	    $temp_read_cluster_number++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$temp_read_cluster_data[$temp_read_cluster_number][$j]=$read_cluster_data[$h][$j];
	    }
	}
	
	@sort_read_cluster_index=();
	for (my $h=1; $h<=$temp_read_cluster_number; $h++) {
	    if ($temp_read_cluster_data[$h][$reduced_cross_column_number+2]=~/\-/) {
		my (@tempID)=split("\-", $temp_read_cluster_data[$h][$reduced_cross_column_number+2]);
		$sort_read_cluster_index[$h][1]=$tempID[0];
	    }
	    else {
		$sort_read_cluster_index[$h][1]=$temp_read_cluster_data[$h][$reduced_cross_column_number+2];
	    }
	    $sort_read_cluster_index[$h][2]=$h;
	}
	
	quicksort_read_cluster_data_by_anchor1_geneExon_after_anchor0_withoutClusterID(1, $temp_read_cluster_number);
	
	#put re-sorted data back into reduced_cross-reduced_anchordata
	for (my $h=1; $h<=$temp_read_cluster_number; $h++) {
	    my $row_index=$sort_read_cluster_index[$h][2];
	    my $true_row_number=$temp_start_pointer+$h-1;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$read_cluster_data[$true_row_number][$j]=$temp_read_cluster_data[$row_index][$j];
	    }
	}
	
	#re-start the pointer;
	$moving_pointer++;
	$temp_start_pointer=$moving_pointer;
	$temp_end_pointer=$moving_pointer;
	if ($temp_start_pointer<=$cluster_row_number) {
	    $gene_ID=$read_cluster_data[$temp_start_pointer][$reduced_cross_column_number+1];
	}
    }
    
}

sub quicksort_read_cluster_data_by_anchor1_geneExon_after_anchor0_withoutClusterID {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1, $temp_2);
    
    if ($right>=$left) {
	$pivot=$sort_read_cluster_index[$right][1];
	# print "pivot= $pivot\n";
	$pointer_left=$left-1;
	$pointer_right=$right;
	do
	{
	    do
	    {
		$pointer_left=$pointer_left+1;
	    } until ($sort_read_cluster_index[$pointer_left][1] ge $pivot);
	    
	    do
	    { 
		$pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_read_cluster_index[$pointer_right][1] le $pivot));
	    
	    $temp_1=$sort_read_cluster_index[$pointer_left][1];
	    $sort_read_cluster_index[$pointer_left][1]=$sort_read_cluster_index[$pointer_right][1];
	    $sort_read_cluster_index[$pointer_right][1]=$temp_1;
	    
	    $temp_2=$sort_read_cluster_index[$pointer_left][2];
	    $sort_read_cluster_index[$pointer_left][2]=$sort_read_cluster_index[$pointer_right][2];
	    $sort_read_cluster_index[$pointer_right][2]=$temp_2;
	    
	} until ($pointer_right<=$pointer_left);
	
	$sort_read_cluster_index[$pointer_right][1]=$sort_read_cluster_index[$pointer_left][1];
	$sort_read_cluster_index[$pointer_left][1]=$sort_read_cluster_index[$right][1];
	$sort_read_cluster_index[$right][1]=$temp_1;
	
	$sort_read_cluster_index[$pointer_right][2]=$sort_read_cluster_index[$pointer_left][2];
	$sort_read_cluster_index[$pointer_left][2]=$sort_read_cluster_index[$right][2];
	$sort_read_cluster_index[$right][2]=$temp_2;
		
	quicksort_read_cluster_data_by_anchor1_geneExon_after_anchor0_withoutClusterID($left, $pointer_left-1);
	quicksort_read_cluster_data_by_anchor1_geneExon_after_anchor0_withoutClusterID($pointer_left+1, $right);
    }
}

sub define_cluster_with_reassign_FirstReadlocation {
    #my (@unit_data, $unit_row_number)=@_;

    #deal with first mate
    my (@temp_gene)=split("\/\/", $unit_data[1][$reduced_cross_column_number+1]);
    my $temp_index=$cDNA_order{$temp_gene[0]};
    my $exonCount=$hgGenome_refGene[$temp_index][1];
    #start to deal exon location
    my $exonlocation=$hgGenome_refGene[$temp_index][2];
    $exonlocation=~s/\|\|//;
    my (@temp_exonlocation)=split(",", $exonlocation);
    
    #put the exon location into an array
    my @true_exonlocation=();
    my @reverse_exonlocation=();
    for (my $i=1; $i<=$exonCount; $i++) {
	$true_exonlocation[$i][1]=$temp_exonlocation[$i-1];
	$true_exonlocation[$i][2]=$temp_exonlocation[$i-1+$exonCount];
	
	$reverse_exonlocation[$i][1]=$temp_exonlocation[$i-1];
	$reverse_exonlocation[$i][2]=$temp_exonlocation[$i-1+$exonCount];
    }

    #reverse (anti-sense) gene
    if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	for (my $i=$exonCount; $i>=1; $i--) {
	    $true_exonlocation[$i][1]=$reverse_exonlocation[$exonCount-$i+1][1];
	    $true_exonlocation[$i][2]=$reverse_exonlocation[$exonCount-$i+1][2];
	}

    }

    #start to deal with the read location
    my (@temp_exonID)=split("\/\/", $unit_data[1][$reduced_cross_column_number+1]);
    my $firstMate_exonID=$temp_exonID[1];
    $firstMate_exonID=~s/exon|intron//;
    for (my $i=2; $i<=$unit_row_number; $i++) {
	#forward read
	if ($unit_data[$i][2]=~/\+/) {
	    if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+1]);
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#bigger exon number
		if ($firstMate_exonID<$temp_firstMate_exonID) {
		    $firstMate_exonID=$temp_firstMate_exonID;
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	    elsif ($hgGenome_refGene[$temp_index][4]=~/\-/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+1]);
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#smaller exon number
		if ($firstMate_exonID>$temp_firstMate_exonID) {
		    $firstMate_exonID=$temp_firstMate_exonID;
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	}
	#reverse read
	elsif ($unit_data[$i][2]=~/\-/) {
	    if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+1]);
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#bigger exon number
		if ($firstMate_exonID<$temp_firstMate_exonID) {
		    $firstMate_exonID=$temp_firstMate_exonID;
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	    elsif ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+1]);
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#smaller exon number
		if ($firstMate_exonID>$temp_firstMate_exonID) {
		    $firstMate_exonID=$temp_firstMate_exonID;
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	    
	} #elsif ($unit_data[$i][2]==~/\-/) {    
    }

    #re-calculate the read location after removing the intron
    for (my $i=1; $i<=$unit_row_number; $i++) {
	#forward read
	if ($unit_data[$i][2]=~/\+/) {
	    if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+1]);
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#bigger exon number
		if ($firstMate_exonID!=$temp_firstMate_exonID) {
		    #intron size between rwo reads
		    my $intron_size=0;
		    for ($h=$temp_firstMate_exonID; $h<$firstMate_exonID; $h++) {
			$intron_size=$intron_size+$true_exonlocation[$h+1][1]-$true_exonlocation[$h][2];
		    }
		    #my $temp_loc=$unit_data[$i][1]+$intron_size;
		    #$unit_data[$i][1]=$unit_data[$i][1]."\/\/".$temp_loc;
                    $unit_data[$i][1]=$unit_data[$i][1]+$intron_size;
		    if ($temp_exonID[1]=~/exon/) {
			$unit_data[$i][$reduced_cross_column_number+1]=$temp_exonID[0]."\/\/exon".
			    $firstMate_exonID;
		    }
		    elsif ($temp_exonID[1]=~/intron/) {
			$unit_data[$i][$reduced_cross_column_number+1]=$temp_exonID[0]."\/\/intron".
			    $firstMate_exonID;
		    }
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	    elsif ($hgGenome_refGene[$temp_index][4]=~/\-/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+1]);
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#smaller exon number
		if ($firstMate_exonID!=$temp_firstMate_exonID) {
		    my $intron_size=0;
		    for ($h=$temp_firstMate_exonID; $h>$firstMate_exonID; $h--) {
			$intron_size=$intron_size+$true_exonlocation[$h-1][1]-$true_exonlocation[$h][2];
		    }
		    #my $temp_loc=$unit_data[$i][1]+$intron_size;
		    #$unit_data[$i][1]=$unit_data[$i][1]."\/\/".$temp_loc;
		    $unit_data[$i][1]=$unit_data[$i][1]+$intron_size;
		    if ($temp_exonID[1]=~/exon/) {
			$unit_data[$i][$reduced_cross_column_number+1]=$temp_exonID[0]."\/\/exon".
			    $firstMate_exonID;
		    }
		    elsif ($temp_exonID[1]=~/intron/) {
			$unit_data[$i][$reduced_cross_column_number+1]=$temp_exonID[0]."\/\/intron".
			    $firstMate_exonID;
		    }
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	}
	#reverse read
	elsif ($unit_data[$i][2]=~/\-/) {
	    if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+1]);
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#bigger exon number
		if ($firstMate_exonID!=$temp_firstMate_exonID) {
		    #intron size between rwo reads
		    my $intron_size=0;
		    for ($h=$temp_firstMate_exonID; $h<$firstMate_exonID; $h++) {
			$intron_size=$intron_size+$true_exonlocation[$h][1]-$true_exonlocation[$h+1][2];
		    }
		    #my $temp_loc=$unit_data[$i][1]-$intron_size;
		    #$unit_data[$i][1]=$unit_data[$i][1]."\/\/".$temp_loc;
		    $unit_data[$i][1]=$unit_data[$i][1]-$intron_size; 
		    if ($temp_exonID[1]=~/exon/) {
			$unit_data[$i][$reduced_cross_column_number+1]=$temp_exonID[0]."\/\/exon".
			    $firstMate_exonID;
		    }
		    elsif ($temp_exonID[1]=~/intron/) {
			$unit_data[$i][$reduced_cross_column_number+1]=$temp_exonID[0]."\/\/intron".
			    $firstMate_exonID;
		    }
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	    elsif ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+1]);
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#smaller exon number
		if ($firstMate_exonID!=$temp_firstMate_exonID) {
		    my $intron_size=0;
		    for ($h=$temp_firstMate_exonID; $h>$firstMate_exonID; $h--) {
			$intron_size=$intron_size+$true_exonlocation[$h][1]-$true_exonlocation[$h-1][2];
		    }
		    #my $temp_loc=$unit_data[$i][1]-$intron_size;
		    #$unit_data[$i][1]=$unit_data[$i][1]."\/\/".$temp_loc;
		    $unit_data[$i][1]=$unit_data[$i][1]-$intron_size;
		    if ($temp_exonID[1]=~/exon/) {
			$unit_data[$i][$reduced_cross_column_number+1]=$temp_exonID[0]."\/\/exon".
			    $firstMate_exonID;
		    }
		    elsif ($temp_exonID[1]=~/intron/) {
			$unit_data[$i][$reduced_cross_column_number+1]=$temp_exonID[0]."\/\/intron".
			    $firstMate_exonID;
		    } 
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	    
	} #elsif ($unit_data[$i][2]==~/\-/) {    
    }

    #return($firstMate_exonID);
   
}

sub define_cluster_with_reassign_SecondReadlocation {
    #my (@unit_data, $unit_row_number)=@_;   
    #deal with second mate
    my (@temp_gene)=split("\/\/", $unit_data[1][$reduced_cross_column_number+2]);
    my $temp_index=$cDNA_order{$temp_gene[0]};
    my $exonCount=$hgGenome_refGene[$temp_index][1];
    #start to deal exon location
    my $exonlocation=$hgGenome_refGene[$temp_index][2];
    $exonlocation=~s/\|\|//;
    my (@temp_exonlocation)=split(",", $exonlocation);
    
    #put the exon location into an array
    my @true_exonlocation=();
    my @reverse_exonlocation=();
    for (my $i=1; $i<=$exonCount; $i++) {
	$true_exonlocation[$i][1]=$temp_exonlocation[$i-1];
	$true_exonlocation[$i][2]=$temp_exonlocation[$i-1+$exonCount];
	
	$reverse_exonlocation[$i][1]=$temp_exonlocation[$i-1];
	$reverse_exonlocation[$i][2]=$temp_exonlocation[$i-1+$exonCount];
    }

    #reverse (anti-sense) gene
    if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	for (my $i=$exonCount; $i>=1; $i--) {
	    $true_exonlocation[$i][1]=$reverse_exonlocation[$exonCount-$i+1][1];
	    $true_exonlocation[$i][2]=$reverse_exonlocation[$exonCount-$i+1][2];
	}

    }

    
    #satrt to deal with the read location
    (@temp_exonID)=split("\/\/", $unit_data[1][$reduced_cross_column_number+2]);
    my $secondMate_exonID=$temp_exonID[1];
    $secondMate_exonID=~s/exon|intron//;
    for (my $i=2; $i<=$unit_row_number; $i++) {
	#forward read
	if ($unit_data[$i][7]=~/\+/) {
	    if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+2]);
		my $temp_secondMate_exonID=$temp_exonID[1];
		$temp_secondMate_exonID=~s/exon|intron//;
		#bigger exon number
		if ($secondMate_exonID<$temp_secondMate_exonID) {
		    $secondMate_exonID=$temp_secondMate_exonID;
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	    elsif ($hgGenome_refGene[$temp_index][4]=~/\-/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+2]);
		my $temp_secondMate_exonID=$temp_exonID[1];
		$temp_secondMate_exonID=~s/exon|intron//;
		#smaller exon number
		if ($secondMate_exonID>$temp_secondMate_exonID) {
		    $secondMate_exonID=$temp_secondMate_exonID;
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	}
	#reverse read
	elsif ($unit_data[$i][7]=~/\-/) {
	    if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+2]);
		my $temp_secondMate_exonID=$temp_exonID[1];
		$temp_secondMate_exonID=~s/exon|intron//;
		#bigger exon number
		if ($secondMate_exonID<$temp_secondMate_exonID) {
		    $secondMate_exonID=$temp_secondMate_exonID;
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	    elsif ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+2]);
		my $temp_secondMate_exonID=$temp_exonID[1];
		$temp_secondMate_exonID=~s/exon|intron//;
		#smaller exon number
		if ($secondMate_exonID>$temp_secondMate_exonID) {
		    $secondMate_exonID=$temp_secondMate_exonID;
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	    
	} #elsif ($unit_data[$i][7]==~/\-/) {    
    }
    
   #re-calculate the read location after removing the intron
    for (my $i=1; $i<=$unit_row_number; $i++) {
	#forward read
	if ($unit_data[$i][7]=~/\+/) {
	    if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+2]);
		my $temp_secondMate_exonID=$temp_exonID[1];
		$temp_secondMate_exonID=~s/exon|intron//;
		#bigger exon number
		if ($secondMate_exonID!=$temp_secondMate_exonID) {
		    #intron size between rwo reads
		    my $intron_size=0;
		    for ($h=$temp_secondMate_exonID; $h<$secondMate_exonID; $h++) {
			$intron_size=$intron_size+$true_exonlocation[$h+1][1]-$true_exonlocation[$h][2];
		    }
		    #my $temp_loc=$unit_data[$i][6]+$intron_size;
		    #$unit_data[$i][6]=$unit_data[$i][6]."\/\/".$temp_loc;
		    $unit_data[$i][6]=$unit_data[$i][6]+$intron_size; 
		    if ($temp_exonID[1]=~/exon/) {
			$unit_data[$i][$reduced_cross_column_number+2]=$temp_exonID[0]."\/\/exon".
			    $secondMate_exonID;
		    }
		    elsif ($temp_exonID[1]=~/intron/) {
			$unit_data[$i][$reduced_cross_column_number+2]=$temp_exonID[0]."\/\/intron".
			    $secondMate_exonID;
		    }
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	    elsif ($hgGenome_refGene[$temp_index][4]=~/\-/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+2]);
		my $temp_secondMate_exonID=$temp_exonID[1];
		$temp_secondMate_exonID=~s/exon|intron//;
		#smaller exon number
		if ($secondMate_exonID!=$temp_secondMate_exonID) {
		    my $intron_size=0;
		    for ($h=$temp_secondMate_exonID; $h>$secondMate_exonID; $h--) {
			$intron_size=$intron_size+$true_exonlocation[$h-1][1]-$true_exonlocation[$h][2];
		    }
		    #my $temp_loc=$unit_data[$i][6]+$intron_size;
		    #$unit_data[$i][6]=$unit_data[$i][6]."\/\/".$temp_loc;
		    $unit_data[$i][6]=$unit_data[$i][6]+$intron_size;
		    if ($temp_exonID[1]=~/exon/) {
			$unit_data[$i][$reduced_cross_column_number+2]=$temp_exonID[0]."\/\/exon".
			    $secondMate_exonID;
		    }
		    elsif ($temp_exonID[1]=~/intron/) {
			$unit_data[$i][$reduced_cross_column_number+2]=$temp_exonID[0]."\/\/intron".
			    $secondMate_exonID;
		    }
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	}
	#reverse read
	elsif ($unit_data[$i][7]=~/\-/) {
	    if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+2]);
		my $temp_secondMate_exonID=$temp_exonID[1];
		$temp_secondMate_exonID=~s/exon|intron//;
		#bigger exon number
		if ($secondMate_exonID!=$temp_secondMate_exonID) {
		    #intron size between rwo reads
		    my $intron_size=0;
		    for ($h=$temp_secondMate_exonID; $h<$secondMate_exonID; $h++) {
			$intron_size=$intron_size+$true_exonlocation[$h][1]-$true_exonlocation[$h+1][2];
		    }
		    #my $temp_loc=$unit_data[$i][6]-$intron_size;
		    #$unit_data[$i][6]=$unit_data[$i][6]."\/\/".$temp_loc;
		    $unit_data[$i][6]=$unit_data[$i][6]-$intron_size; 
		    if ($temp_exonID[1]=~/exon/) {
			$unit_data[$i][$reduced_cross_column_number+2]=$temp_exonID[0]."\/\/exon".
			    $secondMate_exonID;
		    }
		    elsif ($temp_exonID[1]=~/intron/) {
			$unit_data[$i][$reduced_cross_column_number+2]=$temp_exonID[0]."\/\/intron".
			    $secondMate_exonID;
		    }
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\+/) {
	    elsif ($hgGenome_refGene[$temp_index][4]=~/\+/) {
		my (@temp_exonID)=split("\/\/", $unit_data[$i][$reduced_cross_column_number+2]);
		my $temp_secondMate_exonID=$temp_exonID[1];
		$temp_secondMate_exonID=~s/exon|intron//;
		#smaller exon number
		if ($secondMate_exonID!=$temp_secondMate_exonID) {
		    my $intron_size=0;
		    for ($h=$temp_secondMate_exonID; $h>$secondMate_exonID; $h--) {
			$intron_size=$intron_size+$true_exonlocation[$h][1]-$true_exonlocation[$h-1][2];
		    }
		    #my $temp_loc=$unit_data[$i][6]-$intron_size;
		    #$unit_data[$i][6]=$unit_data[$i][6]."\/\/".$temp_loc;
		    $unit_data[$i][6]=$unit_data[$i][6]-$intron_size;
		    if ($temp_exonID[1]=~/exon/) {
			$unit_data[$i][$reduced_cross_column_number+2]=$temp_exonID[0]."\/\/exon".
			    $secondMate_exonID;
		    }
		    elsif ($temp_exonID[1]=~/intron/) {
			$unit_data[$i][$reduced_cross_column_number+2]=$temp_exonID[0]."\/\/intron".
			    $secondMate_exonID;
		    }
		}
	    } #if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	    
	} #elsif ($unit_data[$i][7]==~/\-/) {    
    }
   
    #return($secondMate_exonID);
}

sub unit_to_unit_clustering {
    my ($geneLoc, $readLoc)=@_;
    $MSK_stats_means=$ARGV[3]; #fragment length
    $MSK_stats_stdev=$ARGV[4]; #fragment stdev
    #for in-silico fragment length
    #the half read length=50
    #my $MSK_stats_3stdev_minorhalfread_distance=$MSK_stats_means+3*$MSK_stats_stdev-50;
    my $MSK_stats_3stdev_minorhalfread_distance=$MSK_stats_means+3*$MSK_stats_stdev;
   
    #deal with first mate
    my (@temp_gene)=split("\/\/", $unit_forCluster_data[1][$reduced_cross_column_number+$geneLoc]);
    my $temp_index=$cDNA_order{$temp_gene[0]};
    my $exonCount=$hgGenome_refGene[$temp_index][1];
    #start to deal exon location
    my $exonlocation=$hgGenome_refGene[$temp_index][2];
    $exonlocation=~s/\|\|//;
    my (@temp_exonlocation)=split(",", $exonlocation);
    
    #put the exon location into an array
    my @true_exonlocation=();
    my @reverse_exonlocation=();
    for (my $i=1; $i<=$exonCount; $i++) {
	$true_exonlocation[$i][1]=$temp_exonlocation[$i-1];
	$true_exonlocation[$i][2]=$temp_exonlocation[$i-1+$exonCount];
	
	$reverse_exonlocation[$i][1]=$temp_exonlocation[$i-1];
	$reverse_exonlocation[$i][2]=$temp_exonlocation[$i-1+$exonCount];
    }
    #reverse (anti-sense) gene
    if ($hgGenome_refGene[$temp_index][4]=~/\-/) {
	for (my $i=$exonCount; $i>=1; $i--) {
	    $true_exonlocation[$i][1]=$reverse_exonlocation[$exonCount-$i+1][1];
	    $true_exonlocation[$i][2]=$reverse_exonlocation[$exonCount-$i+1][2];
	}		
    }
  
    #start to deal with both gene sense and read orientation
    if (($hgGenome_refGene[$temp_index][4]=~/\+/) && ($unit_forCluster_data[1][$readLoc+1]=~/\-/)) {
	my $dis_read_clusterID=0;
	my $temp_start_pointer=1;
	my $temp_end_pointer=1;
	my (@exonID)=split("\/\/", $unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+$geneLoc]);    
	my $moving_pointer=0;
	while ($moving_pointer<$unit_forCluster_row_number) {
	    my @grouping_size=();
	    while ($moving_pointer<$unit_forCluster_row_number) {
		#calculate the distance between two exons
		my (@temp_exonID)=split("\/\/", $unit_forCluster_data[$moving_pointer+1][$reduced_cross_column_number+$geneLoc]);
		my $firstMate_exonID=$exonID[1];
		$firstMate_exonID=~s/exon|intron//;
		
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#intron size		
		my $between_exon_size=0;
		for ($h=$temp_firstMate_exonID; $h>$firstMate_exonID; $h--) {
		    $between_exon_size=$between_exon_size+$true_exonlocation[$h][1]-$true_exonlocation[$h-1][2];
		}
		
		#grouping
		$grouping_size[$moving_pointer+1]=
		    ($unit_forCluster_data[$moving_pointer+1][$readLoc]-$between_exon_size)-$unit_forCluster_data[$temp_start_pointer][$readLoc];
		#debugging
		#if (($hgGenome_refGene[$temp_index][7]=~/RARA/) && ($exonID[1]=~/exon/) && ($temp_exonID[1]=~/exon/)) {
		 #   print "Gene=$hgGenome_refGene[$temp_index][7]\t exon1=$exonID[1]\t exon2=$temp_exonID[1]\t groupingsize=$grouping_size[$moving_pointer+1]\t";
		  #  print "Intron=$between_exon_size\t loc1=$unit_forCluster_data[$temp_start_pointer][$readLoc]\t ";
		   # print "loc2=$unit_forCluster_data[$moving_pointer+1][$readLoc]\n";
		#}

		#if (($firstMate_exonID eq $temp_firstMate_exonID) || 
		if (($exonID[1]=~/exon/) && ($temp_exonID[1]=~/exon/)) {
		    if ($grouping_size[$moving_pointer+1]<=$MSK_stats_3stdev_minorhalfread_distance) {
			$moving_pointer++;
		    }
		    else {
			last;
		    }
		}
		elsif (($exonID[1]=~/intron/) && ($temp_exonID[1]=~/intron/) && ($exonID[1] eq $temp_exonID[1])) {
		    if ($grouping_size[$moving_pointer+1]<=$MSK_stats_3stdev_minorhalfread_distance) {
			$moving_pointer++;
		    }
		    else {
			last;
		    }
		}
		else {
		    last;
		}
		
	    }
	    
	    #start a new cluster ID
	    $dis_read_clusterID++;
	    $temp_end_pointer=$moving_pointer;
	    #remove outlier and re-set the moving_pointer at at least 10 pairs
	    if (abs($temp_end_pointer-$temp_start_pointer)>=9) {
		my $beforeCluster_number=0;
		@beforeCluster_data=();
		$grouping_size[$temp_start_pointer]=0;
		for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
		    $beforeCluster_number++;
		    $beforeCluster_data[$beforeCluster_number][1]=$grouping_size[$h];
		}
		outlier_detection_for_beforeClusterdata($beforeCluster_number);
		#re-set the pointer
		my $outlier_number_atBoundary=0;
		my $outlier_number_atBottom=0;
		for (my $h=1; $h<=$beforeCluster_number; $h++) {
		    if ($beforeCluster_data[$h][2]=~/Yes/) {
			$outlier_number_atBoundary++;
		    }
		    elsif ($beforeCluster_data[$h][2]=~/No/) {
			last;
		    }
		}
		for (my $h=$beforeCluster_number; $h>=1; $h--) {
		    if ($beforeCluster_data[$h][2]=~/Yes/) {
			$outlier_number_atBottom++;
		    }
		    elsif ($beforeCluster_data[$h][2]=~/No/) {
			last;
		    }
		}
		#Re-set
		#print "outlier_atBoundary=$outlier_number_atBoundary\toutlier_atBottom=$outlier_number_atBottom\t"; 
		#print "totalNumber=$beforeCluster_number\t gene=$hgGenome_refGene[$temp_index][7]\t\+\-\n";
		if ($outlier_number_atBoundary>=1) {
		    $moving_pointer=$temp_start_pointer+$outlier_number_atBoundary;
		}
		elsif ($outlier_number_atBottom>=1) {
		    $moving_pointer=$moving_pointer-$outlier_number_atBottom;
		    if ($moving_pointer<$temp_start_pointer) { $moving_pointer=$temp_start_pointer; }
		}
	    }
	    
	    #re-set the final pointer
	    $temp_end_pointer=$moving_pointer;	
	    #####unit_data
	    @unit_data=();
	    $unit_row_number=0;
	    for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
		$unit_row_number++;
		for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		    $unit_data[$unit_row_number][$j]=$unit_forCluster_data[$h][$j];
		}
	    }
	    
	    if ($geneLoc==1) {
		define_cluster_with_reassign_FirstReadlocation();
	    }
	    elsif ($geneLoc==2) {
		define_cluster_with_reassign_SecondReadlocation();
	    }
	    #put data back into unit_forCluster after clustering
	    $unit_row_number=0;
	    for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
		$unit_row_number++;
		for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		    $unit_forCluster_data[$h][$j]=$unit_data[$unit_row_number][$j];		    
		}
		#just for anchor0
		for (my $j=($reduced_cross_column_number+$geneLoc); $j<=($reduced_cross_column_number+$geneLoc); $j++) {
		    $unit_forCluster_data[$h][$j]=$unit_data[$unit_row_number][$j]."\-".$dis_read_clusterID;
		}
		
	    }
	    
	    #re-start the pointer;
	    $moving_pointer++;
	    $temp_start_pointer=$moving_pointer;
	    $temp_end_pointer=$moving_pointer;
	    if ($temp_start_pointer<=$unit_forCluster_row_number) {
		(@exonID)=split("\/\/", $unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+$geneLoc]);
	    }
	}
    } #if (($hgGenome_refGene[$temp_index][4]=~/\+/) && if ($unit_forCluster_data[1][2]=~/\-/))

    elsif (($hgGenome_refGene[$temp_index][4]=~/\+/) && ($unit_forCluster_data[1][$readLoc+1]=~/\+/)) {
	my $dis_read_clusterID=0;
	my $temp_start_pointer=$unit_forCluster_row_number;
	my $temp_end_pointer=$unit_forCluster_row_number;
	my (@exonID)=split("\/\/", $unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+$geneLoc]);    
	my $moving_pointer=$unit_forCluster_row_number+1;
	while ($moving_pointer>1) {
	    my @grouping_size=();
	    while ($moving_pointer>1) {
		#calculate the distance between two exons
		my (@temp_exonID)=split("\/\/", $unit_forCluster_data[$moving_pointer-1][$reduced_cross_column_number+$geneLoc]);
		my $firstMate_exonID=$exonID[1];
		$firstMate_exonID=~s/exon|intron//;
		
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		
		#intron size		
		my $between_exon_size=0;
		for ($h=$temp_firstMate_exonID; $h<$firstMate_exonID; $h++) {
		    $between_exon_size=$between_exon_size+$true_exonlocation[$h+1][1]-$true_exonlocation[$h][2];
		}
		
		#grouping
		$grouping_size[$moving_pointer-1]=
		    $unit_forCluster_data[$temp_start_pointer][$readLoc]-($unit_forCluster_data[$moving_pointer-1][$readLoc]+$between_exon_size);
		#if (($hgGenome_refGene[$temp_index][7]=~/RARA/) && ($exonID[1]=~/exon/) && ($temp_exonID[1]=~/exon/)) {
		 #   print "Gene=$hgGenome_refGene[$temp_index][7]\t exon1=$exonID[1]\t exon2=$temp_exonID[1]\t groupingsize=$grouping_size[$moving_pointer-1]\t";
		 #   print "Intron=$between_exon_size\t loc1=$unit_forCluster_data[$temp_start_pointer][$readLoc]\t ";
		 #   print "loc2=$unit_forCluster_data[$moving_pointer-1][$readLoc]\n";
		#}

		if (($exonID[1]=~/exon/) && ($temp_exonID[1]=~/exon/)) {
		    if ($grouping_size[$moving_pointer-1]<=$MSK_stats_3stdev_minorhalfread_distance) {
			$moving_pointer--;
		    }
		    else {
			last;
		    }
		}
		elsif (($exonID[1]=~/intron/) && ($temp_exonID[1]=~/intron/) && ($exonID[1] eq $temp_exonID[1])) {
		    if ($grouping_size[$moving_pointer-1]<=$MSK_stats_3stdev_minorhalfread_distance) {
			$moving_pointer--;
		    }
		    else {
			last;
		    }
		}
		else {
		    last;
		}

	    }
	    
	    #start a new cluster ID
	    $dis_read_clusterID++;	    
	    #start to assign genename the ordinal number
	    $temp_end_pointer=$moving_pointer;
	    #remove outlier and re-set the moving_pointer at at least 10 pairs
	    if (abs($temp_end_pointer-$temp_start_pointer)>=9) {
		my $beforeCluster_number=$temp_start_pointer-$temp_end_pointer+1;
		@beforeCluster_data=();
		$grouping_size[$temp_start_pointer]=0;
		for (my $h=$temp_start_pointer; $h>=$temp_end_pointer; $h--) {
		    $beforeCluster_data[$beforeCluster_number][1]=$grouping_size[$h];
		    $beforeCluster_number--;
		}
		
		$beforeCluster_number=$temp_start_pointer-$temp_end_pointer+1;
		outlier_detection_for_beforeClusterdata($beforeCluster_number);
		#re-set the pointer
		my $outlier_number_atBoundary=0;
		my $outlier_number_atBottom=0;
		for (my $h=$beforeCluster_number; $h>=1; $h--) {
		    if ($beforeCluster_data[$h][2]=~/Yes/) {
			$outlier_number_atBoundary++;
		    }
		    elsif ($beforeCluster_data[$h][2]=~/No/) {
			last;
		    }
		}
		for (my $h=1; $h<=$beforeCluster_number; $h++) {
		    if ($beforeCluster_data[$h][2]=~/Yes/) {
			$outlier_number_atBottom++;
		    }
		    elsif ($beforeCluster_data[$h][2]=~/No/) {
			last;
		    }
		}
		#Re-set
		#print "outlier_atBoundary=$outlier_number_atBoundary\t outlier_atBottom=$outlier_number_atBottom\t";
		#print "totalNumber=$beforeCluster_number\tgene=$hgGenome_refGene[$temp_index][7]\t\+\+\n";
		if ($outlier_number_atBoundary>=1) {
		    #print "outlier=$outlier_number\n";
		    $moving_pointer=$temp_start_pointer-$outlier_number_atBoundary;
		}
		elsif ($outlier_number_atBottom>=1) {
		    #print "outlier=$outlier_number\n";
		    $moving_pointer=$moving_pointer+$outlier_number_atBottom;
		    if ($moving_pointer>$temp_start_pointer) { $moving_pointer=$temp_start_pointer; }
		}
	    }

	    #re-set the final pointer
	    $temp_end_pointer=$moving_pointer;
	
	    #####unit_data
	    @unit_data=();
	    $unit_row_number=0;
	    for (my $h=$temp_start_pointer; $h>=$temp_end_pointer; $h--) {
		$unit_row_number++;
		for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		    $unit_data[$unit_row_number][$j]=$unit_forCluster_data[$h][$j];
		}
	    }
	    
	    if ($geneLoc==1) {
		define_cluster_with_reassign_FirstReadlocation();
	    }
	    elsif ($geneLoc==2) {
		define_cluster_with_reassign_SecondReadlocation();
	    }
	    #put data back into unit_forCluster after clustering
	    $unit_row_number=0;
	    for (my $h=$temp_start_pointer; $h>=$temp_end_pointer; $h--) {
		$unit_row_number++;
		for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		    $unit_forCluster_data[$h][$j]=$unit_data[$unit_row_number][$j];
		}
		#just for anchor0
		for (my $j=($reduced_cross_column_number+$geneLoc); $j<=($reduced_cross_column_number+$geneLoc); $j++) {
		    $unit_forCluster_data[$h][$j]=$unit_data[$unit_row_number][$j]."\-".$dis_read_clusterID;  
		}
	    }
	    
	    #re-start the pointer;
	    $moving_pointer--;
	    $temp_start_pointer=$moving_pointer;
	    $temp_end_pointer=$moving_pointer;
	    if ($temp_start_pointer>=1) {
		(@exonID)=split("\/\/", $unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+$geneLoc])
	    }
	}
    } #if (($hgGenome_refGene[$temp_index][4]=~/\+/) && if ($unit_forCluster_data[1][2]=~/\+/))

    elsif (($hgGenome_refGene[$temp_index][4]=~/\-/) && ($unit_forCluster_data[1][$readLoc+1]=~/\-/)) {
	my $dis_read_clusterID=0;
	my $temp_start_pointer=1;
	my $temp_end_pointer=1;
	my (@exonID)=split("\/\/", $unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+$geneLoc]);    
	my $moving_pointer=0;
	while ($moving_pointer<$unit_forCluster_row_number) {
	    my @grouping_size=();
	    while ($moving_pointer<$unit_forCluster_row_number) {
		#calculate the distance between two exons
		my (@temp_exonID)=split("\/\/", $unit_forCluster_data[$moving_pointer+1][$reduced_cross_column_number+$geneLoc]);
		my $firstMate_exonID=$exonID[1];
		$firstMate_exonID=~s/exon|intron//;
		
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		
		#intron size
		my $between_exon_size=0;
		for ($h=$temp_firstMate_exonID; $h<$firstMate_exonID; $h++) {
		    $between_exon_size=$between_exon_size+$true_exonlocation[$h][1]-$true_exonlocation[$h+1][2];
		}
		
		#grouping
		$grouping_size[$moving_pointer+1]=
		    ($unit_forCluster_data[$moving_pointer+1][$readLoc]-$between_exon_size)-$unit_forCluster_data[$temp_start_pointer][$readLoc];

		#if (($hgGenome_refGene[$temp_index][7]=~/RUNX1/) && ($exonID[1]=~/exon/) && ($temp_exonID[1]=~/exon/)) {
		#    print "Gene=$hgGenome_refGene[$temp_index][7]\t exon1=$exonID[1]\t exon2=$temp_exonID[1]\t groupingsize=$grouping_size[$moving_pointer+1]\t";
		#    print "Intron=$between_exon_size\t loc1=$unit_forCluster_data[$temp_start_pointer][1]\t ";
		#    print "loc2=$unit_forCluster_data[$moving_pointer+1][1]\n";
		#}

		#if (($firstMate_exonID eq $temp_firstMate_exonID) ||
		if (($exonID[1]=~/exon/) && ($temp_exonID[1]=~/exon/)) {
		    if ($grouping_size[$moving_pointer+1]<=$MSK_stats_3stdev_minorhalfread_distance) {
			$moving_pointer++;
		    }
		    else {
			last;
		    }
		}
		elsif (($exonID[1]=~/intron/) && ($temp_exonID[1]=~/intron/) && ($exonID[1] eq $temp_exonID[1])) {
		    if ($grouping_size[$moving_pointer+1]<=$MSK_stats_3stdev_minorhalfread_distance) {
			$moving_pointer++;
		    }
		    else {
			last;
		    }
		}
		else {
		    last;
		}

	    }
	    
	    #start a new cluster ID
	    $dis_read_clusterID++;
	    #start to assign genename the ordinal number
	    $temp_end_pointer=$moving_pointer;	
	    #remove outlier and re-set the moving_pointer with at least 10 pairs
	    if (abs($temp_end_pointer-$temp_start_pointer)>=9) {
		my $beforeCluster_number=0;
		@beforeCluster_data=();
		$grouping_size[$temp_start_pointer]=0;
		for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
		    $beforeCluster_number++;
		    $beforeCluster_data[$beforeCluster_number][1]=$grouping_size[$h];
		}
		outlier_detection_for_beforeClusterdata($beforeCluster_number);
		#re-set the pointer
		my $outlier_number_atBoundary=0;
		my $outlier_number_atBottom=0;
		for (my $h=1; $h<=$beforeCluster_number; $h++) {
		    if ($beforeCluster_data[$h][2]=~/Yes/) {
			$outlier_number_atBoundary++;
		    }
		    elsif ($beforeCluster_data[$h][2]=~/No/) {
			last;
		    }
		}
		for (my $h=$beforeCluster_number; $h>=1; $h--) {
		    if ($beforeCluster_data[$h][2]=~/Yes/) {
			$outlier_number_atBottom++;
		    }
		    elsif ($beforeCluster_data[$h][2]=~/No/) {
			last;
		    }
		}
		#Re-set
		#print "outlier_atBoundary=$outlier_number_atBoundary\t outlier_atBottom=$outlier_number_atBottom\t";
		#print "totalNumber=$beforeCluster_number\tgene=$hgGenome_refGene[$temp_index][7]\t\-\-\n";
		if ($outlier_number_atBoundary>=1) {
		    $moving_pointer=$temp_start_pointer+$outlier_number_atBoundary;
		}
		elsif ($outlier_number_atBottom>=1) {
		    $moving_pointer=$moving_pointer-$outlier_number_atBottom;
		    if ($moving_pointer<$temp_start_pointer) { $moving_pointer=$temp_start_pointer; }
		}
	    }

	    #re-set the final pointer
	    $temp_end_pointer=$moving_pointer;

	    #####unit_data
	    @unit_data=();
	    $unit_row_number=0;
	    for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
		$unit_row_number++;
		for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		    $unit_data[$unit_row_number][$j]=$unit_forCluster_data[$h][$j];
		}
	    }
	    
	    if ($geneLoc==1) {
		define_cluster_with_reassign_FirstReadlocation();
	    }
	    elsif ($geneLoc==2) {
		define_cluster_with_reassign_SecondReadlocation();
	    }
	    #put data back into unit_forCluster after clustering
	    $unit_row_number=0;
	    for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
		$unit_row_number++;
		for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		    $unit_forCluster_data[$h][$j]=$unit_data[$unit_row_number][$j];
		}
		#just for anchor0
		for (my $j=($reduced_cross_column_number+$geneLoc); $j<=($reduced_cross_column_number+$geneLoc); $j++) {
		    $unit_forCluster_data[$h][$j]=$unit_data[$unit_row_number][$j]."\-".$dis_read_clusterID;
		}
	    }
	    
	    #re-start the pointer;
	    $moving_pointer++;
	    $temp_start_pointer=$moving_pointer;
	    $temp_end_pointer=$moving_pointer;
	    if ($temp_start_pointer<=$unit_forCluster_row_number) {
		(@exonID)=split("\/\/", $unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+$geneLoc]);
	    }
	}
    } #if (($hgGenome_refGene[$temp_index][4]=~/\-/) && if ($unit_forCluster_data[1][2]=~/\-))

    elsif (($hgGenome_refGene[$temp_index][4]=~/\-/) && ($unit_forCluster_data[1][$readLoc+1]=~/\+/)) {
	my $dis_read_clusterID=0;
	my $temp_start_pointer=$unit_forCluster_row_number;
	my $temp_end_pointer=$unit_forCluster_row_number;
	my (@exonID)=split("\/\/", $unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+$geneLoc]);    
	my $moving_pointer=$unit_forCluster_row_number+1;
	while ($moving_pointer>1) {
	    my @grouping_size=();
	    while ($moving_pointer>1) {
		#calculate the distance between two exons
		my (@temp_exonID)=split("\/\/", $unit_forCluster_data[$moving_pointer-1][$reduced_cross_column_number+$geneLoc]);
		my $firstMate_exonID=$exonID[1];
		$firstMate_exonID=~s/exon|intron//;
		
		my $temp_firstMate_exonID=$temp_exonID[1];
		$temp_firstMate_exonID=~s/exon|intron//;
		#intron size
		my $between_exon_size=0;
		for ($h=$temp_firstMate_exonID; $h>$firstMate_exonID; $h--) {
		    $between_exon_size=$between_exon_size+$true_exonlocation[$h-1][1]-$true_exonlocation[$h][2];
		}
		
		#grouping
		$grouping_size[$moving_pointer-1]=
		    $unit_forCluster_data[$temp_start_pointer][$readLoc]-($unit_forCluster_data[$moving_pointer-1][$readLoc]+$between_exon_size);
		#if (($firstMate_exonID eq $temp_firstMate_exonID) || 

		if (($exonID[1]=~/exon/) && ($temp_exonID[1]=~/exon/)) {
		    if ($grouping_size[$moving_pointer-1]<=$MSK_stats_3stdev_minorhalfread_distance) {
			$moving_pointer--;
		    }
		    else {
			last;
		    }
		}
		elsif (($exonID[1]=~/intron/) && ($temp_exonID[1]=~/intron/) && ($exonID[1] eq $temp_exonID[1])) {
		    if ($grouping_size[$moving_pointer-1]<=$MSK_stats_3stdev_minorhalfread_distance) {
			$moving_pointer--;
		    }
		    else {
			last;
		    }
		}
		else {
		    last;
		}
	    }

	    #start a new cluster ID
	    $dis_read_clusterID++;	    
	    #start to assign genename the ordinal number
	    $temp_end_pointer=$moving_pointer;	
	    #remove outlier and re-set the moving_pointer at at least 10 pairs
	    if (abs($temp_end_pointer-$temp_start_pointer)>=9) {
		my $beforeCluster_number=$temp_start_pointer-$temp_end_pointer+1;
		@beforeCluster_data=();
		$grouping_size[$temp_start_pointer]=0;
		for (my $h=$temp_start_pointer; $h>=$temp_end_pointer; $h--) {
		    $beforeCluster_data[$beforeCluster_number][1]=$grouping_size[$h];
		    $beforeCluster_number--;
		}
		
		$beforeCluster_number=$temp_start_pointer-$temp_end_pointer+1;
		outlier_detection_for_beforeClusterdata($beforeCluster_number);
		#re-set the pointer
		my $outlier_number_atBoundary=0;
		my $outlier_number_atBottom=0;
		for (my $h=$beforeCluster_number; $h>=1; $h--) {
		    if ($beforeCluster_data[$h][2]=~/Yes/) {
			$outlier_number_atBoundary++;
		    }
		    elsif ($beforeCluster_data[$h][2]=~/No/) {
			last;
		    }
		}
		for (my $h=1; $h<=$beforeCluster_number; $h++) {
		    if ($beforeCluster_data[$h][2]=~/Yes/) {
			$outlier_number_atBottom++;
		    }
		    elsif ($beforeCluster_data[$h][2]=~/No/) {
			last;
		    }
		}
		#Re-set
		#print "outlier_atBoundary=$outlier_number_atBoundary\t outlier_atBottom=$outlier_number_atBottom\t";
		#print "totalNumber=$beforeCluster_number\tgene=$hgGenome_refGene[$temp_index][7]\t\-\+\n";
		if ($outlier_number_atBoundary>=1) {
		    #print "outlier=$outlier_number\n";
		    $moving_pointer=$temp_start_pointer-$outlier_number_atBoundary;
		}
		elsif ($outlier_number_atBottom>=1) {
		    #print "outlier=$outlier_number\n";
		    $moving_pointer=$moving_pointer+$outlier_number_atBottom;
		    if ($moving_pointer>$temp_start_pointer) { $moving_pointer=$temp_start_pointer; }
		}
	    }

	    #re-set the final pointer
	    $temp_end_pointer=$moving_pointer;

	    #####unit_data
	    @unit_data=();
	    $unit_row_number=0;
	    for (my $h=$temp_start_pointer; $h>=$temp_end_pointer; $h--) {
		$unit_row_number++;
		for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		    $unit_data[$unit_row_number][$j]=$unit_forCluster_data[$h][$j];
		}
	    }
	    if ($geneLoc==1) {
		define_cluster_with_reassign_FirstReadlocation();
	    }
	    elsif ($geneLoc==2) {
		define_cluster_with_reassign_SecondReadlocation();
	    }
	    #put data back into unit_forCluster after clustering
	    $unit_row_number=0;
	    for (my $h=$temp_start_pointer; $h>=$temp_end_pointer; $h--) {
		$unit_row_number++;
		for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		    $unit_forCluster_data[$h][$j]=$unit_data[$unit_row_number][$j];
		}
		#just for anchor0
		for (my $j=($reduced_cross_column_number+$geneLoc); $j<=($reduced_cross_column_number+$geneLoc); $j++) {
		    $unit_forCluster_data[$h][$j]=$unit_data[$unit_row_number][$j]."\-".$dis_read_clusterID;   
		}
	    }
	    
	    #re-start the pointer;
	    $moving_pointer--;
	    $temp_start_pointer=$moving_pointer;
	    $temp_end_pointer=$moving_pointer;
	    if ($temp_start_pointer>=1) {
		(@exonID)=split("\/\/", $unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+$geneLoc]);
	    }
	}
    } #if (($hgGenome_refGene[$temp_index][4]=~/\-/) && if ($unit_forCluster_data[1][2]=~/\+))

}

sub sort_unitforClusterdata_by_anchor0_readLoc {
    
    #sort the data by anchor0 location
    @sort_exon_index=();
    for (my $i=1; $i<=$unit_forCluster_row_number; $i++) {
	$sort_exon_index[$i][1]=$unit_forCluster_data[$i][1];
	$sort_exon_index[$i][2]=$i;
    }

    quicksort_unitforClusterdata_by_readLoc(1, $unit_forCluster_row_number);

    #re-order the data by geneName
    my @temp_sort_exon_data=();
    for (my $i=1; $i<=$unit_forCluster_row_number; $i++) {
	my $row_index=$sort_exon_index[$i][2];
	for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
	    $temp_sort_exon_data[$i][$j]=$unit_forCluster_data[$row_index][$j];
	}
    }

    undef(@unit_forCluster_data);
    @unit_forCluster_data=();
    @unit_forCluster_data=@temp_sort_exon_data;  
    undef(@temp_sort_exon_data);

}

sub sort_unitforClusterdata_by_anchor1_readLoc {
    
    #sort the data by anchor0 location
    @sort_exon_index=();
    for (my $i=1; $i<=$unit_forCluster_row_number; $i++) {
	$sort_exon_index[$i][1]=$unit_forCluster_data[$i][6];
	$sort_exon_index[$i][2]=$i;
    }

    quicksort_unitforClusterdata_by_readLoc(1, $unit_forCluster_row_number);

    #re-order the data by geneName
    my @temp_sort_exon_data=();
    for (my $i=1; $i<=$unit_forCluster_row_number; $i++) {
	my $row_index=$sort_exon_index[$i][2];
	for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
	    $temp_sort_exon_data[$i][$j]=$unit_forCluster_data[$row_index][$j];
	}
    }

    undef(@unit_forCluster_data);
    @unit_forCluster_data=();
    @unit_forCluster_data=@temp_sort_exon_data;  
    undef(@temp_sort_exon_data);

}

sub resort_unitforClusterdata_by_anchor1_readLoc {
    
    my $temp_start_pointer=1;
    my $temp_end_pointer=1;
    my $gene_ID=$unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+1];
    
    my $moving_pointer=0;
    while ($moving_pointer<$unit_forCluster_row_number) {
	while (($moving_pointer<$unit_forCluster_row_number) && 
	       ($gene_ID eq $unit_forCluster_data[$moving_pointer+1][$reduced_cross_column_number+1])) {
	    $moving_pointer++;
	}
	
	#start to assign genename the ordinal number
	$temp_end_pointer=$moving_pointer;
	
	
	#re-sort the data by anchor1 genename
	my @temp_anchorGene_data=();
	my $temp_anchorGene_number=0;
	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	    $temp_anchorGene_number++;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$temp_anchorGene_data[$temp_anchorGene_number][$j]=$unit_forCluster_data[$h][$j];
	    }
	}

	@sort_exon_index=();
	for (my $h=1; $h<=$temp_anchorGene_number; $h++) {
	    $sort_exon_index[$h][1]=$temp_anchorGene_data[$h][6];
	    $sort_exon_index[$h][2]=$h;
	}
	quicksort_unitforClusterdata_by_readLoc(1, $temp_anchorGene_number);
	
	#put re-sorted data back into cross-anchordata
	for (my $h=1; $h<=$temp_anchorGene_number; $h++) {
	    my $row_index=$sort_exon_index[$h][2];
	    my $true_row_number=$temp_start_pointer+$h-1;
	    for (my $j=0; $j<=($reduced_cross_column_number+2); $j++) {
		$unit_forCluster_data[$true_row_number][$j]=$temp_anchorGene_data[$row_index][$j];
	    }
	}
	
	#re-start the pointer;
	$moving_pointer++;
	$temp_start_pointer=$moving_pointer;
	$temp_end_pointer=$moving_pointer;
	if ($moving_pointer<=$unit_forCluster_row_number) {
	    $gene_ID=$unit_forCluster_data[$temp_start_pointer][$reduced_cross_column_number+1];
	}	
    }
 
}

sub quicksort_unitforClusterdata_by_readLoc {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1, $temp_2);

    if ($right>=$left) {
	$pivot=$sort_exon_index[$right][1];
	# print "pivot= $pivot\n";
	$pointer_left=$left-1;
	$pointer_right=$right;
	do
	{
	    do
	    {
		$pointer_left=$pointer_left+1;
	    } until ($sort_exon_index[$pointer_left][1] >= $pivot);
	    
	    do
	    { 
		$pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_exon_index[$pointer_right][1] <= $pivot));
	    
	    $temp_1=$sort_exon_index[$pointer_left][1];
	    $sort_exon_index[$pointer_left][1]=$sort_exon_index[$pointer_right][1];
	    $sort_exon_index[$pointer_right][1]=$temp_1;
	    
	    $temp_2=$sort_exon_index[$pointer_left][2];
	    $sort_exon_index[$pointer_left][2]=$sort_exon_index[$pointer_right][2];
	    $sort_exon_index[$pointer_right][2]=$temp_2;
	    
	} until ($pointer_right<=$pointer_left);
	
	$sort_exon_index[$pointer_right][1]=$sort_exon_index[$pointer_left][1];
	$sort_exon_index[$pointer_left][1]=$sort_exon_index[$right][1];
	$sort_exon_index[$right][1]=$temp_1;
	
	$sort_exon_index[$pointer_right][2]=$sort_exon_index[$pointer_left][2];
	$sort_exon_index[$pointer_left][2]=$sort_exon_index[$right][2];
	$sort_exon_index[$right][2]=$temp_2;
	
	
	quicksort_unitforClusterdata_by_readLoc($left, $pointer_left-1);
	quicksort_unitforClusterdata_by_readLoc($pointer_left+1, $right);
    }
}

sub outlier_detection_for_beforeClusterdata {
    my ($beforeCluster_number)=@_;
    
    @max_ESD=();		
    %outlier_index=();
    
    for (my $beforeCluster_order=1; $beforeCluster_order<=$beforeCluster_number; $beforeCluster_order++) {
	$outlier_index{$beforeCluster_order}=0;
    }
    
    $ESD_sig_number=calculate_ESD_by_beforeCluster($beforeCluster_number);
    
    #assign the outlier
    for (my $k=1; $k<=$ESD_sig_number; $k++) {
	$outlier_index{$max_ESD[$k][2]}=1;
    }
    
    
    #generate the indicator for outliers
    for ($h=1; $h<=$beforeCluster_number; $h++) {
	$beforeCluster_data[$h][2]="No";
    }
    #assign the outliers
    if ($ESD_sig_number>=1) {
	for (my $k=1; $k<=$ESD_sig_number; $k++) {
	    if ($outlier_index{$max_ESD[$k][2]}==1) {
		$beforeCluster_data[$max_ESD[$k][2]][2]="Yes";
	    }
	}
	
    } #for (my $j=1; $j<=$affyBeforeCluster_column_number; $j++)
    
}

sub calculate_ESD_by_beforeCluster {

    my ($beforeCluster_number)=@_;
    $max_outlier_number=int($beforeCluster_number/10);
    $ESD_sig_level=0.05;
    
    #my ($beforeCluster_line)=@_;
    
    #initialize
    my %index_buffer_for_ESD=();
    for (my $temp_number=1; $temp_number<=$beforeCluster_number; $temp_number++) {
	$index_buffer_for_ESD{$temp_number}=0;
    }
    
    for (my $k=1; $k<=$max_outlier_number; $k++) {
	
	my @temp_beforeCluster_data=();
	for (my $temp_number=1; $temp_number<=$beforeCluster_number; $temp_number++) {
	    if ($index_buffer_for_ESD{$temp_number}==0) {
		#use SU's formula instead of LI's formula
		$temp_beforeCluster_data[$temp_number]=$beforeCluster_data[$temp_number][1];
	    }
	    elsif ($index_buffer_for_ESD{$temp_number}==1) {
		$temp_beforeCluster_data[$temp_number]="NA"; #excluded data point
	    }
	}
	
	
	#calculate the mean and Stdev  
	my $temp_mean=0;
	my $temp_square=0;
	
	for (my $temp_number=1; $temp_number<=$beforeCluster_number; $temp_number++) {
	    if (!($temp_beforeCluster_data[$temp_number]=~/NA/)) {
		$temp_mean=$temp_mean+$temp_beforeCluster_data[$temp_number];
		$temp_square=$temp_square+$temp_beforeCluster_data[$temp_number]*$temp_beforeCluster_data[$temp_number];
	    }
	}
	
	#calculate ESD
	
	$temp_mean=$temp_mean/($beforeCluster_number-$k+1);
	my $overall_mean=$temp_mean;
	$temp_mean=$temp_mean*$temp_mean;
	my $overall_stdev=1;
	
	my $temp_condition=$temp_square-($beforeCluster_number-$k+1)*$temp_mean;
	if ($temp_condition<=0) {
	    $overall_stdev=1;
	}
	elsif ($temp_condition>0) {
	    $overall_stdev=sqrt($temp_condition/($beforeCluster_number-$k));
	}
	
	@temp_ESD_value=();
	for (my $temp_number=1; $temp_number<=$beforeCluster_number; $temp_number++) {
	    if (!($temp_beforeCluster_data[$temp_number]=~/NA/)) {
		$temp_ESD_value[$temp_number][1]=abs($temp_beforeCluster_data[$temp_number]-$overall_mean)/$overall_stdev;
		$temp_ESD_value[$temp_number][2]=$temp_number;
	    }
	    elsif ($temp_beforeCluster_data[$temp_number]=~/NA/) {
		$temp_ESD_value[$temp_number][1]=0;
		$temp_ESD_value[$temp_number][2]=$temp_number;
	    }
	}
	
    
	quicksort_lookup_table_each_beforeCluster(1, $beforeCluster_number);
	
	$max_ESD[$k][1]=$temp_ESD_value[$beforeCluster_number][1];
	$max_ESD[$k][2]=$temp_ESD_value[$beforeCluster_number][2]; #save the index
	
	$index_buffer_for_ESD{$temp_ESD_value[$beforeCluster_number][2]}=1; #max ESD
	
	undef(@temp_beforeCluster_data);
	undef(@temp_ESD_value);

	#print "k=$k ESD=$max_ESD[$k][1]\n";

    } #for (my $k=1; $k<=$max_outlier_number; $k++)
    
    
    #calculate the p-value for each ESD
    #ESD_distribution();
    for (my $k=1; $k<=$max_outlier_number; $k++) {
	my $n=$beforeCluster_number-$k+1;

#	my $z=$max_ESD[$k][1];
#	if ((($n*($n-2)*$z*$z)/(($n-1)*($n-1)-$n*$z*$z))<=0) {
#	    $z=1;
#	    # print "the beforeCluster beforeCluster line= $beforeCluster_line  z=$z\n";
#	}
#	
#	my $T=sqrt(($n*($n-2)*$z*$z)/(($n-1)*($n-1)-$n*$z*$z));    
#	my $p=2*Statistics::Distributions::tprob($n-2, $T);
#	
#	$max_ESD[$k][3]=$p*$n; #save the pvalue
	$max_ESD[$k][3]=ESD_pvalue_calculation($n, $k, $ESD_sig_level);
	#print "pvalue= $max_ESD[$k][3]\n";
    }
    
    #decision rule for significant beforeClusters
    
    for (my $k=$max_outlier_number; $k>=1; $k--) {
	if ($max_ESD[$k][3]<=$ESD_sig_level) {
	    return($k);
	    # last;
	}
    }
    return(0);
}

sub ESD_pvalue_calculation {
    my ($n, $k, $alpha)=@_;
    if ($alpha>0.01) {
	if ($n==5) {
	    if ($max_ESD[$k][1]<1.72) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==6) {
	    if ($max_ESD[$k][1]<1.89) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==7) {
	    if ($max_ESD[$k][1]<2.02) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==8) {
	    if ($max_ESD[$k][1]<2.13) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==9) {
	    if ($max_ESD[$k][1]<2.21) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==10) {
	    if ($max_ESD[$k][1]<2.29) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==11) {
	    if ($max_ESD[$k][1]<2.36) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==12) {
	    if ($max_ESD[$k][1]<2.41) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==13) {
	    if ($max_ESD[$k][1]<2.46) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==14) {
	    if ($max_ESD[$k][1]<2.51) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==15) {
	    if ($max_ESD[$k][1]<2.55) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==16) {
	    if ($max_ESD[$k][1]<2.59) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==17) {
	    if ($max_ESD[$k][1]<2.62) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==18) {
	    if ($max_ESD[$k][1]<2.65) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==19) {
	    if ($max_ESD[$k][1]<2.68) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==20) {
	    if ($max_ESD[$k][1]<2.71) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==21) {
	    if ($max_ESD[$k][1]<2.73) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==22) {
	    if ($max_ESD[$k][1]<2.76) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==23) {
	    if ($max_ESD[$k][1]<2.78) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==24) {
	    if ($max_ESD[$k][1]<2.80) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==25) {
	    if ($max_ESD[$k][1]<2.82) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==26) {
	    if ($max_ESD[$k][1]<2.84) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==27) {
	    if ($max_ESD[$k][1]<2.86) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==28) {
	    if ($max_ESD[$k][1]<2.88) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==29) {
	    if ($max_ESD[$k][1]<2.89) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n==30) {
	    if ($max_ESD[$k][1]<2.91) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=31) && ($n<=35)) {
	    if ($max_ESD[$k][1]<2.98) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=36) && ($n<=40)) {
	    if ($max_ESD[$k][1]<3.04) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=41) && ($n<=45)) {
	    if ($max_ESD[$k][1]<3.09) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=46) && ($n<=50)) {
	    if ($max_ESD[$k][1]<3.13) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=51) && ($n<=60)) {
	    if ($max_ESD[$k][1]<3.20) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=61) && ($n<=70)) {
	    if ($max_ESD[$k][1]<3.26) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=71) && ($n<=80)) {
	    if ($max_ESD[$k][1]<3.31) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=81) && ($n<=90)) {
	    if ($max_ESD[$k][1]<3.35) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=91) && ($n<=100)) {
	    if ($max_ESD[$k][1]<3.38) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=101) && ($n<=150)) {
	    if ($max_ESD[$k][1]<3.52) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=151) && ($n<=200)) {
	    if ($max_ESD[$k][1]<3.61) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=201) && ($n<=300)) {
	    if ($max_ESD[$k][1]<3.72) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif (($n>=301) && ($n<=400)) {
	    if ($max_ESD[$k][1]<3.80) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	elsif ($n>=401) {
	    if ($max_ESD[$k][1]<(3.86+0.06*int(($n-400)/100))) {
		return(0.1);
	    }
	    else {
		return(0.049);
	    }		
	}
	else {
	    return(0.9);
	}
    }
    
    my @ESD_pvalue=();
    $ESD_pavlue[5]=1.72;
    $ESD_pavlue[6]=1.89;   
    $ESD_pavlue[7]=2.02;
    $ESD_pavlue[8]=2.13;
    $ESD_pavlue[9]=2.21;
    $ESD_pavlue[10]=2.29;
    $ESD_pavlue[11]=2.36;
    $ESD_pavlue[12]=2.41;
    $ESD_pavlue[13]=2.46;
    $ESD_pavlue[14]=2.51;
    $ESD_pavlue[15]=2.55;
    $ESD_pavlue[16]=2.59;
    $ESD_pavlue[17]=2.62;
    $ESD_pavlue[18]=2.65;
    $ESD_pavlue[19]=2.68;
    $ESD_pavlue[20]=2.71;
    $ESD_pavlue[21]=2.73;
    $ESD_pavlue[22]=2.76;
    $ESD_pvalue[23]=2.78;
    $ESD_pvalue[24]=2.80;

    $ESD_pavlue[25]=2.82;
    $ESD_pavlue[26]=2.84;   
    $ESD_pavlue[27]=2.86;
    $ESD_pavlue[28]=2.88;
    $ESD_pavlue[29]=2.89;
    $ESD_pavlue[30]=2.91;
    $ESD_pavlue[35]=2.98;
    $ESD_pavlue[40]=3.04;
    $ESD_pavlue[45]=3.09;
    $ESD_pavlue[50]=3.13;
    $ESD_pavlue[60]=3.20;
    $ESD_pavlue[70]=3.26;
    $ESD_pavlue[80]=3.31;
    $ESD_pavlue[90]=3.35;
    $ESD_pavlue[100]=3.38;
    $ESD_pavlue[150]=3.52;
    $ESD_pavlue[200]=3.61;
    $ESD_pavlue[300]=3.72;
    $ESD_pvalue[400]=3.80;
    $ESD_pvalue[500]=3.86;

}

sub quicksort_lookup_table_each_beforeCluster {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1, $temp_2);
    
    if ($right>=$left) {
	$pivot=$temp_ESD_value[$right][1];
	# print "pivot= $pivot\n";
	$pointer_left=$left-1;
	$pointer_right=$right;
	do
	{
	    do
	    {
		$pointer_left=$pointer_left+1;
	    } until ($temp_ESD_value[$pointer_left][1]>=$pivot);
	    
	    do
	    { 
		$pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($temp_ESD_value[$pointer_right][1]<=$pivot));
	    
	    $temp_1=$temp_ESD_value[$pointer_left][1];
	    $temp_ESD_value[$pointer_left][1]=$temp_ESD_value[$pointer_right][1];
	    $temp_ESD_value[$pointer_right][1]=$temp_1;
	    
	    $temp_2=$temp_ESD_value[$pointer_left][2];
	    $temp_ESD_value[$pointer_left][2]=$temp_ESD_value[$pointer_right][2];
	    $temp_ESD_value[$pointer_right][2]=$temp_2;
	    
	} until ($pointer_right<=$pointer_left);
	
	$temp_ESD_value[$pointer_right][1]=$temp_ESD_value[$pointer_left][1];
	$temp_ESD_value[$pointer_left][1]=$temp_ESD_value[$right][1];
	$temp_ESD_value[$right][1]=$temp_1;
	
	$temp_ESD_value[$pointer_right][2]=$temp_ESD_value[$pointer_left][2];
	$temp_ESD_value[$pointer_left][2]=$temp_ESD_value[$right][2];
	$temp_ESD_value[$right][2]=$temp_2;
	
	
	quicksort_lookup_table_each_beforeCluster($left, $pointer_left-1);
	quicksort_lookup_table_each_beforeCluster($pointer_left+1, $right);
    }
}

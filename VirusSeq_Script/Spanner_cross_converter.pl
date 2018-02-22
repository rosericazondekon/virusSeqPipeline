#! /usr/bin/perl

####! /usr/bin/perl -w
#########The Spanner file MSK.anchors.txt is extemely important
#########the file Spanner_anchor_hg19.txt has a different format from MSK.anchors.txt
#########even though they are basically same
#####chr25 is named for chrVirus in the combined reference hg19OV genome 

# my $in_data_file=$ARGV[0];  #"hg19OV_refGene_RIS.txt";
# my $axt_file=$ARGV[1];      #axt file
# my $out_chip_file=$ARGV[2];   #out crossChip filename  

#start to deal with UCSC refGene data
get_hgGenome_refGene();
regenerate_hgGenome_refGene_from_fiveprime_to_threeprime();


sort_hgGenome_refGene_by_geneName();
#remove the redundant gene symbols
remove_redundant_genesymbols_byExonCount_GeneSize();

get_chromosome_genomic_location_hgGenome_refGene();
generate_hgGenome_refGene_Pointer_by_chrom();


#Spanner data
#######The Spanner file MSK.anchors.txt is extemely important
get_Spanner_anchor_data();

get_spanner_crossdata_data_by_batch();
#swap_length_with_orientation();

#annotate Spanner data
annotate_cross_data_anchor0_with_refGene();
annotate_cross_data_anchor1_with_refGene();

#sort the crossdata by anchor0 genename
sort_crossdata_by_anchor0();

#sort the crossdata by anchor1 genename after anchor0 genename
resort_crossdata_by_anchor();

#aggegate the cross data
aggregate_crossdata_by_anchor0_and_anchor1();

#output_annotate_cross_data();
#output_annotate_aggregate_cross_data();
output_annotate_aggregate_cross_data_exclude_unkn();

#map readID and sequence into cross table
build_hash_for_reduced_cross_data();
load_unique_aligned_sorted_axt_file();
map_location_backto_readID();

#start blat filerting
replace_blat_parsing();
#output_fasta_for_blat_from_reduced_annotate_cross_data();
#call_blat_alignment();
#get_blat_data();
filter_blat_reduced_cross_data_for_chip();
output_reduced_annotate_cross_data();

#start to deal with UCSC refGene data
sub get_hgGenome_refGene {
    my $in_data_file=$ARGV[0];
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
   # }
   # print "\n";
  #}

}

sub regenerate_hgGenome_refGene_from_fiveprime_to_threeprime {
  
  for (my $i=1; $i<=$hgGenome_refGene_row_number; $i++) {
    $hgGenome_refGene[$i][3]=~s/chr//;
    if ($hgGenome_refGene[$i][3] eq "X") { $hgGenome_refGene[$i][3]=23; }
    elsif ($hgGenome_refGene[$i][3] eq "Y") { $hgGenome_refGene[$i][3]=24; }
    elsif ($hgGenome_refGene[$i][3] eq "M") { $hgGenome_refGene[$i][3]=25; }

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
  #  }
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
    #  }
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


sub get_spanner_crossdata_data_by_batch {
    #upload the MSK.*.cross.span.txt by Spanner
    
    #my $VirusOnly_indicator=$ARGV[2]; #detect only Virus or other fusion
    @cross_data=();
    $cross_row_number=0;
    %anchor_indicator=();
    
    my $start_chr=25;
    #if ($VirusOnly_indicator=~/Yes/) { $start_chr=25; }
    
    for (my $chrID=$start_chr; $chrID<=25; $chrID++) {
	my @temp_data=();
	undef(@temp_data);
	@temp_data=();
	my $chr_number=$chrID;
	if ($chrID=~/23/) { $chr_number="X"; }
	elsif ($chrID=~/24/){ $chr_number="Y"; }
	elsif ($chrID=~/25/) { $chr_number="M"; }
	
	$in_data_file="MSK"."\.".$chr_number."\.cross\.span\.txt";	
	open (INPUT, "<$in_data_file")
	    || die "Can't open $in_data_file $!";
	
	
	$temp_data_row_number=0; #genes number
	$line_indicator=0;
	#input data from files
	
	while (<INPUT>) {
	    $text=$_;
	    $text=~s/|\n//g;
	    $text=~s/ //;
	    (@line)=split(" ", $text);
	    
	    if ($text=~/position/) {
		$line_indicator++;
	    }
	    
	    #readin the data
	    if ($line_indicator>0) {
		$temp_data_column_number=0;
		for $read_data(@line) {
		    #$read_data=~s/\"//g;
		    $temp_data[$temp_data_row_number][$temp_data_column_number]=$read_data;
		    $temp_data_column_number++;
		}	
		
		$temp_data_row_number++;
		#print "column=$temp_data_column_number\n";
	    }
	} #while(INPUT)
	
	close(INPUT);	
    
	$temp_data_row_number--;
	$temp_data_column_number--;

	#header
	if ($chrID==$start_chr) {
	    for (my $j=0; $j<=$temp_data_column_number; $j++) {
		$cross_data[0][$j]=$temp_data[0][$j];
	    }
	} #if ($chrID==$start_chr)
	
	#put the data into single array
	#content
	for (my $i=1; $i<=$temp_data_row_number; $i++) {
	    my $anchor_0=$temp_data[$i][0]."\_".$temp_data[$i][1];
	    my $anchor_1=$temp_data[$i][5]."\_".$temp_data[$i][6];
	    my $anchor_0_1=$anchor_0."\_".$anchor_1;
	    my $anchor_1_0=$anchor_1."\_".$anchor_0;
	    $anchor_0_1=~s/ //g;
	    $anchor_1_0=~s/ //g;
	    
	    if ((!(exists($anchor_indicator{$anchor_0_1}))) && (!(exists($anchor_indicator{$anchor_1_0})))) {
		$anchor_indicator{$anchor_0_1}=1;
		$anchor_indicator{$anchor_1_0}=1;
		
		$cross_row_number++;
		for (my $j=0; $j<=$temp_data_column_number; $j++) {	    
		    $cross_data[$cross_row_number][$j]=$temp_data[$i][$j];
		}
	    }
	    #else {
		#print "$anchor_0_1   Indicator\n";
		#print "$anchor_1_0   Indicator\n";
	    #}
	} 
	
	undef(@temp_data);
    } #for (my $crID=1; $chrID<=25; $chrID++)
    
    $cross_column_number=$temp_data_column_number;
    #finishing loading thw multiple tables
    

    #convert anchor number back to chromosome number
    for (my $i=1; $i<=$cross_row_number; $i++) {
	$cross_data[$i][0]=$anchor_to_chr{$cross_data[$i][0]};
	$cross_data[$i][5]=$anchor_to_chr{$cross_data[$i][5]};	
    }
    
    
    #remove random
    $new_row=0;
    @new_cross_data=();
    for (my $i=1; $i<=$cross_row_number; $i++) {   
	#if (!(($cross_data[$i][5]=~/random/) || ($cross_data[$i][5]=~/hap/) || ($cross_data[$i][5]=~/M/))) {
	if (!(($cross_data[$i][5]=~/random/) || ($cross_data[$i][5]=~/hap/))) {
	    $new_row++;
	    for (my $j=0; $j<=$cross_column_number; $j++) {
		$new_cross_data[$new_row][$j]=$cross_data[$i][$j];
	    }
	}
    }
    
    #header
    for (my $j=0; $j<=$cross_column_number; $j++) {
	$new_cross_data[0][$j]=$cross_data[0][$j];
    }

    #change X and Y to 23 24
    for (my $i=1; $i<=$new_row; $i++) {
	if ($new_cross_data[$i][0]=~/X/) {
	    $new_cross_data[$i][0]=23;
	}	
	elsif ($new_cross_data[$i][0]=~/Y/) {
	    $new_cross_data[$i][0]=24;
	}
	elsif ($new_cross_data[$i][0]=~/M/) {
	    $new_cross_data[$i][0]=25;
	}
	
	if ($new_cross_data[$i][5]=~/X/) {
	    $new_cross_data[$i][5]=23;
	}
	elsif ($new_cross_data[$i][5]=~/Y/) {
	    $new_cross_data[$i][5]=24;
	}    
	elsif ($new_cross_data[$i][5]=~/M/) {
	    $new_cross_data[$i][5]=25;
	}    
    }
    
    #re-align the location back to Mosaik location
    for (my $i=1; $i<=$new_row; $i++) {
	$new_cross_data[$i][1]++;
	$new_cross_data[$i][6]++;	
	
    }
    
    #re-assign the data
    $cross_row_number=$new_row;
    undef(@cross_data);
    @cross_data=@new_cross_data;
    undef(@new_cross_data);
    
    
    #for (my $i=0; $i<=$cross_row_number; $i++) {
    #for (my $j=0; $j<=$cross_column_number; $j++) {
	 #  print "$cross_data[$i][$j]\t";
    #}
    #print "\n";
    #}

}


sub get_Spanner_anchor_data {
    my $in_data_file="MSK.anchors.txt";
    
    open (INPUT, "<$in_data_file")
	|| die "Can't open $in_data_file $!";
    
    
    @anchor_data=();
    
    $anchor_row_number=0; #genes number
    $line_indicator=0;
    #input data from files
    
    while (<INPUT>) {
	$text=$_;
	$text=~s/|\n//g;
	$text=~s/\t//;
	(@line)=split("\t", $text);
	
	if ($text=~/0\./) {
	    $line_indicator++;
	}
	
	#readin the data
	if ($line_indicator>0) {
	    $anchor_column_number=0;
	    for $read_data(@line) {
		#$read_data=~s/\"//g;
		$anchor_data[$anchor_row_number][$anchor_column_number]=$read_data;
		$anchor_column_number++;
	    }	
	    
	    $anchor_row_number++;
	    #print "column=$anchor_column_number\n";
	}
    } #while(INPUT)
    
    close(INPUT);
    
    $anchor_row_number--;
    $anchor_column_number--;
    
    %anchor_to_chr=();
    for (my $i=0; $i<$anchor_row_number; $i++) {
	$anchor_data[$i][0]=~s/\.//;
	$anchor_to_chr{$anchor_data[$i][0]}=$anchor_data[$i][1];
	#print "$anchor_data[$i][0]=$anchor_data[$i][1]\n";
    }
}

sub annotate_cross_data_anchor0_with_refGene {

    my $aligned_read_length=60; #read_length=60bp
    for (my $i=1; $i<=$cross_row_number; $i++) {
	my $chromosome=$cross_data[$i][0];
	my $start_location=$cross_data[$i][1]+1;
	#print "chr=$chromosome\n";
	#print "start=$start_location\n";
	
	my $hgGenome_start_pointer=$hgGenome_refGene_pointer[$chromosome][1];
	my $hgGenome_end_pointer=$hgGenome_refGene_pointer[$chromosome][2];
	
	my $final_down_index=detect_start_location_hgGenome($start_location,
							$hgGenome_start_pointer, $hgGenome_end_pointer);
	
	$cross_data[$i][$cross_column_number+1]="Unkn".$chromosome;
	
	if (($start_location>=$hgGenome_refGene[$final_down_index][5]) && 
	    ($start_location<=$hgGenome_refGene[$final_down_index][6])) {
	    $cross_data[$i][$cross_column_number+1]=$hgGenome_refGene[$final_down_index][7];

	    my $first_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location);
	    my $second_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location+$aligned_read_length);
	    my $third_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location-$aligned_read_length);
	    if ($first_location=~/exon/) {
		$cross_data[$i][$cross_column_number+1]=
		    $cross_data[$i][$cross_column_number+1]."\/\/".$first_location;
	    }
	    elsif ($second_location=~/exon/) {
		$cross_data[$i][$cross_column_number+1]=
		    $cross_data[$i][$cross_column_number+1]."\/\/".$second_location;
	    }
	    else {
		$cross_data[$i][$cross_column_number+1]=
		    $cross_data[$i][$cross_column_number+1]."\/\/".$third_location;
	    }
	    #$cross_data[$i][$cross_column_number+3]=
		#detect_location_site_within_gene_byLocation($final_down_index, $start_location);
	}
	
	#move up the pointer by 5
	if ($cross_data[$i][$cross_column_number+1]=~/Unkn/) {
	    my $temp_pointer=$final_down_index-1;
	    while ($temp_pointer>=($final_down_index-5)) {
		if ($temp_pointer<=1) {
		    last;
		    $temp_pointer=2;
		}
		if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		    ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		    $cross_data[$i][$cross_column_number+1]=$hgGenome_refGene[$temp_pointer][7];

		    my $first_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location);
		    my $second_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location+$aligned_read_length);
		    my $third_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location-$aligned_read_length);
		    if ($first_location=~/exon/) {
			$cross_data[$i][$cross_column_number+1]=
			    $cross_data[$i][$cross_column_number+1]."\/\/".$first_location;
		    }
		    elsif ($second_location=~/exon/) {
			$cross_data[$i][$cross_column_number+1]=
			    $cross_data[$i][$cross_column_number+1]."\/\/".$second_location;
		    }
		    else {
			$cross_data[$i][$cross_column_number+1]=
			    $cross_data[$i][$cross_column_number+1]."\/\/".$third_location;
		    }
		    
		    #$cross_data[$i][$cross_column_number+3]=
			#detect_location_site_within_gene_byLocation($temp_pointer, $start_location);
		    last;
		}
		#move the pointer
		$temp_pointer--;
	    }  
	}
	
	#move down the pointer by 5
	if ($cross_data[$i][$cross_column_number+1]=~/Unkn/) {
	    my $temp_pointer=$final_down_index+1;
	    while ($temp_pointer<=($final_down_index+5)) {
		if ($temp_pointer>=$hgGenome_refGene_row_number) {
		    last;
		    $temp_pointer--;
		}
		
		if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		    ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		    $cross_data[$i][$cross_column_number+1]=$hgGenome_refGene[$temp_pointer][7];

		    my $first_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location);
		    my $second_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location+$aligned_read_length);
		    my $third_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location-$aligned_read_length);
		    if ($first_location=~/exon/) {
			$cross_data[$i][$cross_column_number+1]=
			    $cross_data[$i][$cross_column_number+1]."\/\/".$first_location;
		    }
		    elsif ($second_location=~/exon/) {
			$cross_data[$i][$cross_column_number+1]=
			    $cross_data[$i][$cross_column_number+1]."\/\/".$second_location;
		    }
		    else {
			$cross_data[$i][$cross_column_number+1]=
			    $cross_data[$i][$cross_column_number+1]."\/\/".$third_location;
		    }
		    #$cross_data[$i][$cross_column_number+3]=
			#detect_location_site_within_gene_byLocation($temp_pointer, $start_location);
		    last;
		}
		#move the pointer
		$temp_pointer++;
	    }  
	}
	
	#assign a nearby gene
	if ($cross_data[$i][$cross_column_number+1]=~/Unkn/) {
	    if ($final_down_index<1) {
		$final_down_index=1;
	    }
	    if ($final_down_index>$hgGenome_refGene_row_number) {
		$final_down_index=$hgGenome_refGene_row_number;
	    }

	    $cross_data[$i][$cross_column_number+1]=$hgGenome_refGene[$final_down_index][7];
	    my $first_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location);
	    my $second_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location+$aligned_read_length);
	    my $third_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location-$aligned_read_length);
	    if ($first_location=~/exon/) {
		$cross_data[$i][$cross_column_number+1]=
		    $cross_data[$i][$cross_column_number+1]."\/\/".$first_location;
	    }
	    elsif ($second_location=~/exon/) {
		$cross_data[$i][$cross_column_number+1]=
		    $cross_data[$i][$cross_column_number+1]."\/\/".$second_location;
	    }
	    else {
		$cross_data[$i][$cross_column_number+1]=
		    $cross_data[$i][$cross_column_number+1]."\/\/".$third_location;
	    }
	}
	
    } #  for (my $i=1; $i<=$cross_row_number; $i++) {
    
}

sub annotate_cross_data_anchor1_with_refGene {
    my $aligned_read_length=60; #read_length=60bp
    for (my $i=1; $i<=$cross_row_number; $i++) {
	my $chromosome=$cross_data[$i][5];
	my $start_location=$cross_data[$i][6]+1;

	
	my $hgGenome_start_pointer=$hgGenome_refGene_pointer[$chromosome][1];
	my $hgGenome_end_pointer=$hgGenome_refGene_pointer[$chromosome][2];
	
	my $final_down_index=detect_start_location_hgGenome($start_location,
							$hgGenome_start_pointer, $hgGenome_end_pointer);

	$cross_data[$i][$cross_column_number+2]="Unkn".$chromosome;
	
	if (($start_location>=$hgGenome_refGene[$final_down_index][5]) && 
	    ($start_location<=$hgGenome_refGene[$final_down_index][6])) {
	    $cross_data[$i][$cross_column_number+2]=$hgGenome_refGene[$final_down_index][7];

	    my $first_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location);
	    my $second_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location+$aligned_read_length);
	    my $third_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location-$aligned_read_length);
	    if ($first_location=~/exon/) {
		$cross_data[$i][$cross_column_number+2]=
		    $cross_data[$i][$cross_column_number+2]."\/\/".$first_location;
	    }
	    elsif ($second_location=~/exon/) {
		$cross_data[$i][$cross_column_number+2]=
		    $cross_data[$i][$cross_column_number+2]."\/\/".$second_location;
	    }
	    else {
		$cross_data[$i][$cross_column_number+2]=
		    $cross_data[$i][$cross_column_number+2]."\/\/".$third_location;
	    }
	    #$cross_data[$i][$cross_column_number+4]=
		#detect_location_site_within_gene_byLocation($final_down_index, $start_location);
	}

	#move up the pointer by 5
	if ($cross_data[$i][$cross_column_number+2]=~/Unkn/) {
	    my $temp_pointer=$final_down_index-1;
	    while ($temp_pointer>=($final_down_index-5)) {
		if ($temp_pointer<=1) {
		    last;
		    $temp_pointer=2;
		}
		if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		    ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		    $cross_data[$i][$cross_column_number+2]=$hgGenome_refGene[$temp_pointer][7];

		    my $first_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location);
		    my $second_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location+$aligned_read_length);
		    my $third_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location-$aligned_read_length);
		    if ($first_location=~/exon/) {
			$cross_data[$i][$cross_column_number+2]=
			    $cross_data[$i][$cross_column_number+2]."\/\/".$first_location;
		    }
		    elsif ($second_location=~/exon/) {
			$cross_data[$i][$cross_column_number+2]=
			    $cross_data[$i][$cross_column_number+2]."\/\/".$second_location;
		    }
		    else {
			$cross_data[$i][$cross_column_number+2]=
			    $cross_data[$i][$cross_column_number+2]."\/\/".$third_location;
		    }
		    
		    #$cross_data[$i][$cross_column_number+4]=
			#detect_location_site_within_gene_byLocation($temp_pointer, $start_location);
		    last;
		}
		#move the pointer
		$temp_pointer--;
	    }  
	}
	
	#move down the pointer by 5
	if ($cross_data[$i][$cross_column_number+2]=~/Unkn/) {
	    my $temp_pointer=$final_down_index+1;
	    while ($temp_pointer<=($final_down_index+5)) {
		if ($temp_pointer>=$hgGenome_refGene_row_number) {
		    last;
		    $temp_pointer--;
		}
		
		if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		    ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		    $cross_data[$i][$cross_column_number+2]=$hgGenome_refGene[$temp_pointer][7];
		    my $first_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location);
		    my $second_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location+$aligned_read_length);
		    my $third_location=detect_location_site_within_gene_byLocation($temp_pointer, $start_location-$aligned_read_length);
		    if ($first_location=~/exon/) {
			$cross_data[$i][$cross_column_number+2]=
			    $cross_data[$i][$cross_column_number+2]."\/\/".$first_location;
		    }
		    elsif ($second_location=~/exon/) {
			$cross_data[$i][$cross_column_number+2]=
			    $cross_data[$i][$cross_column_number+2]."\/\/".$second_location;
		    }
		    else {
			$cross_data[$i][$cross_column_number+2]=
			    $cross_data[$i][$cross_column_number+2]."\/\/".$third_location;
		    }

		    #$cross_data[$i][$cross_column_number+4]=
			#detect_location_site_within_gene_byLocation($temp_pointer, $start_location);
		    last;
		}
		#move the pointer
		$temp_pointer++;
	    }  
	}
	
	#assign a nearby gene
	if ($cross_data[$i][$cross_column_number+2]=~/Unkn/) {
	    if ($final_down_index<1) {
		$final_down_index=1;
	    }
	    if ($final_down_index>$hgGenome_refGene_row_number) {
		$final_down_index=$hgGenome_refGene_row_number;
	    }
	    $cross_data[$i][$cross_column_number+2]=$hgGenome_refGene[$final_down_index][7];

	    my $first_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location);
	    my $second_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location+$aligned_read_length);
	    my $third_location=detect_location_site_within_gene_byLocation($final_down_index, $start_location-$aligned_read_length);
	    if ($first_location=~/exon/) {
		$cross_data[$i][$cross_column_number+2]=
		    $cross_data[$i][$cross_column_number+2]."\/\/".$first_location;
	    }
	    elsif ($second_location=~/exon/) {
		$cross_data[$i][$cross_column_number+2]=
		    $cross_data[$i][$cross_column_number+2]."\/\/".$second_location;
	    }
	    else {
		$cross_data[$i][$cross_column_number+2]=
		    $cross_data[$i][$cross_column_number+2]."\/\/".$third_location;
	    }
	    
	}
	
    } #  for (my $i=1; $i<=$cross_row_number; $i++) {

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

sub annotate_cross_data_anchor0_with_refGene_Old {

    for (my $i=1; $i<=$cross_row_number; $i++) {
	my $chromosome=$cross_data[$i][0];
	my $start_location=$cross_data[$i][1]+1;
	#print "chr=$chromosome\n";
	#print "start=$start_location\n";
	
	my $hgGenome_start_pointer=$hgGenome_refGene_pointer[$chromosome][1];
	my $hgGenome_end_pointer=$hgGenome_refGene_pointer[$chromosome][2];
	
	my $final_down_index=detect_start_location_hgGenome($start_location,
							$hgGenome_start_pointer, $hgGenome_end_pointer);
	
	$cross_data[$i][$cross_column_number+1]="Unkn".$chromosome;
	
	if (($start_location>=$hgGenome_refGene[$final_down_index][5]) && 
	    ($start_location<=$hgGenome_refGene[$final_down_index][6])) {
	    $cross_data[$i][$cross_column_number+1]=$hgGenome_refGene[$final_down_index][7];
	}
	
	#move up the pointer by 5
	if ($cross_data[$i][$cross_column_number+1]=~/Unkn/) {
	    my $temp_pointer=$final_down_index-1;
	    while ($temp_pointer>=($final_down_index-5)) {
		if ($temp_pointer<=1) {
		    last;
		    $temp_pointer=2;
		}
		if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		    ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		    $cross_data[$i][$cross_column_number+1]=$hgGenome_refGene[$temp_pointer][7];
		    last;
		}
		#move the pointer
		$temp_pointer--;
	    }  
	}
	
	#move down the pointer by 5
	if ($cross_data[$i][$cross_column_number+1]=~/Unkn/) {
	    my $temp_pointer=$final_down_index+1;
	    while ($temp_pointer<=($final_down_index+5)) {
		if ($temp_pointer>=$hgGenome_refGene_row_number) {
		    $temp_pointer--;
		    last;
		}
		
		if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		    ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		    $cross_data[$i][$cross_column_number+1]=$hgGenome_refGene[$temp_pointer][7];
		    last;
		}
		#move the pointer
		$temp_pointer++;
	    }  
	}
	
	#assign a nearby gene
	if ($cross_data[$i][$cross_column_number+1]=~/Unkn/) {
	    if ($final_down_index<1) {
		$final_down_index=1;
	    }
	    if ($final_down_index>$hgGenome_refGene_row_number) {
		$final_down_index=$hgGenome_refGene_row_number;
	    }
	    $cross_data[$i][$cross_column_number+1]=$hgGenome_refGene[$final_down_index][7];
	}
	
    } #  for (my $i=1; $i<=$cross_row_number; $i++) {

}

sub annotate_cross_data_anchor1_with_refGene_Old {

    for (my $i=1; $i<=$cross_row_number; $i++) {
	my $chromosome=$cross_data[$i][5];
	my $start_location=$cross_data[$i][6]+1;

	
	my $hgGenome_start_pointer=$hgGenome_refGene_pointer[$chromosome][1];
	my $hgGenome_end_pointer=$hgGenome_refGene_pointer[$chromosome][2];
	
	my $final_down_index=detect_start_location_hgGenome($start_location,
							$hgGenome_start_pointer, $hgGenome_end_pointer);

	$cross_data[$i][$cross_column_number+2]="Unkn".$chromosome;
	
	if (($start_location>=$hgGenome_refGene[$final_down_index][5]) && 
	    ($start_location<=$hgGenome_refGene[$final_down_index][6])) {
	    $cross_data[$i][$cross_column_number+2]=$hgGenome_refGene[$final_down_index][7];
	}

	#move up the pointer by 5
	if ($cross_data[$i][$cross_column_number+2]=~/Unkn/) {
	    my $temp_pointer=$final_down_index-1;
	    while ($temp_pointer>=($final_down_index-5)) {
		if ($temp_pointer<=1) {
		    last;
		    $temp_pointer=2;
		}
		if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		    ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		    $cross_data[$i][$cross_column_number+2]=$hgGenome_refGene[$temp_pointer][7];
		    last;
		}
		#move the pointer
		$temp_pointer--;
	    }  
	}
	
	#move down the pointer by 5
	if ($cross_data[$i][$cross_column_number+2]=~/Unkn/) {
	    my $temp_pointer=$final_down_index+1;
	    while ($temp_pointer<=($final_down_index+5)) {
		if ($temp_pointer>=$hgGenome_refGene_row_number) {
		    last;
		    $temp_pointer--;
		}
		
		if (($start_location>=$hgGenome_refGene[$temp_pointer][5]) && 
		    ($start_location<=$hgGenome_refGene[$temp_pointer][6])) {		
		    $cross_data[$i][$cross_column_number+2]=$hgGenome_refGene[$temp_pointer][7];
		    last;
		}
		#move the pointer
		$temp_pointer++;
	    }  
	}
	
	
	#assign a nearby gene
	if ($cross_data[$i][$cross_column_number+2]=~/Unkn/) {
	    if ($final_down_index<1) {
		$final_down_index=1;
	    }
	    if ($final_down_index>$hgGenome_refGene_row_number) {
		$final_down_index=$hgGenome_refGene_row_number;
	    }
	    $cross_data[$i][$cross_column_number+2]=$hgGenome_refGene[$final_down_index][7];
	}
	
    } #  for (my $i=1; $i<=$cross_row_number; $i++) {
    
}


sub output_annotate_cross_data {
    #output the data

    $cross_data[0][$cross_column_number+1]="Gene0";
    $cross_data[0][$cross_column_number+2]="Gene1";
    
    for (my $i=0; $i<=$cross_row_number; $i++) {
	for (my $j=0; $j<=($cross_column_number+1); $j++) {
	    print "$cross_data[$i][$j]\t";
	}
	print "$cross_data[$i][$cross_column_number+2]\n";
    }
    
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

sub sort_crossdata_by_anchor0 {
    
    #sort the data by anchor1 geneName
    @sort_cross_index=();
    for (my $i=1; $i<=$cross_row_number; $i++) {
	$sort_cross_index[$i][1]=$cross_data[$i][$cross_column_number+1];
	$sort_cross_index[$i][2]=$i;
    }


    quicksort_crossdata_by_geneName(1, $cross_row_number);

    #re-order the data by geneName
    my @temp_sort_cross_data=();
    for (my $i=1; $i<=$cross_row_number; $i++) {
	my $row_index=$sort_cross_index[$i][2];
	for (my $j=0; $j<=($cross_column_number+2); $j++) {
	    $temp_sort_cross_data[$i][$j]=$cross_data[$row_index][$j];
	}
    }

    #header
    $cross_data[0][$cross_column_number+1]="Gene0";
    $cross_data[0][$cross_column_number+2]="Gene1";
    for (my $j=0; $j<=($cross_column_number+2); $j++) {
	$temp_sort_cross_data[0][$j]=$cross_data[0][$j];
    }

    undef(@cross_data);
    @cross_data=();
    @cross_data=@temp_sort_cross_data;  
    undef(@temp_sort_cross_data);

}

sub quicksort_crossdata_by_geneName {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1, $temp_2);

    if ($right>=$left) {
	$pivot=$sort_cross_index[$right][1];
	# print "pivot= $pivot\n";
	$pointer_left=$left-1;
	$pointer_right=$right;
	do
	{
	    do
	    {
		$pointer_left=$pointer_left+1;
	    } until ($sort_cross_index[$pointer_left][1] ge $pivot);
	    
	    do
	    { 
		$pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_cross_index[$pointer_right][1] le $pivot));
	    
	    $temp_1=$sort_cross_index[$pointer_left][1];
	    $sort_cross_index[$pointer_left][1]=$sort_cross_index[$pointer_right][1];
	    $sort_cross_index[$pointer_right][1]=$temp_1;
	    
	    $temp_2=$sort_cross_index[$pointer_left][2];
	    $sort_cross_index[$pointer_left][2]=$sort_cross_index[$pointer_right][2];
	    $sort_cross_index[$pointer_right][2]=$temp_2;
	    
	} until ($pointer_right<=$pointer_left);
	
	$sort_cross_index[$pointer_right][1]=$sort_cross_index[$pointer_left][1];
	$sort_cross_index[$pointer_left][1]=$sort_cross_index[$right][1];
	$sort_cross_index[$right][1]=$temp_1;
	
	$sort_cross_index[$pointer_right][2]=$sort_cross_index[$pointer_left][2];
	$sort_cross_index[$pointer_left][2]=$sort_cross_index[$right][2];
	$sort_cross_index[$right][2]=$temp_2;
	
	
	quicksort_crossdata_by_geneName($left, $pointer_left-1);
	quicksort_crossdata_by_geneName($pointer_left+1, $right);
    }
}


sub resort_crossdata_by_anchor {
    
    my $temp_start_pointer=1;
    my $temp_end_pointer=1;
    my $gene_ID=$cross_data[$temp_start_pointer][$cross_column_number+1];
    
    my $moving_pointer=0;
    while ($moving_pointer<$cross_row_number) {
	while (($moving_pointer<$cross_row_number) && 
	       ($gene_ID eq $cross_data[$moving_pointer+1][$cross_column_number+1])) {
	    $moving_pointer++;
	}
	
	#start to assign genename the ordinal number
	$temp_end_pointer=$moving_pointer;
	
	
	#re-sort the data by anchor1 genename
	my @temp_anchorGene_data=();
	my $temp_anchorGene_number=0;
	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	    $temp_anchorGene_number++;
	    for (my $j=0; $j<=($cross_column_number+2); $j++) {
		$temp_anchorGene_data[$temp_anchorGene_number][$j]=$cross_data[$h][$j];
	    }
	}

	@sort_anchorGene_index=();
	for (my $h=1; $h<=$temp_anchorGene_number; $h++) {
	    $sort_anchorGene_index[$h][1]=$temp_anchorGene_data[$h][$cross_column_number+2];
	    $sort_anchorGene_index[$h][2]=$h;
	}
	quicksort_by_anchor_geneName(1, $temp_anchorGene_number);
	
	#put re-sorted data back into cross-anchordata
	for (my $h=1; $h<=$temp_anchorGene_number; $h++) {
	    my $row_index=$sort_anchorGene_index[$h][2];
	    my $true_row_number=$temp_start_pointer+$h-1;
	    for (my $j=0; $j<=($cross_column_number+2); $j++) {
		$cross_data[$true_row_number][$j]=$temp_anchorGene_data[$row_index][$j];
	    }
	}
	
	#re-start the pointer;
	$moving_pointer++;
	$temp_start_pointer=$moving_pointer;
	$temp_end_pointer=$moving_pointer;
	$gene_ID=$cross_data[$temp_start_pointer][$cross_column_number+1];
	
    }
 
}

sub quicksort_by_anchor_geneName {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1, $temp_2);
    
    if ($right>=$left) {
	$pivot=$sort_anchorGene_index[$right][1];
	# print "pivot= $pivot\n";
	$pointer_left=$left-1;
	$pointer_right=$right;
	do
	{
	    do
	    {
		$pointer_left=$pointer_left+1;
	    } until ($sort_anchorGene_index[$pointer_left][1] ge $pivot);
	    
	    do
	    { 
		$pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_anchorGene_index[$pointer_right][1] le $pivot));
	    
	    $temp_1=$sort_anchorGene_index[$pointer_left][1];
	    $sort_anchorGene_index[$pointer_left][1]=$sort_anchorGene_index[$pointer_right][1];
	    $sort_anchorGene_index[$pointer_right][1]=$temp_1;
	    
	    $temp_2=$sort_anchorGene_index[$pointer_left][2];
	    $sort_anchorGene_index[$pointer_left][2]=$sort_anchorGene_index[$pointer_right][2];
	    $sort_anchorGene_index[$pointer_right][2]=$temp_2;
	    
	} until ($pointer_right<=$pointer_left);
	
	$sort_anchorGene_index[$pointer_right][1]=$sort_anchorGene_index[$pointer_left][1];
	$sort_anchorGene_index[$pointer_left][1]=$sort_anchorGene_index[$right][1];
	$sort_anchorGene_index[$right][1]=$temp_1;
	
	$sort_anchorGene_index[$pointer_right][2]=$sort_anchorGene_index[$pointer_left][2];
	$sort_anchorGene_index[$pointer_left][2]=$sort_anchorGene_index[$right][2];
	$sort_anchorGene_index[$right][2]=$temp_2;
	
	
	quicksort_by_anchor_geneName($left, $pointer_left-1);
	quicksort_by_anchor_geneName($pointer_left+1, $right);
    }
}

sub aggregate_crossdata_by_anchor0_and_anchor1 {
 
    @aggregate_cross_data=();
    $aggregate_cross_row_number=0;
    $aggregate_cross_column_number=$cross_column_number;

    for (my $j=0; $j<=($cross_column_number+2); $j++) {
	$aggregate_cross_data[0][$j]=$cross_data[0][$j];
    }
    $aggregate_cross_data[0][$cross_column_number+3]="Mate_Number";
    $aggregate_cross_data[0][$cross_column_number+4]="anchor0_D";
    $aggregate_cross_data[0][$cross_column_number+5]="anchor1_D";

    my $temp_start_pointer=1;
    my $temp_end_pointer=1;
    my $gene_ID0=$cross_data[$temp_start_pointer][$cross_column_number+1];
    my $gene_ID1=$cross_data[$temp_start_pointer][$cross_column_number+2];
    
    my $moving_pointer=0;
    while ($moving_pointer<$cross_row_number) {
	while (($moving_pointer<$cross_row_number) && 
	       ($gene_ID0 eq $cross_data[$moving_pointer+1][$cross_column_number+1]) &&
	       ($gene_ID1 eq $cross_data[$moving_pointer+1][$cross_column_number+2])) {
	    $moving_pointer++;
	}
	
	#start to assign genename the ordinal number
	$temp_end_pointer=$moving_pointer;
	

	#aggregate the cross data
	#anchor_0 location
	#my $anchor0_left_location=$cross_data[$temp_start_pointer][1];
	#my $anchor0_right_location=$cross_data[$temp_start_pointer][1];

	#my $anchor1_left_location_index=$temp_start_pointer;
	#my $anchor1_right_location_index=$temp_start_pointer;
	
	#for (my $h=($temp_start_pointer+1); $h<=$temp_end_pointer; $h++) {
	 #   if ($anchor0_left_location>$cross_data[$h][1]) {
	#	$anchor0_left_location=$cross_data[$h][1];
	#	$anchor1_left_location_index=$h;
	#    }

	    #if ($anchor0_right_location<$cross_data[$h][1]) {
	#	$anchor0_right_location=$cross_data[$h][1];
	#       $anchor1_right_location_index=$h;
	#    }
	#}

	#anchor_1 location
	#my $anchor1_left_location=$cross_data[$anchor1_left_location_index][6];
	#my $anchor1_right_location=$cross_data[$anchor1_right_location_index][6];


	#$anchor1_left_location=$cross_data[$temp_start_pointer][6];
	#$anchor1_right_location=$cross_data[$temp_start_pointer][6];
	#for (my $h=($temp_start_pointer+1); $h<=$temp_end_pointer; $h++) {
	 #   if ($anchor1_left_location>$cross_data[$h][6]) {
	#	$anchor1_left_location=$cross_data[$h][6];
	#    }
	#    if ($anchor1_right_location<$cross_data[$h][6]) {
		#$anchor1_right_location=$cross_data[$h][6];
	    #}
	#}
	
	my $temp_anchorGene_index=int(($temp_end_pointer-$temp_start_pointer)/2)+$temp_start_pointer;
	$aggregate_cross_row_number++;
	
	for (my $j=0; $j<=($cross_column_number+2); $j++) {
	    $aggregate_cross_data[$aggregate_cross_row_number][$j]=$cross_data[$temp_anchorGene_index][$j];
	}
	$aggregate_cross_data[$aggregate_cross_row_number][$cross_column_number+3]=$temp_end_pointer-$temp_start_pointer+1;

	#generate the MAD for each clusters
	#anchor0
	@temp_row_data=();
	for (my $i=$temp_start_pointer; $i<=$temp_end_pointer; $i++) {
	    $temp_row_data[$i-$temp_start_pointer]=$cross_data[$i][1];
	}	    
	$aggregate_cross_data[$aggregate_cross_row_number][$cross_column_number+4]=cross_data_MAD($temp_end_pointer-$temp_start_pointer+1);
	#$aggregate_cross_data[$aggregate_cross_row_number][$cross_column_number+6]=cross_data_mean_stdev($temp_end_pointer-$temp_start_pointer+1);

	#anchor1
	@temp_row_data=();
	for (my $i=$temp_start_pointer; $i<=$temp_end_pointer; $i++) {
	    $temp_row_data[$i-$temp_start_pointer]=$cross_data[$i][6];
	}	    
	$aggregate_cross_data[$aggregate_cross_row_number][$cross_column_number+5]=cross_data_MAD($temp_end_pointer-$temp_start_pointer+1);
	#$aggregate_cross_data[$aggregate_cross_row_number][$cross_column_number+7]=cross_data_mean_stdev($temp_end_pointer-$temp_start_pointer+1);

	#$aggregate_cross_data[$aggregate_cross_row_number][$cross_column_number+4]=$anchor0_right_location-$anchor0_left_location;
	#$aggregate_cross_data[$aggregate_cross_row_number][$cross_column_number+5]=$anchor1_right_location-$anchor1_left_location;
	
	#re-start the pointer;
	$moving_pointer++;
	$temp_start_pointer=$moving_pointer;
	$temp_end_pointer=$moving_pointer;
	$gene_ID0=$cross_data[$temp_start_pointer][$cross_column_number+1];
	$gene_ID1=$cross_data[$temp_start_pointer][$cross_column_number+2];
	
    }
 
}


sub cross_data_mean_stdev {

    my ($temp_cross_row_number)=@_;
    
    if ($temp_cross_row_number>=2) {
	my $temp_mean=0;
	my $temp_square=0;
	my $good_array=0;
	
	for (my $j=0; $j<$temp_cross_row_number; $j++) {
	    $good_array++;
	    $temp_mean=$temp_mean+$temp_row_data[$j];
	    $temp_square=$temp_square+$temp_row_data[$j]*$temp_row_data[$j];
	}


	$true_mean=$temp_mean/$good_array;
	$temp_mean=$temp_mean/$good_array;
	
	#$temp_row_data[$i][$affy_column_number+1]=$temp_mean;	
	$temp_mean=$temp_mean*$temp_mean;
	if ((($temp_square-$good_array*$temp_mean)/($good_array-1))>0) {
	    $CV=sqrt(($temp_square-$good_array*$temp_mean)/($good_array-1))/$true_mean;
	}
	else {
	    $CV=0;
	}
    }
    else {
	$CV=0;
    }

    return($CV);
  
} #end of subroutine


sub cross_data_MAD {
    
    my ($temp_cross_row_number)=@_;
    
    if ($temp_cross_row_number>=2) {
	my @sort_row_data=sort {$a<=>$b} @temp_row_data;
	
	my $m1=$sort_row_data[int($temp_cross_row_number/2)];
	my $m2=$sort_row_data[int($temp_cross_row_number/2)-1];
	my $median=($m1+$m2)/2;
	
	my @MAD=();
	for (my $i=0; $i<$temp_cross_row_number; $i++) {	  
	    $MAD[$i]=abs($temp_row_data[$i]-$median);
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
  
} #end of subroutine



sub output_annotate_aggregate_cross_data {
    #output the data
    $aggregate_cross_data[0][$aggregate_cross_column_number+1]="Gene0";
    $aggregate_cross_data[0][$aggregate_cross_column_number+2]="Gene1";

    $aggregate_cross_data[0][$aggregate_cross_column_number+3]="Mate_Number";
    $aggregate_cross_data[0][$aggregate_cross_column_number+4]="anchor0_D";
    $aggregate_cross_data[0][$aggregate_cross_column_number+5]="anchor1_D";

    for (my $i=0; $i<=$aggregate_cross_row_number; $i++) {
	for (my $j=0; $j<=($aggregate_cross_column_number+4); $j++) {
	    print "$aggregate_cross_data[$i][$j]\t";
	}
	print "$aggregate_cross_data[$i][$aggregate_cross_column_number+5]\n";
    }
    
}

sub swap_length_with_orientation {
    #replace the read legnth with read orientation
    for (my $i=1; $i<=$cross_row_number; $i++) {
	$cross_data[$i][2]=$cross_data[$i][3];
	$cross_data[$i][7]=$cross_data[$i][8];
    }
}

sub output_annotate_aggregate_cross_data_exclude_unkn {
    #output the data
    my $min_cluster_size=2; #minimum number of abonormal pairing reads
    $aggregate_cross_data[0][$aggregate_cross_column_number+1]="Gene0";
    $aggregate_cross_data[0][$aggregate_cross_column_number+2]="Gene1";
    $aggregate_cross_data[0][$aggregate_cross_column_number+3]="Mate_Number";
    $aggregate_cross_data[0][$aggregate_cross_column_number+4]="anchor0_MAD";
    $aggregate_cross_data[0][$aggregate_cross_column_number+5]="anchor1_MAD";
    
    #re-aggregate the cross data
    @new_aggregate_data=();
    $new_aggregate_row_number=0;
    
    for (my $i=1; $i<=$aggregate_cross_row_number; $i++) {
	#the MAD has to been larger than 0 and at least two mates
	if (($aggregate_cross_data[$i][$aggregate_cross_column_number+3]>=$min_cluster_size) &&
	    (($aggregate_cross_data[$i][$aggregate_cross_column_number+4]>0) &&
	     ($aggregate_cross_data[$i][$aggregate_cross_column_number+5]>0))) {
	    #if ($aggregate_cross_data[$i][$aggregate_cross_column_number+3]>=$min_cluster_size) {
#	    my $temp_ratio1=$aggregate_cross_data[$i][$aggregate_cross_column_number+4]/
#		$aggregate_cross_data[$i][$aggregate_cross_column_number+5];
#	    my $temp_ratio2=$aggregate_cross_data[$i][$aggregate_cross_column_number+5]/
#		$aggregate_cross_data[$i][$aggregate_cross_column_number+4];
#	    
#	    my $temp_ratio=$temp_ratio1;
#	    if ($temp_ratio<$temp_ratio2) { $temp_ratio=$temp_ratio2; }
#	    
#	    #the ratio has to be smaller than 10 fold
#	    if ($temp_ratio<=20) {
	    
	    #if (!(($aggregate_cross_data[$i][$aggregate_cross_column_number+1]=~/Unkn/) ||
		#  ($aggregate_cross_data[$i][$aggregate_cross_column_number+2]=~/Unkn/))) {
	    $new_aggregate_row_number++;
	    for (my $j=0; $j<=($aggregate_cross_column_number+5); $j++) {
		$new_aggregate_data[$new_aggregate_row_number][$j]=$aggregate_cross_data[$i][$j];
	    }
	}
#	    }
#	}
    }
    #}

    #head
    for (my $j=0; $j<=($aggregate_cross_column_number+5); $j++) {
	$new_aggregate_data[0][$j]=$aggregate_cross_data[0][$j];
    }    
    
    #RE-ASSIGN
    $aggregate_cross_row_number=$new_aggregate_row_number;
    for (my $i=1; $i<=$aggregate_cross_row_number; $i++) {	
	for (my $j=0; $j<=($aggregate_cross_column_number+5); $j++) {
	    $aggregate_cross_data[$i][$j]=$new_aggregate_data[$i][$j];
	}    
    }
    
    #undef(@aggregate_cross_data);
    #@aggregate_cross_data=@new_aggregate_data;
    undef(@new_aggregate_data);

    #output
    #for (my $i=0; $i<=$aggregate_cross_row_number; $i++) {
	#if (!(($aggregate_cross_data[$i][$aggregate_cross_column_number+1]=~/Unkn/) ||
	 #     ($aggregate_cross_data[$i][$aggregate_cross_column_number+2]=~/Unkn/))) {
	  #  for (my $j=0; $j<=($aggregate_cross_column_number+4); $j++) {
		#print "$aggregate_cross_data[$i][$j]\t";
	    #}
	    #print "$aggregate_cross_data[$i][$aggregate_cross_column_number+5]\n";
	#}
    #}
    
}

sub build_hash_for_reduced_cross_data {
 
    #build indicator for putative fusions genes
    %putative_fusions_indicator=();
    for (my $i=1; $i<=$aggregate_cross_row_number; $i++) {
	my $temp_genes=$aggregate_cross_data[$i][$aggregate_cross_column_number+1]."\_".
	    $aggregate_cross_data[$i][$aggregate_cross_column_number+2];
	$putative_fusions_indicator{$temp_genes}=1;
    }
    
    
    #generate a new data array with putative fusions genes
    @reduced_cross_data=();
    $reduced_cross_row_number=0;
    $reduced_cross_column_number=$cross_column_number;

    $cross_data[0][$cross_column_number+1]="Gene0";
    $cross_data[0][$cross_column_number+2]="Gene1";
    #head
    for (my $j=0; $j<=($cross_column_number+2); $j++) {
	$reduced_cross_data[0][$j]=$cross_data[0][$j];
    }

    #content
    for (my $i=1; $i<=$cross_row_number; $i++) {
	my $temp_genes=$cross_data[$i][$cross_column_number+1]."\_".$cross_data[$i][$cross_column_number+2];
	if (exists($putative_fusions_indicator{$temp_genes})) {
	    $reduced_cross_row_number++;
	    for (my $j=0; $j<=($cross_column_number+2); $j++) {
		$reduced_cross_data[$reduced_cross_row_number][$j]=$cross_data[$i][$j];
	    }
	}
    }


    #build the alignment location indicator for putative fusions genes
    %fusions_location_indicator=();
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	my $temp_location_0=$reduced_cross_data[$i][0]."\_".$reduced_cross_data[$i][1];
	my $temp_location_1=$reduced_cross_data[$i][5]."\_".$reduced_cross_data[$i][6];
	$temp_location_0=~s/23\_/X\_/;
	$temp_location_0=~s/24\_/Y\_/;
	$temp_location_0=~s/25\_/M\_/;

	$temp_location_1=~s/23\_/X\_/;
	$temp_location_1=~s/24\_/Y\_/;
	$temp_location_1=~s/25\_/M\_/;

	$fusions_location_indicator{$temp_location_0}=1;
	$fusions_location_indicator{$temp_location_1}=1;
    }

}

sub load_unique_aligned_sorted_axt_file {
    
    my $axt_file=$ARGV[1];   
    open (INPUT, "<$axt_file")
	|| die "Can't open $axt_file $!";
    
    @location_readID=();
    $location_readID_row_number=0;
    $location_readID_indicator=0;

    %aligned_location_pointer1=();
    %aligned_location_pointer2=();

    my $temp_chr_line=0;    
    while (<INPUT>) {
	$text=$_;
	$text=~s/|\n//g;
	
	if ($text=~/chr/) {
	    $temp_chr_line=0;
	    $temp_chr_line++;
	    
	    (@line)=split(" ", $text);
	    my @axt_data=split(" ", $text);
	    my $temp_location=$axt_data[1]."\_".$axt_data[2];
	    $temp_location=~s/chr//;
	    if (exists($fusions_location_indicator{$temp_location})) {
		$location_readID_row_number++;
		$location_readID_indicator=1;
		$location_readID[$location_readID_row_number][0]=$temp_location; #read aligned location
		$location_readID[$location_readID_row_number][1]=$axt_data[7]; #read aligned orientation
		$location_readID[$location_readID_row_number][2]=$axt_data[4]; #readID
		if (!(exists($aligned_location_pointer1{$temp_location}))) {
		    $aligned_location_pointer1{$temp_location}=$location_readID_row_number;
		    $aligned_location_pointer2{$temp_location}=$location_readID_row_number;
		}
		else {
		    $aligned_location_pointer2{$temp_location}=$location_readID_row_number;
		}		    
	    }
	} #if ($text=~/chr/)
	else {
	    $temp_chr_line++;
	}
	
	#assign the sequence to readID
	if (($temp_chr_line==3) && ($location_readID_indicator==1)) {
	    $text=~s/\-//g;
	    $location_readID[$location_readID_row_number][3]=$text;
	    $location_readID_indicator=0;
	}
	
    } #while(INPUT)
    
    close(INPUT);

    #sort the location_readID
    sort_location_readID();
    generate_aligned_location_pointer_by_spanner_location();    
    
}

sub sort_location_readID {
    
    #sort the data by readID geneName
    @sort_readID_index=();
    for (my $i=1; $i<=$location_readID_row_number; $i++) {
	$sort_readID_index[$i][1]=$location_readID[$i][0];
	$sort_readID_index[$i][2]=$i;
    }
    

    quicksort_readID_by_spanner_location(1, $location_readID_row_number);
    
    #re-order the data by location
    my @temp_sort_readID=();
    for (my $i=1; $i<=$location_readID_row_number; $i++) {
	my $row_index=$sort_readID_index[$i][2];
	for (my $j=0; $j<=3; $j++) {
	    $temp_sort_readID[$i][$j]=$location_readID[$row_index][$j];
	}
    }

    #move data back
    for (my $i=1; $i<=$location_readID_row_number; $i++) {
	for (my $j=0; $j<=3; $j++) {
	    $location_readID[$i][$j]=$temp_sort_readID[$i][$j];
	}
    }
    
}

sub quicksort_readID_by_spanner_location {
    my ($left, $right)=@_;
    my ($pointer_left, $pointer_right);
    my ($pivot, $temp_1, $temp_2);

    if ($right>=$left) {
	$pivot=$sort_readID_index[$right][1];
	# print "pivot= $pivot\n";
	$pointer_left=$left-1;
	$pointer_right=$right;
	do
	{
	    do
	    {
		$pointer_left=$pointer_left+1;
	    } until ($sort_readID_index[$pointer_left][1] ge $pivot);
	    
	    do
	    { 
		$pointer_right=$pointer_right-1; 
	    } until (($pointer_right<=$left) || ($sort_readID_index[$pointer_right][1] le $pivot));
	    
	    $temp_1=$sort_readID_index[$pointer_left][1];
	    $sort_readID_index[$pointer_left][1]=$sort_readID_index[$pointer_right][1];
	    $sort_readID_index[$pointer_right][1]=$temp_1;
	    
	    $temp_2=$sort_readID_index[$pointer_left][2];
	    $sort_readID_index[$pointer_left][2]=$sort_readID_index[$pointer_right][2];
	    $sort_readID_index[$pointer_right][2]=$temp_2;
	    
	} until ($pointer_right<=$pointer_left);
	
	$sort_readID_index[$pointer_right][1]=$sort_readID_index[$pointer_left][1];
	$sort_readID_index[$pointer_left][1]=$sort_readID_index[$right][1];
	$sort_readID_index[$right][1]=$temp_1;
	
	$sort_readID_index[$pointer_right][2]=$sort_readID_index[$pointer_left][2];
	$sort_readID_index[$pointer_left][2]=$sort_readID_index[$right][2];
	$sort_readID_index[$right][2]=$temp_2;
	
	
	quicksort_readID_by_spanner_location($left, $pointer_left-1);
	quicksort_readID_by_spanner_location($pointer_left+1, $right);
    }
}

sub generate_aligned_location_pointer_by_spanner_location {
    %aligned_location_pointer1=();  
    %aligned_location_pointer2=();

    my $start_row_index=1;
    my $end_row_index=1;
    my $spanner_location=$location_readID[$start_row_index][0]; #start readID location
    
    my $moving_pointer=1;
    while ($moving_pointer<=$location_readID_row_number) {
	while (($moving_pointer<$location_readID_row_number) &&
	       ($location_readID[$moving_pointer+1][0] eq $spanner_location)) {
	    $moving_pointer++;
	}
	
	#start to assign blat the ordinal number
	$end_row_index=$moving_pointer;
	
	#save the spanner_location index
	$aligned_location_pointer1{$spanner_location}=$start_row_index;
	$aligned_location_pointer2{$spanner_location}=$end_row_index;
	
	#re-start the pointer;
	$moving_pointer++;
	$start_row_index=$moving_pointer;
	$end_row_index=$moving_pointer;
	
	if ($start_row_index<=$location_readID_row_number) {
	    $spanner_location=$location_readID[$start_row_index][0]; #reset the spanner_location
	}
    } #while ($moving_pointer<=$location_readID_row_number)
    
    
}


sub map_location_backto_readID {
    
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	my $temp_location_0=$reduced_cross_data[$i][0]."\_".$reduced_cross_data[$i][1];
	my $temp_location_1=$reduced_cross_data[$i][5]."\_".$reduced_cross_data[$i][6];		

	$temp_location_0=~s/23\_/X\_/;
	$temp_location_0=~s/24\_/Y\_/;
	$temp_location_0=~s/25\_/M\_/;

	$temp_location_1=~s/23\_/X\_/;
	$temp_location_1=~s/24\_/Y\_/;
	$temp_location_1=~s/25\_/M\_/;
	
	if ((exists($aligned_location_pointer1{$temp_location_0})) &&
	    (exists($aligned_location_pointer2{$temp_location_0})) &&
	    (exists($aligned_location_pointer1{$temp_location_1})) &&
	    (exists($aligned_location_pointer2{$temp_location_1}))) {
	    

	    $anchor0_pointer_1=$aligned_location_pointer1{$temp_location_0};
	    $anchor0_pointer_2=$aligned_location_pointer2{$temp_location_0};
	    
	    
	    $anchor1_pointer_1=$aligned_location_pointer1{$temp_location_1};
	    $anchor1_pointer_2=$aligned_location_pointer2{$temp_location_1};
	    
	    
	    for (my $h=$anchor0_pointer_1; $h<=$anchor0_pointer_2; $h++) {
		my @temp_anchor0_readID=split("\/", $location_readID[$h][2]);
		my $last_indicator=0;
		for (my $k=$anchor1_pointer_1; $k<=$anchor1_pointer_2; $k++) {
		    my @temp_anchor1_readID=split("\/", $location_readID[$k][2]);
		    
		    if ($temp_anchor0_readID[0] eq $temp_anchor1_readID[0]) {
			#strand, readID, and sequence
			$reduced_cross_data[$i][2]=$location_readID[$h][1];
			$reduced_cross_data[$i][3]="\>".$location_readID[$h][2];
			#reverse the complementary sequence
			if ($location_readID[$h][1]=~/\+/) {
			    $reduced_cross_data[$i][4]=$location_readID[$h][3];
			}
			else {
			    $reduced_cross_data[$i][4]=
				generate_complementary_sequence($location_readID[$h][3]);
			}

			#$reduced_cross_data[$i][4]=$location_readID[$h][3];

			#strand, readID, and sequence			
			$reduced_cross_data[$i][7]=$location_readID[$k][1];
			$reduced_cross_data[$i][8]="\>".$location_readID[$k][2];
			#reverse the complementary sequence
			if ($location_readID[$k][1]=~/\+/) {
			    $reduced_cross_data[$i][9]=$location_readID[$k][3];
			}
			else {
			    $reduced_cross_data[$i][9]=
				generate_complementary_sequence($location_readID[$k][3]);
			}

			#$reduced_cross_data[$i][9]=$location_readID[$k][3];
			$last_indicator=1;
			last;
			#$last_indicator=1;
		    }
		    
		}
		
		if ($last_indicator==1) {
		    last;
		}
	    }
	}
	
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
    #my $axt_file=$ARGV[1];   
    my $temp_fa_file=$ARGV[2];
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
    my $temp_fa_file=$ARGV[2];
    #my $minScore=$ARGV[5]; #=30; #this one is for Blat search sensitivity
    #$temp_fa_file=~s/Chip/Blat/;
    $temp_fa_file=~s/\.txt/\.fa/;
    my $temp_psl_file=$temp_fa_file;
    $temp_psl_file=~s/\.fa/\.psl/;
    #system("blat ../../RefSeq/hgGenome/hgGenome.2bit temp_blat.fa temp_blat.psl");
    my $genome_version="hg19";
    my $data_path_file="\/RIS\/home\/xsu1\/Mosaik\/hg19_Mosaik\/".$genome_version."\_2bit"."\/".$genome_version."OV\.2bit";
    #system("blat /home2/xsu/Mosaik/RefSeq/hgGenome/hgGenome.2bit $temp_fa_file $temp_psl_file");
    system("/RIS/home/xsu1/bin/blat $data_path_file $temp_fa_file $temp_psl_file");    
}

sub get_blat_data {
    #my $blat_file="temp_blat.psl";
    my $temp_fa_file=$ARGV[2];
    #$temp_fa_file=~s/Chip/Blat/;
    $temp_fa_file=~s/\.txt/\.fa/;
    my $blat_file=$temp_fa_file;
    $blat_file=~s/\.fa/\.psl/;
    #my $blat_file="temp_blat.psl";
    open (INPUT, "<$blat_file")
	|| die "Can't open $blat_file $!";

    my $start_line=0;
    my $refGene_RNA_Genome=$ARGV[0]; #=hg19_refGene_RIS.txt or hg19_TruecDNA_refGene_RIS.txt
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
    
    unlink($temp_fa_file);
    unlink($blat_file);
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

sub filter_blat_reduced_cross_data_for_chip {

    #chip data indicator
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	$reduced_cross_data[$i][$reduced_cross_column_number+6]="No";
    }
    
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
	
	my $temp_BlatHit_yes=0;
	my $temp_BlatHit_no=0;
	my $temp_BlatHit_dup=0;
	for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
	    if ($reduced_cross_data[$h][$reduced_cross_column_number+5]=~/Yes/) {
		$temp_BlatHit_yes++;
	    }
	    else {
		$temp_BlatHit_no++;
	    }
	    if ($reduced_cross_data[$h][$reduced_cross_column_number+5]=~/Duplicate/) {
		$temp_BlatHit_dup++;
	    }
	}
	
	
	if (($temp_BlatHit_yes>($temp_BlatHit_no/2)) && ($temp_BlatHit_dup<=2)) {
	    for (my $h=$temp_start_pointer; $h<=$temp_end_pointer; $h++) {
		if ((($reduced_cross_data[$h][2]=~/\+/) || ($reduced_cross_data[$h][2]=~/\-/)) &&
		    (($reduced_cross_data[$h][7]=~/\+/) || ($reduced_cross_data[$h][7]=~/\-/))) {
		    $reduced_cross_data[$h][$reduced_cross_column_number+6]="Chip";
		}
	    }
	}
	
	#re-start the pointer;
	$moving_pointer++;
	$temp_start_pointer=$moving_pointer;
	$temp_end_pointer=$moving_pointer;
	$gene_ID_0=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+1];
	$gene_ID_1=$reduced_cross_data[$temp_start_pointer][$reduced_cross_column_number+2];
	
    }
    
}

sub output_reduced_annotate_cross_data {
    #output the data

    $temp_chip_file=$ARGV[2];
    open (OUTPUT, ">$temp_chip_file")
	|| die "Can't open $temp_chip_file $!";

    for (my $i=0; $i<=0; $i++) {
	for (my $j=0; $j<=($reduced_cross_column_number+1); $j++) {
	    print OUTPUT "$reduced_cross_data[$i][$j]\t";
	}
	print OUTPUT "$reduced_cross_data[$i][$reduced_cross_column_number+2]\n";
    }
    #content
    for (my $i=1; $i<=$reduced_cross_row_number; $i++) {
	if ($reduced_cross_data[$i][$reduced_cross_column_number+6]=~/Chip/) {
	    for (my $j=0; $j<=($reduced_cross_column_number+1); $j++) {
		print OUTPUT "$reduced_cross_data[$i][$j]\t";
	    }
	    print OUTPUT "$reduced_cross_data[$i][$reduced_cross_column_number+2]\n";
	}
    }
    close(OUTPUT);
}


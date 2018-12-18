#!/usr/bin/perl
##########################################################################
###  snpfx.pl                                                          ###
##########################################################################
###  DETERMINES THE EFFECT OF A SNP ON A REFERENCE SEQUENCE            ###
##########################################################################
##########################################################################


#LOAD BIOPERL MODULES
# Bioperl SeqIO module

use Bio::Perl;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::RefSeq;


##READ SNP FILE NAME FROM COMMANDLINE ARGUMENTS
$snp_file=$ARGV[0];

##READ REFERENCE SEQUENCE FROM COMMAND LINE ARGUMENTS
$filename=$ARGV[1];


#$gb= new Bio::DB::RefSeq;


########################################################################################
###########  LOAD NCBI REFERENCE SEQUENCE     ##########################################
########################################################################################
##### IF SEQUENCE ALREADY EXISTS LOCALLY, THEN WE DON'T NEED TO PULL FROM NCBI #########
##### SO THIS FUNCTION FIRST CHECKS FOR THE FILE BY SEARCHING LOCAL DIR FOR A  #########
##### FOR A FILE WITH THAT NAME.  IF IT EXISTS, WE LOAD THAT.  IF IT DOESN'T   #########
##### THEN WE USE THE BIOPERL REFSEQ MODULE TO PULL FROM THE NCBI              #########
########################################################################################



if (-e $filename) {
   $input_sequence = Bio::SeqIO->new( -file => "<$filename");
	#print "Found Filename\n";
  }
else {
  die "SNP effect infers the consequence of SNPs within coding genes on the gene product.\n
       It takes 2 inputs, 1\) a list of SNPs relative to a reference genome sequence\n
                             format \= position \(n-1\) \"tab\" sample=nt \"tab\" ref=nt\n
                          2\) an annotated reference genome in genbank .gbk format\n
       Usage: snpfx.pl <snp-list.txt> <reference.gbk> >name-outfile.txt\n";
  };

##READ SEQ OBJECT FROM STREAM##
my $seq_object = $input_sequence->next_seq;

##########################################################################################
######### EXTRACT GENE FEATURES FROM DOWNLOADED SEQUENCE  ################################
##########################################################################################

for my $gene_features ($seq_object->get_SeqFeatures) {
  #PULL ONLY THE GENE TAGS FROM FEATURE LIST AND STORE IN ARRAY
    if ($gene_features->primary_tag eq "CDS") {

       #Make sure feature has locus_tag ID
       if ($gene_features->has_tag('locus_tag')) {
          for $locusID ($gene_features->get_tag_values('locus_tag')) {
             #OUTPUT LOCUS ID TAG
             #print "locus ID:  ", $locusID, "\t";
              $feature_array[$loop_counter][1]=$locusID;
           };

if ($gene_features->has_tag('gene')) {
          for $gene_name ($gene_features->get_tag_values('gene')){
               $feature_array[$loop_counter][5]=$gene_name;
            }; }
      else {
	     $feature_array[$loop_counter][5]="--";
           };

    ###RETRIEVE ORF START LOCATION USING BIOPERL LOCATION FUNCTION
       my $gene_start = $gene_features->location->start;
    ###RETRIEVE ORF END LOCATION USING BIOPERL LOCATION FUNCTION
       my $gene_end = $gene_features->location->end;

    ###FIND STRAND ORIENTATION FOR EACH ORF
       my $strand_feature = $gene_features->location->strand;
           #  print "$gene_start  $gene_end  $strand_feature \n";

    ###LOAD FEATURES INTO Nx4 ARRAY
       $feature_array[$loop_counter][2]=$gene_start;
       $feature_array[$loop_counter][3]=$gene_end;
       $feature_array[$loop_counter][4]=$strand_feature;
       $loop_counter++;

          };
            };
};


###OPEN INPUT SNP FILE
open (SNPFILE, $snp_file) or die "Can't find SNP file!";

###########LOAD SNPFILE CONTENTS INTO ARRAY################
foreach $line (<SNPFILE>) {
	@row=(split(/\t/, $line));
        push @SNParray, [@row];
};

#TEST: PRINT OUT ARRAY
#foreach $line(@SNParray) {	
#       foreach $item(@$line){
#	  print "$item\t";
#	 };
#        print "\n";
#};
print "Position\(n-1\)\tGeneID\tName\tType\tPos_in_Gene\tINDEL_Length\tRef_base\tSample_base\tRef_Codon\tSample_Codon\tRef_AA\tSample_AA\tSubstitution\n";
#ITERATE THROUGH FEATURE ARRAY AND LOOK TO SEE IF SNP POSITION IS BETWEEN $GENE_START AND $GENE_END
	foreach $line(@SNParray) {
		$codon="";
		$effect="";
		$amino_acid="";
		$ref_amino_acid="";
		@row=@$line;
		$found_gene=0;
		#$item=$row[0];
		#print "$item\n";
			 $i=0;
			 
			 ##EXTRACT VARIANT SEQUENCE
				  foreach $item(@row){
					chomp $item;
					if($item=~ m/sample/){
						$index=index($item, "=");
						$sample_base= substr($item,($index+1));
					};
					if($item=~ m/ref/){
						$index=index($item, "=");
						$ref_base= substr($item,($index+1));
						$ref_len=(length $ref_base);
					};
				   };
				  
			foreach $feature (@feature_array) {
			   ###TEST IF SNP IS BETWEEN START AND END
			   if((($row[0])+1)>=$feature_array[$i][2] && $row[0]<=$feature_array[$i][3]) {
			       #WRITE SNP POSITION AND NAME OF GENE IT'S INSIDE TO OUTPUT FILE
			      #open (OUTFILE, ">>$outfile");
			       #print $feature_array[$i][1] \t $feature_array[$i][5]; 
				  
				  ##WE ONLY WANT SNPs SO EXCLUDING ANY VARIANT LARGER THAN 1 BP IN LENGTH
				  $len=(length $sample_base);
				  if(($len==1) && ($ref_len==1)){
					  
					###DETERMINE NEW CODON IF GENE IS ON FORWARD STRAND###
					if(($feature_array[$i][4])==1){
						####FIND WHERE VARIANT IS LOCATED IN THE GENE
						 $position_in_gene=(($row[0]+1)-$feature_array[$i][2]);
						   ####FIND WHICH CODON POSITION VARIANT IS LOCATED
						   $codon_position=($position_in_gene % 3);
						if($codon_position==0) {
							 $codon_seq=$seq_object->subseq(($row[0]+1),($row[0]+3));
							  $mut_seq=$sample_base.($seq_object->subseq(($row[0]+2),($row[0]+3)));
							  $codon_obj=Bio::Seq->new(-seq =>$mut_seq, -alphabet => 'dna' );
							  $codon=$codon_obj->seq;
						   };
						  if($codon_position==1) {
							  $codon_seq=$seq_object->subseq(($row[0]),($row[0]+2));
							  $mut_seq=($seq_object->subseq(($row[0]),($row[0]))).($sample_base).($seq_object->subseq(($row[0]+2),($row[0]+2)));
							  $codon_obj=Bio::Seq->new(-seq =>$mut_seq, -alphabet => 'dna' );
							  $codon=$codon_obj->seq;
						  };
						  if($codon_position==2) {
							  $codon_seq=$seq_object->subseq(($row[0]-1),($row[0]+1));
							  $mut_seq=($seq_object->subseq(($row[0]-1),($row[0]))).($sample_base);
							  $codon_obj=Bio::Seq->new(-seq =>$mut_seq, -alphabet => 'dna' );
							  $codon=$codon_obj->seq;
						  };
					};
				####ETERMINE NEW CODON IF GENE IS ON REVERSE STRAND###
					if(($feature_array[$i][4])==-1){
						####FIND WHERE VARIANT IS LOCATED IN THE GENE
						 $position_in_gene=(($feature_array[$i][3])-($row[0]+1));
						   ####FIND WHICH CODON POSITION VARIANT IS LOCATED
						   $codon_position=($position_in_gene % 3);
						 if($codon_position==0) {
							$codon_seq=$seq_object->subseq(($row[0]-1),($row[0]+1));
							  $mut_seq=($seq_object->subseq(($row[0]-1),($row[0]))).($sample_base);
							  $mut_seq_obj=Bio::Seq->new(-seq=>$mut_seq, -alphabet => 'dna' );
							  $codon=$mut_seq_obj->revcom->seq;
							  $ref_codon_obj=Bio::Seq->new(-seq =>$codon_seq, -alphabet => 'dna' );
							  $codon_seq=$ref_codon_obj->revcom->seq
						  };
						 if($codon_position==1) {
							  $codon_seq=$seq_object->subseq(($row[0]),($row[0]+2));
							  $mut_seq=($seq_object->subseq(($row[0]),($row[0]))).($sample_base).($seq_object->subseq(($row[0]+2),($row[0]+2)));
							  $mut_seq_obj=Bio::Seq->new(-seq=>$mut_seq, -alphabet => 'dna' );
							  $codon=$mut_seq_obj->revcom->seq;
							  $ref_codon_obj=Bio::Seq->new(-seq =>$codon_seq, -alphabet => 'dna' );
							  $codon_seq=$ref_codon_obj->revcom->seq
							  
						   };
						  if($codon_position==2) {
							  $codon_seq=$seq_object->subseq(($row[0]+1),($row[0]+3));
							  $mut_seq=$sample_base.($seq_object->subseq(($row[0]+2),($row[0]+3)));
							  $mut_seq_obj=Bio::Seq->new(-seq=>$mut_seq, -alphabet => 'dna' );
							  $codon=$mut_seq_obj->revcom->seq;
							  $ref_codon_obj=Bio::Seq->new(-seq =>$codon_seq, -alphabet => 'dna' );
							  $codon_seq=$ref_codon_obj->revcom->seq
							  
						  };
					};
						###TRANSLATE CODON TO AMINO ACID CODES
						$amino_acid=translate_as_string($codon);
						$ref_amino_acid=translate_as_string($codon_seq);
						if ($amino_acid eq $ref_amino_acid) {
							$effect="synonymous";
						
						}
						else {
							$effect="nonsynonymous";
						};
						
						print "$row[0]\t $feature_array[$i][1] \t $feature_array[$i][5]\t SNP \t $position_in_gene \t --- \t REF=$ref_base \t SAMPLE=$sample_base \t $codon_seq \t $codon \t";
						print "$ref_amino_acid \t $amino_acid \t $effect\n";
						
						}
					else {   
						if(($feature_array[$i][4])==1){
						####FIND WHERE VARIANT IS LOCATED IN THE GENE
						 $position_in_gene=(($row[0]+1)-$feature_array[$i][2]);
						};
						if(($feature_array[$i][4])==-1){
						####FIND WHERE VARIANT IS LOCATED IN THE GENE
						 $position_in_gene=(($feature_array[$i][3])-($row[0]+1));
						};
						$indel_length= ($len) - ($ref_len);
						
						print "$row[0]\t $feature_array[$i][1] \t $feature_array[$i][5]\t INDEL \t $position_in_gene \t $indel_length \t REF=$ref_base \t SAMPLE=$sample_base\n";
						
					};
					 $found_gene=1;
					 };
					$i++;
				};
				if ($found_gene==0){
				print "$row[0]\t IG \n";
			 };	
	};


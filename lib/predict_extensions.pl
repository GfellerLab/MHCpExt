use List::MoreUtils qw(uniq);

##################################################################
# Predict cases of C-terminal extensions from HLA peptodomics data
# Model 1:
# For each 9-mer motif, generate $lg-mer ($lg=10,11,12) peptides based on the 9-mer + bulge model derived from HLA peptidomics 9-mer data. Generate them in the same proportions as the 9-mer MS data.
# Make predictions with the three models (buldge, Cterm extension, Nterm extension) for all motifs.
# Compare the number of MS $lg-mers predicted to form an extension based on the selected motif with the number that you expdect by chance.
# Renormalization is done by the number of $lg-mers predicted to bind to this motif, either bulge, C-term or N-term (this is to consider the fact that differences in peptide presentation may be observed between 9- and $lg-mers)
#
# One limitation is if the contribution of each motif changes and if some data from another motif (with more ligands than in the 9-mers) look like extension of the motif under consideration.
# A second limitation is if the extensions have some specificity that is similar to the anchor residue (impossible to study with sequence-based approaches)
# A third limitation are cases of extension + buldge (for 11+ mers)
# A forth limitation is due to very annecdotal cases that do not pass the Z-score thresholds
###################################################################


%map=("A", 0, "C", 1, "D", 2, "E", 3, "F", 4, "G", 5, "H", 6, "I", 7, "K", 8, "L", 9, "M", 10, "N", 11, "P", 12, "Q", 13, "R", 14, "S", 15, "T", 16, "V", 17, "W", 18, "Y", 19);
@letter=qw(A C D E F G H I K L M N P Q R S T V W Y);
$N=20;
$tl=9;
$pseudo_count=200;

$input_file=$ARGV[0];
$output_dir=$ARGV[1];

$ncluster=$ARGV[2];
$lg=$ARGV[3];			#Length of the extensions

$lib_path=$ARGV[4];
    
$cluster_name=$ARGV[5];

if ($cluster_name ne "") {
    @cluster=split(",", $cluster_name);
} else {
    for ($j=1; $j<=$ncluster; $j++) {
	$cluster[$j-1]="Cluster_$j";
    }
}

print "\n#########\n$lg-mers\n";

$rep=100;	    #Number of repetition for the random model
$Nr=100000;	    #Number of random peptides from the human proteome
$Zs=2;		    #Z-score for the difference with the random model

$Thresh_resp=0;	#Threshold on the responsibilities to build the 9-mer PWMs
$Thresh_score=2; #Threshold on the scores of the models (P1, P2, P3, + 2 C-terminal). This is to select random peptides and to decide which of the MS peptides are following a given model.
$fr_top=0.02;

###########
# Take the human proteome
###########


open IN, "$lib_path/../data/Homo_sapiens_all_RefSeq.txt", or die;
$l=<IN>;
$nprot=0;
while ($l=<IN>) {
    $l=<IN>;
    $l =~ s/\r?\n$//;
    chomp($l);
    push @human_seq, $l;
    $len[$nprot]=length($l);
    $nprot++;
}
close IN;


#################
# Build a list of random peptides
#################

@rseq=();
srand(101);
for ($r=0; $r<$rep; $r++) {
    for ($n=0; $n<$Nr; $n++) {
	$r1=int(rand($nprot));
	$r2=int(rand($len[$r1]-$lg));
	push @{$rseq[$r]}, substr($human_seq[$r1],$r2, $lg); 
    }
    @{$rseq[$r]}=uniq(@{$rseq[$r]});
}


#######################################
# Take the human amino acid frequencies
#######################################

open IN, "$lib_path/../data/humanProteome.txt", or die;
$tot=0;
while ($l=<IN>) {
    $l =~ s/\r?\n$//;
    chomp($l);
    @a=split(' ', $l);
    $fr[$map{$a[0]}]=$a[1];
    $tot=$tot+$a[1];
}
close IN;
for ($i=0; $i<$N; $i++) {
    $fr[$i]=$fr[$i]/$tot;
}


########################
# Load the Blosum Matrix
########################

open IN, "$lib_path/../data/blosum62.txt", or die;
$l=<IN>;
$j=0;
while ($l=<IN>) {
    $l =~ s/\r?\n$//;
    chomp($l);
    @a=split(' ', $l);
    for ($i=1; $i<scalar @a; $i++) {
	$blo[$j][$i-1]=$a[$i];
    }
    $j++;
}
close IN;

############
# Create output directories
############

if ( -d "$output_dir/Cterm/class1_$lg") {
    system("rm -r $output_dir/Cterm/class1_$lg");
}
system("mkdir $output_dir/Cterm/class1_$lg");

if ( -d "$output_dir/Nterm/class1_$lg") {
    system("rm -r $output_dir/Nterm/class1_$lg");
}
system("mkdir $output_dir/Nterm/class1_$lg");

if ( -d "$output_dir/buldge/class1_$lg") {
    system("rm -r $output_dir/buldge/class1_$lg");
}
system("mkdir $output_dir/buldge/class1_$lg");



#############################################
# For each sample, go through longer peptides
#############################################

open OUTC, ">$output_dir/Cterm/Cterm_extensions_$lg\mer.txt";
print OUTC "Cluster\t$lg-mers\tCterm_Extension\tExpected\tSD\tTot_pep_ass\tFraction\tZ-score\n";

open OUTN, ">$output_dir/Nterm/Nterm_extensions_$lg\mer.txt";
print OUTN "Cluster\t$lg-mers\tNterm_Extension\tExpected\tSD\tTot_pep_ass\tFraction\tZ-score\n";

   
#################
# Load the 9-mers and their respective cluster assignment by MixMHCp
#################

$c1=0;
$Npep9=0;

@pep9=([]);


open IN, "$output_dir/MixMHCp/responsibility/resp_$ncluster.txt", or die;
$l=<IN>;
@a=split(' ', $l);
if($a[$ncluster+1] eq "Trash"){
    $trash=1;
} else {
    $trash=0;
}
while ($l=<IN>) {
    $l =~ s/\r?\n$//;
    chomp($l);
    @a=split(' ', $l);
    $max=0;
    for ($cl=0; $cl<$ncluster+$trash; $cl++) { 
	if ($a[$cl+1]>$max) {
	    $cluster=$cl;
	    $max=$a[$cl+1];
	}
    }
    if ($max>$Thresh_resp && $cluster<$ncluster) { # Do not include peptides from the trash
	push @{$pep9[$cluster]}, $a[0]; 
	$c1++;
    }
    $Npep9++;
}
close IN;

################
# Build the PWMs for each motif, including the blosum correction
################

@pwm=([[]]);

for ($cl=0; $cl<$ncluster; $cl++) {
    
    $le=scalar @{$pep9[$cl]};

    for ($s=0; $s<$tl; $s++) {
	for ($j=0; $j<$N; $j++) {
	    $pwm[$cl][$j][$s]=0;
	}
    }
	
    foreach $ps (@{$pep9[$cl]}) {
	@a=split('', $ps);
	if (scalar @a != $tl) {
	    print "Problem: $tl\t$ps\t$cl\n";
	    die;
	}
	for ($s=0; $s<$tl; $s++) {
	    $pwm[$cl][$map{$a[$s]}][$s]++;
	}
    }
	
    for ($s=0; $s<$tl; $s++) {
	for ($j=0; $j<$N; $j++) {
	    $g[$s][$j]=0;
	    #Correct the raw counts by the background freq
	    for ($m=0; $m<$N; $m++) {
		$g[$s][$j]=$g[$s][$j]+$blo[$m][$j]*($pwm[$cl][$m][$s]/$le);
	    }
	}
    }
    for ($s=0; $s<$tl; $s++) {
	for ($j=0; $j<$N; $j++) {
	    $pwm[$cl][$j][$s]=($pwm[$cl][$j][$s]+$pseudo_count*$g[$s][$j])/($le+$pseudo_count);
	}
    }
}



##############################
# Find the optimal threshold for each cluster
##############################

@Thresh=();
for ($j=0; $j<$ncluster; $j++) {

    @score_distr=();
    
    for ($n=0; $n<$Nr; $n++) {
	@seq=split('', $rseq[0][$n]);
	$tscore=0;
	
	for ($p=0; $p<3; $p++) {
	    $s=$map{$seq[$p]};
	    $tscore=$tscore+log($pwm[$j][$s][$p]/$fr[$s]); 
	}
	for ($p=$tl-2; $p<$tl; $p++) {
	    $s=$map{$seq[$p-$tl+$lg]};
	    $tscore=$tscore+log($pwm[$j][$s][$p]/$fr[$s]); 
	}
	push @score_distr, $tscore;
    }

    @score_distr=sort {$b <=> $a} @score_distr;
    $Thresh[$j]=$score_distr[int($Nr*$fr_top)];
    $Thresh[$j]=2;
	
}


###############################
# Select random sets of peptides.
# by running the bulge model only.
###############################

$Npep10=0;
    
open IN, "$input_file", or die;
 
while ($l=<IN>) {
    $Npep10++;
}
close IN;

print "9-mers: $Npep9\t$lg-mers: $Npep10\n";


@rseq_p=([[]]);	
for ($r=0; $r<$rep; $r++) {

    #Go through all random peptides.

    #From all the peptides that passed the threshold for the different bulge motifs built from the 9-mer,
    #select a number corresponding to the expected number if motifs are sampled with the frequency observed in the 9-mers.

    for ($j=0; $j<$ncluster; $j++) {

	@{$rseq_p[$r][$j]}=();
	    
	$n9=scalar @{$pep9[$j]};
	$n10=int( $n9 / $Npep9 * $Npep10 );
	$ct10=0;
	for ($n=0; $n<$Nr && $ct10<$n10; $n++) {
	    @seq=split('', $rseq[$r][$n]);
	    $tscore=0;
		
	    for ($p=0; $p<3; $p++) {
		$s=$map{$seq[$p]};
		$tscore=$tscore+log($pwm[$j][$s][$p]/$fr[$s]); 
	    }
	    for ($p=$tl-2; $p<$tl; $p++) {
		$s=$map{$seq[$p-$tl+$lg]};
		$tscore=$tscore+log($pwm[$j][$s][$p]/$fr[$s]); 
	    }
	    if ($tscore>$Thresh[$j]) {
		push @{$rseq_p[$r][$j]}, $rseq[$r][$n];
		$ct10++;
	    }
	}
    }
}

    

######################
# Make the predictions
######################
	
for ($j=0; $j<$ncluster; $j++) {

    $tj=$j+1;
    
    $map=0;
 	
    # Go through the actual $lg-mer peptides
    # Find those that are better predicted with the 9-mer logo + C-terminal extensions.
	    
    $Cct=0;
    $Nct=0;
    $bct=0;
    $ct=0;
    @Cpep=();
    @Npep=();
    @bpep=();
    $Tot_pep_ass=0;


    @non_bpep=();
    
    open IN, "$input_file", or die;
    while ($l=<IN>) {
	$l =~ s/\r?\n$//;
	chomp($l);
	$l = uc $l;
		
	$ct++;
		
	@res=check_extensions($l, $j); # 1 for C-terminal ext; 1 for N-terminal ext; 1 for bulges; score for C-term; Score for N-term; Score for bulges; 1 for unambiguous allele assignment 

	$max=-1000;
	for ($ii=3; $ii<=5; $ii++) {
	    if ($max<$res[$ii]) {
		$max=$res[$ii];
	    }
	}
	$check=0;

	if ($res[5] < $Thresh[$j]) {
	    push @non_bpep, $l;
	}
	
	#Count cases where the peptide is better described by the Cterm extension
	if ( $res[3] > $Thresh[$j] && $res[0]==1) {
	    $Cct++;
	    push @Cpep, $l;
	    $check=1;
	}
	#Count cases where the peptide is better described by the Nterm extension
	if ($res[4] > $Thresh[$j] && $res[1]==1) {
	    $Nct++;
	    push @Npep, $l;
	    $check=1;
	}
	#Count cases where the peptide is better described by the bulge
	if ($res[5] > $Thresh[$j] && $res[2]==1) {
	    $bct++;
	    push @bpep, $l;
	    $check=1;
	}
	#Count how many peptides are assigned to this allele (either bulge, N- or C-term)
	if ( ($res[3] > $Thresh[$j] ||  $res[4] > $Thresh[$j] || $res[5] > $Thresh[$j]) && $res[6]==1 ) {
	    $Tot_pep_ass++;
	}

	#####
	# Note that longer peptides that would be in the trash of MixMHCp2.0 are automatically excluded, because they score with any model is smaller than 0
	#####
	
    }
    close IN;
    
    print "\n$lg-mers assigned to $cluster[$j]: $Tot_pep_ass\n";

    #die;
    
    #Decide which random model to use depending on evidences of buldge
    if ($Tot_pep_ass >= 5 && $Tot_pep_ass > 0.1*scalar @{$pep9[$j]}/$Npep9*$Npep10 ) {
	#If the bulge is well supported and at least 10% of $lg-mer follow it for this allele. This is useful since some alleles are basically absent in $lg-mers, like HLA-C alleles.
	$model=1;
    } else {
	$model=2;
    }
    
    
    
    @Npep10_sub=();
    for ($r=0; $r<$rep; $r++) {
	
	#If many sequences match the 9-mer + bulge model, consider random sequence generated based on all 9-mer motifs + buldge (test Model 1: The $lg-mers are fully explained by the bulge model).
	if ( $model==1 ) {
	    $rCct_all[$r]=0;
	    $rNct_all[$r]=0;
	    $rbct_all[$r]=0;
	    $rTot_pep_ass[$r]=0;
		    
	    @rseq_all=();
	    for ($j2=0; $j2<$ncluster; $j2++) {
		push @rseq_all, @{$rseq_p[$r][$j2]};
			
	    }
	    @rseq_all=uniq(@rseq_all);
		    
	    foreach $l (@rseq_all) {

		@res=check_extensions($l, $j);
		#Count cases where the peptide is better described by the Cterm extension

		
		if ( $res[3] > $Thresh[$j] && $res[0]==1) {
		    $rCct_all[$r]++;
		}
		
		#Count cases where the peptide is better described by the Nterm extension
		if ($res[4] > $Thresh[$j] && $res[1]==1) {
		    $rNct_all[$r]++;
		}
		#Count cases where the peptide is better described by the bulge
		if ($res[5] > $Thresh[$j] && $res[2]==1) {
		    $rbct_all[$r]++;
		}
		if ( ($res[3] > $Thresh[$j] ||  $res[4] > $Thresh[$j] || $res[5] > $Thresh[$j]) && $res[6]==1 ) {
		    $rTot_pep_ass[$r]++;
		}
		
	    }
	    #Do the normalization so as to have the same number assigned to each allele as in the 10-mer HLA peptidome
	    $rCct_all[$r]=$rCct_all[$r]/$rTot_pep_ass[$r]*$Tot_pep_ass;
	    $rNct_all[$r]=$rNct_all[$r]/$rTot_pep_ass[$r]*$Tot_pep_ass;
	    $rbct_all[$r]=$rbct_all[$r]/$rTot_pep_ass[$r]*$Tot_pep_ass;
	  
	} 
    }

    if ($model==1) {
	
	printf OUTC "$cluster[$j]\t$ct\t$Cct\t";
	if (mean(@rbct_all) > 0 ) {
	    printf OUTC "%.2f\t%.2f", mean(@rCct_all), std(@rCct_all);
	} else {
	    printf OUTC "0\t0";
	}
	if (mean(@rCct_all)>0) {
	    $Z=($Cct-mean(@rCct_all))/std(@rCct_all);
	} else {
	    $Z=100;
	}
	printf OUTC "\t$Tot_pep_ass\t%.3f\t%.3f\n", $Cct/$Tot_pep_ass, $Z;

	printf "C-terminal extensions assigned to $cluster[$j]: $Cct\n";
	printf "Expected number: %.3f (Z-score: %.3f)\n", mean(@rCct_all), $Z;
	

	printf OUTN "$cluster[$j]\t$ct\t$Nct\t";
	if (mean(@rbct_all) > 0 ) {
	    printf OUTN "%.2f\t%.2f", mean(@rNct_all), std(@rNct_all);
	} else {
	    printf OUTN "0\t0";
	}
	if (mean(@rNct_all)>0) {
	    $Z=($Nct-mean(@rNct_all))/std(@rNct_all);
	} else {
	    $Z=100;
	}

	printf "N-terminal extensions assigned to $cluster[$j]: $Nct\n";
	printf "Expected number: %.3f (Z-score: %.3f)\n", mean(@rNct_all), $Z;

	
	printf OUTN "\t$Tot_pep_ass\t%.3f\t%.3f\n", $Nct/$Tot_pep_ass, $Z;
    }
    if ($model==2) {
	printf OUTC "$cluster[$j]\t$ct\tNA\tNA\tNA\tNA\tNA\tNA\n";
	printf OUTN "$cluster[$j]\t$ct\tNA\tNA\tNA\tNA\tNA\tNA\n";
    }

    
    open OUT, ">$output_dir/Cterm/class1_$lg/cluster_$tj.fa";
    foreach $p (@Cpep) {
	print OUT ">pep\n$p\n";
    }
    close OUT;

    open OUT, ">$output_dir/Nterm/class1_$lg/cluster_$tj.fa";
    foreach $p (@Npep) {
	print OUT ">pep\n$p\n";
    }
    close OUT;
	    
    open OUT, ">$output_dir/buldge/class1_$lg/cluster_$tj.fa";
    foreach $p (@bpep) {
	print OUT ">pep\n$p\n";
    }
    close OUT;
}

close OUTC;
close OUTN;

sub mean(){

    my $m=0;
    my $ct=0;
    foreach my $i (@_) {
	$m=$m+$i;
	$ct++;
    }
    if ($ct>0) {
	$m=$m/$ct;
    } else {
	$m=0;
    }
    
    return($m);
  
}

sub std(){

    my $m=0;
    my $ct=0;
    foreach my $i (@_) {
	$m=$m+$i;
	$ct++;
    }
    if ($ct>0) {
	$m=$m/$ct;
    } else {
	$m=0;
    }
    my $sd=0;
    foreach my $i (@_) {
	$sd=$sd+($i-$m)*($i-$m);
    }
    if ($ct>0) {
	$sd=sqrt($sd/$ct);
    } else {
	$ct=0;
    }
    return($sd);
  
}

sub other_gap{
    #Other possibilities of dealing with gaps
    if ($fixed_gap==0) {
	for ($gap=3; $gap<$lg-2; $gap++) { # Gaps at different positions
	    $tscore=0;
	    for ($p=0; $p<$gap; $p++) {
		$s=$map{$seq[$p]};
		$tscore=$tscore+log($pwm[$j2][$s][$p]/$fr[$s]);
	    }
	    for ($p=$gap+1; $p<$lg; $p++) {
		$s=$map{$seq[$p]};
		$tscore=$tscore+log($pwm[$j2][$s][$p-1]/$fr[$s]);
	    }
	    if ($score<$tscore) {
		$score=$tscore;
		$pgap=$gap;
		$pcl=$j2;
	    }
				#Keep track of the score for this cluster
	    if ($j2==$j) {
		if ($score10<$tscore) {
		    $score10=$tscore;
		}
	    }
	}
    } elsif ($fixed_gap > 0) {
	$tscore=0;
	for ($p=0; $p<$fixed_gap; $p++) {
	    $s=$map{$seq[$p]};
	    $tscore=$tscore+log($pwm[$j2][$s][$p]/$fr[$s]);
	}
	for ($p=$fixed_gap+$lg-$tl; $p<$lg; $p++) {
	    $s=$map{$seq[$p]};
	    $tscore=$tscore+log($pwm[$j2][$s][$p-1]/$fr[$s]);
	}
	if ($score<$tscore) {
	    $score=$tscore;
	    $pgap=-1;
	    $pcl=$j2;
	}
	#Keep track of the score for this cluster
	if ($j2==$j) {
	    if ($score10<$tscore) {
		$score10=$tscore;
	    }
	}
    }
}


sub check_extensions(){

    my $l=$_[0];
    my @seq=split('', $l);
    my $j=$_[1];
    
    my @bulge_score=();
    my @Cterm_score=();
    my @Nterm_score=();

    my $j2;
    my $cd;
    my $Max=-100;
		    
    for ($j2=0; $j2<$ncluster; $j2++) {
		
			
	#Make sure it cannot be explained by a bulge on any allele
	#Only take positions P1,P2,P3 and and the two last ones
	$tscore=0;
	for ($p=0; $p<3; $p++) {
	    $s=$map{$seq[$p]};
	    $tscore=$tscore+log($pwm[$j2][$s][$p]/$fr[$s]);
	}
	for ($p=$lg-2; $p<$lg; $p++) {
	    $s=$map{$seq[$p]};
	    $tscore=$tscore+log($pwm[$j2][$s][$p+$tl-$lg]/$fr[$s]);
	    
	}
	push @bulge_score, $tscore;
	if ($j == $j2) {
	    $b_score=$tscore;
	    if ($b_score>$Max) {
		$Max=$b_score;
	    }
	}

	#Or by a C-terminal extension
	$tscore=0; 
	for ($p=0; $p<3; $p++) {
	    $s=$map{$seq[$p]};
	    $tscore=$tscore+log($pwm[$j2][$s][$p]/$fr[$s]);
	}
	for ($p=$tl-2; $p<$tl; $p++) {
	    $s=$map{$seq[$p]};
	    $tscore=$tscore+log($pwm[$j2][$s][$p]/$fr[$s]);
	  
	}
	push @Cterm_score, $tscore;
	if ($j == $j2) {
	    $C_score=$tscore;
	    if ($C_score>$Max) {
		$Max=$C_score;
	    }
	}
	
	#Or by a N-terminal extension
	$tscore=0;
	for ($p=0; $p<3; $p++) {
	    $s=$map{$seq[$p+$lg-$tl]};
	    $tscore=$tscore+log($pwm[$j2][$s][$p]/$fr[$s]);
	}
	for ($p=$tl-2; $p<$tl; $p++) {
	    $s=$map{$seq[$p+$lg-$tl]};
	    $tscore=$tscore+log($pwm[$j2][$s][$p]/$fr[$s]);
	}
	push @Nterm_score, $tscore;
	if ($j == $j2) {
	    $N_score=$tscore;
	    if ($N_score>$Max) {
		$Max=$N_score;
	    }
	}
	
    }
    
    
    $C_count=0;
    $N_count=0;
    $b_count=0;
    
    #Here we compare each score with all other possible explanations.
    #Basically, we count how often the score for one model of the jth motif is less than another score + Thresh_score. This number is at least one.
    #If it is equal to one, then the peptide can be unambiguously assigned. If not, there is ambiguity in the assignment.
    #We only select if the score for one model is significantly larger than all other models with any motif.
    for ($j2=0; $j2<$ncluster; $j2++) {
	if ($C_score < $bulge_score[$j2] + $Thresh_score) {
	    $C_count++;
	}
	if ($N_score < $bulge_score[$j2] + $Thresh_score) {
	    $N_count++;
	}
	if ($b_score < $bulge_score[$j2] + $Thresh_score) {
	    $b_count++;
	}
	if ($C_score < $Cterm_score[$j2] + $Thresh_score) {
	    $C_count++;
	}
	if ($N_score < $Cterm_score[$j2] + $Thresh_score) {
	    $N_count++;
	}
	if ($b_score < $Cterm_score[$j2] + $Thresh_score) {
	    $b_count++;
	}
	if ($C_score < $Nterm_score[$j2] + $Thresh_score) {
	    $C_count++;
	}
	if ($N_score < $Nterm_score[$j2] + $Thresh_score) {
	    $N_count++;
	}
	if ($b_score < $Nterm_score[$j2] + $Thresh_score) {
	    $b_count++;
	}
    }
    
    # Check if the peptide can be at least unambiguously assigned to the jth allele.
    # This includes cases where it could not be distinguished between the bulge, C-term or N-term, but at least assignement to other alleles could be excluded.
    # This will be useful to compare with the null model.
    $cd=1;
    for ($j2=0; $j2<$ncluster; $j2++) {
	if ($j2 != $j) {
	    if ($bulge_score[$j2] > $Max - $Thresh_score || $Cterm_score[$j2] > $Max - $Thresh_score || $Nterm_score[$j2] > $Max - $Thresh_score) {
		$cd=0;
	    }
	}
    }	    
    
    return($C_count, $N_count, $b_count, $C_score, $N_score, $b_score, $cd);
}

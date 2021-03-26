#!/usr/bin/perl
$DB = $ARGV[0];
$OUTPUT = $ARGV[1];
$path = `pwd`;
chomp $path;
$code = "";
for(my $i = 0; $i <= 3; $i++){
    $num = int(rand(10));
    $letter = ('A' .. 'Z')[26 * rand];
    $code .= $num . $letter;
}
# Make temporary Directory
if (!(-d "./tmp$code")){system("mkdir tmp$code")};
# Remove redundancy from DB
print "Remove redundancy from DB\n";
my $DB_NR = &remove_fasta_redundancy($DB);
# Get canonical and isoform sequences
print "Get canonical and isoform sequences\n";
$REF = &parser_canonical_isoform($DB_NR);
$CAN = $REF->{canonical};
$ISO = $REF->{isoform};
# Encode headers
print "Encode headers\n";
$REF = &encode_headers($ISO);
$ISO_encoded = $REF->{coded};
$ISO_table = $REF->{table};
# Trypsin digestion
print "Trypsin digestion\n";
$digested =  $path . "/tmp$code/digested.peptides";
system("digest -seqall $ISO_encoded -outfile $digested -menu 1 -mono no  -unfavoured" );
print "$digested\n";
$parsed_digested = &parser_peptides($digested);
# Uncode headers
print "Uncode headers\n";
$parsed_digested_uncoded = &uncode_headers($parsed_digested, $ISO_table);
# Remove Redundant peptides
print "Remove Redundant peptides\n";
$peptides_NR = &remove_redundant_peptides($parsed_digested_uncoded);
# Find proteotypic peptides
$peptides_NR = "tmp$code/NR_uncoded_parsed_digested.txt";
$CAN = "tmp$code/canonical.fasta";
print "Find proteotypic peptides\n";
$proteotypic = &find_proteotypic_peptides($peptides_NR, $CAN);
# Make OUTPUT
print "Make OUTPUT\n";
system("cat $uniprot_canonical $proteotypic_peptides > $OUTPUT");
# Encode custom database
print "Encode custom database\n";
$REF = &encode_headers($OUTPUT);
$OUTPUT_encoded = $REF->{coded};
$OUTPUT_table = $REF->{table};

####--------------------------------------------------------------------------------####
####--------------------------------------------------------------------------------####
sub remove_fasta_redundancy {
    my $pre_file = shift;
    my $file = "";
    if($pre_file =~ m/.+\/(.+)$/){$file = $1;} else {$file = $pre_file;}
    my $pre_path = `pwd`;
    chomp $pre_path;
    my $path = $pre_path . "/tmp$code/";
    my $count = 0;
    my $count_before = 0;
    my $header = "";
    my $seq = "";
    my %seq_header;
    my $file_non_redundant = "NR_" . $file;
    my $full_non_redundant = $path . $file_non_redundant;

    open(IN,"<",$pre_file) or die "cannot open the file ". $pre_file ."\n";
    while (my $line1 = <IN>) {
	chomp $line1;
	if ($line1 =~ m/>/){
	    $count++;
	    $count_before++;
	    if ($count > 1) {
		$header =~ s/>//;
		$seq_header{$seq}{$header} = 1;
		$seq = "";
	    }
	    $header = $line1;
	} else {
	    $seq .= $line1;
	}
    }
    $header =~ s/>//;
    $seq_header{$seq}{$header} = 1;
    # Printing the non_redundant file
    open(NR,">",$full_non_redundant) or die "cannot open the file! ". $full_non_redundant ."\n";
    foreach my $chave_seq (sort {$a cmp $b} keys %seq_header) {
	my $check = 0;
	my $header_concat = "";
	foreach my $chave_header (sort {$a cmp $b} keys %{$seq_header{$chave_seq}}) {
	    # Very first header or the only one
	    $header = $chave_header;
	    # Headers are concatened in this variable
	    $header_concat .= "\=\=\=" . $header;
	    # Check if that's only one sequence for this header or more
	    $check++;
	}
	# if check > 1 it means that the sequence is redudant (e.g, that are more than one sequence for this header)
	if ($check > 1) {
	    # The first "===" is removed
	    $header_concat =~ s/\=\=\=//;
	    # The sequence with all the respective headers are printed out at this file
	    print NR ">".$header_concat."\n".$chave_seq."\n";
	    # else, if the sequence is unique, the unique header is also printed out at the file
	} else {
	    print NR ">".$header."\n".$chave_seq."\n";
	}
    }
    close(NR);

    # Here is a control structere that shows how many sequences were before and after
    print "\tBefore: |$count_before|\t";
    print "\tAfter: |".scalar(keys(%seq_header))."|\t";
    my $total = ($count_before - scalar(keys(%seq_header)));
    print "\tRedundant Sequences: |$total|\n";
    # File containing the non redudant sequences are returned to the main program
    return($full_non_redundant);
}
####--------------------------------------------------------------------------------####
####--------------------------------------------------------------------------------####
sub remove_redundant_peptides {
    my $full_path = shift;
    $full_path =~  m/.+\/tmp$code\/(.+)/;
    my $file = $1;
    my $count = 0;
    my $count_before = 0;
    my $count_after = 0;
    my $header = "";
    my @peps = ();
    my %seq_header;
    my $pre_path = `pwd`;
    chomp $pre_path;
    my $path = $pre_path . "/tmp$code/";
    my $path_non_redundant = $path . "NR_" . $file;

    open(IN,"<",$full_path) or die "cannot open the file". $full_path."\n";
    while (my $line1 = <IN>) {
	chomp $line1;
	if ($line1 =~ m/>/){
	    $count++;
	    if ($count > 1) {
		$header =~ s/>//;
		foreach my $pep (@peps) {
		    $seq_header{$pep}{$header} = 1;
		}
		@peps = ();
	    }
	    $header =~ s/>//;
	    $header = $line1;
	} else {
	    my @split = split(/\s+/, $line1);
	    $count_before++;
	    push(@peps, $split[6]);
	}
    }

    foreach my $pep (@peps) {
	$seq_header{$pep}{$header} = 1;
    }

    open(NR,">",$path_non_redundant) or die "cannot open the file" . $path_non_redundant . "\n";
    foreach my $chave_seq (sort {$a cmp $b} keys %seq_header) {
	$count_after++;
	$count = 0;
	my $header_concat = "";
	foreach my $chave_header (sort {$a cmp $b} keys %{$seq_header{$chave_seq}}) {
	    $header = $chave_header;
	    $header_concat .= "\:\:\:" . $header;
	    $count++;
	}
	if ($count > 1) {
	    $header_concat =~ s/\:\:\://;
	    print NR ">".$header_concat."\n".$chave_seq."\n";
	} else {
	    print NR ">".$header."\n".$chave_seq."\n";
	}
    }
    close(NR);

    print "\tBefore: |$count_before|\t";
    print "\tAfter: |$count_after|\t";
    my $final = ($count_before)-($count_after);
    print "\tRedundant Peptides: |$final|\n";
    return($path_non_redundant);
}
####--------------------------------------------------------------------------------####
####--------------------------------------------------------------------------------####
sub parser_canonical_isoform {

    my $file = shift;
    my $count = 0;
    my $header = "";
    my $seq = "";
    my $pre_path = `pwd`;
    chomp $pre_path;
    my $path = $pre_path . "/tmp$code/";
    my $canonical = $path . "canonical.fasta";
    my $isoform = $path . "isoform.fasta";

    open(IN,"<", $file) or die "cannot open the file1 ".$file."\n";
    open(CAN,">",$canonical) or die "cannot open the file2 ".$canonical."\n";
    open(ISO,">",$isoform) or die "cannot open the file3 ".$isoform."\n";

    while(<IN>){
	chomp $_;
	if($_ =~ /^>/){
	    $can=0;
	    $iso=0;
	    if($_ =~ /sp\|[A-Z,0-9]+\|/){#if($_ =~ /\]\s+[A-Z,0-9]+\-1|]*\s+[A-Z,0-9]+\-1|sp\|[A-Z,0-9]\|/){
		print CAN "$_\n";
		$can=1;
	    }else{
		print ISO "$_\n";
		$iso=1;
	    }
	}elsif($can){
	    print CAN "$_\n";
	}elsif($iso){
	    print ISO "$_\n";
	}
    }
    close(IN);
    close(CAN);
    close(ISO);
    return({canonical => $canonical, isoform => $isoform});
}
####--------------------------------------------------------------------------------####
####--------------------------------------------------------------------------------####
sub encode_headers {

    my $full_path = shift;
    $full_path =~  m/.+\/(.+)/;
    my $file = $1;
    my $pre_path = `pwd`;
    chomp $pre_path;
    my $path = $pre_path . "/tmp$code/";
    my $coded = $path . "encoded_" . $file;
    my $table = $path . "table_" . $file;
    my $header = "";
    my $seq = "";
    my $count = 0;
    my $check = 0;

    open(IN, "<", $full_path) or die "cannot open the file " . $full_path . "\n";
    open(OUT, ">", $coded) or die "cannot open the file " . $coded . "\n";
    open(TAB, ">", $table) or die "cannot open the file " . $table . "\n";

    while (my $line = <IN>){
	chomp $line;
	if ($line =~ m/^>/) {
	    $check++;
	    if ($check > 1) {
		print OUT ">" . $count . "\n";
		print OUT "$seq\n";
		print TAB "$header!!$count\n";
	    }
	    $count++;
	    $header = $line;
	} else {
	    $seq = $line;
	}
    }

    print OUT ">" . $count . "\n";
    print OUT "$seq\n";
    print TAB "$header!!$count\n";

    close(IN);
    close(OUT);
    close(TAB);

    return( { coded => $coded, table => $table } );

}
####--------------------------------------------------------------------------------####
####--------------------------------------------------------------------------------####
sub uncode_headers {

    my $full_path = $_[0];
    my $table = $_[1];
    $full_path =~  m/.+\/tmp$code\/(.+)/;
    my $pre_path = `pwd`;
    chomp $pre_path;
    my $path = $pre_path . "/tmp$code/";
    my $file = $1;
    my $uncoded = $path . "uncoded_" . $file;
    my $header = "";
    my $header_real = "";
    my $seq = "";
    my $count = 0;
    my $check = 0;
    my %tb_cod_header;
    my $cod = "";
    my $check_pep = 0;
    my @peps = ();

    open(IN, "<", $full_path) or die "cannot open the file " . $full_path . "\n";
    open(OUT, ">", $uncoded) or die "cannot open the file " . $uncoded . "\n";
    open(TAB, "<", $table) or die "cannot open the file " . $table . "\n";

    while (my $line2 = <TAB>) {
	chomp $line2;
	my @split = split (/!!/, $line2);
	$split[0] =~ s/\s//;
	$split[0] =~ s/^>//;
	my $header = $split[0];
	$split[1] =~ s/\s//;
	my $codigo = $split[1];
	$tb_cod_header{$codigo} = $header;
    }

    close(TAB);
    while (my $line = <IN>){
	chomp $line;
	if ($line =~ m/^>/) {
	    $cod = $line;
	    $cod =~ s/^>//;
	    if (exists $tb_cod_header{$cod}) {
		$header_real = $tb_cod_header{$cod};
		print OUT ">$header_real\n";
	    }
	} else {
	    print OUT "$line\n";
	}
    }

    close(IN);
    close(OUT);
    close(TAB);

    return($uncoded);

}
####--------------------------------------------------------------------------------####
####--------------------------------------------------------------------------------####
sub parser_peptides {
    my $full_path = shift;
    my $file = 'digested.peptides';
    my $path = "./tmp$code/";
    my $parsed = $path . "parsed_" . $file;
    open (LOG, "> ./tmp$code/log_parsed_digested.txt");
    open (OUT, ">", $parsed);
    my $count_line = 0;
    my $check_1st_header = 0;
    my $check_2nd_header = 0;
    my $count_end = 0;
    my $check_hit_count = 0;
    my $file_hits = "";

    my $check_end = 0;
    my $header_check = 0;
    my $count_pep_nos_parametros = 0;
    my $count_header = 0;
    my $header = "";
    my @pep_dentro_dos_parametros = ();

    open (IN, "<", $full_path);
    while (<IN>){
	$count_line++;
	chomp($_);
	if (($_ =~ /^\#\={2}/) and ($check_1st_header == 0)){
	    $check_1st_header = 1;
	    $count_pep_nos_parametros = 0;
	}
	if ($check_1st_header == 1){
	    if ($_ =~ /^#\sSequence\:/){
		my @split_seq = split (/from:/);
		$header = $split_seq[0];
		$header =~ s/#//;
		$header =~ s/\s+//g;
		$header =~ s/Sequence\://;
		$count_header++;
	    }
	    if ($_ =~ /^#\sHitCount\:/){
		$file_hits = $_;
		$file_hits =~ s/# HitCount://;
		$file_hits =~ s/\s+//;
	    }
	    if (!(($header =~ /^$/) and ($file_hits =~ /^$/))){
		$header_check = 1;
	    }
	}
	if (($_ =~ /^\#\={2}/) and ($check_2nd_header == 0) and ($check_1st_header == 1) and ($header_check == 1)){
	    $check_2nd_header = 1;
	}
	if (($check_2nd_header == 1) and ($check_1st_header == 1) and ($header_check == 1)){
	    if ($_ =~ /^\s+\d+/){
		$check_hit_count++;
		my @split_pep = split (/\s+/, $_);
		my $length = length($split_pep[6]);
		if (($length >= 6) and ($length <= 24)){

		    push (@pep_dentro_dos_parametros, $_);
		    $count_pep_nos_parametros++;
		}
	    }
	}
	if (($_ =~ /^\#\-{2}/) and ($check_2nd_header == 1) and ($check_1st_header == 1) and ($file_hits == $check_hit_count) and ($header_check == 1)){
	    $count_end++;
	    if ($count_end == 2){
		$check_2nd_header = $check_1st_header = $check_hit_count = $count_end = $header_check = $check_end = 0;
		$file_hits = "";
		if ($count_pep_nos_parametros > 0){
		    print OUT ">$header\n";
		    foreach my $info (@pep_dentro_dos_parametros){
			print OUT "$info\n";
		    }
		}
		@pep_dentro_dos_parametros = ();

		if ($count_pep_nos_parametros == 0){
		    print LOG "Cabecalho $header nao teve peptideos entre 6 e 24 aminacidos\n";
		}
	    }
	}
	if ($count_line == 150){
	}
    }
    close(OUT);
    close(LOG);
    return($parsed);
}
####--------------------------------------------------------------------------------####
####--------------------------------------------------------------------------------####
sub find_proteotypic_peptides {

    my $peptides = $_[0];
    my $uniprot_canonical = $_[1];
    my %uniprot_seq;
    my $header;
    my $pep;
    my $seq;
    my $pre_path = `pwd`;
    chomp $pre_path;
    my $path = $pre_path . "/tmp$code/";
    my $proteotypic = $path . "proteotypic_peptides.fasta";

    open(UNI, "<", $uniprot_canonical) or die "cannot open the file" . $uniprot_canonical . "\n";
    while (my $line1 = <UNI>) {
	chomp $line1;
	unless ($line1 =~ m/^>/){
	    $seq = $line1;
	    $uniprot_seq{$seq}=1;
	}
    }
    close(UNI);

    open(PRO, ">", $proteotypic) or die "cannot open the file" . $proteotypic . "\n";
    open(PEP, "<", $peptides) or die "cannot open the file" . $peptides . "\n";
    while (my $line2 = <PEP>) {
	chomp $line2;
	if ($line2 =~ m/^>/){
	    $header = $line2;
	} else {
	    $pep = $line2;
	    my @pep_grep = grep (/$pep/, %uniprot_seq);
	    if (!(@pep_grep)) {
		print PRO "$header\n";
		print PRO "$pep\n";
	    }
	}
    }
    close(PEP);
    close(PRO);
    return($proteotypic);
}

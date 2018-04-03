#!/user/local/bin/perl -w

###########################################################################
# Date: 02-10-2016                                                        #
# Author: Xiaowei Yan                                                     #
# Run as: perl VariantsFromVCF_run.pl                                     #
# Program Description: 				                       	  #
#   Import VCF file, filter out unwanted loci, keep those known variants, #
#   output to an excel file.                                              #
###########################################################################

use strict;
use warnings;
#use Cwd;
use Spreadsheet::WriteExcel;
#use Data::Dumper;

if (!defined $ARGV[0]) {
	print "\nCould not find <run id> as in the argument\nTry use date format yymmdd as the <run id>\n\n";
	exit 0;
}

#import amino acid names
my $filename = "/media/ngs/YanXiaowei/VCF filtering/table/AminoAcidName.txt";
if (!(-f $filename)) {
	$filename = "N:/YanXiaowei/VCF filtering/table/AminoAcidName.txt";
}
my %AminoAcid = ();
open(FILE, "< $filename") || die "ERROR: Cannot open $filename for reading\n".$!;
while (my $line = <FILE>) {
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	my @fields = split(/\t/, $line);
	$AminoAcid{$fields[0]} = $fields[1];
}
close FILE;

#import artifact info
$filename = "/media/ngs/YanXiaowei/VCF filtering/table/Artifacts_10312015.txt";
if (!(-f $filename)) {
	$filename = "N:/YanXiaowei/VCF filtering/table/Artifacts_10312015.txt";
}
my %Artifact = ();
open(FILE, "< $filename") || die "ERROR: Cannot open $filename for reading\n".$!;
while (my $line = <FILE>) {
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	my @fields = split(/\t/, $line);
	$Artifact{$fields[1]."\_".$fields[2]."\_".$fields[3]."\_".$fields[4]} = 1;
}
close FILE;

my $out_xls = "XLS";
my $out_tsv = "TSV";
my $out_vcf = "VCF";

if (!(-d $out_xls)) {
	mkdir $out_xls;
}
if (!(-d $out_tsv)) {
	mkdir $out_tsv;
}
if (!(-d $out_vcf)) {
	mkdir $out_vcf;
}
if (!(-d "tmp")) {
	mkdir "tmp";
}

my $workbook_m = Spreadsheet::WriteExcel->new("$out_xls/MissedVariants_$ARGV[0].xls");
my $worksheet_m = $workbook_m->add_worksheet("mismatched");
my $row2 = 2;
my $color_flag = 50;
my $format_m = $workbook_m->add_format(bg_color => $color_flag);
$worksheet_m -> write("A1", "Accession No", $format_m);
$worksheet_m -> write("B1", "Biomarker", $format_m);
$worksheet_m -> write("C1", "Variant", $format_m);
$worksheet_m -> write("D1", "Run ID", $format_m);
$worksheet_m -> write("E1", "Notes", $format_m);

my $workbook_s = Spreadsheet::WriteExcel->new("$out_xls/MissedSamples_$ARGV[0].xls");
my $worksheet_s = $workbook_s->add_worksheet("missed samples");
my $row3 = 2;
my $format_s = $workbook_s->add_format(bg_color => $color_flag);
$worksheet_s -> write("A1", "Sample ID", $format_s);
$worksheet_s -> write("B1", "status", $format_s);

###    go through the whole dir

my %NGS_variant = ();
my %NGS_sample = ();

open(OUT_SAMPLES, "> samples_$ARGV[0].txt") || die "ERROR: Cannot open samples.txt for writing\n".$!;

###    go through the whole dir

my @run_dir = ('/media/ngs/MiSeqOne/','/media/ngs/MiSeqOne/2016/','/media/ngs/MiSeqOne/2014/','/media/ngs/MiSeqOne/2015/','/media/ngs/MiSeqTwo/','/media/ngs/MiSeqTwo/2014/','/media/ngs/MiSeqTwo/2015/','/media/ngs/NGS_patient_runs/2014/');
if (!(-d $run_dir[0])) {
	@run_dir = ('N:/MiSeqOne/','N:/MiSeqOne/2016/','N:/MiSeqOne/2014/','N:/MiSeqOne/2015/','N:/MiSeqTwo/','N:/MiSeqTwo/2014/','N:/MiSeqTwo/2015/','N:/NGS_patient_runs/2014/');
}

my $run_found = 0;
for (my $r=0; $r<scalar(@run_dir); $r++) {
	opendir(RUN_DIR, $run_dir[$r]) or die "ERROR: Cannot open $run_dir[$r] for reading\n".$!;
	while (my $run_name = readdir(RUN_DIR)) {
		chomp $run_name;
		next if (!($run_name =~ m/^$ARGV[0]_/));

		$run_found = 1;
		mkdir $out_xls if (!(-d $out_xls));
		mkdir $out_tsv if (!(-d $out_tsv));
		mkdir $out_vcf if (!(-d $out_vcf));
		mkdir "$out_xls/$ARGV[0]" if (!(-d "$out_xls/$ARGV[0]"));
		mkdir "$out_tsv/$ARGV[0]" if (!(-d "$out_tsv/$ARGV[0]"));
		mkdir "$out_vcf/$ARGV[0]" if (!(-d "$out_vcf/$ARGV[0]"));

		my $vcf_dir = $run_dir[$r].$run_name."/Data/Intensities/BaseCalls/Alignment";
		if (-d $vcf_dir) {
		}
		elsif (-d $run_dir[$r].$run_name."/Alignment") {
			$vcf_dir = $run_dir[$r].$run_name."/Alignment";
		}
		elsif (-d $run_dir[$r].$run_name."/Alignment2") {
			$vcf_dir = $run_dir[$r].$run_name."/Alignment2";
		}
		else {
			next;	# We only want directory
		}

		my @run_id = split(/_/, $run_name);
		if ($run_id[0] eq "T") {
			$run_id[0] = "T_".$run_id[1];
		}

######   import SCI case information

		my %SampleNames = ();	#hash of hashes for sample names and run names
		opendir(VCF_DIR, $vcf_dir) or die "ERROR: Cannot open $vcf_dir for reading\n".$!;
		while (my $vcf_file = readdir(VCF_DIR)) {
			chomp $vcf_file;
			next unless (-f "$vcf_dir/$vcf_file");	# We only want files
			next unless ($vcf_file =~ m/\.genome\.vcf$/);	# Use a regular expression to find files ending in .genome.vcf
			next if ($vcf_file =~ m/^TSCA\-Ctrl/ or $vcf_file =~ m/^Water/ or $vcf_file =~ m/^Coriel/);
			my @fields_filename = split("\_", $vcf_file);
			my @fields_samplename = split("-", $fields_filename[0]);
			if (defined $fields_samplename[1]) {
				$fields_samplename[0] .= "-".$fields_samplename[1];
			}
			$SampleNames{$fields_samplename[0]} = 1;
		}
		closedir(VCF_DIR);

		if (!defined $ARGV[0]) {
			print "\nCould not find <run id> as in the argument\nTry use date format yymmdd as the <run id>\n";
			print "Try to run as:  perl VariantsFromVCF_run.pl yymmdd [case_database]\n\n";
			exit 0;
		}

		$filename = "/media/ngs/YanXiaowei/VCF filtering/table/Cases_03282016.txt";
		if (!(-f $filename)) {
			$filename = "N:/YanXiaowei/VCF filtering/table/Cases_03282016.txt";
		}
		if (defined $ARGV[1]) {
			$filename = $ARGV[1];
			if (!(-f $filename)) {
				$filename = "/media/ngs/MiniVCF/doc/".$ARGV[1];
			}
			if (!(-f $filename)) {
				$filename = "/media/ngs/YanXiaowei/VCF filtering/table/".$ARGV[1];
			}
			if (!(-f $filename)) {
				$filename = "N:/MiniVCF/doc/".$ARGV[1];
			}
			if (!(-f $filename)) {
				$filename = "N:/YanXiaowei/VCF filtering/table/".$ARGV[1];
			}
		}

		open(FILE, "< $filename") || die "$!\nPlease Run as: perl VariantsFromVCF_run.pl yymmdd [case_database]\nERROR: Cannot open case database: $filename\n".$!;

		while (my $line = <FILE>) {
			$line =~ s/\r//g;
			$line =~ s/\n//g;
			my @fields = split(/\t/, $line);
			next if (!defined $fields[0]);
			next if (!defined $SampleNames{$fields[0]});

			if (!defined $fields[2]) {
				$NGS_sample{$fields[0]} = 0 if (!defined $NGS_sample{$fields[0]});
				next;
			}
			if ($fields[2] eq "N/A") {
				$NGS_sample{$fields[0]} = 0 if (!defined $NGS_sample{$fields[0]});
				next;
			}

			$fields[2] =~ s/^ //;
			$fields[2] =~ s/ $//;
			if (0 eq index $fields[2], $fields[1]) {
				$fields[2] = substr $fields[2], length($fields[1]);
			}
			$fields[2] =~ s/^\_//;
			$fields[2] =~ s/^ //;

			if ($fields[2] =~ m/fs/) {
				$fields[2] = substr $fields[2], 0, index($fields[2], "fs");
				while (length($fields[2])>0 && ($fields[2] =~ /[^0-9]$/)) {
					$fields[2] = substr($fields[2], 0, length($fields[2])-1);
				}
				$fields[2] .= "fs";
			}
			elsif ($fields[2] =~ m/^splice site /) {
				$fields[2] = substr $fields[2], 12;
			}
			elsif (substr($fields[2],2,8) eq "-UTR_Del") {
				$fields[2] = substr($fields[2],2);
			}
			elsif (substr($fields[2],3,8) eq "-UTR_Del") {
				$fields[2] = substr($fields[2],3);
			}
			elsif ($fields[2] eq "DPYD\*2A") {
				$fields[2] = "1905+1G>A";
			}
			elsif ($fields[2] =~ m/delinsP/) {
				$fields[2] = "L755fs";
			}
			elsif ($fields[2] =~ m/del$/) {
				my @f = split(/_/, $fields[2]);
				if (!defined $f[1]) {
					$fields[2] = substr($f[0],0,-3)."_".substr($f[0],0,-3)."del";
				}
			}
			$NGS_variant{$fields[0]}{$fields[1]}{$fields[2]} = 0;
			if (!defined $NGS_sample{$fields[0]}) {
				$NGS_sample{$fields[0]} = 0;
			}
		}
		close FILE;

################

		my %samples = ();	#hash of hashes for sample names and run names

		opendir(VCF_DIR, $vcf_dir) or die "ERROR: Cannot open $vcf_dir for reading\n".$!;

		my $file_count = 0;
		my $line_count = 0;

		while (my $vcf_file = readdir(VCF_DIR)) {
			chomp $vcf_file;
			next unless (-f "$vcf_dir/$vcf_file");	# We only want files
			next unless ($vcf_file =~ m/\.genome\.vcf$/);	# Use a regular expression to find files ending in .genome.vcf

			next if ($vcf_file =~ m/^TSCA\-Ctrl/ or $vcf_file =~ m/^Water/ or $vcf_file =~ m/^Coriel/);

###########################
#next unless ($vcf_file =~ m/^S16\-083334/);
###########################


			my $run_name = $vcf_file;
			$run_name =~ s/\.genome\.vcf$//;
			my @fields_filename = split("\_", $vcf_file);
			my @fields_samplename = split("-", $fields_filename[0]);
			my $sample_name = $fields_samplename[0];
			if (defined $fields_samplename[1]) {
				$sample_name .= "-".$fields_samplename[1];
			}
			if (defined $NGS_sample{$sample_name}) {
				$NGS_sample{$sample_name} = 1;
			}
			else {
#				next;
			}

			if (defined $samples{$sample_name}) {
				$samples{$sample_name}{"run_count"} ++;
			}
			else {
				$samples{$sample_name}{"run_count"} = 1;
			}
			$samples{$sample_name}{$samples{$sample_name}{"run_count"} - 1} = $run_name;

			print $run_name."\t".$run_dir[$r].$run_id[0]."\n";
			print OUT_SAMPLES $sample_name."\t".$run_name."\t".$run_id[0]."\t".$run_dir[$r]."\n";

			my $field_number = 7;
			my %header_fields = ();

			open(VCF_FILE, "< $vcf_dir/$vcf_file") || die "ERROR: Cannot open $vcf_dir/$vcf_file for reading\n".$!;
			my $output_file = "tmp/temp.vcf";
			open(OUTPUT_FILE, "> $output_file") || die "ERROR: Cannot open $output_file for writing\n".$!;

			while (my $line = <VCF_FILE>) {
				$line =~ s/\r//g;
				$line =~ s/\n//g;
				if ($line =~ m/^\#\#FORMAT=.*/ || $line =~ m/^\#\#INFO=.*/) {
					my @fields_line = split("ID=", $line);
					my @fields_ID = split(",", $fields_line[1]);
					$header_fields{$fields_ID[0]} = $field_number;
					$field_number++;
					print OUTPUT_FILE $line."\n";
					next;
				}
				if ($line =~ m/^\#\#.*/) {
					print OUTPUT_FILE $line."\n";
					next;
				}
				if ($line =~ m/^\#.*/) {
					my @fields_headline = split("\#", $line);
					@fields_headline = split("\t", $fields_headline[1]);
					for (my $i=0; $i<7; $i++) {
						$header_fields{$fields_headline[$i]} = $i;
					}
					print OUTPUT_FILE $line."\n";
					next;
				}

				my @fields_data = split("\t", $line);
				if ($fields_data[4] eq ".") {
				}
				else {
					print OUTPUT_FILE $line."\n";
				}
			}
			close OUTPUT_FILE;
			close VCF_FILE;
			if (-f "/media/ngs/YanXiaoWei/snpEff/snpEff.jar") {
				system("java -Xmx16G -jar /media/ngs/YanXiaoWei/snpEff/snpEff.jar hg19 $output_file > tmp/ann.vcf");
			}
			else {
				system("java -Xmx16G -jar N:/YanXiaoWei/snpEff/snpEff.jar hg19 $output_file > tmp/ann.vcf");
			}
			###    reading ann.vcf   ########

			open(VCF_FILE, "< tmp/ann.vcf") || die "ERROR: Cannot open ann.vcf for reading\n".$!;
			$output_file = "tmp/".$run_name."_variants_".$run_id[0].".txt";
			open(OUTPUT_FILE, "> $output_file") || die "ERROR: Cannot open $output_file for writing\n".$!;

			$field_number = 7;
			for (keys %header_fields) {
				delete $header_fields{$_};
			}
			my @ANN_fields = ();
			my @LOF_fields = ();

			$line_count = 0;
			while (my $line = <VCF_FILE>) {
				$line_count ++;

				$line =~ s/\r//g;
				$line =~ s/\n//g;

				if ($line =~ m/^\#\#FORMAT=.*/) {
					my @fields_line = split("ID=", $line);
					my @fields_ID = split(",", $fields_line[1]);
					$header_fields{$fields_ID[0]} = $field_number;
					$field_number++;
					next;
				}
				elsif ($line =~ m/^\#\#INFO=.*/) {
					my @fields_line = split("ID=", $line);
					@fields_line = split(/'/, $fields_line[1]);
					my @fields_ID = split(/,/, $fields_line[0]);
					if ($fields_ID[0] eq "ANN") {
						$fields_line[1] =~ s/ //g;
						@fields_ID = split(/\|/, $fields_line[1]);
						for (my $j=0; $j<scalar(@fields_ID); $j++) {
							$ANN_fields[$j] = $fields_ID[$j];
							$header_fields{$ANN_fields[$j]} = $field_number;
							$field_number++;
						}
					}
					elsif ($fields_ID[0] eq "LOF") {
						$fields_line[1] =~ s/ //g;
						@fields_ID = split(/\|/, $fields_line[1]);
						for (my $j=0; $j<scalar(@fields_ID); $j++) {
							$LOF_fields[$j] = $fields_ID[$j]."_LOF";
							$header_fields{$LOF_fields[$j]} = $field_number;
							$field_number++;
						}
					}
					elsif ($fields_ID[0] eq "REF") {
						$header_fields{$fields_ID[0]."_INFO"} = $field_number;
						$field_number++;
					}
					else {
						$header_fields{$fields_ID[0]} = $field_number;
						$field_number++;
					}
					next;
				}
				elsif ($line =~ m/^\#\#.*/) {
					next;
				}
				elsif ($line =~ m/^\#.*/) {
					my @fields_headline = split(/\#/, $line);
					@fields_headline = split(/\t/, $fields_headline[1]);
					for (my $i=0; $i<7; $i++) {
						$header_fields{$fields_headline[$i]} = $i;
					}
					my ($key) = grep {0 == $header_fields{$_}} keys %header_fields;
					print OUTPUT_FILE $key;

					for (my $i=1; $i<$field_number; $i++) {
						($key) = grep {$i == $header_fields{$_}} keys %header_fields;
						print OUTPUT_FILE "\t".$key;
					}
					print OUTPUT_FILE "\n";
					next;
				}

				my @fields_data = split(/\t/, $line);
				my @INFO_data = split(/\;/, $fields_data[7]);
				my @FORMAT_names = split(/\:/, $fields_data[8]);
				my @FORMAT_values = split(/\:/, $fields_data[9]);
				my %dataline = ();

				for (my $i=0; $i<$field_number; $i++) {
					my ($k) = grep {$i == $header_fields{$_}} keys %header_fields;
					if ($i<7) {
						$dataline{$k} = $fields_data[$i];
					}
					else {
						$dataline{$k} = "-";
					}
				}

				for (my $i=0; $i<scalar(@INFO_data); $i++) {
					my @INFO_ID = split(/=/, $INFO_data[$i]);
					if ($INFO_ID[0] eq "ANN") {
						my @ANN = split(/\,/, $INFO_ID[1]);
						for (my $k=0; $k<scalar(@ANN); $k++) {
							my @fields_value = split(/\|/, $ANN[$k]);
							for (my $j=0; $j<scalar(@ANN_fields); $j++) {
								if ($k>0) {
									$dataline{$ANN_fields[$j]} = $dataline{$ANN_fields[$j]}."|";
								}
								else {
									$dataline{$ANN_fields[$j]} = "";
								}
								if (defined $fields_value[$j]) {
									$dataline{$ANN_fields[$j]} = $dataline{$ANN_fields[$j]}.$fields_value[$j];
								}
							}
						}
					}
					elsif ($INFO_ID[0] eq "LOF") {
						$INFO_ID[1] =~ s/^\(//;
						$INFO_ID[1] =~ s/\)$//;
						my @fields = split(/\|/, $INFO_ID[1]);
						for (my $j=0; $j<scalar(@LOF_fields); $j++) {
							if (defined $fields[$j]) {
								$dataline{$LOF_fields[$j]} = $fields[$j];
							}
						}
					}
					elsif ($INFO_ID[0] eq "REF") {
						$dataline{$INFO_ID[0]."_INFO"} = $INFO_ID[0];
					}
					elsif (defined $INFO_ID[1]) {
						$dataline{$INFO_ID[0]} = $INFO_ID[1];
					}
					else {
						$dataline{$INFO_ID[0]} = $INFO_ID[0];
					}
				}

				for (my $i=0; $i<scalar(@FORMAT_names); $i++) {
					if ($FORMAT_names[$i] eq "GT") {
						if ($FORMAT_values[$i] eq "0/1") {
							$dataline{$FORMAT_names[$i]} = "het";
						}
						elsif ($FORMAT_values[$i] eq "1/1") {
							$dataline{$FORMAT_names[$i]} = "hom";
						}
						else {
							$dataline{$FORMAT_names[$i]} = $FORMAT_values[$i];
						}
					}
					else {
						$dataline{$FORMAT_names[$i]} = $FORMAT_values[$i];
					}
				}

				my ($k) = grep {0 == $header_fields{$_}} keys %header_fields;
				print OUTPUT_FILE $dataline{$k};

				for (my $n=1; $n<$field_number; $n++) {
					($k) = grep {$n == $header_fields{$_}} keys %header_fields;
					print OUTPUT_FILE "\t".$dataline{$k};
				}
				print OUTPUT_FILE "\n";
			}
			close OUTPUT_FILE;
			close VCF_FILE;
		}
		closedir(VCF_DIR);


		###   calculate variant frequency   ###

		my $workbook = "";
		my $worksheet = "";

		my %variant_freq = ();
		my $total_samples = 0;

		foreach my $sample ( keys %samples ) {
			$total_samples++ if (!($sample =~ m/^Water.*/ || $sample =~ m/^Coriel.*/ || $sample =~ m/^TSCA-Ctrl.*/));

			my %variants = ();	#hash of keys for each unique variant

			for (my $i=0; $i<$samples{$sample}{"run_count"}; $i++) {
				my $line_count = 0;
				my @local_head_fields = ();

				open(VARIANT_FILE, "< tmp/".$samples{$sample}{$i}."_variants_".$run_id[0].".txt") || die "ERROR: Cannot open $samples{$sample}{$i}_variants.txt for reading\n".$!;
				while (my $line = <VARIANT_FILE>) {
					$line_count ++;
					if ($line_count == 1) {
						next;
					}
					$line =~ s/\r//g;
					$line =~ s/\n//g;
					my @data_fields = split("\t", $line);
					my $key = $data_fields[0]."_".$data_fields[1]."_".$data_fields[3]."_".$data_fields[4];
					if (defined $variants{$key}) {
						next;
					}
					else {
						$variants{$key} = 1;
						if (defined $variant_freq{$key}) {
							if ($sample =~ m/^Water.*/ || $sample =~ m/^Coriel.*/ || $sample =~ m/^TSCA-Ctrl.*/) {
							}
							else {
								$variant_freq{$key}{'count'} ++;
							}
							$variant_freq{$key}{'sample'} .= "\|$sample";

						}
						else {
							if ($sample =~ m/^Water.*/ || $sample =~ m/^Coriel.*/ || $sample =~ m/^TSCA-Ctrl.*/) {
								$variant_freq{$key}{'count'} = 0;
							}
							else {
								$variant_freq{$key}{'count'} = 1;
							}
							$variant_freq{$key}{'sample'} = $sample;
						}
					}
				}
			}
		}



		##   output to an excel and tsv file    ##

		foreach my $sample ( keys %samples ) {
			my %header_fields = ();
			my $header_fields_count = 0;
			my %variants = ();	#hash of hashes for each unique variant record

			for (my $n=0; $n<$samples{$sample}{"run_count"}; $n++) {
				foreach my $var ( keys %variants ) {
					$variants{$var}{'flag'} = 0;
				}
				my $line_count = 0;
				my @local_head_fields = ();

				open(VARIANT_FILE, "< tmp/".$samples{$sample}{$n}."_variants_".$run_id[0].".txt") || die "ERROR: Cannot open $samples{$sample}{$n}_variants.txt for reading\n".$!;
				while (my $line = <VARIANT_FILE>) {
					$line_count ++;
					$line =~ s/\r//g;
					$line =~ s/\n//g;
					if ($line_count == 1) {
						@local_head_fields = split(/\t/, $line);
						for (my $j=0; $j<scalar(@local_head_fields); $j++) {
							if (!defined $header_fields{$local_head_fields[$j]}) {
								$header_fields{$local_head_fields[$j]} = $header_fields_count;
								$header_fields_count ++;
							}
						}
						next;
					}

					my @data_fields = split("\t", $line);
					my %data_temp = ();
					for (my $j=0; $j<scalar(@local_head_fields); $j++) {
						$data_temp{$local_head_fields[$j]} = $data_fields[$j];
					}
					my @GI = split(/,/, $data_temp{'GI'});
					for (my $i=1; $i<scalar(@GI); $i++) {
						if ($GI[0] =~ m/-AS/ and defined$GI[$i]) {
							$GI[0] = $GI[$i];
						}
						else {
							last;
						}
					}
					my $key = $data_fields[0]."_".$data_fields[1]."_".$data_fields[3]."_".$data_fields[4];


					##########  match to SCI case  #############

					my $case_matched = 0;
					my $target = "";
					my $target_var = "";
					my @v = ();
					my @HGVSp = split(/\|/, $data_temp{'HGVS.p'});
					for (my $i=0; $i<scalar(@HGVSp); $i++) {
						if ($HGVSp[$i] =~ m/^p\./) {
							@v = split(/\./, $HGVSp[$i]);
							$v[1] =~ s/\*/Ter/g;
							my $start = substr $v[1], 0, 3;
							my $pos = substr $v[1], 3, length($v[1])-6;
							my $end = substr $v[1], -3;
							if ($v[1] =~ m/fs/) {
								$pos = substr $v[1], 3, index($v[1], "fs")-3;
								if (! $pos =~ /^[0-9,.E]+$/) {
									$pos = substr $pos, 0, length($pos)-3;
								}
								$end = "fs";
							}
							elsif ($v[1] =~ m/delinsdel/) {
								$end = "delinsdel";
								$v[1] = substr $v[1], 0, index($v[1], $end);
								my @fields = split(/\_/, $v[1]);
								$end = "del";

								$pos = substr $fields[0], 3;
								my $start_temp = "";
								if (defined $fields[1]) {
									$start_temp = substr $fields[1], 0, 3;
									if (defined $AminoAcid{$start_temp}) {
										$start_temp = $AminoAcid{$start_temp};
									}
									$pos = $pos."_".$start_temp.substr($fields[1],3);
								}
							}
							elsif ($v[1] =~ m/delins/) {
								my @fields = split("delins", $v[1]);
								$end = "";
								my $start_temp = "";
								if (defined $fields[1]) {
									while (length($fields[1])>=3) {
										$start_temp = substr $fields[1], 0, 3;
										if (defined $AminoAcid{$start_temp}) {
											$end = $end.$AminoAcid{$start_temp};
										}
										else {
											$end = $end.$start_temp;
										}
										$fields[1] = substr $fields[1],3;
									}
									if (length($fields[1])>=1) {
										$end = $end.$start_temp;
									}
								}
								$end = "delins".$end;

								@fields = split(/\_/, $fields[0]);
								$pos = substr $fields[0], 3;
								$start_temp = "";
								if (defined $fields[1]) {
									$start_temp = substr $fields[1], 0, 3;
									if (defined $AminoAcid{$start_temp}) {
										$start_temp = $AminoAcid{$start_temp};
									}
									$pos = $pos."_".$start_temp.substr($fields[1], 3);
								}
							}
							elsif ($v[1] =~ m/ins/) {
								my @fields = split("ins", $v[1]);
								$end = "";
								my $start_temp = "";
								if (defined $fields[1]) {
									while (length($fields[1])>=3) {
										$start_temp = substr $fields[1], 0, 3;
										if (defined $AminoAcid{$start_temp}) {
											$end = $end.$AminoAcid{$start_temp};
										}
										else {
											$end = $end.$start_temp;
										}
										$fields[1] = substr $fields[1],3;
									}
									if (length($fields[1])>=1) {
										$end = $end.$start_temp;
									}
								}
								$end = "ins".$end;

								@fields = split(/\_/, $fields[0]);
								$pos = substr $fields[0], 3;
								$start_temp = "";
								if (defined $fields[1]) {
									$start_temp = substr $fields[1], 0, 3;
									if (defined $AminoAcid{$start_temp}) {
										$start_temp = $AminoAcid{$start_temp};
									}
									$pos = $pos."_".$start_temp.substr($fields[1], 3);
								}
							}
							elsif ($end eq "Ter") {
								$end = "\*";
							}
							elsif ($end eq "dup" || $end eq "del") {
								$v[1] = substr $v[1], 0, index($v[1], $end);
								my @fields = split(/\_/, $v[1]);

								$pos = substr $fields[0], 3;
								my $start_temp = "";
								if (defined $fields[1]) {
									$start_temp = substr $fields[1], 0, 3;
									if (defined $AminoAcid{$start_temp}) {
										$start_temp = $AminoAcid{$start_temp};
									}
									$pos = $pos."_".$start_temp.substr($fields[1],3);
								}
								else {
									$start_temp = $start;
									if (defined $AminoAcid{$start_temp}) {
										$start_temp = $AminoAcid{$start_temp};
									}
									$pos = $pos."_".$start_temp.$pos;
								}
							}

							if (defined $AminoAcid{$end}) {
								$end = $AminoAcid{$end};
							}
							if (defined $AminoAcid{$start}) {
								$start = $AminoAcid{$start};
							}

							if (defined $NGS_variant{$sample}{$GI[0]}{$start.$pos.$end}) {
								$case_matched = 1;
								$NGS_variant{$sample}{$GI[0]}{$start.$pos.$end} = 1;
								$target = $GI[0];
								$target_var = $start.$pos.$end;
							}
							elsif (defined $NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$start.$pos.$end}) {
								$case_matched = 1;
								$NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$start.$pos.$end} = 1;
								$GI[0] = $data_temp{'Gene_Name'};
								$target = $GI[0];
								$target_var = $start.$pos.$end;
							}
							elsif ($end eq "fs") {
								$end = "\*";
								if (defined $NGS_variant{$sample}{$GI[0]}{$start.$pos.$end}) {
									$case_matched = 1;
									$NGS_variant{$sample}{$GI[0]}{$start.$pos.$end} = 1;
									$target = $GI[0];
									$target_var = $start.$pos.$end;
								}
								elsif (defined $NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$start.$pos.$end}) {
									$case_matched = 1;
									$NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$start.$pos.$end} = 1;
									$GI[0] = $data_temp{'Gene_Name'};
									$target = $GI[0];
									$target_var = $start.$pos.$end;
								}
							}
						}
						last if ($case_matched==1);
					}

					if ($case_matched==0) {
						my @HGVSc = split(/\|/, $data_temp{'HGVS.c'});
						for (my $i=0; $i<scalar(@HGVSc); $i++) {
							if ($HGVSc[$i] =~ m/c\./) {
								@v = split("c\.", $HGVSc[$i]);
								if (defined $NGS_variant{$sample}{$GI[0]}{$v[1]}) {
									$case_matched = 1;
									$NGS_variant{$sample}{$GI[0]}{$v[1]} = 1;
									$target = $GI[0];
									$target_var = $v[1];
								}
								elsif (defined $NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$v[1]}) {
									$case_matched = 1;
									$NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$v[1]} = 1;
									$GI[0] = $data_temp{'Gene_Name'};
									$target = $GI[0];
									$target_var = $v[1];
								}
							}
						}
					}

					if ($case_matched==0) {
						@v = split(/\,/, $data_temp{'FC'});
						for (my $i = 0; $i<scalar(@v); $i++) {
							my @f = split(/\_/, $v[$i]);
							$f[0] = $f[1] if (defined $f[1]);
							$f[0] =~ s/X/\*/g;
							if (defined $NGS_variant{$sample}{$GI[0]}{$f[0]}) {
								$case_matched = 1;
								$NGS_variant{$sample}{$GI[0]}{$f[0]} = 1;
								$target = $GI[0];
								$target_var = $f[0];
								my @v2 = split(/\,/, $data_temp{'TI'});
								if (defined $v2[$i]) {
									$data_temp{'HGVS.p'} = $v2[$i];
								}
								else {
									$data_temp{'HGVS.p'} = "";
								}
								$data_temp{'HGVS.p'} = $data_temp{'HGVS.p'}."_".$f[0];
							}
							elsif (defined $NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$f[0]}) {
								$case_matched = 1;
								$NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$f[0]} = 1;
								$GI[0] = $data_temp{'Gene_Name'};
								$target = $GI[0];
								$target_var = $f[0];
								my @v2 = split(/\,/, $data_temp{'TI'});
								if (defined $v2[$i]) {
									$data_temp{'HGVS.p'} = $v2[$i];
								}
								else {
									$data_temp{'HGVS.p'} = "";
								}
								$data_temp{'HGVS.p'} = $data_temp{'HGVS.p'}."_".$f[0];
							}
							last if ($case_matched==1);
						}
					}

					if ($case_matched==0) {
						if  ($data_temp{'HGVS.c'} =~ m/del/) {
							my $var = "-UTR_Del".substr($data_temp{'REF'},1);
							if (defined $NGS_variant{$sample}{$GI[0]}{$var}) {
								$case_matched = 1;
								$NGS_variant{$sample}{$GI[0]}{$var} = 1;
								$target = $GI[0];
								$target_var = $var;
							}
							elsif (defined $NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$var}) {
								$case_matched = 1;
								$NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$var} = 1;
								$GI[0] = $data_temp{'Gene_Name'};
								$target = $GI[0];
								$target_var = $var;
							}
							$var = "\*447_\*452del".substr($data_temp{'REF'},1);
							if (defined $NGS_variant{$sample}{$GI[0]}{$var}) {
								$case_matched = 1;
								$NGS_variant{$sample}{$GI[0]}{$var} = 1;
								$target = $GI[0];
								$target_var = $var;
							}
							elsif (defined $NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$var}) {
								$case_matched = 1;
								$NGS_variant{$sample}{$data_temp{'Gene_Name'}}{$var} = 1;
								$GI[0] = $data_temp{'Gene_Name'};
								$target = $GI[0];
								$target_var = $var;
							}
						}
					}

					if (defined $variants{$key}) {
						$variants{$key}{"ConcordanceCount"} ++;
						for (my $j=0; $j<scalar(@local_head_fields); $j++) {
							if ($local_head_fields[$j] eq "QUAL" or
									$local_head_fields[$j] eq "FILTER" or
									$local_head_fields[$j] eq "GQX" or
									$local_head_fields[$j] eq "GT" or
									$local_head_fields[$j] eq "GQ" or
									$local_head_fields[$j] eq "AD" or
									$local_head_fields[$j] eq "VF" or
									$local_head_fields[$j] eq "NL" or
									$local_head_fields[$j] eq "SB" or
									$local_head_fields[$j] eq "DP" ) {
								$variants{$key}{$local_head_fields[$j]} .= "|$data_fields[$j]";
							}
						}
					}
					else {
						$variants{$key}{"ConcordanceCount"} = 1;
						$variants{$key}{'CaseMatched'} = $case_matched;
						for (my $j=0; $j<scalar(@local_head_fields); $j++) {
							$variants{$key}{$local_head_fields[$j]} = "";
							if ($local_head_fields[$j] eq "QUAL" or
									$local_head_fields[$j] eq "FILTER" or
									$local_head_fields[$j] eq "GQX" or
									$local_head_fields[$j] eq "GT" or
									$local_head_fields[$j] eq "GQ" or
									$local_head_fields[$j] eq "VF" or
									$local_head_fields[$j] eq "NL" or
									$local_head_fields[$j] eq "SB" or
									$local_head_fields[$j] eq "DP" ) {
								for (my $k=0; $k<$n; $k++) {
									$variants{$key}{$local_head_fields[$j]} .= "0|";
								}
							}
							elsif ($local_head_fields[$j] eq "AD") {
								for (my $k=0; $k<$n; $k++) {
									$variants{$key}{$local_head_fields[$j]} .= "0,0|";
								}
							}
							$variants{$key}{$local_head_fields[$j]} .= $data_fields[$j];
						}
					}
					$variants{$key}{'freq'} = $variant_freq{$key}{count};
					$variants{$key}{flag} = 1;
					if ($case_matched == 1) {
						$variants{$key}{target} = $target;
						$variants{$key}{var} = $target_var;
					}
				}
				foreach my $var ( keys %variants ) {
					if ($variants{$var}{'flag'} == 0) {
						for (my $j=0; $j<scalar(@local_head_fields); $j++) {
							if ($local_head_fields[$j] eq "QUAL" or
									$local_head_fields[$j] eq "FILTER" or
									$local_head_fields[$j] eq "GQX" or
									$local_head_fields[$j] eq "GT" or
									$local_head_fields[$j] eq "GQ" or
									$local_head_fields[$j] eq "VF" or
									$local_head_fields[$j] eq "NL" or
									$local_head_fields[$j] eq "SB" or
									$local_head_fields[$j] eq "DP" ) {
								$variants{$var}{$local_head_fields[$j]} .= "|0";
							}
							elsif ($local_head_fields[$j] eq "AD") {
								$variants{$var}{$local_head_fields[$j]} .= "|0,0";
							}
						}
					}
				}
				close VARIANT_FILE;
			}



			##########   output variants to xls file

			for (my $n=0; $n<$samples{$sample}{"run_count"}; $n++) {
				my @worksheet_columns = ("A".."ZZ");
				$workbook = Spreadsheet::WriteExcel->new("$out_xls/$run_id[0]/$samples{$sample}{$n}"."_".$run_id[0].".xls");
				$worksheet = $workbook->add_worksheet($samples{$sample}{$n});
				my $color_flag = 1;
				my $format=$workbook->add_format(bg_color => "grey");

				$worksheet -> write("A1", "Gene", $format);
				$worksheet -> write("B1", "Variant", $format);
				$worksheet -> write("C1", "Chr", $format);
				$worksheet -> write("D1", "Coordinate", $format);
				$worksheet -> write("E1", "Type", $format);
				$worksheet -> write("F1", "Genotype", $format);
				$worksheet -> write("G1", "Exonic", $format);
				$worksheet -> write("H1", "Filters", $format);
				$worksheet -> write("I1", "Quality", $format);
				$worksheet -> write("J1", "GQX", $format);
				$worksheet -> write("K1", "Variant Freq", $format);
				$worksheet -> write("L1", "Read Depth", $format);
				$worksheet -> write("M1", "Alt Read Depth", $format);
				$worksheet -> write("N1", "Allelic Depth", $format);
				$worksheet -> write("O1", "dbSNP ID", $format);
				$worksheet -> write("P1", "Consequence", $format);
				$worksheet -> write("Q1", "HGVSc", $format);
				$worksheet -> write("R1", "HGVSp", $format);
				$worksheet -> write("S1", "Reported", $format);

				$worksheet -> freeze_panes(1, 0);

				my $row = 1;
				for my $key (sort {$variants{$b}{ConcordanceCount}<=>$variants{$a}{ConcordanceCount}} keys %variants){
					if (!defined $variants{$key}) {
						next;
					}
					if ($variants{$key}{'freq'}/$total_samples > 0.8) {
						if ($variants{$key}{'CaseMatched'}==0) {
							next;
						}
						else {
							$NGS_variant{$sample}{$variants{$key}{target}}{$variants{$key}{var}} = 2;
						}
					}
					if (defined $Artifact{$key}) {
						if ($variants{$key}{'CaseMatched'}==0) {
							next;
						}
						else {
							$NGS_variant{$sample}{$variants{$key}{target}}{$variants{$key}{var}} = 3;
						}
					}
					if ($variants{$key}{"ConcordanceCount"} < $samples{$sample}{"run_count"} && $samples{$sample}{"run_count"} > 1) {
						if ($variants{$key}{'CaseMatched'}==0) {
							next;
						}
						else {
							$NGS_variant{$sample}{$variants{$key}{target}}{$variants{$key}{var}} = 4;
						}
					}
					my @SB = split(/\|/, $variants{$key}{"SB"});
					my @DP = split(/\|/, $variants{$key}{"DP"});
					my @FILTER = split(/\|/, $variants{$key}{"FILTER"});
					my @QUAL = split(/\|/, $variants{$key}{"QUAL"});
					my @VF = split(/\|/, $variants{$key}{"VF"});
					my @GT = split(/\|/, $variants{$key}{"GT"});
					my @GQ = split(/\|/, $variants{$key}{"GQ"});
					my @GI = split(/,/, $variants{$key}{"GI"});
					for (my $i=1; $i<scalar(@GI); $i++) {
						if ($GI[0] =~ m/-AS/ and defined$GI[$i]) {
							$GI[0] = $GI[$i];
						}
						else {
							last;
						}
					}

					my $fail = 0;
					for (my $i=0; $i<scalar(@SB); $i++) {
						if ($SB[$i] > -100 ) {
							$fail = 1;
							last;
						}
					}
					if ($fail == 1) {
						if ($variants{$key}{'CaseMatched'} == 0) {
							next;
						}
						else {
							$NGS_variant{$sample}{$variants{$key}{target}}{$variants{$key}{var}} = 5;
						}
					}

					$fail = 0;
					for (my $i=0; $i<scalar(@DP); $i++) {
						if ($DP[$i] < 400 ) {
							$fail = 1;
							last;
						}
					}
					if ($fail == 1) {
						if ($variants{$key}{'CaseMatched'} == 0) {
							next;
						}
						else {
							$NGS_variant{$sample}{$variants{$key}{target}}{$variants{$key}{var}} = 6;
						}
					}

					$fail = 0;
					for (my $i=0; $i<scalar(@FILTER); $i++) {
						if ($FILTER[$i] ne "PASS" ) {
							$fail = 1;
							last;
						}
					}
					if ($fail == 1) {
						if ($variants{$key}{'CaseMatched'} == 0) {
							next;
						}
						else {
							$NGS_variant{$sample}{$variants{$key}{target}}{$variants{$key}{var}} = 7;
						}
					}

					$fail = 0;
					for (my $i=0; $i<scalar(@QUAL); $i++) {
						if ($QUAL[$i] < 100 ) {
							$fail = 1;
							last;
						}
					}
					if ($fail == 1) {
						if ($variants{$key}{'CaseMatched'} == 0) {
							next;
						}
						else {
							$NGS_variant{$sample}{$variants{$key}{target}}{$variants{$key}{var}} = 8;
						}
					}

					$fail = 0;
					for (my $i=0; $i<scalar(@VF); $i++) {
						if ($VF[$i] < 0.05 ) {
							$fail = 1;
							last;
						}
					}
					if ($fail == 1) {
						if ($variants{$key}{'CaseMatched'} == 0) {
							next;
						}
						else {
							$NGS_variant{$sample}{$variants{$key}{target}}{$variants{$key}{var}} = 9;
						}
					}

					if ($samples{$sample}{"run_count"} == 1) {
						$color_flag = 43;
					}
					else {
						$color_flag = 22;
					}
					$format->set_bg_color($color_flag);
#					$format->set_format_properties(bg_color => $color_flag);

					$row ++;
					$worksheet->write("A$row", $GI[0], $format);
					$worksheet->write("B$row", $variants{$key}{'REF'}.">".$variants{$key}{'ALT'}, $format);
					$worksheet->write("C$row", $variants{$key}{'CHROM'}, $format);
					$worksheet->write("D$row", $variants{$key}{'POS'}, $format);

					my $Type = "";
					if ($variants{$key}{'HGVS.c'} =~ m/delins/) {
						$Type = "indel";
					}
					elsif ($variants{$key}{'HGVS.c'} =~ m/ins/) {
						$Type = "insertion";
					}
					elsif ($variants{$key}{'HGVS.c'} =~ m/del/) {
						$Type = "deletion";
					}
					elsif ($variants{$key}{'HGVS.c'} =~ m/dup/) {
						$Type = "duplication";
					}
					elsif ($variants{$key}{'HGVS.c'} =~ m/inv/) {
						$Type = "inversion";
					}
					elsif ($variants{$key}{'HGVS.c'} =~ m/>/) {
						$Type = "SNP";
					}
					$worksheet->write("E$row", $Type, $format);

					$worksheet->write("F$row", $GT[$n], $format);

					my $exon = "no";
					if ($variants{$key}{'EXON'} eq 'EXON') {
						$exon = "yes";
					}
					$worksheet->write("G$row", $exon, $format);

					$worksheet->write("H$row", $FILTER[$n], $format);
					$worksheet->write("I$row", $QUAL[$n], $format);
					$worksheet->write("J$row", $GQ[$n], $format);
					$worksheet->write("K$row", $VF[$n], $format);
					$worksheet->write("L$row", $DP[$n], $format);

					my @Allelic_Depth = split(/\|/, $variants{$key}{'AD'});
					my @Alt_Depth = split(/\,/, $Allelic_Depth[$n]);
					$worksheet->write("M$row", $Alt_Depth[1], $format);

					$worksheet->write("N$row", $Allelic_Depth[$n], $format);
					$worksheet->write("O$row", $variants{$key}{'ID'}, $format);
					$worksheet->write("P$row", $variants{$key}{'Annotation'}, $format);

					my @Feature_ID = split(/\|/, $variants{$key}{'Feature_ID'});
					my @HGVSc = split(/\|/, $variants{$key}{'HGVS.c'});
					my $column = "";
					for (my $i=0; $i<scalar(@HGVSc); $i++) {
						$HGVSc[$i] = "" if (!defined $HGVSc[$i]);
						$HGVSc[$i] = $Feature_ID[$i].":".$HGVSc[$i] if (defined $Feature_ID[$i]);
						$column = $column."|" if ($i>0);
						$column = $column.$HGVSc[$i];
					}
					$worksheet->write("Q$row", $column, $format);

					$worksheet->write("R$row", $variants{$key}{'HGVS.p'}, $format);
					if ($variants{$key}{'CaseMatched'} == 0) {
						$worksheet->write("S$row", "no", $format);
					}
					else {
						$worksheet->write("S$row", "yes", $format);
					}
				}

				$workbook->close() or die "Error closing file: $!";



				##### output tsv

				if (defined $NGS_sample{$sample}) {
					my $output_file = "$out_tsv/$run_id[0]/$samples{$sample}{$n}"."_".$run_id[0].".tsv";
					open(OUTPUT_FILE, "> $output_file") || die "ERROR: Cannot open $output_file for writing\n".$!;

					print OUTPUT_FILE "Gene\tVariant\tChr\tCoordinate\tType\tGenotype\tExonic\tFilters\tQuality\tGQX\tVariant Freq\tRead Depth\tAlt Read Depth\tAllelic Depth\tdbSNP ID\tConsequence\tHGVSc\tHGVSp\tReported\n";

					for my $key (sort {$variants{$b}{ConcordanceCount}<=>$variants{$a}{ConcordanceCount}} keys %variants){
						if (!defined $variants{$key}) {
							next;
						}
						if ($variants{$key}{'CaseMatched'}==0 && $variants{$key}{'freq'}/$total_samples > 0.8) {
							next;
						}
						if ($variants{$key}{'CaseMatched'}==0 && defined $Artifact{$key}) {
							next;
						}
						if ($variants{$key}{'CaseMatched'}==0 && $variants{$key}{"ConcordanceCount"} < $samples{$sample}{"run_count"} && $samples{$sample}{"run_count"} > 1) {
							next;
						}
						my @SB = split(/\|/, $variants{$key}{"SB"});
						my @DP = split(/\|/, $variants{$key}{"DP"});
						my @FILTER = split(/\|/, $variants{$key}{"FILTER"});
						my @QUAL = split(/\|/, $variants{$key}{"QUAL"});
						my @VF = split(/\|/, $variants{$key}{"VF"});
						my @GT = split(/\|/, $variants{$key}{"GT"});
						my @GQ = split(/\|/, $variants{$key}{"GQ"});
						my @GI = split(/,/, $variants{$key}{"GI"});
						for (my $i=1; $i<scalar(@GI); $i++) {
							if ($GI[0] =~ m/-AS/ and defined$GI[$i]) {
								$GI[0] = $GI[$i];
							}
							else {
								last;
							}
						}
						if ($variants{$key}{'CaseMatched'} == 0) {
							my $fail = 0;
							for (my $i=0; $i<scalar(@SB); $i++) {
								if ($SB[$i] > -100 ) {
									$fail = 1;
									last;
								}
							}
							if ($fail == 1) {
								next;
							}

							$fail = 0;
							for (my $i=0; $i<scalar(@DP); $i++) {
								if ($DP[$i] < 400 ) {
									$fail = 1;
									last;
								}
							}
							if ($fail == 1) {
								next;
							}

							$fail = 0;
							for (my $i=0; $i<scalar(@FILTER); $i++) {
								if ($FILTER[$i] ne "PASS" ) {
									$fail = 1;
									last;
								}
							}
							if ($fail == 1) {
								next;
							}

							$fail = 0;
							for (my $i=0; $i<scalar(@QUAL); $i++) {
								if ($QUAL[$i] < 100 ) {
									$fail = 1;
									last;
								}
							}
							if ($fail == 1) {
								next;
							}

							$fail = 0;
							for (my $i=0; $i<scalar(@VF); $i++) {
								if ($VF[$i] < 0.05 ) {
									$fail = 1;
									last;
								}
							}
							if ($fail == 1) {
								next;
							}
						}

						print OUTPUT_FILE $GI[0];
						print OUTPUT_FILE "\t".$variants{$key}{'REF'}.">".$variants{$key}{'ALT'};
						print OUTPUT_FILE "\t".$variants{$key}{'CHROM'};
						print OUTPUT_FILE "\t".$variants{$key}{'POS'};

						my $Type = "";
						if ($variants{$key}{'HGVS.c'} =~ m/delins/) {
							$Type = "indel";
						}
						elsif ($variants{$key}{'HGVS.c'} =~ m/ins/) {
							$Type = "insertion";
						}
						elsif ($variants{$key}{'HGVS.c'} =~ m/del/) {
							$Type = "deletion";
						}
						elsif ($variants{$key}{'HGVS.c'} =~ m/dup/) {
							$Type = "duplication";
						}
						elsif ($variants{$key}{'HGVS.c'} =~ m/inv/) {
							$Type = "inversion";
						}
						elsif ($variants{$key}{'HGVS.c'} =~ m/>/) {
							$Type = "SNP";
						}
						print OUTPUT_FILE "\t".$Type;

						print OUTPUT_FILE "\t".$GT[$n];

						my $exon = "no";
						if ($variants{$key}{'EXON'} eq 'EXON') {
							$exon = "yes";
						}
						print OUTPUT_FILE "\t".$exon;

						print OUTPUT_FILE "\t".$FILTER[$n];
						print OUTPUT_FILE "\t".$QUAL[$n];
						print OUTPUT_FILE "\t".$GQ[$n];
						print OUTPUT_FILE "\t".$VF[$n];
						print OUTPUT_FILE "\t".$DP[$n];

						my @Allelic_Depth = split(/\|/, $variants{$key}{'AD'});
						my @Alt_Depth = split(/\,/, $Allelic_Depth[$n]);
						print OUTPUT_FILE "\t".$Alt_Depth[1];

						print OUTPUT_FILE "\t".$Allelic_Depth[$n];
						print OUTPUT_FILE "\t".$variants{$key}{'ID'};
						print OUTPUT_FILE "\t".$variants{$key}{'Annotation'};
						print OUTPUT_FILE "\t".$variants{$key}{'Feature_ID'}.":".$variants{$key}{'HGVS.c'};
						print OUTPUT_FILE "\t".$variants{$key}{'HGVS.p'};
						if ($variants{$key}{'CaseMatched'} == 0) {
							print OUTPUT_FILE "\tno";
						}
						else {
							print OUTPUT_FILE "\tyes";
						}
						print OUTPUT_FILE "\n";
					}
					close OUTPUT_FILE;


					##### convert tsv to vcf

					open(TS, "< $output_file") || die "ERROR: Cannot open $output_file for reading\n$!";

					$output_file = "$out_vcf/$run_id[0]/$samples{$sample}{$n}"."_".$run_id[0].".vcf";
					open(OUT, "> $output_file") || die "ERROR: Cannot open $output_file for writing\n$!";

					my %tsv_header_hash = ();
					my $header_flag = 0;

					while( my $tsv_line = <TS> ){
						chomp $tsv_line;
						$tsv_line =~ s/\"//g;
						## Check if the line is empty
						if($tsv_line =~ m/.*\S.*/){
							## Check if file has header information and print OUT appropriate headers in vcf file
							if($header_flag == 0){
								$header_flag = 1;
								my @temp_data=split("\t",$tsv_line);
								for(my $i=0; $i<scalar(@temp_data);$i++){
									$tsv_header_hash{$temp_data[$i]}=$i;
								}
								my $out_str = &HeadersOfVCF;
								print OUT $out_str, $sample."_".$run_id[0], "\n";
								next;
							}
						}
						else {
							next;
						}

						if(!defined $tsv_header_hash{"Gene"}){
							print OUT "ERROR\: Header absent in file. Please input the correct file and try again.\n";
							print "ERROR\: Header absent in file. Please input the correct file and try again.\n";
							last;
						}
						else {
							my @temp_data = split("\t",$tsv_line);

							# 00_CHROM and 01_POS
							print OUT $temp_data[$tsv_header_hash{"Chr"}],"\t",$temp_data[$tsv_header_hash{"Coordinate"}],"\t";

							# 02_ID
							if($temp_data[$tsv_header_hash{"dbSNP ID"}] =~ m/^\S+/){	### one or more non-whitespace character
								print OUT $temp_data[$tsv_header_hash{"dbSNP ID"}],"\t";
							}
							else {
								print OUT "\.\t";
							}

							# 03_REF and 04_ALT
							if($temp_data[$tsv_header_hash{"Variant"}] =~ m/(\w+)\>(\w+)/){
								print OUT $1,"\t",$2;
							}
							else {
								print OUT "ERROR\:There seems to be an error with the standard variant nomenclature. Please check and try again.\n";
								print "ERROR\:There seems to be an error with the standard variant nomenclature. Please check and try again.\n";
								last;
							}

							# 05_QUAL
							my @fields = split(/\|/, $temp_data[$tsv_header_hash{"Quality"}]);
							print OUT "\t", $fields[0];

							# 06_FILTER
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"Filters"}]);
							print OUT "\t", $fields[0], "\t";

							# 07_INFO: DP, TI, GI, and FC
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"Read Depth"}]);
							print OUT "DP\=",$fields[0];
							@fields = split(/\:/, $temp_data[$tsv_header_hash{"HGVSc"}]);
							print OUT "\;TI\=",$fields[0];
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"Gene"}]);
							print OUT "\;GI\=",$fields[0];
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"Consequence"}]);
							print OUT "\;FC\=",$fields[0];

							# 07_INFO: HGVS.c and HGVS.p
							@fields = split(/\:/, $temp_data[$tsv_header_hash{"HGVSc"}]);
							print OUT "\;HGC\=",$fields[1];
							if (defined $temp_data[$tsv_header_hash{"HGVSp"}]) {
								if ($temp_data[$tsv_header_hash{"HGVSp"}] ne "-") {
									print OUT "\;HGP\=",$temp_data[$tsv_header_hash{"HGVSp"}];
								}
							}

							# 07_INFO: REPORT
							print OUT "\;REPORT\=",$temp_data[$tsv_header_hash{"Reported"}];

							# 08_FORMAT:
							print OUT "\tGT\:GQ\:AD\:VF\:GQX\t";

							#09_SampleINFO
							#0/1:100:2022,65:0.0310:20:-100.0000:100
							#GT
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"Genotype"}]);
							$fields[0] =~ s/het/0\/1/g;
							$fields[0] =~ s/hom/1\/1/g;
							print OUT $fields[0];
							#GQ
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"Quality"}]);
							print OUT "\:",$fields[0];
							#AD
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"Allelic Depth"}]);
							print OUT "\:",$fields[0];
							#VF (This is: Alt Variant Freq_13/100)
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"Variant Freq"}]);
							print OUT "\:",$fields[0];
							#GQX
							@fields = split(/\|/, $temp_data[$tsv_header_hash{"GQX"}]);
							print OUT "\:",$fields[0];
							print OUT "\n";
						}
					}
					close TS;
					close OUT;
				}
			}


			##### output variants missed in NGS_variants

			my $total_NGS_variant = 0;
			my $row = $row2;
			if (defined $NGS_variant{$sample}) {
				for my $k1 (keys %{$NGS_variant{$sample}}) {
					for my $k2 (keys %{$NGS_variant{$sample}{$k1}}) {
						$total_NGS_variant++;
						if ($NGS_variant{$sample}{$k1}{$k2} == 0) {
							$format_m->set_bg_color(22);
							$worksheet_m->write("E$row", "missed", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 1) {
							$format_m->set_bg_color(50);
							$worksheet_m->write("E$row", "found", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 2) {
							$format_m->set_bg_color(66);
							$worksheet_m->write("E$row", "found but freq>0.8", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 3) {
							$format_m->set_bg_color(66);
							$worksheet_m->write("E$row", "found but artifact-like", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 4) {
							$format_m->set_bg_color(66);
							$worksheet_m->write("E$row", "found but discordant", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 5) {
							$format_m->set_bg_color(66);
							$worksheet_m->write("E$row", "found but SB>-100", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 6) {
							$format_m->set_bg_color(66);
							$worksheet_m->write("E$row", "found but DP<400", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 7) {
							$format_m->set_bg_color(66);
							$worksheet_m->write("E$row", "found but not PASS", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 8) {
							$format_m->set_bg_color(66);
							$worksheet_m->write("E$row", "found but QUAL<100", $format_m);
						}
						elsif ($NGS_variant{$sample}{$k1}{$k2} == 9) {
							$format_m->set_bg_color(66);
							$worksheet_m->write("E$row", "found but VF<0.05", $format_m);
						}
						$worksheet_m->write("A$row", $sample, $format_m);
						$worksheet_m->write("B$row", $k1, $format_m);
						$worksheet_m->write("C$row", $k2, $format_m);
						$worksheet_m->write("D$row", $run_id[0], $format_m);
						$row++;
					}
				}
			}
			$row2 = $row;

		}
		last if ($run_found == 1);
	}
	closedir(RUN_DIR);
	last if ($run_found == 1);
}
close OUT_SAMPLES;

$row2 = 2;
for my $key (sort {$NGS_sample{$a}<=>$NGS_sample{$b}} keys %NGS_sample){
	if ($NGS_sample{$key} == 0) {
		$format_s->set_bg_color(22);
		$worksheet_s->write("B$row2", "missed", $format_s);
	}
	else {
		$format_s->set_bg_color(50);
		$worksheet_s->write("B$row2", "found", $format_s);
	}
	$worksheet_s->write("A$row2", $key, $format_s);
	$row2++;
}


$workbook_m->close() or die "Error closing file: $!";
$workbook_s->close() or die "Error closing file: $!";

exit 0;

########### SUBROUTINES #############
sub HeadersOfVCF {
	my $ret_string = "\#\#fileformat\=VCFv4\.1\n";
	$ret_string .= "\#\#FORMAT\=\<ID\=GQX\,Number\=1\,Type\=Integer\,Description\=\"Minimum of \{Genotype quality assuming variant position\,Genotype quality assuming non\-variant position\}\"\>\n";
	$ret_string .= "\#\#FORMAT\=\<ID\=GT\,Number\=1\,Type\=String\,Description\=\"Genotype\"\>\n";
	$ret_string .= "\#\#FORMAT\=\<ID\=GQ\,Number\=1\,Type\=Integer\,Description\=\"Genotype Quality\"\>\n";
	$ret_string .= "\#\#FORMAT\=\<ID\=AD\,Number\=\.\,Type\=Integer\,Description\=\"Allele Depth\"\>\n";
	$ret_string .= "\#\#FORMAT\=\<ID\=VF\,Number\=1\,Type\=Float\,Description\=\"Variant Frequency\"\>\n";
#	$ret_string .= "\#\#FORMAT\=\<ID\=NL\,Number\=1\,Type\=Integer\,Description\=\"Applied BaseCall Noise Level\"\>\n";
#	$ret_string .= "\#\#FORMAT\=\<ID\=SB\,Number\=1\,Type\=Float\,Description\=\"StrandBias Score\"\>\n";
	$ret_string .= "\#\#INFO\=\<ID\=TI\,Number\=\.\,Type\=String\,Description\=\"Transcript ID\"\>\n";
	$ret_string .= "\#\#INFO\=\<ID\=GI\,Number\=\.\,Type\=String\,Description\=\"Gene ID\"\>\n";
	$ret_string .= "\#\#INFO\=\<ID\=EXON\,Number\=0\,Type\=Flag\,Description\=\"Exon Region\"\>\n";
	$ret_string .= "\#\#INFO\=\<ID\=REPORT\,Number\=0\,Type\=Flag\,Description\=\"Reported to SCI\"\>\n";
	$ret_string .= "\#\#INFO\=\<ID\=FC\,Number\=\.\,Type\=String\,Description\=\"Functional Consequence\"\>\n";
	$ret_string .= "\#\#INFO\=\<ID\=DP\,Number\=1\,Type\=Integer\,Description\=\"Total Depth\"\>\n";
	$ret_string .= "\#\#INFO\=\<ID\=HGC\,Number\=1\,Type\=String\,Description\=\"HGVS c.\"\>\n";
	$ret_string .= "\#\#INFO\=\<ID\=HGP\,Number\=1\,Type\=String\,Description\=\"HGVS p.\"\>\n";
	$ret_string .= "\#\#FILTER\=\<ID\=LowVariantFreq\,Description\=\"Low variant frequency \< 0\.01\"\>\n";
	$ret_string .= "\#\#FILTER\=\<ID\=LowGQ\,Description\=\"GQ below \< 30\.00\"\>\n";
	$ret_string .= "\#\#FILTER\=\<ID\=R8\,Description\=\"IndelRepeatLength is greater than 8\"\>\n";
	$ret_string .= "\#\#FILTER\=\<ID\=SB\,Description\=\"Variant strand bias too high\"\>\n";
	$ret_string .= "\#\#source\=CallSomaticVariantsv2\.1\.12\.0\n";
	$ret_string .= "\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	return $ret_string;
}




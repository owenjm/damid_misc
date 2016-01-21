#!/usr/bin/perl -w

# gff_average.pl -- averages multiple gff files
# Copyright © 2015 Owen Marshall

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version. 
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 
# USA

use strict;
use File::Basename;
$|++;

my @scores;
my @logscores;
my @probes;
my @files;
my @out;

my %probes;
my %scores;

my $format = 'bedgraph';


unless (@ARGV) {
	help();
}

my %vars = (
	'strip_chrs' => 0,
);

my %vars_details = (
	'strip_chrs' => 'Remove "chr" prefixes from chromosome numbers',
);

my @in_files;

process_cli();

my $num_samples = @in_files;

foreach (@in_files) {
	load_gff($_);
}

my ($fhead,$dir,$ext) = fileparse($ARGV[0], qr/\.[^.]*/);
$fhead =~ s/-n\d+-/-/; # remove -n1- -n2- type annotations
my $fout = "$fhead.average.$format";

analyse();
sort_array();
output();

print STDERR "All done.\n\n";




sub analyse {
	print STDERR "Averaging ...                     \n";
	my $count=0;
	foreach my $k (keys %probes) {
		
		my @ld = @{ $scores{$k} };
		
		if (@ld < $num_samples) {
			print STDERR "Warning: unmatched probe at $k ... \n";
			next;
		}
		
		my ($chr, $start, $end) = @{ $probes{$k} };
		my ($mean, $sd, $se) = stdev(@ld);
		my $fano = fano($mean, $sd);
		$count++ if $fano > 1;
		
		push @out,  [$chr, '.', '.', $start, $end, $mean, '.', '.', '.'];
		
	}
	my $t = $#out;
	print "  $t probes processed\n  $count with fano > 1\n";
}

sub sort_array {
	print STDERR "Sorting ...\n";
	@out = sort {$a->[3] <=> $b->[3]} @out;
	@out = sort {$a->[0] cmp $b->[0]} @out;
}

sub output {
	print STDERR "Writing output to $fout ...\n";
	open (OUT, ">$fout") || die "Cannot open output file for writing: $!\n";
	
	print OUT qq(track type=bedGraph name="$fhead.average" description="$fhead.average"\n) unless $format eq 'gff';
	
	foreach (@out) {
		my ($chr, $source, $type, $start, $end, $score, $b, $c, $name) = @{$_};
		if ($format eq 'gff') {
			print OUT join("\t", $chr, '.', '.', $start, $end, $score, '.', '.', '.'), "\n";
		} else {
			print OUT join("\t", $chr, $start, $end, $score), "\n";
		}
	}
}

sub load_gff {
	my ($fn, $ar2) = @_;
	push @files, $fn;
	print STDERR "Reading input file: $fn ...\n";
	open (IN, "<$fn") || die "Unable to read $fn: $!\n";
	
	my $p=0;
	foreach (<IN>) {
		next if m/^track/;
		chomp;
		
		my @line = split('\t');
		
		next unless $#line > 2;
		
		my ($chr, $start, $end, $score);
		if ($#line == 3) {
			# bedgraph
			($chr, $start, $end, $score) = @line;
		} else {
			# GFF
			($chr, $start, $end, $score) = @line[0,3,4,5];
			$format = 'gff';
		}
		
		next unless $start;
		
		print STDERR "Read $p lines ...\r" if ($p%10000 == 0);
		
		$chr=~s/chr//g if $vars{'strip_chrs'};
		
		$probes{"$chr-$start"} = [$chr, $start, $end];
		push @{ $scores{"$chr-$start"} }, $score;
		
		$p++;
	}
	close (IN);
}

sub stdev {
	my @a = @_;
	my $n = @a;
	return ($a[0], 0, 0, $a[0]) if $n<2;
	
	my $total;
	foreach (@a) {
		$total+=$_;
	}
	my $mean = $total/$n;
	
	my $sum;
	foreach my $x (@a) {
		$sum+=abs($x-$mean);
	}
	my $sd= $sum/($n-1);
	my $se = $sd/sqrt($n);
	
	return ($mean, $sd, $se);
}

sub fano {
	# the number of probes with a fano factor > 1 is a quick indication of the overall correlation
	my ($mean, $sd) = @_;
	my $variance = $sd**2;
	$mean = 1 unless $mean;
	my $fano = $variance/$mean;
	return $fano;
}


sub process_cli {
	foreach (@ARGV) {
		if (/--(.*)=(.*)/) {
			unless (defined($vars{$1})) {
				print STDERR "Did not understand $_ ...\n";
				help();
			}
			my ($v, $opt) = ($1,$2);
			$vars{$v} = $opt;
			next;
		} elsif (/--h[elp]*/) {
			help();
		} elsif (/--(.*)/) {
			# if no parameter is specified we assume it's a switch ...
			if (defined($vars{$1})) {
				$vars{$1} = 1;
			} else {
				print STDERR "Did not understand $_ ...\n";
				help();
			}
			next;
		}
		push @in_files, $_;
	}
}


sub help {
	print STDOUT <<EOT;

Usage: gff_average [gff files to process]

Options:
EOT
	
	my $opt_len = 0;
	foreach (keys %vars) {
		my $l = length($_);
		$opt_len = $l if $l > $opt_len;
	}
	
	$opt_len+=2;
	
	my $cols= `tput cols` || 80;
	
	my ($v, $val, $def, $def_format);
	my $help_format = "format STDOUT =\n"
		.' '.'^'.'<'x$opt_len . ' '. '^' . '<'x($cols-$opt_len-4) . "\n"
		.'$v, $def_format'."\n"
		.' '.'^'.'<'x$opt_len . '   '. '^' . '<'x($cols-$opt_len-6) . "~~\n"
		.'$v, $def_format'."\n"
		.".\n";
		
	eval $help_format;
	die $@ if $@;
	
	foreach my $k (sort (keys %vars)) {
		($v, $val, $def) = ($k, $vars{$k}, $vars_details{$k});
		$def||="";
		$def_format = $val ? "$def\n\r[Current value: $val]" : $def;
		$v = "--$v";
#		format =
# ^<<<<<<<<<<<<<<<<<<<< ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#$v, $def_format
# ^<<<<<<<<<<<<<<<<<<<<   ^<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ~~
#$v, $def_format
#.

		write();
		
	}
	print STDOUT "\n";
	exit 1;
}
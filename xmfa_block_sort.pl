#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $xmfa_file;
my $print_seq_ids;
my @sort_order = ();
my $fasta_outdir = '.';
my $fasta_postfix = '.sort.fa';
my $fasta_linker;
my $coords_file;
my $sorted_xmfa_file;
my $null_record = 'NA';
my $help;

my %block_seq_order = ();
my %block_seqs = ();
my %block_lens = ();
my %ref_block_starts = ();
my %ref_block_stops = ();
my %block_ref_pos = ();
my @block_order = ();
my %unplaced_blocks = ();
my %seq_id_names = ();
my %seq_id_fasta_fhs = ();

parse_args();
parse_xmfa_header();
my ($coords_fh, $sorted_xmfa_fh) = open_fhs();
parse_blocks();
sort_blocks();
print_blocks();

exit(0);


sub sort_blocks {
	foreach my $block_id (keys %block_seqs) {
		$unplaced_blocks{$block_id}++;
	}

	foreach my $seq_id (@sort_order) {
		if ($seq_id == $sort_order[0]) {
			foreach my $block_start (sort { $a <=> $b } keys %{$ref_block_starts{$seq_id}}) {
				my $block_id = $ref_block_starts{$seq_id}{$block_start};

				if (! exists($unplaced_blocks{$block_id})) {
					next();
				}

				if (! @block_order) {
					push(@block_order, $block_id);
					delete($unplaced_blocks{$block_id});

					next();
				}

				my $last_seq_id_index;

				foreach my $index (0..$#block_order) {
					my $index_block_id = $block_order[$index];

					if (! exists($block_ref_pos{$index_block_id}{$seq_id}{'start'})) {
						next();
					}

					my $index_seq_start = $block_ref_pos{$index_block_id}{$seq_id}{'start'};

					$last_seq_id_index = $index;

					if ($block_start < $index_seq_start) {
						splice(@block_order, $index, 0, $block_id);
						delete($unplaced_blocks{$block_id});

						last();
					}
				}

				if (exists($unplaced_blocks{$block_id}) && defined($last_seq_id_index)) {
					splice(@block_order, $last_seq_id_index + 1, 0, $block_id);
					delete($unplaced_blocks{$block_id});

					next();
				}
			}
		}

		else {
			my @static_block_order = @block_order;

			foreach my $block_id (@static_block_order) {
				my $block_order_index = get_block_order_index($block_id);

				if (defined($block_order_index)) {
					extend_block($block_order_index, $block_id, $seq_id, 'left');
				}

				$block_order_index = get_block_order_index($block_id);

				if (defined($block_order_index)) {
					extend_block($block_order_index, $block_id, $seq_id, 'right');
				}
			}
		}
	}

	return(0);
}


sub print_blocks {
	my %fasta_seq = ();
	my %fasta_pos = ();

	my $block_count = 1;

	if (defined($coords_file)) {
		print($coords_fh join("\t", '#sort_block_pos', 'orig_block_pos', 'block_len'));

		foreach my $seq_id (@sort_order) {
			my $seq_name = $seq_id_names{$seq_id};

			print($coords_fh join("\t", '', "$seq_name.fasta.start", "$seq_name.fasta.stop", "$seq_name.xmfa.start", "$seq_name.xmfa.stop", "$seq_name.xmfa.strand"));
		}

		print($coords_fh "\n");
	}

	foreach my $block_id (@block_order) {
		my $block_len = $block_lens{$block_id};

		if (defined($coords_file)) {
			print($coords_fh "$block_count\t$block_id\t$block_len");
		}

		foreach my $seq_id (@sort_order) {
			my $xmfa_start = $null_record;
			my $xmfa_stop = $null_record;
			my $xmfa_strand = $null_record;
			my $seq;

			if (exists($block_seqs{$block_id}{$seq_id})) {
				$xmfa_start = $block_ref_pos{$block_id}{$seq_id}{'start'};
				$xmfa_stop = $block_ref_pos{$block_id}{$seq_id}{'stop'};
				$xmfa_strand = $block_ref_pos{$block_id}{$seq_id}{'strand'};

				if ($xmfa_stop < $xmfa_start) {
					my $tmp = $xmfa_start;
					$xmfa_start = $xmfa_stop;
					$xmfa_stop = $tmp;
				}

				$seq = $block_seqs{$block_id}{$seq_id};
				$seq =~ s/\n//g;
			}

			else {
				$seq = '-' x $block_len;
			}

			$fasta_seq{$seq_id} .= "$seq";

			if (! exists($fasta_pos{$seq_id})) {
				$fasta_pos{$seq_id} = 0;
			}

			my $fasta_start = $fasta_pos{$seq_id} + 1;
			my $fasta_stop = $fasta_start + $block_len - 1;

			$fasta_pos{$seq_id} = $fasta_stop;

			if (defined($fasta_linker)) {
				$fasta_seq{$seq_id} .= "$fasta_linker";
				$fasta_pos{$seq_id} += length($fasta_linker);
			}

			if (defined($coords_file)) {
				print($coords_fh join("\t", '', $fasta_start, $fasta_stop, $xmfa_start, $xmfa_stop, $xmfa_strand));
			}
		}

		if (defined($sorted_xmfa_file)) {
			foreach my $block_seq_index (sort { $a <=> $b } keys %{$block_seq_order{$block_id}}) {
				my $seq_id = $block_seq_order{$block_id}{$block_seq_index}{'seq_id'};
				my $header = $block_seq_order{$block_id}{$block_seq_index}{'header'};
				my $seq = $block_seqs{$block_id}{$seq_id};

				print($sorted_xmfa_fh "$header\n"."$seq");
			}

			print($sorted_xmfa_fh "=\n");
		}

		if (defined($coords_file)) {
			print($coords_fh "\n");
		}

		$block_count++;
	}

	foreach my $seq_id (keys %seq_id_fasta_fhs) {
		my $fasta_fh = $seq_id_fasta_fhs{$seq_id};
		my $wrapped_seq = $fasta_seq{$seq_id};

		$wrapped_seq =~ s/.{72}\K/\n/g;

		print($fasta_fh "$wrapped_seq");

		close($fasta_fh);
	}


	return(0);
}


sub proc_block_seq {
	my $fasta_seq_ref = shift();
	my $fasta_pos_ref = shift();
	my $seq_id = shift();
	my $start = shift();
	my $stop = shift();
	my $strand = shift();

	my $seq_len = length(abs($start - $stop)) + 1;

	if (! defined($$fasta_pos_ref{$seq_id})) {
		$$fasta_pos_ref{$seq_id} = 0;
	}


	if (defined($fasta_linker)) {
		$$fasta_seq_ref{$seq_id} .= "$fasta_linker";
	}

	return(0);
}


sub extend_block {
	my $block_order_index = shift();
	my $block_id = shift();
	my $seq_id = shift();
	my $dir = shift();

	if (! defined($block_order_index) || ! defined($block_id) || ! defined($seq_id) || ! defined($dir)) {
		return(0);
	}

	my $start = $block_ref_pos{$block_id}{$seq_id}{'start'};
	my $stop = $block_ref_pos{$block_id}{$seq_id}{'stop'};

	if (! defined($start) || ! defined($stop)) {
		return(0);
	}

	foreach my $unplaced_block_id (keys %unplaced_blocks) {
		my $unplaced_start = $block_ref_pos{$unplaced_block_id}{$seq_id}{'start'};
		my $unplaced_stop = $block_ref_pos{$unplaced_block_id}{$seq_id}{'stop'};

		if (! defined($unplaced_start) || ! defined($unplaced_stop)) {
			next();
		}

		my $match = 0;

		if ($dir eq 'left') {
			if ($unplaced_stop == $start - 1) {
				splice(@block_order, $block_order_index, 0, $unplaced_block_id);
				delete($unplaced_blocks{$unplaced_block_id});
				$match = 1;

				extend_block($block_order_index, $unplaced_block_id, $seq_id, $dir);
			}
		}

		elsif ($dir eq 'right') {
			if ($unplaced_start == $stop + 1) {
				splice(@block_order, $block_order_index + 1, 0, $unplaced_block_id);
				delete($unplaced_blocks{$unplaced_block_id});
				$match = 1;

				extend_block($block_order_index + 1, $unplaced_block_id, $seq_id, $dir);
			}
		}

		if ($match == 1) {
			last();
		}
	}

	return(0);
}


sub get_block_order_index {
	my $block_id = shift();

	foreach my $index (0..$#block_order) {
		if ($block_order[$index] == $block_id) {
			return($index);
		}
	}

	return(undef);
}


sub parse_blocks {
	open(XMFA, '<', $xmfa_file) or error("$!");

	my $block_id = 1;
	my $valid_seq = 0;
	my $block_seq_count = 0;
	my $seq_id;
	my $pos;
	my $seq;

	while (my $line = <XMFA>) {
		chomp($line);

		if ($line =~ /^#/) {
			if (defined($sorted_xmfa_file)) {
				print($sorted_xmfa_fh "$line\n");
			}
		}

		elsif ($line =~ /^\>/) {
			my $mod_line = $line;

			$mod_line =~ s/^>\s*//;

			my ($ref_pos, $strand, $ref_seq_file) = split(/\s+/, $mod_line);
			($seq_id, $pos) = split(':', $ref_pos);
			my ($start, $stop) = split('-', $pos);

			if (defined($start) && defined($stop) && $start > $stop) {
				my $tmp = $start;
				$start = $stop;
				$stop = $start;
			}

			if (defined($seq_id) && defined($start) && defined($stop) && defined($strand)) {
				$block_seq_count++;

				$ref_block_starts{$seq_id}{$start} = $block_id;
				$ref_block_stops{$seq_id}{$stop} = $block_id;
				$block_ref_pos{$block_id}{$seq_id}{'start'} = $start;
				$block_ref_pos{$block_id}{$seq_id}{'stop'} = $stop;
				$block_ref_pos{$block_id}{$seq_id}{'strand'} = $strand;
				$block_seq_order{$block_id}{$block_seq_count}{'seq_id'} = $seq_id;
				$block_seq_order{$block_id}{$block_seq_count}{'header'} = $line;

				$valid_seq = 1;
			}

			else {
				$valid_seq = 0;
			}

			$seq = '';
		}

		elsif ($line =~ /^\=/) {
			$block_seq_count = 0;
			$block_lens{$block_id} = length($seq);
			$block_id++;
		}

		else {
			if ($valid_seq == 1) {
				$seq .= "$line";
				$block_seqs{$block_id}{$seq_id} .= "$line\n";
			}
		}
	}

	close(XMFA);

	return(0);
}


sub open_fhs {
	my $coords_fh;
	my $sorted_xmfa_fh;

	if (defined($coords_file)) {
		open($coords_fh, '>', $coords_file) or error("$!");
	}

	if (defined($sorted_xmfa_file)) {
		open($sorted_xmfa_fh, '>', $sorted_xmfa_file) or error("$!");
	}

	return($coords_fh, $sorted_xmfa_fh);
}


sub parse_xmfa_header {
	my %header_seq_ids = ();

	open(XMFA, '<', $xmfa_file) or error("$!");

	while (my $line = <XMFA>) {
		if ($line =~ /^#Sequence(\d+)File/) {
			chomp($line);

			my $seq_id = $1;

			$header_seq_ids{$seq_id}++;

			my ($header, $ref_seq_file) = split(/\t/, $line);

			if (defined($print_seq_ids)) {
				print("id: $seq_id\tseq file: $ref_seq_file\n");
			}

			else {
				my $file_base = $ref_seq_file;

				$file_base =~ s/^.*\///;
				$file_base =~ s/\.f[ast]*a//;

				$seq_id_names{$seq_id} = $file_base;

				my $seq_id_fasta_file = "$fasta_outdir/$file_base"."$fasta_postfix";

				open($seq_id_fasta_fhs{$seq_id}, '>', $seq_id_fasta_file) or error("$!");

				my $fasta_fh = $seq_id_fasta_fhs{$seq_id};

				print($fasta_fh ">$file_base"."$fasta_postfix\n");
			}
		}

		elsif ($line =~ /^>/) {
			last();
		}
	}

	close(XMFA);

	foreach my $seq_id (@sort_order) {
		if (exists($header_seq_ids{$seq_id})) {
			delete($header_seq_ids{$seq_id});
		}

		else {
			error("seq id: $seq_id not found in $xmfa_file header\nrun '$0 -x $xmfa_file -p' to see available seq ids");
		}
	}

	foreach my $seq_id (sort { $a <=> $b } keys %header_seq_ids) {
		push(@sort_order, $seq_id);
	}

	if (defined($print_seq_ids)) {
		exit(0);
	}

	return(0);
}


sub error {
	my $msg = shift();

	if (defined($msg)) {
		print("error: $msg\n");
	}

	exit(0);
}


sub arg_error {
	my $msg = shift();

	if (defined($msg)) {
		print("error: $msg\n\n");
	}

	print("use option -h or --help to display help menu\n");

	exit(0);
}


sub parse_args {
	if ($#ARGV == -1) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	GetOptions ('x|xfma=s' => \$xmfa_file,
				'p|print' => \$print_seq_ids,
				's|sort=s{,}' => \@sort_order,
				'o|outdir=s' => \$fasta_outdir,
				'postfix=s' => \$fasta_postfix,
				'l|linker=s' => \$fasta_linker,
				'c|coords=s' => \$coords_file,
				'xmfasort=s' => \$sorted_xmfa_file,
				'n|null=s' => \$null_record,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	if (! defined($xmfa_file)) {
		arg_error('xmfa file required');
	}

	$fasta_outdir =~ s/\/$//;
}


__END__

=head1 NAME

xmfa_block_sort.pl

=head1 SYNOPSIS

xmfa_block_sort.pl -x xmfa.file [options]

=head1 DESCRIPTION

sort XMFA file by block (LCB) and generate gapped fasta files with fasta/xmfa coordinates

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -x --xmfa     input xmfa file (required)

 -p --print    print reference seq ids, then exit

 -s --sort     sort order
                 provide one or more ref seq ids (from header)
                 default: sort by 1, 2, 3... 
                 example: -s 2 3 1 (sort by 2, then 3, then 1)
                 example: -s 3 (sort by 3, then 1, then 2)

 -c --coords   output fasta coords file

 --xmfasort    output sorted xmfa file

 -o --outdir   output fasta directory
                 default: current working directory

 --postfix     output fasta file postfix
                 default: ".sort.fa"

 -l --linker   output fasta block sequence linker
                 default: no linker
                 example: -l "NNNNNNNNNN"

 -n --null     null record value
                 default: "NA"

 -h --help     display help menu

=cut

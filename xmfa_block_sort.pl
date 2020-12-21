#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my $xmfa_file;
my $print_seq_ids;
my @sort_order = ();
my $help;

my %blocks = ();
my %ref_block_starts = ();
my %ref_block_stops = ();
my %block_ref_pos = ();

parse_args();
parse_xmfa_header();
parse_blocks();
sort_blocks();

exit(0);


sub sort_blocks {
	my %proc_blocks = ();
	my @block_order = ();

	foreach my $seq_id (@sort_order) {
		foreach my $block_start (sort { $a <=> $b } keys %{$ref_block_starts{$seq_id}}) {
			my $block_id = $ref_block_starts{$seq_id}{$block_start};

			if (exists($proc_blocks{$block_id})) {
				next();
			}

			if (! @block_order) {
				push(@block_order, $block_id);
				$proc_blocks{$block_id}++;

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
					$proc_blocks{$block_id}++;

					last();
				}
			}

			if (! exists($proc_blocks{$block_id}) && defined($last_seq_id_index)) {
				splice(@block_order, $last_seq_id_index + 1, 0, $block_id);
				$proc_blocks{$block_id}++;

				next();
			}
		}
	}

	foreach my $block_id (@block_order) {
		print("$blocks{$block_id}=\n");
	}

	return(0);
}


sub parse_blocks {
	open(XMFA, '<', $xmfa_file) or error("$!");

	my $block_id = 1;
	my $valid_seq = 0;

	while (my $line = <XMFA>) {
		chomp($line);

		if ($line =~ /^#/) {
			print("$line\n");
		}

		elsif ($line =~ /^\>/) {
			my $mod_line = $line;

			$mod_line =~ s/^>\s*//;

			my ($ref_pos, $strand, $ref_seq_file) = split(/\s+/, $mod_line);
			my ($seq_id, $pos) = split(':', $ref_pos);
			my ($start, $stop) = split('-', $pos);

			if (defined($start) && defined($stop) && $start > $stop) {
				my $tmp = $start;
				$start = $stop;
				$stop = $start;
			}

			if (defined($seq_id) && defined($start) && defined($stop) && defined($strand)) {
				$ref_block_starts{$seq_id}{$start} = $block_id;
				$ref_block_stops{$seq_id}{$stop} = $block_id;
				$block_ref_pos{$block_id}{$seq_id}{'start'} = $start;
				$block_ref_pos{$block_id}{$seq_id}{'stop'} = $stop;
				$blocks{$block_id} .= "$line\n";

				$valid_seq = 1;
			}

			else {
				$valid_seq = 0;
			}
		}

		elsif ($line =~ /^\=/) {
			$block_id++;
		}

		else {
			if ($valid_seq == 1) {
				$blocks{$block_id} .= "$line\n";
			}
		}
	}

	close(XMFA);

	return(0);
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
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	if (! defined($xmfa_file)) {
		arg_error('xmfa file required');
	}

	if (! defined($print_seq_ids)) {
	}
}


__END__

=head1 NAME

xmfa_block_sort.pl

=head1 SYNOPSIS

xmfa_block_sort.pl [options]

=head1 DESCRIPTION

xmfa_block_sort.pl sorts XMFA files by block (LCB)

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

 -x --xmfa     xmfa file (required)

 -p --print    print reference seq ids, then exit

 -s --sort     sort order
                 provide one or more ref seq ids (from header)
                 default: sort by 1, 2, 3... 
                 example: -s 2 3 1 (sort by 2, then 3, then 1)
                 example: -s 3 (sort by 3, then 1, then 2)

 -h --help     display help menu

=cut

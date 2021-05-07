#!/usr/bin/env perl

use Getopt::Long;
use Pod::Usage;

use strict;
use warnings;

my @xmfa_files = ();
my $xmfa_list_file;
my $gfa_file;
my @pack_files = ();
my $pack_list_file;
my $vg_path;
my $sampling_interval = 50;
my $edge_cov = 0;
my $min_cov = 0;
my $help;

my %pack_names = ();
my @file_bases = ();

parse_args();

print(STDOUT "path_name\tpath_pos\txmfa_col\tnode_id\tpack_name\tcov\n");

foreach my $xmfa_file (@xmfa_files) {
	parse_xmfa_file($xmfa_file);
}

exit(0);


sub parse_xmfa_file {
	my $xmfa_file = shift();

	my %seq_ids = ();
	my %seq_names = ();
	my %seq_id_names = ();
	my %block_seqs = ();
	my %block_seq_info = ();
	my %seq_pos_nodes = ();
	my %xmfa_nodes = ();
	my %pack_covs = ();
	my $gfa_parsed = 0;
	my $use_seq_names = 1;
	my $block_seq_count = 0;
	my $block_id = 1;
	my $valid_seq = 0;
	my $xmfa_col = -1;
	my $seq_id;
	my $pos;

	open(XMFA, '<', $xmfa_file) or error("can't open xmfa file: $!");

	print(STDERR "processing xmfa file: $xmfa_file\n");

	while (my $line = <XMFA>) {
		chomp($line);

		if ($line =~ /^#/) {
			if ($line =~ /^#Sequence(\d+)File/) {
				$seq_id = $1;

				$seq_ids{$seq_id}++;

				my ($header_field, $header_val) = split(/\t/, $line);
				my $seq_name = $header_val;

				$seq_name =~ s/^.*\///;
				$seq_name =~ s/\.f[ast]*a//;
				$seq_name =~ s/\s+/_/g;

				if ($use_seq_names == 1 && exists($seq_names{$seq_name})) {
					print(STDOUT "duplicate sequence name: $seq_name found, use of seq names in output disabled, seq ids (from xmfa header) will be used instead\n");

					$use_seq_names = 0;
				}

				$seq_names{$seq_name}++;
				$seq_id_names{$seq_id} = $seq_name;
			}
		}

		elsif ($line =~ /^\>/) {
			if ($gfa_parsed == 0) {
				parse_gfa_file(\%seq_pos_nodes, \%seq_id_names, \%xmfa_nodes, $use_seq_names);

				$gfa_parsed = 1;

				foreach my $pack_file (@pack_files) {
					parse_pack_file($pack_file, \%xmfa_nodes, \%pack_covs);
				}
			}

			$line =~ s/^>\s*//;

			my ($seq_id_pos, $strand, $seq_file) = split(/\s+/, $line);
			($seq_id, $pos) = split(':', $seq_id_pos);
			my ($start, $stop) = split('-', $pos);

			$valid_seq = 0;

			if (defined($seq_id) && defined($start) && defined($stop) && defined($strand)) {
				$block_seq_count++;
				$block_seq_info{$seq_id}{'start'} = $start;
				$block_seq_info{$seq_id}{'stop'} = $stop;
				$block_seq_info{$seq_id}{'strand'} = $strand;
				$valid_seq = 1;
			}
		}

		elsif ($line =~ /^\=/) {
			my $block_len;
			my %seq_offset = ();

			foreach my $seq_id (keys %block_seqs) {
				$block_len = length($block_seqs{$seq_id});

				last();
			}

			my $block_col = 0;

			foreach my $block_col (0..($block_len - 1)) {
				$xmfa_col++;

				if ($xmfa_col % $sampling_interval != 0) {
					next();
				}

				foreach my $seq_id (sort { $a <=> $b } keys %block_seqs) {
					my $seq_name = $seq_id;

					if ($use_seq_names == 1 && exists($seq_id_names{$seq_id})) {
						$seq_name = $seq_id_names{$seq_id};
					}

					my $sub_seq;

					if ($block_col == 0) {
						$sub_seq = substr($block_seqs{$seq_id}, $block_col, 1);
					}

					else {
						$sub_seq = substr($block_seqs{$seq_id}, $block_col - $sampling_interval + 1, $sampling_interval);
					}

					$sub_seq =~ s/\-//g;

					my $sub_seq_len = length($sub_seq);

					if ($sub_seq_len == 0) {
						next();
					}

					$seq_offset{$seq_id} += $sub_seq_len;

					my $col_base = substr($block_seqs{$seq_id}, $block_col, 1);

					if ($col_base !~ /[acgtACGT]/) {
						next();
					}

					my $seq_offset = $seq_offset{$seq_id};
					my $seq_start = $block_seq_info{$seq_id}{'start'};
					my $seq_stop = $block_seq_info{$seq_id}{'stop'};
					my $seq_strand = $block_seq_info{$seq_id}{'strand'};

					my $seq_pos;

					if ($seq_strand eq '+') {
						$seq_pos = $seq_start + $seq_offset - 1;
					}

					else {
						$seq_pos = $seq_stop - $seq_offset + 1;
					}


					if (! exists($seq_pos_nodes{$seq_name}{$seq_pos})) {
						print(STDERR "node not found for path name: $seq_name, pos: $seq_pos\n");

						next();
					}

					my ($node_id, $node_index);

					if ($edge_cov) {
						$node_id = $seq_pos_nodes{$seq_name}{$seq_pos};
						$node_index = 0;
					}

					else {
						($node_id, $node_index) = split(/\t/, $seq_pos_nodes{$seq_name}{$seq_pos});
					}

					foreach my $file_base (@file_bases) {
						my $cov = 0;

						if (defined($node_id) && exists($pack_covs{$node_id}{$file_base})) {
							$cov = ${$pack_covs{$node_id}{$file_base}}[$node_index];

							if (! defined($cov)) {
								print(STDERR "coverage not defined for node id: $node_id, pack: $file_base, index: $node_index\n");

								next();
							}
						}


						if ($cov < $min_cov) {
							next();
						}

						print(STDOUT "$seq_name\t$seq_pos\t$xmfa_col\t$node_id\t$file_base\t$cov\n");
					}
				}
			}

			$block_seq_count = 0;
			undef %block_seqs;
			undef %block_seq_info;
			$block_id++;
		}

		else {
			if ($valid_seq == 1) {
				$block_seqs{$seq_id} .= "$line";
			}
		}
	}

	close(XMFA);

	return(0);
}


sub parse_gfa_file {
	my $seq_pos_nodes_ref = shift();
	my $seq_id_names_ref = shift();
	my $xmfa_nodes_ref = shift();
	my $use_seq_names = shift();

	my %path_nodes = ();
	my %node_lens = ();
	my %seq_names = ();

	foreach my $seq_id (keys %$seq_id_names_ref) {
		if ($use_seq_names == 1) {
			$seq_names{$$seq_id_names_ref{$seq_id}}++;
		}

		else {
			$seq_names{$seq_id}++;
		}
	}

	open(GFA, '<', $gfa_file) or error("can't open gfa file: $!");

	print(STDERR "processing gfa file: $gfa_file\n");

	while (my $line = <GFA>) {
		if ($line =~ /^P/) {
			chomp($line);

			my ($rec_type, $path_name, $nodes, $overlaps) = split(/\t/, $line);

			if (! exists($seq_names{$path_name})) {
				next();
			}

			my @nodes = split(',', $nodes);

			foreach my $node (@nodes) {
				$node =~ s/[\-\+]$//;
				$node_lens{$node} = 0;
			}

			push(@{$path_nodes{$path_name}}, @nodes);
		}
	}

	close(GFA);


	open(GFA, '<', $gfa_file) or error("can't open gfa file: $!");

	while (my $line = <GFA>) {
		if ($line =~ /^S/) {
			chomp($line);

			my ($rec_type, $node, $seq) = split(/\t/, $line);

			if (exists($node_lens{$node}) && defined($seq)) {
				$node_lens{$node} = length($seq);
			}
		}
	}

	close(GFA);


	foreach my $path_name (keys %path_nodes) {
		my $path_pos = 0;

		foreach my $index (0..$#{$path_nodes{$path_name}}) {
			my $node = $path_nodes{$path_name}[$index];

			$node =~ s/[\-\+]$//;

			if (! exists($node_lens{$node})) {
				print(STDERR "node: $node length not found\n");

				next();
			}

			my $len = $node_lens{$node};

			$$xmfa_nodes_ref{$node}++;

			my $start_pos = $path_pos + 1;
			my $stop_pos = $path_pos + $len;

			if ($stop_pos < $start_pos) {
				next();
			}

			my $index = 0;

			foreach my $pos ($start_pos..$stop_pos) {
				if ($edge_cov) {
					$$seq_pos_nodes_ref{$path_name}{$pos} = "$node";
				}

				else {
					$$seq_pos_nodes_ref{$path_name}{$pos} = "$node\t$index";
					$index++;
				}
			}

			$path_pos = $stop_pos;
		}
	}

	print(STDERR "completed processing gfa file: $gfa_file\n");

	return(0);
}


sub parse_pack_file {
	my $pack_file = shift();
	my $xmfa_nodes_ref = shift();
	my $pack_covs_ref = shift();

	my $file_base = $pack_file;

	my $pack_fh;

	if ($pack_file =~ /\.gz$/) {
		open($pack_fh, '-|', "gzip -dc $pack_file") or error("can't open pack file: $!");
	}

	elsif ($pack_file =~ /\.bz2$/) {
		open($pack_fh, '-|', "bzip2 -dc $pack_file") or error("can't open pack file: $!");
	}

	else {
		open($pack_fh, '<', $pack_file) or error("can't open pack file: $!");
	}

	print(STDERR "processing pack file: $pack_file\n");

	$file_base =~ s/^.*\///;
	$file_base =~ s/[\.\_].*$//;

	if (exists($pack_names{$pack_file})) {
		$file_base = $pack_names{$pack_file};
	}

	push(@file_bases, $file_base);

	my $file_type;

	while (my $line = <$pack_fh>) {
		chomp($line);

		if ($line =~ /^seq\.pos\tnode\.id/) {
			$file_type = 'pos_cov';

			next();
		}

		elsif ($line =~ /^node\.id\tpos\.covs/) {
			$file_type = 'seg_cov';

			next();
		}

		if (! defined($file_type)) {
			error("invalid pack/seg cov file format: $pack_file\n");
		}

		if ($file_type eq 'pos_cov') {
        	my ($seq_pos, $node_id, $node_offset, $cov) = split(/\t/, $line);

			if (! exists($$xmfa_nodes_ref{$node_id})) {
				next();
			}

			if ($edge_cov) {
				if ($node_offset == 0) {
					push(@{$$pack_covs_ref{$node_id}{$file_base}}, $cov);
				}
			}

			else {
				push(@{$$pack_covs_ref{$node_id}{$file_base}}, $cov);
			}
		}

		elsif ($file_type eq 'seg_cov') {
			my ($node_id, $covs) = split(/\t/, $line);

			if (! exists($$xmfa_nodes_ref{$node_id})) {
				next();
			}

			if ($edge_cov) {
				$covs =~ s/\,.*$//;

				push(@{$$pack_covs_ref{$node_id}{$file_base}}, $covs);
			}

			else {
				push(@{$$pack_covs_ref{$node_id}{$file_base}}, split(',', $covs));
			}
		}
	}

	close($pack_fh);

	return(0);
}



sub error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n");
	}

	exit(0);
}


sub arg_error {
	my $msg = shift();

	if (defined($msg)) {
		print(STDERR "error: $msg\n\n");
	}

	print("use option -h or --help to display help menu\n");

	exit(0);
}


sub parse_args {
	if ($#ARGV == -1) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	GetOptions ('x|xmfa=s{,}' => \@xmfa_files,
				'xmfalist=s' => \$xmfa_list_file,
				'g|gfa=s' => \$gfa_file,
				'p|pack=s{,}' => \@pack_files,
				'packlist=s' => \$pack_list_file,
				'interval=i' => \$sampling_interval,
				'e|edge' => \$edge_cov,
				'c|cov=i' => \$min_cov,
				'h|help' => \$help) or error('cannot parse arguments');

	if (defined($help)) {
		pod2usage(-exitval => 0, -verbose => 1);
	}

	if (defined($xmfa_list_file)) {
		open(XMFALIST, '<', $xmfa_list_file) or error("can't open xmfa list file: $!");

		foreach my $line (<XMFALIST>) {
			chomp($line);

			push(@xmfa_files, $line);
		}

		close(XMFALIST);
	}

	if (! @xmfa_files) {
		arg_error('at least 1 xmfa file must be provided');
	}

	if (! defined($gfa_file)) {
		arg_error('gfa file required');
	}

	if (defined($pack_list_file)) {
		open(PACKLIST, '<', $pack_list_file) or error("can't open pack list file: $!");

		foreach my $line (<PACKLIST>) {
			chomp($line);

			my ($file, $name) = split(/\t/, $line);

			push(@pack_files, $file);

			if (defined($name)) {
				$pack_names{$file} = $name;
			}
		}

		close(PACKLIST);
	}

	if (! @pack_files) {
		arg_error('at least 1 pack table file must be provided');
	}

	if ($sampling_interval < 1) {
		arg_error('sample frequency must be a positive integer');
	}

	return(0);
}


__END__

=head1 NAME

xmfa_graph_cov.pl

=head1 SYNOPSIS

xmfa_graph_cov.pl [options] > graph.packs.cov

=head1 DESCRIPTION

xmfa_graph_cov.pl produces graph segment coverages based on
sampled xmfa positions

=head1 AUTHOR

Brian Abernathy

=head1 OPTIONS

=head2 general options

 -x --xmfa      xmfa file(s)
                  ex: -x chr01.xmfa chr02.xmfa

 --xmfalist     text file containing list of xmfa file paths
                  1 xmfa file per line

1 or more xmfa files must be specified using -x/--xmfa and/or
--xmfalist. Multiple xmfa files are processed sequentially, so
it is more efficient to run multiple commands for multiple
chromosome genomes. (1 command for each chromosome) Output for
multiple chromosomes can be concatenated together, if desired.

 -g --gfa       genome gfa file, vg-based (required)

 --interval     xmfa sampling interval (sample coverage every
                  X multiple alignment columns)
                  defautl: 50

 -e --edge      report node leading edge coverage (instead of node
                  position coverage)
                  default: disabled

 -c --cov       minimum node coverage
                  do not report nodes falling below threshold
                  default: 0

 -p --pack      vg pack table or segment coverage file(s)
                  ex: -p sample1.pack.table sample2.seg.cov.gz

 --packlist     text file containing list of pack file paths
                  pack names (found in output header) may optionally
                  be specified in the second field
                  1 pack file (and name) per line, tab-delimited

1 or more pack (table or segment coverage) files must be specified
using -p/--pack and/or --packlist. Pack table files are generated
using the `vg pack -d` command. Pack segment coverage files are
generated using pack_table_to_seg_cov.pl, which is part of the
gfa_var_genotyper project. 
(https://github.com/brianabernathy/gfa_var_genotyper) Pack files
may be uncompressed or compressed using either gzip or bzip2.
(.gz or .bz2 file extension)

Note that pack names (found in output header) are automatically
extracted from pack file names by removing all trailing characters
after the first '.' or '_'. If you prefer to specify pack names,
the --packlist option must be used, with names appearing in the
file's second field.

=head2 help

 -h --help      display help menu

=cut

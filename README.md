# xmfa_tools 

### xmfa_tools.pl
Description:

    sort XMFA file by block (LCB) to generate gapped fasta
    files with fasta/xmfa coordinates and/or vg-based gfa files

    This is accomplished by defining a multi-pass sort using
    sequence IDs found in the XMFA header. The primary sort
    sequence ID is used to establish a backbone of blocks.
    Additional sorting passes recursively extend the backbone
    blocks using the non-primary sequence ID block positions.

    Blocks are oriented so that the primary sort sequence is
    always on the forward strand. Blocks that do not contain
    a primary sort alignment will be oriented based on the
    sequence alignment of highest sort priority. When the
    highest sort priority sequence alignment is on the reverse
    strand, all block sequences are reverse complemented.


Usage:

    xmfa_tools.pl -x xmfa.file [options]

Options:

  required:

     -x --xmfa     input xmfa file (required)

  processing:

     -p --print    print xmfa seq ids, names and files, then exit

     -s --sort     enable block sorting

     -o --order    sort order
                     provide one or more seq ids (from header)
                     default: sort by 1, 2, 3... 
                     example: --order 2 3 1 (sort by 2, then 3, then 1)
                     example: --order 3 (sort by 3, then 1, then 2)

     --nogapfilter disable invalid gap filter
                     default: invalid gap filtering is enabled
                     invalid gap example:
                     ACTAGCTGATG--------CTGACGTAATCGTGATGATCGATGCTGA
                     ACTAGCTGATGCTGACGTA--------ATCGTGATGATCGATGCTGA
                     ACTAGCTGATGCTGACGTA--------ATCGTGATGATCGATGCTGA

     --noseqnames  disable use of seq names in fasta and gfa files
                     By default, seq names (viewed using -p option)
                     are used to name output fasta files and gfa path
                     records. The --noseqnames option will disable
                     the use of seq names and use seq ids instead.

     -i --include  include only specified seq ids in output
                     default: include all
                     example: --include 1 3 7

     -t --threads  number of vg threads to use for gfa processing
                     Current vg versions yield a very modest gain
                     in processing speed with multithreading enabled.
                     default: 1

  output:

     --xmfaout     output xmfa file (requires -s/--sort or
                     -i/--include)

     --xmfasort    output sorted xmfa file
                     deprecated: use --xmfaout with -s/--sort

     -g --gfa      output gfa file (vg-based, similar to v1 spec)

     --gfapostfix  include gfa seq postfix in gfa path records
                     default: no postfix
                     example: --gfapostfix ".chr01"

     -v --vg       path to vg executable (required for gfa processing)
                     default: autodetect in $PATH (if available)

     --catgfapaths concatenate gfa paths with the same gfa name
                     default: output separate path records for each
                     block

     -c --coords   output fasta coords file
                     specifying a fasta coords file will generate
                     a gapped fasta file for each (selected) sequence

     -n --null     output fasta coords null record value
                     default: "NA"

     --fastadir    output fasta directory
                     default: current working directory

     --fapostfix   output gapped fasta file postfix
                     default: ".sort.gapped.fa"

     -l --linker   output fasta block sequence linker
                     default: no linker
                     example: -l "NNNNNNNNNN"

  help:

     -h --help     display help menu

---

### xmfa_graph_cov.pl

Description:

   Step through xmfa file(s) and report graph node/segment coverage for
   each pack file at specified sampling intervals.

Usage:

    xmfa_graph_cov.pl [options] > graph.packs.cov

Options:

  general options:

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

     -c --cov       minimum node coverage
                      do not report nodes falling below threshold
                      default: 0

     -p --pack      vg pack table or segment coverage file(s)
                      ex: -p sample1.pack.table sample2.seg.cov.gz

     --packlist     text file containing list of pack file paths
                      pack names (found in output header) may optionally
                      be specified in the second field
                      1 pack file (and name) per line, tab-delimited

    1 or more pack (table or segment coverage) files must be specified using
    -p/--pack and/or --packlist. Pack table files are generated using the
    `vg pack -d` command. Pack segment coverage files are generated using
    pack_table_to_seg_cov.pl, which is part of the gfa_var_genotyper
    project. (https://github.com/brianabernathy/gfa_var_genotyper) Pack
    files may be uncompressed or compressed using either gzip or bzip2. (.gz
    or .bz2 file extension)

    Note that pack names (found in output header) are automatically
    extracted from pack file names by removing all trailing characters after
    the first '.' or '_'. If you prefer to specify pack names, the
    --packlist option must be used, with names appearing in the file's
    second field.

  help:

     -h --help      display help menu

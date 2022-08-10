<p xmlns:dct="http://purl.org/dc/terms/">
<a rel="license" href="http://creativecommons.org/publicdomain/mark/1.0/">
<img src="http://i.creativecommons.org/p/mark/1.0/88x31.png"
     style="border-style: none;" alt="Public Domain Mark" />
</a>
<br />
This work is free of known copyright restrictions.
</p>

# xmfa_tools

Test data sets, scripts and additional documentation specific to skim-seq processing are also available at the [PanPipes](https://github.com/USDA-ARS-GBRU/PanPipes#demos "PanPipes") repository.

Description:

    xmfa_tools.pl can be used to create GFA sequence graphs,
    gapped fasta files and/or sorted multi-sequence alignments
    from raw XMFA files. Additionally, invalid gaps (see
    --nogapfilter option below) are also removed.

    Sorting is accomplished by defining a multi-pass sort using
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

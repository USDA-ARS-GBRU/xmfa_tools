# xmfa_tools 

### xmfa_tools.pl
Description:

	sort XMFA file by block (LCB) and generate gapped fasta
	files with fasta/xmfa coordinates

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

     -s --sort     enable block sorting

     -o --order    sort order
                     provide one or more seq ids (from header)
                     default: sort by 1, 2, 3... 
                     example: --order 2 3 1 (sort by 2, then 3, then 1)
                     example: --order 3 (sort by 3, then 1, then 2)

     -p --print    print ref seq ids, file and base names, then exit

     --nobasenames disable use of base names in fasta and gfa files
                     By default, base names (viewed using -p option)
                     are used to name output fasta files and gfa path
                     records. The --nobasenames option will disable
                     the use of base names and use seq ids instead.

     -i --include  include only specified seq ids in output
                     default: include all
                     example: --include 1 3 7

     -t --threads  number of vg threads to use for gfa processing
                     Current vg versions yield a very modest gain
                     in processing speed with multithreading enabled.
                     default: 1

  output:

     --xmfasort    output sorted xmfa file

     -g --gfa      output gfa file (vg-based, similar to v1 spec)

     --gfapostfix  include gfa seq postfix in gfa path records
                     default: no postfix
                     example: --gfapostfix ".chr01"

     -v --vg       path to vg executable (required for gfa processing)
                     default: autodetect in $PATH (if available)

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

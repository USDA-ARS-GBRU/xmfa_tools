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

     -p --print    print reference seq ids, then exit

     -i --include  include only specified seq ids in output
                     default: include all
                     example: --include 1 3 7

  output:

     --xmfasort    output sorted xmfa file

     -g --gfa      output gfa (v1) file

     --gfabasename include sequence 'base name' in gfa path records
                     If xmfa header 'Entry' records are present,
                     this option will automatically be disabled.
                     (use -p option to see xmfa seq base names)
                     default: no base name

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

     --fapostfix   output fasta file postfix
                     default: ".sort.fa"

     -l --linker   output fasta block sequence linker
                     default: no linker
                     example: -l "NNNNNNNNNN"

  help:

     -h --help     display help menu

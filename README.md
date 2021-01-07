# xmfa_tools 

### xmfa_block_sort.pl
Description:

	sort XMFA file by block (LCB) and generate gapped fasta
	files with fasta/xmfa coordinates

    This is accomplished by defining a multi-pass sort using
    sequence IDs found in the XMFA header. The primary sort 
    sequence ID is used to establish a backbone of blocks. 
    Additional sorting passes recursively extend the backbone
    blocks using the non-primary sequence ID block positions.

Usage:

    xmfa_block_sort.pl -x xmfa.file [options]

Options:

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

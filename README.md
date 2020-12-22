# xmfa_tools 

### xmfa_block_sort.pl
Description:

    sort XMFA files by block (LCB)

    This is accomplished by defining a multi-pass sort using
    sequence IDs found in the XMFA header. The primary sort 
    sequence ID is used to establish a backbone of blocks. 
    Additional sorting passes recursively extend the backbone
    blocks using the non-primary sequence ID block positions.

Usage:

    xmfa_block_sort.pl [options]

Options:

     -x --xmfa     xmfa file (required)

     -p --print    print reference seq ids, then exit

     -s --sort     sort order
                     provide one or more ref seq ids (from header)
                     default: sort by 1, 2, 3... 
                     example: -s 2 3 1 (sort by 2, then 3, then 1)
                     example: -s 3 (sort by 3, then 1, then 2)

     -h --help     display help menu

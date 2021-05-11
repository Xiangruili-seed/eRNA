import pysam
import sys

### Input varibles to set
# file to split on
unsplit_file = open(sys.argv[1])
# where to place output files
out_dir = open(sys.argv[2])
# variable to hold barcode index
CB_hold = open(sys.argv[3])
itr = 0
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile( unsplit_file, "rb")
for read in samfile.fetch( until_eof=True):
    # barcode itr for current read
    CB_itr = read.get_tag( 'CB')
    # if change in barcode or first line; open new file
    if( CB_itr!=CB_hold or itr==0):
        # close previous split file, only if not first read in file
        if( itr!=0):
            split_file.close()
        CB_hold = CB_itr
        itr+=1
        split_file = pysam.AlignmentFile( out_dir + "CB_{}.bam".format( itr), "wb", template=samfile)

    # write read with same barcode to file
    split_file.write( read)
split_file.close()
samfile.close()

from pyliftover import LiftOver



def convert_coordinated(chr,pos):
    lo = LiftOver('data/hg19ToHg38.over.chain.gz')
    chr_string = 'chr' + str(chr)
    pos_int = int(pos)
    ret_val = lo.convert_coordinate(chr_string, pos_int)
    return ret_val[0][1]



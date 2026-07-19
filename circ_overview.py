import pycircos
import pybigwig
def circos_wg(chrom_sizes):
    Garc=pycircos.Garc
    Gcircle=pycircos.Gcircle
    circle=Gcircle(#browser page size..)
    with open(chrom_sizes,'r') as f:
        for line in f.readline():
            line=line.rstrip().split('\t')
            chrom=line[0]
            length=int(line[-1])
            crom_arc=Garc(arc_id=chrom,size=length,interspace=2,raxis_range=(),label_position=..,label_visible=True)
            circle.add_argc(chrom_arc)
    return circle 

def circos_bw(circle,bwtype,args=()):
    #args=() is bw list like (rhppr_bw_fwd.bw,rhppr_bw_rev.bw, opened with pybigwig.open)
    #add place holder for the circle.lineplot()
    circlines=[]
    bwfwd=pyBigWig.load(args.get(#))
    bwrev=pyBigWig.load(args.get(#))
    #combine them together 
    circle.lineplot()#the bw coverage/read depth as the lines 
    return circle

def 


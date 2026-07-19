class Loci(object):
    def __init__(self,chrom,start,end,strand,class):
        self.chrom=chrom
        self.start=int(start)
        self.end=int(end)
        self.strand=strand
        self.class=class
    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end}:{self.strand}"
    def __iter__(self):
        for i in range(self.start,self.end+1):
            yield i
class PCGene(Loci):#class for protein coding gene (mRNA region)
    def __init__(self,chrom,start,end,strand,cds,gene,transcript_id,proxy):
        super().__init__(chrom=chrom,start=start,end=end,strand=strand)
        self.proxy=proxy
        self.timeout=timeout 
        self.cds=cds
        self.gene=gene 
        self.transcript_id=transcript_id
    def __tx_range__(self):
        tx_info=defaultdict(list)
        domains=namedtuple('domain genome region',['category','type','description','start','end'])
        for i in self.cds.keys():
            curr_prot=

from Bio import SeqIO
def parsefa()
def gen_msa():
    pssm_gen=[]
    pool=Pool(1)
    res=pool.map(search,seqdir)
    pool.close()
    pool.join()



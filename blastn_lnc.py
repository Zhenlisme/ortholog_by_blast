## 2020.05.19
## This scripts was writtern for repeat the nature's work on lncRNA.(https://www.nature.com/articles/s41586-019-1341-x)

import os,subprocess,gzip
from multiprocessing import Pool

class NatureRedone():
    def __init__(self,gtfdict):
        self.gffs=["/".join([gtfdict,f]) for f in os.listdir(gtfdict)]
    def transfertobed12(self,bed12dict):
        for file in self.gffs:
            bed12file=".".join([".".join(file.split("/")[-1].split(".")[0:-1]),"bed12"])
            subprocess.run("gtfToGenePred %s stdout|genePredToBed stdin %s" %
                           (file,"/".join([bed12dict,bed12file])),shell=True)
    def extract_exon(self,genomedict,bed12dict,lncdict):
        for file in os.listdir(genomedict):
            if not file.endswith(".fasta"):
                continue
            bnm=file.split(".")[0]
            genomefile="/".join([genomedict,file])
            bedfile="/".join([bed12dict,".".join([bnm,"lncRNA","bed12"])])
            lncfile="/".join([lncdict,".".join([bnm,"lnc","fasta"])])
            subprocess.run("bedtools getfasta -fi %s -bed %s -split -s -name > %s" %
                          (genomefile,bedfile,lncfile),shell=True)
    def mkblastdb(self,lncdict):
        subprocess.run("cat %s/*.lnc.fasta > %s" % (lncdict,"lncseq.fasta"),shell=True)
        subprocess.run("makeblastdb -in %s -dbtype nucl -out %s" % ("lncseq.fasta","lncdb"),shell=True)
    def seqalign(self,query,db,out):
        subprocess.run("blastn -query %s -out %s -db %s -outfmt 6 -evalue 1e-3 -num_threads 20" %
                       (query,db,out), shell=True)

    def sonfilter(self,inplist):
        opline=[]
        print("subprocess was done!")
        for line in inplist:
            print(line)
            line=line.decode()
            splitline = line.split("\t")
     #       first, second = splitline[0].split(":")[0], splitline[1].split(":")[0]
     #       if not second in self.selected_lnc_list or not first in self.selected_lnc_list:
     #           continue
            if float(splitline[2]) < 0.1 or abs(int(splitline[6]) - int(splitline[7])) < 50 or abs(
                    int(splitline[8]) - int(splitline[9])) < 50:
                continue
            opline.append(line)
        return "".join(opline)
    def filter(self,selected_lnc,align_result,opfile):
        oplines=[]
    #    self.selected_lnc_list=[]
    #    for file in os.listdir(selected_lnc):
    #        with open("/".join([selected_lnc,file]),'r') as F:
    #            F.readline()
    #            self.selected_lnc_list.extend([line.split(",")[2] for line in F])
        print("Begine reading!")
        with gzip.open(align_result, 'r') as F:
            alignlist=F.readlines()
        print("Have read done!")
        persize=int(len(alignlist)/40)
        steps=[i for i in range(0,len(alignlist),persize)]
        steps.append(len(alignlist))
        rangelist=[[i, steps[steps.index(i) + 1]] for i in steps if steps.index(i) != len(steps) - 1]
        mypool=Pool(processes=20)
        for p in rangelist:
            goalist=alignlist[p[0]:p[1]]
            oplines.append(mypool.apply_async(self.sonfilter,(goalist,)).get())
        mypool.close()
        mypool.join()
        with open(opfile,'w') as opf:
            opf.writelines(oplines)

if __name__ == "__main__":
    os.chdir("/home/zhluo/Project/lizhen_lncRNA/colinearity/NatureWorkReDone/")
    NatureWork=NatureRedone(gtfdict="lnc_ano")
    #NatureWork.transfertobed12(bed12dict="bed12dict")
    #NatureWork.extract_exon(genomedict="genomes",bed12dict="bed12dict",lncdict="lncdict")
    #NatureWork.mkblastdb(lncdict="lncdict")
    NatureWork.filter(selected_lnc="lnclist",align_result="lncseq.aln.xlsx.gz",opfile="full_filtered_aln_result.xlsx")


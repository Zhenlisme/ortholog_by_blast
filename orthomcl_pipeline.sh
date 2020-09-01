

cd /home/zhluo/Project/lizhen_lncRNA/colinearity/NatureWorkReDone/lnc_ano/

for i in `ls *.gtf`;do bnm=`basename $i .lncRNA.gtf`;awk -v OFS="\t" '{print $1","$2","$3","$10","$7,$4,$5}' $i|sort -k1,1 -k2,2n > ${bnm}.bmg.bed;done

for i in `ls *.bmg.bed`;do bnm=`basename $i .bmg.bed`;bedtools merge -i $i > ${bnm}.amg.bed;done

for i in `ls *.amg.bed`;do bnm=`basename $i .amg.bed`;cat $i|tr -s "," "\t" |sort -k1,1 -k4,4 -k6,6n -k7,7rn -k3,3r|awk -v OFS="\t" '{print $1,$2,$3,$6,$7,"1000",$5,".", "gene_id "$4" transcript_id "$4}'> ${bnm}.amg.gtf;done

for i in `ls *.amg.gtf`;do bnm=`basename $i .amg.gtf`;gtfToGenePred $i stdout|genePredToBed stdin ${bnm}.amg.bed12;done

for i in `ls *.bed12`;do bnm=`basename $i .amg.bed12`;bedtools getfasta -fi ../genomes/${bnm}.fasta -bed $i -split -s -name > ../lncdict/${bnm}.lnc.fasta;done

cd /home/zhluo/Project/lizhen_lncRNA/colinearity/NatureWorkReDone/lncdict

for i in `ls *.lnc.fasta`;do bnm=`basename $i .lnc.fasta`;orthomclAdjustFasta ${bnm} $i 1;done
for i in `ls *.fasta|grep -v ".lnc.fasta"`;do mv $i ../complaintFasta/;done

cd ../

##过滤掉长度小于50的序列，10000指终止密码子个数小于10000
/home/nazhang/Fish_genome/source/orthomclSoftware-v2.0.9/bin/orthomclFilterFasta complaintFasta/ 50 10000

mv goodProteins.fasta goodlnc.fasta

makeblastdb -in goodlnc.fasta -dbtype nucl -out goodlncdb
blastn -query goodlnc.fasta -out all-all.blastn -db goodlncdb -outfmt 6 -evalue 1e-3 -num_threads 30 &
/home/nazhang/Fish_genome/source/orthomclSoftware-v2.0.9/bin/orthomclBlastParser all-all.blastn complaintFasta/ >similarSequences.txt
/home/nazhang/Fish_genome/source/orthomclSoftware-v2.0.9/bin/orthomclInstallSchema my_orthomcl_dir/orthomcl.config my_orthomcl_dir/install_schema.log
/home/nazhang/Fish_genome/source/orthomclSoftware-v2.0.9/bin/orthomclLoadBlast my_orthomcl_dir/orthomcl.config similarSequences.txt
/home/nazhang/Fish_genome/source/orthomclSoftware-v2.0.9/bin/orthomclPairs my_orthomcl_dir/orthomcl.config orthomcl_pairs.log cleanup=no
/home/nazhang/Fish_genome/source/orthomclSoftware-v2.0.9/bin/orthomclDumpPairsFiles my_orthomcl_dir/orthomcl.config 2>dump.err &
mcl mclInput --abc -I 1.5 -o mclOutput
/home/nazhang/Fish_genome/source/orthomclSoftware-v2.0.9/bin/orthomclMclToGroups cluster- 1000 < mclOutput > groups.txt
/home/nazhang/Fish_genome/source/orthomclSoftware-v2.0.9/bin/orthomclSingletons goodlnc.fasta groups.txt >> singletons.txt

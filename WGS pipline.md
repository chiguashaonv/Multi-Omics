# 全基因组测序数据分析

[参考资料1](https://mp.weixin.qq.com/mp/homepage?__biz=MzAxOTUxOTM0Nw==&hid=1&sn=d945cf61bd86e85724e146df42af5bcc&scene=1&devicetype=android-26&version=26060637&lang=zh_CN&nettype=WIFI&ascene=7&session_us=gh_2942f3f5dbfe&wx_header=1)

[参考资料2](http://www.biotrainee.com/home.php?mod=space&uid=378&do=thread&view=me&type=thread&order=dateline&from=space&page=5)

## 1. 下载数据并且了解数据

### 数据下载

[数据下载](ftp://ftp.kobic.re.kr/pub/KPGP/2015_release_candidate/WGS/KPGP-00001/)
下载下来的数据是一个人的全基因组测序数据，包括６个lane：

### 数据检验

![数据](http://imglf5.nosdn0.126.net/img/SWliemNmRGVaVmw5MGw4eEt4MDZTUG9SU0dvYVlESTQ3TjVHZzgvMnVsVUxWN2RSYm52ZU1nPT0.png?imageView&thumbnail=1680x0&quality=96&stripmeta=0)

md5文件是用来检验在文件传输过程中有没有发生数据的遗失或者损坏等：

![md5检验](http://imglf6.nosdn0.126.net/img/SWliemNmRGVaVmw5MGw4eEt4MDZTQjRWbnNFZFJ6OEl2Vkx0TC8zM2VRL0dYalVHNjFaWGhnPT0.png?imageView&thumbnail=1680x0&quality=96&stripmeta=0)

如上图所示，下载下来的md5文件的md5码和下载下来的fastq.gz文件的md5码一致，说明了文件在下载过程中没有缺失或者损坏

### 了解数据

在我们下载的数据中包括了６个lane，在处理过程中，我们一般是在对比后对bam文件进行合并

## 2. 质控

[数据质控](https://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798281&idx=1&sn=c3448e0e656a38808d0000ac8337e25d&scene=19#wechat_redirect)

```powershell
ls *.gz | while read id; do
	fastqc -t 10 -o ~/project/technique/fastqc --noextract $id
done
```
**multiqc**
```powershell
multiqc ./fastqc/
```
速度非常快,这个数据的[结果](https://github.com/chiguashaonv/Multi-Omics/blob/master/WGS/multiqc_report.html)也非常好．

## 3. 比对

**bwa**

[bwa参考１](http://starsyi.github.io/2016/05/24/BWA-%E5%91%BD%E4%BB%A4%E8%AF%A6%E8%A7%A3/) 
[bwa参考２](https://blog.csdn.net/oxygenjing/article/details/77747750)

### 3.1　建立索引
```powershell
bwa index [-p prefix] [-a algoType] <in.db.fasta>
```
在我们的服务器上已经建立了bwa的索引文件，因此这一步就可以省略

### 3.2 基因组比对
```powershell
#用法：mem Usage: bwa mem [options] ref.fa reads.fq [mates.fq]
for i in $(seq 1 6);do
	bwa mem -t 10 -R "@RG\tID:L${i}\tPL:ILLUMINA\tSM:KPGP-00001" -M \
	~/reference/Genome/BWAIndex/hg19.fa \
	~/project/technique/sequence/KPGP-00001_L${i}_R1.fq.gz ~/project/technique/sequence/KPGP-00001_L${i}_R2.fq.gz > ~/project/technique/bwa.mapping/KPGP-00001_L${i}.sam 2> ./mem-pe.log
done

```
六个lane的数据，每个数据的每个reads大约是5-6G,比对之后生成的６个sam文件大约每个50G左右．生成的sam文件排序为按照原本的reads顺序进行排序．
在阅读资料时，我发现在使用gatk call变异的时候，比对过程中设置的`-R`这个选项是必须的，所以必须在这步中进行设置，这个命令的解释可以参考这篇[文章](https://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798296&idx=1&sn=790d0141eec792b25083c63e87fee14c&scene=19#wechat_redirect)，但是我发现这个问题的时候已经跑完比对过程了，如果重新再跑一次，时间太久，因此需要寻找为SAM/BAM添加Read Groups的[方法](https://www.jianshu.com/p/215b17c12174):
```powershell
TAG='@RG\tID:xzg\tSM:Ebola\tLB:patient_100'　＃一般来讲，ID指lane id，SM指样本id,LB指测序文库的名字，重要性不高
# Add the tags during alignment 在比对是进行添加
bwa mem -R $TAG $REF $R1 $R2 | samtools sort > bwa.bam
samtools index bwa.bam
# Add tags with samtools addreplacerg　在比对后进行添加
samtools addreplacerg -r $TAG bwa_sorted.bam -o bwa_sorted_with_rg.bam
 ```
所以在我们下面转换为bam文件之后，添加rg信息进去



### 3.3 sam文件转换

为了节省储存空间，所以我们将比较大的sam文件转换为二进制的bam文件
```powershell
#bam to sam
samtools view -h abc.bam > abc.sam
#sam to bam
samtools view -@ 10 -b -S abc.sam > abc.bam
```
```powershell
for i in $(seq 1 6);do
	samtools view -@ 10 -b -S KPGP-00001_L${i}.sam > KPGP-00001_L${i}.bam
done
```
在转换为bam文件之后，每个文件的大小变为了15,16G，接下来，我们添加一下上面命令中漏掉的RG信息
```powershell
for i in $(seq 1 6):do
	TAG='@RG\tID:L{i}\tSM:KPGP-00001'
	samtools addreplacerg -r $TAG KPGP-00001_L{i}.bam -o KPGP-00001_L{i}.with_rg.bam
done
```

### 3.4 将bam文件按照位置进行排序
```powershell
for i in $(seq 1 6);do
	samtools sort -@ 10 -O bam -o KPGP-00001_L${i}.sort.bam KPGP-00001_L${i}.bam
done
```

## 4. 去除pcr重复

## 5. 合并不同的lane

## 5. 局部重比对

## 6. 变异检测



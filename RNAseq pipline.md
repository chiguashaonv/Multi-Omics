# RNAseq pepline

## 数据下载与处理
### 将数据从sra格式转化为fastq.gz格式
对于双端测序，用`fastq-dump`命令的`--split-3`选项
```powershell
fastq-dump --split-3 -A xxx.sra \
--gzip \
-O outpath 
```
  
### 检测数据质量QC
```powershell
fastqc -o outpath --noextarct xxx.fastq.gz
```

### 自动化处理数据
可同时处理双端数据
```powershell
fastp -i xxx.1.fastq.gz -o xxx.1.out.fastq.gz -I xxx.2.fastq.gz -O xxx.2.out.fastq.gz
```

## 序列比对

### 收集并了解不同的比对软件
[比较不同的比对软件的优势以及劣势](https://blog.csdn.net/qazplm12_3/article/details/76700045)。常用的比对软件包括：STAR, TopHat和HISAT2
>STAR具有最高比例的在基因组上有唯一比对位置的reads，尤其是对读长为300 nt的MCF7样品也有最高的比对率。与TopHat和HISAT2不同，STAR只保留双端reads都比对到基因组的序列，但对低质量的比对 (允许更多的错配碱基和soft-clip事件) 容忍度高。这一点在长reads >(MCF7-300)样品中的体现更为明显。TopHat则不允许soft-clip事件。

>在比对速度方面，HISAT2比STAR快2.5倍，比TopHat快大约100倍。

>soft-clip事件: 即reads末端存在低质量碱基或接头导致比对不上的, STAR会自动尝试截去未比对部分，只保留比对上的部分。

接下来，我们对这三种软件都安装并且进行尝试：[STAR安装](https://github.com/alexdobin/STAR)，[TopHat2安装](http://ccb.jhu.edu/software/tophat/index.shtml)，[HISAT2安装](https://ccb.jhu.edu/software/hisat2/index.shtml)

在了解这三种软件的创始人以及使用目的时，我发现TopHat2和HISAT2都是约翰霍普金斯大学计算生物学中心发表的软件，而且再TopHat2的首页上也已经提出了：
>Please note that TopHat has entered a low maintenance, low support stage as it is now largely superseded by HISAT2 which provides the same core functionality (i.e. spliced alignment of RNA-Seq reads), in a more accurate and much more efficient way.

因此，我们下面只尝试STAR和HISAT2这两种方法

**STAR**

[使用手册](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

[生信菜鸟团的使用说明](http://www.bio-info-trainee.com/727.html)

[参考阅读１](https://www.jianshu.com/p/eca16bf2824e)
[参考阅读２](http://www.bioinfo-scrounger.com/archives/288)
[参考阅读３](http://starsyi.github.io/2016/05/24/SAM%E6%96%87%E4%BB%B6%E6%A0%BC%E5%BC%8F%E8%AF%B4%E6%98%8E/)

1. 生成genome Index

```powershell
STAR \
 --runThreadN 30	\
 --runMode genomeGenerate	\
 --genomeDir /home/liuke/reference/genomeDir \
 --genomeFastaFiles /home/liuke/reference/hg19.fa	\
 --sjdbGTFfile /home/liuke/reference/hg19.gtf	\
 --sjdbOverhang 99
```

2. 进行序列比对

```powershell
STAR \
--runThreadN 30 \
--genomeDir genomeDir \
--readFilesCommand zcat \
--readFilesIn fastq1 fastq2 \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix ${mapping_dir}/${sample_id}. 
```
输出的文件：![图片](http://imglf5.nosdn0.126.net/img/SWliemNmRGVaVmtnZE1GNnkzREpEV2tmWmkxS0NuZkFvcDVUT2pxMEROamEwNWMwdVNSL0hBPT0.png?imageView&thumbnail=1680x0&quality=96&stripmeta=0)
Aligned.sortedByCoord.out.bam: 比对结果
SJ.out.tab：包含剪接信息

>column 1: chromosome

>column 2: first base of the intron (1-based)

>column 3: last base of the intron (1-based)

>column 4: strand (0: undefined, 1: +, 2: -)

>column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5:AT/AC, 6: GT/AT

>column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)

>column 7: number of uniquely mapping reads crossing the junction

>column 8: number of multi-mapping reads crossing the junction

>column 9: maximum spliced alignment overhang


3. STAR软件可以进行二次比对

>为了发现更加灵敏的new junction，STAR建议使用2-pass mode，其能增加检测到的new junction数目，使得更多的splices reads能mapping到new junction。因此STAR先用一般参数做一遍mapping，收集检测到的junction信息，然后利用这已经annotated junction来做第二次mapping．

3.1 重新建立索引

```powershell
STAR \
--runThreadN 10 \
--runMode genomeDir \
--genomeDir new.path.to.genomeDir \
--genomeFastaFiles /home/liuke/reference/hg19.fa \
--sjdbGTFfile /home/liuke/reference/hg19.gtf \
--sjdbOverhang 99 \
--sjdbFileChrStartEnd all.SJ.out.tab.list \ #相比与第一次建立索引，只增加了一个命令选项，就是把SJ.out.tab文件加入到建立索引中
```

























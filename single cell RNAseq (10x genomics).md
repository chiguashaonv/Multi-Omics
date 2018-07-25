# single cell RNAseq analysis (10x genomics)

课题组开了关于单细胞方向的课题，在正式的拿到自己的测序数据之前，先用了别人的数据熟悉流程，观察别人的思路，因此，下载了2018年发表在nature medicine上的一篇文章的数据，这篇文章的数据存储在EMBL-EBI上，分为两个数据集，E-MTAB-6149和E-MTAB-6653，建库方法是10x genomics，这篇文章中也使用了10x genomics推荐的pipline－－Cell Ranger，下面，就按照文章的思路，对数据进行分析．

## 下载数据

```powershell
for i in `cat download.link.file`;do mwget $i;done
```

## 表达矩阵

[参考文件1](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)

下载下来的数据集首先用md5验证在数据传输过程中是否发生了损伤
```powershell
vi 6149.md5
for i in `ls *.gz`;do md5sum $i >> 6149.md5;done
```
将计算出来的md5与下载的表格中的md5比较，检验是否一致，在我的结果中，完全一致，说明数据在下载时已经下载完全

安装cellranger：
```powershell
wget -O cellranger-2.1.1.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-2.1.1.tar.gz?Expires=1532453768&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTIuMS4xLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTUzMjQ1Mzc2OH19fV19&Signature=CHM31vix7Nm2U6yPntuZX4rjIjcEvVL2y3ORUOF2Elo9AOOgXIFLWRocBGshtnonZcCYHu4A4VoQgRSYGFKEv6XSZEPbvToEt7tLxSq0KWIcWCfWdrJtcAJfOhj1YFVL~Fvjc-WDWdxmLenNYREQSbmFLn175Dn1sCGgKhuAmy3oJBbYwrULNraxinMhGEJKQm-vnr0C-FKRDLCnaMRoRpObIlncyJYYoGKv4SG2yaMh7d9Wqkc3Rt1IqNEDmC4jgLD8HuOYeKY6-fxQF4E-zh10LO9eU37YsHr3X8ng~W3IecfPKaU7hmMNCdZVKyxORIDbC~9p6-Is1bfNWWexxA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

tar -xzvf cellranger-2.1.1.tar.gz
```

下载reference:
```powershell
#下载hg19版本
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-1.2.0.tar.gz

tar -xzvf refdata-cellranger-hg19-1.2.0.tar.gz
```

测序数据采用了两种不同的建库方法，Single Cell 3' v1和Single Cell 3' v2，这两种方法测序的read1,read2的内容不同，具体可以参考[Sequencing Requirements for Single Cell 3'](https://support.10xgenomics.com/single-cell-gene-expression/index/doc/specifications-sequencing-requirements-for-single-cell-3)


在阅读cellranger的使用说明时，我发现对于这个软件来说，输入的fastq文件名字要有特定的格式：`Sample1_S1_L001_I1_001.fastq.gz`，因此，我在使用这个软件的时候，就需要对我们的数据进行改名，才能继续使用，因为我们的数据是一个样本一个fastq(R1,R2,R3)文件，所以在改名之后，每个样本都需要运行一次`cellranger count`拿到每一个样本的表达矩阵，然后用`cellranger aggr`命令进行合并．下面，以`Sample1_S1_L001_R1_001.fastq.gz,Sample1_S1_L001_R2_001.fastq.gz,Sample1_S1_L001_I1_001.fastq.gz`这个名字为例，写出示例代码．

```powershell
cellranger count --id=samplename --sample=Sample1 \
    --transcriptome=/home/liuke/reference/refdata-cellranger-hg19-1.2.0 \
    --fastqs=folder_path.of.fastqs \
    --expect-cells=4000 --localcores=20
```


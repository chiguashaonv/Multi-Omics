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

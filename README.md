## Uniprot API

#### Uniprot API of UniProt website REST API

利用Uniprot官网的API爬取Uniprot数据库

简单例子：[初识Uniprot API](http://www.bioinfo-scrounger.com/archives/417)

具体需求：抓取Uniprot蛋白数据库中的指定蛋白对应数据

脚本：`UniProt_website_REST API_GO.R`，通过REST API获取对应蛋白xml格式的数据并解析，用多线性加快爬取速度



#### Proteins REST API

利用EBI的protein API爬取Uniprot数据库，以及其他多维的蛋白数据（protein and genome information）

教程及例子：[https://www.ebi.ac.uk/proteins/api/doc/](https://www.ebi.ac.uk/proteins/api/doc/)

相比上面UniProt website REST API方法，此方法的抓取速度有非常大的提升；因为该方法可以一次性最多抓取100个蛋白

脚本：`Proteins_REST_API_GO.R`，通过REST API获取对应的多个蛋白的数据，主要代码是将爬取的数据转化为json格式

    d <- GET(url, accept("application/json"))
    mydata <- toJSON(content(d))
    mydata <- fromJSON(mydata)




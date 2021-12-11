library(tidyverse)

gtf <- readr::read_tsv("herv.gtf", col_names=c("chr","src","type","s","e","score","strand","frame","anno"))
parsed_gtf <- gtf %>% transmute(len=e-s, gid=sub("gene_id \"(.*)\"; transcript_id.*", "\\1", anno))
len_distrib <- ggplot(parsed_gtf,aes(x=gid,y=len))+geom_boxplot()+
  geom_jitter(color="blue")+
  labs(x="Gene ID", y="Length", title = "Length distribution of HERV Answer")
ggsave("answer_len_distrib.png",len_distrib)

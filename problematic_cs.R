library(reshape2)
library(ggplot2)
library(dplyr)

#read all seqs
seqs <- c(read_uniprot(paste0(pathway, "sp2010_2015.txt"), euk = TRUE),
          read_uniprot(paste0(pathway, "sp1950_2010.txt"), euk = TRUE))

#train model on all seqs
nokmer_model <- train_hsmm(seqs, aaaggregation)

preds <- predict(nokmer_model, seqs)

df <- pred2df(preds)

#get prediction and divide data - did model recognise properly seqs on whose it was trained?
hard_seqs <- t(sapply(seqs[which(df[["sp.probability"]] < 0.5)], function(single_seq) {
  cut <- attr(single_seq, "sig")[2]
  single_seq[(cut-2):(cut+1)]
}))

easy_seqs <- t(sapply(seqs[which(df[["sp.probability"]] >= 0.5)], function(single_seq) {
  cut <- attr(single_seq, "sig")[2]
  single_seq[(cut-2):(cut+1)]
}))



rownames(tmp) <- NULL

dat <- data.frame(c(rep("hard", 20), rep("easy", 20)), 
                  a()[-1], 
                  rbind(apply(hard_seqs, 2, function(i) 
                    table(factor(i, levels = a()[-1]))),
                    apply(easy_seqs, 2, function(i) 
                      table(factor(i, levels = a()[-1])))))

colnames(dat) <- c("rec", "aa", paste0("P", c(-3:-1, 1)))
mdat <- melt(dat, id.vars = c(1, 2))
colnames(mdat) <- c("rec", "aa", "pos", "count")

mdat <- cbind(mdat, freq = c((mdat %>% filter(rec == "hard") %>% select(count) %>%
     unlist)/nrow(hard_seqs),
  (mdat %>% filter(rec == "easy") %>% select(count) %>%
     unlist)/nrow(easy_seqs)))

ggplot(mdat, aes(x = pos, y = freq, fill = rec)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_wrap(~aa, ncol = 5)

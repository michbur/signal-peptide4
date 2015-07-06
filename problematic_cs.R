library(reshape2)
library(ggplot2)
library(dplyr)
library(biogram)
library(seqinr)
librfary(signalHsmm)

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


dat <- data.frame(c(rep("no", 20), rep("yes", 20)), 
                  a()[-1], 
                  rbind(apply(hard_seqs, 2, function(i) 
                    table(factor(i, levels = a()[-1]))),
                    apply(easy_seqs, 2, function(i) 
                      table(factor(i, levels = a()[-1])))))

colnames(dat) <- c("rec", "aa", paste0("P", c(-3:-1, 1)))

mdat <- melt(dat, id.vars = c(1, 2))
colnames(mdat) <- c("rec", "aa", "pos", "count")

mdat <- cbind(mdat, freq = c((mdat %>% filter(rec == "no") %>% select(count) %>%
                                unlist)/nrow(hard_seqs),
                             (mdat %>% filter(rec == "yes") %>% select(count) %>%
                                unlist)/nrow(easy_seqs)))

cs_unigram <- mdat

#ngrams

give_me_df <- function(u, d, hard_seqs, easy_seqs) {
  dat <- data.frame(c(rep("no", length(u)^2), rep("yes", length(u)^2)), 
                    paste0(create_ngrams(2, u), "_", d), 
                    rbind(apply(seq2ngrams(hard_seqs, 2, d = d), 2, function(i) 
                      table(factor(i, levels = paste0(create_ngrams(2, u), "_", d)))),
                      apply(seq2ngrams(easy_seqs, 2, d = d), 2, function(i) 
                        table(factor(i, levels = paste0(create_ngrams(2, u), "_", d))))))
  colnames(dat) <- c("rec", "aa", paste0("P", c(-3:(-1 - d))))
  
  mdat <- melt(dat, id.vars = c(1, 2))
  colnames(mdat) <- c("rec", "aa", "pos", "count")
  
  mdat <- cbind(mdat, freq = c((mdat %>% filter(rec == "no") %>% select(count) %>%
                                  unlist)/nrow(hard_seqs),
                               (mdat %>% filter(rec == "yes") %>% select(count) %>%
                                  unlist)/nrow(easy_seqs)))
  
  levels(mdat[["aa"]]) <- decode_ngrams(levels(mdat[["aa"]]))
  
  #strange bug with NA
  
  # test1 <- test_features(c(rep(0, nrow(hard_seqs)), rep(1, nrow(easy_seqs))),
  #                        count_ngrams(rbind(hard_seqs, easy_seqs), n = 2, 
  #                                     u = a()[-1], d = 0, pos = TRUE),
  #                        adjust = NULL)
  # imp_ngrams  <- cut(test1, breaks = c(0, 0.05, 1))[[1]]
  # 
  # dec_ngrams <- decode_ngrams(imp_ngrams)
  # positions <- paste0("P-", ngrams2df(imp_ngrams)[, "position"])
  # 
  # lapply(1L:length(positions), function(i)
  #   mdat[mdat[["pos"]] == positions[i] & mdat[["aa"]] == dec_ngrams[i], ])
  
  ordered_ngrams <- order(mdat %>% filter(rec == "no") %>% select(freq) %>% unlist -
                            mdat %>% filter(rec == "yes") %>% select(freq) %>% unlist,
                          decreasing = TRUE)
  
  rbind((mdat %>% filter(rec == "no"))[ordered_ngrams, ][1L:20, ],
        (mdat %>% filter(rec == "yes"))[ordered_ngrams, ][1L:20, ])
}

cs_bigram0 <- give_me_df(a()[-1], 0, hard_seqs = hard_seqs, easy_seqs = easy_seqs)
cs_bigram1 <- give_me_df(a()[-1], 1, hard_seqs = hard_seqs, easy_seqs = easy_seqs)
cs_bigram0d <- give_me_df(1L:4, 0, hard_seqs = degenerate(hard_seqs, aaaggregation), 
                         easy_seqs = degenerate(easy_seqs, aaaggregation))
cs_bigram1d <- give_me_df(1L:4, 1, hard_seqs = degenerate(hard_seqs, aaaggregation), 
                         easy_seqs = degenerate(easy_seqs, aaaggregation))

save(cs_unigram, cs_bigram0, cs_bigram1, 
     cs_bigram0d, cs_bigram1d, file = "sumpres.RData")


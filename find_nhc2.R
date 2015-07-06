library(hmeasure)

find_nhc2 <- function (protein, signal = NULL) {
  protein <- toupper(protein)
  if (is.null(signal)) 
    signal <- attr(protein, "sig")
  sig <- protein[signal[1]:signal[2]]
  start_c <- length(sig) - 2
  noh <- 0
  while (noh < 2 && start_c > 1) {
    start_c <- start_c - 1
    noh <- noh + ifelse(sig[start_c] %in% c("A", "I", "L", 
                                            "F", "W", "V"), 1, -1)
    noh <- ifelse(noh < 0, 0, noh)
  }
  start_c <- start_c + 2
  start_h <- start_c - 6
  if (start_h > 1) {
    nonh <- 0
    noc <- 0
    while (nonh < 3 && noc == 0 && start_h > 1) {
      start_h <- start_h - 1
      nonh <- nonh + ifelse(sig[start_h] %in% c("A", "I", "L", 
                                                "F", "W", "V"), -1, 1)
      nonh <- ifelse(nonh < 0, 0, nonh)
      noc <- ifelse(sig[start_h] %in% c("G", "S", "T", "P", "Y", "R", "H", "K", "E"), 1, 0)
    }
  }
  else {
    start_h <- 1
  }
  prestart_c <- start_c - 1
  noh <- 0
  while (noh == 0 && start_h < prestart_c) {
    start_h <- start_h + 1
    noh <- noh + ifelse(sig[start_h] %in% c("A", "I", "L", 
                                            "F", "W", "V"), 1, 0)
  }
  c(start_n = 1, start_h = start_h, start_c = start_c, cs = signal[2])
}

calc_t2 <- function (list_prots, aa_list) {
  nhc <- t(vapply(list_prots, find_nhc2, rep(0, 4)))
  n_region <- NULL
  h_region <- NULL
  c_region <- NULL
  rest <- NULL
  for (i in 1L:length(list_prots)) {
    region_starts <- nhc[i, ]
    n_region <- c(n_region, list_prots[[i]][1:(region_starts[2] - 
                                                 1)])
    h_region <- c(h_region, list_prots[[i]][region_starts[2]:(region_starts[3] - 
                                                                1)])
    c_region <- c(c_region, list_prots[[i]][region_starts[3]:(region_starts[4] - 
                                                                1)])
    rest <- c(rest, list_prots[[i]][region_starts[4]:length(list_prots[[i]])])
  }
  t1 <- rep(0, length(aa_list))
  temp <- table(degenerate(n_region, aa_list))
  t1[as.numeric(names(temp))] <- temp
  names(t1) <- 1:length(aa_list)
  t2 <- rep(0, length(aa_list))
  temp <- table(degenerate(h_region, aa_list))
  t2[as.numeric(names(temp))] <- temp
  names(t2) <- 1:length(aa_list)
  t3 <- rep(0, length(aa_list))
  temp <- table(degenerate(c_region, aa_list))
  t3[as.numeric(names(temp))] <- temp
  names(t3) <- 1:length(aa_list)
  t4 <- rep(0, length(aa_list))
  temp <- table(degenerate(rest, aa_list))
  t4[as.numeric(names(temp))] <- temp
  names(t4) <- 1:length(aa_list)
  len_c <- nhc[, "cs"] - nhc[, "start_c"]
  len_h <- nhc[, "start_c"] - nhc[, "start_h"]
  len_n <- nhc[, "start_h"] - nhc[, "start_n"]
  lengths <- matrix(c(len_n, len_h, len_c), ncol = 3)
  colnames(lengths) <- c("n", "h", "c")
  list(mean_cs = mean(nhc[, 4]), sd_cs = sd(nhc[, 4]), t1 = t1, 
       t2 = t2, t3 = t3, t4 = t4, lengths = lengths)
}

train_hsmm2 <- function (train_data, aa_group, max_length = 32) {
  train_data <- lapply(train_data, toupper)
  ts <- calc_t2(train_data, aa_group)
  t1 <- ts[["t1"]]
  t2 <- ts[["t2"]]
  t3 <- ts[["t3"]]
  t4 <- ts[["t4"]]
  overall <- t4
  overall.probs <- overall/sum(overall)
  overall_probs_log = log(overall.probs)
  lengths <- ts[["lengths"]]
  params <- apply(lengths, 2, signalHsmm:::measure_region, max_length = max_length)
  params <- cbind(params, rep(1/max_length, max_length))
  additional_margin = 10
  pipar <- c(1, 0, 0, 0)
  tpmpar <- matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 
                     0, 0, 0), 4, byrow = TRUE)
  od <- matrix(c((t1/sum(t1))[1L:4], (t2/sum(t2))[1L:4], (t3/sum(t3))[1L:4], 
                 (t4/sum(t4))[1L:4]), 4, byrow = TRUE)
  res <- list(aa_group = aa_group, pipar = pipar, tpmpar = tpmpar, 
              od = od, overall_probs_log = overall_probs_log, params = params)
  class(res) <- "sighsmm_model"
  res
}

nokmer_model2 <- train_hsmm2(seqs, aaaggregation)
preds2 <- predict(nokmer_model2, seqs)
sum(pred2df(preds2)[["sp.probability"]] < 0.5)

real_labels <- c(rep(1, 214), rep(0, 214))

signalHSMM2010h <- train_hsmm2(read_uniprot(paste0(pathway, "sp1950_2010.txt"), euk = TRUE, what = "signal"), aaaggregation)
preds3 <- predict(nokmer_model2, benchmark_data2)
real_labels <- c(rep(1, 214), rep(0, 214))
pred_vals <- pred2df(preds3)[["sp.probability"]]

library(dplyr)
library(biogram)
library(signalHsmm)
library(seqinr)
library(hmeasure)
library(ggplot2)
library(knitr)
library(foreach)
library(doMC)


signalHSMM_proteins <- read_uniprot(paste0(pathway, "sp1950_2010.txt"), euk = TRUE, what = "signal")

                                    
b <- lapply(signalHSMM_proteins, function(x){
  cut <- attr(x, "sig")[2]
  deg_sample <- as.numeric(degenerate(toupper(x)[1L:min(70,length(x))], aaaggregation  ))
  #remove atypical amino acids
  deg_sample <- na.omit(deg_sample)
  deg_sample[(cut-2):(cut+2)]
})

b <- matrix(unlist(b), ncol=5, byrow = TRUE)
maxSignal <- 45
p4 <- 0.1
p3 <- 0.25/(1-p4)
p2 <- 0.35/(1-p4)/(1-p3)
p1 <- 0.1/(1-p4)/(1-p3)/(1-p2)

sum(b[,1]%in%c(3) & b[,2]%in%c(1) & b[,3] %in% c(2) & b[,4]%in%c(1) & b[,4]%in%c(1,3))/nrow(b)

kMer1 <- list(3, 1, 2, 1, c(1,3))
pState1 <- 3
nState1 <- 4
pTrans1 <- p1

# -3|222
sum(b[,1]%in%c(2,3) & b[,2]%in%c(1,2,3) & b[,3] %in% c(2,4))/nrow(b)
kMer2 <- list(c(2,3),c(1,2), c(2,4))
pState2 <- 3
nState2 <- 4
pTrans2 <- p2

#-3|4(1,2,3,4)(2,4)
sum(b[,1]%in%c(4) & b[,2]%in%c(1,2,3) & b[,3] %in% c(2,4))/nrow(b)
kMer3 <- list(4, c(1,2,3), c(2,4))
pState3 <- 3
nState3 <- 4
pTrans3 <- p3

# kMer4
sum(b[,1]%in%c(2) & b[,2]%in%c(1) & b[,3] %in% c(1,4) & b[,4]%in%c(1,2,3,4))/nrow(b)
kMer4 <- list(2, 1, c(1,4))
pState4 <- 3
nState4 <- 4
pTrans4 <- p4

lastState1 = 4+1+length(kMer1)

registerDoMC(cores=7)
n <- length(benchmark_data2)
results <- foreach(i=1:n) %dopar% {
  prot <- benchmark_data2[[i]]
  deg_sample <- na.omit(as.numeric(degenerate(toupper(prot)[1L:min(maxSignal, length(prot))], aaaggregation  )))
  
  viterbi_res <- duration_viterbi(deg_sample-1, signalHSMM2010h$pipar, signalHSMM2010h$tpmpar, 
                                  signalHSMM2010h$od, signalHSMM2010h$params)
  viterbi_path <- viterbi_res[["path"]]+1
  c_site <- ifelse(any(viterbi_path == 4), 
                   max(which(viterbi_path == 3)), 
                   length(deg_sample))
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  
  viterbi_res5 <- duration_viterbi(deg_sample-1, parametersSet5$pipar, parametersSet5$tpmpar, 
                                   parametersSet5$od, parametersSet5$params)
  viterbi_path5 <- viterbi_res5[["path"]]+1
  c_site5 <- ifelse(any(viterbi_path5 == 4), 
                    min(which(viterbi_path5 == 4))-1, 
                    length(deg_sample))
  c_site5 <- ifelse(any(viterbi_path5 == lastState1), 
                    min(which(viterbi_path5 == lastState1)), 
                    c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState2), 
                    min(which(viterbi_path5 == lastState2)), 
                    c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState3), 
                    min(which(viterbi_path5 == lastState3)), 
                    c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState4), 
                    min(which(viterbi_path5 == lastState4)), 
                    c_site5)
  #get probabilities of signal peptide model
  prob.signal5 <- viterbi_res5[["viterbi"]][c_site5, viterbi_path5[c_site5]]
  #get probabilities of no signal peptide model
  prob.non5 <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site5], 0)
  prob.total5 <- exp(prob.signal5 - prob.non5)
  
  c(1 - 1/(1 + prob.total), 1 - 1/(1 + prob.total5), c_site, c_site5)
}
probs <- t(sapply(results, function (x) head(x, 2)))
cut <- t(sapply(results, function (x) tail(x, 2)))
lastState2 = lastState1+1+length(kMer2)
lastState3 = lastState2+1+length(kMer3)
lastState4 = lastState3+1+length(kMer4)

parametersSet2 <- add_k_mer_state(kMer1, signalHSMM2010h$pipar, signalHSMM2010h$tpmpar, signalHSMM2010h$od, 
                                  signalHSMM2010h$params, pState1, nState1, pTrans1, d=3)
parametersSet3 <- add_k_mer_state(kMer2, pipar = parametersSet2$pipar, tpmpar = parametersSet2$tpmpar, 
                                  od = parametersSet2$od, params = parametersSet2$params,
                                  pState2, nState2, pTrans2, d=3)
parametersSet4 <- add_k_mer_state(kMer3, pipar = parametersSet3$pipar, tpmpar = parametersSet3$tpmpar, 
                                  od = parametersSet3$od, params = parametersSet3$params, 
                                  pState3, nState3, pTrans3, d=3)
parametersSet5 <- add_k_mer_state(kMer4, pipar = parametersSet4$pipar, tpmpar = parametersSet4$tpmpar, 
                                  od = parametersSet4$od, params = parametersSet4$params, 
                                  pState4, nState4, pTrans4, d=3)
overall_probs_log <- signalHSMM2010h$overall_probs_log
maxSignal <- 45


registerDoMC(cores=7)
n <- length(benchmark_data2)
results <- foreach(i=1:n) %dopar% {
  prot <- benchmark_data2[[i]]
  deg_sample <- na.omit(as.numeric(degenerate(toupper(prot)[1L:min(maxSignal, length(prot))], aaaggregation)))
  
  viterbi_res <- duration_viterbi(deg_sample-1, signalHSMM2010h$pipar, signalHSMM2010h$tpmpar, 
                                  signalHSMM2010h$od, signalHSMM2010h$params)
  viterbi_path <- viterbi_res[["path"]]+1
  c_site <- ifelse(any(viterbi_path == 4), 
                   max(which(viterbi_path == 3)), 
                   length(deg_sample))
  prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
  prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
  prob.total <- exp(prob.signal - prob.non)
  
  viterbi_res5 <- duration_viterbi(deg_sample-1, parametersSet5$pipar, parametersSet5$tpmpar, 
                                   parametersSet5$od, parametersSet5$params)
  viterbi_path5 <- viterbi_res5[["path"]]+1
  c_site5 <- ifelse(any(viterbi_path5 == 4), 
                    min(which(viterbi_path5 == 4))-1, 
                    length(deg_sample))
  c_site5 <- ifelse(any(viterbi_path5 == lastState1), 
                    min(which(viterbi_path5 == lastState1)), 
                    c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState2), 
                    min(which(viterbi_path5 == lastState2)), 
                    c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState3), 
                    min(which(viterbi_path5 == lastState3)), 
                    c_site5)
  c_site5 <- ifelse(any(viterbi_path5 == lastState4), 
                    min(which(viterbi_path5 == lastState4)), 
                    c_site5)
  #get probabilities of signal peptide model
  prob.signal5 <- viterbi_res5[["viterbi"]][c_site5, viterbi_path5[c_site5]]
  #get probabilities of no signal peptide model
  prob.non5 <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site5], 0)
  prob.total5 <- exp(prob.signal5 - prob.non5)
  
  c(1 - 1/(1 + prob.total), 1 - 1/(1 + prob.total5), c_site, c_site5)
}
probs <- t(sapply(results, function (x) head(x, 2)))
cut <- t(sapply(results, function (x) tail(x, 2)))

ggplot(data.frame(noMer=probs[,1],kMer=probs[,2], sig = factor(real_labels)), aes(x=noMer,y=kMer, color=sig)) +
  geom_point() + xlab("No k-mer") + ylab("with k-mer") + 
  scale_color_discrete("Signal peptide") +geom_abline(slope=1)

# bench_metrics2 <- rbind(bench_metrics2, HMeasure(real_labels, data.frame(signalHsmmH = pred_vals, signalHsmmHKmer = probs[,2]),
#                                threshold = c(rep(0.5, 5), 0.1, 0.1, 0.5))[["metrics"]])
# save(plot_dat, bench_metrics2, file = "ngram_extension.RData")

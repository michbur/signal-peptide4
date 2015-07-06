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
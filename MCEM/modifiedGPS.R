modifiedGPS = function (DATABASE, RR0 = 1, MIN.n11 = 1, DECISION = 1, DECISION.THRES = 0.05, 
                RANKSTAT = 1, TRONC = FALSE, TRONC.THRES = 1, PRIOR.INIT = c(alpha1 = 0.2, beta1 = 0.06, alpha2 = 1.4, beta2 = 1.8, w = 0.1), PRIOR.PARAM = NULL)
{
  DATA <- DATABASE$data
  N <- DATABASE$N
  L <- DATABASE$L
  n11 <- DATA[, 1]
  n1. <- DATA[, 2]
  n.1 <- DATA[, 3]
  E <- DATA[, 2] * DATA[, 3]/N
  P_OUT <- TRUE
  if (is.null(PRIOR.PARAM)) {
    P_OUT <- FALSE
    if (TRONC == FALSE) {
      data_cont <- xtabs(DATA[, 1] ~ L[, 1] + L[, 2])
      n1._mat <- apply(data_cont, 1, sum)
      n.1_mat <- apply(data_cont, 2, sum)
      n1._c <- rep(n1._mat, times = length(n.1_mat))
      n.1_c <- rep(n.1_mat, each = length(n1._mat))
      E_c <- n1._c * n.1_c/N
      n11_c <- as.vector(data_cont)
      p_out <- suppressWarnings(nlm(.lik2NB, p = PRIOR.INIT, 
                                    n11 = n11_c, E = E_c, iterlim = 500))
    }
    if (TRONC == TRUE) {
      tronc <- TRONC.THRES - 1
      p_out <- suppressWarnings(nlm(.likTronc2NB, p = PRIOR.INIT, 
                                    n11 = n11[n11 >= TRONC.THRES], E = E[n11 >= TRONC.THRES], 
                                    tronc, iterlim = 500))
    }
    PRIOR.PARAM <- p_out$estimate 
    code.convergence <- p_out$code
  }
  if (MIN.n11 > 1) {
    E <- E[n11 >= MIN.n11]
    n1. <- n1.[n11 >= MIN.n11]
    n.1 <- n.1[n11 >= MIN.n11]
    LL <- data.frame(drugs = L[, 1], events = L[, 2], n11)
    LL1 <- LL[, 1][n11 >= MIN.n11]
    LL2 <- LL[, 2][n11 >= MIN.n11]
    rm(list = "L")
    L <- data.frame(LL1, LL2)
    n11 <- n11[n11 >= MIN.n11]
  }
  Nb.Cell <- length(n11)
  post.H0 <- vector(length = Nb.Cell)
  Q <- PRIOR.PARAM[5] * dnbinom(n11, size = PRIOR.PARAM[1], 
                                prob = PRIOR.PARAM[2]/(PRIOR.PARAM[2] + E))/(PRIOR.PARAM[5] * 
                                                                               dnbinom(n11, size = PRIOR.PARAM[1], prob = PRIOR.PARAM[2]/(PRIOR.PARAM[2] + 
                                                                                                                                            E)) + (1 - PRIOR.PARAM[5]) * dnbinom(n11, size = PRIOR.PARAM[3], 
                                                                                                                                                                                 prob = PRIOR.PARAM[4]/(PRIOR.PARAM[4] + E)))
  post.H0 <- Q * pgamma(RR0, PRIOR.PARAM[1] + n11, PRIOR.PARAM[2] + 
                          E) + (1 - Q) * pgamma(RR0, PRIOR.PARAM[3] + n11, PRIOR.PARAM[4] + 
                                                  E)
  # Posterior Expectation of log2(lambda) formula 10 in the paper
  postE <- log(2)^(-1) * (Q * (digamma(PRIOR.PARAM[1] + n11) - 
                                 log(PRIOR.PARAM[2] + E)) + (1 - Q) * (digamma(PRIOR.PARAM[3] + 
                                                                                 n11) - log(PRIOR.PARAM[4] + E)))
  LB <- .QuantileDuMouchel(0.05, Q, PRIOR.PARAM[1] + n11, PRIOR.PARAM[2] + 
                             E, PRIOR.PARAM[3] + n11, PRIOR.PARAM[4] + E)
  
  alpha1 <- PRIOR.PARAM[1] + n11
  beta1 <- PRIOR.PARAM[2] + E
  alpha2 <- PRIOR.PARAM[3] + n11
  beta2 <- PRIOR.PARAM[4] + E
  var <- Q*(1-Q)*(alpha1/beta1-alpha2/beta2)^2+Q*alpha1/(beta1^2)+(1-Q)*alpha2/(beta2^2)
  var <- var*1.0/N # s.e.
  
  if (RANKSTAT == 1) 
    RankStat <- post.H0
  if (RANKSTAT == 2) 
    RankStat <- LB
  if (RANKSTAT == 3) 
    RankStat <- postE
  if (RANKSTAT == 1) {
    FDR <- (cumsum(post.H0[order(RankStat)])/(1:length(post.H0)))
    FNR <- rev(cumsum((1 - post.H0)[order(1 - RankStat)]))/(Nb.Cell - 
                                                              1:length(post.H0))
    Se <- cumsum((1 - post.H0)[order(RankStat)])/(sum(1 - 
                                                        post.H0))
    Sp <- rev(cumsum(post.H0[order(1 - RankStat)]))/(Nb.Cell - 
                                                       sum(1 - post.H0))
  }
  if (RANKSTAT == 2 | RANKSTAT == 3) {
    FDR <- (cumsum(post.H0[order(RankStat, decreasing = TRUE)])/(1:length(post.H0)))
    FNR <- rev(cumsum((1 - post.H0)[order(1 - RankStat, decreasing = TRUE)]))/(Nb.Cell - 
                                                                                 1:length(post.H0))
    Se <- cumsum((1 - post.H0)[order(RankStat, decreasing = TRUE)])/(sum(1 - 
                                                                           post.H0))
    Sp <- rev(cumsum(post.H0[order(1 - RankStat, decreasing = TRUE)]))/(Nb.Cell - 
                                                                          sum(1 - post.H0))
  }
  if (DECISION == 1) 
    Nb.signaux <- sum(FDR <= DECISION.THRES)
  if (DECISION == 2) 
    Nb.signaux <- min(DECISION.THRES, Nb.Cell)
  if (DECISION == 3) {
    if (RANKSTAT == 1) 
      Nb.signaux <- sum(RankStat <= DECISION.THRES, na.rm = TRUE)
    if (RANKSTAT == 2 | RANKSTAT == 3) 
      Nb.signaux <- sum(RankStat >= DECISION.THRES, na.rm = TRUE)
  }
  Q_func <- Q
  Expect_Q <- postE
  RES <- vector(mode = "list")
  RES$INPUT.PARAM <- data.frame(RR0, MIN.n11, DECISION, DECISION.THRES, 
                                RANKSTAT, TRONC, TRONC.THRES)
  # RES$STOPPING <- data.frame(Q_func, Expect_Q)
  RES$PARAM <- vector(mode = "list")
  if (P_OUT == TRUE) 
    RES$PARAM$PRIOR.PARAM <- data.frame(PRIOR.PARAM)
  if (P_OUT == FALSE) {
    RES$PARAM$PRIOR.INIT <- data.frame(PRIOR.INIT)
    RES$PARAM$PRIOR.PARAM <- PRIOR.PARAM
    RES$PARAM$CONVERGENCE <- code.convergence
  }
  if (RANKSTAT == 1) {
    RES$ALLSIGNALS <- data.frame(L[, 1][order(RankStat)], 
                                 L[, 2][order(RankStat)], n11[order(RankStat)], E[order(RankStat)], 
                                 RankStat[order(RankStat)], (n11/E)[order(RankStat)], 
                                 n1.[order(RankStat)], n.1[order(RankStat)], FDR, 
                                 FNR, Se, Sp,var)
    colnames(RES$ALLSIGNALS) <- c("drug", "event", "count", 
                                  "expected count", "postH0", "n11/E", "drug margin", 
                                  "event margin", "FDR", "FNR", "Se", "Sp","var")
  }
  if (RANKSTAT == 2 | RANKSTAT == 3) {
    RES$ALLSIGNALS <- data.frame(L[, 1][order(RankStat, decreasing = TRUE)], 
                                 L[, 2][order(RankStat, decreasing = TRUE)], n11[order(RankStat, 
                                                                                       decreasing = TRUE)], E[order(RankStat, decreasing = TRUE)], 
                                 RankStat[order(RankStat, decreasing = TRUE)], (n11/E)[order(RankStat, 
                                                                                             decreasing = TRUE)], n1.[order(RankStat, decreasing = TRUE)], 
                                 n.1[order(RankStat, decreasing = TRUE)], FDR, FNR, 
                                 Se, Sp, post.H0[order(RankStat, decreasing = TRUE)], var)
    if (RANKSTAT == 2) 
      colnames(RES$ALLSIGNALS) <- c("drug", "event", "count", 
                                    "expected count", "Q_0.05(lambda)", "n11/E", 
                                    "drug margin", "event margin", "FDR", "FNR", 
                                    "Se", "Sp", "postH0", "var")
    if (RANKSTAT == 3) 
      colnames(RES$ALLSIGNALS) <- c("drug", "event", "count", 
                                    "expected count", "post E(Lambda)", "n11/E", 
                                    "drug margin", "event margin", "FDR", "FNR", 
                                    "Se", "Sp", "postH0", "var")
  }
  RES$ALLSIGNALS$Q_func = Q_func; RES$ALLSIGNALS$Expect_Q = Expect_Q
  RES$SIGNALS <- RES$ALLSIGNALS[1:Nb.signaux, ]
  RES$NB.SIGNALS <- Nb.signaux
  return(RES)
}

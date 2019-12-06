library(PhViD)
library(dplyr)
source("modifiedGPS.R")
MCEM = function(ADR, Drug, RankStat = "Q_0.05(lambda)", MaxIter = 100)
{
  # Create count table
  CT = Count_table(ADR, Drug)
  CT.count = as.PhViD(CT$Count_table)
  # initiate GPS
  GPS.Result = modifiedGPS(CT.count, RANKSTAT = 2)
  Signals = GPS.Result$ALLSIGNALS[,c(RankStat, "Se")]
  Prior = GPS.Result$PARAM$PRIOR.PARAM
  # first iteration
  # -- MC STEP
  NextLongTable = MCStep(CT$Long_table, Signals)
  NextLongTable = rbind(CT$Long_table, NextLongTable)
  # -- Max Step
  NextGPSResult = MaxStep(NextLongTable, Prior)
  NextGPSResultRetained = retainGPS(GPS.Result, NextGPSResult)
  # iterate
  epsilon = 10^-6 # stopping criterion
  for (i in 2:MaxIter){
    PrevGPSResult = NextGPSResult
    Signals = NextGPSResult$ALLSIGNALS[,c(RankStat, "Se")]
    PrevLongTable = NextLongTable
    PrevGPSResultRetained = NextGPSResultRetained
    # -- MC STEP
    NextLongTable = MCStep(PrevLongTable, Signals)
    # NextLongTable = rbind(PrevLongTable, NextLongTable)
    # -- Max Step
    NextGPSResult = MaxStep(NextLongTable, Prior)
    NextGPSResultRetained = retainGPS(PrevGPSResultRetained, NextGPSResult)
    # -- check convergence
    PrevQ = PrevGPSResultRetained$ALLSIGNALS$Q_func
    NextQ = NextGPSResultRetained$ALLSIGNALS$Q_func
    cat(i, " ", mean(abs(NextQ-PrevQ)), "\n")
    if (mean(abs(NextQ-PrevQ)) < epsilon){
      Convergence = TRUE
      break
    }else{
      Convergence = FALSE
    }
  }
  return(list(GPS = NextGPSResultRetained, Convergence = Convergence, iters = i))
}

#' Convert binary data to count table
#' @param Drug binary matrix of drugs
#' @param ADR binary matrix of ADR
library(dplyr)
Count_table = function(ADR, Drug){
  # reverse-crosstab ADR
  ADR_Long = data.frame(id = 1:nrow(ADR), stack(ADR))
  ADR_Long = ADR_Long %>% 
    filter(values == 1) %>%
    rename(ADR = ind) %>%
    dplyr::select(id, ADR)
  # reverse-crosstab Drug
  Drug_Long = data.frame(id = 1:nrow(Drug), stack(Drug))
  Drug_Long = Drug_Long %>%
    filter(values == 1)%>%
    rename(Drug = ind) %>%
    dplyr::select(id, Drug) 
  # join
  combined = inner_join(ADR_Long, Drug_Long, by = "id")
  
  Count_table = combined %>%
    group_by(Drug, ADR) %>%
    summarise(Count = n())
  return(list(Count_table = as.data.frame(Count_table), Long_table = as.data.frame(combined)))
}


MCStep = function(Long_table, Signals){
  AllCases = Long_table %>% mutate(DrugADR = paste(Drug, ADR),
                                   lambda = ifelse(is.na(Signals[DrugADR, 1]), 0, Signals[DrugADR, 1]))
  SelectedCases = AllCases %>%
    group_by(id, ADR) %>%
    summarise(Drug = Drug[which(rmultinom(1, 1, lambda)>0)])
  return(as.data.frame(SelectedCases))
}
MaxStep = function(Long_table, prior){
  ContingencyTable = Long_table %>%
    group_by(Drug, ADR) %>%
    summarise(Count = n()) %>%
    filter(Count > 0) %>%
    as.data.frame()
  PhViDdata <- as.PhViD(ContingencyTable)
  GPS.Result = modifiedGPS(PhViDdata, RANKSTAT = 2, PRIOR.INIT = prior)
  return(GPS.Result)
}

retainGPS = function(PrevGPSResult, NextGPSResult){
  Prev = PrevGPSResult$ALLSIGNALS
  Next = NextGPSResult$ALLSIGNALS
  joined = left_join(  
    Prev %>%
      dplyr::select(drug, event),
    Next,
    by = c("drug", "event"))
  joined[is.na(joined)] = 0
  NextGPSResult$ALLSIGNALS = joined
  return(NextGPSResult)
}

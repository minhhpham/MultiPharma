sample_data = read.csv("sample input.csv")
Drug = sample_data[,2:4]
ADR=  sample_data[,5:7]

set.seed(1)
MCEM_result = MCEM(ADR, Drug)
MCEM_result$GPS

RGPS_result = RGPS(ADR1 ~ DRUG1 + DRUG2 + DRUG3, data = sample_data)
RGPS_result

# MultiPharma
MCEM publication:
  https://www.nature.com/articles/s41598-018-19979-7
  
RGPS publication:
  https://pdfs.semanticscholar.org/d90f/698d4c6aefeb039bc92826414ba28fc0248a.pdf

Example Usage:
```
sample_data = read.csv("sample input.csv")
Drug = sample_data[,2:4]
ADR=  sample_data[,5:7]

set.seed(1)
MCEM_result = MCEM(ADR, Drug)
MCEM_result$GPS

set.seed(1)
RGPS_result = RGPS(ADR1 ~ DRUG1 + DRUG2 + DRUG3, data = sample_data)
RGPS_result
```

Both `MCEM_result$GPS` and `RGPS_result` have the same object class as returned by `PhViD::GPS()`. Please refer to `PhViD` documentation to understand the results.

# Future Developments
I can make this repository an R package when I can find free time in my schedule. If more people show interest in using my work, I will give this a higher priority. For any question, please feel free to contact or open an issue on Github.
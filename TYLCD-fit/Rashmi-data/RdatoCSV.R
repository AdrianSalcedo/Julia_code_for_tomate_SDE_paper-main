library(dplyr)
library(tidyverse)
dataFrame <- data.frame(df_fit_CIS)
dataFrame[,1:3] <- dataFrame[,1:3]/1000
dataFrameReal <- data.frame(df_sample_N_real)
dataFrameReal[,1] <- dataFrameReal[,1]/1000
write_csv(dataFrameReal,"/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model6/estimation/Rashmi_data/Rashmi_real_data.csv")
write_csv(dataFrame,"/home/gabrielsalcedo/Documentos/TYLCVD-fit/Model6/estimation/Rashmi_data/NegBin_Rashmi_data.csv")

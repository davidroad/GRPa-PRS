# newmetaSWE should include at least the prs you calculated and the diagnosis indicator. 

subgroup_newmetaSWE <- newmetaSWE %>% 
  mutate(
    
    prsSCZ = case_when(case_control == 1 ~ prs),
    prsCtr = case_when(case_control == 0 ~ prs),
    
    #20
    TB20all = case_when( prs > quantile( prs, 0.775) &  prs < quantile( prs, 0.975) ~ 1,
                         prs < quantile( prs, 0.225) &  prs > quantile( prs, 0.025) ~ 0),
    TB20SCZ = case_when( prsSCZ > quantile( prsSCZ[!is.na( prsSCZ)], 0.775) &  prsSCZ < quantile( prsSCZ[!is.na( prsSCZ)], 0.975)~ 1,
                         prsSCZ < quantile( prsSCZ[!is.na( prsSCZ)], 0.225) &  prsSCZ > quantile( prsSCZ[!is.na( prsSCZ)], 0.025) ~ 0),
    TB20Ctr = case_when( prsCtr > quantile( prsCtr[!is.na( prsCtr)], 0.775) &  prsCtr < quantile( prsCtr[!is.na( prsCtr)], 0.975)~ 1,
                         prsCtr < quantile( prsCtr[!is.na( prsCtr)], 0.225) &  prsCtr > quantile( prsCtr[!is.na( prsCtr)], 0.025) ~ 0),
    
    
  )

#  
prsSCZ = sort(subgroup_newmetaSWE$prsSCZ[!is.na(subgroup_newmetaSWE$prsSCZ)])
prsCtr = sort(subgroup_newmetaSWE$prsCtr[!is.na(subgroup_newmetaSWE$prsCtr)])
sprs = sort(c(prsSCZ,prsCtr))

BL = quantile(prsSCZ, 0.025)
BU20 = sprs[which(abs(sprs - BL) == min(abs(sprs - BL))) + round(length(sprs)*0.20)]

TU = quantile(prsCtr, 0.975)
TL20 = sprs[which(abs(sprs - TU) == min(abs(sprs - TU))) - round(length(sprs)*0.20)]


subgroup_newmetaSWE = subgroup_newmetaSWE %>% 
  mutate(# 20
         T20 = case_when(prs >= TL20 & prs <= TU & case_control == 1 ~ 1,
                         prs >= TL20 & prs <= TU & case_control == 0 ~ 0),
         B20 = case_when(prs >= BL & prs <= BU20 & case_control == 1 ~ 1,
                         prs >= BL & prs <= BU20 & case_control == 0 ~ 0)
  )


write.table(subgroup_newmetaSWE, file = "/data2/xli37/data/TWAS/subgrp_SWEmeta0217.txt", sep = "\t",col.names = TRUE, row.names = FALSE, quote = FALSE)


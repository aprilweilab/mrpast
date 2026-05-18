//Number of population samples (demes)  
5 
//Population effective sizes (number of genes)  
AFR_N0  
EUR_N0
ASIA_N0
NAT_N0
ADMIX_N0
//Sample sizes  
200
200
200
200
200
//Growth rates : negative growth implies population expansion  
0
0 
0 
0 
0 
//Number of migration matrices : 0 implies no migration between demes  
0 
//historical event: time, source, sink, migrants, new size, new growth rate, migr. matrix  
9 historical event  
T1 4 0 AM1 1 0 0
T1 4 1 AMPAR2 1 EUR_NG 0
T1 2 2 0 ASIA_N1 ASIA_NG 0
T1 4 3 1 1 0 0
T2 1 1 0 EUR_N1 0 0
T2 3 2 1 ASIA_N1 0 0
T3 2 1 1 EUR_N2 0 0
T4 1 0 1 1 0 0  
T5 0 0 0 AFR_N1 0 0  
//Number of independent loci [chromosome]  
1 0  
//Per chromosome: Number of linkage blocks  
1 
//per Block: data type, num loci, rec. rate and mut rate + optional parameters  
FREQ 1 1e-8 2.35e-8

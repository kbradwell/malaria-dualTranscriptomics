1. use malipaxgene_de_v2 to obtain differential expression over timepoints. Change line 79 to qlf <- glmQLFTest(fit,coef=2:3) to test differential expression across individuals
2. to conduct a similar test as in point #1, while also including a covariate such as parasite developmental stage, see malipaxgene_de_v2_covariate

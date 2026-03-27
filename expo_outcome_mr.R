#### MR ####
library(VariantAnnotation)
library(gwasglue)
library(TwoSampleMR)
load("ebi-a-GCST90020053exposure_dat.rda")
exposure_dat$exposure="Frailty index"
dat_list_res=exposure_dat
dat_list_res1=dat_list_res[which(dat_list_res$pval.exposure<5e-8),]
dat=dat_list_res1
N=dat[1,"samplesize.exposure"]    
dat=transform(dat,R2=2*((beta.exposure)^2)*eaf.exposure*(1-eaf.exposure))
dat=transform(dat,R2=((beta.exposure)^2/((beta.exposure)^2+((se.exposure)^2)*N))) #计算R2
dat_list_res1=dat
save(dat_list_res1,file = "expo_list_final.F.rda")
dat_list_res1=dat_list_res1[which(dat_list_res1$F>10),]
table(dat_list_res1$id.exposure)
length(unique(dat_list_res1$id.exposure))
exp_dat=dat_list_res1
exposure_idall1 <- unique(exp_dat$id.exposure)
exposure_idall=exposure_idall1
outcomeID <- c('G6_SLEEPAPNO')
outcomeName <- c('Sleep apnoea')
dat_list <- list()
res_single_list=list()
res_single_OR_list=list()
for (exposureID in exposure_idall) {
 
  exposure_dat=exp_dat[which(exp_dat$id.exposure==exposureID),]
  outcomeTab<-merge(exposure_dat, outcomeData1, by.x="SNP", by.y="SNP")
  if(dim(outcomeTab)[1]>2){
    write.csv(outcomeTab[,-(2:ncol(exposure_dat))], file="outcome.csv")
    outcome_dat<-read_outcome_data(snps=exposure_dat$SNP,
                                   filename="outcome.csv", sep = ",")
    outcome_dat$outcome=outcomeName
    outcome_dat$id.outcome=outcomeID
    dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exposure_dat,
      outcome_dat = outcome_dat
    )
    res <- mr(dat)
    if(dim(res)[1]>2){
      res_single_list[[exposureID]] <- res
      OR <-generate_odds_ratios(res)
      res_single_OR_list[[exposureID]] <- OR
      dat_list[[exposureID]] <- dat
    }
  }
}
save(dat_list,file = "exposure_outcome_harmonise_all.rda")
save(res_single_list,file = "exposure_outcome_res_list.rda")
save(res_single_OR_list,file = "exposure_outcome_res_OR_list.rda")
res_single_list1=do.call(rbind,res_single_list)
res_single_OR_list1=do.call(rbind,res_single_OR_list)
dat_list1=do.call(rbind,dat_list)
write.csv(dat_list1,"exposure_outcome_harmonise.csv",row.names = FALSE)
write.csv(res_single_list1,"result_MR.csv",row.names = FALSE)
write.csv(res_single_OR_list1,"result_MR_OR.csv",row.names = FALSE)

id=unique(res_single_list1$id.exposure)
exposure_idall=id
library(MRPRESSO)
pleio_list <- list()
hete_list <- list()
presso_list=list()
presso_list_main=list()
for (exposureID in exposure_idall) {
  dat <- dat_list[[exposureID]]
  het <- mr_heterogeneity(dat)
  hete_list[[exposureID]] <- het
  pleio <- mr_pleiotropy_test(dat)
  if(dim(dat)[1]>3){
    presso <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                        OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat, NbDistribution = 1000,  
                        SignifThreshold = 0.05)
    presso$`Main MR results`$result=exposureID
    presso_list[[exposureID]]=presso
    presso_list_main[[exposureID]]=presso$`Main MR results`
    
    pleio$PRESSO_P <- presso$`MR-PRESSO results`$`Global Test`$Pvalue
    pleio$RSSobs <- presso$`MR-PRESSO results`$`Global Test`$RSSobs
  }else{
    pleio$PRESSO_P <-""
    pleio$RSSobs <-""
  }
  pleio_list[[exposureID]] <- pleio
}


pleio_merge <- do.call(rbind,pleio_list)
hete_merge <- do.call(rbind,hete_list)
presso_list_main_merge=do.call(rbind,presso_list_main)
write.csv(pleio_merge,"pleiotropy_test.csv",row.names = FALSE)
write.csv(hete_merge,"heterogeneity_test.csv",row.names = FALSE)
write.csv(presso_list_main_merge,"mr_presso_main.csv")
save(presso_list,file = "mr_presso_result.rda")

for (exposureID in exposure_idall){

  library(patchwork)
  library(ggplot2)
  single <- mr_leaveoneout(dat=dat_list[[exposureID]])
  write.csv(single,paste0("leaveone_",exposureID,".csv"),row.names = FALSE)
  pa <- mr_leaveoneout_plot(single)[[1]]

  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'none')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste(dat_list[[exposureID]]$exposure,paste0(' on ',dat_list[[exposureID]]$outcome),sep = "\n"))
  pdf(paste0("leaveone_",exposureID,".pdf"))
  print(pa)
  dev.off()

  
  
  pa <- mr_scatter_plot(dat=dat_list[[exposureID]],mr_results = res_single_list[[exposureID]])[[1]]
  pa$labels$x <- unlist(strsplit(pa$labels$x,split = ' \\|\\| '))[1]

  pa$labels$y <- paste0("SNP effect on ",dat_list[[exposureID]]$outcome)
  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'top')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste0("SNP effect on ",dat_list[[exposureID]]$exposure))
  pdf(paste0("scatter_",exposureID,".pdf"),height = 8,width = 8)
  print(pa)
  dev.off()

  
  res_single_2 <- mr_singlesnp(dat=dat_list[[exposureID]])
  write.csv(res_single_2,paste0("forest_",exposureID,".csv"),row.names = FALSE)
  pa <- mr_forest_plot(res_single_2)[[1]]

  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'none')+
    update_geom_defaults("line", list(size = 5))+
    labs(x=paste(dat_list[[exposureID]]$exposure,paste0(' on ',dat_list[[exposureID]]$outcome),sep = "\n"))
  pdf(paste0("forest_",exposureID,".pdf"),width = 8)
  print(pa)
  dev.off()
  
  res_single_2 <- mr_singlesnp(dat=dat_list[[exposureID]])
  nrow(res_single_2)
  res_single_2=res_single_2[-nrow(res_single_2),]
  pa <- mr_funnel_plot(res_single_2)[[1]]

  pa <- pa+scale_color_brewer(palette = 'Set1')+
    theme_classic(base_size = 18)+
    theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
          axis.ticks = element_line(size = 1),
          legend.position = 'top')+
    update_geom_defaults("line", list(size = 5))
  pdf(paste0("funnel_",exposureID,".pdf"),height = 8,width = 8)
  print(pa)
  dev.off()
  
}





load("../../czf/1.shenrentangniaobing/medium_effect/cis-eQTLs.pvalue0.05-ProbeLevel_new.rda")

s=which(snp.p$SNP%in%dat_list1$SNP)
snp.p1=snp.p[s,]
write.csv(snp.p1,"exposure_snp.gene.p5e-8.csv")
table(snp.p1$GeneSymbol)

expo.snp=read.csv("old/result/ebi-a-GCST90020053/leaveone_ebi-a-GCST90020053.csv",header = T,stringsAsFactors = F,check.names = F)
s=which(snp.p$SNP%in%expo.snp$SNP)
snp.p1=snp.p[s,]
write.csv(snp.p1,"exposure_snp.gene.p5e-6.csv")
length(unique(snp.p1$GeneSymbol))









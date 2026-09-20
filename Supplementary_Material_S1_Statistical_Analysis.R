# ==============================================================================
# Supplementary R script
# In vitro germination of polyembryonic mango (Mangifera indica L.) cv. 'Chulucanas'
# Tables 1-3 | Figures 1-5 | Model diagnostics | Data exports
# ==============================================================================

# 1. SETUP ---------------------------------------------------------------------
rm(list=ls()); graphics.off(); gc(); set.seed(123)
options(contrasts=c("contr.sum","contr.poly"),scipen=999)

required_packages <- c("googlesheets4","tidyverse","janitor","car","glmmTMB","DHARMa",
                       "emmeans","brglm2","multcomp","multcompView","FactoMineR","Hmisc",
                       "ComplexHeatmap","circlize","cowplot","scales","writexl")
missing_packages <- required_packages[!vapply(required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)>0) stop("Install required packages before running the script: ",paste(missing_packages,collapse=", "))
invisible(lapply(required_packages,library,character.only=TRUE))

dir.create("Results_FINAL",showWarnings=FALSE,recursive=TRUE)
FONT <- "serif"
treatment_levels <- paste0("T",1:12)
treatment_colors <- c(T1="#1B9E77",T2="#D95F02",T3="#7570B3",T4="#E7298A",T5="#66A61E",T6="#E6AB02",
                      T7="#A6761D",T8="#666666",T9="#1F78B4",T10="#B2DF8A",T11="#FB9A99",T12="#CAB2D6")

theme_article <- function(base_size=12){
  ggplot2::theme_bw(base_size=base_size,base_family=FONT)+
    ggplot2::theme(text=ggplot2::element_text(family=FONT,color="black"),
                   axis.title=ggplot2::element_text(size=base_size+1),
                   axis.text=ggplot2::element_text(size=base_size-1,color="black"),
                   axis.text.x=ggplot2::element_text(angle=0,hjust=.5),
                   panel.border=ggplot2::element_rect(color="black",fill=NA,linewidth=.7),
                   panel.grid.major=ggplot2::element_line(linewidth=.3,color="grey85"),
                   panel.grid.minor=ggplot2::element_blank(),
                   legend.title=ggplot2::element_text(face="bold"),
                   strip.text=ggplot2::element_text(face="bold"),
                   plot.margin=ggplot2::margin(8,8,8,8))
}
ggplot2::theme_set(theme_article())

# 2. DATA IMPORT AND CLEANING ---------------------------------------------------
SHEET_URL <- "https://docs.google.com/spreadsheets/d/14KitmqBSxjkpdRTOGGMeaQZM_E0FU8PFmquQUVeLj3Q/edit"
googlesheets4::gs4_deauth()
data_raw <- googlesheets4::read_sheet(ss=SHEET_URL,sheet="Hoja 1",show_col_types=FALSE)
if(ncol(data_raw)<17) stop("The input dataset must contain at least 17 columns.")

data_seed <- data_raw %>% janitor::clean_names()
names(data_seed)[1:17] <- c("ID","Description","Rep","Treatment_raw","Seed","Fruit_length","Fruit_width","Fruit_weight",
                            "TSS","Firmness","Seed_weight","Embryos","Shoots","Leaves","Phenolics","Contamination","Establishment")
numeric_vars <- c("Fruit_length","Fruit_width","Fruit_weight","TSS","Firmness","Seed_weight","Embryos","Shoots","Leaves",
                  "Phenolics","Contamination","Establishment")
data_seed <- data_seed %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(numeric_vars),~suppressWarnings(as.numeric(unlist(.x,use.names=FALSE)))),
                Treatment_num=as.numeric(stringr::str_extract(as.character(Treatment_raw),"\\d+")),
                Treatment=factor(paste0("T",Treatment_num),levels=treatment_levels),
                Maturity=factor(ifelse(Treatment_num<=6,"PM","CM"),levels=c("PM","CM")),
                CW=factor(dplyr::case_when(Treatment_num %in% c(1,2,3,7,8,9)~"0",Treatment_num %in% c(4,5,6,10,11,12)~"20"),levels=c("0","20")),
                LJ=factor(dplyr::case_when(Treatment_num %in% c(1,4,7,10)~"0",Treatment_num %in% c(2,5,8,11)~"1.5",
                                           Treatment_num %in% c(3,6,9,12)~"3"),levels=c("0","1.5","3")),
                Rep_ID=interaction(Treatment,Rep,drop=TRUE),Phenolics_yes=as.integer(Phenolics>0),
                Contamination_yes=as.integer(Contamination>0),Establishment_yes=as.integer(Establishment>0))
if(anyNA(data_seed$Treatment)) stop("Treatment coding could not be resolved for all observations.")
if(!all(table(data_seed$Treatment)==20)) stop("Each treatment must contain 20 seeds.")

rep_data <- data_seed %>%
  dplyr::group_by(Treatment,Maturity,CW,LJ,Rep,Rep_ID) %>%
  dplyr::summarise(n_seed=dplyr::n(),Fruit_length=mean(Fruit_length,na.rm=TRUE),Fruit_width=mean(Fruit_width,na.rm=TRUE),
                   Fruit_weight=mean(Fruit_weight,na.rm=TRUE),TSS=mean(TSS,na.rm=TRUE),Firmness=mean(Firmness,na.rm=TRUE),
                   Seed_weight=mean(Seed_weight,na.rm=TRUE),Embryos_total=sum(Embryos,na.rm=TRUE),Shoots_total=sum(Shoots,na.rm=TRUE),
                   Embryos_mean=mean(Embryos,na.rm=TRUE),Shoots_mean=mean(Shoots,na.rm=TRUE),Positive_seed_n=sum(Shoots>0,na.rm=TRUE),
                   Leaves_conditional=ifelse(any(Shoots>0,na.rm=TRUE),mean(Leaves[Shoots>0],na.rm=TRUE),NA_real_),
                   Phenolics_yes=sum(Phenolics_yes,na.rm=TRUE),Contamination_yes=sum(Contamination_yes,na.rm=TRUE),
                   Establishment_yes=sum(Establishment_yes,na.rm=TRUE),Phenolics_prop=Phenolics_yes/n_seed,
                   Contamination_prop=Contamination_yes/n_seed,Establishment_prop=Establishment_yes/n_seed,.groups="drop")
if(nrow(rep_data)!=48) stop("Expected 48 experimental replicates (12 treatments x 4 replicates).")

# 3. TABLE 1: FRUIT AND SEED CHARACTERIZATION ----------------------------------
char_vars <- c(Fruit_length="Fruit length (mm)",Fruit_width="Fruit width (mm)",Fruit_weight="Fruit weight (g)",
               TSS="Total soluble solids (°Brix)",Firmness="Fruit firmness (kg cm⁻²)",Seed_weight="Seed weight without endocarp (g)")
char_abbr <- c(Fruit_length="FL",Fruit_width="FW",Fruit_weight="FWT",TSS="TSS",Firmness="FF",Seed_weight="SW")
Table1 <- dplyr::bind_rows(lapply(names(char_vars),function(v){
  df <- data_seed %>% dplyr::select(Maturity,dplyr::all_of(v)) %>% dplyr::filter(!is.na(.data[[v]]))
  tt <- stats::t.test(df[[v]]~df$Maturity,var.equal=FALSE)
  sm <- df %>% dplyr::group_by(Maturity) %>% dplyr::summarise(Mean=mean(.data[[v]],na.rm=TRUE),SE=sd(.data[[v]],na.rm=TRUE)/sqrt(dplyr::n()),.groups="drop")
  tibble::tibble(Variable=char_vars[[v]],Abbreviation=char_abbr[[v]],
                 PM=sprintf("%.2f ± %.2f",sm$Mean[sm$Maturity=="PM"],sm$SE[sm$Maturity=="PM"]),
                 CM=sprintf("%.2f ± %.2f",sm$Mean[sm$Maturity=="CM"],sm$SE[sm$Maturity=="CM"]),p_value=tt$p.value)
})) %>% dplyr::mutate(p=dplyr::if_else(p_value<.001,"<0.001",sprintf("%.3f",p_value))) %>%
  dplyr::select(Variable,Abbreviation,PM,CM,p)

# 4. TABLE 2A: FACTORIAL MODELS ------------------------------------------------
m_DE <- glmmTMB::glmmTMB(Embryos_total~Maturity*CW*LJ,family=glmmTMB::compois(link="log"),data=rep_data)
m_DP <- glmmTMB::glmmTMB(Shoots_total~Maturity*CW*LJ,family=stats::poisson(link="log"),data=rep_data)
fit_binary <- function(success){
  f <- stats::as.formula(paste0("cbind(",success,",n_seed-",success,")~Maturity*CW*LJ"))
  stats::glm(f,family=stats::binomial("logit"),data=rep_data,method=brglm2::brglmFit,type="AS_mean")
}
m_PE <- fit_binary("Phenolics_yes"); m_MC <- fit_binary("Contamination_yes"); m_IVE <- fit_binary("Establishment_yes")
anova_DE <- car::Anova(m_DE,type=3); anova_DP <- car::Anova(m_DP,type=3)
anova_PE <- car::Anova(m_PE,type=3); anova_MC <- car::Anova(m_MC,type=3); anova_IVE <- car::Anova(m_IVE,type=3)
effect_order <- c("Maturity","CW","LJ","Maturity:CW","Maturity:LJ","CW:LJ","Maturity:CW:LJ")
effect_labels <- c("Maturity stage (MS)","Coconut water (CW)","Lemon juice (LJ)","MS × CW","MS × LJ","CW × LJ","MS × CW × LJ")
extract_wald <- function(x){
  d <- as.data.frame(x); d$Effect <- rownames(d); sc <- grep("Chisq",names(d),value=TRUE)[1]; pc <- grep("Pr",names(d),value=TRUE)[1]
  d %>% dplyr::transmute(Effect,df=Df,Statistic=.data[[sc]],p=.data[[pc]]) %>% dplyr::filter(Effect %in% effect_order) %>%
    dplyr::mutate(Effect=factor(Effect,levels=effect_order)) %>% dplyr::arrange(Effect)
}
format_wald <- function(x) x %>% dplyr::mutate(Result=paste0(sprintf("%.3f",Statistic)," (",ifelse(p<.001,"<0.001",sprintf("%.3f",p)),")",
                                                               dplyr::case_when(p<.001~"***",p<.01~"**",p<.05~"*",TRUE~"")))
DEf <- format_wald(extract_wald(anova_DE)); DPf <- format_wald(extract_wald(anova_DP)); PEf <- format_wald(extract_wald(anova_PE))
MCf <- format_wald(extract_wald(anova_MC)); IVEf <- format_wald(extract_wald(anova_IVE))
Table2A <- tibble::tibble(Effect=effect_labels,df=DEf$df,Developed_embryos=DEf$Result,Developed_plumules=DPf$Result,
                          Phenolic_exudation=PEf$Result,Microbial_contamination=MCf$Result,In_vitro_establishment=IVEf$Result)

# 5. TABLE 2B: CONDITIONAL LEAF DEVELOPMENT ------------------------------------
leaf_data <- rep_data %>% dplyr::filter(!is.na(Leaves_conditional)) %>% droplevels()
m_LDS <- stats::lm(Leaves_conditional~Maturity+CW+LJ,data=leaf_data)
anova_LDS <- car::Anova(m_LDS,type=2); shapiro_LDS <- stats::shapiro.test(stats::residuals(m_LDS))
levene_LDS <- car::leveneTest(Leaves_conditional~Treatment,data=leaf_data,center=median)
LDS_aov <- as.data.frame(anova_LDS)
Table2B <- tibble::tibble(Effect=c("Maturity stage (MS)","Coconut water (CW)","Lemon juice (LJ)"),df=LDS_aov$Df[1:3],
                          F=round(LDS_aov$`F value`[1:3],3),p=round(LDS_aov$`Pr(>F)`[1:3],3))
emm_LDS_MS <- emmeans::emmeans(m_LDS,~Maturity); emm_LDS_CW <- emmeans::emmeans(m_LDS,~CW); emm_LDS_LJ <- emmeans::emmeans(m_LDS,~LJ)

# 6. TABLE 3: MODEL-ESTIMATED PE, MC AND IVE -----------------------------------
add_treatment <- function(df){
  df %>% dplyr::mutate(Treatment=dplyr::case_when(
    Maturity=="PM"&CW=="0"&LJ=="0"~"T1",Maturity=="PM"&CW=="0"&LJ=="1.5"~"T2",Maturity=="PM"&CW=="0"&LJ=="3"~"T3",
    Maturity=="PM"&CW=="20"&LJ=="0"~"T4",Maturity=="PM"&CW=="20"&LJ=="1.5"~"T5",Maturity=="PM"&CW=="20"&LJ=="3"~"T6",
    Maturity=="CM"&CW=="0"&LJ=="0"~"T7",Maturity=="CM"&CW=="0"&LJ=="1.5"~"T8",Maturity=="CM"&CW=="0"&LJ=="3"~"T9",
    Maturity=="CM"&CW=="20"&LJ=="0"~"T10",Maturity=="CM"&CW=="20"&LJ=="1.5"~"T11",Maturity=="CM"&CW=="20"&LJ=="3"~"T12"),
    Treatment=factor(Treatment,levels=treatment_levels))
}
get_binary_emm <- function(model){
  e <- emmeans::emmeans(model,~Maturity*CW*LJ,type="response")
  z <- as.data.frame(multcomp::cld(e,adjust="tukey",Letters=letters,reversed=TRUE)) %>% add_treatment() %>%
    dplyr::mutate(Mean=prob,Lower_CI=asymp.LCL,Upper_CI=asymp.UCL,Letters=stringr::str_trim(.group)) %>% dplyr::arrange(Treatment)
  if(anyNA(z$Mean)) stop("NA detected in binary-model estimated marginal means."); z
}
PE_emm <- get_binary_emm(m_PE); MC_emm <- get_binary_emm(m_MC); IVE_emm <- get_binary_emm(m_IVE)
fmt_prob <- function(x,name) x %>% dplyr::transmute(Treatment,!!name:=sprintf("%.1f (%.1f–%.1f)",Mean*100,Lower_CI*100,Upper_CI*100))
Treatment_information <- rep_data %>% dplyr::distinct(Treatment,Maturity,CW,LJ) %>% dplyr::arrange(Treatment)
Table3 <- Treatment_information %>% dplyr::left_join(fmt_prob(PE_emm,"PE"),by="Treatment") %>%
  dplyr::left_join(fmt_prob(MC_emm,"MC"),by="Treatment") %>% dplyr::left_join(fmt_prob(IVE_emm,"IVE"),by="Treatment")

# 7. MODEL DIAGNOSTICS ----------------------------------------------------------
set.seed(123)
sim_DE <- DHARMa::simulateResiduals(m_DE,plot=FALSE); sim_DP <- DHARMa::simulateResiduals(m_DP,plot=FALSE)
Diagnostics <- tibble::tibble(Model=c("DE_COMPoisson","DP_Poisson"),
  Dispersion_p=c(DHARMa::testDispersion(sim_DE)$p.value,DHARMa::testDispersion(sim_DP)$p.value),
  Zero_inflation_p=c(DHARMa::testZeroInflation(sim_DE)$p.value,DHARMa::testZeroInflation(sim_DP)$p.value))

# 8. FIGURE 1: PE, MC AND IVE --------------------------------------------------
draw_binary <- function(df,y_label){
  ymax <- max(1.08,max(df$Upper_CI,na.rm=TRUE)+.12)
  ggplot2::ggplot(df,ggplot2::aes(Treatment,Mean,fill=Treatment))+
    ggplot2::geom_col(width=.65,color="black",linewidth=.5,alpha=.8)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower_CI,ymax=Upper_CI),width=.18,linewidth=.4)+
    ggplot2::geom_text(ggplot2::aes(y=Upper_CI+.018,label=Letters),family=FONT,size=4,vjust=0)+
    ggplot2::scale_fill_manual(values=treatment_colors,drop=FALSE)+ggplot2::scale_x_discrete(drop=FALSE)+
    ggplot2::scale_y_continuous(labels=scales::percent_format(accuracy=1),breaks=c(0,.25,.50,.75,1),limits=c(0,ymax),expand=c(0,0))+
    ggplot2::labs(x="Treatment",y=y_label)+ggplot2::guides(fill="none")+theme_article(12)+ggplot2::coord_cartesian(clip="off")
}
Fig1A <- draw_binary(PE_emm,"PE (%)")+ggplot2::labs(tag="A"); Fig1B <- draw_binary(MC_emm,"MC (%)")+ggplot2::labs(tag="B")
Fig1C <- draw_binary(IVE_emm,"IVE (%)")+ggplot2::labs(tag="C")
Figure1 <- cowplot::plot_grid(cowplot::plot_grid(Fig1A,Fig1B,ncol=2,align="hv",axis="tblr"),
                              cowplot::plot_grid(NULL,Fig1C,NULL,ncol=3,rel_widths=c(.5,1,.5)),ncol=1)
ggplot2::ggsave("Results_FINAL/Figure1_PE_MC_IVE.jpg",Figure1,width=12,height=9,units="in",dpi=600,bg="white")

# 9. FIGURE 2: DE, DP AND LDS --------------------------------------------------
mean_se <- function(data,variable) data %>% dplyr::group_by(Treatment) %>%
  dplyr::summarise(n=sum(!is.na(.data[[variable]])),Mean=mean(.data[[variable]],na.rm=TRUE),
                   SE=stats::sd(.data[[variable]],na.rm=TRUE)/sqrt(n),.groups="drop")
DE_summary <- mean_se(rep_data,"Embryos_mean"); DP_summary <- mean_se(rep_data,"Shoots_mean")
LDS_summary <- rep_data %>% dplyr::group_by(Treatment) %>% dplyr::summarise(n=sum(!is.na(Leaves_conditional)),
  Mean=ifelse(n>0,mean(Leaves_conditional,na.rm=TRUE),NA_real_),SE=ifelse(n>1,sd(Leaves_conditional,na.rm=TRUE)/sqrt(n),NA_real_),.groups="drop")
draw_mean_se <- function(df,y_label){
  ymax <- max(df$Mean+ifelse(is.na(df$SE),0,df$SE),na.rm=TRUE)*1.25
  ggplot2::ggplot(df,ggplot2::aes(Treatment,Mean,fill=Treatment))+ggplot2::geom_col(width=.65,color="black",linewidth=.5,alpha=.85)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=pmax(0,Mean-SE),ymax=Mean+SE),width=.18,linewidth=.4)+
    ggplot2::geom_text(ggplot2::aes(y=Mean+SE+ymax*.04,label=sprintf("%.2f",Mean)),family=FONT,size=3.8)+
    ggplot2::scale_fill_manual(values=treatment_colors,drop=FALSE)+ggplot2::scale_x_discrete(drop=FALSE)+
    ggplot2::scale_y_continuous(limits=c(0,ymax),expand=c(0,0),labels=function(x)sprintf("%.2f",x))+
    ggplot2::labs(x="Treatment",y=y_label)+ggplot2::guides(fill="none")+theme_article(12)+ggplot2::coord_cartesian(clip="off")
}
Fig2A <- draw_mean_se(DE_summary,"DE")+ggplot2::labs(tag="A"); Fig2B <- draw_mean_se(DP_summary,"DP")+ggplot2::labs(tag="B")
LDS_ymax <- max(LDS_summary$Mean+ifelse(is.na(LDS_summary$SE),0,LDS_summary$SE),na.rm=TRUE)*1.25
Fig2C <- ggplot2::ggplot(LDS_summary,ggplot2::aes(Treatment,Mean,fill=Treatment))+
  ggplot2::geom_col(data=LDS_summary %>% dplyr::filter(!is.na(Mean)),width=.65,color="black",linewidth=.5,alpha=.85)+
  ggplot2::geom_errorbar(data=LDS_summary %>% dplyr::filter(!is.na(SE)),ggplot2::aes(ymin=pmax(0,Mean-SE),ymax=Mean+SE),width=.18,linewidth=.4)+
  ggplot2::geom_text(data=LDS_summary %>% dplyr::filter(!is.na(Mean)),ggplot2::aes(y=Mean+ifelse(is.na(SE),0,SE)+LDS_ymax*.04,label=sprintf("%.2f",Mean)),family=FONT,size=3.8)+
  ggplot2::geom_text(data=LDS_summary %>% dplyr::filter(is.na(Mean)),ggplot2::aes(y=LDS_ymax*.06,label="NE"),family=FONT,fontface="italic",size=3.8)+
  ggplot2::scale_fill_manual(values=treatment_colors,drop=FALSE)+ggplot2::scale_x_discrete(drop=FALSE)+
  ggplot2::scale_y_continuous(limits=c(0,LDS_ymax),expand=c(0,0),labels=function(x)sprintf("%.2f",x))+
  ggplot2::labs(x="Treatment",y="LDS",tag="C")+ggplot2::guides(fill="none")+theme_article(12)+ggplot2::coord_cartesian(clip="off")
Figure2 <- cowplot::plot_grid(cowplot::plot_grid(Fig2A,Fig2B,ncol=2,align="hv",axis="tblr"),
                              cowplot::plot_grid(NULL,Fig2C,NULL,ncol=3,rel_widths=c(.5,1,.5)),ncol=1)
ggplot2::ggsave("Results_FINAL/Figure2_DE_DP_LDS.jpg",Figure2,width=12,height=9,units="in",dpi=600,bg="white")

# 10. FIGURE 3: SPEARMAN CORRELATION MATRICES ----------------------------------
corr_vars <- c("Fruit_length","Fruit_width","Fruit_weight","TSS","Firmness","Seed_weight","Embryos_mean","Shoots_mean",
               "Leaves_conditional","Phenolics_prop","Contamination_prop","Establishment_prop")
pretty_names_corr <- c("FL","FW","FWT","TSS","FF","SW","DE","DP","LDS","PE","MC","IVE")
var_order <- pretty_names_corr; y_levels_plot <- rev(var_order); x_levels_plot <- var_order[-length(var_order)]
compute_corr_long <- function(df){
  corr_df <- df %>% dplyr::select(dplyr::all_of(corr_vars)); colnames(corr_df) <- pretty_names_corr
  res <- Hmisc::rcorr(as.matrix(corr_df),type="spearman")
  r_long <- as.data.frame(res$r) %>% tibble::rownames_to_column("Var1") %>% tidyr::pivot_longer(-Var1,names_to="Var2",values_to="r")
  p_long <- as.data.frame(res$P) %>% tibble::rownames_to_column("Var1") %>% tidyr::pivot_longer(-Var1,names_to="Var2",values_to="p")
  dplyr::left_join(r_long,p_long,by=c("Var1","Var2")) %>% dplyr::mutate(Var1=factor(Var1,levels=var_order),Var2=factor(Var2,levels=var_order)) %>%
    dplyr::filter(as.numeric(Var2)<as.numeric(Var1)) %>% dplyr::mutate(sig=dplyr::case_when(p<.001~"***",p<.01~"**",p<.05~"*",TRUE~""),
      label=paste0(sprintf("%.2f",r),sig),Var1=factor(Var1,levels=y_levels_plot),Var2=factor(Var2,levels=x_levels_plot))
}
build_corr_plot <- function(corr_long,tag_label){
  ggplot2::ggplot(corr_long,ggplot2::aes(Var2,Var1,fill=r))+ggplot2::geom_tile(color="grey90",linewidth=.4)+
    ggplot2::geom_text(ggplot2::aes(label=label),family=FONT,size=3.5)+
    ggplot2::scale_fill_gradient2(low="firebrick3",mid="white",high="steelblue4",midpoint=0,limits=c(-1,1),breaks=seq(-1,1,.5),name=expression(rho))+
    ggplot2::scale_x_discrete(position="bottom",drop=FALSE)+ggplot2::scale_y_discrete(drop=FALSE)+ggplot2::coord_fixed(clip="off")+
    ggplot2::labs(x=NULL,y=NULL,tag=tag_label)+ggplot2::theme_minimal(base_family=FONT)+
    ggplot2::theme(panel.grid=ggplot2::element_blank(),axis.text=ggplot2::element_text(color="black",size=10),plot.tag=ggplot2::element_text(size=14,face="bold"))
}
Fig3A <- build_corr_plot(compute_corr_long(rep_data %>% dplyr::filter(Maturity=="PM")),"A")
Fig3B <- build_corr_plot(compute_corr_long(rep_data %>% dplyr::filter(Maturity=="CM")),"B")
Figure3 <- cowplot::plot_grid(Fig3A,Fig3B,ncol=1,align="v",axis="lr")
ggplot2::ggsave("Results_FINAL/Figure3_Spearman.jpg",Figure3,width=13,height=20,units="in",dpi=600,bg="white")

# 11. FIGURE 4: PCA -------------------------------------------------------------
pca_vars <- c("Fruit_length","Fruit_width","Fruit_weight","TSS","Firmness","Seed_weight","Embryos_mean","Shoots_mean",
              "Phenolics_prop","Contamination_prop","Establishment_prop")
pretty_names_pca <- c("FL","FW","FWT","TSS","FF","SW","DE","DP","PE","MC","IVE")
pca_df_clean <- rep_data %>% dplyr::select(dplyr::all_of(c(pca_vars,"Treatment","Rep_ID"))) %>% stats::na.omit()
pca_matrix <- pca_df_clean %>% dplyr::select(dplyr::all_of(pca_vars)); colnames(pca_matrix) <- pretty_names_pca
res_pca <- FactoMineR::PCA(pca_matrix,scale.unit=TRUE,graph=FALSE)
PCA_loadings <- as.data.frame(res_pca$var$coord); PCA_contributions <- as.data.frame(res_pca$var$contrib)
ind_coords <- as.data.frame(res_pca$ind$coord[,1:2]); colnames(ind_coords) <- c("Dim1","Dim2")
ind_coords$Treatment <- factor(pca_df_clean$Treatment,levels=treatment_levels)
var_coords <- as.data.frame(res_pca$var$coord[,1:2]); colnames(var_coords) <- c("Dim1","Dim2"); var_coords$Variable <- rownames(var_coords)
pct_dim1 <- round(res_pca$eig[1,2],1); pct_dim2 <- round(res_pca$eig[2,2],1)
var_coords <- var_coords %>% dplyr::mutate(Dim1=Dim1*3,Dim2=Dim2*3,label_x=Dim1+sign(Dim1)*.15,label_y=Dim2+sign(Dim2)*.15,
  label_hjust=dplyr::case_when(Dim1>0~0,Dim1<0~1,TRUE~.5),label_vjust=dplyr::case_when(Dim2>0~0,Dim2<0~1,TRUE~.5))
max_val <- max(abs(c(ind_coords$Dim1,ind_coords$Dim2,var_coords$Dim1,var_coords$Dim2))); panel_lim <- ceiling(max_val)+.3
Figure4 <- ggplot2::ggplot()+ggplot2::geom_hline(yintercept=0,linetype="dashed",linewidth=.4)+ggplot2::geom_vline(xintercept=0,linetype="dashed",linewidth=.4)+
  ggplot2::geom_point(data=ind_coords,ggplot2::aes(Dim1,Dim2,color=Treatment),size=2.5,alpha=.85)+
  ggplot2::geom_segment(data=var_coords,ggplot2::aes(x=0,y=0,xend=Dim1,yend=Dim2),arrow=grid::arrow(length=grid::unit(.2,"cm")),linewidth=.6)+
  ggplot2::geom_text(data=var_coords,ggplot2::aes(label_x,label_y,label=Variable,hjust=label_hjust,vjust=label_vjust),family=FONT,size=4)+
  ggplot2::scale_color_manual(values=treatment_colors,limits=treatment_levels,drop=FALSE)+ggplot2::coord_cartesian(xlim=c(-panel_lim,panel_lim),ylim=c(-panel_lim,panel_lim),clip="off")+
  ggplot2::labs(x=paste0("PC1 (",pct_dim1,"%)"),y=paste0("PC2 (",pct_dim2,"%)"),color="Treatment")+theme_article(12)+
  ggplot2::theme(legend.position="bottom",legend.direction="horizontal")+ggplot2::guides(color=ggplot2::guide_legend(nrow=2,byrow=TRUE))
ggplot2::ggsave("Results_FINAL/Figure4_PCA.jpg",Figure4,width=11,height=9.5,units="in",dpi=600,bg="white",limitsize=FALSE)

# 12. FIGURE 5: Z-SCORE HEATMAP ------------------------------------------------
heat_raw <- rep_data %>% dplyr::group_by(Treatment) %>% dplyr::summarise(FL=mean(Fruit_length,na.rm=TRUE),FW=mean(Fruit_width,na.rm=TRUE),
  FWT=mean(Fruit_weight,na.rm=TRUE),TSS=mean(TSS,na.rm=TRUE),FF=mean(Firmness,na.rm=TRUE),SW=mean(Seed_weight,na.rm=TRUE),
  DE=mean(Embryos_mean,na.rm=TRUE),DP=mean(Shoots_mean,na.rm=TRUE),PE=mean(Phenolics_prop,na.rm=TRUE)*100,
  MC=mean(Contamination_prop,na.rm=TRUE)*100,IVE=mean(Establishment_prop,na.rm=TRUE)*100,.groups="drop") %>% tibble::column_to_rownames("Treatment")
if(anyNA(heat_raw)) stop("Missing values remain in heatmap data.")
heat_z_matrix <- scale(heat_raw); if(any(!is.finite(heat_z_matrix))) stop("Non-finite Z-scores detected.")
z_text_matrix <- matrix(sprintf("%.2f",heat_z_matrix),nrow=nrow(heat_z_matrix),ncol=ncol(heat_z_matrix),dimnames=dimnames(heat_z_matrix))
ht <- ComplexHeatmap::Heatmap(heat_z_matrix,name="Z-score",col=grDevices::colorRampPalette(c("#B2182B","#F7F7F7","#2166AC"))(100),
  cluster_rows=TRUE,cluster_columns=TRUE,rect_gp=grid::gpar(col="white",lwd=1.2),
  heatmap_legend_param=list(title="Z-score",title_position="leftcenter-rot",at=c(-4,-2,0,2,4),labels=c("−4","−2","0","2","4")),
  cell_fun=function(j,i,x,y,w,h,fill) grid::grid.text(z_text_matrix[i,j],x,y,gp=grid::gpar(fontsize=9,fontfamily=FONT,col="black")),
  row_names_gp=grid::gpar(fontsize=10,fontfamily=FONT),column_names_gp=grid::gpar(fontsize=10,fontfamily=FONT),column_names_rot=0,border="white")
grDevices::jpeg("Results_FINAL/Figure5_Heatmap.jpg",width=11,height=9,units="in",res=600,quality=95)
ComplexHeatmap::draw(ht,heatmap_legend_side="left",padding=grid::unit(c(2,2,2,2),"mm")); grDevices::dev.off()

# 13. EXPORT TABLES AND REPRODUCIBILITY INFORMATION -----------------------------
writexl::write_xlsx(list(Table1=Table1,Table2A=Table2A,Table2B=Table2B,Table3=Table3,Diagnostics=Diagnostics,
                         PCA_loadings=PCA_loadings,PCA_contributions=PCA_contributions),"Results_FINAL/Statistical_results.xlsx")
utils::write.csv(rep_data,"Results_FINAL/Replicate_level_data.csv",row.names=FALSE,fileEncoding="UTF-8")
utils::write.csv(PE_emm,"Results_FINAL/PE_estimated_marginal_means.csv",row.names=FALSE,fileEncoding="UTF-8")
utils::write.csv(MC_emm,"Results_FINAL/MC_estimated_marginal_means.csv",row.names=FALSE,fileEncoding="UTF-8")
utils::write.csv(IVE_emm,"Results_FINAL/IVE_estimated_marginal_means.csv",row.names=FALSE,fileEncoding="UTF-8")
writeLines(capture.output(sessionInfo()),"Results_FINAL/sessionInfo.txt")
cat("\nAnalysis completed successfully. Outputs saved in 'Results_FINAL'.\n")

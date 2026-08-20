#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survcomp)
  library(jsonlite)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: analyze_luad_figure5.R FIG4_ROOT PDXE_TSV_DIR OUTPUT_ROOT")
fig4 <- normalizePath(args[[1]], mustWork = TRUE)
pdxe_dir <- normalizePath(args[[2]], mustWork = TRUE)
out <- normalizePath(args[[3]], mustWork = TRUE)
tables <- file.path(out, "results", "tables")
figures <- file.path(out, "results", "figures")
audit <- file.path(out, "audit")
dir.create(tables, recursive = TRUE, showWarnings = FALSE)
dir.create(figures, recursive = TRUE, showWarnings = FALSE)
dir.create(audit, recursive = TRUE, showWarnings = FALSE)
set.seed(20260820)

norm_key <- function(x) toupper(gsub("[^A-Z0-9]", "", fifelse(is.na(x), "", as.character(x))))
read_activity <- function(path) {
  x <- fread(path, check.names = FALSE); ids <- x[[1]]; x[[1]] <- NULL
  m <- as.matrix(x); storage.mode(m) <- "double"; rownames(m) <- ids; m
}
safe_ci <- function(x, outcome_time) {
  keep <- is.finite(x) & is.finite(outcome_time)
  x <- x[keep]; outcome_time <- outcome_time[keep]
  if (length(x) < 10L || sd(x) == 0 || sd(outcome_time) == 0) {
    return(list(c.index=NA_real_, p.value=NA_real_, lower=NA_real_, upper=NA_real_))
  }
  tryCatch({
    fit <- survcomp::concordance.index(
      x=x, surv.time=outcome_time, surv.event=rep(1, length(x)),
      method="noether", na.rm=TRUE
    )
    list(c.index=unname(fit$c.index), p.value=unname(fit$p.value),
         lower=unname(fit$lower), upper=unname(fit$upper))
  }, error=function(e) list(c.index=NA_real_, p.value=NA_real_, lower=NA_real_, upper=NA_real_))
}

inputs <- readRDS(file.path(out, "data", "processed", "figure5_cellline_inputs.rds"))
activity <- inputs$activity
hc <- inputs$hc_tfs
if (length(hc) != 31L || !all(hc %in% rownames(activity))) stop("Frozen 31 HC-TF input mismatch")

test_one_drug <- function(dataset, response) {
  response <- response[, .(aac=median(aac), drug_name=first(na.omit(drug_name)),
                           max_aac=max(aac)), by=.(canonical_drug, depmap_id)]
  response <- response[depmap_id %in% colnames(activity)]
  if (uniqueN(response$depmap_id) < 10L || max(response$aac) <= 0.2) return(NULL)
  tf_results <- rbindlist(lapply(hc, function(tf) {
    values <- activity[tf, response$depmap_id]
    fit <- safe_ci(values, 1-response$aac)
    data.table(dataset=dataset, canonical_drug=response$canonical_drug[[1]],
               drug_name=response$drug_name[[1]], TF=tf, n=nrow(response),
               max_AAC=max(response$aac), c_index=fit$c.index,
               p_value=fit$p.value, lower95=fit$lower, upper95=fit$upper,
               spearman_rho=unname(cor(values,response$aac,method="spearman")))
  }))
  tf_results
}

tasks <- unlist(lapply(names(inputs$datasets), function(dataset) {
  response <- copy(inputs$datasets[[dataset]]$paired)
  response <- response[eligible_luad == TRUE & depmap_id %in% colnames(activity)]
  split(response, by="canonical_drug", keep.by=TRUE)
}), recursive=FALSE)
task_dataset <- unlist(lapply(names(inputs$datasets), function(dataset) {
  response <- inputs$datasets[[dataset]]$paired
  response <- response[eligible_luad == TRUE & depmap_id %in% colnames(activity)]
  rep(dataset, length(unique(response$canonical_drug)))
}), use.names=FALSE)
if (length(tasks) != length(task_dataset)) stop("Cell-line task construction mismatch")

cores <- min(12L, max(1L, detectCores()-2L))
cell_results <- mclapply(seq_along(tasks), function(i) test_one_drug(task_dataset[[i]], as.data.table(tasks[[i]])), mc.cores=cores)
cell_results <- rbindlist(cell_results, fill=TRUE)
cell_results <- cell_results[is.finite(c_index) & is.finite(p_value)]
cell_results[, FDR := p.adjust(p_value, method="BH"), by=dataset]
cell_results[, direction := fifelse(c_index > 0.5, "higher TF activity -> sensitivity",
                                    "higher TF activity -> resistance")]
cell_results[, significant := FDR <= 0.05]
cell_results[, effect_abs := abs(c_index-0.5)]
setorder(cell_results, dataset, FDR, -effect_abs)
fwrite(cell_results, file.path(tables, "Figure5ABC_cellline_concordance_all.tsv.gz"), sep="\t")

sig <- cell_results[significant == TRUE]
pair_repl <- sig[, .(
  datasets=paste(sort(unique(dataset)),collapse=";"),
  dataset_n=uniqueN(dataset),
  direction_n=uniqueN(direction),
  direction=first(direction),
  min_FDR=min(FDR),
  mean_c_index=mean(c_index),
  drug_name=first(na.omit(drug_name))
), by=.(canonical_drug,TF)]
pair_repl[, replicated := dataset_n >= 2L & direction_n == 1L]
fwrite(pair_repl, file.path(tables, "Figure5D_replicated_drug_TF_pairs.tsv"), sep="\t")

membership <- sig[, .(membership=paste(sort(unique(dataset)),collapse=" + ")),
                  by=.(canonical_drug,TF)]
membership_counts <- membership[, .N, by=membership][order(-N,membership)]
fwrite(membership_counts, file.path(tables, "Figure5D_significant_pair_membership.tsv"), sep="\t")

cox <- fread(file.path(fig4, "results", "tables", "Figure4_all_Cox_models.tsv"))
multi <- cox[model == "multivariable" & p_value <= 0.05]
prog <- multi[, .(
  prognosis_calls=paste(sort(unique(fifelse(HR < 1,"favorable","adverse"))),collapse=";"),
  significant_endpoint_n=.N,
  best_prognostic_p=min(p_value)
), by=TF]
replicated <- pair_repl[replicated == TRUE]
tf_integration <- replicated[, .(
  sensitivity_pairs=sum(direction == "higher TF activity -> sensitivity"),
  resistance_pairs=sum(direction == "higher TF activity -> resistance"),
  replicated_pairs=.N
), by=TF]
tf_integration <- merge(data.table(TF=hc), tf_integration, by="TF", all.x=TRUE)
for (column in c("sensitivity_pairs","resistance_pairs","replicated_pairs")) set(tf_integration, which(is.na(tf_integration[[column]])), column, 0L)
tf_integration <- merge(tf_integration, prog, by="TF", all.x=TRUE)
tf_integration[is.na(prognosis_calls), prognosis_calls := "not nominally prognostic"]
fwrite(tf_integration, file.path(tables, "Figure5EF_prognosis_drug_integration.tsv"), sep="\t")

# PDXE v2: the source labels lung models as NSCLC. Only the independently
# classified LUAD-like subset is used for the anchor-style response validation.
pdxe_scores <- fread(file.path(tables, "Figure5_PDXE_v2_lung_state_scores.tsv"))
pdxe_activity <- read_activity(file.path(out,"results","pdxe_projection","PDXE_V2_viper_activity.tsv.gz"))
luad_like <- pdxe_scores[prediction == "LUAD-like", sample_id]
model <- fread(file.path(pdxe_dir,"model.tsv"), na.strings=c("","NA"))
response <- fread(file.path(pdxe_dir,"sensitivity_model.tsv"), na.strings=c("","NA"))
pdx <- merge(model[tissue == "NSCLC"], response, by="model.id", all.x=FALSE)
pdx <- pdx[patient.id %in% luad_like & is.finite(best.average.response)]
pdx <- pdx[!grepl("\\+",drug) & !grepl("control|vehicle|untreated",drug,ignore.case=TRUE)]

cell_drug_names <- rbindlist(lapply(inputs$datasets, function(z) z$drugs[,.(canonical_drug,drug_name)]), fill=TRUE)
cell_drug_names[, drug_key := norm_key(drug_name)]
name_map <- cell_drug_names[drug_key != "", .(
  canonical_n=uniqueN(canonical_drug), canonical_drug=first(canonical_drug),
  cell_drug_name=first(na.omit(drug_name))
), by=drug_key][canonical_n == 1L]
pdx[, drug_key := norm_key(drug)]
pdx <- merge(pdx, name_map[,.(drug_key,canonical_drug,cell_drug_name)], by="drug_key", all.x=TRUE)
pdx <- pdx[, .(
  best_average_response=median(best.average.response),
  mRECIST={value <- na.omit(mRECIST); if(length(value)) first(value) else NA_character_},
  model_id=first(model.id)
), by=.(patient_id=patient.id,drug,drug_key,canonical_drug)]
coverage <- pdx[, .(PDX_n=uniqueN(patient_id), mapped_to_cellline_drug=sum(!is.na(canonical_drug))), by=drug][order(-PDX_n)]
fwrite(coverage, file.path(tables,"SupplementaryFigure9A_PDXE_drug_coverage.tsv"), sep="\t")

pdx[, PDX_n := .N, by=drug]
pdx_tasks <- split(pdx[PDX_n >= 10L], by="drug")
pdx_tests <- rbindlist(lapply(pdx_tasks, function(z) {
  rbindlist(lapply(hc, function(tf) {
    values <- pdxe_activity[tf,z$patient_id]
    fit <- safe_ci(values,z$best_average_response)
    data.table(drug=z$drug[[1]], drug_key=z$drug_key[[1]], canonical_drug=z$canonical_drug[[1]],
               TF=tf,n=nrow(z),c_index=fit$c.index,p_value=fit$p.value,
               lower95=fit$lower,upper95=fit$upper,
               spearman_sensitivity_rho=unname(cor(values,-z$best_average_response,method="spearman")))
  }))
}), fill=TRUE)
if (nrow(pdx_tests)) {
  pdx_tests <- pdx_tests[is.finite(c_index) & is.finite(p_value)]
  pdx_tests[, FDR := p.adjust(p_value,method="BH")]
  pdx_tests[, direction := fifelse(c_index > 0.5,"higher TF activity -> sensitivity","higher TF activity -> resistance")]
  pdx_tests <- merge(pdx_tests, replicated[,.(canonical_drug,TF,cellline_direction=direction,dataset_n)],
                     by=c("canonical_drug","TF"),all.x=TRUE)
  pdx_tests[, cellline_replicated_same_direction := !is.na(dataset_n) & direction == cellline_direction]
} else {
  pdx_tests <- data.table(drug=character(),drug_key=character(),canonical_drug=character(),TF=character(),
                          n=integer(),c_index=numeric(),p_value=numeric(),lower95=numeric(),upper95=numeric(),
                          spearman_sensitivity_rho=numeric(),FDR=numeric(),direction=character(),
                          cellline_direction=character(),dataset_n=integer(),cellline_replicated_same_direction=logical())
}
fwrite(pdx_tests, file.path(tables,"SupplementaryFigure9B_PDXE_all_concordance_tests.tsv"), sep="\t")
validated <- pdx_tests[cellline_replicated_same_direction == TRUE & FDR <= 0.05 & n >= 10L]
setorder(validated,FDR,-n,drug,TF)
selected <- if (nrow(validated)) validated[1] else validated
fwrite(validated,file.path(tables,"Figure5GH_PDX_validated_pairs.tsv"),sep="\t")

theme_paper <- theme_bw(base_size=10) + theme(panel.grid.minor=element_blank(),plot.title=element_text(face="bold",size=11))
save_plot <- function(plot,name,width=6.3,height=5.0) {
  ggsave(file.path(figures,paste0(name,".pdf")),plot,width=width,height=height,device=cairo_pdf)
  ggsave(file.path(figures,paste0(name,".png")),plot,width=width,height=height,dpi=220)
}
for (i in seq_along(c("GDSC2","CTRPv2","PRISM"))) {
  dataset_name <- c("GDSC2","CTRPv2","PRISM")[[i]]
  panel <- c("A","B","C")[[i]]
  d <- cell_results[dataset == dataset_name]
  d[, plot_FDR := pmax(FDR,1e-300)]
  d[, class := fifelse(FDR <= 0.05 & c_index > 0.5,"Sensitivity",
                       fifelse(FDR <= 0.05,"Resistance","Not significant"))]
  labels <- head(d[FDR <= 0.05][order(FDR)],8)
  p <- ggplot(d,aes(c_index-0.5,-log10(plot_FDR),color=class)) +
    geom_point(alpha=.55,size=.8) + geom_vline(xintercept=0,linetype=2,color="grey45") +
    geom_hline(yintercept=-log10(.05),linetype=3,color="grey45") +
    geom_text(data=labels,aes(label=paste0(drug_name,"–",TF)),size=2.2,check_overlap=TRUE,vjust=-.5,show.legend=FALSE) +
    scale_color_manual(values=c(Sensitivity="#B2182B",Resistance="#2166AC",`Not significant`="#BDBDBD")) +
    labs(title=sprintf("%s  %s LUAD drug–TF associations",panel,dataset_name),x="Concordance index − 0.5",y="−log10(FDR)",color=NULL) + theme_paper
  save_plot(p,paste0("Figure5",panel,"_",dataset_name,"_volcano"))
}

pD <- ggplot(membership_counts,aes(x=reorder(membership,N),y=N,fill=grepl("\\+",membership))) +
  geom_col(width=.72) + coord_flip() +
  scale_fill_manual(values=c(`TRUE`="#6A3D9A",`FALSE`="#8DA0CB"),guide="none") +
  labs(title="D  Cross-dataset membership of significant drug–TF pairs",x=NULL,y="Unique drug–TF pairs") + theme_paper
save_plot(pD,"Figure5D_significant_pair_membership",6.3,4.7)

if (nrow(replicated)) {
  top_pairs <- replicated[order(min_FDR)][1:min(.N,30),paste0(drug_name," | ",TF)]
  heat <- copy(cell_results[paste0(drug_name," | ",TF) %in% top_pairs])
  heat[, pair := paste0(drug_name," | ",TF)]
  heat[, pair := factor(pair,levels=rev(unique(top_pairs)))]
  pE <- ggplot(heat,aes(dataset,pair,fill=c_index-0.5)) + geom_tile(color="white") +
    scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-.5,.5),oob=scales::squish) +
    labs(title="E  Replicated drug–TF associations",x=NULL,y=NULL,fill="CI − 0.5") + theme_paper +
    theme(axis.text.y=element_text(size=6.5))
} else {
  pE <- ggplot()+annotate("text",x=0,y=0,label="No drug–TF pair replicated at FDR ≤ 0.05\nin at least two datasets",size=5)+
    xlim(-1,1)+ylim(-1,1)+labs(title="E  Replicated drug–TF associations")+theme_void()
}
save_plot(pE,"Figure5E_replicated_association_heatmap",7.0,6.1)

plot_tf <- tf_integration[replicated_pairs>0][order(replicated_pairs)][1:min(.N,20)]
if (nrow(plot_tf)) {
  long <- melt(plot_tf,id.vars=c("TF","prognosis_calls"),measure.vars=c("sensitivity_pairs","resistance_pairs"),
               variable.name="association",value.name="pairs")
  long[association=="resistance_pairs",pairs:=-pairs]
  pF <- ggplot(long,aes(reorder(TF,abs(pairs),FUN=sum),pairs,fill=association)) + geom_col() + coord_flip() +
    scale_fill_manual(values=c(sensitivity_pairs="#B2182B",resistance_pairs="#2166AC"),
                      labels=c("Sensitivity","Resistance")) +
    geom_text(data=unique(plot_tf[,.(TF,prognosis_calls)]),aes(x=TF,y=0,label=ifelse(prognosis_calls=="not nominally prognostic","",prognosis_calls)),
              inherit.aes=FALSE,size=2.2,hjust=-.1) +
    labs(title="F  HC-TFs linked to replicated treatment response",x=NULL,y="Replicated pairs (resistance shown negative)",fill=NULL)+theme_paper
} else {
  pF <- ggplot()+annotate("text",x=0,y=0,label="No cross-dataset replicated associations",size=5)+xlim(-1,1)+ylim(-1,1)+
    labs(title="F  HC-TFs linked to treatment response")+theme_void()
}
save_plot(pF,"Figure5F_TF_drug_integration",6.5,5.4)

if (nrow(selected)) {
  z <- pdx[drug == selected$drug]
  z[, TF_activity := pdxe_activity[selected$TF,patient_id]]
  z[, activity_group := ifelse(TF_activity >= median(TF_activity),"High TF activity","Low TF activity")]
  pG <- ggplot(z,aes(TF_activity,best_average_response)) + geom_point(size=2.5,aes(color=activity_group)) +
    geom_smooth(method="lm",se=TRUE,color="#333333",linewidth=.7) +
    scale_color_manual(values=c(`High TF activity`="#B2182B",`Low TF activity`="#2166AC")) +
    labs(title=sprintf("G  %s activity and %s response",selected$TF,selected$drug),
         subtitle=sprintf("LUAD-like PDX n=%d; CI=%.2f; FDR=%.3g",selected$n,selected$c_index,selected$FDR),
         x=paste(selected$TF,"VIPER activity"),y="Best average response (lower = greater sensitivity)",color=NULL)+theme_paper
  z <- z[order(best_average_response)]; z[,patient_id:=factor(patient_id,levels=patient_id)]
  pH <- ggplot(z,aes(patient_id,best_average_response,fill=activity_group)) + geom_col() +
    coord_flip() + scale_fill_manual(values=c(`High TF activity`="#B2182B",`Low TF activity`="#2166AC")) +
    labs(title=sprintf("H  %s response across LUAD-like PDXs",selected$drug),x=NULL,y="Best average response",fill=NULL)+theme_paper
} else {
  message_text <- "No cell-line-replicated drug–TF pair reached PDX FDR ≤ 0.05\nwith at least 10 computationally LUAD-like PDXs"
  pG <- ggplot()+annotate("text",x=0,y=0,label=message_text,size=4.5)+xlim(-1,1)+ylim(-1,1)+
    labs(title="G  Prespecified PDX validation")+theme_void()
  pH <- ggplot()+annotate("text",x=0,y=0,label="PDX waterfall plot not generated because\nthe validation criterion was not met",size=4.5)+
    xlim(-1,1)+ylim(-1,1)+labs(title="H  PDX response distribution")+theme_void()
}
save_plot(pG,"Figure5G_PDX_validation",6.4,5.0)
save_plot(pH,"Figure5H_PDX_waterfall",6.4,5.0)

# Supplementary audit panels.
dataset_audit <- fread(file.path(tables,"Figure5_dataset_audit.tsv"))
pS8A <- ggplot(dataset_audit,aes(dataset,mapped_LUAD_lines,fill=dataset))+geom_col()+geom_text(aes(label=mapped_LUAD_lines),vjust=-.3)+
  labs(title="A  LUAD cell-line mapping",x=NULL,y="Mapped LUAD models")+theme_paper+theme(legend.position="none")
save_plot(pS8A,"SupplementaryFigure8A_cellline_mapping",5.2,4.2)
pS8B <- ggplot(dataset_audit,aes(dataset,eligible_drugs_n10_AACgt0.2,fill=dataset))+geom_col()+
  geom_text(aes(label=eligible_drugs_n10_AACgt0.2),vjust=-.3)+labs(title="B  Drugs entering association tests",x=NULL,y="Eligible drugs")+
  theme_paper+theme(legend.position="none")
save_plot(pS8B,"SupplementaryFigure8B_drug_coverage",5.2,4.2)
classifier_scores <- rbindlist(list(
  fread(file.path(tables,"SupplementaryFigure8D_TCGA_classifier_scores.tsv"))[,cohort:="TCGA"],
  fread(file.path(tables,"SupplementaryFigure8E_PDMR_classifier_scores.tsv"))[,cohort:="PDMR"],
  fread(file.path(tables,"SupplementaryFigure8F_DepMap_classifier_scores.tsv"))[,cohort:="DepMap"]
),fill=TRUE)
classifier_threshold <- fromJSON(file.path(audit,"pdxe_lung_classifier_receipt.json"))$frozen_threshold
pS8C <- ggplot(classifier_scores[truth %in% c("LUAD","LUSC")],aes(state_score,fill=truth,color=truth))+
  geom_density(alpha=.22)+facet_wrap(~cohort,scales="free_y")+geom_vline(xintercept=classifier_threshold,linetype=2)+
  scale_fill_manual(values=c(LUAD="#B2182B",LUSC="#2166AC"))+scale_color_manual(values=c(LUAD="#B2182B",LUSC="#2166AC"))+
  labs(title="C  External validation of the frozen lung-state classifier",x="LUAD-minus-LUSC activity score",y="Density",fill=NULL,color=NULL)+theme_paper
save_plot(pS8C,"SupplementaryFigure8C_classifier_validation",8.0,4.3)

pS9A <- ggplot(head(coverage,20),aes(reorder(drug,PDX_n),PDX_n,fill=mapped_to_cellline_drug>0))+geom_col()+coord_flip()+
  scale_fill_manual(values=c(`TRUE`="#6A3D9A",`FALSE`="#BDBDBD"),name="Mapped to\ncell-line drug")+
  labs(title="A  PDXE v2 single-agent coverage in LUAD-like PDXs",x=NULL,y="PDX models")+theme_paper
save_plot(pS9A,"SupplementaryFigure9A_PDX_drug_coverage",6.5,5.5)
if (nrow(pdx_tests)) {
  pdx_tests[,plot_FDR:=pmax(FDR,1e-300)]
  pS9B <- ggplot(pdx_tests,aes(c_index-0.5,-log10(plot_FDR),color=cellline_replicated_same_direction))+
    geom_point(alpha=.7)+geom_hline(yintercept=-log10(.05),linetype=2)+geom_vline(xintercept=0,linetype=3)+
    scale_color_manual(values=c(`TRUE`="#B2182B",`FALSE`="#9E9E9E"),name="Cell-line replicated\nand same direction")+
    labs(title="B  Exploratory PDX drug–TF concordance screen",x="Concordance index − 0.5",y="−log10(FDR)")+theme_paper
} else {
  pS9B <- ggplot()+annotate("text",x=0,y=0,label="No PDX drug had at least 10 LUAD-like models",size=5)+xlim(-1,1)+ylim(-1,1)+theme_void()
}
save_plot(pS9B,"SupplementaryFigure9B_PDX_association_screen",6.5,5.2)

receipt <- list(
  status="passed", frozen_HC_TFs=length(hc), cellline_tests=nrow(cell_results),
  significant_cellline_pairs=nrow(sig), replicated_same_direction_pairs=nrow(replicated),
  LUAD_like_PDXs=length(luad_like), PDX_drugs_n_ge_10=uniqueN(pdx_tests$drug),
  PDX_tests=nrow(pdx_tests), PDX_validated_pairs=nrow(validated),
  selected_PDX_pair=if(nrow(selected)) as.list(selected[1]) else NULL,
  PDX_interpretation=if(nrow(selected)) "anchor-style validation criterion met" else "Figure 5G-H validation criterion not met",
  c_index_orientation=list(cell_lines="CI > 0.5: higher TF activity associates with higher AAC",
                           PDX="CI > 0.5: higher TF activity associates with lower best average response")
)
write_json(receipt,file.path(audit,"figure5_analysis_receipt.json"),pretty=TRUE,auto_unbox=TRUE,na="null")
cat(toJSON(receipt,auto_unbox=TRUE,na="null"),"\n")

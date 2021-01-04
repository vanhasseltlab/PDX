###Figures and metrics for paper

# libraries
library(tidyverse)
library(gridExtra)
library(cowplot)
library(egg)
library(ggnetwork)
library(network)
library(sna)
library(RxODE)
library(jpeg)
library(grid)
library(pdftools)

# import data
#From script 2
dirs <- sort(list.dirs(path = "results/NONMEM/", full.names = FALSE, recursive = FALSE), decreasing = TRUE)
runName <- grep("KG", dirs, value = TRUE)[1]
untreated <- read.table(paste0("results/NONMEM/", runName, "/m0_natural_growth.tab"), header = 1, skip = 1)
#used in Figure 2

#From script 5
load("data/clean/outcomes_allData.Rdata")
#used in Figure 2, 3, S1

#From script 7a 
load("results/MVLasso_predictions/KGKD - lambda_repeatedCV.Rdata")    
load("results/MVLasso_predictions/KGKD - predictions_repeatedCV.Rdata")
load("results/MVLasso_predictions/KGKD - sABC_over_time.Rdata")
#used in Figure 3

#From script 7b
load("results/beta_grplasso.Rdata") 
#used in Figure 4 and 5

#extract chosen treatments
treatments <- as.character(unique(lambda_mat$treatment))
chosen_treatments <- treatments[sapply(smooth_curves, function(x) any(x[-1, "cverror"] < x[1, "cverror"]))]
chosen_treatments <-  chosen_treatments[!chosen_treatments %in% c("TAS266", "LFW527 + binimetinib")]
#used in Figure 2, 3, 4, 5

#Figure 2: P1 – Nonmem – Observations and estimated curves #----------------------------------------------####
#prepare data
outcomes$ID <- sapply(strsplit(outcomes$ID_TREAT, "_"), function(x) x[1])
untreated <- untreated %>% filter(EVID == 0) %>% 
  mutate(ID = as.character(ID)) %>% 
  left_join(unique(outcomes[, c("ID", "ID_name")])) %>% 
  filter(!is.na(ID_name))
chosen_curves <- matrix(c("X-3483", "X-3483", "X-3483", "X-3483", 
                          "LDK378", "BKM120", "encorafenib + binimetinib", "encorafenib"), ncol = 2)
untreated_merge <- data.frame()
for (i in 1:nrow(chosen_curves)) {
  untreated_merge <- rbind.data.frame(untreated_merge, data.frame(untreated[untreated$ID_name == chosen_curves[i, 1], ], Treat_name = chosen_curves[i, 2]))
}
untreated_merge$natural_growth <- untreated_merge$IPRED

df_plot <- allData %>% 
  mutate(ID_Treat_key = paste(Treat_name, ID_name, sep = "_")) %>% 
  filter(ID_Treat_key %in% paste(chosen_curves[, 2], chosen_curves[, 1], sep = "_")) %>% 
  mutate(Treat_name = factor(Treat_name, levels = chosen_curves[, 2]))

summ_example <- df_plot %>% group_by(ID_Treat_key, ID_name, Treat_name) %>% 
  summarize(kg = kg[1], kd = kd[1], kr = kr[1]) %>% ungroup() %>% 
  pivot_longer(c(kg, kd, kr), values_to = "value", names_to = "Parameter") %>% 
  mutate(Parameter = factor(Parameter, levels = c("kg", "kd", "kr")))

colouring_blue <- colorRampPalette(c("#89A1F5", "black"))(70)[c(1, 1:6*5, (1:33) + 30)]
colouring_orange <- colorRampPalette(c("#F46B2D", "black"))(70)[c(1, 1:6*5, (1:33) + 30)]
colouring_black <- colorRampPalette(c("white", "black"))(70)[c(1:6*5, (1:34) + 30)+5]

dat_2a <- outcomes %>% filter(Treat_name %in% chosen_treatments) %>% 
  pivot_longer(cols = c(kg, kd, kr), names_to = "Parameter", values_to = "value") %>% 
  filter(!((value == 0 | value > 0.5) & Parameter == "kr")) %>% 
  mutate(Parameter = factor(Parameter, levels = c("kg", "kd", "kr")))
levels(dat_2a$Parameter) <- c(expression(italic("k"["g"])*" (day"^-1*")"),
                              expression(italic("k"["d"])*" (day"^-1*")"),
                              expression(italic("k"["r"])*" (day"^-1*")"))

mapping_param <- summ_example %>% 
  mutate(ID_Treat_key = ifelse(Parameter == "kg", "LDK378_X-3483", as.character(ID_Treat_key))) %>% 
  mutate(ID_Treat_key = ifelse(Parameter == "kr" & Treat_name != "encorafenib", "LDK378_X-3483", as.character(ID_Treat_key))) %>% 
  mutate(x = as.numeric(Parameter)*65 + 50, y = 1400)

mapping_param$colour <- colouring_blue[17]
mapping_param$colour[mapping_param$Parameter == "kd" & mapping_param$Treat_name == "LDK378"] <- colouring_black[1]
mapping_param$colour[mapping_param$Parameter == "kd" & mapping_param$Treat_name == "BKM120"] <- colouring_black[3]
mapping_param$colour[mapping_param$Parameter == "kd" & mapping_param$Treat_name == "encorafenib + binimetinib"] <- colouring_black[8]
mapping_param$colour[mapping_param$Parameter == "kd" & mapping_param$Treat_name == "encorafenib"] <- colouring_black[12]
mapping_param$colour[mapping_param$Parameter == "kr"] <-"#F47339"
mapping_param$colour[mapping_param$Parameter == "kr" & mapping_param$Treat_name == "encorafenib"] <- colouring_orange[3]

levels(mapping_param$Parameter) <- c(expression(italic("k"["g"])), expression(italic("k"["d"])),
                                     expression(italic("k"["r"])))

#Figure 2A 
figure2a <- dat_2a %>% 
  ggplot(aes(x = value)) +
  geom_histogram(aes(y = ..count../sum(..count..), fill = value), bins = 40, position = "dodge", 
                 fill = c(colouring_blue, colouring_black, colouring_orange)) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  facet_wrap(~ Parameter, strip.position = "bottom", scales = "free", labeller = label_parsed) +
  labs(y = "Proportion", x = NULL, tag = "A") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.placement = "outside")

#Figure 2B
figure2b <- ggplot(df_plot, aes(x = TIME)) + 
  geom_line(data = untreated_merge, aes(y = natural_growth, linetype = 'dashed'), colour = 'dark red') +
  geom_point(aes(y = DV, size = DOSE), shape = 1, colour = 'black') +
  geom_line(aes(y = IPRED, linetype = 'solid'), colour = 'dark red') +
  facet_wrap(~ Treat_name, nrow = 1) +
  labs(x = "Time (days)", y = bquote('Volume ('~mm^3~')'), tag = "B") +
  scale_linetype_manual(name = NULL, values = c('solid' = 'solid', 'dashed' = 'dashed'), 
                        labels = c('Modeled untreated','Modeled treated')) +
  scale_size_continuous(name = NULL, breaks = 1, labels = "Observed treated", range = c(1.5,1.5)) +
  scale_y_continuous(limits = c(0, 1500), expand = expansion(mult = c(0.01, 0))) +
  
  geom_label(data = mapping_param, aes(x = x, y = y-200, label = sprintf("%0.2f", round(value, digits = 2))), 
             fill = mapping_param$colour, colour = "white", show.legend = F) +
  geom_text(data = mapping_param, aes(x = x, y = y , label = Parameter), parse = T) +
  
  theme_bw() +
  theme(legend.spacing.x = unit(0.2, "cm"), text = element_text(size = 10), 
        legend.position="bottom", strip.background =  element_rect(fill =  "white"))


pdf(file = "paper/figures/Figure2_parameters.pdf", width = 8, height = 6)
grid.arrange(figure2a, figure2b, layout_matrix = matrix(c(2, 3, 2, 3, 2, 3), ncol = 3))
dev.off()



#Figure 3: P2 - Predictions - sABC #-----------------------------------------------------------------------------------------####

#prepare data
mod1 <- RxODE({
  d/dt(V) <- kg*V - kd*exp(-kr*t)*V
})
days <- seq(0, 100, by = 0.5)
ev <- eventTable(amount.units = 'mm^3', time.units = 'days')%>%
  add.sampling(days)
names(predictions) <- treatments

M <- 20
cv_pred <- predictions[["encorafenib"]]
ids <- paste0("X-", c("2861", "4644", "1441", "1906"))
sABC <- sABC_treat_total %>% filter(time == 56 & Treat_name == "encorafenib" & ID_name %in% ids)

data_abc <- data.frame()

#plot_avgcurve <- list()
for (id in ids) {
  if (ncol(cv_pred) == 2) {
    cv_pred <- abind::abind(cv_pred, array(0, replace(dim(cv_pred), 2, 1)), along = 2)
    colnames(cv_pred)[3] <- "kr"
  }
  pred_id <- cv_pred[id, ,]
  outcome_id <- outcomes[outcomes$ID_name == id & outcomes$Treat_name == "encorafenib", ]
  y_base <- outcome_id$base
  average_curve <- rep(0, length(days))
  
  out_truth <- mod1$solve(outcome_id[, c("kg", "kd", "kr")], ev, y_base)
  colnames(out_truth)[2] <- "truth"
  
  for (m in 1:M){
    input_pars_pred <- c(pred_id[1, m], pred_id[2, m], pred_id[3, m])
    out_pred <- mod1$solve(input_pars_pred, ev, y_base)
    colnames(out_pred) <- c("time", "predicted")
    average_curve <- average_curve + out_pred[, "predicted"]*(1/M)
  }
  
  average_curve <- data.frame(time = out_pred[, "time"], volume = average_curve)
  data_abc <- rbind.data.frame(data_abc, data.frame(merge(average_curve, out_truth), ID_name = id, limit_y = min(7000, max(c(out_truth[, "truth"], average_curve$volume)))))
}

observed <- allData %>% filter(Treat_name == "encorafenib" & ID_name %in% ids & TIME <= 56) %>% 
  left_join(sABC[, c("ID_name", "sABC_model")]) %>% 
  mutate(ID_and_sABC = factor(paste("sABC =", round(sABC_model, 3)), 
                              levels = unique(paste("sABC =", round(sABC_model, 3))[order(sABC_model)])))

#table with values for this plot
sABC_treats <- sABC_treat_total %>% filter(time == 56) %>%
  group_by(Treat_name) %>% summarize(mean_sABC = mean(sABC_model), sd_sABC = sd(sABC_model),
                                     median_sABC = median(sABC_model), min_sABC = min(sABC_model),
                                     max_sABC = max(sABC_model), percentile75 = quantile(sABC_model, 0.75))
order_treats <- as.character(sABC_treats$Treat_name[order(sABC_treats$median_sABC)])

sABC_treats <- sABC_treats %>% filter(Treat_name %in% chosen_treatments)

#Figure 3A
figure3a <- data_abc %>% left_join(sABC[, c("sABC_model","ID_name")]) %>% 
  mutate(ID_and_sABC = factor(paste("sABC =", round(sABC_model, 3)), 
                              levels = unique(paste("sABC =", round(sABC_model, 3))[order(sABC_model)]))) %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = volume, colour = "pred")) +
  geom_line(aes(y = truth, colour = "est")) +
  geom_ribbon(aes(ymin = 0, ymax = truth), fill = "orange", alpha = "0.3") +
  geom_ribbon(aes(ymin = truth, ymax = pmin(2500, volume)), fill = "black", alpha = "0.3") +
  labs(x = "Time (days)", y = bquote('Volume ('~mm^3~')'), tag = "A") +
  scale_y_continuous(limits = c(0, min(2500, max(c(out_truth[, "truth"], average_curve$volume)))),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(limits = c(0, 56),expand = expansion(mult = c(0, 0))) +
  facet_wrap(~ ID_and_sABC, nrow = 1) +
  scale_colour_manual(name = NULL, values = c('est' = '#f46e32', 'pred' = 'black'), 
                      labels = c('Estimation','Prediction from omics')) +
  theme_minimal() +
  theme(plot.title = element_text(size = 10), plot.margin = unit(c(5.5, 7.5, 5.5, 5.5), "pt"),
        legend.position = "top", panel.spacing = unit(1, "lines"), legend.margin = margin(0, 0, 0, 0),
        legend.box.margin=margin(-10, -10, -10, -10))

#Figure 3B
figure3b <- sABC_treat_total %>% filter(time == 56) %>% filter(Treat_name %in% chosen_treatments) %>% 
  mutate(Treat_name = fct_relevel(Treat_name, order_treats)) %>% 
  ggplot(aes(y = sABC_model, x = Treat_name)) +
  geom_boxplot(outlier.shape = NA, fill = "#FFE3B2", colour = "black") +
  geom_hline(yintercept = 0, colour = "#B2B2B2", linetype = 1) +
  labs(x = "Treatment", y = bquote(sABC["t = 56"]), tag = "B") +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.25), limits = c(0, max(sABC_treats$percentile75))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), legend.position = "top",
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "#B2B2B2"), 
        panel.grid.minor.y = element_line(color = "#B2B2B2"))

pdf(file = "paper/figures/Figure3_sABC.pdf", width = 8, height = 7)
grid.arrange(figure3a, figure3b, layout_matrix = matrix(c(1, 2, 2)))
dev.off()



#Figure 4: P3 - Pathway selection #-----------------------------------------------------------------------------------------####
#select treatments from cross-validation

PlotPathwaySelection <- function(beta_) {
  treat_by_pathway_df <- unique(beta_ %>% filter(treatment %in% chosen_treatments) %>% 
                                  mutate(outcome = as.character(outcome), treatment = tools::toTitleCase(as.character(treatment)) ) %>%  
                                  select(treatment, group, outcome))
  ind_both <- duplicated(treat_by_pathway_df[, c("treatment", "group")]) | rev(duplicated(rev(treat_by_pathway_df[, c("treatment", "group")])))
  treat_by_pathway_df$outcome[ind_both] <- "kd and kr"
  
  treat_by_pathway_df$outcome <- factor(treat_by_pathway_df$outcome, levels = c("kd", "kr", "kd and kr"))
  levels(treat_by_pathway_df$outcome) <- c("k[d]", "k[r]", "k[d] and k[r]")
  
  treat_by_pathway <- treat_by_pathway_df %>% 
    ggplot(aes(y = group, x = treatment, fill = outcome)) +
    geom_tile() +
    scale_fill_manual(values = c("black", "#F46B2D", "#89A1F5"), 
                      labels = c(parse(text = "k[d]"), parse(text = "k[r]"), parse(text = "k[d]~and~k[r]"))) +
    labs(x = "Treatment", y = "Pathway", fill = "Response") +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 80)) +
    scale_x_discrete(position = "bottom") +
    theme_bw() + 
    geom_vline(xintercept = seq(1.5, length(chosen_treatments) - 0.5, 1), colour = "grey85") +
    geom_hline(yintercept = seq(1.5, length(unique(treat_by_pathway_df$group)) - 0.5, 1), colour = "grey85") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust=0.2), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), legend.text.align = 0)
  
  
  treats_plot <- intersect(treat_by_pathway_df$treatment, tools::toTitleCase(chosen_treatments))
  
  histogram_treats <- treat_by_pathway_df %>% filter(treatment %in% treats_plot) %>% 
    ggplot(aes(x = treatment, fill = outcome)) +
    scale_fill_manual(values = c("black", "#F46B2D", "#89A1F5")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)), 
                       breaks = 0:max(table(treat_by_pathway_df$treatment) + 1)) +
    scale_x_discrete(breaks = NULL) +
    geom_bar(stat = "count") +
    geom_vline(xintercept = seq(1.5, length(treats_plot) - 0.5, 1), colour = "grey85") +
    labs(x = NULL, y = "Count") +
    theme_bw() +
    theme(axis.text.x = element_blank(), legend.position = "none", panel.grid.minor = element_blank())
  
  
  
  ggarrange(histogram_treats, treat_by_pathway, 
            nrow = 2, ncol = 1, widths = c(1), heights = c(1, 5))
}

pdf("paper/figures/Figure4_selected_pathways_Wiki.pdf", width = 9, height = 14, onefile = FALSE)
PlotPathwaySelection(beta_cn_wiki)
dev.off()


#Figure 5: P3 - Network pathways #-----------------------------------------------------------------------------------------####
#prepare data
beta_ <- beta_cn_wiki %>% 
  filter(treatment %in% chosen_treatments) %>% 
  mutate(pathway = str_sub(group,-6,-1), treatment = tools::toTitleCase(as.character(treatment)))
edges = "treatment"
nodes = "pathway"

test <- unique(beta_[, c(edges, nodes, "outcome")])

test_short <- unique(test) %>% group_by(treatment, pathway) %>% summarize(omic_type_n = n(), outcome = as.character(outcome[1]))
test_short$outcome[test_short$omic_type_n == 2] <- "kd and kr"
test_short <- as.data.frame(test_short)
test_short <- apply(test_short[, c("treatment", "pathway", "outcome")], 2, as.character)
treats <- unique(test_short[, nodes])
paths <- unique(test_short[, edges])

adj_mat <- matrix(0, ncol = length(treats), nrow = length(paths), dimnames = list(paths, treats))
edge_omic <- character(length = nrow(test_short))
for (i in 1:nrow(test_short)) {
  adj_mat[test_short[i, 1], test_short[i, 2]] <- 1
  edge_omic[i] <- test_short[i, 3]
}

network_ <- network(adj_mat, bipartite = T)
set.edge.attribute(network_, "outcome", edge_omic)
set.seed(2019)
ggnet_ <- ggnetwork(network_)
ggnet_$type <- ifelse(ggnet_$vertex.names %in% treats, nodes, edges)

ggnet_$outcome <- factor(ggnet_$outcome, levels = c("kd", "kr", "kd and kr"))
length(grep("WP", unique(ggnet_$vertex.names)))
colour_extra <- c(rep("white", 71), rep("black", 19))


#Figure 5
figure5 <- ggplot(ggnet_, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(color = outcome, linetype = outcome), curvature = 0) +
  geom_nodelabel(aes(label = vertex.names, fill = type), fontface = "bold", size = 2) +
  scale_fill_manual(values = c("#FFE3B2", "white")) +
  scale_color_manual(values = c("black", "#EFA98A", "#9BAEF3")) +
  lims(x = c(-0.05, 1.05)) +
  theme_blank() +
  theme(legend.title = element_blank())


pdf("paper/figures/Figure5_network_pathways.pdf", width = 8, height = 5)
print(figure5)
dev.off()

#Figure S1: P1 - Boxplot residuals NONMEM #--------------------------------------------------------------------------####
#prepare data
plots1_data <- allData
medians_order <- plots1_data %>% group_by(Treat_name) %>% 
  summarize(median_CWRES = median(CWRES))

plots1_data$Treat_name <- factor(plots1_data$Treat_name, levels = medians_order$Treat_name[order(medians_order$median_CWRES)])

#Figure S1
figures1 <- ggplot(plots1_data, aes(y = CWRES, x = Treat_name)) +
  geom_hline(yintercept = c(-2, 0, 2), colour = "grey35", linetype = c(2, 1, 2)) +
  geom_boxplot() +
  labs(x = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.2))


pdf(file = "paper/figures/FigureS1_nonmem_GOF.pdf", width = 8, height = 7)
print(figures1)
dev.off()


#Figure S2: P1 - Visual check NONMEM estimation #-----------------------------------------------------------####
PlotGrowthCurves <- function(allData, untreated, sorted = TRUE, treatments = NULL) {
  allData$ID_kd <- paste0(allData$ID_name, ", kg: ", allData$kg, "\nkd: ", allData$kd, "\nkr: ", allData$kr)
  
  predict_col <- "#F46E32"
  untreated_col <- "#5CB1EB"
  observation_col <- "#001158"
  
  if (is.null(treatments)) {
    treatments <- sort(unique(allData$Treat_name))
  }
  
  p_list <- list()
  for (treat in treatments) {
    dat1 <- allData[allData$Treat_name == treat, ]
    dat2 <- untreated[untreated$ID %in% dat1$ID, ]
    if (treat == "untreated") {
      dat1$ID_kd <- as.factor(dat1$ID_name)
      dat2 <- merge(dat2, unique(dat1[c("ID", "ID_kd")]), by = "ID")
      plot_col <- untreated_col
    } else {
      dat2 <- merge(dat2, unique(dat1[, c("ID", "ID_kd")]), by = "ID")
      plot_col <- predict_col
    }
    if (is.null(allData$converged)){
      title_plot <- paste0("Treatment ", treat)
    } else {
      title_plot <- paste0("Treatment ", treat, c(", NOT CONVERGED", ", CONVERGED")[dat1$converged[1] + 1])
    }
    p_list[[treat]] <- ggplot(dat1, aes(x = TIME, group = ID_kd)) +
      geom_point(aes(y = DV), colour = observation_col) +
      geom_line(data = dat2, aes(y = IPRED_NAT), colour = untreated_col, linetype = "dashed") +
      geom_line(aes(y = IPRED), colour = plot_col) +
      facet_wrap(~ ID_kd) +
      ggtitle(title_plot) +
      ylim(0, max(1200, max(dat1$DV))) +
      theme_bw() +
      theme(plot.title = element_text(size = 16, face = "bold"))
  }
  return(p_list)
}

untreated_plot <- untreated
names(untreated_plot)[11] <- "IPRED_NAT"
plotsChosen <- PlotGrowthCurves(allData, untreated_plot, sorted = TRUE)

pdf(file = "paper/figures/FigureS2_estimated_growth_curves.pdf", width = 20, height = 16)
for (treat in names(plotsChosen)) {
  print(plotsChosen[[treat]])
}
dev.off()

caption_S2 <- "Set of figures with the non-linear mixed effect model fits for all PDXs. Within every treatment, the volume (DV) is plotted against time (TIME) for every PDX. The modeled 
natural growth curves are shown (blue dashed line) and the model (red solid line) fit to the observed values (dark blue points). Tumor growth curves from treatment TAS266 
show a bad model fit."


pdf('paper/figures/FigureS2_caption.pdf',width = 20, height = 3)
old_mar <- par()$mar
par(mar= rep(0,4))
plot(NA, xlim=c(0,2), ylim=c(-0.5,0.8), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(0, 0.5, "Supplementary Figure S2:", font = 2, adj = 0, cex = 1.5)
text(0, 0,caption_S2, adj = 0, font = 1, cex = 1.5)
par(mar = old_mar)
dev.off()

pdf_combine(c("paper/figures/FigureS2_caption.pdf", "paper/figures/FigureS2_estimated_growth_curves.pdf"), 
            output = "paper/figures/FigureS2_estimated_growth_curves_with_caption.pdf")

#Figure S3: P3 - prediction errors with and without KR #----------------------------------------------------####

model_parameters <- c("KGKDKR", "KGKD")
sABC_models <- NULL
all_results_df <- NULL
long_outcomes <- outcomes %>% ungroup() %>% select(c(ID_name, Treat_name, kg, kd, kr)) %>% 
  pivot_longer(-c(ID_name, Treat_name), names_to = "parameter", values_to = "estimated")

for (m_p in model_parameters) {
  prefix_output <- paste0(m_p, " - ")
  #load data
  load(paste0("results/MVLasso_predictions/", prefix_output, "sABC_over_time.Rdata"))
  load(paste0("results/MVLasso_predictions/", prefix_output, "predictions_repeatedCV.Rdata"))# predictions smooth curves
  
  chosen_treatments <- unique(treatments[sapply(smooth_curves, function(x) any(x[-1, "cverror"] < x[1, "cverror"]))])
  chosen_treatments <-  chosen_treatments[!chosen_treatments %in% c("TAS266", "LFW527 + binimetinib")]
  
  sABC_models <- rbind(sABC_models, sABC_treat_total %>% filter(time == 56) %>% 
                         select(sABC_intercept, sABC_model, time, ID_name, Treat_name)%>%
                         mutate(included_parameters = m_p,
                                chosen = Treat_name %in% chosen_treatments))


  load(paste0("results/MVLasso_predictions/", prefix_output, "predictions_repeatedCV.Rdata"))
  
  if (is.null(names(predictions))) {
    names(predictions) <- treatments
  }
  dat_all_t <- NULL
  for (treat in treatments) {
    dat_t <- data.frame(kg = rowMeans(predictions[[treat]][, "kg", ]), Treat_name = treat, stringsAsFactors = F,
                        kd = rowMeans(predictions[[treat]][, "kd", ]),
                        kr = ifelse("kr" %in% colnames(predictions[[treat]]), rowMeans(predictions[[treat]][, "kr", ], na.rm = T), 0))
    dat_t$ID_name <- rownames(dat_t)
    dat_all_t <- rbind(dat_all_t, dat_t)
  }
  all_results_df <- rbind(all_results_df, dat_all_t %>% 
                            pivot_longer(-c(ID_name, Treat_name), names_to = "parameter", values_to = "predicted") %>% 
                            mutate(included_parameters = m_p) %>% right_join(long_outcomes))
}
all_results_df <- all_results_df %>% left_join(sABC_models)

#select treatments with non-zero kr
kr_treats <- NULL
load("data/clean/lasso_data_filtered.Rdata")
for (i in 1:length(treatments)) {
  dat_Y <- outcomes[outcomes$Treat_name == treatments[i], ]
  rownames(dat_Y) <- dat_Y$ID_name
  ids <- intersect(rownames(X_cn), dat_Y$ID_name)
  n <- length(ids)
  outcome <- c("kg", "kd", "kr")

    #Model with or without KR
  if (var(dat_Y[ids, "kr"]) > 0 & length(unique(unlist(dat_Y[ids, "kr"]))) > 2) {
    kr_treats <- c(kr_treats, treatments[i])
  }
}

#plot prediction error Kg and Kd
figures3a <- all_results_df %>%
  filter(Treat_name != "untreated" & ID_name != "X-4215") %>% 
  filter(Treat_name %in% kr_treats & parameter != "kr") %>% 
  mutate(Treat_name = fct_relevel(Treat_name, order_treats)) %>% 
  mutate(prediction_error = estimated - predicted) %>% 
  ggplot(aes(y = abs(prediction_error), x = Treat_name, fill = included_parameters)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, colour = "#B2B2B2", linetype = 1) +
  labs(y = "Absolute prediction error (parameter)", tag = "A") +
  scale_y_continuous(breaks = seq(-2, 2, by = 0.05), limits = c(0, 0.155)) +
  scale_fill_manual(values = c("grey55", "white")) +
  facet_wrap(~ parameter, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), legend.position = "top",
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "#B2B2B2"), 
        panel.grid.minor.y = element_line(color = "#B2B2B2"), axis.title.x = element_blank())

#plot prediction error sABC
figures3b <- all_results_df %>%
  filter(Treat_name != "untreated" & ID_name != "X-4215") %>% 
  filter(Treat_name %in% kr_treats & parameter == "kg") %>% 
  mutate(Treat_name = fct_relevel(Treat_name, order_treats)) %>% 
  ggplot(aes(y = sABC_model, x = Treat_name, fill = included_parameters)) +
  geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, colour = "#B2B2B2", linetype = 1) +
  labs(x = "Treatment", y = bquote(sABC["t = 56"]), tag = "B") +
  scale_y_continuous(breaks = seq(0, 2, by = 0.25), limits = c(0, 1.75)) +
  scale_fill_manual(values = c("grey55", "white"))  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2), legend.position = "none",
        panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "#B2B2B2"), 
        panel.grid.minor.y = element_line(color = "#B2B2B2"))

pdf(file = "paper/figures/FigureS3_prediction_errors.pdf", width = 8, height = 11)
grid.arrange(figures3a, figures3b, layout_matrix = matrix(c(1, 1, 1, 2, 2)))
dev.off()


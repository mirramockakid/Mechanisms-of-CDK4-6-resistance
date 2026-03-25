design <- model.matrix(~ 0 + Group + Pool, data = meta)
v <- limma::voom(dge, design)

fit <- limma::lmFit(v, design)

contrast.matrix <- limma::makeContrasts(
  "MCF7 Abemaciclib 0.5" = GroupX0.5.uM.Abemaciclib.resistant.MCF7 - GroupCDK4.6i.sensitive.MCF7,
  "MCF7 Abemaciclib 1" = GroupX1.0.uM.Abemaciclib.resistant.MCF7 - GroupCDK4.6i.sensitive.MCF7,
  "MCF7 Abemaciclib 1.5" = GroupX1.5.uM.Abemaciclib.resistant.MCF7 - GroupCDK4.6i.sensitive.MCF7,
  "MCF7 Palbociclib 1.2" = GroupX1.2.uM.Palbociclib.resistant.MCF7 - GroupCDK4.6i.sensitive.MCF7,
  "MCF7 Palbociclib 2.4" = GroupX2.4.uM.Palbociclib.resistant.MCF7 - GroupCDK4.6i.sensitive.MCF7,
  "MCF7 Palbociclib 3.6" = GroupX3.6.uM.Palbociclib.resistant.MCF7 - GroupCDK4.6i.sensitive.MCF7,
  "MCF7 Palbociclib 4.8" = GroupX4.8.uM.Palbociclib.resistant.MCF7 - GroupCDK4.6i.sensitive.MCF7,
  "T47D Palbociclib 1.2" = GroupX1.2.uM.Palbociclib.resistant.T47D - GroupCDK4.6i.sensitive.T47D,
  "T47D Palbociclib 2.4" = GroupX2.4.uM.Palbociclib.resistant.T47D - GroupCDK4.6i.sensitive.T47D,
  "T47D Palbociclib 3.6" = GroupX3.6.uM.Palbociclib.resistant.T47D - GroupCDK4.6i.sensitive.T47D,
  "T47D Palbociclib 4.8" = GroupX4.8.uM.Palbociclib.resistant.T47D - GroupCDK4.6i.sensitive.T47D,
  "MCF7 vs T47D" = GroupCDK4.6i.sensitive.MCF7 - GroupCDK4.6i.sensitive.T47D,

  levels = design
)

contrast <- limma::contrasts.fit(fit, contrast.matrix)
contrast <- limma::eBayes(contrast)

results_list <- lapply(colnames(contrast), function(contrast_name) {
  res <- limma::topTable(
    contrast,
    coef = contrast_name,
    number = Inf,
    sort.by = "P"
  ) |> dplyr::select(Gene = gene_id, logFC, P.Value, adj.P.Val)
}) |> setNames(colnames(contrast))


lapply(names(results_list), function(contrast_name) {

  res <- results_list[[contrast_name]]
  res <- res[order(res$adj.P.Val), ]

  # results table
  write.csv(res, 
  file = paste0("tables/", contrast_name, "_dea_results.csv"),
  row.names = FALSE
  )

  plot_volcano(results_list[[contrast_name]], file_name = paste0(contrast_name, "_volcano.pdf"), title = paste0(contrast_name))

})

lapply(names(results_list), function(contrast_name) {
  plot_heatmap(expr_mat = v$E, metadata = meta, results = results_list[[contrast_name]], gene_names = v$genes$gene_id, file_name = paste0(contrast_name, "_heatmap.pdf"), annotation = c("Resistance", "Dose uM", "Cell line"))
})


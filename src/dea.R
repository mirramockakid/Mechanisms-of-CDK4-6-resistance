design <- model.matrix(~ 0 + Group, data = meta)
v <- limma::voom(dge, design)

fit <- limma::lmFit(v, design)

contrast.matrix <- makeContrasts(
  "MCF7 vs T47D" = GroupCDK4.6i.sensitive.MCF7 - GroupCDK4.6i.sensitive.T47D,
  levels = design
)

contrast <- limma::contrasts.fit(fit, contrast.matrix)
contrast <- limma::eBayes(contrast)
res <- topTable(contrast, coef = "MCF7 vs T47D", number = Inf)

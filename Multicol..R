numeric_vars <- c("Temperature", "Rain", "Sun", "Wind", "Raining on sampling day")
corr_numeric <- corr[, numeric_vars]
cor_matrix <- cor(corr_numeric, method = "spearman", use = "pairwise.complete.obs")
corrplot(cor_matrix, method = "color", type = "upper",
         order = "hclust", tl.col = "black", tl.srt = 45,
         addCoef.col = "black", number.cex = 0.8)
tiff('Spearman rank corr.tiff', units="in", width=7, height=5, res=300)
plot(corrplot(cor_matrix, method = "color", type = "upper",
              order = "hclust", tl.col = "black", tl.srt = 45,
              addCoef.col = "black", number.cex = 0.8))
dev.off()
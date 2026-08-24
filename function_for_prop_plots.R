save_prop_plots <- function(x, output_dir, gender_label = "none", facet = "none", ...) {
  library(ggplot2)
  library(reshape2)
  library(plyr)
  
  cd <- data.frame(mice::complete(x, "long", include = TRUE))
  cd$.imp <- factor(cd$.imp)
  
  r <- as.data.frame(is.na(x$data))
  impcat <- x$meth != "" & sapply(x$data, is.factor)
  vnames <- names(impcat)[impcat]
  
  # Add dep_anx if it exists and isn't already included
  if (!"dep_anx" %in% vnames && "dep_anx" %in% names(x$data)) {
    vnames <- c(vnames, "dep_anx")
  }
  
  for (xvar in vnames) {
    cd_tmp <- cd
    select <- cd_tmp$.imp != 0 & !r[[xvar]]
    cd_tmp[select, xvar] <- NA
    
    # Convert dep_anx to factor for plotting if needed
    if (xvar == "dep_anx") {
      cd_tmp[[xvar]] <- as.factor(cd_tmp[[xvar]])
    }
    
    meltDF <- melt(cd_tmp[, c(".imp", xvar)], id.vars = ".imp", variable.name = "variable")
    meltDF <- meltDF[!is.na(meltDF$value), ]
    
    a <- ddply(meltDF, c(".imp", "variable", "value"), summarize, count = length(value))
    b <- ddply(meltDF, c(".imp", "variable"), summarize, tot = length(value))
    mdf <- merge(a, b)
    mdf$prop <- mdf$count / mdf$tot
    
    plotDF <- merge(unique(meltDF), mdf)
    plotDF$value <- factor(plotDF$value,
                           levels = levels(x$data[[xvar]]),
                           ordered = TRUE)
    
    filename <- if (tolower(gender_label) == "none") {
      file.path(output_dir, "multiple_imputation", "mi_diagnostics", "prop_plot_3int_072025", paste0("prop_", xvar, ".png"))
    } else {
      file.path(output_dir, "multiple_imputation", "mi_diagnostics", "prop_plot_3int_072025", paste0("prop_", xvar, "_", gender_label, ".png"))
    }
    
    p <- ggplot(plotDF, aes(x = value, fill = .imp, y = prop)) +
      geom_bar(position = "dodge", stat = "identity") +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.position = "none") +
      ylab("Proportion") +
      xlab(xvar) +
      scale_fill_manual(
        name = "",
        values = c("black", colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(x$m + 3)[1:(x$m + 3)])
      ) +
      guides(fill = guide_legend(nrow = 1)) +
      ggtitle(paste("Observed vs Imputed Proportions -", xvar))
    
    ggsave(
      filename = filename,
      plot = p,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
}

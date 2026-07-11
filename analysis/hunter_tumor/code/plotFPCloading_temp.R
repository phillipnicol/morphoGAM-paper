#Temporary version of plotFPCloading

plotFPCloading_temp <- 
function(mgam_object,
         curve.fit,
         L = 1,
         num_genes = 5,
         type = "t",
         sign = "+") {

  if(type == "t") {
    if(sign == "-") {
      top5 <- order(mgam_object$fpca.t$u[, L], decreasing = FALSE)[1:num_genes]
    } else {
      top5 <- order(mgam_object$fpca.t$u[, L], decreasing = TRUE)[1:num_genes]
    }

    if(sum(mgam_object$fxs.t[top5, ] %*% mgam_object$fpca.t$v[, L]) < 0) {
      mgam_object$fpca.t$v[, L] <- -1 * mgam_object$fpca.t$v[, L]
    }

    mat <- matrix(0, nrow = ncol(mgam_object$fxs.t), ncol = num_genes + 1)
    mat[, 1] <- mgam_object$fpca.t$d[L] * mgam_object$fpca.t$v[, L]
    colnames(mat) <- rep("Eigenfn", num_genes + 1)

    for(i in 1:num_genes) {
      mat[, i + 1] <- mgam_object$fxs.t[top5[i], ]
      colnames(mat)[i + 1] <- rownames(mgam_object$fxs.t)[top5[i]]
    }

    mat <- mat[, -1]
    df <- as.data.frame(mat)
    df$t <- curve.fit$xyt$t
    df <- df |> pivot_longer(cols = -t)

    base_colors <- brewer.pal(9, "Set1")[-6]

    if (num_genes <= length(base_colors)) {
      custom_colors <- base_colors[1:num_genes]
    } else {
      custom_colors <- colorRampPalette(base_colors)(num_genes)
    }

    names(custom_colors) <- colnames(mat)

    label_data <- df %>%
      group_by(name) %>%
      slice_max(value, with_ties = FALSE, n = 1)

    p <- ggplot(data = df, aes(x = t, y = value,
                               group = name,
                               color = name)) +
      geom_line(aes(size = ifelse(name == "Eigenfn", 2, 0.75)),
                show.legend = FALSE) +
      scale_color_manual(values = custom_colors) +
      scale_size_identity() +
      theme_bw() +
      geom_text_repel(data = label_data,
                      aes(label = name),
                      size = 4,
                      nudge_x = 0.1,
                      hjust = 0,
                      show.legend = FALSE,
                      max.overlaps = Inf) +
      ylab("Log FC from baseline") +
      geom_abline(slope = 0, intercept = 0, color = "grey", linetype = "dashed") +
      xlab("t")

    return(p)

  } else {

    if(sign == "-") {
      top5 <- order(mgam_object$fpca.r$u[, L], decreasing = FALSE)[1:num_genes]
    } else {
      top5 <- order(mgam_object$fpca.r$u[, L], decreasing = TRUE)[1:num_genes]
    }

    if(sum(mgam_object$fxs.r[top5, ] %*% mgam_object$fpca.r$v[, L]) < 0) {
      mgam_object$fpca.r$v[, L] <- -1 * mgam_object$fpca.r$v[, L]
    }

    mat <- matrix(0, nrow = ncol(mgam_object$fxs.r), ncol = num_genes + 1)
    mat[, 1] <- mgam_object$fpca.r$d[L] * mgam_object$fpca.r$v[, L]
    colnames(mat) <- rep("Eigenfn", num_genes + 1)

    for(i in 1:num_genes) {
      mat[, i + 1] <- mgam_object$fxs.r[top5[i], ]
      colnames(mat)[i + 1] <- rownames(mgam_object$fxs.r)[top5[i]]
    }

    mat <- mat[, -1]
    df <- as.data.frame(mat)
    df$r <- curve.fit$xyt$r
    df <- df |> pivot_longer(cols = -r)

    base_colors <- brewer.pal(9, "Set1")[-6]

    if (num_genes <= length(base_colors)) {
      custom_colors <- base_colors[1:num_genes]
    } else {
      custom_colors <- colorRampPalette(base_colors)(num_genes)
    }

    names(custom_colors) <- colnames(mat)

    label_data <- df %>%
      group_by(name) %>%
      slice_max(value, with_ties = FALSE, n = 1)

    p <- ggplot(data = df, aes(x = r, y = value,
                               group = name,
                               color = name)) +
      geom_line(aes(size = ifelse(name == "Eigenfn", 2, 0.75)),
                show.legend = FALSE) +
      scale_color_manual(values = custom_colors) +
      scale_size_identity() +
      theme_bw() +
      geom_text_repel(data = label_data,
                      aes(label = name),
                      size = 4,
                      nudge_x = 0.1,
                      hjust = 0,
                      show.legend = FALSE,
                      max.overlaps = Inf) +
      ylab("Log FC from baseline") +
      xlab("r")

    return(p)
  }
}
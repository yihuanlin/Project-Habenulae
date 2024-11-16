#!/usr/bin/env RScript
rm(list = ls(all = TRUE)) # Remove previous objects from globe environment
tryCatch(assign("last.warning", NULL, envir = baseenv()), error = function(e) {}) # Remove previous warnings
tryCatch(dev.off(dev.list()["RStudioGD"]), error = function(e) {}) # Remove previous plots in RStudio

loadPackage <- function(x, ...) {
  x <- c(x, ...)
  un <- x[!(x %in% installed.packages()[, "Package"])]
  if (length(un)) {
    invisible(readline(prompt = paste("Press [Enter] to install package \"", paste(un, collapse = ", "), "\"", sep = "")))
    install.packages(un, dependencies = TRUE)
  }
  suppressPackageStartupMessages(sapply(x, require, character.only = TRUE))
}

posthocTGH <- function(y, x, method = c("games-howell", "tukey"),
                       conf.level = 0.95, digits = 2,
                       p.adjust = "none", formatPvalue = TRUE) {
  ### Based on http://www.psych.yorku.ca/cribbie/6130/games_howell.R and userfriendlyscience package (version 0.7.2)
  method <- tolower(method)
  tryCatch(method <- match.arg(method), error = function(err) {
    stop("Argument for 'method' not valid!")
  })

  res <- list(input = as.list(environment()))

  res$intermediate <- list(
    x = factor(x[complete.cases(x, y)]),
    y = y[complete.cases(x, y)]
  )
  res$intermediate$n <- tapply(y, x, length)
  res$intermediate$groups <- length(res$intermediate$n)
  res$intermediate$df <- sum(res$intermediate$n) - res$intermediate$groups
  res$intermediate$means <- tapply(y, x, mean)
  res$intermediate$variances <- tapply(y, x, var)
  res$intermediate$names <- levels(res$intermediate$x)
  res$intermediate$pairNames <- combn(res$intermediate$groups, 2, function(ij) {
    paste0(rev(res$intermediate$names[ij]), collapse = "-")
  })

  res$intermediate$descriptives <- cbind(
    res$intermediate$n,
    res$intermediate$means,
    res$intermediate$variances
  )
  rownames(res$intermediate$descriptives) <- levels(res$intermediate$x)
  colnames(res$intermediate$descriptives) <- c("n", "means", "variances")

  ### Start on Tukey
  res$intermediate$errorVariance <-
    sum((res$intermediate$n - 1) * res$intermediate$variances) /
      res$intermediate$df
  res$intermediate$se <- combn(res$intermediate$groups, 2, function(ij) {
    sqrt(res$intermediate$errorVariance * sum(1 / res$intermediate$n[ij]))
  })
  res$intermediate$dmeans <- combn(res$intermediate$groups, 2, function(ij) {
    diff(res$intermediate$means[ij])
  })
  res$intermediate$t <- abs(res$intermediate$dmeans) / res$intermediate$se
  res$intermediate$p.tukey <- ptukey(res$intermediate$t * sqrt(2),
    res$intermediate$groups,
    res$intermediate$df,
    lower.tail = FALSE
  )
  res$intermediate$alpha <- (1 - conf.level)
  res$intermediate$qcrit <- qtukey(res$intermediate$alpha,
    res$intermediate$groups,
    res$intermediate$df,
    lower.tail = FALSE
  ) / sqrt(2)
  res$intermediate$tukey.low <- res$intermediate$dmeans - (res$intermediate$qcrit * res$intermediate$se)
  res$intermediate$tukey.high <- res$intermediate$dmeans + (res$intermediate$qcrit * res$intermediate$se)
  res$output <- list()
  res$output$tukey <- data.frame(
    res$intermediate$dmeans,
    res$intermediate$tukey.low,
    res$intermediate$tukey.high,
    res$intermediate$t,
    res$intermediate$df,
    res$intermediate$p.tukey
  )
  columnNames <- c("diff", "ci.lo", "ci.hi", "t", "df", "p")
  if (p.adjust != "none") {
    res$output$tukey$p.tukey.adjusted <- p.adjust(res$intermediate$p.tukey,
      method = p.adjust
    )
    columnNames <- c(columnNames, "p.adjusted")
  }

  rownames(res$output$tukey) <- res$intermediate$pairNames
  colnames(res$output$tukey) <- columnNames

  ### Start on Games-Howell
  res$intermediate$df.corrected <- combn(res$intermediate$groups, 2, function(ij) {
    sum(res$intermediate$variances[ij] /
      res$intermediate$n[ij])^2 /
      sum((res$intermediate$variances[ij] /
        res$intermediate$n[ij])^2 /
        (res$intermediate$n[ij] - 1))
  })
  res$intermediate$se.corrected <- combn(res$intermediate$groups, 2, function(ij) {
    sqrt(sum(res$intermediate$variances[ij] / res$intermediate$n[ij]))
  })
  res$intermediate$t.corrected <- abs(res$intermediate$dmeans) / res$intermediate$se.corrected

  res$intermediate$qcrit.corrected <-
    qtukey(res$intermediate$alpha,
      res$intermediate$groups,
      res$intermediate$df.corrected,
      lower.tail = FALSE
    ) / sqrt(2)

  res$intermediate$gh.low <- res$intermediate$dmeans -
    res$intermediate$qcrit.corrected * res$intermediate$se.corrected
  res$intermediate$gh.high <- res$intermediate$dmeans +
    res$intermediate$qcrit.corrected * res$intermediate$se.corrected


  res$intermediate$p.gameshowell <- ptukey(res$intermediate$t.corrected * sqrt(2),
    res$intermediate$groups,
    res$intermediate$df.corrected,
    lower.tail = FALSE
  )
  res$output$games.howell <- data.frame(
    res$intermediate$dmeans,
    res$intermediate$gh.low,
    res$intermediate$gh.high,
    res$intermediate$t.corrected,
    res$intermediate$df.corrected,
    res$intermediate$p.gameshowell
  )
  columnNames <- c("diff", "ci.lo", "ci.hi", "t", "df", "p")
  if (p.adjust != "none") {
    res$output$games.howell$p.gameshowell.adjusted <- p.adjust(res$intermediate$p.gameshowell,
      method = p.adjust
    )
    columnNames <- c(columnNames, "p.adjusted")
  }
  rownames(res$output$games.howell) <- res$intermediate$pairNames
  colnames(res$output$games.howell) <- columnNames

  ### Set class and return object
  class(res) <- "posthocTGH"
  return(res)
}

print.posthocTGH <- function(x, digits = x$input$digits, ...) {
  print(x$intermediate$descriptives, digits = digits)
  cat("\n")
  if (x$input$method == "tukey") {
    dat <- x$output$tukey
  } else if (x$input$method == "games-howell") {
    dat <- x$output$games.howell
  }
  dat[, 6:ncol(dat)] <- sapply(dat[, 6:ncol(dat)],
    format.pval,
    digits = digits,
    includeP = FALSE
  )
  print(dat, digits = digits)
}

dCohens <- function(x, y = FALSE, mu = FALSE) {
  x <- na.omit(x)
  if (is.numeric(mu)) {
    cd <- abs(mean(x) - mu) / sd(x)
    return(signif(cd, digits = 4))
  } else {
    y <- na.omit(y)
    lx <- length(x) - 1
    ly <- length(y) - 1
    md <- abs(mean(x) - mean(y))
    csd <- lx * var(x) + ly * var(y)
    csd <- csd / (lx + ly)
    csd <- sqrt(csd)
    cd <- md / csd
    return(signif(cd, digits = 4))
  }
}

fCohens <- function(anova) {
  su <- summary(anova)[[1]]
  l <- length(su[, "Sum Sq"])
  if (l > 2) {
    for (i in 1:(l - 1)) {
      eta <- su[i, "Sum Sq"] / (su[i, "Sum Sq"] + su[l, "Sum Sq"])
      f <- sqrt(eta^2 / (1 - eta^2))
      cat(paste("Cohen's F for", rownames(su)[i], "is", signif(f, digits = 4), "(Partial)", "\n"))
    }
  }
  eta <- sum(su[1:l - 1, "Sum Sq"]) / sum(su[, "Sum Sq"])
  f <- sqrt(eta^2 / (1 - eta^2))
  cat(paste("One way Cohen's F:", signif(f, digits = 4), "\n"))
}

wCohens <- function(df) {
  chi <- suppressWarnings(chisq.test(df))
  w <- sqrt(chi$statistic / sum(chi$observed))
  cat(paste("Cohen's W:", signif(w, digits = 4), "\n"))
}

mean.getConfidenceInterval <- function(x, SD, n) {

}

mean.getDistribution <- function(x, dist = "CI") {
  x <- x[complete.cases(x)]
  mean <- mean(x)
  N <- length(x)
  SD <- sd(x)
  SE <- SD / sqrt(N)
  ME <- qt(0.975, N - 1) * SE # Get 95% confidence interval by default
  if (dist == "CI") {
    dist <- ME
  } else if (dist == "SE") {
    dist <- SE
  } else if (dist == "SD") {
    dist <- SD
  }
  return(paste(
    sprintf("%.4f", mean), "±",
    sprintf("%.4f", dist), "; N = ",
    sprintf("%.1f", N), "; from ",
    sprintf("%.4f", mean - ME), " to ",
    sprintf("%.4f", mean + ME),
    sep = ""
  ))
}

plotPieChart <- function(df) {
  title <- sub("\\..*", "", deparse(substitute(df)))
  labels <- c("Normal", "Mild", "Severe", "Dead")
  labels <- labels[!df$value == 0]
  colours <- c("#F9F5EB", "#EAE3D2", "#81C6E8", "#607EAA")
  colours <- colours[!df$value == 0]
  df <- df[!df$value == "0", ]
  df <- df %>%
    mutate(
      csum = rev(cumsum(rev(value))),
      pos = value / 2 + lead(csum, 1),
      pos = if_else(is.na(pos), value / 2, pos),
      percentage = value / sum(value),
      label = paste0(round(value, 1), " (", round(percentage * 100), "%", ")")
    )
  if (title == "crtl") {
    title.a <- bquote("Uninjected" ~ "(" * italic("N") ~ "=" ~ .(paste0(sum(df$value), ")")))
  } else {
    title.a <- bquote(italic(.(title)) ~ "(" * italic("N") ~ "=" ~ .(paste0(sum(df$value), ")")))
  }
  plot <- ggplot(df, aes(x = "", y = value, fill = Severity)) +
    geom_col(color = "#000000") +
    coord_polar(theta = "y") +
    geom_label_repel(data = df, aes(y = pos, label = label), size = 4.5, nudge_x = 0, show.legend = FALSE, min.segment.length = 10) +
    scale_fill_manual(labels = labels, values = colours) +
    ggtitle(title.a) +
    theme_void() +
    theme(plot.title.position = "plot", plot.title = element_text(hjust = 0.34))
  ggsave(plot = plot, filename = paste0(title, ".png"), width = 4, height = 4)
  # plot
}

getLateralityDF <- function(df) {
  name <- sub("\\..*", "", deparse(substitute(df)))
  df <- df[, -c(2:5)]
  df <- df[rowSums(df[, sapply(df, is.numeric)]) != 0, ]
  df <- data.frame(Treatment = name, df)
  return(df)
}

getLateralityCT <- function(df) {
  name <- sub("\\..*", "", deparse(substitute(df)))
  df <- data.frame(
    WT.Laterality = sum(df$WT.Laterality),
    dHb.Vicera.Inverted = sum(df$dHb.Vicera.Inverted),
    Vicera.Inverted = sum(df$Vicera.Inverted),
    dHb.Inverted = sum(df$dHb.Inverted)
  )
  rownames(df) <- name
  return(df)
}

addSignificanceMarks <- function(df, y = FALSE, mu = FALSE) {
  for (col in colnames(df)) {
    if (is.logical(y)) {
      p <- t.test(df[, col], mu = 0)$p.value
    } else {
      p <- t.test(df[, col], y[, col])$p.value
    }
    if (p > 0.05) {
      labels <- "" # ns
    } else if (p < 0.05 && p >= 0.01) {
      labels <- "*"
    } else if (p < 0.01 && p >= 0.001) {
      labels <- "**"
    } else if (p < 0.001) {
      labels <- "***"
    }
    if (is.logical(mu)) {
      p <- 0
    } else {
      p <- t.test(df[, col], mu = mu)$p.value
    }
    if (p > 0.05) {
      text(x = which(colnames(df) == col), y = max(df) * 1.002, labels = labels, pos = 3, xpd = TRUE, col = "red")
    } else {
      text(x = which(colnames(df) == col), y = max(df) * 1.002, labels = labels, pos = 3, xpd = TRUE)
    }
  }
}

normaliseColumns <- function(df) {
  for (col in names(df)) {
    df[[col]] <- df[[col]] / sum(df[[col]]) * length(df[[col]])
  }
  return(df)
}

df.mean.getDistribution <- function(df) {
  for (col in names(df)) {
    print(col)
    print(mean.getDistribution(df[[col]]))
    cat("\n")
  }
}

loadPackage("ggplot2", "ggrepel", "dplyr", "readxl", "ggsignif")
if (requireNamespace("rstudioapi", quietly = TRUE)) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}


toxicity <- read.csv("./toxicity and laterality.csv")
toxicity[is.na(toxicity)] <- 0
toxicity <- toxicity[, !colnames(toxicity) == "X"]
crtl.tox <- subset(toxicity, (Treatment == "Uninjected") | (Treatment == "Uninjected.Practice") | (Treatment == "Uninjected.Aaron"))[, -1]
crtl.tox.sum <- c(sum(crtl.tox$Normal), sum(crtl.tox$Mild), sum(crtl.tox$Severe), sum(crtl.tox$Dead))
crtl.tox.sum <- data.frame(value = crtl.tox.sum, Severity = c("1", "2", "3", "4"))
dbx1b.tox <- subset(toxicity, Treatment == "dbx1b")[, -1]
dbx1b.tox.sum <- c(sum(dbx1b.tox$Normal), sum(dbx1b.tox$Mild), sum(dbx1b.tox$Severe), sum(dbx1b.tox$Dead))
dbx1b.tox.sum <- data.frame(value = dbx1b.tox.sum, Severity = c("1", "2", "3", "4"))
irx7.tox <- subset(toxicity, Treatment == "irx7")[, -1]
irx7.tox.sum <- c(sum(irx7.tox$Normal), sum(irx7.tox$Mild), sum(irx7.tox$Severe), sum(irx7.tox$Dead))
irx7.tox.sum <- data.frame(value = irx7.tox.sum, Severity = c("1", "2", "3", "4"))
lrp5.tox <- subset(toxicity, Treatment == "lrp5")[, -1]
lrp5.tox.sum <- c(sum(lrp5.tox$Normal), sum(lrp5.tox$Mild), sum(lrp5.tox$Severe), sum(lrp5.tox$Dead))
lrp5.tox.sum <- data.frame(value = lrp5.tox.sum, Severity = c("1", "2", "3", "4"))
GFP.tox <- subset(toxicity, Treatment == "GFP")[, -1]
GFP.tox.sum <- c(sum(GFP.tox$Normal), sum(GFP.tox$Mild), sum(GFP.tox$Severe), sum(GFP.tox$Dead))
GFP.tox.sum <- data.frame(value = GFP.tox.sum, Severity = c("1", "2", "3", "4"))
plotPieChart(irx7.tox.sum)
plotPieChart(lrp5.tox.sum)
plotPieChart(GFP.tox.sum)
plotPieChart(crtl.tox.sum)
plotPieChart(dbx1b.tox.sum)
crtl.laterality <- getLateralityDF(crtl.tox)
dbx1b.laterality <- getLateralityDF(dbx1b.tox)
irx7.laterality <- getLateralityDF(irx7.tox)
GFP.laterality <- getLateralityDF(GFP.tox)
laterality <- rbind(getLateralityCT(crtl.laterality), getLateralityCT(dbx1b.laterality), getLateralityCT(irx7.laterality), getLateralityCT(GFP.laterality))
laterality.dHb <- data.frame(WT = laterality$WT.Laterality + laterality$Vicera.Inverted, Inverted = laterality$dHb.Vicera.Inverted + laterality$dHb.Inverted)
rownames(laterality.dHb) <- rownames(laterality)
laterality.vicera <- data.frame(WT = laterality$WT.Laterality + laterality$dHb.Inverted, Inverted = laterality$dHb.Vicera.Inverted + laterality$Vicera.Inverted)
rownames(laterality.vicera) <- rownames(laterality)
laterality.dHb[3, 2] <- laterality.dHb[3, 2] + 2 # Include 2 samples where vicera organs not stained but dHb inverted
# laterality.dHb[3, 1] <- 32 # Use gata2a:GFP screen instead of ISH
# laterality.dHb[3, 2] <- 16 # Use gata2a:GFP screen instead of ISH
# fisher.test(laterality[c(1, 3), ]) # irx7 vs crtl
# fisher.test(laterality[c(3, 4), ]) # irx7 vs GFP
fisher.test(laterality.vicera[c(3, 4), ]) # irx7 vs GFP
wCohens(laterality.vicera[c(3, 4), ])
fisher.test(laterality.dHb[c(3, 4), ]) # irx7 vs GFP
wCohens(laterality.dHb[c(3, 4), ])
fisher.test(laterality.vicera[c(2, 4), ]) # dbx1b vs GFP
wCohens(laterality.vicera[c(2, 4), ])
fisher.test(laterality.dHb[c(2, 4), ]) # dbx1b vs GFP
wCohens(laterality.dHb[c(2, 4), ])
suppressWarnings(chisq.test(laterality[c(3, 4), ])$expected)
laterality.prop <- round(laterality / rowSums(laterality), 3) * 100


suppressWarnings(MiSeq <- read_excel("MiSeq.xlsx", col_types = c("text", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "skip", "skip")))
g1.modified <- MiSeq$`%Modified`[grepl("^F", MiSeq$Well)]
g2.modified <- MiSeq$`%Modified`[grepl("^G", MiSeq$Well)]
g3.modified <- MiSeq$`%Modified`[grepl("^H", MiSeq$Well)]
g1.frameshift <- MiSeq$`%Frameshift`[grepl("^F", MiSeq$Well)]
g2.frameshift <- MiSeq$`%Frameshift`[grepl("^G", MiSeq$Well)]
g3.frameshift <- MiSeq$`%Frameshift`[grepl("^H", MiSeq$Well)]
total.frameshift <- MiSeq$`%Estimated Total Frameshift`
total.frameshift <- total.frameshift[!is.na(total.frameshift)]
dbx1b.modified <- data.frame(g1 = g1.modified[1:6], g2 = g2.modified[1:6], g3 = g3.modified[1:6])
control.modified <- data.frame(g1 = g1.modified[7:9], g2 = g2.modified[7:9], g3 = g3.modified[7:9])
dbx1b.frameshift <- data.frame(g1 = g1.frameshift[1:6], g2 = g2.frameshift[1:6], g3 = g3.frameshift[1:6], bi = total.frameshift[1:6]^2)
control.frameshift <- data.frame(g1 = g1.frameshift[7:9], g2 = g2.frameshift[7:9], g3 = g3.frameshift[7:9], bi = total.frameshift[7:9]^2)
colour <- c("#00a3ff", "#5dd934", "#f9bb00", "#ff2501", "#d41876", "#929292", "#b0643c", "#2994ac", "#3f7cad")
weights <- sapply(seq(1, nrow(MiSeq), by = 3), function(i) sum(MiSeq$Reads[i:(i + 2)]))
weights <- data.frame(
  g1 = MiSeq$Reads[grepl("^F", MiSeq$Well)],
  g2 = MiSeq$Reads[grepl("^G", MiSeq$Well)],
  g3 = MiSeq$Reads[grepl("^H", MiSeq$Well)],
  bi = weights
)
dbx1b.weights <- normaliseColumns(weights[1:6, ])
control.weights <- normaliseColumns(weights[7:9, ])
dbx1b.frameshift.weighted <- dbx1b.frameshift * dbx1b.weights
dbx1b.modified.weighted <- dbx1b.modified * dbx1b.weights[, -4]
df.mean.getDistribution(dbx1b.modified)
df.mean.getDistribution(dbx1b.frameshift)
t.test(dbx1b.frameshift$bi, control.frameshift$bi)
dCohens(dbx1b.frameshift$bi, control.frameshift$bi)
t.test(dbx1b.frameshift$g2, control.frameshift$g2)
dCohens(dbx1b.frameshift$g2, control.frameshift$g2)
t.test(dbx1b.modified$g2, control.modified$g2)
dCohens(dbx1b.modified$g2, control.modified$g2)


png("dbx1b.modified.png", width = 960, height = 1440, res = 300)
boxplot(dbx1b.modified, col = colour, xlab = "CRISPR guide", ylab = "Proportion of modified haplotypes")
addSignificanceMarks(dbx1b.modified, control.modified)
dev.off()

png("dbx1b.frameshift.png", width = 960, height = 1440, res = 300)
boxplot(dbx1b.frameshift[, -4], col = colour, xlab = "CRISPR guide", ylab = "Proportion of frameshifted haplotypes")
addSignificanceMarks(dbx1b.frameshift[, -4], control.frameshift[, -4])
dev.off()

png("dbx1b.bi.png", width = 720, height = 1440, res = 300)
boxplot(dbx1b.frameshift[, 4], xlab = "", ylab = "Probability of biallelic kouckout") # col = colour[4],
# text(x = 1, y = 1.015, labels = "ns", pos = 3, xpd = TRUE)
dev.off()

# Connected plots to see correlation between efficacy among guides and embryo
png("dbx1b.modified.png", width = 960, height = 1440, res = 300)
plot(c(1, 2, 3), dbx1b.modified[1, -4], col = colour[1], type = "b", ylim = c(0, 1), xlab = "CRISPR guide", ylab = "Proportion of modified haplotypes", xaxt = "n")
axis(side = 1, at = c(1, 2, 3), labels = c("g1", "g2", "g3"))
lines(c(1, 2, 3), dbx1b.modified[2, -4], col = colour[2], type = "b")
lines(c(1, 2, 3), dbx1b.modified[3, -4], col = colour[3], type = "b")
lines(c(1, 2, 3), dbx1b.modified[4, -4], col = colour[4], type = "b")
lines(c(1, 2, 3), dbx1b.modified[5, -4], col = colour[6], type = "b")
lines(c(1, 2, 3), dbx1b.modified[6, -4], col = colour[8], type = "b")
dev.off()

png("dbx1b.frameshift.png", width = 960, height = 1440, res = 300)
plot(c(1, 2, 3), dbx1b.frameshift[1, -4], col = colour[1], type = "b", ylim = c(0, 1), xlab = "CRISPR guide", ylab = "Proportion of frameshifted haplotypes", xaxt = "n")
axis(side = 1, at = c(1, 2, 3), labels = c("g1", "g2", "g3"))
lines(c(1, 2, 3), dbx1b.frameshift[2, -4], col = colour[2], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[3, -4], col = colour[3], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[4, -4], col = colour[4], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[5, -4], col = colour[6], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[6, -4], col = colour[8], type = "b")
dev.off()

png("dbx1b.modified.png", width = 960, height = 1440, res = 300)
boxplot(dbx1b.modified, xlab = "CRISPR guide", ylab = "Proportion of modified haplotypes", xata = "n")
addSignificanceMarks(dbx1b.modified, control.modified)
axis(side = 1, at = c(1, 2, 3), labels = c("g1", "g2", "g3"))
lines(c(1, 2, 3), dbx1b.modified[1, -4], col = colour[1], type = "b")
lines(c(1, 2, 3), dbx1b.modified[2, -4], col = colour[2], type = "b")
lines(c(1, 2, 3), dbx1b.modified[3, -4], col = colour[3], type = "b")
lines(c(1, 2, 3), dbx1b.modified[4, -4], col = colour[4], type = "b")
lines(c(1, 2, 3), dbx1b.modified[5, -4], col = colour[6], type = "b")
lines(c(1, 2, 3), dbx1b.modified[6, -4], col = colour[8], type = "b")
dev.off()

png("dbx1b.frameshift.png", width = 960, height = 1440, res = 300)
boxplot(dbx1b.frameshift[, -4], xlab = "CRISPR guide", ylab = "Proportion of frameshifted haplotypes", xata = "n")
addSignificanceMarks(dbx1b.frameshift[, -4], control.frameshift[, -4])
axis(side = 1, at = c(1, 2, 3), labels = c("g1", "g2", "g3"))
lines(c(1, 2, 3), dbx1b.frameshift[1, -4], col = colour[1], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[2, -4], col = colour[2], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[3, -4], col = colour[3], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[4, -4], col = colour[4], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[5, -4], col = colour[6], type = "b")
lines(c(1, 2, 3), dbx1b.frameshift[6, -4], col = colour[8], type = "b")
dev.off()

exp1 <- read.csv("./preliminary.csv") # First experiment aHuC/D aGFP immuno, dbx1b and lrp5 CRISPR, 96 hpf
crtl.exp1.l <- exp1$Volume.L[exp1$Test == "Uninjected"]
crtl.exp1.r <- exp1$Volume.R[exp1$Test == "Uninjected"]
crtl.exp1.tot <- crtl.exp1.l + crtl.exp1.r
crtl.exp1 <- cbind(L = crtl.exp1.l, R = crtl.exp1.r)
dbx1b.exp1.l <- exp1$Volume.L[exp1$Test == "dbx1b"]
dbx1b.exp1.r <- exp1$Volume.R[exp1$Test == "dbx1b"]
dbx1b.exp1.tot <- dbx1b.exp1.l + dbx1b.exp1.r
lrp5.exp1.l <- exp1$Volume.L[exp1$Test == "lrp5"]
lrp5.exp1.r <- exp1$Volume.R[exp1$Test == "lrp5"]
lrp5.exp1.tot <- cbind(L = lrp5.exp1.l, R = lrp5.exp1.r)
mean.getDistribution(dbx1b.exp1.tot)
mean.getDistribution(crtl.exp1.tot)
mean(lrp5.exp1.l + lrp5.exp1.r)
t.test(dbx1b.exp1.tot, mu = lrp5.exp1.l + lrp5.exp1.r)

dHb <- read_excel("dHb.volume.xlsx") # aCachd1 aGFP immuno, dbx1b CRISPR, 48/72/96 hpf
dHb$AI <- (dHb$Volume.L - dHb$Volume.R) / (dHb$Volume.L + dHb$Volume.R)
dHb$Total <- dHb$Volume.L + dHb$Volume.R
crtl.4dpf.l <- dHb$Volume.L[dHb$Embryo == "Uninjected.4dpf"]
crtl.4dpf.r <- dHb$Volume.R[dHb$Embryo == "Uninjected.4dpf"]
crtl.4dpf.tot <- crtl.4dpf.l + crtl.4dpf.r
crtl.4dpf <- cbind(L = crtl.4dpf.l, R = crtl.4dpf.r)
dbx1b.4dpf.l <- dHb$Volume.L[dHb$Embryo == "dbx1b.4dpf"]
dbx1b.4dpf.r <- dHb$Volume.R[dHb$Embryo == "dbx1b.4dpf"]
dbx1b.4dpf.tot <- dbx1b.4dpf.l + dbx1b.4dpf.r
dbx1b.4dpf <- cbind(L = dbx1b.4dpf.l, R = dbx1b.4dpf.r)

crtl.3dpf.l <- dHb$Volume.L[dHb$Embryo == "Uninjected.3dpf"]
crtl.3dpf.r <- dHb$Volume.R[dHb$Embryo == "Uninjected.3dpf"]
crtl.3dpf.tot <- crtl.3dpf.l + crtl.3dpf.r
crtl.3dpf <- cbind(L = crtl.3dpf.l, R = crtl.3dpf.r)
dbx1b.3dpf.l <- dHb$Volume.L[dHb$Embryo == "dbx1b.3dpf"]
dbx1b.3dpf.r <- dHb$Volume.R[dHb$Embryo == "dbx1b.3dpf"]
dbx1b.3dpf.tot <- dbx1b.3dpf.l + dbx1b.3dpf.r
dbx1b.3dpf <- cbind(L = dbx1b.3dpf.l, R = dbx1b.3dpf.r)

crtl.2dpf.l <- dHb$Volume.L[dHb$Embryo == "Uninjected.2dpf"]
crtl.2dpf.r <- dHb$Volume.R[dHb$Embryo == "Uninjected.2dpf"]
crtl.2dpf.tot <- crtl.2dpf.l + crtl.2dpf.r
crtl.2dpf <- cbind(L = crtl.2dpf.l, R = crtl.2dpf.r)
dbx1b.2dpf.l <- dHb$Volume.L[dHb$Embryo == "dbx1b.2dpf"]
dbx1b.2dpf.r <- dHb$Volume.R[dHb$Embryo == "dbx1b.2dpf"]
dbx1b.2dpf.tot <- dbx1b.2dpf.l + dbx1b.2dpf.r
dbx1b.2dpf <- cbind(L = dbx1b.2dpf.l, R = dbx1b.2dpf.r)

lov <- read_excel("leftover.volume.xlsx")
lov$AI <- (lov$Volume.L - lov$Volume.R) / (lov$Volume.L + lov$Volume.R)
lov$Total <- lov$Volume.L + lov$Volume.R
crtl.lov.4dpf.l <- lov$Volume.L[lov$Embryo == "Uninjected.4dpf"]
crtl.lov.4dpf.r <- lov$Volume.R[lov$Embryo == "Uninjected.4dpf"]
crtl.lov.4dpf <- data.frame(L = crtl.lov.4dpf.l, R = crtl.lov.4dpf.r, Embryo = "Uninjected")
crtl.lov.4dpf.tot <- crtl.lov.4dpf.l + crtl.lov.4dpf.r
crtl.lov.4dpf.a <- (crtl.lov.4dpf.l - crtl.lov.4dpf.r) / (crtl.lov.4dpf.l + crtl.lov.4dpf.r)
dbx1b.lov.4dpf.l <- lov$Volume.L[lov$Embryo == "dbx1b.4dpf"]
dbx1b.lov.4dpf.r <- lov$Volume.R[lov$Embryo == "dbx1b.4dpf"]
dbx1b.lov.4dpf <- data.frame(L = dbx1b.lov.4dpf.l, R = dbx1b.lov.4dpf.r, Embryo = "dbx1b")
dbx1b.lov.4dpf.tot <- dbx1b.lov.4dpf.l + dbx1b.lov.4dpf.r
dbx1b.lov.4dpf.a <- (dbx1b.lov.4dpf.l - dbx1b.lov.4dpf.r) / (dbx1b.lov.4dpf.l + dbx1b.lov.4dpf.r)
lov.4dpf.length <- max(c(length(crtl.lov.4dpf.a), length(dbx1b.lov.4dpf.a)))
length(crtl.lov.4dpf.a) <- lov.4dpf.length
length(dbx1b.lov.4dpf.a) <- lov.4dpf.length
lov.4dpf.a <- cbind(Uninjected = crtl.lov.4dpf.a, dbx1b = dbx1b.lov.4dpf.a)
lov.4dpf <- rbind(crtl.lov.4dpf, dbx1b.lov.4dpf)
t.test(dbx1b.lov.4dpf.a, mu = 0) # Test if dbx1b dHbL is still asymmetric baised to the left
# png("box.4dpf.png", width = 960, height = 1440, res = 300)
# boxplot(lov.4dpf.a, col = colour, main = "Asymmetry Index")
# text(x = 1.5, y = max(lov.4dpf.a, na.rm = TRUE) * 1.002, labels = "***", pos = 3, xpd = TRUE)
# boxplot(AI ~ Embryo, lov[grepl("4dpf", lov$Embryo), ], col = colour, main = "Asymmetry Index")
# dev.off()
png("sp.4dpf.png", width = 1920, height = 960, res = 300)
plot(crtl.lov.4dpf.l, crtl.lov.4dpf.r, col = colour[1], xlim = c(0, max(crtl.lov.4dpf.l)), ylim = c(0, max(crtl.lov.4dpf.r)), xlab = expression("Volume of left dHb"[L]), ylab = expression("Volume of right dHb"[L]), pch = 16)
points(dbx1b.lov.4dpf.l, dbx1b.lov.4dpf.r, col = colour[2], pch = 16)
lines(seq(0, max(crtl.lov.4dpf.l, crtl.lov.4dpf.r), length.out = 100), seq(0, max(crtl.lov.4dpf.l, crtl.lov.4dpf.r), length.out = 100), col = "black")
# abline(lm(R ~ L, data = crtl.lov.4dpf), col = colour[1])
# abline(lm(R ~ L, data = dbx1b.lov.4dpf), col = colour[2])
legend("bottomright", fill = colour, legend = c("Uninjected", "dbx1b"), bty = "n", bg = "white")
dev.off()

crtl.lov.3dpf.l <- lov$Volume.L[lov$Embryo == "Uninjected.3dpf"]
crtl.lov.3dpf.r <- lov$Volume.R[lov$Embryo == "Uninjected.3dpf"]
crtl.lov.3dpf <- data.frame(L = crtl.lov.3dpf.l, R = crtl.lov.3dpf.r, Embryo = "Uninjected")
crtl.lov.3dpf.tot <- crtl.lov.3dpf.l + crtl.lov.3dpf.r
crtl.lov.3dpf.a <- (crtl.lov.3dpf.l - crtl.lov.3dpf.r) / (crtl.lov.3dpf.l + crtl.lov.3dpf.r)
dbx1b.lov.3dpf.l <- lov$Volume.L[lov$Embryo == "dbx1b.3dpf"]
dbx1b.lov.3dpf.r <- lov$Volume.R[lov$Embryo == "dbx1b.3dpf"]
dbx1b.lov.3dpf <- data.frame(L = dbx1b.lov.3dpf.l, R = dbx1b.lov.3dpf.r, Embryo = "Uninjected")
dbx1b.lov.3dpf.tot <- dbx1b.lov.3dpf.l + dbx1b.lov.3dpf.r
dbx1b.lov.3dpf.a <- (dbx1b.lov.3dpf.l - dbx1b.lov.3dpf.r) / (dbx1b.lov.3dpf.l + dbx1b.lov.3dpf.r)
lov.3dpf.length <- max(c(length(crtl.lov.3dpf.a), length(dbx1b.lov.3dpf.a)))
length(crtl.lov.3dpf.a) <- lov.3dpf.length
length(dbx1b.lov.3dpf.a) <- lov.3dpf.length
lov.3dpf.a <- cbind(Uninjected = crtl.lov.3dpf.a, dbx1b = dbx1b.lov.3dpf.a)
lov.3dpf <- rbind(crtl.lov.3dpf, dbx1b.lov.3dpf)
png("sp.3dpf.png", width = 1920, height = 960, res = 300) # 1440, 1200
plot(crtl.lov.3dpf.l, crtl.lov.3dpf.r, col = colour[1], xlim = c(0, max(crtl.lov.3dpf.l)), ylim = c(0, max(crtl.lov.3dpf.r)), xlab = expression("Volume of left dHb"[L]), ylab = expression("Volume of right dHb"[L]), pch = 16)
points(dbx1b.lov.3dpf.l, dbx1b.lov.3dpf.r, col = colour[2], pch = 16)
lines(seq(0, max(crtl.lov.3dpf.l, crtl.lov.3dpf.r), length.out = 100), seq(0, max(crtl.lov.3dpf.l, crtl.lov.3dpf.r), length.out = 100), col = "black")
# abline(lm(R ~ L, data = crtl.lov.3dpf), col = colour[1])
# abline(lm(R ~ L, data = dbx1b.lov.3dpf), col = colour[2])
# legend("top", fill = colour, legend = c("Uninjected", "dbx1b"), bty = "n", bg = "white")
dev.off()

crtl.lov.2dpf.l <- lov$Volume.L[lov$Embryo == "Uninjected.2dpf"]
crtl.lov.2dpf.r <- lov$Volume.R[lov$Embryo == "Uninjected.2dpf"]
crtl.lov.2dpf <- data.frame(L = crtl.lov.2dpf.l, R = crtl.lov.2dpf.r, Embryo = "Uninjected")
crtl.lov.2dpf.tot <- crtl.lov.2dpf.l + crtl.lov.2dpf.r
crtl.lov.2dpf.a <- (crtl.lov.2dpf.l - crtl.lov.2dpf.r) / (crtl.lov.2dpf.l + crtl.lov.2dpf.r)
dbx1b.lov.2dpf.l <- lov$Volume.L[lov$Embryo == "dbx1b.2dpf"]
dbx1b.lov.2dpf.r <- lov$Volume.R[lov$Embryo == "dbx1b.2dpf"]
dbx1b.lov.2dpf <- data.frame(L = dbx1b.lov.2dpf.l, R = dbx1b.lov.2dpf.r, Embryo = "Uninjected")
dbx1b.lov.2dpf.tot <- dbx1b.lov.2dpf.l + dbx1b.lov.2dpf.r
dbx1b.lov.2dpf.a <- (dbx1b.lov.2dpf.l - dbx1b.lov.2dpf.r) / (dbx1b.lov.2dpf.l + dbx1b.lov.2dpf.r)
lov.2dpf.length <- max(c(length(crtl.lov.2dpf.a), length(dbx1b.lov.2dpf.a)))
length(crtl.lov.2dpf.a) <- lov.2dpf.length
length(dbx1b.lov.2dpf.a) <- lov.2dpf.length
lov.2dpf.a <- cbind(Uninjected = crtl.lov.2dpf.a, dbx1b = dbx1b.lov.2dpf.a)
lov.2dpf <- rbind(crtl.lov.2dpf, dbx1b.lov.2dpf)
png("sp.2dpf.png", width = 1920, height = 960, res = 300)
plot(crtl.lov.2dpf.l, crtl.lov.2dpf.r, col = colour[1], xlim = c(0, max(crtl.lov.2dpf.l)), ylim = c(0, max(crtl.lov.2dpf.r)), xlab = expression("Volume of left dHb"[L]), ylab = expression("Volume of right dHb"[L]), pch = 16)
points(dbx1b.lov.2dpf.l, dbx1b.lov.2dpf.r, col = colour[2], pch = 16)
lines(seq(0, max(crtl.lov.2dpf.l, crtl.lov.2dpf.r), length.out = 100), seq(0, max(crtl.lov.2dpf.l, crtl.lov.2dpf.r), length.out = 100), col = "black")
# abline(lm(R ~ L, data = crtl.lov.2dpf), col = colour[1])
# abline(lm(R ~ L, data = dbx1b.lov.2dpf), col = colour[2])
# legend("top", fill = colour, legend = c("Uninjected", "dbx1b"), bty = "n", bg = "white")
dev.off()

png("sp.png", width = 1440, height = 1920, res = 300)
plot(crtl.lov.4dpf.l, crtl.lov.4dpf.r, col = colour[1], xlim = c(0, max(crtl.lov.4dpf.l)), ylim = c(0, max(crtl.lov.4dpf.r)), xlab = expression("Volume of left dHb"[L]), ylab = expression("Volume of right dHb"[L]), pch = 0)
points(dbx1b.lov.4dpf.l, dbx1b.lov.4dpf.r, col = colour[2], pch = 0)
lines(seq(0, max(crtl.lov.4dpf.l, crtl.lov.4dpf.r), length.out = 100), seq(0, max(crtl.lov.4dpf.l, crtl.lov.4dpf.r), length.out = 100), col = "black")
points(crtl.lov.3dpf.l, crtl.lov.3dpf.r, col = colour[1], pch = 1)
points(dbx1b.lov.3dpf.l, dbx1b.lov.3dpf.r, col = colour[2], pch = 1)
points(crtl.lov.2dpf.l, crtl.lov.2dpf.r, col = colour[1], pch = 2) # #99DBFF colour[1] #006AA3
points(dbx1b.lov.2dpf.l, dbx1b.lov.2dpf.r, col = colour[2], pch = 2) # #BBEFA9 colour[2] #36891A
points(62000, 1500, col = colour[2], pch = 15)
points(62000, 2000, col = colour[2], pch = 16)
points(62000, 2500, col = colour[2], pch = 17)
text(x = 40500, y = 1500, labels = "96 hpf dbx1b", pos = 4, xpd = TRUE)
text(x = 40500, y = 2000, labels = "72 hpf dbx1b", pos = 4, xpd = TRUE)
text(x = 40500, y = 2500, labels = "48 hpf dbx1b", pos = 4, xpd = TRUE)
points(62000, 0, col = colour[1], pch = 16)
points(62000, 500, col = colour[1], pch = 16)
points(62000, 1000, col = colour[1], pch = 17)
text(x = 34000, y = 0, labels = "96 hpf Uninjected", pos = 4, xpd = TRUE)
text(x = 34000, y = 500, labels = "72 hpf Uninjected", pos = 4, xpd = TRUE)
text(x = 34000, y = 1000, labels = "48 hpf Uninjected", pos = 4, xpd = TRUE)
points(mean(crtl.lov.4dpf.l), mean(crtl.lov.4dpf.r), col = colour[1], pch = 15, cex = 2)
points(mean(crtl.lov.3dpf.l), mean(crtl.lov.3dpf.r), col = colour[1], pch = 16, cex = 2)
points(mean(crtl.lov.2dpf.l), mean(crtl.lov.2dpf.r), col = colour[1], pch = 17, cex = 2)
points(mean(dbx1b.lov.4dpf.l), mean(dbx1b.lov.4dpf.r), col = colour[2], pch = 15, cex = 2)
points(mean(dbx1b.lov.3dpf.l), mean(dbx1b.lov.3dpf.r), col = colour[2], pch = 16, cex = 2)
points(mean(dbx1b.lov.2dpf.l), mean(dbx1b.lov.2dpf.r), col = colour[2], pch = 17, cex = 2)
dev.off()

lov.a <- lov[lov$Volume.L != 0, ]
lov.a$Embryo <- as.character(lov.a$Embryo)
lov.a$Day <- substr(lov.a$Embryo, nchar(lov.a$Embryo) - 3, nchar(lov.a$Embryo))
lov.a$Embryo <- substr(lov.a$Embryo, 1, nchar(lov.a$Embryo) - 5)
lov.a$Embryo <- factor(lov.a$Embryo, levels = c("Uninjected", "dbx1b"))
# bartlett.test(lov$AI ~ lov$Embryo) # Test for homogeneity of variances
lov.a.aov <- aov(lov.a$AI ~ lov.a$Embryo * lov.a$Day)
fCohens(lov.a.aov)
anova(lov.a.aov)
posthocTGH(y = lov[lov$Volume.L != 0, ]$AI, x = lov[lov$Volume.L != 0, ]$Embryo)
t.test(dbx1b.lov.3dpf.a, crtl.lov.3dpf.a)
dCohens(dbx1b.lov.3dpf.a, crtl.lov.3dpf.a)
t.test(dbx1b.lov.4dpf.a, crtl.lov.4dpf.a)
dCohens(dbx1b.lov.4dpf.a, crtl.lov.4dpf.a)

# oneway.test(lov.a$AI ~ lov.a$Embryo * lov.a$Day, var.equal = FALSE)
# write.csv(lov.a, "lov.a.csv", row.names = FALSE)

lov.size <- table(lov$Embryo)[c(4, 1, 5, 2, 6, 3)]
lov.size <- paste0("(", lov.size, ")")
lov.t <- lov
lov.t$Embryo <- as.character(lov.t$Embryo)
lov.t$Day <- substr(lov.t$Embryo, nchar(lov.t$Embryo) - 3, nchar(lov.t$Embryo))
lov.t$Embryo <- substr(lov.t$Embryo, 1, nchar(lov.t$Embryo) - 5)
lov.t$Embryo <- factor(lov.t$Embryo, levels = c("Uninjected", "dbx1b"))
lov.tot.aov <- aov(lov.t$Total ~ lov.t$Embryo * lov.t$Day)
fCohens(lov.tot.aov)
anova(lov.tot.aov)
posthocTGH(y = lov$Total, x = lov$Embryo)
t.test(dbx1b.lov.3dpf.tot, dbx1b.lov.4dpf.tot)
dCohens(dbx1b.lov.3dpf.tot, dbx1b.lov.4dpf.tot)
png("dHbL.png", width = 1440, height = 1440, res = 300)
boxplot(Total ~ Embryo + Day, lov.a, col = colour[1:2], xlab = "Age of embryo", ylab = expression("Total volume of dHb"[L]), at = c(0, 1, 2.5, 3.5, 5, 6), xaxt = "n", ylim = c(0, 90000))
axis(side = 1, at = c(0.5, 3, 5.5), labels = c("48 hpf", "72 hpf", "96 hpf"))
# legend("topleft", fill = colour, legend = c("Uninjected", "dbx1b"), bty = "n")
lines(c(0, 1), rep(80000, 2), lwd = 2, col = colour[5])
text(x = 0.5, y = 76000, labels = "***", pos = 3, xpd = TRUE, col = colour[5])
lines(c(2.5, 3.5), rep(80000, 2), lwd = 2, col = colour[5])
text(x = 3, y = 76000, labels = "***", pos = 3, xpd = TRUE, col = colour[5])
lines(c(5, 6), rep(80000, 2), lwd = 2, col = colour[5])
text(x = 5.5, y = 76000, labels = "***", pos = 3, xpd = TRUE, col = colour[5])
lines(c(2.5, 5), rep(84500, 2), lwd = 2, col = colour[1])
text(x = 3.75, y = 80500, labels = "**", pos = 3, xpd = TRUE, col = colour[1])
lines(c(3.5, 6), rep(89000, 2), lwd = 2, col = colour[2])
text(x = 4.75, y = 87000, labels = "ns", pos = 3, xpd = TRUE, col = colour[2])
lines(c(0, 2.5), rep(85500, 2), lwd = 2, col = colour[1])
text(x = 1.25, y = 81500, labels = "***", pos = 3, xpd = TRUE, col = colour[1])
lines(c(1, 3.5), rep(90000, 2), lwd = 2, col = colour[2])
text(x = 2.25, y = 86000, labels = "***", pos = 3, xpd = TRUE, col = colour[2])
text(c(0, 1, 2.5, 3.5, 5, 6), par("usr")[3] * 0.8, lov.size, pos = 1, xpd = TRUE)
text(c(3), par("usr")[3] * 4.2, "(Sample size)", pos = 1, xpd = TRUE)
dev.off()

png("asymmetry.png", width = 1440, height = 1440, res = 300)
boxplot(AI ~ Embryo + Day, lov.a, col = colour[1:2], xlab = "Age of embryo", ylab = expression("Asymmetry index of dHb"[L]), at = c(0, 1, 2.5, 3.5, 5, 6), xaxt = "n", ylim = c(-0.04, 1.14))
axis(side = 1, at = c(0.5, 3, 5.5), labels = c("48 hpf", "72 hpf", "96 hpf"))
legend("bottomleft", fill = colour, legend = c("Uninjected", "dbx1b"), bty = "n")
lines(c(0, 1), rep(1.03, 2), lwd = 2, col = colour[5])
text(x = 0.5, y = 0.975, labels = "**", pos = 3, xpd = TRUE, col = colour[5])
lines(c(2.5, 3.5), rep(1.03, 2), lwd = 2, col = colour[5])
text(x = 3, y = 0.975, labels = "**", pos = 3, xpd = TRUE, col = colour[5])
lines(c(5, 6), rep(1.03, 2), lwd = 2, col = colour[5])
text(x = 5.5, y = 0.975, labels = "***", pos = 3, xpd = TRUE, col = colour[5])
lines(c(2.5, 5), rep(1.08, 2), lwd = 2, col = colour[1])
text(x = 3.75, y = 1.025, labels = "**", pos = 3, xpd = TRUE, col = colour[1])
lines(c(3.5, 6), rep(1.13, 2), lwd = 2, col = colour[2])
text(x = 4.75, y = 1.095, labels = "ns", pos = 3, xpd = TRUE, col = colour[2])
lines(c(0, 2.5), rep(1.09, 2), lwd = 2, col = colour[1])
text(x = 1.25, y = 1.055, labels = "ns", pos = 3, xpd = TRUE, col = colour[1])
lines(c(1, 3.5), rep(1.14, 2), lwd = 2, col = colour[2])
text(x = 2.25, y = 1.085, labels = "***", pos = 3, xpd = TRUE, col = colour[2])
# points(AI ~ 0, lov.a[lov$Embryo == "Uninjected.4dpf", ], col = colour[1], pch = 16)
text(c(0, 1, 2.5, 3.5, 5, 6), par("usr")[3] * 0.8, lov.size, pos = 1, xpd = TRUE)
text(c(3), par("usr")[3] * 2.7, "(Sample size)", pos = 1, xpd = TRUE)
dev.off()

dHb.size <- table(dHb$Embryo)[c(4, 1, 5, 2, 6, 3)]
dHb.size <- paste0("(", dHb.size, ")")
dHb.v <- dHb
dHb.v$Embryo <- as.character(dHb.v$Embryo)
dHb.v$Day <- substr(dHb.v$Embryo, nchar(dHb.v$Embryo) - 3, nchar(dHb.v$Embryo))
dHb.v$Embryo <- substr(dHb.v$Embryo, 1, nchar(dHb.v$Embryo) - 5)
dHb.v$Embryo <- factor(dHb.v$Embryo, levels = c("Uninjected", "dbx1b"))
# bartlett.test(dHb$Total ~ dHb$Embryo) # Test for homogeneity of variances
oneway.test(dHb.v$Total ~ dHb.v$Embryo * dHb.v$Day, var.equal = FALSE)
dHb.aov <- aov(dHb.v$Total ~ dHb.v$Embryo * dHb.v$Day)
fCohens(dHb.aov)
anova(dHb.aov)
posthocTGH(y = dHb$Total, x = dHb$Embryo)
mean.getDistribution(dbx1b.2dpf.tot)
mean.getDistribution(crtl.2dpf.tot)
mean.getDistribution(dbx1b.3dpf.tot)
mean.getDistribution(crtl.3dpf.tot)
mean.getDistribution(dbx1b.4dpf.tot)
mean.getDistribution(crtl.4dpf.tot)
t.test(dbx1b.2dpf.tot, crtl.2dpf.tot)
dCohens(dbx1b.2dpf.tot, crtl.2dpf.tot)
t.test(dbx1b.3dpf.tot, dbx1b.4dpf.tot)
dCohens(dbx1b.3dpf.tot, dbx1b.4dpf.tot)
t.test(crtl.3dpf.tot, crtl.4dpf.tot)
dCohens(crtl.3dpf.tot, crtl.4dpf.tot)
oneway.test(dHb.v$AI ~ dHb.v$Embryo * dHb.v$Day, var.equal = FALSE)
dHb.a.aov <- aov(dHb.v$AI ~ dHb.v$Embryo * dHb.v$Day)
fCohens(dHb.a.aov)
anova(dHb.a.aov)
posthocTGH(y = dHb$AI, x = dHb$Embryo)

png("dHb.png", width = 1440, height = 1440, res = 300)
boxplot(Total ~ Embryo + Day, dHb.v, col = colour[1:2], xlab = "Age of embryo", ylab = expression("Total volume of dHb (μm³)"), at = c(0, 1, 2.5, 3.5, 5, 6), xaxt = "n", ylim = c(0, max(dHb.v$Total) + 26000))
axis(side = 1, at = c(0.5, 3, 5.5), labels = c("48 hpf", "72 hpf", "96 hpf"))
# legend("topleft", fill = colour, legend = c("Uninjected", "dbx1b"), bty = "n")
lines(c(0, 1), rep(151000, 2), lwd = 2, col = colour[5])
text(x = 0.5, y = 143000, labels = "***", pos = 3, xpd = TRUE, col = colour[5])
lines(c(2.5, 3.5), rep(151000, 2), lwd = 2, col = colour[5])
text(x = 3, y = 143000, labels = "***", pos = 3, xpd = TRUE, col = colour[5])
lines(c(5, 6), rep(151000, 2), lwd = 2, col = colour[5])
text(x = 5.5, y = 143000, labels = "***", pos = 3, xpd = TRUE, col = colour[5])
lines(c(2.5, 5), rep(158000, 2), lwd = 2, col = colour[1])
text(x = 3.75, y = 150000, labels = "***", pos = 3, xpd = TRUE, col = colour[1])
lines(c(3.5, 6), rep(165000, 2), lwd = 2, col = colour[2])
text(x = 4.75, y = 161000, labels = "ns", pos = 3, xpd = TRUE, col = colour[2])
lines(c(0, 2.5), rep(160000, 2), lwd = 2, col = colour[1])
text(x = 1.25, y = 152000, labels = "***", pos = 3, xpd = TRUE, col = colour[1])
lines(c(1, 3.5), rep(167000, 2), lwd = 2, col = colour[2])
text(x = 2.25, y = 159000, labels = "***", pos = 3, xpd = TRUE, col = colour[2])
text(c(0, 1, 2.5, 3.5, 5, 6), par("usr")[3] * 0.8, dHb.size, pos = 1, xpd = TRUE)
text(c(3), par("usr")[3] * 4.2, "(Sample size)", pos = 1, xpd = TRUE)
dev.off()

png("asymmetry.dHb.png", width = 1440, height = 1440, res = 300)
boxplot(AI ~ Embryo + Day, dHb.v, col = colour[1:2], xlab = "Age of embryo", ylab = expression("Asymmetry index of dHb"), at = c(0, 1, 2.5, 3.5, 5, 6), xaxt = "n")
axis(side = 1, at = c(0.5, 3, 5.5), labels = c("48 hpf", "72 hpf", "96 hpf"))
legend("topright", fill = colour, legend = c("Uninjected", "dbx1b"), bty = "n")
text(c(0, 1, 2.5, 3.5, 5, 6), par("usr")[3] * 0.98, dHb.size, pos = 1, xpd = TRUE)
text(c(3), par("usr")[3] * 1.35, "(Sample size)", pos = 1, xpd = TRUE)
dev.off()

suppressWarnings(pHH3 <- read_excel("pHH3.xlsx", col_types = c("text", "skip", "numeric", "numeric")))
pHH3$Embryo <- factor(pHH3$Embryo, levels = c("Uninjected", "dbx1b"))
pHH3$tot <- pHH3$Count.L + pHH3$Count.R
pHH3.size <- table(pHH3$Embryo)
pHH3.size <- paste0("(", pHH3.size, ")")
crtl.pHH3.tot <- pHH3$tot[pHH3$Embryo == "Uninjected"]
dbx1b.pHH3.tot <- pHH3$tot[pHH3$Embryo == "dbx1b"]
t.test(dbx1b.pHH3.tot, crtl.pHH3.tot)
dCohens(dbx1b.pHH3.tot, crtl.pHH3.tot)
png("pHH3.png", width = 1440, height = 1440, res = 300)
boxplot(tot ~ Embryo, pHH3, col = colour, xlab = "Embryo", ylab = "Number of pHH3 positive cells")
legend("topleft", fill = colour, legend = c("Uninjected", "dbx1b"), bty = "n")
text(c(1, 2), par("usr")[3] * 0.98, pHH3.size, pos = 1, xpd = TRUE)
text(c(1.5), par("usr")[3] * 0.65, "(Sample size)", pos = 1, xpd = TRUE)
dev.off()

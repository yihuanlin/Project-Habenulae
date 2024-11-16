#!/usr/bin/env RScript
rm(list = ls(all = TRUE)) # Remove previous objects from globe environment
tryCatch(assign("last.warning", NULL, envir = baseenv()), error = function(e) {}) # Remove previous warnings # nolint: line_length_linter.
tryCatch(dev.off(dev.list()["RStudioGD"]), error = function(e) {}) # Remove previous plots in RStudio

translate <- function(dna_sequence) {
  codons <- strsplit(dna_sequence, "(?<=.{3})", perl = TRUE)[[1]]
  amino_acids <- c()

  for (codon in codons) {
    amino_acid <- switch(codon,
      "TTT" = "F",
      "TTC" = "F",
      "TTA" = "L",
      "TTG" = "L",
      "CTT" = "L",
      "CTC" = "L",
      "CTA" = "L",
      "CTG" = "L",
      "ATT" = "I",
      "ATC" = "I",
      "ATA" = "I",
      "ATG" = "M",
      "GTT" = "V",
      "GTC" = "V",
      "GTA" = "V",
      "GTG" = "V",
      "TCT" = "S",
      "TCC" = "S",
      "TCA" = "S",
      "TCG" = "S",
      "CCT" = "P",
      "CCC" = "P",
      "CCA" = "P",
      "CCG" = "P",
      "ACT" = "T",
      "ACC" = "T",
      "ACA" = "T",
      "ACG" = "T",
      "GCT" = "A",
      "GCC" = "A",
      "GCA" = "A",
      "GCG" = "A",
      "TAT" = "Y",
      "TAC" = "Y",
      "TAA" = "*",
      "TAG" = "*",
      "CAT" = "H",
      "CAC" = "H",
      "CAA" = "Q",
      "CAG" = "Q",
      "AAT" = "N",
      "AAC" = "N",
      "AAA" = "K",
      "AAG" = "K",
      "GAT" = "D",
      "GAC" = "D",
      "GAA" = "E",
      "GAG" = "E",
      "TGT" = "C",
      "TGC" = "C",
      "TGA" = "*",
      "TGG" = "W",
      "CGT" = "R",
      "CGC" = "R",
      "CGA" = "R",
      "CGG" = "R",
      "AGT" = "S",
      "AGC" = "S",
      "AGA" = "R",
      "AGG" = "R",
      "GGT" = "G",
      "GGC" = "G",
      "GGA" = "G",
      "GGG" = "G"
    )
    amino_acids <- c(amino_acids, amino_acid)
  }

  return(paste(amino_acids, collapse = ""))
}

highlight <- function(seq_mutated, cat) {
  cat(paste0("\033[32m", cat, "\033[0m\n"))
  highlighted_seq <- ""
  stop <- FALSE
  count <- 0

  for (i in 1:nchar(seq)) {
    if (stop) {
      count <- count + 1
      highlighted_seq <- paste0(highlighted_seq, "\033[31m", substr(seq, i, i), "\033[0m")
    } else if (substr(seq_mutated, i, i) == substr(seq, i, i)) {
      highlighted_seq <- paste0(highlighted_seq, substr(seq_mutated, i, i))
    } else if (substr(seq_mutated, i, i) == "*") {
      stop <- TRUE
      highlighted_seq <- paste0(highlighted_seq, "\033[34m", substr(seq, i, i), "(", i, ")->", substr(seq_mutated, i, i), "\033[0m")
    } else {
      highlighted_seq <- paste0(highlighted_seq, "\033[34m", substr(seq, i, i), "(", i, ")->", substr(seq_mutated, i, i), "\033[0m")
    }
  }

  cat(highlighted_seq, "\n")
  if (count > 5) {
    cat("Truncated:", count, "\n")
  }
}

if (requireNamespace("rstudioapi", quietly = TRUE)) {
  is_rstudio <- TRUE
} else {
  is_rstudio <- FALSE
}
if (is_rstudio) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

seq <- readLines("./Danio_rerio_dbx1b_sequence.fa", warn = FALSE) # Load coding sequence of Dbx1b
seq <- paste(seq[-1], collapse = "")

seq_mutated <- gsub("CGTCGCCATGTCTACTTCTT", "CGTCGCCATGTCTACTGCTT", seq)
seq_mutated <- translate(seq_mutated)

seq_mutated_1 <- gsub("TCCGAGCTTTCTGCGGCCCTCGTCGGCTCT", "TCCGAGCTTTCTGCGGCCCGAGCTTTGCTGCGTCGGCTCT", seq)
seq_mutated_1 <- gsub("TCACCACCACCTCCCGTGCAAGGAATGAATGCAAA", "TCACCACCACCTCCCAAGGAATGAATGCAAA", seq_mutated_1)
seq_mutated_1 <- translate(seq_mutated_1)

seq_mutated_2 <- gsub("TCCGAGCTTTCTGCGGCCCTCGTCGGCTCTTTCCTTGCC", "TCCGAGCTTTCTGCGGCCCTTCGTCGGCTCTTTCCTTGCC", seq)
seq_mutated_2 <- gsub("TCACCACCACCTCCCGTGCAAGGAATGAATGCAAA", "TCACCACCACCTCCCAAGGAATGAATGCAAA", seq_mutated_2)
seq_mutated_2 <- translate(seq_mutated_2)

seq_mutated_3 <- gsub("TCCGAGCTTTCTGCGGCCCTCGTCGGCTCTTTCCTTGCCT", "TCCGAGCTTTCTGCGGCTCTTTCCTTGCCT", seq)
seq_mutated_3 <- gsub("TCACCACCACCTCCCGTGCAAGGAATGAATGCAAA", "TCACCACCACCTCCAAGGAATGAATGCAAA", seq_mutated_3)
seq_mutated_3 <- translate(seq_mutated_3)

seq_mutated_4 <- gsub("TCCGAGCTTTCTGCGGCCCTCGTCGGCTCTTTCCTTGCCT", "TCCGAGCTTTCTGCGGCTCTTTCCTTGCCT", seq)
seq_mutated_4 <- gsub("TCACCACCACCTCCCGTGCAAGGAATGAATGCAAA", "TCAAGGAATGAATGCAAA", seq_mutated_4)
seq_mutated_4 <- translate(seq_mutated_4)

seq_mutated_6 <- gsub("TCCGAGCTTTCTGCGGCCCTCGTCGGCTCTTTCCTTGC", "TCCGAGCTTTCTGCGGCCCTTTCGTCGGCTCTTTCCTTGC", seq)
seq_mutated_6A <- gsub("TCACCACCACCTCCCGTGCAAGGAATGAATGCAAA", "TCACCACCACCTGCAAGGAATGAATGCAAA", seq_mutated_6)
seq_mutated_6B <- gsub("CGTCGCCATGTCTACTTCTTTGG", "CGTCGCCATGTCTACTGG", seq_mutated_6)
seq_mutated_6A <- translate(seq_mutated_6A)
seq_mutated_6B <- translate(seq_mutated_6B)

seq <- translate(seq)

highlight(seq_mutated_1, "H1+F1: (8 sg1 Shift +1 X 52 sg2 Shift -1 not effective)")
highlight(seq_mutated_2, "H2+F2: (1 sg1 Shift +1 X 9 sg2 Shift -1 not effective)")
highlight(seq_mutated_3, "H3+F3: (sg1 Shift -1 effective)")
highlight(seq_mutated_4, "H4+F4: (sg1 Shift -1 effective)")
cat(paste0("\033[32m", "H5+F5: sg3 no frameshift", "\033[0m\n"))
highlight(seq_mutated_6A, "H6+F6: (sg1 Shift -1 + sg3 Shift +1 effective)")
highlight(seq_mutated_6B, "H6+F6: (sg1 Shift -1 + sg2 Shift +1 effective)")

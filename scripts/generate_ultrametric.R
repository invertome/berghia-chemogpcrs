#!/usr/bin/env Rscript
# generate_ultrametric.R
# Purpose: Convert a ML phylogenetic tree to ultrametric using penalized likelihood.
# Inputs: ML tree file (Newick format), output file path
# Outputs: Ultrametric tree file (Newick format)
# Method: Uses chronos from ape package with penalized likelihood
# Author: Jorge L. Perez-Moreno, Ph.D., Katz Lab, University of Massachusetts, Amherst
#
# Usage:
#   Rscript generate_ultrametric.R input_tree.tre output_ultrametric.tre [model] [lambda]
#
# Arguments:
#   input_tree.tre       - Input ML tree (Newick format, e.g., from IQ-TREE)
#   output_ultrametric.tre - Output ultrametric tree
#   model                - Clock model: "correlated" (default), "relaxed", "discrete"
#   lambda               - Smoothing parameter (default: 1)
#
# Notes:
#   - IMPORTANT: Use BUSCO species tree, NOT GPCR gene tree (to avoid bias)
#   - Default scales root to age 1 (sufficient for CAFE5 relative rates)
#   - No fossil calibration required for relative dating

suppressPackageStartupMessages({
    library(ape)
})

# --- Parse Command Line Arguments ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
    cat("Usage: Rscript generate_ultrametric.R input_tree.tre output_tree.tre [model] [lambda]\n")
    cat("\nModels: correlated (default), relaxed, discrete\n")
    cat("Lambda: smoothing parameter (default: 1)\n")
    quit(status = 1)
}

input_file <- args[1]
output_file <- args[2]
model <- ifelse(length(args) >= 3, args[3], "correlated")
lambda <- ifelse(length(args) >= 4, as.numeric(args[4]), 1)

# --- Validate Input ---
if (!file.exists(input_file)) {
    stop(paste("Input file not found:", input_file))
}

cat("=== Ultrametric Tree Generation ===\n")
cat("Input:", input_file, "\n")
cat("Output:", output_file, "\n")
cat("Model:", model, "\n")
cat("Lambda:", lambda, "\n\n")

# --- Read Tree ---
cat("Reading input tree...\n")
tree <- read.tree(input_file)

# Print tree info
cat("  Tips:", Ntip(tree), "\n")
cat("  Nodes:", Nnode(tree), "\n")
cat("  Rooted:", is.rooted(tree), "\n")

# Check if tree is rooted
if (!is.rooted(tree)) {
    cat("  WARNING: Tree is unrooted. Attempting to root at midpoint...\n")
    tree <- midpoint.root(tree)
    if (!is.rooted(tree)) {
        stop("Failed to root tree. Please provide a rooted tree.")
    }
    cat("  Tree rooted at midpoint.\n")
}

# Check for zero-length branches (can cause issues)
if (any(tree$edge.length == 0)) {
    cat("  WARNING: Tree has zero-length branches. Adding small values...\n")
    tree$edge.length[tree$edge.length == 0] <- 1e-8
}

# Check for negative branches
if (any(tree$edge.length < 0)) {
    cat("  WARNING: Tree has negative branch lengths. Setting to small positive...\n")
    tree$edge.length[tree$edge.length < 0] <- 1e-8
}

# --- Check if Already Ultrametric ---
is_ultra <- is.ultrametric(tree, tol = 1e-6)
cat("  Already ultrametric:", is_ultra, "\n\n")

if (is_ultra) {
    cat("Tree is already ultrametric. Copying to output.\n")
    write.tree(tree, output_file)
    cat("Done.\n")
    quit(status = 0)
}

# --- Convert to Ultrametric using chronos ---
cat("Converting to ultrametric using chronos...\n")
cat("  (This may take a moment for large trees)\n\n")

# Set up calibration: default scales root to age 1
# No fossil calibration needed for relative rate analysis
calibration <- makeChronosCalib(tree)

# Run chronos with specified model
tryCatch({
    if (model == "discrete") {
        # For strict clock, use discrete model with 1 rate category
        control <- chronos.control(nb.rate.cat = 1)
        ultra_tree <- chronos(tree, lambda = lambda, model = "discrete",
                              calibration = calibration, control = control)
    } else {
        # Correlated or relaxed clock
        ultra_tree <- chronos(tree, lambda = lambda, model = model,
                              calibration = calibration)
    }
}, error = function(e) {
    cat("ERROR during chronos:", conditionMessage(e), "\n")
    cat("Attempting with relaxed model...\n")
    ultra_tree <<- chronos(tree, lambda = lambda, model = "relaxed",
                           calibration = calibration)
})

# --- Verify Result ---
cat("\nVerifying ultrametric conversion...\n")

# Check if result is ultrametric
if (!is.ultrametric(ultra_tree, tol = 1e-6)) {
    cat("  WARNING: Result may not be perfectly ultrametric.\n")
    cat("  Attempting to force ultrametric...\n")

    # Force ultrametric by adjusting terminal branches
    root_to_tip <- node.depth.edgelength(ultra_tree)
    max_depth <- max(root_to_tip[1:Ntip(ultra_tree)])

    for (i in 1:Ntip(ultra_tree)) {
        tip_depth <- root_to_tip[i]
        if (tip_depth < max_depth) {
            # Find the edge leading to this tip
            edge_idx <- which(ultra_tree$edge[,2] == i)
            if (length(edge_idx) > 0) {
                ultra_tree$edge.length[edge_idx] <-
                    ultra_tree$edge.length[edge_idx] + (max_depth - tip_depth)
            }
        }
    }
}

# Final check
is_ultra_final <- is.ultrametric(ultra_tree, tol = 1e-6)
cat("  Final ultrametric check:", is_ultra_final, "\n")

# --- Print Statistics ---
cat("\nTree Statistics:\n")
cat("  Tips:", Ntip(ultra_tree), "\n")
cat("  Root age:", max(node.depth.edgelength(ultra_tree)), "\n")

# Branch length summary
cat("  Branch lengths:\n")
cat("    Min:", min(ultra_tree$edge.length), "\n")
cat("    Max:", max(ultra_tree$edge.length), "\n")
cat("    Mean:", mean(ultra_tree$edge.length), "\n")

# --- Write Output ---
cat("\nWriting ultrametric tree to:", output_file, "\n")
write.tree(ultra_tree, output_file)

# Also write in NEXUS format (useful for some downstream tools)
nexus_file <- sub("\\.tre$", ".nex", output_file)
if (nexus_file != output_file) {
    cat("Writing NEXUS format to:", nexus_file, "\n")
    write.nexus(ultra_tree, file = nexus_file)
}

# --- Write Log File ---
log_file <- sub("\\.tre$", "_chronos.log", output_file)
if (log_file != output_file) {
    cat("\nWriting log to:", log_file, "\n")
    sink(log_file)
    cat("Chronos Ultrametric Conversion Log\n")
    cat("===================================\n\n")
    cat("Input:", input_file, "\n")
    cat("Output:", output_file, "\n")
    cat("Model:", model, "\n")
    cat("Lambda:", lambda, "\n")
    cat("Tips:", Ntip(ultra_tree), "\n")
    cat("Is ultrametric:", is_ultra_final, "\n")
    cat("Root age:", max(node.depth.edgelength(ultra_tree)), "\n")
    cat("\nBranch length summary:\n")
    print(summary(ultra_tree$edge.length))
    cat("\nPhylogenetic diversity:", sum(ultra_tree$edge.length), "\n")
    sink()
}

cat("\n=== Conversion Complete ===\n")

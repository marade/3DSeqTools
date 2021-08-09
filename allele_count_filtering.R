
library(evobiR)

args <- commandArgs(trailingOnly = TRUE)
inpath <- args[1]
outfile <- args[2]
if (is.na(inpath) || is.na(outfile)) {
    message("Syntax: Rscript allele_count_filtering.R <input-directory> <output-file>")
    quit()
}

## Processing of 3D-seq replicate data after generation of .tab files by initial
## sequence read processing.
## Overall workflow (additional comments below):
## 1.  get set of .tab data, combine
## 2.  filter for low coverage in individual replicate samples (part of getting .tab data)
## 3.  filter snps (universal or individual)
## 4.  calc mean aaf per site (with minXwYr filtering)
## 5.  filter out isolated sites
## 6.  calc sliding windows

## Get list of 'excluded sites' (e.g., phage hypervariable sites)
excluded_sites <- read.delim('./excluded_sites.txt', header = F)[, 1]

tabulate_all_samples <- function(tabfile_set,
                                 short_ids,
                                 min_coverage_per_sample = 15) {
  
  ## Generates a combined table of ref and alt allele counts for all replicate
  ## samples in a set
  
  ## PARAMETERS:
  ## tabfile_set - list of the .tab files (with paths) for the replicates in
  ##   the set (these are the .tab output files from the preliminary processing)
  ## short_ids - list of short identifiers for the same replicates
  
  out <- NULL
  
  for (i in 1:length(tabfile_set)) {
    print(paste0(" Reading file ", i, " of ", length(tabfile_set)))
    
    s <- read.delim(tabfile_set[i])
    
    # Filter
    s <- s[which(s[, "refseq"] %in% c("TC", "GA")), ]
    s <- s[which(!s[, "position"] %in% excluded_sites), ]
    
    nxt <- s[, c("ref_allele_count", "alt_allele_count")]
    
    # Set ref and alt counts to 0 for sites with total coverage less than min_coverage_per_sample
    low_cov_sites <- apply(nxt, 1, sum) < min_coverage_per_sample
    nxt[low_cov_sites, 'ref_allele_count'] <- 0
    nxt[low_cov_sites, 'alt_allele_count'] <- 0
    
    colnames(nxt) <- c(paste0(short_ids[i], "_ref"), paste0(short_ids[i], "_alt"))
    
    if (is.null(out)) {
      pos <- s[, "position"]
      out <- cbind(pos, nxt)
      colnames(out)[1] <- "position"
    } else {
      out <- cbind(out, nxt)
    }
  }

  return (out)
}

filter_snps <- function(tab,
                        filter_method = "individual",
                        aaf_threshold = 0.95) {
  
  ## Removes "SNPs".  I.e., sets alt and ref counts to 0 at sites where aaf >
  ## aaf_threshold (aaf = alt allele frequency = alt/(ref+alt))
  
  ## PARAMETERS:
  ## tab - output table from tabulate_all_samples()
  ## filter_method - either "individual", in which SNPs are removed
  ##   independently from each replicate sample, or "universal", in which SNPs
  ##   are removed only if the site is classified as a SNP in all replicate samples
  ##   in the set.
  ## aaf_threshold - defines "SNP"
  
  if (filter_method == "individual") {
    print("  Filtering SNPs in individual samples.")
    
    ref_col_nos <- grep('ref', colnames(tab))
    
    for (rcn in ref_col_nos) {
      acn = rcn + 1
      aaf = tab[,acn] / (tab[,rcn] + tab[,acn])
      sample_snps <- which(aaf > aaf_threshold)
      tab[sample_snps, acn] <- 0
      tab[sample_snps, rcn] <- 0
    }
    
  } else if (filter_method == "universal") {
    print("  Filtering universal SNPs.")
    
    ref_col_nos <- grep('ref', colnames(tab))
    alt_col_nos <- ref_col_nos + 1
    
    all_samples_are_snps <- function(tab_row) {
      aafs <- tab_row[alt_col_nos] / (tab_row[ref_col_nos] + tab_row[alt_col_nos])
      return (sum(aafs > aaf_threshold, na.rm=T) == sum(!is.na(aafs)) &
                sum(!is.na(aafs)) > 0)
    }
    
    universal_snps <- which(apply(tab, 1, all_samples_are_snps))
    tab[universal_snps, c(ref_col_nos, alt_col_nos)] <- 0
    
  } else {
    
    print("Set filter_method to 'individual' or 'universal'")
    return(FALSE)
  }
  
  return(tab)
}

calculate_mean_aafs <- function(tab,
                                min_alt_reads_per_rep = 1,
                                min_reps_w_min_alt_reads = 3) {
  
  ## Calculate the mean alt allele frequency from the individual sample aafs.
  ## Set mean aaf to 0 if fewer than min_reps_w_min_alt_reads individual
  ##   replicate samples have at least min_alt_reads_per_rep alt reads.
  ## Return a table with just positions and mean aafs
  
  ## PARAMETERS:
  ## tab - output table from tabulate_all_samples() or filter_snps()

  ref_col_nos <- grep('ref', colnames(tab))
  alt_col_nos <- ref_col_nos + 1
  
  maaf <- function(r) {
    ## r is a row from tab
    ## Return aaf value
    ## Return 0 if not at least min_alt_reads_per_rep alt reads in at least
    ##   min_reps_w_min_alt_reads reps
    
    refs = r[ref_col_nos]
    alts = r[alt_col_nos]
    
    if (sum(alts >= min_alt_reads_per_rep, na.rm=T) >= min_reps_w_min_alt_reads) {
      aafs <- alts / (refs + alts)
      return(mean(aafs, na.rm=T))
    } else {
      return(0)
    }
  }
  
  mean_aaf <- apply(tab, 1, maaf)
  mean_aaf_table <- cbind(tab[,"position"], mean_aaf)
  colnames(mean_aaf_table) <- c("position", "mean_aaf")
  
  return(mean_aaf_table)
}

filter_isolated_sites <- function(maaf_tab,
                                  isolation_bp = 100,
                                  max_sites_within_window = 1) {
  
  ## Remove isolated sites.  I.e., set maaf to 0 for each site with maaf>0 if
  ##   there are not more than max_sites_within_window maaf>0 sites (including the
  ##   site itself) within isolation_bp on either side of the site.
  
  ## PARAMETERS: maaf_tab - a two-column table with position and maaf (or aaf),
  ##   such as an output table from calculate_mean_aafs()
  
  ## Returns the same table but with mean_aaf set to 0 for isolated sites
  
  pos <- maaf_tab[,1]
  maaf <- maaf_tab[,2]
  
  np <- nrow(maaf_tab)
  isol <- rep(F, np)
  
  for (i in which(maaf > 0)) {
    nt <- pos[i]
    # To save processing time, analyze a sub-portion of the full table that
    # includes the site and is at least as large as the isolation span to query.
    # The sub-portion is larger than necessary but much smaller than the full table.
    subportion <- max(c(0, i - isolation_bp)) : min(c(np, i + isolation_bp))
    psub <- pos[subportion]
    asub <- maaf[subportion]
    # calculate if the site is isolated
    isol[i] <- 
      sum(asub[psub > nt - isolation_bp & psub < nt + isolation_bp] > 0) <= max_sites_within_window
  }
  
  maaf[isol] <- 0
  maaf_tab <- cbind(pos, maaf)
  
  return(maaf_tab)
}

calculate_moving_average <- function(maaf_tab,
                                     sliding_window_bp = 75) {
  
  ## Calculate the moving average maaf (or aaf).  I.e., for each position in a
  ##   maaf_tab table, calculate that mean aaf of all positions within a window of
  ##   size sliding_window_bp centered on that position.
  
  ## PARAMETERS: 
  ## maaf_tab - a two-column table with position and maaf (or aaf), such as an
  ##   output table from calculate_mean_aafs() or filter_isolated_sites()
  ## sliding_window_bp - the size of the sliding window. Use an ODD value so
  ##   that the window can be centered on a given position!
  
  ## Returns a table like the input table but with a new column of the moving average
  
  pos <- maaf_tab[,1]
  maaf <- maaf_tab[,2]
  maxpos <- max(pos)
  
  print(paste0("   Calculating moving avg with ", sliding_window_bp, "bp window size"))
  
  half_wind = floor(sliding_window_bp/2)
  if (half_wind == sliding_window_bp/2) { 
    print('WARNING - use odd no of bp for sliding window size!') 
  }
  
  allbases <- rep(NA, maxpos)
  allbases[pos] <- maaf
  allbases_wrapped <- c(tail(allbases, half_wind),
                        allbases,
                        head(allbases, half_wind),
                        NA)  # extra NA so sliding window vector is the correct length
  
  f <- function(x) {
    mean(x, na.rm = T)
  }
  mov_avg <- SlidingWindow(f, allbases_wrapped, sliding_window_bp, 1)
  
  newcolname <- paste0('mov_avg_', sliding_window_bp, 'bp')
  out_colnames <- colnames(maaf_tab)
  out <- cbind(maaf_tab, mov_avg[pos])
  colnames(out) <- c(out_colnames, newcolname)
  
  return(out)
}

## See text of publication for parameters used for different analyses.

tfset <- list.files(inpath), full.names=TRUE)

sids <- c('rep1', 'rep2', 'rep3', 'rep4')

## Combine the replicate data into a single table and minimally filter to remove
## hypervariable sites and sites with low coverage.
combined_table <- tabulate_all_samples(tfset, sids,
                                       min_coverage_per_sample = 15)
## (Save and visualize as needed)

## Filter out transitions that likely represent SNPs in the parent strain(s).
## Filter using "individual" method to remove SNPs from each replicate
## independently.  Filter using "universal" method to remove only SNPs that are
## apparent in all replicates.
table_wo_snps <- filter_snps(combined_table,
                             filter_method = "universal",
                             aaf_threshold = 0.95)

## Filter for transitions that are reproduced in multiple replicates (e.g., 3 of
## the 4) with a minimum read count threshold for mutation events (e.g., 3
## alternate-allele reads). After filtering calculate average transition
## frequencies (mean_aaf).
mean_aaf_table <- calculate_mean_aafs(table_wo_snps,
                                      min_alt_reads_per_rep = 3,
                                      min_reps_w_min_alt_reads = 3)

## Eliminate positions lacking a neighboring transition event within 100 bp.
mean_aaf_table_no_isol <- filter_isolated_sites(mean_aaf_table,
                                                isolation_bp = 100,
                                                max_sites_within_window = 1)

## Calculate a moving average of the mean_aaf values using a 75-bp window.
mean_aaf_w_mov_avg <- calculate_moving_average(mean_aaf_table_no_isol,
                                               sliding_window_bp = 75)
## (Save and visualize as needed)
message("Writing output...")
message(output)
write.table(mean_aaf_w_mov_avg, file=outfile, col.names=NA, quote=FALSE, sep='\t')

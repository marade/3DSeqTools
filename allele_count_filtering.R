
library(evobiR)

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



## Example of processing a data set
## Strain:  gcsR-dddA  Tn7::Para-dddI  delta-ung
## Condition:  no ara, 3 passages

## See text of publication for parameters used for different analyses.

tfset <- c('29.mq30.deduped.raw.tab',
           '30.mq30.deduped.raw.tab',
           '31.mq30.deduped.raw.tab',
           '32.mq30.deduped.raw.tab')
tfset <- paste0('./sample_tab_data/', tfset)

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





#### TEMP
plot_tbl <- function(dat,
                     xlim = c(0, 6264398),
                     ylim = c(0,0.5),
                     use_local_ylim = F,
                     title = "",
                     subtitle = "",
                     site_name = "",
                     span_show_all_zeros = 10000) {
  # (show all 0s for xlim span of span_show_all_zeros (or less); for larger
  # spans, show only a proportional subset of the 0s)
  probcut <- span_show_all_zeros / (xlim[2] - xlim[1])
  subdat <- dat[dat[,"mean_aaf"] > 0 | runif(nrow(dat)) < probcut,]
  pos <- subdat[, "position"]
  maaf <- subdat[ ,"mean_aaf"]
  if (use_local_ylim) {
    ylim <- c(0, max(0.1, max(maaf[pos >= xlim[1] & pos <= xlim[2]], 
                              na.rm=T), 
                     na.rm=T))
  }
  plot(pos, maaf,
       xlab = 'genome position',
       ylab = 'maaf',
       xlim = xlim,
       ylim = ylim, 
       main = title,
       sub = subtitle)
  text(xlim[1] + 0.001 * (xlim[2] - xlim[1]), 
       y = ylim[1] + 0.98 * (ylim[2] - ylim[1]),
       labels = site_name,
       adj = c(0, 1), cex = 1.4)
}

pos <- combined_table[,'position']
dat <- mean_aaf_w_mov_avg
cent <- site_info['GcsR1', 'center']  #2747027
xlim <- cent + c(-48, 48)
plot_tbl(dat, xlim=xlim, ylim = c(0, 0.2))
abline(v=cent, col="Pink")

combined_table[pos >= xlim[1] & pos <= xlim[2],]
combined_table[pos >= xlim[1] & pos <= xlim[2] &
                 mean_aaf_table[,'mean_aaf']>0,]
# sites to use for qPCR:
candidate_sites <- c(2747024, 2747109, 2747205, 2747245)
points(candidate_sites, dat[dat[,'position'] %in% candidate_sites, 'mean_aaf'],
       col="Red")
combined_table[pos %in% candidate_sites,]


## TEMP2 (FleQ)
tfiles <- c('57.mq30.deduped.raw.TC.tab',
            '58.mq30.deduped.raw.TC.tab',
            '17.mq30.deduped.raw.TC.tab',
            '18.mq30.deduped.raw.TC.tab')
tpaths <- paste0('../seq_context_data/', c('NextSeq_Run_03/', 'NextSeq_Run_03/',
                                           'NextSeq_Run_04/', 'NextSeq_Run_04/'))
tfset <- paste0(tpaths, tfiles)

sids <- c('57', '58', '17', '18')

combined_table <- tabulate_all_samples(tfset, sids,
                                       min_coverage_per_sample = 15)

table_wo_snps <- filter_snps(combined_table,
                             filter_method = "individual",
                             aaf_threshold = 0.95)

mean_aaf_table <- calculate_mean_aafs(table_wo_snps,
                                      min_alt_reads_per_rep = 1,
                                      min_reps_w_min_alt_reads = 3)

mean_aaf_w_mov_avg <- calculate_moving_average(mean_aaf_table,
                                               sliding_window_bp = 75)

pos <- combined_table[,'position']
dat <- mean_aaf_w_mov_avg
cent <- 2453403  # pslA
xlim <- cent + c(-270, 270)
plot_tbl(dat, xlim=xlim, ylim = c(0, 0.2))
abline(v=cent, col="Pink")

sub_set <- pos >= xlim[1] & pos <= xlim[2] &
  mean_aaf_table[,'mean_aaf']>-1
combined_table[sub_set,]
points(pos[sub_set], mean_aaf_w_mov_avg[sub_set, 'mean_aaf'], col="Red")

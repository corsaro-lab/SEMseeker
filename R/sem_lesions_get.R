# LESIONS clustering by GENOMIC distance (bp) rather than matrix
# distance (probe count). The legacy sliding_window_size parameter is gone;
# the only knob is now LESIONS_BP, registered in core_init_env() and exposed via
# semseeker(LESIONS_BP = ...).
#
# Semantics of LESIONS_BP: maximum bp distance for two probes to be in the
# same enrichment window. For probe i, ENRICHMENT[i] = sum(MUTATIONS[j]) over
# j with |START[j] - START[i]| <= LESIONS_BP. This makes the LESIONS callset
# comparable across array densities (450K / EPIC / EPICv2) and consistent
# with biology defined in bp (CpG islands ~500-2000bp, DMRs, enhancers).

#' @importFrom dplyr %>%
#' @importFrom rlang .data
sem_lesions_get <- function(grouping_column, mutation_annotated_sorted)
{

  ssEnv <- core_get_session_info()

  if( is.null(mutation_annotated_sorted))
    return (mutation_annotated_sorted)

  if(nrow(mutation_annotated_sorted) == 0)
    return (mutation_annotated_sorted)

  mutationAnnotatedSortedLocal <- mutation_annotated_sorted

  summed <- stats::aggregate(mutationAnnotatedSortedLocal$MUTATIONS, by = list(mutationAnnotatedSortedLocal[,grouping_column]), FUN = sum)
  colnames(summed) <- c(grouping_column,"MUTATIONS_COUNT")
  counted <- stats::aggregate(mutationAnnotatedSortedLocal$MUTATIONS, by = list(mutationAnnotatedSortedLocal[,grouping_column]), FUN = length)
  colnames(counted) <- c(grouping_column,"PROBES_COUNT")
  mutationAnnotatedSortedLocal <- merge(mutationAnnotatedSortedLocal,summed, by = grouping_column)
  mutationAnnotatedSortedLocal <- merge(mutationAnnotatedSortedLocal,counted, by = grouping_column)
  rm(counted)
  rm(summed)

  lesions_bp <- as.integer(ssEnv$LESIONS_BP)
  if (is.na(lesions_bp) || lesions_bp < 0L)
    stop("LESIONS_BP must be a non-negative integer (bp distance), got: ",
         ssEnv$LESIONS_BP)

  # Per-row enrichment / window-size / window-span within a single grouping
  # unit. start, mutations: numeric vectors of equal length (one row per
  # probe). Returns a list with 3 vectors of the same length:
  #   enrichment[i]  = sum(mutations[j]) for j with |start[j]-start[i]|<=lesions_bp
  #   window_size[i] = count of probes within the bp window centered on i
  #   span[i]        = max(start[j]) - min(start[j]) over those j
  # O(N log N) via order + findInterval + cumsum on the SORTED vector.
  .bp_window_stats <- function(start, mutations, lesions_bp) {
    n <- length(start)
    if (n == 0L)
      return(list(enrichment = numeric(0), window_size = integer(0),
                  span = numeric(0)))
    o <- order(start)
    s <- start[o]
    m <- mutations[o]
    # findInterval(x, vec) returns max j s.t. vec[j] <= x (left-continuous).
    # left  = first j with s[j] >= s[i] - lesions_bp
    # right = last  j with s[j] <= s[i] + lesions_bp
    left  <- findInterval(s - lesions_bp - 0.5, s) + 1L
    right <- findInterval(s + lesions_bp,        s)
    cs    <- c(0, cumsum(m))
    enr_s <- cs[right + 1L] - cs[left]
    ws_s  <- right - left + 1L
    span_s <- s[right] - s[left]
    enr  <- numeric(n); ws <- integer(n); sp <- numeric(n)
    enr[o] <- enr_s; ws[o] <- ws_s; sp[o] <- span_s
    list(enrichment = enr, window_size = ws, span = sp)
  }

  # Vectorise per grouping unit. We process each gene/chr group independently;
  # within a group the window calculation is O(N log N).
  idx_by_grp <- split(seq_len(nrow(mutationAnnotatedSortedLocal)),
                      mutationAnnotatedSortedLocal[[grouping_column]])
  enrichment     <- numeric(nrow(mutationAnnotatedSortedLocal))
  window_size    <- integer(nrow(mutationAnnotatedSortedLocal))
  basepair_count <- numeric(nrow(mutationAnnotatedSortedLocal))
  for (idx in idx_by_grp) {
    if (length(idx) == 0L) next
    st  <- mutationAnnotatedSortedLocal$START[idx]
    mu  <- mutationAnnotatedSortedLocal$MUTATIONS[idx]
    stats_grp <- .bp_window_stats(st, mu, lesions_bp)
    enrichment[idx]     <- stats_grp$enrichment
    window_size[idx]    <- stats_grp$window_size
    basepair_count[idx] <- stats_grp$span
  }
  mutationAnnotatedSortedLocal$ENRICHMENT     <- enrichment
  mutationAnnotatedSortedLocal$WINDOW_SIZE    <- window_size
  mutationAnnotatedSortedLocal$BASEPAIR_COUNT <- basepair_count

  mutationAnnotatedSortedLocal$ENRICHMENT[is.na(mutationAnnotatedSortedLocal$ENRICHMENT)] <- 0

  # H0: ENRICHMENT ~ Binomial(WINDOW_SIZE, p0)
  # WINDOW_SIZE is the row-specific count of probes in the bp window (replaces
  # the legacy constant size = sliding_window_size). p0 is the empirical
  # background mutation rate of the grouping unit.
  # p-value = P(X >= ENRICHMENT) = pbinom(ENRICHMENT - 1, size, p0, lower.tail = FALSE)
  p0 <- mutationAnnotatedSortedLocal$MUTATIONS_COUNT / mutationAnnotatedSortedLocal$PROBES_COUNT

  lesionpValue <- stats::pbinom(
    mutationAnnotatedSortedLocal$ENRICHMENT - 1L,
    size       = as.integer(mutationAnnotatedSortedLocal$WINDOW_SIZE),
    prob       = p0,
    lower.tail = FALSE
  )

  lesionpValue[is.nan(lesionpValue)] <- 1
  lesionpValue[is.na(lesionpValue)]  <- 1

  tt <- data.frame(mutationAnnotatedSortedLocal, lesionpValue)

  # Bonferroni weighted by log10(BASEPAIR_COUNT + 9). The "+9" floor keeps
  # log10 >= 1 even when the window-span collapses to 0 (singleton-window
  # rows, e.g. LESIONS_BP=0 or isolated probes with no neighbours within the
  # bp threshold), preserving the spirit of the legacy weighting where the
  # divisor never reached zero on a non-trivial window.
  # The same correction as sem_lesions_get_bulk(): the number of probes times
  # log10 of the nominal window, 2 * LESIONS_BP. It used to be log10 of the span
  # the probes of each window actually cover, which is a lighter correction on a
  # sparse array and made this path call a few lesions the bulk path did not.
  # Which of the two the method should use is an open question; until it is
  # settled the two engines apply the same one.
  lesionWeighted <- (tt$lesionpValue) < (as.numeric(ssEnv$bonferroni_threshold) / (length(tt$PROBES_COUNT) * log10(2 * lesions_bp)))
  rm(tt)

  lesionWeighted <- data.frame(as.data.frame(mutationAnnotatedSortedLocal), "LESIONS" = lesionWeighted)

  lesionWeighted <- anno_sort_by_chr_and_start(lesionWeighted)
  lesionWeighted <- subset(lesionWeighted, lesionWeighted$LESIONS == TRUE)[, c("CHR", "START", "END")]

  core_log_event("DEBUG: ", format(Sys.time(), "%a %b %d %X %Y"), " Got lesions for sample !")
  return(lesionWeighted)

}

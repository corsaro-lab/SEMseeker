#' sem_lesions_get_bulk - multi-sample bp-window LESIONS computation
#'
#' Versione **multi-sample (bulk)** di sem_lesions_get, usata dal path
#' `sem_analyze_population_bulk()` per processare N sample in parallelo via
#' polars+R-side matrix ops. La versione single-sample (per il path legacy
#' `sem_analyze_population()` per-sample loop) e' in `R/lesions_get.R`. Pattern:
#' `<funzione_single_sample>` + `<funzione_single_sample>_bulk` come
#' altri sites del package (io_signal_save, sem_deltaX_get).
#'
#' Usa una **physical window in bp** centrata su ogni position p, con radius
#' `LESIONS_BP` letto da ssEnv (default 5000 bp).
#' Per ogni probe i su chr c, considera le sonde j con
#'   |START[j] - START[i]| <= LESIONS_BP   AND   CHR[j] == CHR[i]
#' Per ogni sample s, ENRICHMENT[i,s] = sum_{j in window} MUTATIONS[j,s];
#' size del binomial test = WINDOW_SIZE[i] (variabile per riga).
#'
#' Conseguenza biologica:
#' - Regioni dense (CpG island, ~50bp/sonda): la window contiene molte sonde
#'   (es. 100-200 in 5 kbp) -> binomial test sensibile a cluster locali
#' - Regioni sparse intergeniche: la window contiene poche sonde (es. 2-5 in
#'   5 kbp) -> binomial test richiede ENRICHMENT proporzionalmente piu' alto
#'   per essere significativo
#' - **Il concetto "aggregati mono-direzionali localizzati" diventa rigoroso**:
#'   il radius e' fisico, non logico.
#'
#' @param mut_df data.frame con colonne CHR, START, END + colonne sample
#'   contenenti i flag MUTATIONS (0/1). Devono essere ordinati per (CHR, START).
#' @param sample_cols vector di character: nomi delle colonne sample in mut_df.
#' @param LESIONS_BP integer: radius della finestra in bp. Default NULL -> letto
#'   da ssEnv$LESIONS_BP. Fallback finale 5000L se ssEnv vuoto.
#' @param bonf_threshold numeric: soglia Bonferroni base (default 0.05).
#'
#' @return matrice integer di 0/1, n_rows x length(sample_cols), con LESIONS
#'   per ogni (probe, sample).
#'
#' @keywords internal
#' @noRd
sem_lesions_get_bulk <- function(mut_df, sample_cols,
                             LESIONS_BP = NULL,
                             bonf_threshold = 0.05) {

  if (!all(c("CHR", "START") %in% colnames(mut_df))) {
    stop("[sem_lesions_get_bulk] mut_df must contain CHR + START columns")
  }
  if (length(sample_cols) == 0L) return(matrix(0L, nrow = nrow(mut_df), ncol = 0L))

  if (is.null(LESIONS_BP)) {
    ssEnv <- tryCatch(core_get_session_info(), error = function(e) NULL)
    LESIONS_BP <- if (!is.null(ssEnv) && !is.null(ssEnv$LESIONS_BP))
      as.integer(ssEnv$LESIONS_BP) else 5000L
  }
  window_bp <- as.numeric(LESIONS_BP)
  if (is.na(window_bp) || window_bp < 0)
    stop("[sem_lesions_get_bulk] LESIONS_BP must be >= 0 (bp), got: ", LESIONS_BP)

  # ----- Step 1: precompute window boundaries (left_idx, right_idx) per row -----
  # Per ogni riga i, trovare il range [j_left, j_right] tale che:
  #   mut_df$CHR[j] == mut_df$CHR[i]
  #   mut_df$START[i] - window_bp <= mut_df$START[j] <= mut_df$START[i] + window_bp
  # Implementazione: two-pointer per chr-block.
  chr <- as.character(mut_df$CHR)
  pos <- as.integer(mut_df$START)
  n   <- nrow(mut_df)

  left_idx  <- integer(n)
  right_idx <- integer(n)

  # Process by chromosome block (run-length)
  rle_chr <- rle(chr)
  ends   <- cumsum(rle_chr$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)

  for (blk in seq_along(rle_chr$values)) {
    bs <- starts[blk]; be <- ends[blk]
    block_pos <- pos[bs:be]
    # Two-pointer sweep within sorted block
    j_left <- 1L; j_right <- 1L
    n_block <- length(block_pos)
    for (i in seq_len(n_block)) {
      lo <- block_pos[i] - window_bp
      hi <- block_pos[i] + window_bp
      while (j_left <= n_block && block_pos[j_left] < lo) j_left <- j_left + 1L
      while (j_right <= n_block && block_pos[j_right] <= hi) j_right <- j_right + 1L
      left_idx[bs + i - 1L]  <- bs + j_left - 1L
      right_idx[bs + i - 1L] <- bs + j_right - 2L  # j_right pointed PAST last in-range
    }
  }

  # window_size_per_row = numero sonde nella finestra (varia per riga)
  window_size <- right_idx - left_idx + 1L

  # Bonferroni-weighted threshold: la finestra fisica ha larghezza definita
  # dal parametro semseeker(LESIONS_BP), quindi il peso log10(2*W) e'
  # **configurabile via parametro ma costante per la run**. Funzione monotona
  # di LESIONS_BP: aumentando il radius -> Bonferroni piu' permissivo (cluster
  # piu' larghi sono ammissibili). Riducendo -> piu' stringente (cluster solo
  # compatti).
  log_bp_window <- log10(2 * window_bp)
  n_probes <- n
  bonf_per_row <- bonf_threshold / (n_probes * log_bp_window)

  # ----- Step 2: ENRICHMENT + binomial test per sample (column-wise) -----
  les_mat <- matrix(0L, nrow = n, ncol = length(sample_cols))
  colnames(les_mat) <- sample_cols

  # Per ogni sample column: rolling sum via cumsum vectorized
  for (s_idx in seq_along(sample_cols)) {
    s <- sample_cols[s_idx]
    mut_s <- as.integer(mut_df[[s]])

    # ENRICHMENT_i = sum(mut_s[left_idx[i] : right_idx[i]])
    # Vettorizzato via cumsum: enrichment_i = cum[right_idx[i]] - cum[left_idx[i] - 1]
    cum <- cumsum(mut_s)
    cum_zero <- c(0L, cum)  # cum_zero[k+1] == cum[k]
    enrichment <- cum_zero[right_idx + 1L] - cum_zero[left_idx]

    # p0 per-CHR: rate di outlier in ciascun cromosoma
    p0_per_chr_block <- numeric(n)
    for (blk in seq_along(rle_chr$values)) {
      bs <- starts[blk]; be <- ends[blk]
      m_block <- sum(mut_s[bs:be])
      n_block <- be - bs + 1L
      p0_block <- if (m_block <= 0 || m_block >= n_block) NA_real_ else m_block / n_block
      p0_per_chr_block[bs:be] <- p0_block
    }

    # Binomial test: P(X >= ENRICHMENT | size = window_size_per_row, prob = p0_per_chr)
    valid <- !is.na(p0_per_chr_block) & p0_per_chr_block > 0 & p0_per_chr_block < 1 &
             window_size > 0L
    pvals <- rep(1, n)
    if (any(valid)) {
      pvals[valid] <- stats::pbinom(enrichment[valid] - 1L,
                                    size = window_size[valid],
                                    prob = p0_per_chr_block[valid],
                                    lower.tail = FALSE)
      pvals[is.na(pvals) | is.nan(pvals)] <- 1
    }
    les_mat[, s_idx] <- as.integer(pvals < bonf_per_row)
  }

  les_mat
}

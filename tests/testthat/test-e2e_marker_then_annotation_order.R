# Meta-test: genomic annotation must run AFTER derived-marker computation.
#
# Invariant (richiesto esplicitamente): le annotazioni genomiche (Body, TSS,
# Island shores/shelves, ...) prodotte da anno_annotate_position_pivots() devono
# essere SEMPRE successive al calcolo dei marker derivati DELTAQ/DELTARQ/
# DELTAP/DELTARP fatto da sem_deltaX_get(). Inoltre, dove i POSITION pivot di base
# vengono materializzati con anno_create_position_pivots(), quella chiamata deve
# precedere sem_deltaX_get() (sem_deltaX_get legge i pivot DELTAS/DELTAR base).
#
# anno_annotate_position_pivots() NON ricalcola i marker: legge il POSITION/WHOLE
# pivot di ciascun marker (inclusi i derivati) e lo proietta su AREA/SUBAREA.
# Se qualcuno invertisse l'ordine, le annotazioni dei marker derivati
# leggerebbero pivot inesistenti o stantii. Questo test blocca l'inversione.
#
# I corpi delle funzioni sono letti dal NAMESPACE installato (getFromNamespace +
# body()), NON dai file R/*.R: sotto R CMD check il pacchetto e' installato e la
# cartella sorgente R/ non e' presente, quindi leggere da disco fallirebbe.

# Walk del parse-tree del corpo di una funzione: restituisce, per ciascun
# nome-di-chiamata richiesto, l'indice di prima comparsa in ordine sorgente
# (NA se assente). L'indice e' un contatore globale incrementato a ogni nodo,
# quindi rispetta l'ordinamento testuale degli statement.
.first_call_orders <- function(fun_body, call_names) {
  counter <- 0L
  first <- stats::setNames(rep(NA_integer_, length(call_names)), call_names)
  walk <- function(node) {
    counter <<- counter + 1L
    if (is.call(node)) {
      head <- node[[1L]]
      if (is.symbol(head)) {
        nm <- as.character(head)
        if (nm %in% call_names && is.na(first[[nm]])) {
          first[[nm]] <<- counter
        }
      }
      for (i in seq_along(node)) {
        el <- node[[i]]
        if (!missing(el) && (is.call(el) || is.pairlist(el))) walk(el)
      }
    } else if (is.pairlist(node)) {
      for (i in seq_along(node)) {
        el <- node[[i]]
        if (!missing(el) && (is.call(el) || is.pairlist(el))) walk(el)
      }
    }
    invisible(NULL)
  }
  walk(fun_body)
  first
}

test_that("sem_deltaX_get() precedes anno_annotate_position_pivots() in every caller", {
  callers <- c("sem_core", "core_recover", "association_analysis")

  for (fun_name in callers) {
    fn <- tryCatch(getFromNamespace(fun_name, "SEMseeker"),
                   error = function(e) NULL)
    expect_false(is.null(fn),
      info = paste("could not find function", fun_name, "in SEMseeker namespace"))
    if (is.null(fn)) next

    ord <- .first_call_orders(
      body(fn),
      c("anno_create_position_pivots", "sem_deltaX_get", "anno_annotate_position_pivots")
    )

    # sem_deltaX_get() and anno_annotate_position_pivots() must BOTH appear and be ordered.
    expect_false(is.na(ord[["sem_deltaX_get"]]),
      info = paste(fun_name, "must call sem_deltaX_get()"))
    expect_false(is.na(ord[["anno_annotate_position_pivots"]]),
      info = paste(fun_name, "must call anno_annotate_position_pivots()"))
    expect_lt(ord[["sem_deltaX_get"]], ord[["anno_annotate_position_pivots"]])

    # Where anno_create_position_pivots() is present, it must precede sem_deltaX_get()
    # (sem_deltaX_get reads the base DELTAS/DELTAR POSITION pivots it writes).
    if (!is.na(ord[["anno_create_position_pivots"]])) {
      expect_lt(ord[["anno_create_position_pivots"]], ord[["sem_deltaX_get"]])
    }
  }
})

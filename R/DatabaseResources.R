.spamtp_db_cache <- new.env(parent = emptyenv())

.spamtp_db_legacy_names <- c(
  chem_props = "chem_props",
  source_df = "source_df",
  analyte = "analyte",
  analytehaspathway = "analytehaspathway",
  pathway = "pathway",
  ramp_db_metadata = "ramp_db_metadata",
  ramp_hmdb = "RAMP_hmdb",
  ramp_kegg = "RAMP_kegg",
  ramp_reactome = "RAMP_Reactome",
  ramp_wikipathway = "RAMP_wikipathway",
  hmdb_db = "HMDB_db",
  chebi_db = "Chebi_db",
  lipidmaps_db = "Lipidmaps_db",
  gnps_db = "GNPS_db",
  filtered_fmp10 = "filtered_fmp10"
)

# Structure features are intentionally separate from the pruned chemical table.
.spamtp_db_legacy_names <- c(
  .spamtp_db_legacy_names,
  smiles_features = "smiles_features"
)

.spamtp_db_normalise_resources <- function(resources) {
  resources <- tolower(trimws(as.character(resources)))
  resources <- unique(resources[nzchar(resources)])
  unknown <- setdiff(resources, names(.spamtp_db_legacy_names))
  if (length(unknown)) {
    stop(
      "Unknown SpaMTP database resource(s): ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  resources
}

.spamtp_db_cache_key <- function(resource, version, source, local_dir, verify) {
  local_key <- if (is.null(local_dir)) "<default>" else {
    normalizePath(local_dir, mustWork = FALSE)
  }
  paste(
    resource,
    version %||% "latest",
    source,
    local_key,
    verify,
    sep = "\r"
  )
}

.spamtp_db_label <- function(value, fallback = "user-supplied database") {
  metadata <- attr(value, "spamtp_database", exact = TRUE)
  if (!is.list(metadata)) return(fallback)
  version <- metadata$version %||% "unknown"
  paste0(metadata$source %||% "SpaMTP", " RaMP ", version, " ", metadata$resource)
}

.spamtp_db_resource <- function(resource,
                                version = "latest",
                                source = c("auto", "spamtpdb"),
                                local_dir = NULL,
                                hub = NULL,
                                offline = FALSE,
                                refresh = FALSE,
                                verify = TRUE) {
  resource <- .spamtp_db_normalise_resources(resource)
  if (length(resource) != 1L) {
    stop("resource must identify exactly one database resource.", call. = FALSE)
  }
  source <- match.arg(source)
  if (!is.logical(verify) || length(verify) != 1L || is.na(verify)) {
    stop("verify must be TRUE or FALSE.", call. = FALSE)
  }
  if (is.null(local_dir)) {
    configured <- getOption("SpaMTPdb.resource_dir", "")
    if (!nzchar(configured)) configured <- Sys.getenv("SPAMTPDB_RESOURCE_DIR", "")
    if (nzchar(configured)) local_dir <- configured
  }
  key <- .spamtp_db_cache_key(resource, version, source, local_dir, verify)
  if (!isTRUE(refresh) && exists(key, envir = .spamtp_db_cache, inherits = FALSE)) {
    return(get(key, envir = .spamtp_db_cache, inherits = FALSE))
  }

  if (!requireNamespace("SpaMTPdb", quietly = TRUE)) {
    stop(
      "Loading versioned annotation resources requires the SpaMTPdb package. ",
      "Install SpaMTPdb or supply a named custom database bundle.",
      call. = FALSE
    )
  }

  value <- SpaMTPdb::spaMTPdbResource(
    resource = resource,
    version = version,
    local_dir = local_dir,
    hub = hub,
    offline = offline,
    verify = verify
  )
  resource_metadata <- SpaMTPdb::spaMTPdbResource(
    resource = resource,
    version = version,
    metadata = TRUE
  )

  attr(value, "spamtp_database") <- list(
    resource = resource,
    version = as.character(resource_metadata$version[[1L]]),
    source = "spamtpdb",
    local_dir = local_dir,
    offline = offline,
    verify_local = verify
  )
  assign(key, value, envir = .spamtp_db_cache)
  value
}

.spamtp_db_bundle <- function(resources,
                              database = NULL,
                              version = "latest",
                              source = c("auto", "spamtpdb"),
                              local_dir = NULL,
                              hub = NULL,
                              offline = FALSE,
                              refresh = FALSE,
                              verify = TRUE) {
  resources <- .spamtp_db_normalise_resources(resources)
  source <- match.arg(source)
  if (!is.null(database)) {
    if (!is.list(database) || is.data.frame(database) || is.null(names(database))) {
      stop("database must be a named list of SpaMTP database resources.", call. = FALSE)
    }
    missing <- setdiff(resources, tolower(names(database)))
    if (length(missing)) {
      stop(
        "database is missing resource(s): ", paste(missing, collapse = ", "),
        ".", call. = FALSE
      )
    }
    names(database) <- tolower(names(database))
    return(database[resources])
  }

  values <- lapply(resources, function(resource) {
    .spamtp_db_resource(
      resource = resource,
      version = version,
      source = source,
      local_dir = local_dir,
      hub = hub,
      offline = offline,
      refresh = refresh,
      verify = verify
    )
  })
  stats::setNames(values, resources)
}

#' Load versioned SpaMTP annotation resources
#'
#' Loads a coherent group of annotation resources from [SpaMTPdb]. Retrieved
#' resources are cached for the current R session. A named custom bundle can be
#' supplied for offline, testing, or user-curated workflows.
#'
#' @param resources Character vector of resource names. Use
#'   [spaMTPDatabaseInfo()] to list valid names.
#' @param version SpaMTPdb/RaMP resource version, or `"latest"`.
#' @param source Database source. `"auto"` and `"spamtpdb"` resolve versioned
#'   resources through SpaMTPdb.
#' @param database Optional named list containing the requested resources. When
#'   supplied, no Hub lookup is performed.
#' @param local_dir Optional directory containing staged SpaMTPdb `.rds` files.
#' @param hub Optional pre-created `AnnotationHub` object passed to SpaMTPdb.
#' @param offline If `TRUE`, do not query AnnotationHub.
#' @param refresh If `TRUE`, bypass SpaMTP's in-session resource cache.
#' @param verify Verify local resource files against the SpaMTPdb registry.
#'   Defaults to `TRUE`. Set to `FALSE` only for intentional development
#'   fixtures; use `database` for a named custom bundle instead of presenting
#'   modified resources as an official snapshot.
#'
#' @return A named list containing the requested resources.
#' @export
#'
#' @examples
#' utils::str(formals(loadSpaMTPDatabase))
#' example_database <- list(
#'   ramp_db_metadata = list(ramp_version = "example")
#' )
#' database <- loadSpaMTPDatabase(
#'   "ramp_db_metadata",
#'   database = example_database
#' )
#' names(database)
loadSpaMTPDatabase <- function(
    resources = c(
      "chem_props", "source_df", "analyte", "analytehaspathway", "pathway"
    ),
    version = "latest",
    source = c("auto", "spamtpdb"),
    database = NULL,
    local_dir = NULL,
    hub = NULL,
    offline = FALSE,
    refresh = FALSE,
    verify = TRUE) {
  .spamtp_db_bundle(
    resources = resources,
    database = database,
    version = version,
    source = match.arg(source),
    local_dir = local_dir,
    hub = hub,
    offline = offline,
    refresh = refresh,
    verify = verify
  )
}

#' Inspect SpaMTP database resources
#'
#' @param version Optional SpaMTPdb resource version. `NULL` lists every
#'   available external version.
#'
#' @return A data frame describing resources in SpaMTPdb.
#' @export
#'
#' @examples
#' utils::str(formals(spaMTPDatabaseInfo))
#' spaMTPDatabaseInfo()
spaMTPDatabaseInfo <- function(version = NULL) {
  if (!requireNamespace("SpaMTPdb", quietly = TRUE)) {
    stop(
      "spaMTPDatabaseInfo() requires the SpaMTPdb package.",
      call. = FALSE
    )
  }
  SpaMTPdb::spaMTPdbResources(version = version)
}

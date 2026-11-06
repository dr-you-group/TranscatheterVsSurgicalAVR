# Copyright 2020 Observational Health Data Sciences and Informatics
#
# This file is part of TavrSavr
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# The TAVR cohorts analysed by default: TAVR_ob0 (874) and TAVR_less80_ob0 (876).
udiDiTavrCohortIds <- c(874, 876)

# The tables written by runUdiDiAnalysis(), in outputFolder/udiDi, one set per
# TAVR cohort: tab_udi_products_874.csv, tab_udi_products_876.csv, ...
udiDiResultTables <- c("tab_udi_products",
                       "tab_udi_reintervention_by_manufacturer",
                       "tab_udi_reintervention_by_generation",
                       "tab_udi_reintervention_by_size",
                       "tab_udi_reintervention_by_size_pooled")
udiDiFileName <- function(table, tavrCohortId) {
  paste0(table, "_", tavrCohortId, ".csv")
}
udiDiResultFiles <- as.vector(outer(udiDiResultTables, udiDiTavrCohortIds, udiDiFileName))

#' Run the UDI-DI analysis
#'
#' @description
#' Resolves the index valve of each TAVR cohort (by default 874 and 876) to a product - model
#' generation and valve size - through its UDI-DI, and describes aortic valve
#' reintervention (cohort 797) by manufacturer, generation and size. The SAVR
#' comparator (782) has no UDI-DI and is not analysed.
#'
#' @details
#' Read-only: a few SELECT queries, aggregated in R. Needs only the cohort
#' table, not the CohortMethod output. Writes five aggregate CSVs per TAVR
#' cohort to \code{outputFolder/udiDi}, the cohort ID appended to the file
#' name; see \code{udiDiResultFiles}. Only the files of the default cohorts
#' (\code{udiDiTavrCohortIds}) are picked up by \code{exportResults}. Hospital item codes
#' (device_source_value) and raw unique_device_id values are never written.
#'
#' The outcome tables are DESCRIPTIVE: no propensity score, no adjusted
#' estimate. Reintervention is rare, so the generation and size cells are a
#' feasibility check rather than a result.
#'
#' The outcome cohort is First-limited (one row per person), so a patient whose
#' row falls on or before the index date can never have a post-index event.
#' Those patients are excluded from the at-risk set rather than counted as
#' denominator.
#'
#' @param connectionDetails    An object of type \code{connectionDetails} as created using the
#'                             \code{\link[DatabaseConnector]{createConnectionDetails}} function.
#' @param cdmDatabaseSchema    Schema name where your patient-level data in OMOP CDM format resides.
#' @param cohortDatabaseSchema Schema name where the study cohort table resides.
#' @param cohortTable          The name of the study cohort table.
#' @param oracleTempSchema     Should be used in Oracle to specify a schema where the user has write
#'                             priviliges for storing temporary tables.
#' @param outputFolder         Name of local folder where the results are written, in a sub-folder
#'                             'udiDi'.
#' @param databaseId           A short string for identifying the database.
#' @param deviceExposureTable  Fully qualified name of the table holding the device exposures with
#'                             \code{unique_device_id}. NULL uses the CDM \code{device_exposure}.
#' @param tavrCohortIds        The TAVR cohorts, each analysed separately.
#' @param outcomeId            The reintervention outcome cohort.
#' @param valveConceptId       The standard device concept every TAVR valve maps to.
#' @param riskWindowStart      Start of the risk window in days after index, as in analysis 112.
#' @param riskWindowEnd        End of the risk window in days after index, as in analysis 112.
#'
#' @return
#' Invisibly, the path of the folder the results were written to.
#'
#' @export
runUdiDiAnalysis <- function(connectionDetails,
                             cdmDatabaseSchema,
                             cohortDatabaseSchema = cdmDatabaseSchema,
                             cohortTable = "cohort",
                             oracleTempSchema = NULL,
                             outputFolder,
                             databaseId,
                             deviceExposureTable = NULL,
                             tavrCohortIds = udiDiTavrCohortIds,
                             outcomeId = 797,
                             valveConceptId = 3661561,
                             riskWindowStart = 1,
                             riskWindowEnd = 3285) {
  if (is.null(deviceExposureTable)) {
    deviceExposureTable <- paste0(cdmDatabaseSchema, ".device_exposure")
  }
  udiFolder <- file.path(outputFolder, "udiDi")
  if (!file.exists(udiFolder)) {
    dir.create(udiFolder, recursive = TRUE)
  }
  udiReference <- loadUdiReference()

  # Queries --------------------------------------------------------------------
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))
  runQuery <- function(sqlFilename, tavrCohortId, ...) {
    sql <- SqlRender::loadRenderTranslateSql(sqlFilename = sqlFilename,
                                             packageName = "TavrSavr",
                                             dbms = attr(connection, "dbms"),
                                             oracleTempSchema = oracleTempSchema,
                                             cohort_database_schema = cohortDatabaseSchema,
                                             cohort_table = cohortTable,
                                             cohort_ids = tavrCohortId,
                                             ...)
    result <- DatabaseConnector::querySql(connection, sql)
    names(result) <- tolower(names(result))
    result
  }
  # One set of tables per TAVR cohort ---------------------------------------------
  for (tavrCohortId in tavrCohortIds) {
    raw <- runQuery("UdiIndexDevice.sql",
                    tavrCohortId = tavrCohortId,
                    cdm_database_schema = cdmDatabaseSchema,
                    device_exposure_table = deviceExposureTable,
                    valve_concept_id = valveConceptId)
    outcomeRows <- runQuery("UdiOutcomeRows.sql",
                            tavrCohortId = tavrCohortId,
                            outcome_cohort_ids = outcomeId)

    # Resolve each patient to one product ----------------------------------------
    raw$matchedDi <- matchUdiDi(raw$unique_device_id, udiReference$udiDi)
    person <- raw %>%
      group_by(.data$subject_id) %>%
      summarise(cohortStartDate = min(.data$cohort_start_date),
                opEndDate = max(.data$observation_period_end_date),
                nDistinctDi = dplyr::n_distinct(.data$matchedDi[!is.na(.data$matchedDi)]),
                matchedDi = dplyr::first(.data$matchedDi[!is.na(.data$matchedDi)]),
                .groups = "drop")
    # Patients with no matching UDI-DI, or with two different products on the
    # index date, cannot be assigned to one product.
    mapped <- person %>%
      filter(.data$nDistinctDi == 1) %>%
      inner_join(udiReference, by = c("matchedDi" = "udiDi"))
    ParallelLogger::logInfo("UDI-DI, TAVR cohort ", tavrCohortId, ": ", nrow(person),
                            " patients with an index valve exposure, ",
                            nrow(mapped), " resolved to a single product")

    # Patients per product, every reference product listed ------------------------
    products <- udiReference %>%
      left_join(mapped %>% count(matchedDi = .data$matchedDi, name = "nPatients"),
                by = c("udiDi" = "matchedDi")) %>%
      mutate(nPatients = dplyr::coalesce(.data$nPatients, 0L)) %>%
      arrange(.data$manufacturer, .data$generation, .data$sizeMm) %>%
      select("manufacturer", "generation", "sizeMm", "modelName", "udiDi", "nPatients")
    writeUdiTable(products, udiFolder, udiDiFileName("tab_udi_products", tavrCohortId), databaseId)

    # Reintervention ---------------------------------------------------------------
    atRisk <- buildUdiAtRisk(outcomeRows = outcomeRows,
                             mapped = mapped,
                             riskWindowStart = riskWindowStart,
                             riskWindowEnd = riskWindowEnd)
    groupings <- list(by_manufacturer = "manufacturer",
                      by_generation = c("manufacturer", "generation"),
                      by_size = c("manufacturer", "sizeMm"),
                      by_size_pooled = "sizeMm")
    for (name in names(groupings)) {
      table <- summariseUdiOutcome(atRisk, groupings[[name]])
      fileName <- udiDiFileName(paste0("tab_udi_reintervention_", name), tavrCohortId)
      writeUdiTable(table, udiFolder, fileName, databaseId)
    }
  }

  ParallelLogger::logInfo("UDI-DI results written to ", udiFolder)
  invisible(udiFolder)
}


# The UDI-DI lookup. A UDI-DI is a GS1 GTIN and carries no model or size
# information on its face, so this table is what the analysis needs. Rows are
# in release order; generation is made a factor in that order so that the
# tables sort Evolut R -> Pro -> Pro+ rather than alphabetically.
loadUdiReference <- function() {
  pathToCsv <- system.file("settings", "UdiDiReference.csv", package = "TavrSavr")
  udiReference <- read.csv(pathToCsv, stringsAsFactors = FALSE, colClasses = c(udiDi = "character"))
  if (any(duplicated(udiReference$udiDi))) {
    stop("UdiDiReference.csv contains duplicated UDI-DIs")
  }
  udiReference$generation <- factor(udiReference$generation, levels = unique(udiReference$generation))
  udiReference
}


writeUdiTable <- function(df, udiFolder, fileName, databaseId) {
  df <- cbind(databaseId = databaseId, as.data.frame(df), stringsAsFactors = FALSE)
  write.csv(df, file.path(udiFolder, fileName), row.names = FALSE, na = "")
}


# Match raw unique_device_id values against the reference DIs.
#
# The field may hold the bare GTIN, the GTIN with its GS1 application
# identifier - "(01)00763000017699" - or a full UDI with production data
# appended. Stripping to digits and searching for the reference DI inside covers
# all three. A value matching more than one reference DI is returned as NA.
matchUdiDi <- function(raw, reference) {
  digits <- gsub("[^0-9]", "", ifelse(is.na(raw), "", as.character(raw)))
  vapply(digits, function(d) {
    hits <- reference[vapply(reference, function(di) grepl(di, d, fixed = TRUE), logical(1))]
    if (nchar(d) > 0 && length(hits) == 1) hits else NA_character_
  }, character(1), USE.NAMES = FALSE)
}


# Time at risk opens riskWindowStart days after index and closes at the end of
# observation or riskWindowEnd. Patients whose (only) outcome row is on or
# before index are excluded.
buildUdiAtRisk <- function(outcomeRows, mapped, riskWindowStart, riskWindowEnd) {
  perPerson <- outcomeRows %>%
    inner_join(mapped %>% select("subject_id", "cohortStartDate"), by = "subject_id") %>%
    mutate(daysFromIndex = as.numeric(.data$outcome_date - .data$cohortStartDate)) %>%
    group_by(.data$subject_id) %>%
    summarise(priorOrSameDay = any(.data$daysFromIndex <= 0),
              daysToEvent = suppressWarnings(min(.data$daysFromIndex[.data$daysFromIndex > 0])),
              .groups = "drop")

  atRisk <- mapped %>%
    left_join(perPerson, by = "subject_id") %>%
    filter(!dplyr::coalesce(.data$priorOrSameDay, FALSE)) %>%
    mutate(daysToEnd = pmin(as.numeric(.data$opEndDate - .data$cohortStartDate), riskWindowEnd),
           event = is.finite(.data$daysToEvent) &
             .data$daysToEvent >= riskWindowStart &
             .data$daysToEvent <= .data$daysToEnd,
           timeAtRisk = ifelse(.data$event, .data$daysToEvent, .data$daysToEnd)) %>%
    filter(.data$timeAtRisk >= riskWindowStart)

  ParallelLogger::logInfo("UDI-DI: ", nrow(mapped) - nrow(atRisk), " of ", nrow(mapped),
                          " patients excluded from the at-risk set (outcome on or before index, ",
                          "or no time at risk); ", sum(atRisk$event), " reinterventions")
  atRisk
}


# Kaplan-Meier cumulative risk at `atDays`, as a percentage.
kmUdiRisk <- function(time, event, atDays) {
  if (sum(event) == 0) return(0)
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  round(100 * (1 - summary(fit, times = atDays, extend = TRUE)$surv[1]), 1)
}


# One descriptive row per group: counts, person-time, crude rate, and the
# Kaplan-Meier cumulative risk at one and five years.
summariseUdiOutcome <- function(atRisk, groupColumns) {
  # summarise() returns the groups sorted, generation in release order.
  atRisk %>%
    group_by(dplyr::across(dplyr::all_of(groupColumns))) %>%
    summarise(nPatients = dplyr::n(),
              nEvents = sum(.data$event),
              personYears = round(sum(.data$timeAtRisk) / 365.25, 1),
              medianFollowUpDays = stats::median(.data$timeAtRisk),
              km1yPct = kmUdiRisk(.data$timeAtRisk, .data$event, 365),
              km5yPct = kmUdiRisk(.data$timeAtRisk, .data$event, 1826),
              .groups = "drop") %>%
    mutate(ratePer100Py = round(100 * .data$nEvents / .data$personYears, 2)) %>%
    select(dplyr::all_of(groupColumns), "nPatients", "nEvents", "personYears", "ratePer100Py",
           "medianFollowUpDays", "km1yPct", "km5yPct")
}

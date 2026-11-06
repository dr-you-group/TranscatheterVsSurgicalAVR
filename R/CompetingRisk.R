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

# Non-standard evaluation used by survival::coxph(weights = fgwt):
utils::globalVariables(c("fgwt"))

#' Run the competing risk complementary analysis
#'
#' @description
#' Complementary analysis accounting for the competing risk of death for the
#' non-fatal outcomes of the TAVR versus SAVR comparative effectiveness study.
#'
#' @details
#' Severe aortic stenosis carries substantial mortality. All-cause death is a
#' competing risk for the non-fatal outcomes (re-intervention, stroke, atrial
#' fibrillation, myocardial infarction, new permanent pacemaker implantation,
#' heart failure, aortic valve hospitalisation): once a patient dies these
#' events can no longer occur. The complement of the Kaplan-Meier estimator
#' (1-KM) overestimates absolute risk in this setting regardless of whether
#' death is informative. This function therefore estimates absolute risks via
#' the cumulative incidence function (CIF) and reports two complementary
#' hazard models:
#'
#' \itemize{
#'   \item Cause-specific Cox PH: \code{csHR} = rate of the event among those
#'         still alive and event-free (etiologic question).
#'   \item Fine-Gray subdistribution hazard: \code{sHR} = summary measure that
#'         corresponds directly to the CIF (prognostic question). Its risk set
#'         retains individuals who have already died, so the \code{sHR} should
#'         not be interpreted mechanistically.
#' }
#'
#' Outcomes that are not subject to the competing risk of death are excluded
#' via \code{excludedOutcomeIds}. By default this drops 856 (all-cause
#' mortality, which IS the competing event) together with the two composite
#' endpoints that already contain death as a component: 854 (all-cause
#' mortality or disabling stroke) and 857 (all-cause mortality, inpatient
#' stroke or aortic valve hospitalisation). Those three are analysed by the
#' main CohortMethod pipeline (standard Cox PH + KM), where 1-KM is
#' appropriate.
#'
#' Each model is fitted in the same population as the corresponding main
#' result (the matched or PS-stratified population) and is stratified if and
#' only if that analysis sets \code{stratified = TRUE} in
#' \code{cmAnalysisList.json}. Note that the 1:1 matched analyses carry one
#' \code{stratumId} per matched pair but are fitted unstratified, so the
#' setting cannot be inferred from the population alone.
#'
#' The target/comparator pairs are read from \code{TcosOfInterest.csv} and
#' every pair is analysed for every requested analysis ID, so results do not
#' overwrite one another. Results are written per (target, comparator,
#' analysis) triplet and additionally pooled into one file per table. Every
#' level of the loop is wrapped in \code{tryCatch} so that a single failing
#' model, outcome, or pair cannot abort \code{\link{execute}}.
#'
#' Prerequisites: \code{runAnalyses} in \code{\link{execute}} must have
#' completed, and the \code{cmprsk} package must be installed (it is recorded
#' in \code{renv.lock}).
#'
#' All outputs are aggregated; no patient-level data is written. Counts are
#' censored below \code{minCellCount} by \code{enforceMinCellValue}, and the
#' exported incidence curves are truncated at the 90th percentile of follow-up
#' and reported on a fixed time grid so that no individual's exact event day
#' is disclosed. Four CSV files are written: \code{cr_estimates} (the csHR and
#' sHR with counts), \code{cr_cif} (the curves, with the number still at
#' risk), \code{cr_gray_test} and \code{cr_ph_test}.
#'
#' @param outputFolder         Name of local folder where the results were generated; make sure to
#'                             use forward slashes (/). The same \code{outputFolder} used in
#'                             \code{\link{execute}}.
#' @param databaseId           A short string for identifying the database (e.g. 'Synpuf').
#' @param analysisIds          The CohortMethod analysis IDs to analyse. Defaults to the four
#'                             estimation analyses (9-year and 5-year landmark, pathway-reduced
#'                             covariate PS, 1:1 matching and stratification). Use \code{NULL} to
#'                             run every analysis in \code{cmAnalysisList.json}; note that this
#'                             multiplies runtime by roughly five.
#' @param deathOutcomeId       Cohort ID of the competing event, all-cause mortality.
#' @param excludedOutcomeIds   Outcome IDs that are NOT subject to the competing risk of death and
#'                             are therefore skipped. See details.
#' @param targetComparatorPairs Optional \code{data.frame} with columns \code{targetId} and
#'                             \code{comparatorId} restricting the analysis to those pairs. Defaults
#'                             to all pairs in \code{TcosOfInterest.csv}.
#' @param cifTimeGridDays      Spacing, in days, of the time grid the exported cumulative
#'                             incidence curves are evaluated on. Writing a row per event time
#'                             instead would make cr_cif.csv hundreds of megabytes and would
#'                             disclose the exact event day of individual subjects.
#' @param minCellCount         The minimum number of subjects contributing to a count before it can
#'                             be included in the results.
#'
#' @return
#' Invisibly, the path of the folder the results were written to.
#'
#' @export
runCompetingRiskAnalysis <- function(outputFolder,
                                     databaseId,
                                     analysisIds = c(111, 112, 511, 512),
                                     deathOutcomeId = 856,
                                     excludedOutcomeIds = c(854, 856, 857),
                                     targetComparatorPairs = NULL,
                                     cifTimeGridDays = 30,
                                     minCellCount = 5) {
  if (!is_installed("cmprsk")) {
    stop("The 'cmprsk' package is required for the competing risk analysis. ",
         "It is recorded in renv.lock; run renv::restore(), or ",
         "install.packages('cmprsk').")
  }

  cmOutputFolder <- file.path(outputFolder, "cmOutput")
  referenceFile <- file.path(cmOutputFolder, "outcomeModelReference.rds")
  if (!file.exists(referenceFile)) {
    stop("CohortMethod output not found at ",
         referenceFile,
         ". Run execute() with runAnalyses = TRUE first.")
  }
  reference <- readRDS(referenceFile)

  crFolder <- file.path(outputFolder, "competingRisk")
  if (!file.exists(crFolder)) {
    dir.create(crFolder, recursive = TRUE)
  }

  if (is.null(analysisIds)) {
    analysisIds <- sort(unique(reference$analysisId))
  }
  tcos <- getCompetingRiskTcos(excludedOutcomeIds = excludedOutcomeIds,
                              targetComparatorPairs = targetComparatorPairs)
  getCohortName <- createCohortNameLookup()
  stratifiedByAnalysisId <- getStratifiedByAnalysisId()

  ParallelLogger::logInfo("Running competing risk analysis for ",
                          nrow(tcos),
                          " target-comparator pairs and ",
                          length(analysisIds),
                          " analyses")

  pooled <- list(estimates = list(), cif = list(), gray = list(), ph = list())

  for (i in 1:nrow(tcos)) {
    targetId <- tcos$targetId[i]
    comparatorId <- tcos$comparatorId[i]
    outcomeIds <- tcos$outcomeIds[[i]]
    if (length(outcomeIds) == 0) {
      next
    }
    for (analysisId in analysisIds) {
      label <- sprintf("t%d, c%d, a%d", targetId, comparatorId, analysisId)
      stratified <- stratifiedByAnalysisId[[as.character(analysisId)]]
      if (is.null(stratified)) {
        ParallelLogger::logWarn("Analysis ", analysisId, " is not in cmAnalysisList.json. Skipping.")
        next
      }
      result <- tryCatch({
        runCompetingRiskForTco(reference = reference,
                               cmOutputFolder = cmOutputFolder,
                               targetId = targetId,
                               comparatorId = comparatorId,
                               analysisId = analysisId,
                               outcomeIds = outcomeIds,
                               deathOutcomeId = deathOutcomeId,
                               stratified = stratified,
                               cifTimeGridDays = cifTimeGridDays,
                               getCohortName = getCohortName,
                               databaseId = databaseId)
      }, error = function(e) {
        ParallelLogger::logWarn("Competing risk analysis failed for ", label, ": ", e$message)
        NULL
      })
      if (is.null(result)) {
        next
      }
      exportCompetingRiskTables(result = result,
                                crFolder = crFolder,
                                fileSuffix = sprintf("t%d_c%d_a%d", targetId, comparatorId, analysisId),
                                minCellCount = minCellCount,
                                silent = TRUE)
      for (table in names(pooled)) {
        pooled[[table]] <- c(pooled[[table]], result[[table]])
      }
    }
  }

  ParallelLogger::logInfo("Exporting pooled competing risk results")
  exportCompetingRiskTables(result = pooled,
                            crFolder = crFolder,
                            fileSuffix = NULL,
                            minCellCount = minCellCount,
                            silent = FALSE)

  ParallelLogger::logInfo("Competing risk results written to ", crFolder)
  invisible(crFolder)
}


# Determine which target-comparator pairs and outcomes are subject to the
# competing risk of death. Mirrors the TcosOfInterest.csv parsing in
# createTcos(), but keeps only the outcomes of interest: negative and
# synthesised positive controls are not part of this complementary analysis.
getCompetingRiskTcos <- function(excludedOutcomeIds, targetComparatorPairs) {
  pathToCsv <- system.file("settings", "TcosOfInterest.csv", package = "TavrSavr")
  tcosOfInterest <- read.csv(pathToCsv, stringsAsFactors = FALSE)
  tcs <- unique(tcosOfInterest[, c("targetId", "comparatorId")])
  if (!is.null(targetComparatorPairs)) {
    tcs <- merge(tcs, unique(targetComparatorPairs[, c("targetId", "comparatorId")]))
    if (nrow(tcs) == 0) {
      stop("None of the requested target-comparator pairs occur in TcosOfInterest.csv")
    }
  }
  tcs$outcomeIds <- lapply(1:nrow(tcs), function(i) {
    outcomeIds <- as.character(tcosOfInterest$outcomeIds[tcosOfInterest$targetId == tcs$targetId[i] &
                                                           tcosOfInterest$comparatorId == tcs$comparatorId[i]])
    outcomeIds <- unique(as.numeric(strsplit(outcomeIds, split = ";")[[1]]))
    outcomeIds[!outcomeIds %in% excludedOutcomeIds]
  })
  return(tcs)
}


# Whether each analysis fits a stratified outcome model, keyed by analysis ID.
# This must come from cmAnalysisList.json rather than being inferred from the
# presence of stratumId: the 1:1 matched analyses carry one stratumId per
# matched pair but CohortMethod fits them with stratified = FALSE, so
# stratifying here would stop the competing risk estimates from being
# comparable to the main results.
getStratifiedByAnalysisId <- function() {
  cmAnalysisListFile <- system.file("settings", "cmAnalysisList.json", package = "TavrSavr")
  cmAnalysisList <- CohortMethod::loadCmAnalysisList(cmAnalysisListFile)
  stratified <- lapply(cmAnalysisList, function(cmAnalysis) {
    isTRUE(cmAnalysis$fitOutcomeModelArgs$stratified)
  })
  names(stratified) <- sapply(cmAnalysisList, function(cmAnalysis) as.character(cmAnalysis$analysisId))
  return(stratified)
}


# Readable cohort names for the result tables. The ATLAS names
# carry a '[26HIRAMDV] ' study prefix that only adds noise here.
createCohortNameLookup <- function() {
  pathToCsv <- system.file("settings", "CohortsToCreate.csv", package = "TavrSavr")
  cohortsToCreate <- read.csv(pathToCsv, stringsAsFactors = FALSE)
  function(cohortId) {
    name <- cohortsToCreate$atlasName[cohortsToCreate$cohortId == cohortId]
    if (length(name) == 0) {
      return(paste0("Cohort_", cohortId))
    }
    trimws(sub("^\\[[^]]*\\]", "", name[1]))
  }
}


# Run every competing risk model for one (target, comparator, analysis)
# triplet. Returns a list of lists of data frames, one entry per outcome and
# model, ready to be row-bound.
runCompetingRiskForTco <- function(reference,
                                   cmOutputFolder,
                                   targetId,
                                   comparatorId,
                                   analysisId,
                                   outcomeIds,
                                   deathOutcomeId,
                                   stratified,
                                   cifTimeGridDays,
                                   getCohortName,
                                   databaseId) {
  result <- list(estimates = list(), cif = list(), gray = list(), ph = list())

  # Shared across every table, so that pooled results stay identifiable:
  keys <- data.frame(databaseId = databaseId,
                     targetId = targetId,
                     comparatorId = comparatorId,
                     analysisId = analysisId,
                     stratified = stratified,
                     stringsAsFactors = FALSE)

  for (outcomeId in outcomeIds) {
    outcomeName <- getCohortName(outcomeId)
    ParallelLogger::logInfo(sprintf("  t%d, c%d, a%d, outcome %d: %s",
                                    targetId, comparatorId, analysisId, outcomeId, outcomeName))

    crData <- tryCatch({
      loadAndMergeCrPopulations(reference = reference,
                                cmOutputFolder = cmOutputFolder,
                                targetId = targetId,
                                comparatorId = comparatorId,
                                analysisId = analysisId,
                                outcomeId = outcomeId,
                                deathOutcomeId = deathOutcomeId)
    }, error = function(e) {
      ParallelLogger::logWarn("    Could not assemble competing risk data: ", e$message)
      NULL
    })
    if (is.null(crData) || nrow(crData) == 0) {
      next
    }

    nEvents <- sum(crData$status == 1)
    nDeaths <- sum(crData$status == 2)
    if (nEvents == 0) {
      ParallelLogger::logInfo("    No events of interest observed. Skipping.")
      next
    }
    ParallelLogger::logInfo(sprintf("    n = %d, events = %d, competing deaths = %d, censored = %d",
                                    nrow(crData), nEvents, nDeaths, sum(crData$status == 0)))
    if (nDeaths == 0) {
      ParallelLogger::logInfo("    No competing deaths observed; the subdistribution and ",
                              "cause-specific hazards coincide for this outcome.")
    }

    counts <- cbind(keys,
                    data.frame(outcomeId = outcomeId,
                               outcomeName = outcomeName,
                               nTarget = sum(crData$treatment == 1),
                               nComparator = sum(crData$treatment == 0),
                               eventsTarget = sum(crData$status == 1 & crData$treatment == 1),
                               eventsComparator = sum(crData$status == 1 & crData$treatment == 0),
                               competingEventsTarget = sum(crData$status == 2 & crData$treatment == 1),
                               competingEventsComparator = sum(crData$status == 2 & crData$treatment == 0),
                               stringsAsFactors = FALSE))

    models <- list(cause_specific_cox = fitCauseSpecificCox, fine_gray = fitFineGrayModel)
    for (model in names(models)) {
      fit <- tryCatch(models[[model]](crData, stratified), error = function(e) {
        ParallelLogger::logWarn("    The ", model, " model failed: ", e$message)
        NULL
      })
      if (is.null(fit)) {
        next
      }
      result$estimates[[length(result$estimates) + 1]] <-
        cbind(counts,
              data.frame(model = model,
                         hr = fit$hr,
                         ci95Lb = fit$ci95Lb,
                         ci95Ub = fit$ci95Ub,
                         logHr = fit$logHr,
                         seLogHr = fit$seLogHr,
                         stringsAsFactors = FALSE))
      result$ph[[length(result$ph) + 1]] <-
        cbind(keys,
              data.frame(outcomeId = outcomeId,
                         outcomeName = outcomeName,
                         model = model,
                         chiSq = fit$phChiSq,
                         df = fit$phDf,
                         pValue = fit$phP,
                         stringsAsFactors = FALSE))
    }

    cif <- tryCatch(computeAndTestCif(crData, cifTimeGridDays), error = function(e) {
      ParallelLogger::logWarn("    CIF computation failed: ", e$message)
      NULL
    })
    if (is.null(cif)) {
      next
    }
    result$cif[[length(result$cif) + 1]] <-
      cbind(keys,
            data.frame(outcomeId = outcomeId, stringsAsFactors = FALSE),
            cif$cifData)
    result$gray[[length(result$gray) + 1]] <-
      cbind(keys,
            data.frame(outcomeId = outcomeId,
                       outcomeName = outcomeName,
                       testStatistic = cif$grayTestStat,
                       df = cif$grayDf,
                       pValue = cif$grayP,
                       stringsAsFactors = FALSE))
  }
  return(result)
}


# Merge the event-of-interest population with the all-cause death population
# to build a composite competing risk status:
#   1 = event of interest, 2 = competing death, 0 = censored.
loadAndMergeCrPopulations <- function(reference,
                                      cmOutputFolder,
                                      targetId,
                                      comparatorId,
                                      analysisId,
                                      outcomeId,
                                      deathOutcomeId) {
  getPopFile <- function(id) {
    ref <- reference[reference$targetId == targetId &
                       reference$comparatorId == comparatorId &
                       reference$outcomeId == id &
                       reference$analysisId == analysisId, ]
    if (nrow(ref) == 0) {
      stop(sprintf("No CohortMethod reference row for outcome %d", id))
    }
    # Use the stratified/matched population when the analysis produced one, so
    # that the competing risk models are fitted in the same population as the
    # main results.
    popFile <- if (ref$strataFile[1] != "") ref$strataFile[1] else ref$studyPopFile[1]
    if (popFile == "") {
      stop(sprintf("No study population file for outcome %d", id))
    }
    return(file.path(cmOutputFolder, popFile))
  }

  eventPop <- readRDS(getPopFile(outcomeId))
  deathPop <- readRDS(getPopFile(deathOutcomeId))

  eventCols <- intersect(c("rowId", "treatment", "survivalTime", "outcomeCount", "stratumId"),
                         names(eventPop))
  merged <- merge(eventPop[, eventCols],
                  deathPop[, c("rowId", "survivalTime", "outcomeCount")],
                  by = "rowId",
                  suffixes = c("", "Death"),
                  all.x = TRUE)

  # Subjects present in the outcome population but not in the death population
  # (the two populations are filtered per outcome, so a handful can differ).
  # They contribute follow-up but no competing event.
  unmatched <- is.na(merged$survivalTimeDeath)
  if (any(unmatched)) {
    ParallelLogger::logInfo("    ",
                            sum(unmatched),
                            " subjects absent from the death population; treated as death-free")
    merged$survivalTimeDeath[unmatched] <- merged$survivalTime[unmatched]
    merged$outcomeCountDeath[unmatched] <- 0
  }

  # The event of interest wins a tie: CohortMethod censors the outcome
  # population at the outcome, so survivalTime is already the time to the
  # event when outcomeCount > 0.
  merged$status <- ifelse(merged$outcomeCount > 0,
                          1L,
                          ifelse(merged$outcomeCountDeath > 0, 2L, 0L))
  merged$time <- ifelse(merged$status == 2, merged$survivalTimeDeath, merged$survivalTime)

  merged$outcomeCount <- NULL
  merged$outcomeCountDeath <- NULL
  merged$survivalTime <- NULL
  merged$survivalTimeDeath <- NULL

  # Zero follow-up carries no information and breaks finegray().
  merged <- merged[!is.na(merged$time) & merged$time > 0, ]
  return(merged)
}


# Should this fit be stratified? Only when the corresponding CohortMethod
# analysis stratifies its outcome model AND the population actually carries
# more than one stratum.
hasCrStrata <- function(crData, stratified) {
  stratified &&
    "stratumId" %in% names(crData) &&
    length(unique(crData$stratumId)) > 1
}


# Extract the PH assumption test for the treatment term via Schoenfeld
# residuals. Returns NAs rather than failing: for a weighted Fine-Gray fit
# cox.zph() is not always computable.
testCrProportionalHazards <- function(fit) {
  tryCatch({
    zph <- cox.zph(fit)
    list(chiSq = zph$table["treatment", "chisq"],
         df = zph$table["treatment", "df"],
         p = zph$table["treatment", "p"])
  }, error = function(e) {
    list(chiSq = NA_real_, df = NA_real_, p = NA_real_)
  })
}


# Pull the treatment effect out of a fitted coxph object.
summarizeCrFit <- function(fit) {
  summaryFit <- summary(fit)
  coefRow <- summaryFit$coefficients["treatment", , drop = FALSE]
  ciRow <- summaryFit$conf.int["treatment", , drop = FALSE]
  phTest <- testCrProportionalHazards(fit)
  list(hr = as.numeric(ciRow[, "exp(coef)"]),
       ci95Lb = as.numeric(ciRow[, "lower .95"]),
       ci95Ub = as.numeric(ciRow[, "upper .95"]),
       logHr = as.numeric(coefRow[, "coef"]),
       seLogHr = as.numeric(coefRow[, "se(coef)"]),
       phChiSq = phTest$chiSq,
       phDf = phTest$df,
       phP = phTest$p)
}


# Fit a cause-specific Cox model, censoring competing deaths. The csHR
# estimates the rate of the event among those still alive and event-free
# (etiologic question).
fitCauseSpecificCox <- function(crData, stratified) {
  crData$csStatus <- as.integer(crData$status == 1)
  if (hasCrStrata(crData, stratified)) {
    fit <- coxph(Surv(time, csStatus) ~ treatment + strata(stratumId), data = crData)
  } else {
    fit <- coxph(Surv(time, csStatus) ~ treatment, data = crData)
  }
  return(summarizeCrFit(fit))
}


# Fit a Fine-Gray model via finegray() + weighted coxph(). The sHR
# corresponds to the CIF and so addresses the prognostic question of absolute
# risk in the presence of death.
fitFineGrayModel <- function(crData, stratified) {
  crData$fstatus <- factor(crData$status,
                           levels = 0:2,
                           labels = c("censor", "event", "competingDeath"))
  useStrata <- hasCrStrata(crData, stratified)
  if (useStrata) {
    fgInput <- crData[, c("time", "fstatus", "treatment", "stratumId")]
  } else {
    fgInput <- crData[, c("time", "fstatus", "treatment")]
  }
  fgData <- finegray(Surv(time, fstatus) ~ ., data = fgInput, etype = "event")
  if (useStrata) {
    fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ treatment + strata(stratumId),
                 data = fgData,
                 weights = fgwt)
  } else {
    fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ treatment,
                 data = fgData,
                 weights = fgwt)
  }
  return(summarizeCrFit(fit))
}


# Compute the cumulative incidence functions and Gray's test for the
# between-group comparison of the CIF of the event of interest.
#
# Disclosure control on the exported curves follows prepareKaplanMeier():
# the curves are truncated at the 90th percentile of follow-up, because the
# tail beyond that is driven by a handful of remaining subjects and each of
# its steps is effectively a single patient. On top of that the curves are
# reported on a fixed time grid rather than at the observed event times, so
# no exact event day of any individual is written out.
computeAndTestCif <- function(crData, cifTimeGridDays) {
  cifObj <- cmprsk::cuminc(ftime = crData$time,
                           fstatus = crData$status,
                           group = crData$treatment)
  if (is.null(cifObj$Tests) || !"1" %in% rownames(cifObj$Tests)) {
    stop("Gray's test is not available; both treatment groups must contribute events")
  }
  grayRow <- cifObj$Tests["1", , drop = TRUE]

  dataCutoff <- 0.9
  cutoff <- as.numeric(quantile(crData$time, dataCutoff))
  grid <- seq(from = 0, to = cutoff, by = cifTimeGridDays)
  grid <- grid[grid > 0]

  extractCif <- function(groupLabel, causeLabel) {
    key <- paste(groupLabel, causeLabel)
    group <- if (groupLabel == "1") "target" else "comparator"
    cause <- if (causeLabel == "1") "eventOfInterest" else "competingDeath"
    treatment <- if (groupLabel == "1") 1 else 0
    # Number still at risk at each time, as a count of follow-up times at or
    # beyond it. left.open = TRUE makes findInterval() count the times strictly
    # below, which keeps this O(n log n) rather than a scan per grid point.
    followUp <- sort(crData$time[crData$treatment == treatment])
    countAtRisk <- function(times) {
      length(followUp) - findInterval(times, followUp, left.open = TRUE)
    }
    # Every curve starts at the origin, which cuminc() does not return.
    origin <- data.frame(time = 0,
                         cif = 0,
                         seCif = 0,
                         atRisk = countAtRisk(0),
                         group = group,
                         cause = cause,
                         stringsAsFactors = FALSE)
    if (is.null(cifObj[[key]]) || length(grid) == 0) {
      return(origin)
    }
    curve <- cifObj[[key]]
    rbind(origin,
          data.frame(time = grid,
                     cif = evaluateStepFunction(curve$time, curve$est, grid),
                     seCif = sqrt(evaluateStepFunction(curve$time, curve$var, grid)),
                     atRisk = countAtRisk(grid),
                     group = group,
                     cause = cause,
                     stringsAsFactors = FALSE))
  }

  cifData <- rbind(extractCif("1", "1"),
                   extractCif("1", "2"),
                   extractCif("0", "1"),
                   extractCif("0", "2"))

  list(cifData = cifData,
       grayTestStat = as.numeric(grayRow["stat"]),
       grayDf = as.numeric(grayRow["df"]),
       grayP = as.numeric(grayRow["pv"]))
}


# Evaluate a right-continuous step function, which is what a cumulative
# incidence estimate is, on an arbitrary time grid. Grid points before the
# first jump take the value 0.
evaluateStepFunction <- function(times, values, grid) {
  index <- findInterval(grid, times)
  ifelse(index == 0, 0, values[pmax(index, 1)])
}


# Row-bind and write the accumulated results. fileSuffix identifies the
# (target, comparator, analysis) triplet; NULL writes the pooled files.
exportCompetingRiskTables <- function(result, crFolder, fileSuffix, minCellCount, silent) {
  countFields <- list(estimates = c("nTarget", "nComparator",
                                    "eventsTarget", "eventsComparator",
                                    "competingEventsTarget", "competingEventsComparator"),
                      cif = "atRisk")
  for (table in names(result)) {
    if (length(result[[table]]) == 0) {
      next
    }
    data <- do.call(rbind, result[[table]])
    for (field in countFields[[table]]) {
      data <- enforceMinCellValue(data, field, minCellCount, silent = silent)
    }
    stem <- switch(table,
                   estimates = "cr_estimates",
                   cif = "cr_cif",
                   gray = "cr_gray_test",
                   ph = "cr_ph_test")
    fileName <- if (is.null(fileSuffix)) {
      sprintf("%s.csv", stem)
    } else {
      sprintf("%s_%s.csv", stem, fileSuffix)
    }
    readr::write_csv(data, file.path(crFolder, fileName))
  }
  return(invisible(NULL))
}

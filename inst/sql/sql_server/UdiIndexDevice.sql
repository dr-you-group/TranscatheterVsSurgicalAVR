/************************************************************************
Copyright 2020 Observational Health Data Sciences and Informatics

This file is part of TavrSavr

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
************************************************************************/
{DEFAULT @cdm_database_schema = "cdm"}
{DEFAULT @cohort_database_schema = "cdm"}
{DEFAULT @cohort_table = "cohort"}
{DEFAULT @device_exposure_table = "cdm.device_exposure"}
{DEFAULT @valve_concept_id = 3661561}

-- The index device of each cohort entry: the exposure of the valve concept on
-- the cohort start date. A patient with more than one such record on that date
-- is returned as several rows and resolved in R.
SELECT c.cohort_definition_id,
	c.subject_id,
	c.cohort_start_date,
	de.device_concept_id,
	de.unique_device_id,
	op.observation_period_end_date
FROM @cohort_database_schema.@cohort_table c
INNER JOIN @device_exposure_table de
	ON de.person_id = c.subject_id
		AND de.device_exposure_start_date = c.cohort_start_date
		AND de.device_concept_id = @valve_concept_id
INNER JOIN @cdm_database_schema.observation_period op
	ON op.person_id = c.subject_id
		AND c.cohort_start_date >= op.observation_period_start_date
		AND c.cohort_start_date <= op.observation_period_end_date
WHERE c.cohort_definition_id IN (@cohort_ids);

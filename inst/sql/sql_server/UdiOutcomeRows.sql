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
{DEFAULT @cohort_database_schema = "cdm"}
{DEFAULT @cohort_table = "cohort"}

-- Every outcome cohort row of every person in any of the arm cohorts. Not
-- aggregated on purpose: the outcome cohorts are First-limited, so a row before
-- the index date is the person's ONLY row, and that has to be visible in R.
SELECT o.subject_id,
	o.cohort_definition_id,
	o.cohort_start_date AS outcome_date
FROM @cohort_database_schema.@cohort_table o
WHERE o.cohort_definition_id IN (@outcome_cohort_ids)
	AND o.subject_id IN (
		SELECT subject_id
		FROM @cohort_database_schema.@cohort_table
		WHERE cohort_definition_id IN (@cohort_ids)
		);

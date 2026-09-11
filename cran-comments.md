## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Notes for the reviewer

* The package has no compiled code and no system requirements.
* `Description` intentionally names the estimators it implements (Harrell's
  C, Uno's IPCW C). These are estimator names from the cited literature, not
  software names, so they are not single-quoted.
* Examples run in well under 5 seconds each; the bootstrap-bearing examples
  use small `n_boot` values for that reason.

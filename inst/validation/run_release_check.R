# MIMSTEM installed-package release check
# Run after library(MIMSTEM) or source directly from system.file().

stopifnot(requireNamespace("MIMSTEM", quietly = TRUE))

pkg_ver <- as.character(utils::packageVersion("MIMSTEM"))
cat("MIMSTEM version:", pkg_ver, "\n")

exports <- sort(getNamespaceExports("MIMSTEM"))
expected_exports <- sort(c("mims_tem", "mims_s2tem_ge", "mims_s2tem_pe"))
stopifnot(identical(exports, expected_exports))
cat("Public API check: passed\n")

ext <- system.file("extdata", package = "MIMSTEM")
stopifnot(nzchar(ext))
for (f in c("mims_example_tem.csv", "mims_example_s2tem.csv")) {
  stopifnot(file.exists(file.path(ext, f)))
}
cat("Installed example-data check: passed\n")

cit <- utils::citation("MIMSTEM")
stopifnot(length(cit) >= 1L)
cat("Citation check: passed\n")

vig <- utils::vignette(package = "MIMSTEM")
if (NROW(vig$results) > 0L) {
  stopifnot(any(vig$results[, "Item"] == "MIMSTEM-workflow"))
  cat("Vignette registration check: passed\n")
} else {
  warning("No installed vignette index was found in this installation.")
}

smoke <- system.file("validation", "run_smoke_test.R", package = "MIMSTEM")
stopifnot(nzchar(smoke), file.exists(smoke))
cat("Smoke-test script installed: passed\n")

cat("\nSession information:\n")
print(utils::sessionInfo())
cat("\nRelease metadata/API checks passed.\n")

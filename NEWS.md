# trapdiff 0.3.2

* Use DESeq2 normalized counts for ratio comparison
* Export ratio t-test and normalized counts

# trapdiff 0.3.1

* Added a `NEWS.md` file to track changes to the package.
* Added `filter_regex` options that allows to remove genes by their name
* Added a library size vs size factor debug plot that allows to assess problems with library size and size factors not correlating properly
* The interactive plot as part of the trapdiff standard output displays CPMs instead of TPMs and is filtered for all DE genes, no matter the comparison
* Export CPMs

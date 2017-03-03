ratesbydist <- read.csv("../results/ratebydist_Caenorhabditis.tsv", header=TRUE, sep="\t")
plot(
  ratesbydist$distance,
  ratesbydist$rate,
  main="Gene families in Caenorhabditis",
  xlab="Time in MYR since gene duplication",
  ylab="Per-site substitution rate",
  pch=19
)

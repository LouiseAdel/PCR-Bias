#Commands for Calculations for: Kruskal Wallis and Dunn's Test - Coverage
#Libraries ----
library(rstatix)
library(FSA)

#Coverage Statistics: ----
#Coverage Data - Median of Controls
#Protocol = P
data_BRCA2<-data.frame(
  Protocol = c(
    rep("P1", 10),
    rep("P2", 10),
    rep("P3", 10),
    rep("P4", 10),
    rep("P5", 10),
    rep("P6", 10),
    rep("P7", 2),
    rep("P8", 10)
  ),
  Coverage=c(
    41.98,65.75,70.33,47.86,68.99,77.94,57.31,70.34,53.27,76.94,
    27.23,28.22,11.54,13.05,25.97,39.17,38.49,26.37,9.43,11.56,
    16.46,28.46,26.24,14.61,40.03,20.98,89.34,39.14,44.56,54.91,
    24016.80,30919.10,17803.40,21585.40,20685.80,13299.50,19384.40,13024.00,26662.80,27946.40,
    6254.06,6240.18,6914.06,6531.28,6598.14,5960.34,6744.59,6934.77,7182.01,7358.05,
    1530.08,1165.17,1332.74,3221.22,1421.69,1174.82,1054.39,4584.47,13093.50,5583.94,
    145.27,124.23,
    37.37,44.60,57.41,67.90,31.60,57.04,112.61,84.68,59.28,84.81
  )
)

data_BRCA1<-data.frame(
  Protocol = c(
    rep("P1", 10),
    rep("P2", 10),
    rep("P3", 10),
    rep("P4", 10),
    rep("P5", 10),
    rep("P6", 10),
    rep("P7", 2),
    rep("P8", 10)
  ),
  Coverage=c(
    154.80,194.41,251.15,146.34,195.45,292.25,213.23,235.08,187.53,257.60,
    154.26,89.80,58.89,43.05,109.71,198.23,140.48,113.96,65.30,52.63,
    57.83,52.54,66.97,42.92,90.59,58.18,112.16,97.85,120.21,96.10,
    44526.10,196.04,835.30,27519.30,53445.30,44075.10,23005.40,667.49,37630.80,33398.10,
    4020.08,3308.39,3245.74,3870.37,3352.71,3861.98,3764.09,3882.39,4863.65,2940.89,
    3381.79,1752.89,3586.46,5670.82,3778.59,2662.45,1471.28,11009.80,16259.20,14112.80,
    443.11,426.59,
    53.15,79.17,57.39,101.22,34.87,59.49,114.15,88.80,62.41,83.51
  ) 
)

#Kruskal Wallis Test ----
# Perform Kruskal-Wallis test
kruskal_result_BRCA1 <- kruskal.test(Coverage ~ Protocol, data = data_BRCA1)

# Display results
print(kruskal_result_BRCA1)
#Result: Kruskal-Wallis chi-squared = 60.003, df = 7, p-value = 1.508e-10

kruskal_result_BRCA2 <- kruskal.test(Coverage ~ Protocol, data = data_BRCA2)

# Display results
print(kruskal_result_BRCA2)
#Result: Kruskal-Wallis chi-squared = 64.59, df = 7, p-value = 1.818e-11
#Conclusion: Need Dunn's Test for both BRCA1 and BRCA2

#Dunn's Test ----
dunn_BRCA1 <- dunnTest(Coverage ~ Protocol, data = data_BRCA1, method = "holm")
print(dunn_BRCA1)

dunn_BRCA2 <- dunnTest(Coverage ~ Protocol, data = data_BRCA2, method = "holm")
print(dunn_BRCA2)

# BRCA1 significant differences
cat("=== BRCA1 SIGNIFICANT DIFFERENCES (P.adj < 0.05) ===\n")
brca1_significant <- dunn_BRCA1$res[dunn_BRCA1$res$P.adj < 0.05, ]
if(nrow(brca1_significant) > 0) {
  print(brca1_significant)
  cat("\nNumber of significant comparisons for BRCA1:", nrow(brca1_significant), "\n")
} else {
  cat("No significant differences found for BRCA1\n")
}

cat("\n")

# BRCA2 significant differences
cat("=== BRCA2 SIGNIFICANT DIFFERENCES (P.adj < 0.05) ===\n")
brca2_significant <- dunn_BRCA2$res[dunn_BRCA2$res$P.adj < 0.05, ]
if(nrow(brca2_significant) > 0) {
  print(brca2_significant)
  cat("\nNumber of significant comparisons for BRCA2:", nrow(brca2_significant), "\n")
} else {
  cat("No significant differences found for BRCA2\n")
}

# Combined
all_significant <- rbind(
  data.frame(Gene = "BRCA1", brca1_significant),
  data.frame(Gene = "BRCA2", brca2_significant)
)

cat("\n=== COMBINED SIGNIFICANT RESULTS ===\n")
print(all_significant)

#Result: Significant differences: 
# Gene Comparison         Z      P.unadj        P.adj
#5   BRCA1    P2 - P4 -4.444688 8.801952e-06 2.112469e-04
#6   BRCA1    P3 - P4 -4.904115 9.384957e-07 2.533938e-05
#8   BRCA1    P2 - P5 -3.857049 1.147642e-04 2.295284e-03
#9   BRCA1    P3 - P5 -4.316476 1.585400e-05 3.487881e-04
#12  BRCA1    P2 - P6 -3.921155 8.812555e-05 1.850637e-03
#13  BRCA1    P3 - P6 -4.380582 1.183628e-05 2.722344e-04
#25  BRCA1    P4 - P8  5.064380 4.097315e-07 1.147248e-05
#26  BRCA1    P5 - P8  4.476741 7.579110e-06 1.894778e-04
#27  BRCA1    P6 - P8  4.540847 5.602866e-06 1.456745e-04
#4   BRCA2    P1 - P4 -4.060052 4.906189e-05 1.079362e-03
#51  BRCA2    P2 - P4 -6.293080 3.112282e-10 8.714388e-09
#61  BRCA2    P3 - P4 -5.481070 4.227623e-08 1.141458e-06
#81  BRCA2    P2 - P5 -5.128486 2.920814e-07 7.594117e-06
#91  BRCA2    P3 - P5 -4.316476 1.585400e-05 3.963501e-04
#121 BRCA2    P2 - P6 -4.284423 1.832142e-05 4.397141e-04
#131 BRCA2    P3 - P6 -3.472412 5.158032e-04 1.083187e-02
#251 BRCA2    P4 - P8  4.209632 2.557866e-05 5.883091e-04
#261 BRCA2    P5 - P8  3.045039 2.326504e-03 4.653008e-02



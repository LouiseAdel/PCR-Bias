#Commands for Calculation of Fisher's Exact Test and Log2 Fold-Change - Isoform Rediscovery, PCR Bias

#Comparison of Δ10q levels in samples (S3 and S8) and controls (Sample 15 and 19) 
#from Protocols (P) to levels detected by Direct RNA-seq (P7)

#Libraries ----
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)

#PCR Bias Statistics ----
#Individual Protocols ----
#Read data file
data_controls <- read_excel("PCRBias_delta10q_Levels_controlsonly.xlsx")

#Combine control counts
combined_counts <- data_controls %>%
  group_by(Protocol) %>%
  summarise(
    Delta10q_Reads = sum(Delta10q_Reads),
    FullLength_Reads = sum(FullLength_Reads),
    .groups = "drop"
  ) %>%
  mutate(
    Delta10q_Percentage = 100 * Delta10q_Reads / (Delta10q_Reads + FullLength_Reads)
  )

# Separate out Protocol 7 (Direct RNA-seq)
direct_rna <- combined_counts %>%
  filter(Protocol == 7) %>%
  select(Delta10q_Reads, FullLength_Reads, Delta10q_Percentage)

#Fisher's Exact Test

fisher_results <- combined_counts %>%
  filter(Protocol != 7) %>%
  rowwise() %>%
  mutate(
    Fisher_p = fisher.test(matrix(
      c(Delta10q_Reads, FullLength_Reads,
        direct_rna$Delta10q_Reads, direct_rna$FullLength_Reads),
      nrow = 2, byrow = TRUE
    ))$p.value
  ) %>%
  ungroup()


#Log2 Fold-Change
logfc_data <- combined_counts %>%
  filter(Protocol != 7) %>%
  mutate(
    Log2FC = log2((Delta10q_Percentage + 1e-6) / (direct_rna$Delta10q_Percentage + 1e-6))
  )


#Comnbined Log2FC and Fisher
results_individual <- logfc_data %>%
  left_join(fisher_results %>% select(Protocol, Fisher_p), by = "Protocol")

# Optional: Print results
print(results_individual)

#Result - Individual Protocols:
#Protocol Delta10q_Reads FullLength_Reads Delta10q_Percentage  Log2FC Fisher_p
#1           5926            73204                7.49  0.0158 1   e+ 0
#2             14              106               11.7   0.655  1.79e- 1
#3             70              373               15.8   1.09   6.13e- 4
#4          29521            66809               30.6   2.05   1.89e-22
#5           3723             8165               31.3   2.08   6.61e-23
#6           2921            23487               11.1   0.578  4.97e- 2
#8              8              696                1.14 -2.70   6.68e- 7


#Grouped by amount of PCR:  ----

data_controls_group <- data_controls %>%
  mutate(
    Protocol_Group = case_when(
      Protocol %in% c(1, 2, 3, 6) ~ "PCR_cDNA_only",
      Protocol %in% c(4, 5) ~ "PCR_cDNA_plus_enrichment",
      Protocol == 7 ~ "Direct_RNA",
      TRUE ~ "Other"
    )
  )

# Combine control counts per protocol group
combined_groups <- data_controls_group %>%
  filter(Protocol_Group != "Other") %>%
  group_by(Protocol_Group) %>%
  summarise(
    Delta10q_Reads = sum(Delta10q_Reads),
    FullLength_Reads = sum(FullLength_Reads),
    .groups = "drop"
  ) %>%
  mutate(
    Delta10q_Percentage = 100 * Delta10q_Reads / (Delta10q_Reads + FullLength_Reads)
  )

# Extract Direct RNA group counts
direct_rna_group <- combined_groups %>%
  filter(Protocol_Group == "Direct_RNA") %>%
  select(Delta10q_Reads, FullLength_Reads, Delta10q_Percentage)

# Log2 Fold-Change for grouped protocols vs Direct RNA
logfc_data_group <- combined_groups %>%
  filter(Protocol_Group != "Direct_RNA") %>%
  mutate(
    Log2FC = log2((Delta10q_Percentage + 1e-6) / (direct_rna_group$Delta10q_Percentage + 1e-6))
  )

# Fisher's Exact Test for grouped protocols vs Direct RNA
fisher_results_group <- combined_groups %>%
  filter(Protocol_Group != "Direct_RNA") %>%
  rowwise() %>%
  mutate(
    Fisher_p = fisher.test(matrix(
      c(Delta10q_Reads, FullLength_Reads,
        direct_rna_group$Delta10q_Reads, direct_rna_group$FullLength_Reads),
      nrow = 2, byrow = TRUE
    ))$p.value
  ) %>%
  ungroup()

# Combine Log2FC and Fisher p-values for grouped data
results_grouped <- logfc_data_group %>%
  left_join(fisher_results_group %>% select(Protocol_Group, Fisher_p), by = "Protocol_Group")

# Print grouped protocol results
print(results_grouped)

#Result - Grouped by amount of PCR:
#Protocol_Group           Delta10q_Reads FullLength_Reads Delta10q_Percentage Log2FC Fisher_p
#PCR_cDNA_only                      8931            97170                8.42  0.184 6.01e- 1
#PCR_cDNA_plus_enrichment          33244            74974               30.7   2.05  1.79e-22



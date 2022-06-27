# https://www.tidyverse.org/
library(tidyverse)

# Set the directory holding input files
# Try file.choose() to find yours
input_dir <- "/Users/tkeitt/Desktop/BioPopWork/inputTXT"

# File with jay IDs
# The assumption is that the order of ids in this file matches the order of rows in the reads file
# This turns out not to be true, so don't rely on this output
# We need a file that maps barcodes to jay ids
jay_ids <- read_table(file.path(input_dir, "header-order-genotype-names.txt"), col_names = FALSE)

# This uses piping; the output of each function is piped to the next
# Note the read stats file in the repo has some duplicated rows, so we deduplicate
read_table(file.path(input_dir, "MappedReads.txt"), col_names = FALSE) %>%
  filter(!duplicated(.)) %>%
  select(!X1) %>%
  rename(total = X2, mapped = X3, unmapped = X4, good = X5) %>%
  mutate(id = unlist(jay_ids)) %>%
  pivot_longer(!id, names_to = "quantity", values_to = "count") -> readcount

# Bar chart with ggplot
ggplot(readcount, aes(x = id, fill = quantity, y = count)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90))


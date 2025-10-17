# Install if needed
if (!requireNamespace("KEGGREST", quietly = TRUE)) BiocManager::install("KEGGREST")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("httr", quietly = TRUE)) BiocManager::install("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) BiocManager::install("jsonlite")
if (!requireNamespace("tidyr", quietly = TRUE)) BiocManager::install("tidyr")
if (!requireNamespace("dplyr", quietly = TRUE)) BiocManager::install("dplyr")
if (!requireNamespace("gt", quietly = TRUE)) install.packages("gt")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")

library(KEGGREST)
library(readxl)
library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(gt)
library(scales)
library(viridis)

# Read Excel file
data <- read_excel("~/Documents/UM/Network Biology/Project/reactions.xlsx")   # Replace with your file name

# View first rows
head(data)

syn3a <- data$syn3A
pneumoniae <- data$pneumoniae

get_bigg_info <- function(bigg_id) {
  # Construct the API URL
  url <- paste0("http://bigg.ucsd.edu/api/v2/universal/reactions/", bigg_id)
  
  # Make the GET request
  response <- GET(url)
  
  # Check if the request was successful
  if (status_code(response) == 200) {
    # Parse the JSON content
    content <- fromJSON(rawToChar(response$content), flatten = TRUE)
    
    # Extract name with safe handling
    name <- ifelse(!is.null(content$name), content$name, NA_character_)
    
    # Extract KEGG reaction ID with better handling
    kegg_reaction_id <- NA_character_
    
    if (!is.null(content$database_links)) {
      # Check for KEGG Reaction links
      if (!is.null(content$database_links$KEGG) && !is.null(content$database_links$KEGG$reaction)) {
        kegg_reaction_id <- paste(content$database_links$KEGG$reaction, collapse = "; ")
      }
      # Alternative: check for "KEGG Reaction" (with space)
      else if (!is.null(content$database_links$`KEGG Reaction`)) {
        kegg_reaction_id <- paste(content$database_links$`KEGG Reaction`, collapse = "; ")
      }
      # Alternative: sometimes it's in the identifiers
      else if (!is.null(content$database_links$KEGG) && !is.null(content$database_links$KEGG$id)) {
        kegg_reaction_id <- paste(content$database_links$KEGG$id, collapse = "; ")
      }
    }
    
    return(data.frame(reaction_id = bigg_id,
                      reaction_name = name,
                      kegg_reaction_id = kegg_reaction_id,
                      stringsAsFactors = FALSE))
  } else {
    # If the ID was not found, return a row with NAs
    warning("Reaction ID '", bigg_id, "' not found in BiGG database.")
    return(data.frame(reaction_id = bigg_id,
                      reaction_name = NA_character_,
                      kegg_reaction_id = NA_character_,
                      stringsAsFactors = FALSE))
  }
}

# Use lapply to run the function for each ID and combine the results
syn3a_rxn_list <- lapply(syn3a[1:338], get_bigg_info)
pneumoniae_rxn_list <- lapply(pneumoniae, get_bigg_info)

# Combine the list of dataframes into one
syn3a_rxn_df <- do.call(rbind, syn3a_rxn_list)
pneumoniae_rxn_df <- do.call(rbind, pneumoniae_rxn_list)

syn3a_rxn_df$kegg_reaction_id <- sapply(syn3a_rxn_df$kegg_reaction_id, function(x) strsplit(x, ";")[[1]][1])
pneumoniae_rxn_df$kegg_reaction_id <- sapply(pneumoniae_rxn_df$kegg_reaction_id, function(x) strsplit(x, ";")[[1]][1])

sum(!is.na(syn3a_rxn_df$reaction_name[1:338]))
sum(!is.na(pneumoniae_rxn_df$reaction_name))

sum(!is.na(syn3a_rxn_df$kegg_reaction_id[1:338]))
sum(!is.na(pneumoniae_rxn_df$kegg_reaction_id))




syn3a_rxn_df <- syn3a_rxn_df[!is.na(syn3a_rxn_df$reaction_name),]
pneumoniae_rxn_df <- pneumoniae_rxn_df[!is.na(pneumoniae_rxn_df$reaction_name),]

pneumoniae_rxn_df$kegg_reaction_id[207] <- "R00132,R10092"

# Split it into two rows
pneumoniae_rxn_df <- pneumoniae_rxn_df %>%
  separate_rows(kegg_reaction_id, sep = ",") %>%
  mutate(kegg_reaction_id = trimws(kegg_reaction_id))  # remove whitespace


# Remove na from reactions
syn3a_reactions <- syn3a_rxn_df$kegg_reaction_id
syn3a_reactions <- syn3a_reactions[!is.na(syn3a_reactions)]

pneumoniae_reactions <- pneumoniae_rxn_df$kegg_reaction_id
pneumoniae_reactions <- pneumoniae_reactions[!is.na(pneumoniae_reactions)]

library(KEGGREST)

pathway_id <- function(rxns){
  # Map reactions to pathways
  reaction_to_pathway <- sapply(rxns, function(rxn) {
    pathways <- keggLink("pathway", rxn)
    pathways <- sub("path:", "", pathways)  # clean the KEGG IDs
    return(pathways)
  })
  
  # Convert list to data frame (reaction - pathway pairs)
  rxn_pathway_df <- stack(reaction_to_pathway)
  colnames(rxn_pathway_df) <- c("pathway_id", "reaction_id")
  
  # Remove empty entries
  rxn_pathway_df <- rxn_pathway_df[rxn_pathway_df$pathway_id != "", ]
  
  # Keep only the "map" pathways (remove "rn" duplicates)
  return(rxn_pathway_df <- rxn_pathway_df[grepl("^map", rxn_pathway_df$pathway_id), ])
}

syn3a_pathway_df <-pathway_id(syn3a_reactions)
pneumoniae_pathway_df <-pathway_id(pneumoniae_reactions)


library(dplyr)
library(ggplot2)

# Step 1: Count total reactions per pathway in pneumoniae
pneumoniae_counts <- pneumoniae_pathway_df %>%
  group_by(pathway_id) %>%
  summarize(total_rxn = n_distinct(reaction_id))

# Step 2: Count distinct Syn3A reactions per pathway
syn3a_counts <- syn3a_pathway_df %>%
  group_by(pathway_id) %>%
  summarize(syn3a_rxn = n_distinct(reaction_id))

# Step 3: Merge and calculate fraction retained
pathway_retention <- pneumoniae_counts %>%
  left_join(syn3a_counts, by = "pathway_id") %>%
  mutate(syn3a_rxn = ifelse(is.na(syn3a_rxn), 0, syn3a_rxn),
         fraction_retained = pmin(syn3a_rxn / total_rxn, 1))

# Step 4: Sort pathways by fraction_retained for plotting
pathway_retention <- pathway_retention %>%
  arrange(fraction_retained) %>%
  mutate(pathway_id = factor(pathway_id, levels = pathway_id))


ggplot(pathway_retention, aes(x = fraction_retained, y = pathway_id)) +
  geom_segment(aes(x = 0, xend = fraction_retained, yend = pathway_id), color = "gray70") +
  geom_point(aes(color = fraction_retained), size = 3) +
  scale_color_viridis(option = "viridis", direction = -1) +  # color-blind friendly
  theme_minimal(base_size = 12) +
  labs(
    x = "Fraction of reactions retained in JCVI-syn3A",
    y = "KEGG Pathway ID",
    title = "Pathway Retention Between JCVI-syn3A and M. pneumoniae"
  ) +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 13),
    legend.position = "none"
  )










library(gt)
library(dplyr)
library(scales)



# Separate top and bottom 5 pathways
top_pathways <- pathway_retention %>%
  arrange(desc(fraction_retained), desc(total_rxn)) %>%
  slice_head(n = 5) %>%
  mutate(Category = "Most Retained")

bottom_pathways <- pathway_retention %>%
  arrange(fraction_retained, desc(total_rxn))  %>%
  slice_head(n = 5) %>%
  mutate(Category = "Least Retained")

pathway_names <- data.frame(
  pathway_id = c(top_pathways$pathway_id, bottom_pathways$pathway_id),
  pathway_name = c("Aminoacyl-tRNA biosynthesis",
                   "Carbon fixation by Calvin cycle",
                   "Pentose phosphate pathway",
                   "Galactose metabolism",
                   "Amino sugar and nucleotide sugar metabolism",
                   "Ascorbate and aldarate metabolism",
                   "Glycerophospholipid metabolism",
                   "Arginine biosynthesis",
                   "Pantothenate and CoA biosynthesis",
                   "Teichoic acid biosynthesis"
  )
)



# Combine and join pathway names
pathway_table <- bind_rows(top_pathways, bottom_pathways) %>%
  left_join(pathway_names, by = "pathway_id") %>%
  select(Category, pathway_id, pathway_name, total_rxn, syn3a_rxn) %>%
  rename(
    `KEGG Pathway ID` = pathway_id,
    `Pathway Name` = pathway_name,
    `Total Reactions` = total_rxn,
    `Retained Reactions` = syn3a_rxn
  )

# Build gt table
pathway_table %>%
  gt(groupname_col = "Category") %>%
  tab_header(
    title = md("**Retained Pathways Between JCVI-syn3A and *M. pneumoniae***")  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_title(groups = "title")
  ) %>%
  tab_options(
    table.font.size = 12,
    data_row.padding = px(4),
    column_labels.font.weight = "bold",
    heading.title.font.weight = "bold",
    heading.align = "center",
    row_group.as_column = TRUE
  )


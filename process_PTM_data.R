
# 1. Setup

# Example command
# Rscript /relquant/pp.R
# -d data_package_number
# -o data/test_data_global

library(optparse)
library(MSnID)
library(PlexedPiper)

option_list <- list(
  make_option(c("-d", "--data_package_number"), type="character", default=NULL, 
              help="Data package number", metavar="character"),
  make_option(c("-g", "--global_results_ratio"), type="character", default=NULL, 
              help="Data package number", metavar="character"),
  make_option(c("-o", "--output_folder"), type="character", default=NULL, 
              help="PlexedPiper output folder (Crosstabs)", metavar="character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data_package_number) |
    is.null(opt$global_results_ratio) |
    is.null(opt$output_folder)) {
  print_help(opt_parser)
  stop("3 arguments are required", call.=FALSE)
}

data_package_number <- opt$data_package_number
path_to_global_results_ratio <- opt$global_results_ratio
output_folder <- opt$output_folder

# To DEBUG ---------------------------------------------------------------------
# # Global
# data_package_number = 3625
# path_to_global_results_ratio = "data/test_data_global/results_ratio.txt"
# output_folder = "data/test_data_global"
# To DEBUG ---------------------------------------------------------------------

if(!dir.exists(file.path(output_folder))){
  dir.create(file.path(output_folder), recursive = TRUE)
}

message("- Beginning processing of global data package ", data_package_number, "...")

# 2. Read study design files

message("- Fetch study design tables")

study_design <- get_study_design_by_dataset_package(data_package_number)

fractions  <- study_design$fractions
samples    <- study_design$samples
references <- study_design$references

write.table(fractions, file=file.path(output_folder, "study_design/fractions.txt"),
            quote=F, sep="\t", row.names=F)
write.table(samples, file=file.path(output_folder, "study_design/samples.txt"),
            quote=F, sep="\t", row.names=F)
write.table(references, file=file.path(output_folder, "study_design/references.txt"),
            quote=F, sep="\t", row.names=F)

# 3. Read MS-GF+ data

message("- Prepare MS/MS IDs")
message("   + Read the MS-GF+ output\n")

path_to_MSGF_data <- file.path(output_folder, "msgfData_original.RData")

if (file.exists(path_to_MSGF_data)) {
  load(path_to_MSGF_data)
  # check job records match dataset
  job_records <- get_job_records_by_dataset_package(data_package_number)
  if (!setequal(job_records$Dataset, msnid$Dataset)) {
    warning("Datasets in MS-GF+ data and DMS do not match!")
  }
} else {
  msnid <- read_msms_data_from_DMS(data_package_number)
  save(msnid, file=path_to_MSGF_data)
}

message(show(msnid))


message("   + Read Ascore output")
ascore <- get_AScore_results(data_package_number)

message("   + Select best PTM location by AScore")
msnid <- best_PTM_location_by_ascore(msnid, ascore)

if (proteomics == "ph") {
  msnid <- apply_filter(msnid, "grepl(\"\\\\*\", peptide)")
} else if(proteomics %in% c("ac", "ub")) {
  msnid <- apply_filter(msnid, "grepl(\"\\\\#\", peptide)")
} else {
  stop("proteomics variable not supported")
}

message("   + FDR filter")
msnid <- filter_msgf_data_peptide_level(msnid, 0.01)

message("   + Remove decoy sequences")
msnid <- apply_filter(msnid, "!isDecoy")

message("   + Concatenating redundant RefSeq matches")
msnid <- assess_redundant_protein_matches(msnid)

message("   + Assessing non-inferable proteins")
msnid <- assess_noninferable_proteins(msnid)

message("   + Inference of parsimonius set")
global_results_ratio <- read.table(path_to_global_results_ratio, header=T, sep="\t")
global_protein_ids <- unique(global_results_ratio$protein_id)
msnid <- infer_parsimonious_accessions(msnid, prior=global_protein_ids)

message("   + Mapping sites to protein sequence")
suppressMessages(
  fst <- Biostrings::readAAStringSet(fasta_file)
)
names(fst) <- sub("^(\\S*)\\s.*", "\\1", names(fst))

if(proteomics == "ph") {
  msnid <- map_mod_sites(msnid, 
                         fst, 
                         accession_col = "accession", 
                         peptide_mod_col = "peptide", 
                         mod_char = "*",
                         site_delimiter = "lower")
} else if (proteomics %in% c("ac", "ub")) {
  msnid <- map_mod_sites(msnid, 
                         fst, 
                         accession_col = "accession", 
                         peptide_mod_col = "peptide", 
                         mod_char = "#",
                         site_delimiter = "lower")
} else {
  stop("proteomics variable not supported")
}


message("   + Map flanking sequences")
msnid <- extract_sequence_window(msnid, fst)

path_to_MASIC_data <- file.path(output_folder, "masicData_filtered.RData")
save(masic_data, file=path_to_MASIC_data)

# 4. Read MASIC data

message("- Prepare reporter ion intensities")
message("   + Read MASIC ouput\n")

path_to_MASIC_data <- file.path(output_folder, "masicData_original.RData")
if (file.exists(path_to_MASIC_data)) {
  load(path_to_MASIC_data)
  # check job records match dataset
  job_records <- get_job_records_by_dataset_package(data_package_number)
  if (!setequal(job_records$Dataset, masic_data$Dataset)) {
    warning("Datasets in MASIC data and DMS do not match!")
  }
} else {
  masic_data <- read_masic_data_from_DMS(data_package_number,
                                         interference_score = TRUE)
  save(masic_data, file=path_to_MASIC_data)
}

if (!setequal(fractions$Dataset, masic_data$Dataset)) {
  warning("Datasets in MASIC output and 'fractions.txt' do not match!")
}

message("   + Filtering MASIC data")
masic_data <- filter_masic_data(masic_data, 0.5, 0)

path_to_MASIC_data <- file.path(output_folder, "masicData_filtered.RData")
save(masic_data, file=path_to_MASIC_data)

# 5. Create crosstabs

message("- Create Reporter Ion Intensity Results")
rii_peptide <- make_rii_peptide_ph(msnid = msnid, 
                                   masic_data = masic_data, 
                                   fractions = fractions, 
                                   samples = samples, 
                                   references = references, 
                                   org_name = "Rattus norvegicus")

message("- Create Ratio Results")
results_ratio <- make_results_ratio_ph(msnid =  msnid, 
                                       masic_data = masic_data, 
                                       fractions = fractions, 
                                       samples = samples, 
                                       references = references, 
                                       org_name = "Rattus norvegicus")

# 6. Check results and clean

message("- Save results")

path_to_rii_peptide <- file.path(output_folder, "results_RII-peptide.txt")
message("   + Writing ", path_to_rii_peptide)
write.table(rii_peptide,
            file = path_to_rii_peptide,
            sep="\t",
            row.names = FALSE,
            quote = FALSE)

path_to_results_ratio <- file.path(output_folder, "results_ratio.txt")
message("   + Writing ", path_to_results_ratio)
write.table(results_ratio,
            file = path_to_results_ratio,
            sep="\t",
            row.names = FALSE,
            quote = FALSE)

message("- Done!")

unlink(".Rcache", recursive=TRUE)

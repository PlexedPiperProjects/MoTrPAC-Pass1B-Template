
# 1. Setup

# Example command
# Rscript /relquant/pp.R
# -d data_package_number
# -o data/test_data_global

library(MSnID)
library(PlexedPiper)

# To process from command line

library(optparse)

option_list <- list(
  make_option(c("-d", "--data_package_number"), type="character", default=NULL, 
              help="Data package number", metavar="character"),
  make_option(c("-o", "--output_folder"), type="character", default=NULL, 
              help="PlexedPiper output folder (Crosstabs)", metavar="character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$data_package_number) |
    is.null(opt$output_folder)) {
  print_help(opt_parser)
  stop("2 arguments are required", call.=FALSE)
}

data_package_number <- opt$data_package_number
output_folder <- opt$output_folder

# To process by hand
# data_package_number = 3606
# output_folder = "data/test_data_global"

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

if (!setequal(fractions$Dataset, msnid$Dataset)) {
  warning("Datasets in MS-GF+ output and 'fractions.txt' do not match!")
}

message("   + Correct for isotope selection error\n")
msnid <- correct_peak_selection(msnid)

message("   + MS/MS ID filter and peptide level\n")
msnid <- filter_msgf_data_peptide_level(msnid, 0.01)

message(show(msnid))

message("   + Read FASTA file\n")
path_to_FASTA <- path_to_FASTA_used_by_DMS(data_package_number)
path_to_FASTA_copy <- file.path(output_folder, basename(path_to_FASTA))
if (!file.exists(path_to_FASTA_copy)) file.copy(path_to_FASTA, output_folder)

message("   + MS/MS ID filter at protein level\n")
msnid <- compute_num_peptides_per_1000aa(msnid, path_to_FASTA)
msnid <- filter_msgf_data_protein_level(msnid, 0.01)

message(show(msnid))

message("   + Remove decoy accessions\n")
msnid <- apply_filter(msnid, "!isDecoy")

message(show(msnid))

message("   + Concatenating redundant RefSeq matches\n")
msnid <- assess_redundant_protein_matches(msnid)

message("   + Assessing non-inferable proteins\n")
msnid <- assess_noninferable_proteins(msnid)

message("   + Inference of parsimonious protein set\n")
msnid <- infer_parsimonious_accessions(msnid)

message(show(msnid))

message("   + Compute accsesion coverage\n")
suppressMessages(
  fst <- Biostrings::readAAStringSet(path_to_FASTA)
)
names(fst) <- sub("^(\\S*)\\s.*", "\\1", names(fst))
msnid <- compute_accession_coverage(msnid, fst)

path_to_MSGF_data <- file.path(output_folder, "msgfData_filtered.RData")
save(msnid, file=path_to_MSGF_data)

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
rii_peptide <- make_rii_peptide_gl(msnid = msnid, 
                                   masic_data = masic_data, 
                                   fractions = fractions, 
                                   samples = samples, 
                                   references = references, 
                                   org_name = "Rattus norvegicus")

message("- Create Ratio Results")
results_ratio <- make_results_ratio_gl(msnid =  msnid, 
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

# Finally, copy self into output folder


message("- Done!")

unlink(".Rcache", recursive=TRUE)



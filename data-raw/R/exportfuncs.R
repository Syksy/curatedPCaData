###
#
# Variable functions helping exporting of data/information from curatedPCaData
#
###

#' Export ExperimentHub-friendly metadata csv prototype and MAE-slot Rds files
#'
#' This function greps over MAE-objects available via curatedPCaData and starts constructing ExperimentHub-compatible metadata.csv, first as a data.frame
#' Further, individual MAE-slots are exported as Rds files if
#'
#' @param timestamp Time stamp string in format "YYYYMMDD", and if missing extracted from Sys.time()
#' @param path Path to work in
#' @param export Should also the .Rds files be exported in addition to the metadata.csv
#' @param ... Additional parameters
#'
#' @examples
#' # metadata.csv-file
#' metadat <- curatedPCaData:::export_data(timestamp = "20230215", export = FALSE)
#' write.csv(metadat, file = "metadata.csv", quote = TRUE, row.names = FALSE)
#' # Rds-files
#' curatedPCaData:::export_data(timestamp = "20230215", export = TRUE)
#'
#' @noRd
#' @keywords internal
export_data <- function(timestamp, # Time stamp to append in file names
                        path = getwd(), # Location on drive to save the files into
                        export = TRUE, # Whether the actual .Rds objects should be exported; if FALSE, the relevant metadata is created as if the Rds hada been created
                        ...) {
  # If user does not provide a time stamp, construct one from the system time ('%Y%m%d'-format)
  if (missing(timestamp)) {
    timestamp <- gsub("-", "", substr(Sys.time(), 0, 10))
  }

  # List MAEs available in curatedPCaData
  # Get available maes in curatedPCaData package and specify a path to the folder where to save
  maes <- grep("mae", utils::data(package = "curatedPCaData")$result, value = TRUE)
  # Load MAEs
  data(list = maes, package = "curatedPCaData")

  # Metadata fields that are pre-known; reference list from https://bioconductor.org/packages/3.17/bioc/vignettes/HubPub/inst/doc/CreateAHubPackage.html#example-metadata.csv-file-and-more-information
  ll <- list(
    # "genome: 'character(1)'. Genome. Can be NA.
    genome = NA_character_,
    # BiocVersion: 'character(1)'. The first Bioconductor version the resource was made available for. Unless removed from the hub, the resource will be available for all versions greater than or equal to this field. This generally is the current devel version of Bioconductor.
    biocv = "3.17", # BioConductor version
    # Species: 'character(1)'. Species. For help on valid species see 'getSpeciesList, validSpecies, or suggestSpecies.' Can be NA.
    species = "Homo sapiens", # Human samples only
    # TaxonomyId: 'character(1)'. Taxonomy ID. There are checks for valid taxonomyId given the Species which produce warnings. See GenomeInfoDb::loadTaxonomyDb() for full validation table. Can be NA.
    taxonomyid = 9606,
    # Coordinate_1_based: 'logical'. TRUE if data are 1-based. Can be NA.
    coordinate_1_based = NA,
    # Maintainer: 'character(1)'. Maintainer name and email in the following format: Maintainer Name username@address.
    maintainer = "Teemu Daniel Laajala <teelaa@utu.fi>" # Desired maintainer field
  )

  # Listing for description fields, used with exact names for data slots
  descriptions <- c(
    # CNA copy number alteration fields with differing levels of information
    cna.gistic = "Copy number alteration GISTIC",
    cna.logr = "Copy number alteration log-ratios",
    # GEX - Gene expression fields with different processing pipelines
    gex.relz = "Gene expression (Relative z-score)",
    gex.logq = "Gene expression (Log-quantile norm)",
    gex.rma = "Gene expression (RMA)",
    gex.logr = "Gene expression (Log-ratios)",
    gex.rsem.log = "Gene expression (Log-RSEM)",
    # MUT - Mutation data
    mut = "Mutation",
    # Derived variables
    cibersort = "Deconvolution using CIBERSORTx",
    xcell = "Deconvolution using xCell",
    epic = "Deconvolution using EPIC",
    quantiseq = "Deconvolution using quanTIseq",
    mcp = "Deconvolution using MCP counter",
    estimate = "Deconvolution using ESTIMATE",
    scores = "Various gene expression risk and marker scores"
  )
  # Format an empty df with required fields
  metadata <- data.frame()

  for (i in 1:length(maes)) {
    mae <- maes[i] # Get specific mae & separate the study name for future naming
    study <- gsub("mae_", "", mae)
    ex <- names(MultiAssayExperiment::experiments(get(mae))) # Get the names of object ("experiments") in mae

    # Study specific conditional fields
    ll_study <- list(
      sourcetype = "",
      sourceurl = "",
      sourceversion = "",
      dataprovider = ""
    )

    # Go through all "experiments" and save the object into its own .RDs file
    for (x in ex) {
      oname <- paste0(study, "_", x, "_", timestamp) # object name
      fname <- paste0(oname, ".rds") # file name

      assign(oname, get(mae)[[x]]) # Assign the object to the specific name
      if (export) {
        for (index in 1:length(oname)) {
          saveRDS(get(oname[index]), file = paste0(path, "/", fname[index]))
        }
      }

      # Add info to metadata
      descrp <- paste(oname, descriptions[x], "data of", study, "cohort in curatedPCaData package", sep = " ")

      # Temporary vector to rbind into the full metadata data.frame
      tmp <- c(
        Title = oname,
        Description = descrp,
        BiocVersion = ll$biocv,
        Genome = ll$genome,
        SourceType = "",
        SourceUrl = "",
        SourceVersion = "",
        Species = ll$species,
        TaxonomyId = ll$taxonomyid,
        Coordinate_1_based = ll$coordinate_1_based,
        DataProvider = "",
        Maintainer = ll$maintainer,
        RDataClass = class(get(oname))[1],
        DispatchClass = "Rds",
        ResourceName = fname,
        RDataPath = paste0("curatedPCaData/", fname),
        Tags = x
      )

      tmp[c("SourceType", "SourceUrl", "SourceVersion", "DataProvider")] <- switch(study,
        abida = {
          c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_su2c_2019", SourceVersion = NA, DataProvider = "MSKCC")
        },
        baca = {
          c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_broad_2013", SourceVersion = NA, DataProvider = "Broad/Cornell")
        },
        barbieri = {
          c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_broad", SourceVersion = NA, DataProvider = "Broad/Cornell")
        },
        barwick = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18655", SourceVersion = NA, DataProvider = "Emory University")
        },
        chandran = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6919", SourceVersion = NA, DataProvider = "University of Pittsburgh")
        },
        friedrich = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134051", SourceVersion = NA, DataProvider = "Fraunhofer Institute")
        },
        hieronymus = {
          c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_mskcc_2014", SourceVersion = NA, DataProvider = "MSKCC")
        },
        icgcca = {
          c(SourceType = "tar.gz", SourceUrl = "https://dcc.icgc.org/releases/current/Projects/PRAD-CA", SourceVersion = "Release28", DataProvider = "ICGC")
        },
        igc = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109", SourceVersion = NA, DataProvider = "IGC")
        },
        kim = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119616", SourceVersion = NA, DataProvider = "GenomeDx")
        },
        kunderfranco = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14206", SourceVersion = NA, DataProvider = "Oncology Institute of Southern Switzerland")
        },
        ren = {
          c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_eururol_2017", SourceVersion = NA, DataProvider = "Shanghai Changhai Hospital")
        },
        sun = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136", SourceVersion = NA, DataProvider = "M. D. Anderson Cancer Center")
        },
        taylor = {
          c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21032,https://www.cbioportal.org/study/summary?id=prad_mskcc", SourceVersion = NA, DataProvider = "MSKCC")
        },
        tcga = {
          c(SourceType = "tar.gz", SourceUrl = "https://xenabrowser.net/datapages/?cohort=TCGA%20Prostate%20Cancer%20(PRAD)", SourceVersion = NA, DataProvider = "XenaBrowser")
        },
        true = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5132", SourceVersion = NA, DataProvider = "Fred Hutchinson Cancer Research Center")
        },
        wallace = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6956", SourceVersion = NA, DataProvider = "NCI")
        },
        wang = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8218", SourceVersion = NA, DataProvider = "University of California")
        },
        weiner = {
          c(SourceType = "tar.gz", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157548", SourceVersion = NA, DataProvider = "Decipher Bioscience")
        },
        # Default case
        {
          c(SourceType = NA, SourceUrl = NA, SourceVersion = NA, DataProvider = NA)
        }
      )

      metadata <- rbind(
        metadata,
        tmp
      )
    }

    # Save clinical info
    oname <- paste0(study, "_colData_", timestamp) # object name
    fname <- paste0(oname, ".Rds") # file name

    assign(oname, MultiAssayExperiment::colData(get(mae))) # Assign the object to the spesific name
    if (export) saveRDS(get(oname), file = paste0(path, "/", fname))

    # Add info to metadata
    descrp <- paste(oname, "Clinical metadata (colData-slot) of", study, "cohort in curatedPCaData package", sep = " ")

    tmp_clin <- c(
      Title = oname,
      Description = descrp,
      BiocVersion = ll$biocv,
      Genome = ll$genome,
      SourceType = "",
      SourceUrl = "",
      SourceVersion = "",
      Species = ll$species,
      TaxonomyId = ll$taxonomyid,
      Coordinate_1_based = ll$coordinate_1_based,
      DataProvider = "",
      Maintainer = ll$maintainer,
      RDataClass = class(get(oname))[1],
      DispatchClass = "Rds",
      ResourceName = fname,
      RDataPath = paste0("curatedPCaData/", fname),
      Tags = "clinical"
    )

    tmp_clin[c("SourceType", "SourceUrl", "SourceVersion", "DataProvider")] <- switch(study,
      abida = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_su2c_2019", SourceVersion = NA, DataProvider = "MSKCC")
      },
      baca = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_broad_2013", SourceVersion = NA, DataProvider = "Broad/Cornell")
      },
      barbieri = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_broad", SourceVersion = NA, DataProvider = "Broad/Cornell")
      },
      barwick = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18655", SourceVersion = NA, DataProvider = "Emory University")
      },
      chandran = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6919", SourceVersion = NA, DataProvider = "University of Pittsburgh")
      },
      friedrich = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134051", SourceVersion = NA, DataProvider = "Fraunhofer Institute")
      },
      hieronymus = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_mskcc_2014", SourceVersion = NA, DataProvider = "MSKCC")
      },
      icgcca = {
        c(SourceType = "TXT", SourceUrl = "https://dcc.icgc.org/releases/current/Projects/PRAD-CA", SourceVersion = "Release28", DataProvider = "ICGC")
      },
      igc = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109", SourceVersion = NA, DataProvider = "IGC")
      },
      kim = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119616", SourceVersion = NA, DataProvider = "GenomeDx")
      },
      kunderfranco = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14206", SourceVersion = NA, DataProvider = "Oncology Institute of Southern Switzerland")
      },
      ren = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_eururol_2017", SourceVersion = NA, DataProvider = "Shanghai Changhai Hospital")
      },
      sun = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136", SourceVersion = NA, DataProvider = "M. D. Anderson Cancer Center")
      },
      taylor = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21032,https://www.cbioportal.org/study/summary?id=prad_mskcc", SourceVersion = NA, DataProvider = "MSKCC")
      },
      tcga = {
        c(SourceType = "TXT", SourceUrl = "https://xenabrowser.net/datapages/?cohort=TCGA%20Prostate%20Cancer%20(PRAD)", SourceVersion = NA, DataProvider = "XenaBrowser")
      },
      true = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5132", SourceVersion = NA, DataProvider = "Fred Hutchinson Cancer Research Center")
      },
      wallace = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6956", SourceVersion = NA, DataProvider = "NCI")
      },
      wang = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8218", SourceVersion = NA, DataProvider = "University of California")
      },
      weiner = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157548", SourceVersion = NA, DataProvider = "Decipher Bioscience")
      },
      # Default case
      {
        c(SourceType = NA, SourceUrl = NA, SourceVersion = NA, DataProvider = NA)
      }
    )

    # Row bind
    metadata <- rbind(
      metadata,
      tmp_clin
    )

    # Save sampleMap
    oname <- paste0(study, "_sampleMap_", timestamp) # object name
    fname <- paste0(oname, ".Rds") # file name

    assign(oname, MultiAssayExperiment::sampleMap(get(mae))) # Assign the object to the spesific name
    if (export) saveRDS(get(oname), file = paste0(path, "/", fname))

    # Add info to metadata
    descrp <- paste(oname, "MAE-object sampleMap of", study, "cohort in curatedPCaData package", sep = " ")

    tmp_sampmap <- c(
      Title = oname,
      Description = descrp,
      BiocVersion = ll$biocv,
      Genome = ll$genome,
      SourceType = "",
      SourceUrl = "",
      SourceVersion = NA,
      Species = ll$species,
      TaxonomyId = ll$taxonomyid,
      Coordinate_1_based = ll$coordinate_1_based,
      DataProvider = "",
      Maintainer = ll$maintainer,
      RDataClass = class(get(oname))[1],
      DispatchClass = "Rds",
      ResourceName = fname,
      RDataPath = paste0("curatedPCaData/", fname),
      Tags = "sampleMap"
    )

    tmp_sampmap[c("SourceType", "SourceUrl", "SourceVersion", "DataProvider")] <- switch(study,
      abida = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_su2c_2019", SourceVersion = NA, DataProvider = "MSKCC")
      },
      baca = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_broad_2013", SourceVersion = NA, DataProvider = "Broad/Cornell")
      },
      barbieri = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_broad", SourceVersion = NA, DataProvider = "Broad/Cornell")
      },
      barwick = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18655", SourceVersion = NA, DataProvider = "Emory University")
      },
      chandran = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6919", SourceVersion = NA, DataProvider = "University of Pittsburgh")
      },
      friedrich = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134051", SourceVersion = NA, DataProvider = "Fraunhofer Institute")
      },
      hieronymus = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_mskcc_2014", SourceVersion = NA, DataProvider = "MSKCC")
      },
      icgcca = {
        c(SourceType = "TXT", SourceUrl = "https://dcc.icgc.org/releases/current/Projects/PRAD-CA", SourceVersion = "Release28", DataProvider = "ICGC")
      },
      igc = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2109", SourceVersion = NA, DataProvider = "IGC")
      },
      kim = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119616", SourceVersion = NA, DataProvider = "GenomeDx")
      },
      kunderfranco = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14206", SourceVersion = NA, DataProvider = "Oncology Institute of Southern Switzerland")
      },
      ren = {
        c(SourceType = "TXT", SourceUrl = "https://www.cbioportal.org/study/summary?id=prad_eururol_2017", SourceVersion = NA, DataProvider = "Shanghai Changhai Hospital")
      },
      sun = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25136", SourceVersion = NA, DataProvider = "M. D. Anderson Cancer Center")
      },
      taylor = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE21032,https://www.cbioportal.org/study/summary?id=prad_mskcc", SourceVersion = NA, DataProvider = "MSKCC")
      },
      tcga = {
        c(SourceType = "TXT", SourceUrl = "https://xenabrowser.net/datapages/?cohort=TCGA%20Prostate%20Cancer%20(PRAD)", SourceVersion = NA, DataProvider = "XenaBrowser")
      },
      true = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5132", SourceVersion = NA, DataProvider = "Fred Hutchinson Cancer Research Center")
      },
      wallace = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6956", SourceVersion = NA, DataProvider = "NCI")
      },
      wang = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8218", SourceVersion = NA, DataProvider = "University of California")
      },
      weiner = {
        c(SourceType = "TXT", SourceUrl = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157548", SourceVersion = NA, DataProvider = "Decipher Bioscience")
      },
      # Default case
      {
        c(SourceType = NA, SourceUrl = NA, SourceVersion = NA, DataProvider = NA)
      }
    )

    # Row bind
    metadata <- rbind(
      metadata,
      tmp_sampmap
    )
  }
  # Set column names properly
  colnames(metadata) <- c("Title", "Description", "BiocVersion", "Genome", "SourceType", "SourceUrl", "SourceVersion", "Species", "TaxonomyId", "Coordinate_1_based", "DataProvider", "Maintainer", "RDataClass", "DispatchClass", "ResourceName", "RDataPath", "Tags")

  metadata
}

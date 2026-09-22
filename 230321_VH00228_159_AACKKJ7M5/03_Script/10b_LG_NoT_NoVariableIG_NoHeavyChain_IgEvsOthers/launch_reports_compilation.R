# ####################################################################
# This script launch the compilation of both reports (one per sample)
# and rename them accordingly with the sample name
# ####################################################################

library( knitr)
library( rmarkdown)
library( funr)


# Boolean to activate overwriting security (change at your own risk !)
preventOverwrite = FALSE;

# Boolean to set whether a multi-thread proof strategy must be used (same report
# '.Rmd' name rendered in different folder/context mixes up output).
# Could result in "tempRmdCopy_*.Rmd" files to appear in the 'Rmd' location, and
# require a manual cleanup if the rendering fails, and the script is not able to 
# clean it on exit automatically.
multirenderSafe = TRUE


### Define working folder (contains R/Rmd file for current sample, parent contains global project files)
library( funr)
if( exists( "snakemake"))
{
  cat("\nWorking in Snakemake mode\n")
  WORKING_DIR = snakemake@scriptdir
}else
{
  cat("\nWorking in local script mode\n")
  WORKING_DIR = tryCatch(
    {
      dirname( sys.script())
    },
    error=function(cond) 
    {
      cat("\nWorking in local R session mode\n")
      return( getwd())
    }
  )
}

### Load parameters
# Define an environment that will contain parameters
paramsEnv = new.env();

# Assign the WORKING_DIR to the paramsEnv
assign( "WORKING_DIR" , WORKING_DIR , pos = paramsEnv )

# Load file defining global parameters
globalParamsFilePath = file.path( WORKING_DIR, "../globalParams.R");
if(file.exists(globalParamsFilePath)) 
{
  source( globalParamsFilePath, local = paramsEnv);
} else 
{
  warning("The file 'globalParamsFilePath.R' containing global parameters is missing.");
}

# Load file defining sample parameters
sampleParamsFilePath = file.path( WORKING_DIR, "../sampleParams.R");
if(file.exists(sampleParamsFilePath))
{
  source( sampleParamsFilePath, local = paramsEnv);
} else
{
  warning("The file 'sampleParamsFilePath.R' containing sample parameters is missing.");
}

# Load file defining analysis parameters
analysisParamsFilePath = file.path( WORKING_DIR, "analysisParams.R");
if(file.exists(analysisParamsFilePath)) 
{
  source( analysisParamsFilePath, local = paramsEnv);
} else 
{
  warning("The file 'analysisParamsFilePath.R' containing analysis-specific parameters is missing.");
}

### Compile the HTML report

# Record the original output folder as specified in loaded parameters, for an
# eventual alteration before launch (e.g. generate several reports in subfolders
# to try multiple values for a parameter)
ORIGINAL_PATH_ANALYSIS_OUTPUT = paramsEnv[["PATH_ANALYSIS_OUTPUT"]]


#FINDCLUSTERS_RESOLUTION=0.8
# Eventual loop on a parameter to be tested (override value before rendering)
#for( FINDCLUSTERS_RESOLUTION in c( 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))
{
#  # Assign the new value (and output path) to the environment used for rendering
#  assign( "FINDCLUSTERS_RESOLUTION" , 
#          FINDCLUSTERS_RESOLUTION , 
#          pos = paramsEnv )
#  assign( "PATH_ANALYSIS_OUTPUT" , 
#          file.path( ORIGINAL_PATH_ANALYSIS_OUTPUT, paste0("ClusteringRes_", FINDCLUSTERS_RESOLUTION)), 
#          pos = paramsEnv )
  
  # Create output file name and check if already exists to prevent overwriting
  reportOutputFilename = paste0( paramsEnv[["SCIENTIFIC_PROJECT_NAME"]], "_",
                                 paramsEnv[["EXPERIMENT_PROJECT_NAME"]], "_",
                                 paramsEnv[["ANALYSIS_STEP_NAME"]], ".html");
  
  alreadyExists = file.exists( file.path( paramsEnv[["PATH_ANALYSIS_OUTPUT"]], reportOutputFilename));
  if( alreadyExists && preventOverwrite) stop( paste( "Report file already exists:", reportOutputFilename));
   
  ### Prevent result mixing when parallel rendering
  # Create a copy of original rmd file with unique name.
  # Prevents 'render' to overwrite 'md' file created from name of 'Rmd' file,
  # and scrambled results when parallel-rendering 'Rmd' files with identical name.
  # Specifying 'tmp' folders for 'intermediates_dir' and 'knit_root_dir' prevent
  # usage of WORKING_DIR in reports so we must go for this alternative.
  # Final output 'html' name provided to render anyway (not based on 'Rmd' name).
  rmdCopyName=NULL
  if(multirenderSafe)
  {
    # Create a unique temporary file name
    rmdCopyName = tempfile( pattern = "tempRmdCopy_", tmpdir = WORKING_DIR, fileext = ".Rmd");
    # Copy 'Rmd' file to it
    stopifnot( all( file.copy( from = file.path( WORKING_DIR, "Report.Rmd"),
                               to   = rmdCopyName)));
  }
  
  ?rlang::env_clone
  
  #cat("\ninput:", rmdCopyName)
  #cat("\noutput_dir:", paramsEnv[["PATH_ANALYSIS_OUTPUT"]])
  #cat("\noutput_file:", reportOutputFilename)
  #cat("\n")
  
  rmarkdown::render( input = if(multirenderSafe) rmdCopyName else file.path( WORKING_DIR, "Report.Rmd"),
                     output_dir = paramsEnv[["PATH_ANALYSIS_OUTPUT"]],
                     output_file  = I( reportOutputFilename),
                     envir = rlang::env_clone(paramsEnv), # Clone so the second runs don't get values from previous ones 
                     quiet = FALSE)
  
  # Remove temporary 'Rmd' copy
  if(multirenderSafe)
  {
    if(! file.remove( rmdCopyName)) warning( paste0( "Temporary file '", rmdCopyName, "' could not be removed..."));
  }
}



# TRIP.workflow.R -- an explanatory example workflow

TRIP.PATH <- "~/TRIP_Screen"
# load the scripts
source( file.path( TRIP.PATH, "mutantAligner.R"))

# don't use too many cores on Lenovo to save memory (as each sample will use multiple cores under the hood)
# (if having problems, debugging errors, etc., set to 1)
multicore.setup(4)

# run the full pipeline:  aligning;  resolving mate pairs;  pie charts;  then lastly model the growth defects
do.all( "SampleKey.Example.txt", sharedDayZero="NoDrug")

# the above almost always fails a few times as all issues get fixed.
# Once it works, you will get a plot of QC overview called "AlignmentSuccessOverview.pdf" that shows number of
# genes and the percentage of mate pair reads that were valid for each sample.  
# Samples should cluster in the upper right of the plot.  
# Those landing in the lower left should get flagged for "Exclude" in the sample key file.

# Once all looks good, do the final modeling
model.All.TRIP.Samples( "SampleKey.Example.txt", sharedDayZero="NoDrug")

# all final result files (both .TXT and .HTML) will be in the "Results/Example/" subfolder

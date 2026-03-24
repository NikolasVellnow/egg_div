
# Read the parameters and fst data handed over in the shell

args <- commandArgs(trailingOnly = TRUE)

global_file  <- args[1]
window_file  <- args[2]
output_file  <- args[3]
pop1_name    <- args[4]
pop2_name    <- args[5]
region       <- args[6]
window_size  <- as.numeric(args[7])
step_size    <- as.numeric(args[8])

print(args)

# Read the header line manually
header_line <- readLines(window_file, n = 1)

# Split into column names
header_names <- strsplit(header_line, "\t")[[1]]

# Add missing column name
header_names <- c(header_names, "fst")

# Read data without header and without auto row names
data <- read.delim(
  window_file,
  header = FALSE,
  skip = 1,
  sep = "\t",
  row.names = NULL,
  stringsAsFactors = FALSE
)


title_str <- paste0(
  pop1_name, " vs. ", pop2_name,
  " (window size = ", window_size/1000, " kb, step = ",
  step_size/1000, " kb)"
)

global_fst <- scan(global_file, quiet=T)[2]


# Assign correct column names
colnames(data) <- header_names

# write function to categorize window
is_sex_diff <- function(fst_value, global_fst_value){
  if (is.na(fst_value)) {
    NA
  } else if (fst_value > 2*global_fst_value){
    "sex_diff"
  } else {
    "not_sex_diff"
  }
} 

# vectorize function
v_is_sex_diff <- Vectorize(is_sex_diff)

data$sex_diff<- factor(v_is_sex_diff(data$fst,global_fst) )

max_fst <- max(data$fst)
chr_name <- region

pdf(file = output_file)
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0,1,0,0))


plot(data$midPos, data$fst,
     ylim = c(0, max_fst+0.05),
     #xlim = c(min(chr_data$midPos), max(chr_data$midPos)),
     cex.axis=1.0,
     cex.lab=1.4,
     type = "b",
     pch = 18,
     col=data$sex_diff,
     ylab = expression('F'['ST']), xlab = paste("Position on chromosome", chr_name),
     main = paste(title_str))

abline(h=global_fst, col="red")

dev.off()


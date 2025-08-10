\name{gene_pair_matrix}
\alias{gene_pair_matrix}
\title{
Construct Gene Pair Matrices Using Parallel Processing
}
\description{
This function constructs a matrix of gene pair comparisons using parallel processing.
It takes a dataset and a combination matrix specifying gene pairs to compare, and applies the gene_pair function to each pair in parallel, using multiple cores.
}
\usage{
gene_pair_matrix(dat, combination, ncore)
}
\arguments{
  \item{dat}{
A data frame containing gene expression data, where each column represents a gene and each row represents a sample.
}
 \item{combination}{
A matrix where each row contains a pair of gene names (as column names found in dat) to be compared.
}
  \item{ncore}{
The number of cores to use for parallel processing.
}
}
\details{
The function uses the \code{makeCluster} and \code{parLapply} functions from the \code{parallel} package to distribute the computation of gene pair comparisons across multiple cores. After computation, it combines the results into a single matrix and returns it.
}
\value{
A matrix where each column represents the comparison results of a gene pair across all samples.
}
\references{
% Add references to the literature or websites here if necessary
}
\author{
Bicheng Ye
}
\note{
This function is designed for efficient computation of gene pair comparisons, particularly useful in genomic studies where large datasets are common.
}

\seealso{
% Use \code{\link{other_function}} to link to related functions or documentation
% Example: \code{\link{gene_pair}}
}
\examples{
# Set a random seed for reproducibility
set.seed(123)

# Create a sample dataset with gene expression data
example_data <- data.frame(
  GeneA = rnorm(10, mean = 5, sd = 2), # Randomly generated data for GeneA
  GeneB = rnorm(10, mean = 3, sd = 1.5), # Randomly generated data for GeneB
  GeneC = rnorm(10, mean = 2, sd = 1), # Randomly generated data for GeneC
  GeneD = rnorm(10, mean = 4, sd = 2) # Randomly generated data for GeneD
)

# Name the rows of the dataframe to represent samples
rownames(example_data) <- paste("Sample", 1:10)

# Define a combination of gene pairs to compare
gene_pairs <- combn(colnames(example_data), 2)

# Use the gene_pair_matrix function to compare the expression values of the gene pairs
# Assume we are using 4 cores for parallel processing
result_matrix <- gene_pair_matrix(dat = example_data, combination = gene_pairs, ncore = 4)

# Print the resulting matrix
print(result_matrix)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{methods}

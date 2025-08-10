\name{gene_pair}
\alias{gene_pair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Construct Gene Pairs
%%  ~~function to do ... ~~
}
\description{
The function is used to extract the expression values of two specified genes from a given dataset, compare the expression levels of these two genes, and construct a comparison matrix for the gene pairs.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gene_pair(dat, gene1, gene2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{A data frame containing gene expression data, where each column represents a gene and each row represents a sample.
%%     ~~Describe \code{dat} here~~
}
  \item{gene1}{The name of the first gene.
%%     ~~Describe \code{gene1} here~~
}
  \item{gene2}{The name of the second gene.
%%     ~~Describe \code{gene2} here~~
}
}
\details{
The function will compare the expression values of the two genes in each sample and generate a binary matrix, where 1 indicates that the expression level of gene1 is higher than gene2, and 0 indicates the opposite.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix containing the comparison results of the two genes in each sample.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Bicheng Ye
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Set a random seed for reproducibility
set.seed(123)
# Create a sample dataset with gene expression data
example_data <- data.frame(
  GeneA = rnorm(3, mean = 5, sd = 2), # Randomly generated data for GeneA
  GeneB = rnorm(3, mean = 3, sd = 1.5), # Randomly generated data for GeneB
  GeneC = rnorm(3, mean = 2, sd = 1), # Randomly generated data for GeneC
  GeneD = rnorm(3, mean = 4, sd = 2) # Randomly generated data for GeneD
)

# Name the rows of the dataframe to represent samples
rownames(example_data) <- paste("Sample", 1:3)

# Use the gene_pair function to compare the expression values of GeneA and GeneB
result_matrix <- gene_pair(example_data, "GeneA", "GeneB")

# Print the resulting matrix
print(result_matrix)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

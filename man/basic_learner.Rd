\name{basic_learner}
\alias{basic_learner}
\title{
Construct Basic Learners Using Various Machine Learning Methods
}
\description{
This function constructs predictive models using various machine learning methods, including Elastic Net Regression (ERN), Random Forest (RF), Support Vector Machines (SVM) with different kernels, and K-Nearest Neighbors (KNN). It utilizes parallel processing for efficiency.
}
\usage{
basic_learner(trainingset, CVnumber, Cvrepeats, ncore)
}
\arguments{
  \item{trainingset}{
A data frame containing the training data. The first column should be the response variable, and the remaining columns should be predictor variables.
}
  \item{CVnumber}{
An integer specifying the number of folds to use in cross-validation.
}
  \item{Cvrepeats}{
An integer specifying the number of times to repeat the cross-validation.
}
  \item{ncore}{
An integer specifying the number of cores to use for parallel processing.
}
}
\details{
The function uses the \code{trainControl} function from the \code{caret} package to set up repeated cross-validation with SMOTE sampling. It fits multiple models including Elastic Net Regression, Random Forest, Support Vector Machines (with linear, radial, and polynomial kernels), and K-Nearest Neighbors. Results are stored in lists for each model type.
}
\value{
A list containing fitted models for each machine learning method:
\item{ern_fit}{A list of fitted Elastic Net models with different alpha values.}
\item{KNN_fit}{A list of fitted KNN models with different k values and kernels.}
\item{rf_fit}{A list of fitted Random Forest models with different depths and number of trees.}
\item{svmLinear_fit}{A list of fitted SVM linear models with different C parameters.}
\item{svmPoly_fit}{A list of fitted SVM polynomial models with different C parameters.}
\item{svmRadial_fit}{A list of fitted SVM radial models with different sigma and C parameters.}
}
\references{
Refer to the documentation for the \code{caret} package for more information on model training and evaluation.
}
\author{
Bicheng Ye
}
\note{
This function is designed for users who want to quickly build and evaluate multiple machine learning models on a given dataset using parallel processing for efficiency.
}

\seealso{
\code{\link{caret}} for model training functions.
}

\examples{
# Example usage:
# Set a random seed for reproducibility
set.seed(123)

# Create a sample dataset
example_data <- data.frame(
  res = factor(sample(c(0, 1), 100, replace = TRUE)), # Response variable
  var1 = rnorm(100),
  var2 = rnorm(100),
  var3 = rnorm(100)
)

# Run basic_learner function
results <- basic_learner(trainingset = example_data, CVnumber = 5, Cvrepeats = 3, ncore = 4)

# Print results
print(results)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ ~machine learning }

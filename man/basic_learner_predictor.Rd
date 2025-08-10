\name{basic_learner_predictor}
\alias{basic_learner_predictor}
\title{
Predict Outcomes Using Trained Learners
}
\description{
This function generates predictions for a given test dataset using a list of trained machine learning models, including Elastic Net Regression, K-Nearest Neighbors, Random Forest, and Support Vector Machines.
}
\usage{
basic_learner_predictor(set, learnerlist)
}
\arguments{
  \item{set}{
A data frame containing the test data. The first column should be the response variable (though it may not be used directly for prediction), and the remaining columns should be predictor variables.
}
  \item{learnerlist}{
A list containing the fitted models from the training process, including models for Elastic Net Regression, KNN, Random Forest, and SVM.
}
}
\details{
The function uses the \code{predict} method to generate probabilities or class predictions from the specified models in \code{learnerlist}. It combines the predictions into a single data frame with row names corresponding to the test set samples.
}
\value{
A data frame containing predicted scores for each model:
\item{ern_score}{Predicted probabilities from Elastic Net Regression.}
\item{KNN_score}{Predicted probabilities from K-Nearest Neighbors.}
\item{rf_score}{Predicted probabilities from Random Forest.}
\item{svmLinear_score}{Predicted class labels from SVM with linear kernel.}
\item{svmPoly_score}{Predicted class labels from SVM with polynomial kernel.}
\item{svmRadial_score}{Predicted class labels from SVM with radial kernel.}
}
\references{
Refer to the documentation for the \code{caret} package and other relevant literature on machine learning methods used.
}
\author{
Bicheng Ye
}
\note{
This function is designed to facilitate predictions using multiple machine learning models trained on a given training dataset. It leverages parallel processing capabilities for efficiency.
}

\seealso{
\code{\link{basic_learner}} for training models.
}

\examples{
# Example usage:
# Set a random seed for reproducibility
set.seed(123)

# Create a sample test dataset
test_data <- data.frame(
  var1 = rnorm(10),
  var2 = rnorm(10),
  var3 = rnorm(10)
)

# Assume 'learnerlist' is already created from basic_learner function
# results <- basic_learner_predictor(set = test_data, learnerlist = learnerlist)

# Print predicted scores
# print(results)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ ~prediction }

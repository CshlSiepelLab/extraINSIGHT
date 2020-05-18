library(data.table)
library(gtools)

stripGlmLR = function(cm) {
  cm$y = c()
  cm$model = c()
  cm$residuals = c()
  cm$fitted.values = c()
  cm$effects = c()
  cm$qr$qr = c() 
  cm$linear.predictors = c()
  cm$weights = c()
  cm$prior.weights = c()
  cm$data = c()
  cm$family$variance = c()
  cm$family$dev.resids = c()
  cm$family$aic = c()
  cm$family$validmu = c()
  cm$family$simulate = c()
  attr(cm$terms,".Environment") = c()
  attr(cm$formula,".Environment") = c()

  return(cm)
}

args = commandArgs(trailingOnly=TRUE)
## Read in data
d <- fread(paste("zcat", args[1]), sep="\t")
## Rename columns for greater script readibility
setnames(d, colnames(d), c("chrom","start","end","context_mutation","mutation_flag","coverage","logit_mutation_rate","gc_content","cpg_island"))
## Fit global glm model
model <- glm(mutation_flag ~log(coverage) + logit_mutation_rate + gc_content + cpg_island, family=binomial, data=d, model=FALSE, x=FALSE, y=FALSE)
## Save full model
save(model, file=args[2])
## Save stripped down model with basically everything except the coefficients removed
model <- stripGlmLR(model)
save(model, file=args[3])

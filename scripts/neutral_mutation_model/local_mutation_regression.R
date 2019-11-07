library(data.table)
args = commandArgs(trailingOnly=TRUE)

# load global model
load(args[1])
# 
min_n_mutation = as.numeric(args[2])
# Read in local data from stdin
tryCatch(d<-fread('cat /dev/stdin'), error=function(e) NULL)
# print(head(d))

# if (exists("d")){
if (exists("d") && sum(d$V5 == 1 & d$V11 != '.') >= min_n_mutation){
    d$V14 <- predict(model, newdata=d)
    d_neut <- d[d$V11 != ".", ]

    local_model <- glm(V5~V14, data=d_neut, family="binomial")
    d$V15 <- predict(local_model, newdata=d, type="response")
    write.table(d[, c(1:6, 15)], stdout(), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
}

data <- demo.clustered
head(data)
tail(data)

mock.data <- data[data$Group == 'Mock', ]
unique.sample <- unique(mock.data$Sample)

train.data <- data.frame(matrix(ncol = ncol(mock.data), nrow = 0))
unlabelled.data <- data.frame(matrix(ncol = ncol(mock.data), nrow = 0))
colnames(train.data) <- colnames(mock.data)
colnames(unlabelled.data) <- colnames(mock.data)

for (samp.name in unique.sample) {
  mock.data.sample <- mock.data[mock.data$Sample == samp.name, ]
  meta.cluster <- unique(mock.data.sample$FlowSOM_metacluster)
  
  for (m in c(1:5)) {
    mock.data.cluster <- mock.data.sample[mock.data.sample$FlowSOM_metacluster == m, ]
    split.data <- split(mock.data.cluster, sample(rep(1:2, )))
    
    train.data <- rbind(train.data, split.data[["1"]])
    unlabelled.data <- rbind(unlabelled.data, split.data[["2"]])
    
    colnames(train.data) <- colnames(mock.data)
    colnames(unlabelled.data) <- colnames(mock.data)
    
  }
  
}

split.training <- split(train.data, sample(rep(1:2,)))
write.csv(split.training[["1"]], file = "Training1.csv", row.names = FALSE)
write.csv(split.training[["2"]], file = "Training2.csv", row.names = FALSE)

split.training <- split(unlabelled.data, sample(rep(1:2,)))
write.csv(split.training[["1"]], file = "Labelled1.csv", row.names = FALSE)
write.csv(split.training[["2"]], file = "Labelled2.csv", row.names = FALSE)

unlabelled.data.no.label <- unlabelled.data[, c(1:37)]
split.training <- split(unlabelled.data.no.label, sample(rep(1:2,)))
write.csv(split.training[["1"]], file = "Unlabelled1.csv", row.names = FALSE)
write.csv(split.training[["2"]], file = "Unlabelled2.csv", row.names = FALSE)

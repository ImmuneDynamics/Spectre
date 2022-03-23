library(ggplot2)

# create a dataset
specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
value <- abs(rnorm(12 , 0 , 15))
data <- data.frame(specie,condition,value)

data

# Stacked + percent
dev.off()
ggplot(data, aes(fill=condition, y=value, x=specie)) + 
  geom_bar(position="fill", stat="identity")


new.list

new.list$Group <- c(rep('Mock', 288), rep('WNV', 288))
new.list



dev.off()
ggplot(new.list, aes(fill=Boolean, y=value, x=Group)) + 
  geom_bar(position="fill", stat="identity")

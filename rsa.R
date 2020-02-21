source('final.R')

# back=0 No protection, back=1 Link disjoint, back=2 Node & link disjoint

result = sacred(adjin = 'Adjacency_list.txt',demands = 'Demands.txt',back = 0)

# The Adjacency list and Demands text files must be in the local storage (for Rstudio cloud it is the Files section)

View(result)
View(result$results)
View(result$spectrum)
View(result$blrate)
par(mfrow=c(1,2))
plot(result$blrate$V2, type="l", main="Blocking Rate", xlab="Demands", ylab="Blocking Rate")
plot(result$blrate$V3, type="l", main="Slot Utilisation", xlab="Demands", ylab="Slots Utilised")
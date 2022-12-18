#### Libraries ####

library(argparser)
library(dplyr)
library(tidyr)


######## Arguments ##########
p <- arg_parser("change the format of nucleotide table (melt)")
p <- add_argument(p, "-i", help="input")
# p <- add_argument(p, "-o", help="output")

argv <- parse_args(p)

nuc_table <- read.table(argv$i, header = TRUE)

removeFirstCharacter <- function(x) {
  return(as.numeric(substr(x, 2, nchar(x))))
}

remove_column <- function(data, col_name) {
data <- data[, !names(data) %in% col_name]
return(data)
}

dt_organized <- nuc_table %>% gather(Position, count, 2:ncol(nuc_table))
colnames(dt_organized) <- c("nucleotide", "position", "count")
summary <- dt_organized %>% group_by(position) %>% summarize("total"=sum(count))
d <- merge(dt_organized,summary,by=c("position"))
d$freq <- 100*d$count/d$total
d$position <- sapply(d$position, removeFirstCharacter)
d$readLength <- max(d$position)
# d <- remove_column(d, "count")
d <- remove_column(d, "total")
write.table(d, file="", sep='\t',  quote = F, row.names=F, col.names=F)

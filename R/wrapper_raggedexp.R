wrapper_raggedexp<-function(ragexp){
I=sparseAssay(ragexp,"Variant_Classification")
I=as.data.frame(I)
I$gene=rownames(I)
#I=I[,c(44,1:43)]
I$gene=gsub("\\..*","",I$gene)


v=I %>% 
  dplyr::group_by(gene) %>% 
  dplyr::summarise_all(toString)

v=as.data.frame(v)

b <- data.frame(lapply(v, function(x) {
  gsub("NA, NA, NA", "NA", x)
}))
D<-data.frame(lapply(b, function(x) {
  gsub("NA, NA", "NA", x)
}))
E<-data.frame(lapply(D, function(x) {
  gsub("NA,", "", x)
}))
G<-data.frame(lapply(E, function(x) {
  gsub(", NA", "", x)
}))

rownames(G)=G$gene
G <- G[ , ! names(G) %in% "gene"]

return(G)
}
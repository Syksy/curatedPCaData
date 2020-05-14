
sink(file="all_phenotypes.txt", append=T)
for (i in zz[3:35]){
	
	dd = read.csv(i, stringsAsFactors=F)
	
	capture.output(i, file="all_pheno.txt", append=T)
    capture.output(apply(dd, 2, table), file="all_pheno.txt", append=T)
	
}


##Age of genes
ageDroso <- read.table("ageDroso.txt", sep="\t")
colnames(ageDroso) <- c("Ensembl.Gene.ID", "clade")
levels(ageDroso$clade)
#levels : Bilateria, Drosophila melanogaster, Ecdysozoa, Opisthokonta

ageXeno <- read.table("ageXeno.txt", sep="\t")
colnames(ageXeno) <- c("Ensembl.Gene.ID", "clade")
levels(ageXeno$clade)
# [1] "Bilateria"          "Chordata"           "Euteleostomi"       
# "Opisthokonta"       "Sarcopterygii"     
# [6] "Tetrapoda"          "Vertebrata"         "Xenopus tropicalis"

ageCaeno <- read.table("ageCaeno.txt", sep="\t")
colnames(ageCaeno) <- c("Ensembl.Gene.ID", "clade")
levels(ageCaeno$clade)
# [1] "Bilateria"              "Caenorhabditis elegans" "Ecdysozoa"              "Opisthokonta"    

ageMus <- read.table("ageMus.txt", sep="\t")
colnames(ageMus) <- c("Ensembl.Gene.ID", "clade")
levels(ageMus$clade)
# [1] "Amniota"          "Bilateria"        "Boreoeutheria"    "Chordata"  
# "Euarchontoglires" "Euteleostomi"    
# [7] "Eutheria"         "Glires"           "Mammalia"         "Murinae"       
# "Mus musculus"     "Opisthokonta"    
# [13] "Rodentia"         "Sarcopterygii"    "Sciurognathi"     "Tetrapoda"    
# "Theria"           "Vertebrata" 

ageDanio <- read.table("ageDanio.txt", sep="\t")
colnames(ageDanio) <- c("Ensembl.Gene.ID", "clade")
levels(ageDanio$clade)
# [1] "Bilateria"     "Chordata"      "Clupeocephala" "Danio rerio"   
# "Euteleostomi"  "Neopterygii"   "Opisthokonta" 
# [8] "Otophysa"      "Vertebrata" 

reuni <- function(z) Reduce('union', z)

clades <- reuni(list(levels(ageDroso$clade), levels(ageCaeno$clade),levels(ageXeno$clade),
      levels(ageDanio$clade), levels(ageMus$clade))) #length 26
      

          # clades2 <- union(levels(ageDroso$clade), levels(ageCaeno$clade))
          # clades2 <- union(clades2, levels(ageXeno$clade))
          # clades2 <- union(clades2, levels(ageMus$clade))
          # clades2 <- union(clades2, levels(ageDanio$clade)) #length 26

          #clades2 %in% clades => all TRUE



# [1] "Bilateria"       936             "Drosophila melanogaster" 
# "Ecdysozoa"           936             "Opisthokonta"    1215          
# [5] "Caenorhabditis elegans"          "Chordata"        722                
# "Euteleostomi"        441             "Sarcopterygii"   414       
# [9] "Tetrapoda"       371             "Vertebrata"      535
# "Xenopus tropicalis"                  "Clupeocephala"   265       
# [13] "Danio rerio"                    "Neopterygii"     333 
# "Otophysa"            152              "Amniota"         296                
# [17] "Boreoeutheria"  100             "Euarchontoglires"  92    
# "Eutheria"            104              "Glires"           86       
# [21] "Mammalia"       167             "Murinae"           25      "
# Mus musculus"                          "Rodentia"         77      
# [25] "Sciurognathi"   74               "Theria"           162

age <- c(937, 1, 936, 1215, 1, 722, 441, 414, 371, 535, 1, 265, 1, 333, 152, 296,
         100, 92, 104, 86, 167, 25, 1, 77, 74, 162)

ageTable <- data.frame(Taxon=clades, Age=age)

write.table(ageTable, file="/home/user/Bureau/stage15/ageTable.txt")


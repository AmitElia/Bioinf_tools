rm(list = ls()) #clear environment

print("Hello Bioinformatics")

#Sources: https://stackoverflow.com/questions/21263636/read-fasta-into-a-dataframe-and-extract-subsequences-of-fasta-file
#https://www.biostars.org/p/274312
#https://www.edureka.co/community/2091/how-to-import-text-file-as-a-single-character-string

library(Biostrings)
library(dplyr)
library(stringr)

wdir = "/Users/Amit Elia/Documents/BENG182/datasets/"

#4)

#Read txt file into XStringSet format.
fastaF <- readAAStringSet(paste0(wdir, "datafile.txt"), format = "fasta")

#Converts x to a character vector of the same length as x. 
#The use.names argument controls whether or not names(x) should be propagated to the names of the returned vector.
sequences <- as.character(fastaF, use.names= FALSE)
sequences

# description or comment for each element in x
seq_names = names(fastaF)
seq_names

# A vector of non-negative integers containing the number of letters for each element in x.
seq_lengths <- width(fastaF)
seq_lengths


#creating dataframe with output column
df <- data.frame(seq_names, sequences, seq_lengths)
df["Out"] <- paste(df$seq_names, df$seq_lengths)
View(df)

#print output
df$Out

#5)

df_mouse_rat <- filter(df, grepl('Mus musculus|Rattus norvegicus|_RAT', seq_names))
View(df_mouse_rat)

pattern_list <- data.frame(name = c(df_mouse_rat$seq_names))
pattern_list
fasta_mouse_rat <- fastaF[pattern_list$name]

writeXStringSet(fasta_mouse_rat, "~/BENG182/fasta_mouse_rat.fasta")

#6)

#The first file data.seq contains the concatenation of all of the sequences from each file with no headers, and no newline symbols.
writeLines(sequences, "~/BENG182/data.seq", sep = "@")
df$gi <- str_extract(df$seq_names, "\\|.*?\\|")
df$gi <- str_remove_all(df$gi, "[\\|\\|]")

index = 0
offsets <- c(0)
for (x in seq_lengths){
  offsets <- c(offsets, x + 1 + offsets[length(offsets)])
  index <- index + 1
}
offsets <- head(offsets, -1)
offsets

df_data_in <- data.frame(df$gi, offsets)
df_data_in
write.table(df_data_in, file= "~/BENG182/data.in")

#7)

getSeqFunc <- function(query) {
  seq_data <- readChar("~/BENG182/data.seq", file.info("~/BENG182/data.seq")$size)
  seq_df <- read.table("~/BENG182/data.in")
  View(seq_data)
  View(seq_df)
  locations <- gregexpr(pattern = query, seq_data)
  View(locations)
  
  output <- c()
  idx = 1
  for(x in head(seq_df$offsets, -1)){
    for(loc in locations){
      #print(c(seq_df$offsets[idx], seq_df$offsets[idx+1], loc))
      
      if(seq_df$offsets[idx]< loc & seq_df$offsets[idx+1] > loc){
        output <- c(output, seq_df$df.gi[idx])
      }
    }
    idx <- idx +1
  }
  return(output)
  
}

output <- getSeqFunc("MHIQITDFGTAKVLSPDS")
output

write.table(output, file= "~/BENG182/output.csv")

# To install: devtools::install_github("nacnudus/tidyxl")
# library(tidyxl)
inputFile = "RawData/AxonOutgrowth/2018-10-26_OutgrowthAssay.xlsx"

# Get a denormalized table of cells 
cells <- xlsx_cells(inputFile)

# Get formats with addresses
formats = xlsx_formats(inputFile)

# Parse out colours; make a factored version
COLOURS = unique(formats$local$fill$patternFill$fgColor$rgb[!is.na(formats$local$fill$patternFill$bgColor$rgb)])
COLOURS.factor <- as.factor(COLOURS)

# Add day by colour
cells$Day <- -42 # Make a constant? default day ID.
lapply(1:length(COLOURS), 
       FUN=function(i){
         
         cells$Day[cells$local_format_id %in%
                      which(formats$local$fill$patternFill$fgColor$rgb == COLOURS[i])
                    ] <<- COLOURS.factor[i]
       })

for (day in unique(cells$Day)){
  if (day == -42) { next }
  headers = cells[cells$Day == day, ]
  cells[cells$col %in% headers$col, "Day"] <- day 
}

# Identify variants
variant.columns.character <- unique(
  cells$character[ !is.na(cells$character) ]   
)

variant.characters <- cells$character[cells$character %in% variant.columns.character]
variant.columns <- cells$col[cells$character %in% variant.columns.character]
variant.rows <- cells$row[cells$character %in% variant.columns.character]
variant.batchs <- cells$Day[cells$character %in% variant.columns.character]

# variant.characters[variant.characters %in% c("GFP-9", "T167N")]
# variant.batchs[variant.characters %in% c("GFP-9", "T167N")]

axon.tidy <- do.call(rbind,
                    lapply(1:length(variant.columns),
                           FUN = function(IDX){
                             print(IDX)
                             xcol = variant.columns[IDX]
                             xrow = variant.rows[IDX]
                             variant = variant.characters[IDX]
                             
                             # Variants
                             variants = cells[
                                 cells$col == xcol &
                                 cells$row > xrow &
                                 !is.na(cells$numeric), c("numeric", "Day")
                               ]
                             #variant.molten = melt(table(variant.vialcounts))
                             plates = variants$Day[1:(length(variants$Day))]
                             plates[T] <- -42
                             
                             if (variant == "C124S-B1" ) {
                               current.dishes = c(
                                 1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,3
                                )
                               plates[1:length(current.dishes)] <- current.dishes
                             } else if (variant == "Human WT-B1" ) {
                               current.dishes  = c(
                                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3
                                 )
                               plates[1:length(current.dishes)] <- current.dishes
                               
                             } else if (variant == "D92N-B1" ) {
                               current.dishes <- c(
                                 1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3
                                 )
                               plates[1:length(current.dishes)] <- current.dishes
                               
                             } else if (variant == "G44D-B1" ) {
                               current.dishes <- c(
                                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3
                                 )
                               plates[1:length(current.dishes)] <- current.dishes
                               
                             } else if (variant == "GFP-B1" ) {
                               current.dishes <- c(
                                 1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3
                                 )
                               plates[1:length(current.dishes)] <- current.dishes
                               
                             } else if (variant == "C124S-B2" ) {
                               current.dishes <- c(
                                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3
                                 )
                               plates[1:length(current.dishes)] <- current.dishes
                               
                             } else if (variant == "Human WT-B2" ) {
                               current.dishes <- c(
                                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3
                                 )
                               plates[1:length(current.dishes)] <- current.dishes
                             } else if (variant == "Q171E-B2" ) {
                               current.dishes <- c(
                                 1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3
                                 )
                               plates[1:length(current.dishes)] <- current.dishes
                             } else if (variant == "G44D-B2" ) {
                               current.dishes <- c(
                                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3
                               )
                               plates[1:length(current.dishes)] <- current.dishes
                             } else if (variant == "GFP-B2" ) {
                               current.dishes <- c(
                                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3
                               )
                               plates[1:length(current.dishes)] <- current.dishes
                               
                             } else  if (length(plates) == 30){
                               plates[1:10] <- 1
                               plates[11:20] <- 2
                               plates[21:30] <- 3
                               
                             } else if (length(plates) == 50){
                               plates[1:17] <- 1
                               plates[18:35] <- 2
                               plates[36:50] <- 3
                             } else if (variant == "PTEN-12"){
                               plates[1:8] <- 1
                               plates[9:16] <- 2
                               plates[17:25] <- 3  
                             } else if (variant == "H93R-2"){
                               plates[1:8] <- 1
                               plates[9:15] <- 2
                               plates[16:22] <- 3
                             } else if (variant == "T131I" & xcol == 55){
                               plates[1:14] <- 1
                               plates[15:29] <- 2
                               plates[30:43] <- 3
                             } else {
                               # Weird cases.
                               T1 = ceiling(length(plates)/3)
                               
                               plates[T] <- 3
                               plates[1:T1] <- 1
                               plates[ (T1+1):(T1 + T1)] <- 2
                             }
                             
                             ret.df = data.frame(variant=variant,
                                          day=as.factor(variants$Day[1:(length(variants$Day))] ),
                                          plate=as.factor(plates), 
                                          axon=as.factor(variants$numeric[1:(length(variants$numeric))] )
                                          )
                             return(ret.df)
                           })
)

summary(axon.tidy)
unique(axon.tidy$variant)

axon.tidy$variant[grep(pattern = "-", x = axon.tidy$variant) ]

axon.tidy$variant <- gsub(pattern = "-.*", replacement ="", x = axon.tidy$variant)
axon.tidy$variant <- toupper(axon.tidy$variant)
axon.tidy$variant <- gsub(pattern = "HUMAN ", fixed = T, replacement ="", x = axon.tidy$variant)
axon.tidy$variant <- gsub(pattern = "PTEN", fixed = T, replacement ="WT", x = axon.tidy$variant)
axon.tidy$variant <- gsub(pattern = "G123D", fixed = T, replacement ="G132D", x = axon.tidy$variant)

table(axon.tidy$variant)


################ Handle 4A here!#####################
# # SPECIAL CASE

write.table(axon.tidy,
            file = "Data/PTEN-Axon.tsv", sep = '\t', row.names = F)



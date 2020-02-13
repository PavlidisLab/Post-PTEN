# To install: devtools::install_github("nacnudus/tidyxl")
source("utils.R")
inputFile = "RawData/FlyEclosion/fly_rawdata_2019-06-14.xlsx"

# Get a denormalized table of cells 
cells <- xlsx_cells(inputFile)

# Get formats with addresses
formats = xlsx_formats(inputFile)

# Parse out colours; make a factored version
COLOURS = unique(formats$local$fill$patternFill$fgColor$rgb[!is.na(formats$local$fill$patternFill$bgColor$rgb)])
COLOURS.factor <- as.factor(COLOURS)

# Inspect sheets
table(cells$sheet)

# Add vials by colour
cells$Vial <- -42 # Make a constant? default vial ID.
lapply(1:length(COLOURS), 
       FUN=function(i){
         
         cells$Vial[cells$local_format_id %in%
                  which(formats$local$fill$patternFill$fgColor$rgb == COLOURS[i])
                  ] <<- COLOURS.factor[i]
       })
table(cells$Vial)


# Keep relevant sheet
cells.eclosion <- cells[grepl( pattern="Eclosion", x=cells$sheet ),]
table(cells.eclosion$sheet)


# Identify variants
variant.columns.character <- unique(
  #cells.eclosion$character[grepl(pattern="PTEN|TDGY|attp2|TDDGY|Y65C|G36E|A309S|P357S|L320X|R173H|M1I|R14G|M134I", cells.eclosion$character) & !is.na(cells.eclosion$character) ]
  cells.eclosion$character[ grepl(pattern="PTEN|attp2|Y65C|G36E|A309S|P357S|L320X|R173H|M1I|R14G|M134I|M198I|K6I|C124R|N340H|N340D|K342N|V343L|L345V|I400V|K402N|Y68X|T240X", cells.eclosion$character) &
                            !grepl(pattern="TDGY|TDDGY|tdgy", cells.eclosion$character) &
                            !is.na(cells.eclosion$character) 
                           ]   
)

variant.characters <- cells.eclosion$character[cells.eclosion$character %in% variant.columns.character]
variant.columns <- cells.eclosion$col[cells.eclosion$character %in% variant.columns.character]
variant.rows <- cells.eclosion$row[cells.eclosion$character %in% variant.columns.character]
variant.sheets <- cells.eclosion$sheet[cells.eclosion$character %in% variant.columns.character]
variant.vials <- cells.eclosion$Vial[cells.eclosion$character %in% variant.columns.character]

fly.tidy <- do.call(rbind,
                    lapply(1:length(variant.columns),
                    FUN = function(IDX){
                      print(IDX)
                      xcol = variant.columns[IDX]
                      xrow = variant.rows[IDX]
                      xsheet = variant.sheets[IDX]
                      variant = variant.characters[IDX]
            
                      # Variants
                      variant.vialcounts = cells.eclosion[
                        cells.eclosion$sheet == xsheet &
                        cells.eclosion$col == xcol &
                        cells.eclosion$row > xrow &
                        !is.na(cells.eclosion$numeric), c("numeric", "Vial")
                      ]
                      variant.molten = melt(table(variant.vialcounts))
                      # Internal control
                      control.vialcounts = cells.eclosion[
                        cells.eclosion$sheet == xsheet & 
                        cells.eclosion$col == xcol + 1 &
                          cells.eclosion$row > xrow & 
                          !is.na(cells.eclosion$numeric), c("numeric", "Vial")
                      ]
                      control.molten = melt(table(control.vialcounts))
                      
                      ret.df = rbind(
                        data.frame(eclosion=xsheet,
                                 variant=variant,
                                 vial=as.factor(variant.molten$Vial),
                                 time=as.numeric(variant.molten$numeric),
                                 count=as.numeric(variant.molten$value),
                                 condition = "Experiment"),
                        
                        data.frame(eclosion=xsheet,
                                 variant=variant,
                                 vial=as.factor(control.molten$Vial),
                                 time=as.numeric(control.molten$numeric),
                                 count=as.numeric(control.molten$value),
                                 condition = "Control")
                        )
                      ret.df = ret.df[ret.df$count > 0,]
                      return(ret.df)
                    })
              )

################ Handle 4A here!#####################
# # SPECIAL CASE
# cells.eclosion$character[grep("4ala", cells.eclosion$character ) ]
# "Eclosion 5 & 6"

head(fly.tidy)
fly.tidy[fly.tidy$variant == "4 Ala",]
fly.tidy[fly.tidy$variant == "4A",]
table(fly.tidy$variant)

## Clean up eclosions
fly.tidy$eclosion = as.factor(as.numeric(fly.tidy$eclosion))

## Clean up 
fly.tidy$variant = as.character(fly.tidy$variant )
fly.tidy$variant = gsub(pattern = "PTEN ", replacement = "", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "C124S 4A", replacement = "C124S-4A", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "C124S(6)", fixed = T, replacement = "C124S", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "C124S(7)", fixed = T, replacement = "C124S", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "C136Mfs", replacement = "C136M", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "C136M", replacement = "I135fs", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "G36E/38E", replacement = "G36E", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "G36E/G38E", replacement = "G36E", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "I101 T", replacement = "I101T", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "WT(NEW)", fixed = T, replacement = "WT", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = " Y68N", fixed = T, replacement = "Y68N", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "G285X", fixed = T, replacement = "E285X", fly.tidy$variant)
fly.tidy$variant = gsub(pattern = "T240X", fixed = T, replacement = "Y240X", fly.tidy$variant)

table(fly.tidy$variant)
# Swap R130L and D252G .
## unique(fly.tidy$variant[grep(x=fly.tidy$variant, pattern = "R130")]) # Inspect
## unique(fly.tidy$variant[grep(x=fly.tidy$variant, pattern = "D252")]) # Inspect
fly.tidy$variant[grep(x=fly.tidy$variant, pattern = "R130L")] <-  "placeholder_R130L"
fly.tidy$variant[grep(x=fly.tidy$variant, pattern = "D252G")] <-  "R130L"
fly.tidy$variant[grep(x=fly.tidy$variant, pattern = "placeholder_R130L")] <-  "D252G"

# View(unique(fly.tidy$variant))
# head(fly.tidy)
# View(fly.tidy[order(fly.tidy$time),])
# summary(fly.tidy)

### DROP VARIANTS THAT ARE SEQUENCING ERRORS
SEQUENCING_ERRORS <-
c(
  #"G132D", # Confirmed OK
  # "F279L", # Variant construct wrong ## Confirmed OK afterall
  "A79T", # Variant construct wrong
  # "L295V", # Confirmed OK
  "E307Q" # Variant construct wrong; appears to be WT sequence
  #"A309S",  # Confirmed  OK
  # #"Q298E"  # Confirmed OK,
  # "K402N" # N too small; excluded.
  )
fly.tidy = fly.tidy[fly.tidy$variant %!in% SEQUENCING_ERRORS,]

write.table(fly.tidy,
          file = "Data/PTEN-Fly.tsv", sep = '\t', row.names = F)

# Write out 'identitfy` version of the fly data
fly.identity = data.frame(fly.tidy[
                            rep(seq_len(dim(fly.tidy)[1]), fly.tidy$count), ], 
                          row.names=NULL)
table(fly.identity$time)

# Sanity checks
fly.identity[fly.identity$time %% 12 != 0 ,]
setdiff(unique(fly.identity$variant),
        unique(fly.identity[fly.identity$time == 0,"variant"]))

write.table(fly.identity,
            file = "Data/PTEN-Fly.identity.tsv", sep = '\t', row.names = F)


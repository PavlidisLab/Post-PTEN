Source data files that are not immediately ready to be analysed; requires processing/clean-up.

## Corrections applied to raw data

Below is a description of changes applies to the raw data. To obtain corrected results, used data from the "Data/" folder instead of "RawData/", or produce "Data/" using the provided R scripts.

### Yeast

Below describes some correction applied to the raw data for the yeast assays, where:
    batch #1: RawData/YeastActivity/miniarray_v2
    batch #2: RawData/YeastActivity/180603_PTENMini_MB

Known typos:

    In the yeast batch #1: 
         mistakes   : "A125P", "G285X", "G38E",  "N356Q", "Q346R", "R357S", "V343M", "C136fs",  "C136M", "T240X", "T68X"
         correctionS: "A126P", "E285X", "G36E",  "N356D", "Q396R", "P357S", "V343L", "I135fs", "I135fs", "Y240X", "Y68X"
    
    In the yeast batch #2: A309S was mislabelled as E307Q and corrected before analysis.

Known exclusions:
    In the yeast batch #1: Y180H, E285X and Q396R were excluded because the construct was wrong and failed validation.
    In the yeast batch #2: M35V and E157G failed to mate and were unsucessful.

### Fly

In the fly eclosion assay, the following corrections were applied before analysis:

    "PTEN " prefixes removed from variant names.
    C124S 4A was relabelled to C124S-4A
    C124S(6) was relabelled to 	C124S
    C124S(7)  was relabelled to C124S
    C136Mfs  was relabelled to C136M
    C136M  was relabelled to I135fs
    G36E/38E  was relabelled to G36E
    G36E/G38E was relabelled to G36E
    I101 T was relabelled to I101T
    WT(NEW)  was relabelled to WT
    Y68N  was relabelled to Y68N
    G285X  was relabelled to E285X
    T240X  was relabelled to Y240X
    R130L  was relabelled to D252G
    D252G  was relabelled to R130L

The following varaints were dropped:
    A79T (Wrong construct)
    E307Q (Wrong construct, appears to be WT sequence)
    
### Axon outgrowth

In the axon outgrowth raw data:

    "HUMAN" was stripped from variant names
    PTEN was relabelled WT
    G123D was a typo and fixed to G132D

The batch information was not part of the source data is appended in the "axon.tidy.R"  script.

### Rat

In the "Rat/" data, the "_OE" suffix in variant names was removed.

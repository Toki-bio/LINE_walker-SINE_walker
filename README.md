## Extended Flank Extraction Optimization

The `--extended-flank` and `--try-extended-first` parameters enhance the alignment process by optimizing the flank extraction technique. 

- **--extended-flank**: This parameter triggers the use of bedtools slop to expand the coordinates of previous hits, allowing for the extraction of larger flanks.

- **--try-extended-first**: When enabled, the tool attempts to extract flanks using the extended hit coordinates first before resorting back to a standard genome search.

This optimization significantly improves the alignment speed and accuracy, particularly in challenging genomic regions where standard flanking might not suffice.

When using these parameters, users can expect a more efficient processing time without compromising on the quality of the genomic extraction.
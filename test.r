library(PRIME)
library(GenomicRanges)
library(rtracklayer)
dm_bed <- import("/Users/natsudanav/Documents/Dm_pooled_PRIMEloci_pred_0_75_Pooled_combinedcoreovlwith-d.bed", format="BED")
dm_ctss <- readRDS("/Users/natsudanav/Documents/Dm_all_ctss_rse.rds")

# Ensure all regions have the correct width
if (all(GenomicRanges::width(dm_bed) != 401)) {
  msg <- paste("⚠️ All regions in the object (GRanges) must have width",
               401,
               " : extend 401 bp from thick if existed")
  tc_gr <- extend_fromthick(tc_gr = dm_bed,
                            ext_dis = 200)
} else {
  msg <- paste("✅ All regions in the object (GRanges) have width", len_vec) # nolint: line_length_linter.
}
plc_message(msg)
validate_tc <- plc_validate_tc_object(tc_gr, dm_ctss, ext_dis = 200)

# Pre-processing
current_region_gr <- convert_strand_to_nostrand_gr(tc_gr)
current_region_gr <- remove_metadata_and_duplicates(current_region_gr)

library(dplyr)
ctss_gr <- cast_rse_to_granges(dm_ctss, assay = "counts", coln_assay = 1)
xx <- heatmapData(current_region_gr, ctss_gr, sparse=TRUE)
yy <- plc_batch_heatmapData(current_region_gr, ctss_gr, sparse=TRUE)


check_valid_profile_rownames(xx$`*`$`+`, chr_name)
check_valid_profile_rownames(xx$`*`$`-`, chr_name)


class(xx$`*`$`+`)
class(yy$`*`$`+`)

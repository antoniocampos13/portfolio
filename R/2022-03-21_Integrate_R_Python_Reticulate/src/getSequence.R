getSequence <- function(chr, start, end) {
  gr <- GenomicRanges::GRanges(glue::glue("{chr}:{start}-{end}"))
  refBase <- BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, gr)
  refBase <- as.data.frame(refBase)$x[1]
  
  return(refBase)
}

score_genome_gkmers_parallel<-function(seqnames, seqs, starts=1, target_kmers, kmer_scores, mc.cores=3, gapmin = NULL, gapmax = NULL)
{
  ## score very long sequences with gapped kmers and their enrichments, generating RleLists corresponding to .bw files
  ## target_kmers in the form of "xxxx0nxxxx"
  seqnames=as.character(seqnames)
  seqs=as.character(seqs)
  if(length(seqnames)==1) seqnames=paste0(seqnames,1:length(seqs))
  if(length(starts)==1) starts=rep(starts,length(seqs))
  if(any(is.null(gapmax),is.null(gapmin))){
    gaps= target_kmers %>% str_extract("\\d*(?=n)") %>% as.integer()
    if(is.null(gapmax)) gapmax=max(gaps)
    if(is.null(gapmin)) gapmin=min(gaps)
    print(glue::glue("determined gapmax is {gapmax}, gapmin is {gapmin}"))
  }

  scores=mclapply(1:length(seqs),function(ind){
    score_genome_gkmers_cpp(seqnames=seqnames[ind], seqs=seqs[ind], starts=starts[ind], target_kmers=target_kmers, kmer_scores=kmer_scores, gapmin=gapmin, gapmax=gapmax)
  },mc.cores = mc.cores) %>% set_names(seqnames)
  IRanges::RleList(scores)
}

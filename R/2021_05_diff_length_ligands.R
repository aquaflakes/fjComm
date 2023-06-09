
length_adjust<-function(seqs,output_length=40L,seq_in_middle=FALSE, use_rand_for_N=FALSE)
  # adjust all ligands to the same length, out put seq only
{
  seqs %<>% str_sub(1,output_length)
  Ns=output_length-nchar(seqs)
  if(use_rand_for_N){topaste=stringi::stri_rand_strings(length(seqs),Ns,"[ACGT]")}else{topaste=strrep("N",Ns)}
  if(seq_in_middle){
    leftNs=(Ns/2) %>% as.integer()
    rightNs=Ns-leftNs
    if(use_rand_for_N){leftNs=stringi::stri_rand_strings(length(seqs),leftNs,"[ACGT]")}else{leftNs=strrep("N",leftNs)}
    if(use_rand_for_N){rightNs=stringi::stri_rand_strings(length(seqs),rightNs,"[ACGT]")}else{rightNs=strrep("N",rightNs)}
    seqs=paste0(leftNs,seqs,rightNs)
  }else seqs=paste0(seqs,topaste)
  seqs
}

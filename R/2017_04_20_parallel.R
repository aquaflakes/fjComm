parallel <- function(list_or_vect= 1:10, fun= function(x){x}, workers=15, ...)
{
  pacman::p_load(BiocParallel)
  curr_workers= workers
    max_workers=parallel::detectCores()-2
    if (curr_workers>max_workers) curr_workers=max_workers

    # register(SerialParam())
    # register(MulticoreParam())
  parallelParam= MulticoreParam(workers = curr_workers)
  bplapply(list_or_vect, FUN = fun, BPPARAM = parallelParam,...)
}


parallel_cmd_gen<-function(cmd="samtools view -c {}",params=c("file1","file2","file3"),threads=8)
  glue::glue("parallel -j {threads} -k {cmd} ::: ") %>%paste0(paste0(params,collapse = " "))


mclapply<-function(params,FUN,mc.cores=2,loopMax=10,curr_loop=1,...){
  # robust parall lapply, will try multiple times
  print("mclapply loop: {curr_loop}" %>% glue::glue())
  result=parallel::mclapply(params,FUN,mc.cores = mc.cores,...)
  if(curr_loop>=loopMax) {return(result)} else {curr_loop=curr_loop+1}
  failed_Ind= map(result,~.x[1] %>% stringr::str_detect("Error")) %>% unlist() %>% which
  if(any(failed_Ind)){
    result_= mclapply(params[failed_Ind],FUN,mc.cores = mc.cores,loopMax=loopMax,curr_loop=curr_loop,...)
    for (Ind in seq_along(result_)) {result[[failed_Ind[Ind]]]=result_[[Ind]]}
  }
  return(result)
}

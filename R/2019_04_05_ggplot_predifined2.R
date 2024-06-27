gg_distMDS_2Dplot <- function(distMat, labels=NULL, color="black", fig_anno="", textsize=3, textalpha=0.5, pointsize=2,pointalpha=1)
{
  # visualize distance matrix with 2D xyplot
  # if(is.null(labels))labels=rownames(distMat)
  p_load(ggrepel)
  fit=cmdscale(distMat,eig=TRUE, k=2)
  plotdf=data.frame(x=fit$points[,1], y=fit$points[,2])#,color=color,labels=labels)
  p=ggplot(plotdf,aes(x,y))+geom_point(aes(color=I(color)),size=pointsize,alpha=pointalpha) +xlab("")+ylab("")+gg_anno_grob(x = 0.1,fig_anno)+ coord_fixed(ratio = 1)
  if(!is.null(labels)) p=p+geom_text_repel(aes(label=labels,color=I(color)),size=textsize,alpha=textalpha)
  p
}


gg_PCA_2Dplot<-function(Mat, labels=NULL, color="black", fig_anno="", textsize=3, textalpha=0.5, pointsize=2,pointalpha=1)
{
  # visualize PCA dim1,2 with 2D xyplot
  # if(is.null(labels))labels=rownames(distMat)
  p_load(ggrepel,FactoMineR)
  pca = PCA(Mat, graph = F,ncp = 5)
  pcaplot2=pca$ind$coord %>% as.data.frame() #%>%  mutate(text=labels,color=color)
  pca_eig <- pca$eig %>% round() %>% as.data.frame()
  p=ggplot(pcaplot2)+geom_point(aes(x=Dim.1, y=Dim.2,color=I(color)),size=pointsize
                                ,alpha=pointalpha) +xlab(paste0("PC1 (", pca_eig$`percentage of variance`[1], "%)"))+ylab(paste0("PC2 (", pca_eig$`percentage of variance`[2], "%)"))+gg_anno_grob(x = 0.1,fig_anno)+ coord_fixed(ratio = 1)
  if(!is.null(labels)) p=p+geom_text_repel(aes(x=Dim.1,y=Dim.2,label=labels,color=I(color)),size=textsize,alpha=textalpha)
  p
}

gg_hist <- function(Arr,binwidth=1,fill="white",color="grey30",...)
{
  force(Arr)
  ggplot()+geom_histogram(aes(Arr),binwidth = binwidth,fill=fill,color=color,...)
}

gg_histly <- function(Arr,binwidth=1,fill="white",color="grey30",...)
{
  pacman::p_load(plotly)
  (ggplot()+geom_histogram(aes(Arr),binwidth = binwidth,fill=fill,color=color,...)) %>% ggplotly()
}

gg_2Dheat_dotly <-function(df,scale_dot_size=TRUE)
  # visualize 2D MI heatmap top kmer combs
{
  dot_sizes= if(scale_dot_size) {df$topMIsum %>% {(.-min(.))/max(.)*1.2} } else{1}
  (ggplot(df)+geom_point(aes(pos1,pos2,color=topk1,text=paste0(topk1,pos2-pos1-3,topk2),size=I(dot_sizes) ))) %>% ggplotly()
}


gg_multi_ggoutput<-function(plotlist,ncol=1,col_to_row=FALSE,magnify=100)
{
  plotNum=length(plotlist)
  cols=ncol
  rows=(plotNum/cols) %>% base::ceiling()
  plotNum_=cols*rows

  options(warn=-1)
  layoutMat=matrix(1:plotNum_,ncol=ncol,byrow = TRUE);
  if(col_to_row) {layoutMat %<>% t; cols=rows; rows=ncol}
  plotlist %<>% purrr::map(~ggplotGrob(.))

  placeholder= createDummy(5,5) %>% ggplotGrob()
  plotlist= c(plotlist, vector("list", plotNum_-plotNum) %>% purrr::map(~placeholder))
  curr_p= createDummy((cols)*magnify,(rows)*magnify)+theme(plot.margin = margin())

  for(row in 1:nrow(layoutMat) )
  {
    for (col in 1:ncol(layoutMat))
    {
      curr_p= curr_p+annotation_custom(grob = plotlist[[layoutMat[row,col]]], xmin = (col-1)*magnify, xmax = (col)*magnify, ymin =(rows-row)*magnify, ymax = (rows-row+1)*magnify)
      # browser()
    }
  }
  options(warn=0)
  curr_p
}

complexheatmap_template=function(){cat("color_= circlize::colorRamp2(breaks = seq(-2,6,0.1), colors = colorspace::diverge_hcl(120,palette = 'Blue-Red3',rev = F)[40:120])
ComplexHeatmap::Heatmap(enrich_all,col = color_, cluster_rows = T, show_row_dend = F,cluster_columns = F,column_labels =ChIP_files$Symbol_Name[1:8], row_labels = all_motif_names,
    row_names_gp = gpar(fontsize=7),column_names_gp = gpar(fontsize=8),heatmap_legend_param =list(title = 'Log2 enrichment', at = c(-2, 0, 2,4,6),
    title_gp=gpar(fontsize=8),labels_gp=gpar(fontsize=7),grid_width = unit(3, 'mm') ))")}

dingkun_TF_logo <- function(pwmlist){
  pwmlist %<>% universalmotif::convert_motifs(class = "TFBSTools-PWMatrix")
  #'@ pwmlist  PWMatrixList
  motif_pics=BiocParallel::bplapply(1:length(pwmlist),
                                    function(x){
                                      t_pwmlist <- universalmotif::convert_motifs(pwmlist[[x]])
                                      motif_pic=suppressMessages(universalmotif::view_motifs(t_pwmlist,show.positions = F)+
                                                                   ylab("")+scale_y_continuous() +
                                                                   theme(axis.ticks.y =element_blank(),axis.line.y = element_blank(),axis.text.y = element_blank())  )
                                      ggplotGrob(motif_pic)
                                    },BPPARAM = BiocParallel::MulticoreParam(workers = 30))
  blank_plot <- ggplot() +
    theme_void()+
    xlim(0, 3) +
    ylim(0, 55)#)#+xlim(0,2)

  anno_plots <- c()
  anno_texts <- c()
  for (i in 1:length(pwmlist)) {
    cmd <- glue("annotation_custom(motif_pics[[{i}]],xmin =1,xmax={1+
                length(pwmlist[[i]]@profileMatrix[1,])*0.1},ymin={55-i},ymax={55-(i-1.5)})")
    cmd2 <- glue("annotate('text', x = 0.7, y = {55-(i-1.5)-0.5}, label = '{pwmlist[[i]]@name}',size = 4)")
    anno_plots <- append(anno_plots,cmd)
    anno_texts <- append(anno_texts,cmd2)
  }
  anno_plots %<>% paste0(collapse = "+")
  anno_texts %<>% paste0(collapse = "+")
  plot_text=paste0("blank_plot+ ",anno_plots,"+",anno_texts)
  plot_1 <- eval(parse(text = plot_text))
  plot_1
}

bin_point_density_to_alpha<-function(in_vect){
  # when drawing hex or bin geom, map the density of each bin to alpha, in order to simulate dot plot
  # lowest 4 values to different alpha, others to alpha=1
  ### example ###   ggplot() + geom_hex(aes(x,y,fill=I(bin_point_density_to_alpha(..density..))))
  vals=table(in_vect) %>% sort(decreasing = T) %>% names() %>% as.numeric(); vals=vals*1.1
  case_when(in_vect<vals[1] ~ 0.3,
            in_vect>vals[1]&in_vect<vals[2] ~ 0.5,
            in_vect>vals[2]&in_vect<vals[3] ~ 0.6,
            in_vect>vals[3]&in_vect<vals[4] ~ 0.8,
            .default = 1)
}


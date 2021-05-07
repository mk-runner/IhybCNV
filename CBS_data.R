library(DNAcopy)
"CBS_data"<-function(data,path){
  #data = read.table("/mnt/hgfs/2020.06.08--CNV_LOF/scalRD")
  #data = matrix(data[,1],ncol=96,byrow=F)

  #print(data)
  #data = data[-1, -(1:2)]
  #data = log2(data/2)
  head = matrix(0, nrow(data), 3)
  head[,1] = 1
  head[,2] = 1:nrow(data)
  head[,3] = 1:nrow(data)
  
  chrom <- rep(1, nrow(data))
  maploc <- 1:nrow(data)
  seg.file_g = matrix(0,1,6)
  seg.file_g_one = matrix(0,1,6)
  seg.file = matrix(0, nrow(data), 1)
  
  stac_amp = matrix(0, 1, nrow(data))
  stac_amp[1,] = 1:nrow(data)
  stac_amp_one = matrix(0, 1, nrow(data))
  
  stac_del = matrix(0, 1, nrow(data))
  stac_del[1,] = 1:nrow(data)
  stac_del_one = matrix(0, 1, nrow(data))
  
  
  for (j in 1:ncol(data)){
    
    #cat("sampl No,",j,"\n")
    #data <- CNA(data[,j],chrom,maploc)
    seg<- segment(CNA(data[,j],chrom,maploc))
    for (k in 1:length(seg$output$loc.start)){
      seg.file_g_one[1,1]=j
      seg.file_g_one[1,2]=1
      seg.file_g_one[1,3]=seg$output$loc.start[k]
      seg.file_g_one[1,4]=seg$output$loc.end[k]
      seg.file_g_one[1,5]=seg$output$num.mark[k]
      seg.file_g_one[1,6]=seg$output$seg.mean[k]
      seg.file_g=rbind(seg.file_g,seg.file_g_one)
      seg.file_g_one=matrix(0,1,6)
      
    }
    
    
  }
  seg.file_g = seg.file_g[-1,]
  write.table(seg.file_g,file=path,row.names=F,col.names=F,quote=F,sep="\t")
  
}




  



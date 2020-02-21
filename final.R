sacred <- function(adjin,demands,back=0){
  
  # STARTING COUNTER FOR PERFORMANCE EVALUATION
  timestart=Sys.time() 
  
  # LOADING ADJACENCY MATRIX
  adj=read.csv(adjin,sep = "",stringsAsFactors = F,header = F)
  adj=adj[-1,]
  adj=as.data.frame(adj)
  adj=adj[!duplicated(adj), ]
  adj$V3<-as.numeric(adj$V3)
  rownames(adj)<- NULL
  
  # LOADING DEMANDS
  dem=read.csv(demands,sep = "",stringsAsFactors = F,header = F)
  dem=dem[-1,]
  dem=as.data.frame(dem)
  rownames(dem)<- NULL
  
  # SPECTRUM MATRIX
  adjuniq=adj
  spectr <- as.data.frame(matrix(0,nr=nrow(adjuniq),nc=320))
  rownames(spectr)=paste(adjuniq$V1,adjuniq$V2,sep='-')
  
  # SAVING DIJKSTRA'S PATHS
  results=list()
  for(i in 1:(nrow(dem))){
    temp=dijkstra(adj = adj, start = dem$V1[i],end=dem$V2[i],b=dem$V3[i],backup = back)
    results[[i]]=temp
  }
  
  # SAVING ESTABLISHED DEMANDS 
  epit=as.data.frame(matrix(NA,nr=(nrow(dem)),nc=3))
  epit$V1='Rejected'
  # SPECTRUM ALLOCATION
  count=0
  count1=0
  slotut=0
  for(i in 1:length(results)){
    
    alarm=0
    alarm1=0
    
    count=count+1
    temp=results[[i]]$op_path
    temp=as.numeric(unlist(strsplit(temp,split=' -> ', fixed=TRUE)))
    arxes=c()
    grammes=c()
    for(j in 1:(length(temp)-1)){
      z=which(rownames(spectr) %in% c(paste(temp[j],temp[j+1],sep = '-') ))
      grammes=append(grammes,z)
      z2=which(spectr[z,]==0)
      sl=results[[i]]$slots
      
      for(k in 1:(length(z2)-(sl))){
        z22=z2[(k+sl)]-z2[k]
        if(z22==sl) arxes=append(arxes,z2[k])
      }
    }
    arxes=as.data.frame(table(arxes))
    arxes$arxes=as.numeric(as.character(arxes$arxes))
    arxes=arxes[arxes$Freq == (length(temp)-1),]
    if(dim(arxes)[1]!=0){
      alarm=1
      
    } else {alarm=2}
    
    # ADDING BACKUP PATH TO SPECTRUM
    if(results[[i]]$backup[1]=='None'){
      if(alarm==1){
        arxes=arxes[arxes$arxes == min(arxes$arxes),]
        z1=arxes$arxes
        spectr[grammes,(z1:(z1+sl-1))]=count
        
        slotut=slotut+((results[[i]]$slots)*(length(temp)-1))
        epit$V3[i]=slotut
        
        count=count+1
        
        epit$V2[i]=count1/(nrow(dem))
        epit$V1[i]='Established P'
        next
      } else{
        count=count+1
        
        count1=count1+1
        epit$V2[i]=count1/(nrow(dem))
        epit$V3[i]=slotut
        next
      }
    } else {
      tempz=results[[i]]$backup$path
      tempz=as.numeric(unlist(strsplit(tempz,split=' -> ', fixed=TRUE)))
      arxes1=c()
      grammes1=c()
      sl=results[[i]]$backup$slots
      
      for(j in 1:(length(tempz)-1)){
        z=which(rownames(spectr) %in% c(paste(tempz[j],tempz[j+1],sep = '-') ))
        grammes1=append(grammes1,z)
        z2=which(spectr[z,]==0)
        
        for(k in 1:(length(z2)-(sl))){
          z22=z2[(k+sl)]-z2[k]
          if(z22==sl) arxes1=append(arxes1,z2[k])
        }
      }
      arxes1=as.data.frame(table(arxes1))
      arxes1$arxes1=as.numeric(as.character(arxes1$arxes1))
      arxes1=arxes1[arxes1$Freq == (length(tempz)-1),]
      if(dim(arxes1)[1]!=0){
        alarm1=1
        
      } else {alarm1=2}
    }
    if(alarm==1 & alarm1==1){
      
      arxes=arxes[arxes$arxes == min(arxes$arxes),]
      z1=arxes$arxes
      spectr[grammes,(z1:(z1+sl-1))]=count
      
      count=count+1
      arxes1=arxes1[arxes1$arxes == min(arxes1$arxes),]
      z1=arxes1$arxes
      spectr[grammes1,(z1:(z1+sl-1))]=count
      
      slotut=slotut+((results[[i]]$slots)*(length(temp)-1))+((results[[i]]$backup$slots)*(length(tempz)-1))
      epit$V3[i]=slotut
      
      epit$V2[i]=count1/(nrow(dem))
      epit$V1[i]='Established B'
    } else if(alarm==1 & alarm1!=1){
      
      arxes=arxes[arxes$arxes == min(arxes$arxes),]
      z1=arxes$arxes
      spectr[grammes,(z1:(z1+sl-1))]=count
      
      count=count+1
      slotut=slotut+((results[[i]]$slots)*(length(temp)-1))
      epit$V3[i]=slotut
      epit$V2[i]=count1/(nrow(dem))
      epit$V1[i]='Established P'
    } 
    else {
      count=count+1
      
      count1=count1+1
      epit$V2[i]=count1/(nrow(dem))
      epit$V3[i]=slotut
    }
  }
  
  # STOPING COUNTER FOR PERFORMANCE EVALUATION
  timend=Sys.time()
  ttm=timend-timestart
  outup=list('results'=results,'spectrum'=spectr,'evaluation'=ttm,'blrate'=epit)
  return(outup)
}

### DIJKSTRA'S ALGORITHM FOR SHORTEST PATH CALCULATION ###

dijkstra <- function(adj,start,end,b,backup=0){
  
  ### WARNINGS ###
  if ( !(start %in% unique(append(unique(adj[,1]),unique(adj[,2]))) )){
    return('Wrong start node!')
  }
  
  if ( !(end %in% unique(append(unique(adj[,1]),unique(adj[,2]))) )){
    return('Wrong end node!')
  }
  
  if(class(adj)!="data.frame"){
    return('Adj should be a Data.frame.')
  }
  if(ncol(adj)!=3){
    return('Adj should have 3 columns.')
  }
  
  if(class(adj[,3]) != "numeric"){
    return('Adj 3rd column should contain numeric values.')
  }
  ### END OF WARNINGS ###
  
  # MODULATION MATRIX
  modul=as.data.frame(matrix(NA,nr=4,nc=4))
  modul$V1=c('16QAM','8QAM','QPSK','BPSK')
  modul$V2=c(0,801,1701,4601)
  modul$V3=c(800,1700,4600,9300)
  modul$V4=c(4,3,2,1)
  
  # ADJ CLONE WITH ALL THE DIRECTIONS
  adj1=adj
  for (row in 1:nrow(adj)){
    adj1[(nrow(adj1)+1),]=c(adj$V2[row],adj$V1[row],adj$V3[row])
  }
  adj1=adj1[!duplicated(adj1), ]
  adj=adj1
  
  # DISTANCES MATRIX
  dist=matrix(c(unique(append(unique(adj[,1]),unique(adj[,2])))),nr=length(unique(append(unique(adj[,1]),unique(adj[,2])))),nc=3)
  dist[,2]=Inf
  dist=as.data.frame(dist)
  dist$V2[dist$V1==start]=0
  dist$V3=start
  dist$V2=as.numeric(dist$V2)
  
  # EMPTY VECTOR TO HOLD USED NODES
  P=c() 
  
  # WHILE LOOP (RUNS UNTIL END NODE IS REACHED)
  while (end %in% dist$V1){
    z=dist[which.min(dist[,2]),1]
    n=adj$V1[adj$V2==z]
    
    if (length(P) != 0){
      n=n[-which(n %in% P)]
    }
    for (i in n){
      if (dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1] < dist$V2[dist$V1==i][1]){
        dist$V2[dist$V1==i] = dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1]
        
        l1=length(dist$V3[dist$V1==i])
        l2=length(dist$V3[dist$V1==z])
        l3=min(l1,l2)
        
        for (m in l3){
          dist$V3[dist$V1==i][m] =paste(dist$V3[dist$V1==z][m],i,sep = ' -> ')
        }
        
        if(length(dist$V3[dist$V1==z]) > l3){
          c=dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1]
          
          for (m in (l3 +1):length(dist$V3[dist$V1==z])){
            dist[(nrow(dist)+1),]=NA
            dist[(nrow(dist)),1]=i
            dist[(nrow(dist)),2]=c
            dist[(nrow(dist)),3]=paste(dist$V3[dist$V1==z][m],i,sep = ' -> ')
          }
        }
      } else if (dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1] == dist$V2[dist$V1==i][1]){
        c=dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1]
        for (m in 1:length(dist$V3[dist$V1==z])){
          dist[(nrow(dist)+1),]=NA
          dist[(nrow(dist)),1]=i
          dist[(nrow(dist)),2]=c
          dist[(nrow(dist)),3]=paste(dist$V3[dist$V1==z][m],i,sep = ' -> ')
        }
      }
    }
    L=dist[dist$V1==end,]
    dist=dist[-(which(dist$V1==z)),]
    P=append(P,z)
  }
  L=as.data.frame(L)
  
  # LINK DISJOINT BACKUP PATH
  adj1=adj
  for(m in 1:(length(L$V3))){
    temp=L$V3[m]
    temp=as.numeric(unlist(strsplit(temp,split=' -> ', fixed=TRUE)))
    for (w in 1:(length(temp)-1)) {
      if(length(which(adj1$V1==temp[w] & adj1$V2==temp[w+1])) != 0){
        adj1=adj1[-which(adj1$V1==temp[w] & adj1$V2==temp[w+1]),]
      }
    }
  }
  # NODE & LINK DISJOINT BACKUP PATH
  adj2=adj
  
  for(m in 1:(length(L$V3))){
    temp=L$V3[m]
    temp=as.numeric(unlist(strsplit(temp,split=' -> ', fixed=TRUE)))
    if(length(temp)==2){
      
      for(m in 1:(length(L$V3))){
        temp=L$V3[m]
        temp=as.numeric(unlist(strsplit(temp,split=' -> ', fixed=TRUE)))
        for (w in 1:(length(temp)-1)) {
          if(length(which(adj2$V1==temp[w] & adj2$V2==temp[w+1])) != 0){
            adj2=adj2[-which(adj2$V1==temp[w] & adj2$V2==temp[w+1]),]
          }
        }
      }
    } else{
      for (w in 2:(length(temp)-1)) {
        if(length(which(adj2$V1==temp[w] | adj2$V2==temp[w])) != 0){
          adj2=adj2[-which(adj2$V1==temp[w] | adj2$V2==temp[w]),]
        }
      }  
    }
  }
  
  # KIND OF PROTECTION
  if(backup==0){      # No protection
    for(l in 1:1){
      B=modul$V1[L$V2[l]>=modul$V2 & L$V2[l]<=modul$V3]  
      bi=ceiling(b/(10.7*(modul$V4[modul$V1==B] )))
    }
    retlist=list('op_path'=L$V3[l],'backup'='None','distance'=L$V2[l],'modulation'=B,'slots'=bi)
    return(retlist)
    
  } else if(backup==1){    # Link Disjoint
    for(l in 1:1){
      
      B=modul$V1[L$V2[l]>=modul$V2 & L$V2[l]<=modul$V3]  
      bi=ceiling(b/(10.7*(modul$V4[modul$V1==B] )))
      #print(paste('Optimal Primary path is:',L$V3[l],'with distance:',L$V2[l],'Modulation format:',B,'slots/link:',bi, sep = ' '))
    }
    retlist=list('op_path'=L$V3[l],'backup'=dijkstralt(adj = adj1,start,end,b),'distance'=L$V2[l],'modulation'=B,'slots'=bi)
    return(retlist)
    
  } else if(backup==2){   # Node & Link Disjoint
    for(l in 1:1){
      B=modul$V1[L$V2[l]>=modul$V2 & L$V2[l]<=modul$V3]  
      bi=ceiling(b/(10.7*(modul$V4[modul$V1==B] )))
    }
    retlist=list('op_path'=L$V3[l],'backup'=dijkstralt(adj = adj2,start,end,b),'distance'=L$V2[l],'modulation'=B,'slots'=bi)
    return(retlist)
  }
}

### SECONDARY DIJKSTRA'S FOR PROTECTION PATH (RUNS INSIDE dijkstra FUNCTION) ###

dijkstralt <- function(adj, start,end,b){
  
  ### WARNINGS ###
  if ( !(start %in% unique(append(unique(adj[,1]),unique(adj[,2]))) )){
    return('None')
  }
  
  if ( !(end %in% unique(append(unique(adj[,1]),unique(adj[,2]))) )){
    return('No backup path with node-link disjoint protection')
  }
  
  if(class(adj)!="data.frame"){
    return('Adj should be a Data.frame.')
  }
  if(ncol(adj)!=3){
    return('Adj should have 3 columns.')
  }
  
  if(class(adj[,3]) != "numeric"){
    return('Adj 3rd column should contain numeric values.')
  }
  ### END OF WARNINGS ###
  
  # MODULATION MATRIX
  modul=as.data.frame(matrix(NA,nr=4,nc=4))
  modul$V1=c('16QAM','8QAM','QPSK','BPSK')
  modul$V2=c(0,801,1701,4601)
  modul$V3=c(800,1700,4600,9300)
  modul$V4=c(4,3,2,1)
  
  # DISTANCES MATRIX
  dist=matrix(c(unique(append(unique(adj[,1]),unique(adj[,2])))),nr=length(unique(append(unique(adj[,1]),unique(adj[,2])))),nc=3)
  dist[,2]=Inf
  dist=as.data.frame(dist)
  dist$V2[dist$V1==start]=0
  dist$V3=start
  dist$V2=as.numeric(dist$V2)
  
  P=c()
  
  while (end %in% dist$V1){
    z=dist[which.min(dist[,2]),1]
    n=adj$V2[adj$V1==z]
    
    if (length(P) != 0){
      n=n[-which(n %in% P)]
    }
    for (i in n){
      if (dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1] < dist$V2[dist$V1==i][1]){
        dist$V2[dist$V1==i] = dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1]
        
        l1=length(dist$V3[dist$V1==i])
        l2=length(dist$V3[dist$V1==z])
        l3=min(l1,l2)
        
        for (m in l3){
          dist$V3[dist$V1==i][m] =paste(dist$V3[dist$V1==z][m],i,sep = ' -> ')
        }
        if(length(dist$V3[dist$V1==z]) > l3){
          c=dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)]
          
          for (m in (l3 +1):length(dist$V3[dist$V1==z])){
            dist[(nrow(dist)+1),]=NA
            dist[(nrow(dist)),1]=i
            dist[(nrow(dist)),2]=c[1]
            dist[(nrow(dist)),3]=paste(dist$V3[dist$V1==z][m],i,sep = ' -> ')
          }
        }
      } else if (dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1] == dist$V2[dist$V1==i][1]){
        c=dist$V2[dist$V1==z][1] + adj$V3[(adj$V1==i & adj$V2==z) | (adj$V2==i & adj$V1==z)][1]
        for (m in 1:length(dist$V3[dist$V1==z])){
          dist[(nrow(dist)+1),]=NA
          dist[(nrow(dist)),1]=i
          dist[(nrow(dist)),2]=c
          dist[(nrow(dist)),3]=paste(dist$V3[dist$V1==z][m],i,sep = ' -> ')
        }
      }
    }
    L=dist[dist$V1==end,]
    dist=dist[-(which(dist$V1==z)),]
    P=append(P,z)
  }
  L=as.data.frame(L)
  
  if(L$V2[1]=='Inf'){
    return('None')
  } else{
    for(l in 1:1){
      
      B=modul$V1[L$V2[l]>=modul$V2 & L$V2[l]<=modul$V3]  
      bi=ceiling(b/(10.7*(modul$V4[modul$V1==B] )))
    }
    retlist=list('path'=L$V3[l],'distance'=L$V2[l],'modulation'=B,'slots'=bi)
    return(retlist)
  }
}
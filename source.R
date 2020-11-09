cellphonedb <- function(pval, smeans, metadata){
  ###############
  ## Load Data ##
  ###############
  
  pval <- read.table(pval, sep="\t", stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
  smeans <- read.table(smeans, sep="\t", stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
  
  meta <- read.table(metadata, sep="\t", stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
  
  # Set rownames
  rownames(pval) <- pval$id_cp_interaction
  rownames(smeans) <- smeans$id_cp_interaction
  
  # Work with significant interactions
  df <- smeans[,1:12]
  
  #Create short version of smeans without NAs
  sm.short <- smeans[!apply(is.na(smeans[,13:ncol(smeans)]), 1, all),]
  
  #Create dfs of just pvalues
  smeans <- smeans[,13:ncol(smeans)]
  pval <- pval[,12:ncol(pval)]
  
  # Check if rownames are sample names
  if (all(rownames(head(meta)) != c(1,2,3,4,5,6))){
    meta <- cbind(rownames(meta), meta)
  } 
  
  colnames(meta)[1] <- "sample"
  
  ########################
  ## Count Interactions ##
  ########################
  
  grs <- unique(unlist(sapply(colnames(pval), function(z) strsplit(z, split="\\|")[[1]],simplify = FALSE)))
  
  ### Counts of interaction types between populations #####
  cnt_interactions <- function(xx, grs) {
    CNT <- colSums(xx< 0.05)
    MAT <- matrix(0, ncol=length(grs), nrow=length(grs), dimnames=list(grs,grs))
    for(case in names(CNT)) {
      cnt <- CNT[case]
      gr <- strsplit(case, split="\\|")[[1]]
      MAT[gr[1], gr[2]] <- MAT[gr[1], gr[2]] + cnt
      MAT[gr[2], gr[1]] <- MAT[gr[2], gr[1]] + cnt
    }
    return(MAT)
  }
  
  MAT <- cnt_interactions(xx=pval, grs)
  
  #--- Directional through ligands to receptor (removing receptor x receptor and ligand x ligand)
  cnts_sub1 <- colSums(pval[df$receptor_a=="True" & df$receptor_b=="False",] < 0.05)
  cnts_sub2 <- colSums(pval[df$receptor_a=="False" & df$receptor_b=="True",] < 0.05)
  
  # We have to revert their names to be sum with the previous
  names(cnts_sub2) <- sapply(names(cnts_sub2),function(z) paste(rev(strsplit(z,split="\\|")[[1]]),collapse = "|"))
  cnts_sub2 <- cnts_sub2[names(cnts_sub1)] # same order
  # Add both directionals
  cnts <- cnts_sub1 + cnts_sub2
  RECxLIG <- matrix(NA, ncol=length(grs), nrow=length(grs), dimnames=list(grs,grs))
  for(case in names(cnts)) {
    cnt <- cnts[case]
    gr <- strsplit(case, split="\\|")[[1]]
    RECxLIG[gr[1], gr[2]] <- cnt
  }
  rm(cnt, cnts, cnts_sub1, cnts_sub2, gr, case)
  
  ##########################################
  ## Count celltype senders and recievers ##
  ##########################################
  
  #Count senders and recievers
  cnts_senders <- colSums(pval[df$receptor_a=="True" & df$receptor_b=="False",] < 0.05) %>% as.data.frame()
  cnts_recievers <- colSums(pval[df$receptor_a=="False" & df$receptor_b=="True",] < 0.05) %>% as.data.frame()
  
  #Set column names
  colnames(cnts_senders) <- "count.send"
  colnames(cnts_recievers) <- "count.rec"
  
  #select senders and recievers from rownames
  cnts_senders$sender <- sapply(strsplit(row.names(cnts_senders),"\\|"), `[`, 1)
  cnts_recievers$recievers <- sapply(strsplit(row.names(cnts_recievers),"\\|"), `[`, 2)
  
  cnts_senders <- cnts_senders %>% group_by(sender) %>% mutate(count=sum(count.send)) %>% select(sender, count) %>% unique()
  cnts_recievers <- cnts_recievers %>% group_by(recievers) %>% mutate(count=sum(count.rec)) %>% select(recievers, count) %>% unique()
  
  #Order from largest to smallest
  cnts_senders <- cnts_senders[order(cnts_senders$count, decreasing = FALSE),]
  cnts_senders$sender <- factor(cnts_senders$sender, levels=cnts_senders$sender)
  
  cnts_recievers <- cnts_recievers[order(cnts_recievers$count, decreasing = FALSE),]
  cnts_recievers$recievers <- factor(cnts_recievers$recievers, levels=cnts_recievers$recievers)
  
  #Bar plots of senders and recievers
  plot.senders <- ggplot(cnts_senders, aes(x=sender, y=count)) +
    geom_bar(stat="identity") +
    coord_flip() +
    ylab("Number of Interactions") +
    xlab("Sender") +
    theme(text = element_text(size=15))
  
  plot.recievers <- ggplot(cnts_recievers, aes(x=recievers, y=count)) +
    geom_bar(stat="identity") +
    coord_flip() +
    ylab("Number of Interactions") +
    xlab("Reciever") +
    theme(text = element_text(size=15))
  
  ####################################
  ## Create network of interactions ##
  ####################################
  
  #create an igraph network object from the weighted adjacency matrix stored in pc
  net = igraph::graph_from_adjacency_matrix(MAT, weighted = TRUE)
  #remove multiple edges (meaning keep only one connection between each two cell clusters)
  
  net = igraph::simplify(net, edge.attr.comb = "max")
  
  #Some information to use in our plots
  Num.Interactions = E(net)$weight
  
  Strength = strength(net)
  
  #plot network with ggraph
  set.seed(113)
  lay = ggraph::create_layout(net, layout = "fr")
  network <- ggraph(lay) + 
    geom_edge_link(aes(color=Num.Interactions)) + 
    scale_edge_colour_gradient(  low = "#cacfcf",
                                 high = "#b52002",
                                 space = "Lab",
                                 na.value = "grey50",
                                 guide = "edge_colourbar") +
    geom_node_point(aes(size = Strength)) + 
    geom_node_text(aes(label = grs), repel=TRUE, size=6) +
    theme(panel.background = element_blank(),
          legend.text = element_text(size=14),
          legend.title = element_text(size=14))
  
  ########################
  ## Create Sankey plot ##
  ########################
  
  #Reformat data
  sk <- as.data.frame(RECxLIG) %>% mutate(Sender=row.names(.)) %>% gather(Reciever, Value, -one_of("Sender")) %>% filter(Value > 0)
  sk$Sender <- paste(sk$Sender, ".S", sep="")
  sk$Reciever <- paste(sk$Reciever, ".R", sep="")
  
  sk.nodes <- data.frame(node=unique(c(sk$Sender, sk$Reciever)))
  
  rp <- setNames(sk.nodes$node, as.numeric(row.names(sk.nodes))-1)
  
  sk$sender.n <-  as.numeric(names(rp)[match(sk$Sender, rp)])
  sk$reciever.n <-  as.numeric(names(rp)[match(sk$Reciever, rp)])
  
  sankey <- sankeyNetwork(Links = sk, Nodes = sk.nodes, Source = "sender.n",
                          Target = "reciever.n", Value = "Value", NodeID = "node",
                          units = "Interactions", fontSize = 12)
  
  out <- list(network=network, plot.senders=plot.senders, plot.recievers=plot.recievers, sankey=sankey)
  
  return(out)
  
}
  
  #############################################
  ## Create Sankey plot for gene of interest ##
  #############################################
sk.goi <- function(smeans, metadata, goi){  
  goi <- toupper(goi)
  
  smeans.goi <- read.table(smeans, sep="\t", stringsAsFactors = FALSE, header=TRUE, check.names=FALSE)
  
  smeans.goi <- smeans.goi[smeans.goi$gene_a == goi | smeans.goi$gene_b == goi,]
  
  smeans.goi <- smeans.goi %>% gather(Interaction, Value, 13:ncol(smeans.goi)) %>% filter(!is.na(Value))
  
  smeans.goi$Sender <- sapply(strsplit(smeans.goi$Interaction,"\\|"), `[`, 1)
  smeans.goi$Reciever <- sapply(strsplit(smeans.goi$Interaction,"\\|"), `[`, 2)
  
  #Create dfs for each link in sankey
  ln1 <- data.frame(Sender=paste(smeans.goi$Sender, ".s", sep=""), Reciever=smeans.goi$gene_a, value=1)
  ln2 <- data.frame(Sender=smeans.goi$gene_a, Reciever=smeans.goi$gene_b, value=smeans.goi$Value)
  ln3 <- data.frame(Sender=smeans.goi$gene_b, Reciever=paste(smeans.goi$Reciever, ".r", sep=""), value=1)
  
  lns <- rbind(ln1, ln2, ln3)
  
  #Reformat data
  lns.nodes <- data.frame(node=unique(c(lns$Sender, lns$Reciever)))
  
  rp <- setNames(lns.nodes$node, as.numeric(row.names(lns.nodes))-1)
  
  lns$sender.n <-  as.numeric(names(rp)[match(lns$Sender, rp)])
  lns$reciever.n <-  as.numeric(names(rp)[match(lns$Reciever, rp)])
  
  sankey.goi <- sankeyNetwork(Links = lns, Nodes = lns.nodes, Source = "sender.n",
                          Target = "reciever.n", Value = "value", NodeID = "node",
                          units = "Interactions", fontSize = 12)
  
  return(sankey.goi)
  
}
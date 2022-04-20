#' Set of functions to either simulate or form the binomial distribution of Strand-seq pools for a total of a 1000 individuals.
#' 
#'
#' The defined functions can either simulate the steps of the Pooling process, or generate the 
#' Binomial Distribution of the density of the number of cells successfully generated from each sample. 
#' wells.
#'
#'
#'
#' @name pooling_experiment
#' @author Celia Tsapalou, Wolfram Hoeps
#'
#' Function \code{generate_random_cells} simulates the amount of cells successfully extracted from each pool, assuming 10000 cells in each sample.
#' @param all_samples the total number of samples used in the pooling experiment.
#' @param samples_per_pool number of samples to be used in each one of the pools.
#' @param well define capacity of the cell culture dish.
#' @export
#'


generate_random_cells <- function(all_samples, samples_per_pool, well) {

  pools = all_samples / samples_per_pool
  dropout<-c(1/2)
  cells<-10000
  
  
  # make an empty list to name all the cell counts from all the samples
  
  cell_list <- matrix(0, samples_per_pool*cells, 3)
  colnames(cell_list)<-c('sample', 'cell', 'name')
  
  
  #Generate the list of names for the 10.000 cells in each sample
  
  for (i in 1:samples_per_pool) {
    for (j in 1:cells){ #10000 cells per pool
      cell_list[cells*(i-1)+j, 1]<-i #column 1
      cell_list[cells*(i-1)+j, 2]<-j #column 2
      cell_list[cells*(i-1)+j, 3]<-paste0('sample',i,'_cell',j) #column 3
    }
    
  }
  
  
  random <- data.frame(matrix(vector(),0,1, dimnames=list(c(), c('cell_count'))), stringsAsFactors = F)
  
  #Run it for all the pools
  
  for (p in 1:pools) {
      
    #Make a list with the cells chosen after Random Sampling from the total samples (samples_per_pool*10.000)
    #sample(x, size, replace = FALSE, prob = NULL) 
      
    cell_list_sort<-cell_list[sample(c(1:(samples_per_pool*cells)), (well)),]
    
    
    # Assume Dropout and keep only the alive cells from the well 
    cell_list_sort_live <- cell_list_sort[sample(c(1:well), well*dropout),]
    
    
    # Create empty matrix for plot only for the pooled samples
    cell_number_bar <- matrix(0, samples_per_pool, 1) #matrix with 0 values, 40 rows, 1 column
    rownames(cell_number_bar) <- paste0('pool', p, '_sample_', c(1:samples_per_pool))
    colnames(cell_number_bar) <- 'cell_count'
     
     
    #Alive cells in the samples_per_pool samples
    for (k in 1:samples_per_pool){
      cell_number_bar[k,1] <-sum(cell_list_sort_live[,1]==k)
    }
  
    random <-rbind(random, cell_number_bar) 
    
    }
    
    x<- as.data.frame(rbind(table(random)))
    x_seq = 0:(length(as.numeric(x[1,]))-1)
    counts <- data.frame(x=x_seq, y=as.numeric(x[1,]))
    return(counts)
  

}



#Size=number of cells
#Number of successes : from 0 to as many cells in the well*dropout
#Probability= 1/total number of samples in the pool

#' Function \code{generate_distribution} generates the binomial distribution of the probability of a certain number of successfully extracted cells (x) in a certain number of trials (size) where the probability of success on each trial is fixed (prob).
#' @param all_samples the total number of samples used in the pooling experiment.
#' @param samples_per_pool number of samples to be used in each one of the pools.
#' @param well define capacity of the cell culture dish.
#' @export
#'

generate_distribution <- function(all_samples, samples_per_pool, well) {
  
  pools = all_samples / samples_per_pool
  dropout<-c(1/2)
  
  new_random <- dbinom(x=0:(well*dropout), size= (well * dropout) , prob=1/samples_per_pool)
  new_random <-as.data.frame(new_random)
  new_random <- (new_random$new_random * all_samples)
  new_random <-as.data.frame(new_random)
  
  
  x_seq = 0:(well * dropout)
  counts <- data.frame(x=x_seq, y=new_random$new_random)
  return(counts)
  

  
}

#' Function \code{plot_cells} plots and saves the output of the simulation case and the binomial distribution.
#' @param all_samples the total number of samples used in the pooling experiment.
#' @param samples_per_pool number of samples to be used in each one of the pools.
#' @param well define capacity of the cell culture dish.
#' @param outfile define name of the output plot.
#' @export
#'

plot_cells<-function(samples_per_pool, well, counts, outfile) {

  p<-ggplot(counts, aes(x=x, y=y)) +
    geom_col(color="black", fill="pink")+
    labs(title=paste0(samples_per_pool," samples per Pool. ", well,"-well Plate"),x="Cell Number", y = "Counts") +
    scale_x_continuous(breaks=seq(-1,20), limits=c(-1,20))+
    scale_y_continuous(limits=c(0,500))+
    theme_bw()
  
  
  ggsave(filename = paste0(samples_per_pool," samples per Pool. ", well,"-well Plate", outfile),
                  plot = p, 
                  width = 20, 
                  height = 20, 
                  units = 'cm',
                  dpi = 300, device= "pdf")
  
  print(paste0(samples_per_pool," samples per Pool. ", well,"-well Plate", outfile))
  
  return(p)

}


















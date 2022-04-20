#' Set of functions to either simulate or form the binomial distribution of Strand-seq pools for a total of a 1000 individuals.
#' 
#'
#' The defined functions can either simulate the steps of the Pooling process, or generate the 
#' Binomial Distribution of the density of the number of cells successfully generated from each sample. 
#' wells.
#'
#'
#'
#' @name Pooling_Experiment
#' @author Celia Tsapalou, Wolfram Hoeps
#'
#' Function \code{generate_random_cells} simulates the amount of cells successfully extracted from each pool, assuming 10000 cells in each sample.
#' @param all_samples The total number of samples used in the pooling experiment.
#' @param samples_per_pool Number of samples to be divided in each one of the pools.
#' @param well The capacity of the cell culture dish.
#' @export
#'


generate_random_cells <- function(all_samples, samples_per_pool, well) {

  pools = all_samples / samples_per_pool
  
  # Dropout Rate: 50% of the cells are not successfully pooled
  dropout<-c(1/2)
 
  # Total number of cells for each individual of the experiment
  cells<-10000
  
  
  # Make an empty list to name all the Cell Counts for the Total Number of Samples
  
  cell_list <- matrix(0, samples_per_pool*cells, 3)
  colnames(cell_list)<-c('sample', 'cell', 'name')
  
  
  # Generate the list of names for the 10.000 cells  each sample
  
  for (i in 1:samples_per_pool) {
    for (j in 1:cells){ #10000 cells per pool
      cell_list[cells*(i-1)+j, 1]<-i #column 1
      cell_list[cells*(i-1)+j, 2]<-j #column 2
      cell_list[cells*(i-1)+j, 3]<-paste0('sample',i,'_cell',j) #column 3
    }
    
  }
  
  
  random <- data.frame(matrix(vector(),0,1, dimnames=list(c(), c('cell_count'))), stringsAsFactors = F)
  
  # Run it for the total number of pools
  
  for (p in 1:pools) {
      
    # Make a list with the cells extracted from each individual after Random Sampling (samples_per_pool * 10.000)
    # sample(x, size, replace = FALSE, prob = NULL) 
      
    cell_list_sort<-cell_list[sample(c(1:(samples_per_pool*cells)), (well)),]
    
    
    # Assume 50% dropout rate and keep only the alive cells from each well for the next steps 
    
    cell_list_sort_live <- cell_list_sort[sample(c(1:well), well*dropout),]
    
    
    # Create empty matrix 
    
    cell_number_bar <- matrix(0, samples_per_pool, 1) 
    rownames(cell_number_bar) <- paste0('pool', p, '_sample_', c(1:samples_per_pool))
    colnames(cell_number_bar) <- 'cell_count'
     
     
    # Alive cells in the samples_per_pool samples
    
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



# Size = Total Number of cells
# Number of Successes = from 0 up to the total number of alive cells in each well after dropout
# Probability = 1/Total Number of Samples in the Pool

#' Function \code{generate_distribution} generates the binomial distribution of the probability of a certain number of successfully extracted cells (x) in a certain number of trials (size) where the probability of success on each trial is fixed (prob).
#' @param all_samples The total number of samples used in the pooling experiment.
#' @param samples_per_pool Number of samples to be divided in each one of the pools.
#' @param well The capacity of the cell culture dish.
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
#' @param all_samples The total number of samples used in the pooling experiment.
#' @param samples_per_pool Number of samples to be divided in each one of the pools.
#' @param well The capacity of the cell culture dish.
#' @param outfile Define the name of the output plot.
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


















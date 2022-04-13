library(ggplot2)
library(scales)

#FUNCTION NO1: SIMULATION
generate_random_cells <- function(all_samples, samples_per_pool, well) {


pools = all_samples / samples_per_pool
dropout<-c(1/2)


# make an empty list to name all the cell counts from all the samples
cell_list <- matrix(0, samples_per_pool*10000, 3)
colnames(cell_list)<-c('sample', 'cell', 'name')


#Generate the list of names for the 10.000 cells in each sample
for (i in 1:samples_per_pool){
  for (j in 1:10000){ #10000 cells per pool
    cell_list[10000*(i-1)+j, 1]<-i #column 1
    cell_list[10000*(i-1)+j, 2]<-j #column 2
    cell_list[10000*(i-1)+j, 3]<-paste0('sample',i,'_cell',j) #column 3
  }
}


random <- data.frame(matrix(vector(),0,1, dimnames=list(c(), c('cell_count'))), stringsAsFactors = F)

#Run it for all the pools
for (p in 1:pools) {
    
#Make a list with the cells after Random Sampling (samples_per_pool*10.000)
#sample(x, size, replace = FALSE, prob = NULL) 
  
cell_list_sort<-cell_list[sample(c(1:(samples_per_pool*10000)), (well)),]


#Dropout or keep only the alive cells from the well (96 or 386)
cell_list_sort_live <- cell_list_sort[sample(c(1:(well)), well*dropout),]


#Create empty matrix for plot only for the pooled samples
cell_number_bar <- matrix(0, samples_per_pool, 1) #matrix with 0 values, 40 rows, 1 column
rownames(cell_number_bar) <- paste0('pool', p, '_sample_', c(1:samples_per_pool))
colnames(cell_number_bar) <- 'cell_count'
 
 
#Alive cells in the samples_per_pool samples
  for (k in 1:samples_per_pool){
    cell_number_bar[k,1] <-sum(cell_list_sort_live[,1]==k)
  }


random <-rbind(random, cell_number_bar) 
x<- as.data.frame(rbind(table(random)))
x_seq = 0:(length(as.numeric(x[1,]))-1)
counts <<- as.data.frame(x=x_seq, y=as.numeric(x[1,]))
return(counts)

}
}

result_random <- generate_random_cells(1000,40,96)


##########################################################

#FUNCTION NO2: PLOT DISTRIBUTION

##### Function for generating Distribution ####
distribution_plot<-function(samples_per_pool, well, counts) {

p<-ggplot(counts, aes(x=x, y=y)) +
  geom_col(color="black", fill="pink")+
  labs(title=paste0(samples_per_pool," samples per Pool. ", well*2,"-well Plate"),x="Cell Number", y = "Counts") +
  scale_x_continuous(breaks=seq(-1,20), limits=c(-1,20))+
  scale_y_continuous(limits=c(0,500))+
  theme_bw()
p

}

distribution_plot(40,96, result_random)



###########################################################################################

############################## FUNCTION NO3 ###############################################

#Binomial Distribution

generate_distribution <- function(all_samples, samples_per_pool, well) {
pools = all_samples / samples_per_pool
dropout<-c(1/2)
for (p in 1:pools) {

new_random <<- dbinom(x=0:(well*dropout), size= (well * dropout) , prob=1/samples_per_pool)
new_random<<-as.data.frame(new_random)
new_random <<- (new_random$new_random * all_samples)
new_random<<-as.data.frame(new_random)

return(new_random)


x<<-as.data.frame(rbind(table(new_random)))
x_seq=0:(length(as.numeric(x[1,]))-1)
counts<<-data.frame(x=x_seq, y=as.numeric(x[1,]))

}

}


#Size=number of cells
#Number of successes : from 0 to as many cells in the well*dropout
#Probability= 1/total number of samples in the pool
result_binomial_distribution <-generate_distribution(1000, 40, 96)

#plot(c(0:48), dbinom(x=0:(well*dropout), size= (well * dropout) , prob=1/samples_per_pool))


##### Function for generating Distribution ######
distribution_plot<-function(samples_per_pool, well, new_random) {
  

p<-ggplot(counts, aes(x=x, y=y)) +
    geom_col(color="black", fill="pink")+
    labs(title=paste0(samples_per_pool," samples per Pool. ", well*2,"-well Plate"), x="Cell Number", y = "Counts") +
    scale_x_continuous(breaks=seq(-1,20), limits=c(-1,20))+
    scale_y_continuous(limits=c(0,1))+
    theme_bw()
p
  
}

distribution_plot(40,96, result_binomial_distribution)

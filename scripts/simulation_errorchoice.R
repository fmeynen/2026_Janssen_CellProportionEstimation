#######################
# Idea: Simulate proportions of 10 cell types and identify their errors. See how often they exceed a certain treshold

# Step 1: choose the proportions
  #* Consider a monotone decreasing curve (for example a Beta(alpha, 1))
  #* Take 10 equidistant points on this curve, their values are unscaled proportions
  #* Rescale the 10 proportions such that their sum is equal to 1.
  #* The result is a vector of proportions

# Step 2: Simulate proportions
  #* For the vector of proportions, take a sample from a multinomial distribution, giving us 10 counts that sum to n
  #* Calculate the observed proportions by dividing the counts by n
  #* Calculate the error metrics:
    #* absolute error AE: abs(observed prop - actual prop)
    #* absolute relative error ARE: abs(observed prop - actual prop) / p
  #* This results in a vector of errors.
  #* Note the largest of these errors and note whether it is above or below a certain treshold tau
  #* This results in a binary value for each error metric
  #* Output is the binary values as well as the proportion that resulted in the largest error.

# Step 3: Repeat step 2
 #* Repeat step 2 B times for different values of tau

# Step 4: Graph results
  #* a line graph with on the x-axis the treshold tau and on y the success rate (how much simulations are below the treshold)
  #* same as above, but with both the AE and the ARE
  #* a histogram that shows which of the proportions resulted in the largest error
  #* same as above, but with both the AE and the ARE

#Possible future changes:
  #* Allow for overdispersion in the multinomial distribution
  #* Allow for positive correlations in the multinomial distribution
  #* Consider other monotone decreasing or increasing curves to take proportions from
  #* Consider other error metrics
  #* Extract the largest error to graph the distribution of largest errors.
###########
# These scripts characterize endosperm teM genes expression pattern and structural features

## 1. Explore the structure features 
> 1. The number of exons (Fig 3C)
>
> 2. The total length of CDS, intons, UTRs. (Fig 3D)

## 2. Explore the expression features 
> 1. The maximum expression level in endosperm (Fig 3A)
> 
> 2. Expression fold change between endosperm and other tissues (Fig 3B)
> 
> 3. Transcipts propotion of endosperm teM genes among all the expressed core genes (Fig 3E)

## 3. Perform K-means clustering on endosperm teM genes 
> 1. Calculate the Z-scores using the algorithm:
     $$Z = \frac{X - \mu}{\sigma}$$
    
   Where:
- \( X \) = observed value  
- \( \mu \) = mean of the population  
- \( \sigma \) = standard deviation  
> 3. Using the Elbow Method and Silhouette Method to find the optimal K (Fig S2A-S2B)
> 4. Perform K-means clustering to divide the endosperm teM genes into two groups based on the expression pattern by time (Fig 3G)



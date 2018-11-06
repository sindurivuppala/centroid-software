# centroid-software adapted from Arun et al.; i-xxx, International Journal of Biochemistry Research & Review, Supplementary material
centroid <- vector("numeric")
parent_peptide_mass <- 2136.3
install.packages("readxl")
library(readxl)
data <- read_xlsx("C:/Users/Courses/R statistics and analytics/Braf_digest_complex.xlsx")
data

z <- 1
colnames(data) <- c("mass","intensity")

ram <- which(data$mass >= parent_peptide_mass)
rama <- data[ram[1],]
right <- rama[1,]+(1/(2*z))
left <- rama[1,]-(1/(2*z))
right_vect <- which(data$mass >= right[1,1])
rit <- right_vect[1]
left_vect <- which(data$mass <= left[1,1])
lft <- left_vect[length(left_vect)]
interval <- data[lft:rit,]
interval


kri <- max(interval$intensity)
bal <- which(interval$intensity == kri)
kal <- interval[bal,]
kal
ppm_exp <- kal[1,1]
ppm_exp

# Assign the experimental parent peptide mass value to another variable 'peak_mid_grab' 
peak_mid_grab <- as.numeric(ppm_exp) 
peak_mid_grab

# Perform a loop to describe the peaks 
  a <- data[1]
  b <- data[1]
peaks <- vector("numeric") 

for(a in data[1]){ 
  
  ram <- which(data$mass >= peak_mid_grab) 
  ram
  rama <- data[ram[1],] 
  right <-rama[1] + (1/(2*z)) 
  left <- rama[1] - (1/(2*z)) 
  right_vect <- which(data$mass >= right[1,1]) 
  rit <- right_vect[1] 
  left_vect <- which(data$mass <= left[1,1]) 
  lft <- left_vect[length(left_vect)] 
  interval <- data[lft:rit,] 

  kri <- max(interval$intensity) 
  bal <- which(interval$intensity == kri) 
  kal <- interval[bal,] 
  ppm_mid <- kal[1,1] 

## Trace the left half of distribution 
  position_mbyz <- which(data$mass >= as.numeric(ppm_mid))
  click <- position_mbyz
  i <- c(data$mass)
  t <- c(data$mass)
  for(i in t){ 
    del1 <- as.numeric(data$intensity[click]) 
    del2 <- as.numeric(data$intensity[click - 1]) 
    if(del2>del1){ 
      break 
      } 
    click <- click - 1 
    } 

## Trace the right half of distribution 
  right_click <- position_mbyz 
  i <- c(data$mass)
  t <- c(data$mass)
  for(i in t){ 
    acer1 <- as.numeric(data$intensity[right_click]) 
    acer2 <- as.numeric(data$intensity[right_click + 1] )
    if(acer2 > acer1){ 
      break } 
    right_click <- right_click + 1 
    } 

  ## Subset the entire datapoints corresponding to the single peak 
  peak_dist <- data[click:right_click,] 
  
  #plot(peak_dist) 
  #windows() 
  
  ## Calculate maximium of peak intensity and check if it is increasing again 
  peak_max <- max(peak_dist$intensity) 
  peaks <- c(peaks, peak_max) 
  
  ## Increment the mid mbyz value 
  peak_mid_grab <- peak_mid_grab + 1 
}

  #View(peaks, title = "entire peak list") 


# Type1 : The envelope intensity is like - 50,40,30,20,10. 
# Type2 : The envelope intensity is like - 20,30,50,40,30,20,10. 
# Type3 : The envelope intensity is like - 20,20,30,50,40,30,20,10. 
#                                        or 50,50,40,30,20,10. 
# Create a major loop to decide on the path and to retrieve the end peak position of the envelope as 'empire' 

u <- 1:length(peaks) 
p <- 1:(length(peaks)-1) 

if(peaks[1] > peaks[2]){ 
  for(p in u){
    l <- 1 
    o <- l+1  
    if(peaks[l] > peaks[o]){
      empire <- l+1 
      l <- l+1 
      o <- o+1 
    } 
    else{ 
      empire <- l 
    }   
  } 
}


# Type2 
u <- 1:length(peaks) 
p <- 1:(length(peaks)-1) 

if(peaks[1] < peaks[2]){
  for(p in u){ 
    l <- 1 
    o <- l+1
    if(peaks[l] < peaks[o]){ 
      empire <- l+1 
      l <- l+1 
      o <- o+1 

    }
  }

  if(peaks[1] > peaks[2]){
  for(p in u){ 
    if(peaks[l] > peaks[o]){
      empire <- l+1 
      l <- l+1 
      o <- o+1 
        }
      }
    }
  }


# Type3 
  if(peaks[1] == peaks[2]){
    if(peaks[2] < peaks[3]){
      l <- 2
      o <- l+1
      for(p in u){
        if(peaks[l] < peaks[o]){
          empire <- l+1
          l <- l+1
          o <- o+1
        }
      }
    }
    else{ 
      if(peaks[2] > peaks[3])
      { l <- 2 
        o <- l+1
        for(p in u){
          if(peaks[l] > peaks[o]){
            empire <- l+1 
            l <- l+1
            o <- o+1 
          }
        }
      }
    }
  }


print(empire) 

# Construct the envelope data.frame 
## Find the starting point using the exp max of parent peptide mass 
position_mbyz <- which(data$mass == as.numeric(ppm_exp)) 
click <- as.numeric(position_mbyz )
i <- data[1]
t <- data[1]
for(i in t){
  del1 <- as.numeric(data$intensity[click])
  del2 <- as.numeric(data$intensity[click - 1])
  if(del2>del1){
    break 
  } 
  click <- click - 1 
}


## Find the ending point using the max peak intensity corresponding to 'empire' value 
## Find the times of increment in mbyz value to reach the point
end_incre <- empire - 1
midddle_mbyz_of_end <- as.numeric(ppm_exp + (end_incre * (1/z)))
rama <- midddle_mbyz_of_end 
right <-rama + (1/(2*z)) 
left <- rama - (1/(2*z)) 
right_vect <- which(data$mass >= right) 
rit <- right_vect[1] 
left_vect <- which(data$mass <= left) 
lft <- left_vect[length(left_vect)] 
interval <- data[lft:rit,] 

kri <- max(interval$intensity) 
bal <- which(interval$intensity == kri) 
kal <- interval[bal,] 
midddle_mbyz_of_end_exp <- kal[1,1]
position_mbyz <- which (data$mass == as.numeric(midddle_mbyz_of_end_exp)) 
right_click <- position_mbyz 
i <- data[1] 
t <- data[1]
for (i in t){ 
  acer1 <- as.numeric(data$intensity[right_click]) 
  acer2 <- as.numeric(data$intensity[right_click + 1]) 
  if(acer2 > acer1){ 
    break 
  } 
  right_click <- right_click + 1 
}



## name the new data.frame as 'frame1' 
frame1 <- data [click:right_click,] 

# Find the base peak intensity from 'frame1' 
bp <- max(frame1$intensity) 
bp_mbyz <- frame1$mass[which(frame1$intensity == bp)[1]] 
cut_off <- 0.20 * bp 

## Refine the envelope and store it in 'frame2' 
### Perform a loop to describe the peaks 
peak_mid_grab <- ppm_exp 
vibe <- data.frame() 
peak_max1 <- min(data$intensity) 
a <- data[1]
b <- data[1]
peaks <- vector("numeric") 
peaksr <- vector("numeric") 
for(a in b){
  
  ### Define the distribution of the peak according to condition 1 : successive point intensity 
  ### should be less than the preceeding intensity 
  ### Grab the exact middle mbyz based on maximum intensity 
  
ram <- which(data$mass >= as.numeric(peak_mid_grab)) 
rama <- data[ram[1],] 
right <-rama[1] + (1/(2*z)) 
left <- rama[1] - (1/(2*z)) 
right_vect <- which(data$mass >= right[1,1]) 
rit <- right_vect[1] 
left_vect <- which (data$mass <= left[1,1]) 
lft <- left_vect[length(left_vect)] 
interval <- data[lft:rit,] 

kri <- max(interval$intensity) 
bal <- which(interval$intensity == kri) 
kal <- interval[bal,] 
ppm_mid <- kal[1,1] 

#### Trace the left half of distribution 
position_mbyz <- which (data$mass >= as.numeric(ppm_mid) 
click <- position_mbyz 
i <- 1:1000 
t <- 1:1000 
for(i in t){
  del1 <- data$intensity[click] 
  del2 <- data$intensity[click - 1] 
  if(del2>del1){ 
    break 
  } 
  click <- click - 1 
}


#### Trace the right half of distribution 
right_click <- position_mbyz 
i <- 1:1000 
t <- 1:1000 
for(i in t){
  acer1 <- data$intensity[right_click] 
  acer2 <- data$intensity[right_click + 1] 
  if(acer2 > acer1){
    break 
  } 
  right_click <- right_click + 1 
}


### Subset the entire datapoints corresponding to the single peak 
peak_dist <- data[click:right_click,] 
#plot(peak_dist) 
#windows() 
### Calculate maximium of peak intensity and check if it is increasing again 
peak_max <- max(peak_dist$intensity) 
vibe_mbyz <- peak_dist$mass[peak_dist$intensity == peak_max] 
flash_peak <- data [data$mass == vibe_mbyz[1],] 
if(peak_dist$mass[1] > bp_mbyz & (peak_max1 < peak_max | peak_max < cut_off)){ 
  break 
} 
peaks <- c(peaks, peak_max) 
vibe <- rbind(vibe, flash_peak) 
### Increment the mid mbyz value 
peak_mid_grab <- peak_mid_grab + (1/z) 
peak_max1 <- peak_max 
} 

## ## Trace the starting point of 'frame2' 
######$ 
re_vibe <- data.frame() 
peak_max1 <- data$intensity[data$mass == bp_mbyz] 
a <- 1:1000 
b <- 1:1000 
peak_mid_grab <- bp_mbyz 
for(a in b){ 

### Grab the exact middle mbyz based on maximum intensity 
  ram <- which(data$mass >= peak_mid_grab) 
  rama <- data[ram[1],] 
  right <-rama[1] + (1/(2*z)) 
  left <- rama[1] - (1/(2*z)) 
  right_vect <- which(data$mass >= right[1,1]) 
  rit <- right_vect[1] 
  left_vect <- which (data$mass <= left[1,1]) 
  lft <- left_vect[length(left_vect)] interval <- data[lft:rit,] 
  
  kri <- max(interval$intensity) 
  bal <- which(interval$intensity == kri) 
  kal <- interval[bal,] 
  ppm_mid <- kal[1,1] 
  
  # Should stop if it extends beyond low mbyz of parent peptide 
  if(ppm_mid < ppm_exp){ 
    break
  }
  
  
  #### Trace the left half of distribution 
  position_mbyz <- which (data$mass >= ppm_mid) 
  click <- position_mbyz 
  i <- 1:1000 
  t <- 1:1000 
  for(i in t){
    del1 <- data$intensity[click] 
    del2 <- data$intensity[click - 1] 
    if(del2>del1){ 
      break 
    } 
    click <- click - 1 
  }
  
  
  #### Trace the right half of distribution 
  right_click <- position_mbyz 
  i <- 1:1000 
  t <- 1:1000 
  for(i in t){
    acer1 <- data$intensity[right_click] 
    acer2 <- data$intensity[right_click + 1] 
    if(acer2 > acer1){
      break 
    } 
    right_click <- right_click + 1 
  } 
  
  
  ### Subset the entire datapoints corresponding to the single peak 
  peak_distr <- data[click:right_click,] 
  #plot(peak_dist) 
  #windows() 
  ### Calculate maximum of peak intensity and check if it is increasing again 
  peak_max <- max(peak_distr$intensity) 
  re_vibe_mbyz <- peak_distr$mass[peak_distr$intensity == peak_max] 
  reverse_flash_peak <- data [data$mass == re_vibe_mbyz[1],] 
  
  if(peak_distr$mass[1] < bp_mbyz & (peak_max1 < peak_max | peak_max < cut_off)){ 
    if(peak_max == bp){ 
      peaksr <- c(peaksr, peak_max) 
      break 
    }
    break 
  }
  
  peaksr <- c(peaksr, peak_max) 
  re_vibe <- rbind(re_vibe, reverse_flash_peak) 
  ### Decrement the mid mbyz value 
  peak_mid_grab <- peak_mid_grab - (1/z) 
  peak_max1 <- peak_max 
  }

  
  ## Trace the end point of frame2 
  wall <- length(peaks) 
  fence <- vibe$mass[wall] 
  f2_end_pk <- which (data$mass == fence) 
  right_click <- f2_end_pk 
  i <- 1:1000 
  t <- 1:1000 
  for(i in t){
    acer1 <- data$intensity[right_click] 
    acer2 <- data$intensity[right_click + 1] 
    if(acer2 > acer1){
      break 
    } 
    right_click <- right_click + 1 
  }
  
  
  ## Trace the exact starting point for 'frame2' 
  wallr <- length(peaksr) 
  fencer <- re_vibe$mass[wallr] 
  f2_start_pk <- which(data$mass == fencer) 
  click <- f2_start_pk 
  i <- 1:1000 
  t <- 1:1000 
  for(i in t){ 
    del1 <- data$intensity[click] 
    del2 <- data$intensity[click -1] 
    if(del2>del1){ 
      break
    }
    click <- click - 1 
  }
  
  
  ######$ 
  
  frame2 <- data[click:right_click,] 
  
  # Do the fitting 
  # Assign the mbyz value of the first peak in frame2 to peak_mid_grab 
  peak_mid_grab <- data$mass[data$mass == fencer] 
  
  # Perform a loop to describe the peaks 
  a <- 1:1000 
  b <- 1:1000 
  peaks <- vector("numeric") 
  iso_peak <- vector("numeric") 
  iso_mass <- vector("numeric") 
  opt_sigma <- vector("numeric") 
  for(a in b){
    ram <- which(data$mass >= peak_mid_grab) 
    rama <- data[ram[1],] 
    right <-rama[1] + (1/(2*z)) 
    left <- rama[1] - (1/(2*z)) 
    right_vect <- which(data$mass >= right[1,1]) 
    rit <- right_vect[1] 
    left_vect <- which(data$mass <= left[1,1]) 
    lft <- left_vect[length(left_vect)] 
    interval <- data[lft:rit,] 
    
    kri <- max(interval$intensity) 
    bal <- which(interval$intensity == kri) 
    kal <- interval[bal,] 
    ppm_mid <- kal[1,1] 
    
    ## Trace the left half of distribution 
    position_mbyz <- which(data$mass == ppm_mid) 
    click <- position_mbyz 
    i <- 1:1000 
    t <- 1:1000 for(i in t){
      del1 <- data$intensity[click] 
      del2 <- data$intensity[click - 1] 
      if(del2>del1){ 
        break 
      } 
      click <- click - 1 
    }
    
    
    ## Trace the right half of distribution 
    right_click <- position_mbyz 
    i <- 1:1000 
    t <- 1:1000 
    for(i in t){
      acer1 <- data$intensity[right_click] 
      acer2 <- data$intensity[right_click + 1] 
      if(acer2 > acer1){ 
        break 
      } 
      right_click <- right_click + 1
    }
    
    
    ## Subset the entire datapoints corresponding to the single peak 
    peak_dist <- data[click:right_click,] 
    
    #plot(peak_dist) 
    #windows() 
    
    # Fit it in the gaussian curve to find the maximum value 
    tab1 <- data.frame(x=peak_dist$mass, r=peak_dist$intensity) 
    xmax <- tab1$x[tab1$r == max(tab1$r)]  
    r_at_fwtm <- 0.3 * (max(tab1$r)) 
    tab <- tab1[which(tab1$r >= r_at_fwtm),] 
    plot(tab$x, tab$r, xlab = "mass", ylab = "intensity") 
    b <- round(mean(tab$x)) 
    a <- max(tab$r) 
    xmax <- tab$x[tab$r == max(tab$r)] 
    x1 <- tab$x[tab$x < xmax][which.min(abs(tab$r[tab$x < xmax] - max(tab$r)/2))] 
    x2 <- tab$x[tab$x > xmax][which.min(abs(tab$r[tab$x > xmax] - max(tab$r)/2))] 
    c <- ((x2 - x1)/2.3) 
    print(c) 
    
    fg <- 1:1000 
    for(h in fg){
      (res <- nls( r ~ k*exp(-1/2*(x-mu)^2/sigma^2),start=c(mu=b,sigma=c,k=a) , data = tab, trace = TRUE, control = nls.control(printEval = TRUE, warnOnly = TRUE))) 
      c <- c + 0.1 
      if(res$convInfo$stopCode == 0){
        break 
        }
      } 
    print(summary(res)) 
    v <- summary(res)$parameters[,"Estimate"]
    
    
    plot(r~x, data=tab) 
    plot(function(x) v[3]*exp(-1/2*(x-v[1])^2/v[2]^2),col=2,add=T, xlim = range(tab$x)) 
    windows() 
    
    iso_mass1 <- summary(res)$parameters[1,"Estimate"] 
    iso_peak1 <- summary(res)$parameters[3,"Estimate"] 
    sigma1 <- summary(res)$parameters[2,"Estimate"] 
    
    iso_peak <- c(iso_peak, iso_peak1) 
    iso_mass <- c(iso_mass, iso_mass1) 
    opt_sigma <- c(opt_sigma, sigma1) 
    
    print(tab) 
    
    ## Calculate maximium of peak intensity and check if it is increasing again 
    peak_max <- max(peak_dist$intensity) 
    peaks <- c(peaks, peak_max) 
    
    ## Include the stop point corresponding to the extent of frame2 
    jolly <- nrow(frame2)  
    if(tab1$x[nrow(tab1)] == frame2$mass[jolly]){
      break
    }
    
    
    ## Increment the mid mbyz value 
    peak_mid_grab <- peak_mid_grab + (1/z) 
  }
  
    
    final <- data.frame(iso_mass, iso_peak) 
    names(final) <- c("mass","intensity") 
    
    # Include the virtual points to balance the envelope 
    ## Trace the left virtual point 
    tom <- which(frame2$mass<=final[1,1]) 
    hunter <- frame2[1:length(tom),] 
    ###select A point 
    AA <- which(hunter$intensity <= cut_off) 
    A <- hunter[length(AA),]$mass 
    A1 <- hunter[length(AA),]$intensity 
    ###select B point 
    B <- final[1,1] 
    B1 <- final[1,2] 
    print(A) 
    print(A1) 
    
    print(B) 
    print(B1) #
    ##find the left virtual point 
    walley <- approx(x=c(A1,B1), y=c(A,B), method="linear",xout=cut_off) 
    vir_lft_mbyz <- walley$y 
    vir_lft_intensity <- walley$x 
    if(!is.na(vir_lft_mbyz) & !is.na(vir_lft_intensity)){
      start_point_vir <- c("mass" = vir_lft_mbyz,"intensity" = vir_lft_intensity) 
      final <- rbind(start_point_vir,final) 
    }
    
    
    
    ## Trace the right virtual point 
    jerry <- which(frame2$mass>=final[nrow(final),1]) 
    thawne <- frame2[jerry[1]:jerry[length(jerry)],] 
    ###select A point
    A <- final[length(final$mass), 1]
    A1 <- final[length(final$mass), 2]
    ###select B point
    BB <- which(thawne$intensity <= cut_off)
    B <- thawne[BB[1],]$mass
    B1 <- thawne[BB[1],]$intensity 
    print(A) 
    print(A1)
    print(B)
    print(B1)
    ###find the right virtual point 
    logan <- approx(x=c(A1,B1), y=c(A,B), method="linear",xout=cut_off) 
    vir_rit_mbyz <- logan$y 
    vir_rit_intensity <- logan$x 
    if(!is.na(vir_rit_mbyz) & !is.na(vir_rit_intensity)){
      end_point_vir <- c("mass" = vir_rit_mbyz, "intensity" = vir_rit_intensity)
      final <- rbind(final,end_point_vir) 
    }
    
    
    
    View(final, title = "Exact peak values") 
    
    # Calculate centroid mass 
    mass_vect <- final$mass 
    intensity_vect <- final$intensity 
    denominator <- sum(intensity_vect) 
    numerator <- sum(mass_vect * intensity_vect)
    centroid_mass <- numerator / denominator 
    print("centroid_mass is") 
    print(centroid_mass) 
    
    

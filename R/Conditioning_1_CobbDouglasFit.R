#
# Dorleta Garcia
# 03/03/2011 14:48:33
# changed: 06/03/2011 21:18:44
#-------------------------------------------------------------------------------


CobbDouglasFit <- function(C, B, E, sig.level = 0.05){
    
    # C independent of B => beta = 0
    if( is.null(B)){
        print('Case 0: C independent of B -> beta = 0')
        y <- log(C/E)
        m0 <- lm(y ~ log(E))
        m0.aov  <- anova(m0)
        m0.coef <- unname(coef(m0))
        
        if(m0.coef[2] < -1 | m0.aov[1,5] > sig.level){
            print('Case 0.2: alpha non-significant OR nonsense => alpha = 1')
            m1 <-  lm(y ~ 1)
            m1.coef <- unname(coef(m1))
            return(list(params = c(q = exp(m1.coef[1]), alpha = 1, beta = 0), model = m1))
        }
        # else => alpha > -1 and significant.
       print('Case 0.1: alpha significant and > 0')
       return(list(params = c(q = exp(m0.coef[1]), alpha = m0.coef[2] + 1, beta = 0), model = m0))     
    }
    
    y <- log(C/(B*E))
    
    # If B is constant we don't have information to estimate 'beta' => beta = 1.
    if(var(B) == 0){
        print('Case 1: B constant -> beta = 1')
        m0 <- lm(y ~ log(E))
        m0.aov  <- anova(m0)
        m0.coef <- unname(coef(m0))
        
        if(m0.coef[2] < -1 | m0.aov[1,5] > sig.level){
            print('Case 1.2: alpha non-significant OR nonsense => alpha = 1')
            m1 <-  lm(y ~ 1)
            m1.coef <- unname(coef(m1))
            return(list(params = c(q = exp(m1.coef[1]), alpha = 1, beta = 1), model = m1))
        }
        # else => alpha > -1 and significant.
        print('Case 1.1: alpha significant and > 0')
       return(list(params = c(q = exp(m0.coef[1]), alpha = m0.coef[2] + 1, beta = 1), model = m0))     
    }
    
    # Apply lm
    print('Estimate alpha and beta')
    m0      <- lm(y ~ log(E) + log(B))
    m0.aov  <- anova(m0)
    m0.coef <- unname(coef(m0))

    # Any coefficient == NA??
    if(any(is.na(m0.coef[2:3]))){
        if(is.na(m0.coef[2])){ # => alpha = 1
            print('Case 2.1: alpha = NA -> alpha = 1')
            m0 <- lm(y ~ log(B))
            m0.aov  <- anova(m0)
            m0.coef <- unname(coef(m0))
        
            if(m0.coef[2] < -1 | m0.aov[1,5] > sig.level){
                print('Case 2.1.b: beta non-significant OR nonsense => beta = 1')
                m1 <-  lm(y ~ 1)
                m1.coef <- unname(coef(m1))
                return(list(params = c(q = exp(m1.coef[1]), alpha = 1, beta = 1), model = m1))
            }
            # else => beta > -1 and significant.
            print('Case 2.1.a: beta significant and > 0')
            return(list(params = c(q = exp(m0.coef[1]), alpha = 1, beta = m0.coef[2] + 1), model = m0))     
        }
        if(is.na(m0.coef[3])){ # => beta = 1
            print('Case 2.2: beta = NA -> beta = 1')
            m0 <- lm(y ~ log(E))
            m0.aov  <- anova(m0)
            m0.coef <- unname(coef(m0))
        
            if(m0.coef[2] < -1 | m0.aov[1,5] > sig.level){
                print('Case 2.2.b: alpha non-significant OR nonsense => alpha = 1')
                m1 <-  lm(y ~ 1)
                m1.coef <- unname(coef(m1))
                return(list(params = c(q = exp(m1.coef[1]), alpha = 1, beta = 1), model = m1))
            }
            # else => beta > -1 and significant.
            print('Case 2.2.a: alpha significant and > 0')
            return(list(params = c(q = exp(m0.coef[1]), alpha = m0.coef[2] + 1, beta = 1), model = m0))     
        }
    }
    
    # the coefficients are significant and sensible?
    if(all(m0.aov[,5][1:2] < sig.level) & all(m0.coef[2:3] > -1)){
        print('Case 3: alpha & beta are both significant and sensible')
        return(list(params = c(q = exp(m0.coef[1]), alpha = m0.coef[2] + 1, beta = m0.coef[3] + 1), model = m0))
    }
    
    # alpha AND/OR beta not significant BUT coefficients are sensible
    if(any(m0.aov[,5][1:2] > sig.level) & all(m0.coef[2:3] > -1)){
        sig <- c('log(E)', 'log(B)')[which(m0.aov[,5][1:2] <= sig.level)]
        sig <- ifelse(length(sig) == 0, 1,sig)
        
        if(sig == 1){  # coef non-sensible
            print('Case 4.1: alpha and beta are non-significant => alpha = beta = 1')
            m2 <- lm(y ~ 1) 
            return(list(params = c(q = exp(coef(m2)[1]), alpha = 1, beta = 1), model = m2))   
        }  
        
        form <- as.formula(paste('y ~', sig))
        m1      <- lm(form)
        m1.coef <- coef(m1)
        
        # coef in m1, wil be significant but is it sensible??
        if(m1.coef[2] < -1){  # coef non-sensible
            print('Case 4.1: alpha and beta are non-significant => alpha = beta = 1')
            m2 <- lm(y ~ 1) 
            return(list(params = c(q = exp(coef(m2)[1]), alpha = 1, beta = 1), model = m2))   
        }        
        
        q <- exp(m1.coef[1])
        alpha <- unname(ifelse('log(E)' %in% names(m1.coef), m1.coef['log(E)'] + 1, 1))
        beta  <- unname(ifelse('log(B)' %in% names(m1.coef), m1.coef['log(B)'] + 1, 1))
  
        print('Case 4: alpha OR beta are non-significant => alpha OR beta = 1')
        return(list(params = c(q = q, alpha = alpha, beta = beta), model = m1))
    }
    
    # => alpha OR beta non sensibles =< -1   (=> alpha OR beta = 0, base case)
    nosens <- which(m0.coef[-1] <= -1)

    if(length(nosens) == 2){ #=> alpha = beta = 1
        
        print('Case 5: In the initial fit both alpha AND beta are nonsense => fit the models for alpha and beta individually')
        
        m21    <- lm(y ~ log(B)) # alpha = 1
        m22    <- lm(y ~ log(E)) # beta = 1
        
        m21.coef <- unname(coef(m21))
        m22.coef <- unname(coef(m22))
        
        m21.aov <- anova(m21)
        m22.aov <- anova(m21)
        
       if(all(!((c(m21.coef[2], m22.coef[2]) > -1)*(c(m21.aov[,5][1], m22.aov[,5][1]) < sig.level)))){  # both non significant or nonsense
            m3 <- lm(y~1)
            print('Case 5.3: In the individual fits both are nonsense OR non-significant alpha = beta = 1') 
            return(list(params = c(q = exp(unname(coef(m3))), alpha = 1, beta = 1), model = m3))
        }
        else{ # alpha AND/OR beta > -1  and alpha AND/OR beta significant. (there is AT LEAST one that is at the same time is significant and sensible)
            if(all(((c(m21.coef[2], m22.coef[2]) > -1)*(c(m21.aov[,5][1], m22.aov[,5][1]) < sig.level)))){  # BOTH SENSIBLE AND SIGNIFICANT
                best <- which.min(AIC(m21, m22)[,2])
                print('Case 5.1: In the individual fits both are sensible and significant')
                
                if(best == 1){
                    print('Case 5.1.a: Biomass model is better, alpha = 1 & beta != 1') 
                    return(list(params = c(q = exp(m21.coef[1]), alpha = 1, beta = m21.coef[2] + 1), model = m21))
                }
                else {
                    print('Case 5.1.b: Effort model is better, alpha != 1 & beta = 1')  
                    return(list(params = c(q = exp(m22.coef[1]), alpha = m22.coef[2] + 1, beta = 1), model = m22))    # (best == 2)
                }
            }
            else{ # There is only one valid model
            print('Case 5.2: In the individual fits, there is only 1 valid model')
                if(m21.coef[2] < -1){   # m22.coef[2] > -1.
                    print('Case 5.2.a: alpha != 1, beta = 1')
                    return(list(params = c(q = exp(m22.coef[1]), alpha = m22.coef[2] + 1, beta = 1), model = m22))
                }
                else{      # m21.coef[2] > -1.
                    print('Case 5.2.b: alpha = 1, beta != 1')
                    return(list(params= c(q = exp(m21.coef[1]), alpha = 1, beta = m21.coef[2] + 1), model = m21))
                }
            }
        }

    }
    else{  # alpha _OR_ beta < -1.
  #  print(nosens)
        print('Case 6: In the initial fit there is 1 coefficient that is nonsense')
        if(nosens == 1){   # => alpha < -1   => _alpha = 1_
            m2   <- lm(y ~ log(B))
            m2.coef <- unname(coef(m2))
            m2.aov <- anova(m2)
            
            print('Case 6.1: alpha < 0 => alpha = 1')
            
            if(m2.coef[2] > -1 & m2.aov[1,5] <= sig.level){
                print('Case 6.1.a: beta sensible and significant')
                return(list(params = c(q = exp(m2.coef)[1], alpha = 1, beta = m2.coef[2] + 1), model = m2))
            }
            else{ # alpha nosense OR no significant.
                m3 <- lm(y~1)
                print('Case 6.1.b: beta nonsense OR non-significant => alpha = beta = 1')
                return(list(params = c(q = exp(unname(coef(m3))), alpha = 1, beta = 1), model = m3))
            }
        }
        else{ #  nosens == 2 # => beta < -1   => _beta = 1_
            m2   <- lm(y ~ log(E))
            m2.coef <- unname(coef(m2))
            m2.aov <- anova(m2)
            print('Case 6.2: beta < 0 => beta = 1')
            if(m2.coef[2] > -1 & m2.aov[1,5] <= sig.level){
                print('Case 6.2.a: alpha sensible and significant')
                return(list(params = c(q = exp(m2.coef)[1], alpha = m2.coef[2] + 1, beta = 1), model = m2) )
            }
            else{ # alpha nosense OR no significant.
                m3 <- lm(y~1)
                print('Case 6.1.b: alpha nonsense OR non-significant => alpha = beta = 1')
                return(list(params = c(q = exp(unname(coef(m3))), alpha = 1, beta = 1), model = m3))
            }
        }
    }
}
    


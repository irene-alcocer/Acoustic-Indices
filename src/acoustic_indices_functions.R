# Functions to prepare acoustic indices data

# n_used - N chosen for analysis
tidy_data <- function(df, n_used = "n_adjusted", rm_id = 54, 
                      new_dataset = TRUE,
                      independent_days = TRUE, rm_lt_4 = FALSE){
    
    # Select relevant columns
    df <- select(df, id = ID, Year, Impact_factor = Impact.factor, index = Acoustic.Indices, 
                 taxa = Taxonomic.group, environ = Environment, bio = Biodiversity.parameter, 
                 diversity_source = Diversity.source, 
                 n_original = sample.size, 
                 pseudoreplication = Pseudoreplication, 
                 pseudoreplication_type = Pseudoreplication.type, 
                 n_adjusted = Adjusted.sample.size, 
                 intermediate_group = Independence, Statistical_test, r, t)
    
    # Correct spearman correlations
    spearman <- select(df[df$Statistical_test == "spearman", ], r, .data[[n_used]])
    spearman <- spearman[complete.cases(spearman), ]
    spearman_to_correct <- spearman[spearman[[n_used]] < 90, ]$r
    spearman_correction <- 2 * sin((pi * spearman_to_correct)/6)
    df[which(df$Statistical_test == "spearman" &
                 df[[n_used]] < 90), ]$r <- spearman_correction
    
    # Remove rows without effect sizes
    df <- df[!(is.na(df$r) & is.na(df$t)), ]
    
    if(!rm_lt_4){
        # Put 4 in every n <= 3
        df[[n_used]][which(df[[n_used]] <= 3)] <- 4    
    }else{
        df <- df[- which(df[[n_used]] <= 3), ]
    }
    
    # Put 0.99 in r = 1
    df$r[which(df$r == 1)] <- 0.99
    
    # Compute effect sizes
    noNA_r <- df[!is.na(df$r), ]
    cor2z <- res(r = noNA_r$r, n = noNA_r[[n_used]], id = noNA_r$id, 
                 dig = 4, verbose = FALSE)[c("fisher.z", "var.z")]
    all(is.finite(cor2z$fisher.z)) & all(is.finite(cor2z$var.z))
    
    noNA_t <- df[!is.na(df$t), ]
    n <- floor(noNA_t[[n_used]] / 2)
    t2z <- tes(t = noNA_t$t, n.1 = n, n.2 = n, id = noNA_t$id, 
               dig = 4, verbose = FALSE)[c("fisher.z", "var.z")]
    all(is.finite(t2z$fisher.z)) & all(is.finite(t2z$var.z))
    # Add effect size column
    df[c("z", "var_z")] <- 0
    df[!is.na(df$r), ][c("z", "var_z")] <- cor2z
    df[!is.na(df$t), ][c("z", "var_z")] <- t2z

    # Remove specified rows
    if(!is.null(rm_id)){ 
        df <- df[df$id != rm_id, ] 
        cat("Removed study id\n\t", paste(rm_id, collapse = ", "), "\n")  
    }
    
    # Aggregate studies by intermediate group factor
    old_col <- nrow(df)
    df <- df %>% group_by(intermediate_group, year = Year, 
                          impact_factor = Impact_factor, index, taxa, environ, 
                          bio, id, diversity_source, 
                          pseudoreplication, n = get(n_used)) %>%
                 summarise(z = mean(z), var = composite_var(var_z)) %>%
                 as.data.frame() 
    df <- df[order(df$id), ]
    df <- cbind(entry = rownames(df), df)
    df <- select(df, entry, id, everything(), -intermediate_group) 
    cat(sprintf("Dataframe aggregated from %.0f to %.0f entries", old_col, nrow(df)),
        "\n")
    
    invisible(df)
    
}


clear_moderators <- function(df, 
			     mods = c("taxa", "bio", "diversity_source", "index"),
                             lower_nr = 5){
    
    cat("Levels dropped from dataframe:")
    for(mod in mods){
        if(!is.factor(df[[mod]])){
          df[[mod]] <- as.factor(df[[mod]])
        }
        cat(sprintf("\n\tModerator %s", mod, "\n\t"))
        nr_entries <- apply(table(df[[mod]], df$id), 1, 
                            function(x) length(which(x > 0)))
        less_5 <- names(which(nr_entries < lower_nr))
        df <- filter(df, !(.data[[mod]] %in% less_5))
        df[[mod]] <- droplevels(df[[mod]])
        cat("\n\t\t", less_5)
        
        # Relevel some factors
        if(mod == "diversity_source"){
            df[[mod]] <- relevel(df[[mod]], ref = "no_acoustic")
        }
        if(mod == "environ"){
            df[[mod]] <- relevel(df[[mod]], ref = "T")
        }
        if(mod == "bio"){
            df[[mod]] <- relevel(df[[mod]], ref = "richness")
        }
    }
    cat("\n\n")
    invisible(df)
}

# Calculate composite variance given a vector of variance.
# Formula according to Borenstein 2009 formula 24.5
composite_var <- function(vars, r = 0.5){
    # Do not calculate composite for scalar vectors
    if(length(vars) < 2) return(vars)
    # Calculate covariance with a correlation of r
    cov <- sum(combn(vars, 2, function(x) 2 * r * prod(sqrt(x))))
    # Correct by number of variance estimates
    composite_var <- (sum(vars) + cov) / length(vars)**2
    
    return(composite_var)
}


# Get number of studies per level combination of two factors.
# df must have an id column with study ids
n_studies_level <- function(df, factors){
    if(length(factors) != 2) stop("Exactly two factors must be given")
    # Select second factor to be column
    col_names <- unique(df[[ factors[2] ]]) 
    level_combn <- table(df$id, df[[ factors[1] ]], df[[ factors[2] ]])
    level_combn_single <- lapply(col_names, 
                                 function(x){
                                          colSums(
                                              ifelse(
                                                  level_combn[, , x] > 0, 1, 0)
                                              )
                                 })
                                
    level_combn <- do.call(cbind, level_combn_single)
    
    colnames(level_combn) <- col_names    
    
    return(level_combn)
}

# Get nr studies and nr entries for a given factor
n_studies_entries <- function(df, factor){
    nentries <- rowSums(table(df[[factor]], df$id))
    nstudies <- rowSums(ifelse(table(df[[factor]], df$id) > 0, 1, 0))
    return(cbind(Entries = nentries, Studies = nstudies))
}

# Transform Fisher Z to r
z2r <- function(z){
    num <- exp(2*z) - 1
    denom <- exp(2*z) + 1
    
    return(num/denom)
    
}

multilevel_I <- function(rma_mv_res){
    
    sampling_error <- (sum(1/rma_mv_res$vi) * (rma_mv_res$k-1)) / 
                        (sum(1/rma_mv_res$vi)^2 - sum((1/rma_mv_res$vi)^2))
    
    total_var <- sum(rma_mv_res$sigma2) + sampling_error
    
    # Unexplained variation
    Is <- rma_mv_res$sigma2 / total_var
    
    return(Is)
}

get_predictions <- function(res_rma, intercept = TRUE, format_table = FALSE,
                            clean_labels = FALSE){
    row_names <- rownames(res_rma$b)
    n_params <- nrow(res_rma$b)
    names_columns = c("coef", "estimate", "se", "ci.lb", "ci.ub")
    df_pred <- as.data.frame(matrix(nrow = n_params, 
                                    ncol = length(names_columns)))
    colnames(df_pred) <- names_columns
    for(i in 1:n_params){
        lev <- row_names[i]
        spot_lev <- vector(mode = "numeric", length = n_params)
        spot_lev[grep(lev, rownames(res_rma$b))] <- 1
        if(intercept){ spot_lev <- spot_lev[-1]}
        df_pred[i,] <- c(lev, unlist(predict(res_rma, spot_lev), 
                                     intercept)[1:4])
    }
    df_pred <- mutate_at(df_pred, 2:ncol(df_pred), as.numeric) %>%
                    mutate_if(is.numeric, z2r) %>%    
                    mutate_if(is.numeric, round, 3)
    df_pred$coef <- factor(df_pred$coef,
                           levels = rev(rownames(res_rma$b)))
    # Some formatting to print table
    if(!format_table){
        return(df_pred)
    }
    y_labels <- c("Intercept", "Index ADI", "Index AEI", "Index AR", 
                  "Index BIO", "Index H", "Index NDSI", 
                  "Bio abundance", "Bio diversity",
                  "Bio sound abundance", "Environment Aquatic", "Source acoustic")
    moderators <- str_match(y_labels, "(\\w*)")[, 2]
    moderators[1] <- "Intercept"
    # Get predictions
    df_pred$moderators <- factor(moderators, levels = unique(moderators))
    df_pred <- df_pred %>%
                    mutate(CI = paste0("[", ci.lb,"]", " ",
                                       "[", ci.ub, "]")) %>%
                    select(Moderators = moderators, 
                           Coefficients = coef,
                           Estimate = estimate,
                           SE = se, CI, -ci.lb, -ci.ub,)
    # Smaller labels to try new aesthetics
    small_labels <- str_match(y_labels, "\\w\\s(.*)")[, 2]
    # Capitalize first letter
    small_labels <- sapply(small_labels, 
                           function(x){
                               caps <- toupper(x)
                               caps <- str_sub(caps, 1, 1)
                               new_x <- paste0(caps, str_sub(x, 2))
                           })
    small_labels[1] <- "Intercept"
    names(small_labels) <- NULL
    if(clean_labels){
        df_pred$Coefficients <- small_labels    
    }
    
    return(df_pred)
}

compare_moderators <- function(df, res, moderator, ref_level){
    levs <- paste0(moderator, levels(df[[moderator]])[-1])
    comparisons_cols <- c("Compared", "Estimate", "SE", "CI.lb", "CI.up", "QM", "p")
    df_comparisons <- matrix(nrow = length(levs), ncol = 7,
                               dimnames = list(NULL, comparisons_cols))
    df_comparisons <- as.data.frame(df_comparisons)
    for(i in 1:length(levs)){
        lev2 <- paste0(moderator, ref_level)
        test_levs <- vector(mode = "numeric", length = nrow(res$b))
        test_levs[grep(levs[i], rownames(res$b))] <- -1
        test_levs[grep(lev2, rownames(res$b))] <- 1
        comp <- anova(res, L = test_levs)
        CI.lb <- comp$Xb[,1] - 1.96 * comp$se
        CI.ub <- comp$Xb[,1] + 1.96 * comp$se
        lev1 <- ifelse(levs[i] == lev2, levels(df[[moderator]])[1], levs[i])
        comparison <- str_remove_all(paste0(lev2, " - ", lev1), moderator)
        values <- c(comparison, comp$Xb[,1], comp$se, 
                    CI.lb, CI.ub, comp$QM, comp$QMp)
        df_comparisons[i, ] <- values
        
    }
    df_comparisons <- df_comparisons %>%
                        mutate_at(2:ncol(df_comparisons), as.numeric) 
    
    return(df_comparisons)
}


# This function was taken from dmetar package, as it does a nice plot of I2
# It includes a few minor modifications
mlm.variance.distribution <- function(x){
    
    m = x
    
    # Check class
    if (!(class(m)[1] %in% c("rma.mv", "rma"))){
        stop("x must be of class 'rma.mv'.")
    }
    
    # Check for three level model
    if (m$sigma2s != 2){
        stop("The model you provided does not seem to be a three-level model. This function can only be used for three-level models.")
    }
    
    
    # Get variance diagonal and calculate total variance
    n = m$k.eff
    vector.inv.var = 1/(diag(m$V))
    sum.inv.var = sum(vector.inv.var)
    sum.sq.inv.var = (sum.inv.var)^2
    vector.inv.var.sq = 1/(diag(m$V)^2)
    sum.inv.var.sq = sum(vector.inv.var.sq)
    num = (n-1)*sum.inv.var
    den = sum.sq.inv.var - sum.inv.var.sq
    est.samp.var = num/den
    
    # Calculate variance proportions
    level1=((est.samp.var)/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
    level2=((m$sigma2[1])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
    level3=((m$sigma2[2])/(m$sigma2[1]+m$sigma2[2]+est.samp.var)*100)
    
    # Prepare df for return
    Level=c("Level 1", "Level 2", "Level 3")
    Variance=c(level1, level2, level3)
    df.res=data.frame(Variance)
    colnames(df.res) = c("% of total variance")
    rownames(df.res) = Level
    I2 = c("---", round(Variance[2:3], 2))
    df.res = as.data.frame(cbind(df.res, I2))
    
    totalI2 = Variance[2] + Variance[3]
    
    
    # Generate plot
    df1 = data.frame("Level" = c("Sampling Error", "Total Heterogeneity"),
                     "Variance" = c(df.res[1,1], df.res[2,1]+df.res[3,1]),
                     "Type" = rep(1,2))
    
    df2 = data.frame("Level" = rownames(df.res),
                     "Variance" = df.res[,1],
                     "Type" = rep(2,3))
    
    df = as.data.frame(rbind(df1, df2))
    
    
    g = ggplot(df, aes(fill=Level, y=Variance, x=as.factor(Type))) +
        coord_cartesian(ylim = c(0,1), clip = "off") +
        geom_bar(stat="identity", position="fill", width = 1, color="black") +
        scale_y_continuous(labels = scales::percent)+
        theme(axis.title.x=element_blank(),
              axis.text.y = element_text(color="black"),
              axis.line.y = element_blank(),
              axis.title.y=element_blank(),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.y = element_line(lineend = "round"),
              legend.position = "none",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              legend.background = element_rect(linetype="solid",
                                               colour ="black"),
              legend.title = element_blank(),
              legend.key.size = unit(0.75,"cm"),
              axis.ticks.length=unit(.25, "cm"),
              plot.margin = unit(c(1,3,1,1), "lines")) +
        scale_fill_manual(values = c("deepskyblue1", "darkseagreen1", "darkseagreen2",
                                     "deepskyblue2", "darkseagreen3")) +
        
        # Add Annotation
        
        # Total Variance
        annotate("text", x = 1.5, y = 1.05,
                 label = paste("Total Variance:",
                               round(m$sigma2[1]+m$sigma2[2]+est.samp.var, 3))) +
        
        # Sampling Error
        annotate("text", x = 1, y = (df[1,2]/2+df[2,2])/100,
                 label = paste("Sampling Error Variance: \n", round(est.samp.var, 3)), size = 3) +
        
        # Total I2
        annotate("text", x = 1, y = ((df[2,2])/100)/2-0.02,
                 label = bquote("Total"~italic(I)^2*":"~.(round(df[2,2],2))*"%"), size = 3) +
        annotate("text", x = 1, y = ((df[2,2])/100)/2+0.05,
                 label = paste("Variance not attributable \n to sampling error: ", round(m$sigma2[1]+m$sigma2[2],3)), size = 3) +
        
        # Level 1
        annotate("text", x = 2, y = (df[1,2]/2+df[2,2])/100, label = paste("Level 1: \n",
                                                                           round(df$Variance[3],2), "%", sep=""), size = 3) +
        
        # Level 2
        annotate("text", x = 2, y = (df[5,2]+(df[4,2]/2))/100,
                 label = bquote(italic(I)[Level2]^2*":"~.(round(df[4,2],2))*"%"), size = 3) +
        
        # Level 3
        annotate("text", x = 2, y = (df[5,2]/2)/100,
                 label = bquote(italic(I)[Level3]^2*":"~.(round(df[5,2],2))*"%"), size = 3)

    suppressWarnings(print(g))
    invisible(df.res)
}

final_mod_names <- function(names_vec, shorter = FALSE){
    if(shorter){
        sp <- "Sp."
        abund <- "Abund."
        diver <- "Div."
    }else{
        sp <- "Species"
        abund <- "Abundance of"
        diver <- "Diversity of"
    }
    names_vec <- names_vec %>% 
        str_replace("^(a|A)bundance", paste(sp, "abundance")) %>%
        str_replace("^(d|D)iversity", paste(sp, "diversity")) %>%
        str_replace("^(s|S)ound(_| )abundance", paste(abund, "sounds")) %>%
        str_replace("^(s|S)ound(_| )richness", paste(diver, "sounds")) %>%
        str_replace("A$", "Aquatic") %>%
        str_replace("^acoustic", "Acoustic") %>%
        str_replace("no_acoustic", "Non acoustic") %>%

    return(names_vec)
}

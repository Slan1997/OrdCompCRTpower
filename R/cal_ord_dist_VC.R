library(pacman)
pacman::p_load(tidyr, dplyr,readr,stringr,ggplot2,
               ggpubr,gridExtra,scales,RColorBrewer,
               magrittr,viridis,
               VGAM, # for pbinorm
               VineCopula ,# for BiCopTau2Par
               network, # for plot.RVineMatrix
               polycor
               )

######################################## following from VC_sim_functions.R
odds = function(x){ x/(1-x)}
get_p_trt = function(beta_X,p_ct){
  exp(beta_X)*odds(p_ct)/( exp(beta_X)*odds(p_ct) + 1  )
}
########################################

######################################## following from VCOrdinal_PMF_functions518.R
## get all the setup variables
VCO_Initialize = function(nComp,prev,
                          vine_type = 'D', # plan to accommodate C, D, and maybe R vine
                          method_list, # for later use
                          rho_list=NULL,
                          param_list = NULL){
  # initialize
  nTree = nComp-1
  qs = 1-prev
  
  # eta_list if we only know rho and need to find eta from rho
  if (!is.null(rho_list)){
    eta_list = vector('list',nTree)
  }else{eta_list = NULL}
  
  if (vine_type == 'D'){   # assume the variables are already in order
    # pair list
    pair_list = vector('list',nTree)
    names(pair_list) = paste0('T',2:nComp)
    for (pr in 1:nTree){
      if (pr==1){
        pair_list[[pr]] = matrix(c(1:nTree,2:nComp),byrow=T,nrow=2)
      }else{
        pair_list[[pr]] = matrix(c(1:(2*(nComp-pr))),nrow=2)
      }
    }
    
    # conditioning_joint_node 
    conditioning_joint_node = vector('list',nTree)
    names(conditioning_joint_node) = paste0('T',2:nComp)
    for (l in 1:nTree){
      tree_idx = l+1
      temp = matrix(NA,nrow=tree_idx,ncol=(nComp-l))
      for (r in 1:tree_idx){
        temp[r,] = r:(nTree-l+r)
      }
      cont_num = apply(temp,2,function(x) paste0(x,collapse=''))
      n_c = length(cont_num)
      if (n_c>2){
        final_num = c(cont_num[1],
                      rep(cont_num[2:(n_c-1)],each=2),
                      cont_num[n_c])
      }else{
        final_num = cont_num
      }
      conditioning_joint_node[[l]] = paste0('V',final_num)
    }
    
    # conditioned_on_nodes
    conditioned_on_nodes = vector('list',(nTree-1))
    names(conditioned_on_nodes) = paste0('T',2:(nComp-1))
    # joint_nodes_needed
    joint_nodes_needed = vector('list',(nTree-1))
    names(joint_nodes_needed) = paste0('T',3:nComp)
    for (l in 1:(nTree-1)){
      tree_idx = l+1
      temp = matrix(NA,nrow=l,ncol=(nComp-tree_idx))
      for (r in 1:l){
        temp[r,] = (r+1):(nTree-l+r)
      }
      cont_num = apply(temp,2,function(x) paste0(x,collapse=''))
      joint_nodes_needed[[l]] = paste0('V',cont_num)
      conditioned_on_nodes[[l]] = paste0('V',rep(cont_num,each=2))
    }
    
    
  }else if (vine_type == 'C'){
    # pair list 
    pair_list = vector('list',nTree)
    names(pair_list) = paste0('T',2:nComp)
    
    for (pr in 1:nTree){
      pair_list[[pr]] = matrix(c(rep(1,(nComp-pr)),2:(nComp-pr+1)),byrow=T,nrow=2)
    }
    
    # conditioning_joint_node 
    conditioning_joint_node = vector('list',nTree)
    names(conditioning_joint_node) = paste0('T',2:nComp)
    for (l in 1:nTree){
      conditioning_joint_node[[l]] = paste0('V',paste0(paste0(1:l,collapse=''), (l+1):nComp))
    }
    #conditioning_joint_node
    
    # conditioned_on_nodes
    conditioned_on_nodes = vector('list',(nTree-1))
    names(conditioned_on_nodes) = paste0('T',2:(nComp-1))
    # joint_nodes_needed
    joint_nodes_needed = vector('list',(nTree-1))
    names(joint_nodes_needed) = paste0('T',3:nComp)
    for (l in 1:(nTree-1)){
      conditioned_on_nodes[[l]] = paste0('V',rep( paste0(1:l,collapse=''), nComp-l))
      joint_nodes_needed[[l]] = paste0('V',rep( paste0(1:l,collapse=''), nComp-l-1))
    }
  }
  
  # list_joint_pmf
  V_single = matrix(c(qs,prev),nrow=2,byrow = T) 
  colnames(V_single) = paste0('V',1:nComp)
  rownames(V_single) = paste0('X',0:1)
  
  list_joint_pmf = list(V_single)
  for (ti in 2:(nTree+1)){
    Tidx = paste0('T',ti)
    joints = unique(conditioning_joint_node[[Tidx]])
    temp = matrix(NA,nrow=2^ti,ncol=length(joints))
    colnames(temp) = joints
    rownames(temp) = sort(paste0('X',expand.grid(rep(list(0:1),ti)) %>% 
                                   apply(1,function(x) paste0(x,collapse=''))))
    list_joint_pmf = c(list_joint_pmf,list(temp))
  }
  names(list_joint_pmf) = paste0('T',1:nComp)
  
  return(list(nComp=nComp,
              nTree=nTree,
              prev=prev,
              qs=qs,
              rho_list=rho_list,
              param_list=param_list,
              eta_list=eta_list,
              method_list=method_list,
              pair_list=pair_list,
              conditioning_joint_node=conditioning_joint_node,
              conditioned_on_nodes=conditioned_on_nodes,
              joint_nodes_needed=joint_nodes_needed,
              list_joint_pmf =list_joint_pmf))
}


# Functions for bivariate copula: 
## now have 5 different methods: Independent, Gaussian, Clayton, Frank, Gumbel.
## Note: Independent equals to Gaussian with 0 correlation.
## A good website to visualize copula contour: 
# https://copulatheque.shinyapps.io/copulas/ 

copula_family_map <- c(
  Indep    = 0,
  Gaussian = 1,
  Clayton  = 3,
  Gumbel   = 4,
  Frank    = 5,
  Joe      = 6
)

# name_to_family <- function(method) {
#   if (!method %in% names(copula_family_map)) {
#     stop("Unsupported copula method: ", method)
#   }
#   unname(copula_family_map[[method]])
# }

name_to_family <- function(method) {
  family_map <- c(Indep = 0, Gaussian = 1, Clayton = 3, 
                  Gumbel = 4, Frank = 5, Joe = 6)
  return(unname(family_map[method]))
}

### ### C_fun is a wrapper for VineCopula::BiCopCDF
C_fun <- function(u1, u2,
                  specific_method = "Indep",
                  parC = NULL){
  
  family <- copula_family_map[[specific_method]]
  
  if (is.null(family)) {
    stop(paste0("Unsupported copula method: ", specific_method))
  }
  
  if (family != 0 && is.null(parC)) {
    stop(
      paste0(
        "Missing parameter. Copula '", specific_method,
        "' requires dependence parameter."
      )
    )
  }
  
  VineCopula::BiCopCDF(
    u1 = u1,
    u2 = u2,
    family = family,
    par = parC
  )
}

# C_fun(.8, .9, "Indep")
# 
# C_fun(.8, .9, "Gaussian", .3)
# 
# C_fun(.8, .9, "Clayton", 1)
# 
# C_fun(.8, .9, "Frank", 3)
# 
# C_fun(.8, .9, "Gumbel", 1.5)
# 
# C_fun(.8, .9, "Joe", 1.5)

## compare to VineCopula::BiCopCDF, works the same. Can use BiCopCDF when need more options.
# https://github.com/tnagler/VineCopula/blob/main/R/BiCopCDF.R
# library(VineCopula)
# BiCopCDF(u1=.8,u2=.9,family=0,par=NULL)
# BiCopCDF(u1=.8,u2=.9,family=1,par=.1)
# BiCopCDF(u1=.8,u2=.9,family=3,par=1)
# BiCopCDF(u1=.8,u2=.9,family=5,par=3)
# BiCopCDF(u1=.8,u2=.9,family=4,par=1.5)


### This function is now only valid for Gaussian copula, not sure about the other method yet.
# for other method, we can just assume parameter first.
### if we know the correlation, we can find parameter for copula
find_parC = function(p1,p2,rho,specific_method
){
  # print(paste('find_parC:',specific_method))
  q1 = 1-p1
  q2 = 1-p2
  ## target C(p1,p2,parC), need to find parC by using eq.(2) from lin2021
  target = rho*sqrt(p1*q1*p2*q2)+p1*p2   
  #if (specific_method=='Gaussian'){
  parC_list = (1:1e3)*.001
  #}
  try_get_target = sapply(parC_list, 
                          function(x) C_fun(p1,p2,
                                            specific_method=specific_method,
                                            parC=x)) 
  parC = parC_list[which.min(abs(try_get_target-target ))]
  if (parC == 1) parC = .99
  parC
}

# # pair copula for two binary components
# BiCop_PMF = function(p1,p2,parC,specific_method,print=T){
#   parC0=parC
#   # if (print) print(specific_method)
#   # if (print&(specific_method!='Indep')) print(paste0('parameter input is: ',parC0))
#   biPMF = -1
#   q1 = 1-p1
#   q2 = 1-p2
#   ct = 0
#   while( ! all(biPMF>=0) & ct<10){
#     ct=ct+1
#     key = C_fun(u1=q1,u2=q2,specific_method=specific_method,parC=parC) # parC=parC or find_parC(p1,p2,rho)
#     #print(key)
#     biPMF = c(key,q1-key, q2-key,1-q1-q2+key)
#     #print(biPMF)
#     if (! all(biPMF>=0) ){
#       parC = parC-.001
#       print(biPMF)
#     }
#   }
#   if (! all(biPMF>=0) ) print('biPMF has term less than 0!')
#   ## sometimes the output pmf have some term less than 0
#   if (print & (parC!=parC0)&(specific_method!='Indep')) print(paste0('parameter being used is: ',parC))
#   biPMF
# }

# pair copula for two binary components
check_param <- function(method, par) {
  if (method == "Indep") return(TRUE)
  if (method == "Gaussian") return(abs(par) < 1)
  if (method == "Gumbel" && par < 1) return(FALSE)
  if (method == "Joe" && par < 1) return(FALSE)
  if (method == "Clayton" && par <= 0) return(FALSE)
  TRUE
}

BiCop_PMF <- function(p1, p2, parC, specific_method, print = FALSE) {
  
  if (print) print(specific_method)
  if (print && specific_method != "Indep") {
    print(paste0("parameter input is: ", parC))
  }
  
  q1 <- 1 - p1
  q2 <- 1 - p2
  
  if (!check_param(specific_method, parC)) {
    message(
      "Invalid parameter detected:",
      "\nFamily: ", specific_method,
      "\nParameter: ", parC,
      "\np1: ", p1,
      "\np2: ", p2
    )
    
    stop(
      paste0(
        "Parameter ", round(parC,4),
        " is not valid for copula ", specific_method
      ),
      call. = FALSE
    )
  }
  
  key <- C_fun(
    u1 = q1,
    u2 = q2,
    specific_method = specific_method,
    parC = parC
  )
  
  biPMF <- c(
    key,
    q1 - key,
    q2 - key,
    1 - q1 - q2 + key
  )
  
  tol <- 1e-10
  
  if (all(biPMF >= -tol)) {
    biPMF[biPMF < 0] <- 0
    biPMF <- biPMF / sum(biPMF)
    return(biPMF)
  }
  
  if (print) print(biPMF)
  
  stop(
    paste0(
      "Invalid bivariate PMF for copula family ", specific_method,
      " with parameter ", round(parC,4),
      ". Try another family or weaker dependence."
    ),
    call. = FALSE
  )
}

# for now: rho is only accepted for Gaussian copula
get_temp_bi_conditional_pmf = function(prevalences,pair,
                                       rho=NULL, # should be a vector, need rho for each pair
                                       # currently rho is same for all pair, so it is only a numeric value.
                                       param1=NULL, # should be a vector, need param for each pair
                                       method1 # should be a vector
){
  
  param0 = param1
  temp_bi_conditional_pmf = matrix(NA,nrow=4,ncol=ncol(pair))
  colnames(temp_bi_conditional_pmf) = apply(pair,2,function(x)
    paste0('V',paste0(x,collapse='')))
  #if (method=='Gaussian'){ # if not independent, either need parameter or rho.
  #if (is.null(param)){
  # stopifnot(!is.null(rho))
  # if only know rho, but don't know the parameter for copula, need to find the parameter by eq.2)
  eta_vec = rep(NA,ncol(pair))
  if (!is.null(rho) & length(rho) == 1) rho = rep(rho,ncol(pair))
  if (is.null(param0)) param1 =rep(NA,ncol(pair))
  if (!is.null(param0) &length(param0) == 1) param1 = rep(param1,ncol(pair))
  #print(param1)
  if (length(method1) == 1) method1 = rep(method1,ncol(pair))
  #print(method1)
  #print(eta_vec)
  #   }else{
  #     eta_vec = param
  #   }
  # }else{
  #   eta_vec = param
  # }
  #
  for (i in 1:ncol(pair)){
    bi_idx = paste0(c('V',pair[,i]),collapse = '')
    v1_idx = pair[1,i]
    v2_idx = pair[2,i]
    p1 = as.numeric(prevalences[v1_idx])
    p2 = as.numeric(prevalences[v2_idx])
    if (method1[i]=='Gaussian'&is.null(param0[i])){
      ## may change later when we have param not no rho
      #print(rho[i])
      eta_vec[i] = find_parC(p1,p2,rho[i],
                             specific_method=method1[i])
      param1[i] = eta_vec[i]
      #print(paste0("param1:",param1))
    }
    #bi_pmf = BiCop_PMF(p1,p2,parC=eta_vec[i])
    temp_bi_conditional_pmf[,bi_idx] =
      BiCop_PMF(p1,p2,parC=param1[i],specific_method=method1[i])
  }
  
  #print(eta_vec)
  return(list(temp_bi_conditional_pmf=temp_bi_conditional_pmf,
              eta_vec =eta_vec 
  ))
}

gen_jointpmf = function(SetupInfo,...){
  # Initialize
  nComp = SetupInfo$nComp
  nTree = SetupInfo$nTree
  prevalences = SetupInfo$prev
  rho_list = SetupInfo$rho_list
  param_list = SetupInfo$param_list
  eta_list = SetupInfo$eta_list
  method_list = SetupInfo$method_list
  pair_list = SetupInfo$pair_list
  conditioning_joint_node = SetupInfo$conditioning_joint_node
  conditioned_on_nodes = SetupInfo$conditioned_on_nodes
  joint_nodes_needed = SetupInfo$joint_nodes_needed
  list_joint_pmf = SetupInfo$list_joint_pmf 
  
  for (idx in 1:nTree){
    # print(idx)
    if (!is.null(rho_list)){
      if (is.list(rho_list)){
        rho = rho_list[[idx]]
      }else{ 
        ## rho_list can also be a vector if assume correlation are the same for each layer.
        rho = rho_list[idx] # get correlation if there is any input
      }
    }else{rho = NULL}
    #when idx=1, rho: 0.2 0.1 0.3 0.1 0.2
    
    if (is.list(method_list)){
      method = method_list[[idx]]
    }else{
      ### method_list can also be a vector if assume pair copula method are the same for each layer.
      method = method_list[idx] 
    }
    
    
    if (!is.null(param_list)){
      if (is.list(param_list)){
        param = param_list[[idx]] 
      }else{
        ## param_list can also be a vector if assume parameters are the same for each layer's bivariate copula models.
        param = param_list[idx] # get parameter if there is any input
      }
    }else{
      param = NULL
    }
    
    T_idx = paste0('T',idx+1) # T2, index for tree in text
    T_idx_prev = paste0('T',idx) 
    # get bivariate pmd from each pair
    pair = pair_list[[T_idx]]
    
    if (idx==1){
      #### get bivariate copula pmf (3/28: updated simulation/cluster_VCO/VCO/VCOrdinal_PMF_functions.R get_temp_bi_conditional_pmf to handle vector rho)
      get_bi_result = get_temp_bi_conditional_pmf(prevalences=prevalences,
                                                  pair=pair,
                                                  rho=rho,
                                                  param1=param,
                                                  method1=method)
      
      temp_bi_conditional_pmf = get_bi_result[[1]]
      rownames(temp_bi_conditional_pmf) = paste0(rep(c('X0','X1'),each=2),c('0','1'))
      if (is.null(param)& all(method=='Gaussian')) eta_list[[idx]] = get_bi_result[[2]]
      ## output1: joint pmf, no need to compute, directly extract
      list_joint_pmf[[T_idx]] = temp_bi_conditional_pmf
    }else{
      # prev_row indicate which conditioned prob we are using
      for (prev_row in 1:nrow(prevalences)){ #### for each conditioned value, need to get the joint pmf
        #### get bivariate copula pmf
        get_bi_result = get_temp_bi_conditional_pmf(prevalences[prev_row,],
                                                    pair=pair,
                                                    rho=rho,
                                                    param1=param,
                                                    method1=method)
        temp_bi_conditional_pmf = get_bi_result[[1]]
        rownames(temp_bi_conditional_pmf) = paste0(rep(c('X0','X1'),each=2),c('0','1'))
        if (is.null(param)& all(method=='Gaussian')) eta_list[[idx]] = rbind(eta_list[[idx]],get_bi_result[[2]])
        
        ################## following part is different from Tree 2: need to compute joint pmf now
        ## output1: joint pmf 
        # get the variable conditioned on, needed it for getting joint distribution.
        cond_Vs =joint_nodes_needed[[T_idx]]
        cond_prob = list_joint_pmf[[(idx-1)]][prev_row,cond_Vs] # prev_row corresponding to which numbers are conditioned on
        joint_prob = temp_bi_conditional_pmf
        joint_name = colnames( list_joint_pmf[[T_idx]])
        cond_on = gsub('X','',rownames(prevalences)[prev_row])
        
        for (k in 1:length(joint_name)){
          joint_var = joint_name[k]
          cond_V = gsub('V','',cond_Vs[k])
          # get the position of cond_V
          ## cond_V can be more than 1 digit, e.g., V23, but it has to be in ascending order
          cond_V_digits = strsplit(cond_V,split='')[[1]]
          
          pos = as.vector(sapply(cond_V_digits,function(x) grep(x,strsplit(joint_var,split='')[[1]])))
          
          # get the row for cond_V being 0
          target_rows = sapply(rownames(list_joint_pmf[[T_idx]]),
                               function(x) paste0(strsplit(x,split='')[[1]][pos],
                                                  collapse = '') == cond_on)
          list_joint_pmf[[T_idx]][target_rows,joint_var] = 
            joint_prob[,k]*cond_prob[k]
        }
      }
    }
    
    if (idx != nTree){ # don't need these part for last idx.
      ## output2: get conditional variables' prevalence (new variable, reindexed) to prepare for next Tree 
      ## grab conditioning joint nodes
      cjn = conditioning_joint_node[[T_idx]]
      # e.g. Tree 3: [1]  "V123" "V234" "V234" "V345" "V345" "V456"
      ## grab conditioned nodes
      cn = conditioned_on_nodes[[T_idx]]
      # e.g. Tree 3: [1] "V23" "V23" "V34" "V34" "V45" "V45"
      len = length(cn)
      
      ## find out the position of the digit of non-conded_prob
      # e.g.,position of digit 1 in 12, position of digit 3 in 23, position of digit 2 in 23
      separate_cjn = sapply(gsub('V','',cjn),function(x) strsplit(x,split='')[[1]])
      separate_cn = sapply(gsub('V','',cn),function(x) strsplit(x,split='')[[1]])
      non_conded_digit = rep(NA,len)
      non_conded_digit_pos = rep(NA,len)
      for (s in 1:len ){
        if (is.null(dim(separate_cn))) sub = separate_cn[s]
        else sub = separate_cn[,s]
        non_conded_digit[s] = separate_cjn[!separate_cjn[,s] %in%sub,s]
        non_conded_digit_pos[s] = which(separate_cjn[,s] == non_conded_digit[s])
      }
      # non_conded_digit_pos  # need to let this digit to be 1 in the joint dist when find the conditional pmf
      
      prevalences = matrix(NA,ncol=len,nrow=2^idx)
      colnames(prevalences) = paste(cjn,cn,sep='|')
      rownames(prevalences) = sort(paste0('X',expand.grid(rep(list(0:1),idx)) %>% 
                                            apply(1,function(x) paste0(x,collapse=''))))
      # prevalence = p(V123=1|V23=rowname)
      # prevalences
      
      # for conditioning on 00,01,10,11 (T3)
      for (r in rownames(prevalences)){
        # r =  rownames(prevalences)[2] 
        v_names = matrix(NA,nrow=idx+1,ncol=len)
        for (v in 1:len){
          v_names[non_conded_digit_pos[v],v]='1' # make sure the non_conded_digit always be 1
          v_names[-non_conded_digit_pos[v],v] = strsplit(gsub('X','',r), 
                                                         split='')[[1]]
          # make sure the conditioned digits corresponding to rowname r e.g. X01
        }
        # v_names 
        # target joint PMF we needed for each conditional prevalence
        rjoint = paste0('X',apply( v_names,2,function(x) paste0(x,collapse='')))
        
        ## extract the corresponding conditioned probability
        conded_prob = list_joint_pmf[[T_idx_prev]][r,cn] 
        ## extract the corresponding joint conditioning probability
        conding_prob = apply(rbind(rjoint,cjn),2,
                             function(x) list_joint_pmf[[T_idx]][x[1],x[2]])
        prevalences[r,] = conding_prob /conded_prob 
      }
      # print(prevalences)
    }else{ break }
  }
  return(list(list_joint_pmf=list_joint_pmf,
              eta_list=eta_list,
              jointpmf = list_joint_pmf[[length(list_joint_pmf)]]))
}




Ordinal_PMF = function(jointpmf,method,severity_score_increasing=NULL){
  nComp = log2(nrow(jointpmf)) # find how many components in the model.
  entry = gsub('X','',rownames(jointpmf) )
  if (method=='sum'){
    nc = nchar(gsub('0','',entry )) 
    ## use string pattern to create PMF of ordinal variable.
    Ordinal_sum_PMF = matrix(NA,ncol=(nComp+1),nrow=1)
    colnames(Ordinal_sum_PMF) = paste0('X',0:nComp)
    Ordinal_sum_PMF[1,] = sapply(0:nComp,function(x) sum(jointpmf[nc==x]))
    return(Ordinal_sum_PMF)
  }else{
    # when method=='order'
    # check if severity_score_increasing has been specified
    stopifnot(is.vector(severity_score_increasing))
    entry_new = sapply(entry,function(x) 
      paste0(strsplit(x,split='')[[1]][severity_score_increasing],
             collapse = ''))
    # unless severity_score_increasing = 1:6, otherwise need to reorder the digits
    # after reorder, we are basically assigning severity score 1-6 to component 1-6
    jointpmf_cp = jointpmf
    rownames(jointpmf_cp) = entry_new
    
    Ordinal_ranking_PMF = matrix(NA,ncol=(nComp+1),nrow=1)
    colnames(Ordinal_ranking_PMF) = paste0('X',0:nComp)
    for (j in 0:nComp){
      if (j==0){
        Ordinal_ranking_PMF[,'X0'] = jointpmf_cp[paste0(rep(0,nComp),collapse = ''),]
      }else{
        # put 1 at position j, and 0's for positions > j (there are n-j 0's)
        items = entry_new[grep(paste0(1,paste0(rep(0,nComp-j),collapse = ''),'$',
                                      collapse = ''),entry_new)] 
        # e.g., when j=2, items = "010000", "110000"
        Ordinal_ranking_PMF[,paste0('X',j,collapse = '')] = sum(jointpmf_cp[items,])
      }
    }
    return(Ordinal_ranking_PMF)
  }
}

# Define the choices globally for easier access
copula_choices <- c("Indep" = "Indep", "Gaussian" = "Gaussian", "Clayton" = "Clayton", 
                    "Gumbel" = "Gumbel", "Frank" = "Frank", "Joe" = "Joe")

corr_choices <- c("High (tau=0.6)" = "High", "Moderate (tau=0.4)" = "Moderate", 
                  "Low (tau=0.2)" = "Low", "Very low (tau=0.05)" = "Very low", 
                  "Indep (tau=0)" = "Indep")

# ##################################################################
# 
# ##### Binary components: e.g., adverse events
# #### Calculate the joint distribution given VC configuration
# ## input:
# # 1. number of binary components
# # 2. [control prevalence & treatment prevalence for each] or [control or treatment or overall averaged prevalence for each components + Odds ratio for each components]
# ### maybe add these to parameters output?
# 
# ## choice of vine copula structure: 
# ## default D
# # if just 3 components, then don't need to choose
# # if > 3 components, C or D 
# 
# ## correlation tau (we only need to specify those pairs needed in the vine copula)
# 
# ## want rshiny to plot vinecopula structure based on the vine copula structure, bicopula family, and parameters
# # e.g., if it's D
# 
# 
# #### shiny app tab: Composite
# # to derive the distribution of the ordinal composite outcome consisting of binary components, with assumed prevalence and correlations.
# 
# 
# ###### Design panel
# n_components = 6 # number of binary components
#  
# ## number of categories in the ordinal composite would be determined based on the user input
# ## composition type 
# # 1. "order" Rank by severity (number of categories of ordinal composite = number of binary components)
# # 2. "sum" Summation (number of categories of ordinal composite = number of binary components+1)
# # 3. "covid" COVID-19 Ordinal Scale Example (number of categories of ordinal composite: 6 if not add death, 7 if add death...)
# # Later add: 4. "custom" Custom Composite (number of categories of ordinal composite: user define)
# 
# # •	Summation → fully connected only
# # •	Rank by severity → fully connected only
# # •	Example: COVID-19 Ordinal Scale → fixed complex structure
# # •	Custom composite → user-defined dependence groups and custom logic
# 
# composition_type = "order" 
# 
# if (composition_type =="order"){
#   # need to specify severity_score
#   severity_score = 1:6
# }
# 
# 
# beta_X = rep(log(.8),n_components) # ask user if they want use all same beta? if not, then keep specifying.
# p_ct = c(.2,.3,.2,.3,.2,.3) # prevalence of control group for each binary component. does not need to sum to 1
# 
# ntree = n_components-1
# #beta_X;p_ct;tau;true_copula
# 
# p_trt = get_p_trt(beta_X=beta_X,p_ct=p_ct)
# vinecop_type = "D"  # if n_components==3, then users don't need to choose.
# # => Rshiny app shows vine copula structure. C-vine (at each tree there is a node that acts as a hub and is connected to all the other nodes) 
# ## or D-vine (at each tree there is a path stringing up all the nodes,)
# add_decay_rate = T ### user can decide if they want a decay rate if the pair is on later tree.
# if (add_decay_rate){
#   decay_rate = 2/3 # if add, they user can define, default = 2/3
# }else{
#   decay_rate = 1
# }
# # decay_rate in (0,1], decides the conditional dependence at deeper trees... how to explain this to user doesn't know much at vine copula
# last_idx = ntree - 1
# ## Correlation Assumption button: high, moderate, low, or very low, later can add customize, not for now.
# tau_values = c(High = .6,
#                Moderate = .4,
#                Low = .2,
#                `Very low` = .05,
#                Indep = 0)
# 
# ### add tau sign 
# tau_list = list(High = round(tau_values[["High"]]*decay_rate^(0:last_idx),2),  # 0.60 0.40 0.27 0.18 0.12
#                 Moderate = round(tau_values[["Moderate"]]*decay_rate^(0:last_idx),2),
#                 Low = round(tau_values[["Low"]]*decay_rate^(0:last_idx),2), # 0.20 0.13 0.09 0.06 0.04
#                 `Very low` = round(tau_values[["Very low"]]*decay_rate^(0:last_idx),2)) 
# corr_assump_level = "Low"
# ##### we are more interested in copula on two binary components. (not continuous...)
# ###### => simulate binary components given the vine structure and correlation assumptions and visualize their distributions?
# copula_family = "Gaussian"
# # copula_family_map <- c(
# #   Indep    = 0,
# #   Gaussian = 1,
# #   Clayton  = 3,
# #   Gumbel   = 4,
# #   Frank    = 5,
# #   Joe      = 6
# # )
# ## 1. C1 with all Gaussian,
# ## 2. D1 with all Gaussian,
# copula_idx = copula_family_map[[copula_family]]
# 
# param  = BiCopTau2Par(family = rep(copula_idx, ntree), tau = tau_list[[corr_assump_level]])
# # 0.30901699 0.20278730 0.14090123 0.09410831 0.06279052
# 
# 
# 
# 
# SetupInfo_trt = VCO_Initialize(nComp = n_components,
#                                prev = p_trt,
#                                vine_type = vinecop_type,
#                                method_list=rep(copula_family,ntree),
#                                rho_list = NULL,
#                                param_list = param) # not use the exact parameters simulating data.
#                              
# SetupInfo_ct = VCO_Initialize(nComp = n_components,
#                               prev = p_ct,
#                               vine_type = vinecop_type,
#                               method_list=rep(copula_family,5),
#                               rho_list = NULL,
#                               param_list = param)
# 
# jointpmf_trt = gen_jointpmf(SetupInfo_trt)$jointpmf
# jointpmf_ct = gen_jointpmf(SetupInfo_ct)$jointpmf
# jointpmf = (jointpmf_trt+jointpmf_ct)/2 # this only for allocation ratio 1?
# 
# # e.g., ord_pmf = Ordinal_PMF(jointpmf,method='sum')
# 
# # component 6 least severe, 1 most severe
# ### user specify severity score among components: 
# # severity_score = 1:6
# ord_pmf = Ordinal_PMF(jointpmf,method=composition_type, severity_score_increasing=severity_score)
# # component 1 least severe, 6 most severe
# 
# #### want to plot the ordinal pmf as histogram, compare with naive e.g., 1/(number of categories of the ordinal composite)
# 
# 
# 
# ###### COVID-19 Ordinal Scale Example
# composition_type = "covid"
# ### preload component's name and vine copula structure
# ## User can modify the prevalence of these components
# # Y1 - Hospitalization; 
# # Y2 - Mechanical ventilation or ECMO; 
# # Y3 - Supplemental oxygen; 
# # Y4 - Symptom; 
# # Y5 - Limitation in activity.
# ### default values:
# preval_all = c(0.0140, 0.0028, 0.0084, 0.2885, 0.1513)
# names(preval_all) = c("hospital","vent","oxygen","symptoms_any","usual_activity")
# 
# 
# ### subgroup 1: hospital---vent
# # assume correlation
# ntree_grp1 = 1
# corr_assump_level = "Low"
# copula_family = "Gaussian"
# copula_idx = copula_family_map[[copula_family]]
# 
# param  = BiCopTau2Par(family = copula_idx, tau =tau_values[[corr_assump_level]])
# # 0.30901699 
# 
# p1 = as.numeric(preval_all['hospital'])
# p2 = as.numeric(preval_all['vent'])   
# 
# pmf1 = matrix(BiCop_PMF(p1,p2,parC= param, #find_parC(p1,p2,rho=corr_mat[1,5]),
#                         specific_method=copula_family,print=T
# ),nrow=1)
# colnames(pmf1) = sort(paste0('X',expand.grid(rep(list(0:1),2)) %>% 
#                                apply(1,function(x) paste0(x,collapse=''))))
# pmf1
# 
# ord_pmf_cat1 = pmf1[,'X11']
# 
# ### subgroup 2: hospital---oxygen
# ntree_grp2 = 1
# corr_assump_level = "High"
# copula_family = "Gaussian"
# copula_idx = copula_family_map[[copula_family]]
# 
# param  = BiCopTau2Par(family = copula_idx, tau = tau_values[[corr_assump_level]])
# # 0.30901699 0.20278730 0.14090123 0.09410831 0.06279052
# 
# p1 = as.numeric(preval_all['hospital'])
# p2 = as.numeric(preval_all['oxygen'])   
# 
# pmf2 = matrix(BiCop_PMF(p1,p2,parC=param,#find_parC(p1,p2,rho=corr_mat[1,4]),
#                         specific_method=copula_family),nrow=1)
# colnames(pmf2) = sort(paste0('X',expand.grid(rep(list(0:1),2)) %>% 
#                                apply(1,function(x) paste0(x,collapse=''))))
# pmf2
# 
# ord_pmf_cat2 = pmf2[,'X11']
# ord_pmf_cat3 = pmf2[,'X10']
# 
# 
# ### subgroup 3: hospital---symptoms---activity
# # low corr between hospital---symptoms
# # high corr between symptoms and activity
# # indep between hospital and activity, conditioning on symptoms
# ntree_grp3 = 2 # = number of binary components in this group - 1
# corr_assump_level_tree1 = c("Low","High")
# corr_assump_level_tree2 = "Indep"
# 
# copula_family = "Gaussian"
# copula_idx = copula_family_map[[copula_family]]
# 
# param_tree1 = BiCopTau2Par(family = rep(copula_idx,2), 
#                      tau = as.numeric(tau_values[corr_assump_level_tree1]))
# # 0.309017 0.809017
# param_tree2 = BiCopTau2Par(family = rep(copula_idx,1), 
#                            tau = tau_values[[corr_assump_level_tree2]])
#   
# prev_grp3 = as.numeric(preval_all[c('hospital',
#                                     'symptoms_any',
#                                     'usual_activity')])
# 
# 
# SetupInfo = VCO_Initialize(nComp = 3,
#                            prev = prev_grp3,
#                            vine_type = 'D', # dosen't matter
#                            method_list=rep(copula_family,ntree_grp3), 
#                            rho_list = NULL,
#                            param_list = list(param_tree1,
#                                              param_tree2))
# 
# output = gen_jointpmf(SetupInfo)
# list_joint_pmf = output$list_joint_pmf
# jointpmf = output$jointpmf
# jointpmf
# 
# ord_pmf_cat4 = jointpmf['X011',]
# 
# ord_pmf_cat5 = jointpmf['X010',]
# 
# ord_pmf_cat6 = jointpmf['X000',]
# 
# 
# ord_pmf_all = as.numeric(c(ord_pmf_cat1,ord_pmf_cat2,ord_pmf_cat3,
#                 ord_pmf_cat4,ord_pmf_cat5,ord_pmf_cat6))
# 
# 
# ###### Later: Custom composite (Full Flexibility)
# 
# 
# 
# 
# 
# # BiCopTau2Par, family can be:
# # 0 = independence copula
# # 1 = Gaussian copula
# # 2 = Student t copula (Here only the first parameter can be computed)
# # 3 = Clayton copula
# # 4 = Gumbel copula
# # 5 = Frank copula
# # 6 = Joe copula
# # 13 = rotated Clayton copula (180 degrees; survival Clayton'') \cr `14` = rotated Gumbel copula (180 degrees; survival Gumbel'')
# # 16 = rotated Joe copula (180 degrees; ``survival Joe'')
# # 23 = rotated Clayton copula (90 degrees)
# # `24` = rotated Gumbel copula (90 degrees)
# # `26` = rotated Joe copula (90 degrees)
# # `33` = rotated Clayton copula (270 degrees)
# # `34` = rotated Gumbel copula (270 degrees)
# # `36` = rotated Joe copula (270 degrees)






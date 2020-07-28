
###################################### FUNCTIONS  ################################
library(tidyverse)

### n x n square Table
squareTable <- function(x,y) {
  x <- factor(x)
  y <- factor(y)
  
  commonLevels <- sort(unique(c(levels(x), levels(y))))
  
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  
  table(x,y)
}

### set variable from varlist into factor format
set_factor<- function(df,varlist){
  for (v in varlist){
    df[,v]=factor(df[,v], ordered = FALSE)
  }
  return(df)
}

### get dataframe with variables from a keeplist
get_df_w_var <- function(keeplist, df){
  keep=df[keeplist]
  return(keep)
}

### get dataframe without variables from a droplist
get_df_wo_var <- function(droplist, df){
  df_remain=df[,!names(df) %in% droplist]
  return(df_remain)
}


### reset the imputed variable back to NA
reset_miss <- function(var,df.true, df.misind ){
  rec=as.data.frame(cbind( df.true$KEY, df.misind[[var]]) )
  rec[[var]]=df.true[[var]]
  rec[[var]][rec[["V2"]]==1  ]=NA
  names(rec)=c("KEY",  paste(var,"_misind",sep = ''), paste(var))
  return(rec)
}


imp_rf_chain=function(tar_var, df, seedn, tree_setting = tree_setting, hyp_grid = hyp_grid, dim0=dim0){
  ###
  ### impute tar_var
  ###
  library(h2o)
  h2o.init(nthreads = -1)
  
  print(paste('impute var=', tar_var))
  print(paste('use seed num: ',seedn))
  mis_tar_ind<-df[,c("KEY",tar_var)]
  mis_tar_ind$ind<-missind[[tar_var]]
  
  train<-df[!is.na(df[[tar_var]] ),]
  test<-df[is.na(df[[tar_var]] ),]
  
  pos_KEY_in_train=which(names(train)=='KEY')
  training<-train[,-c(pos_KEY_in_train)] # create training without KEY 
  
  testKEY=test[,c('KEY', tar_var)]# KEY for test data
  
  h2o_train=as.h2o(training)
  
  #h2o_split=h2o.splitFrame(h2o_train, ratios = .8, seed = seedn)
  #h2o_df_train=h2o_split[[1]]
  #h2o_df_valid=h2o_split[[2]]
  
  ### get default tree number
  if (is.null(tree_setting)==FALSE){
    ntree_setting=tree_setting
  } else if (is.null(tree_setting)) {
    dim=dim0
    ntree_setting=10*dim
  }
  
  # hyper-parameters to be considered
  mtries_opt=floor(c(sqrt(dim), dim/3, dim*2/3, dim/4))
  min_rows_opt = c(1, 3, 5, 10)
  max_depth_opt = c(10, 20, 30)
  sample_rate_opt = c(.55, .632, .80)
  
  if (is.null(hyp_grid)){
    hyper_params <- list(mtries = mtries_opt,
                         min_rows = min_rows_opt,
                         max_depth = max_depth_opt,
                         sample_rate = sample_rate_opt
                         
    ) # end hyper_params list
    
    search_criteria <- list(strategy = "RandomDiscrete", seed=seedn,
                            stopping_tolerance = 0.001,   # stop if improvement is < 0.1%
                            stopping_rounds = 10,         # over the last 10 models
                            max_runtime_secs = 60*2)
    
    ### need to complete select by different conditions
    select_by=c('accuracy','r2','logloss','err','rmse')
    
    gl=list()
    
    gl[[1]]=hyper_params
    gl[[2]]=search_criteria
    gl[[3]]=select_by
    
    dl_grid <- h2o.grid("randomForest", y = tar_var,
                        grid_id = "dl_grid",
                        training_frame = h2o_train,
                        #validation_frame = h2o_df_valid,
                        nfolds=5,
                        seed = seedn,
                        ntrees = ntree_setting,
                        #hidden = c(10,10),
                        balance_classes= TRUE,
                        hyper_params = gl[[1]],
                        search_criteria = gl[[2]])
    
    
  } else if (is.null(hyp_grid)==FALSE) {
    
    hyp_grid=hyp_grid
  }
  
  
  typeof(df[[tar_var]])
  is.factor(df[['ASOURCE']])
  is.numeric(df[[tar_var]])
  is.factor(df[['FEMALE']])
  
  if (is.factor(df[[tar_var]])){
    
    #### consider 5 metrics to select best model
    lmod1=list()
    lmod2=list()
    lsort1=c('accuracy','r2') # larger the better
    lsort2=c('logloss','err','rmse') # smaller the better
    lsort=gl[[3]]
    
    for (ls in 1:length(lsort)){
      
      lsd=lsort[ls]
      if (lsd %in% lsort1) {
        d_opt=TRUE
      } else if (lsd %in% lsort2){
        d_opt=FALSE
      }
      
      dl_gridperf <- h2o.getGrid(grid_id = "dl_grid"
                                 ,sort_by = lsd
                                 ,decreasing = d_opt) # smaller the better
      lmod1[[ls]]=dl_gridperf@model_ids[1]
      lmod2[[ls]]=dl_gridperf@model_ids[2]
      
    }# end loop ls
    
    
    ### scoring the model from different metrics with weight 1 for the 1st, 0.5 for the 2nd 
    cmod1=as.data.frame(table(unlist(lmod1)))
    cmod2=as.data.frame(table(unlist(lmod2)))
    cmod1$sg=cmod1$Freq
    cmod2$sg=cmod2$Freq*0.5
    
    cmod=merge(cmod1,cmod2,by = 'Var1',all = TRUE)
    cmod$sg.x= cmod$sg.x %>% replace_na( 0)
    cmod$sg.y= cmod$sg.y %>% replace_na( 0)
    cmod$sg=(cmod$sg.x+cmod$sg.y)
    
    cmod=cmod[order(-cmod$sg.x,-cmod$sg.y),]
    cmod
    cmod$Var1[1]
    
    #cmod[which(cmod$sg==max(cmod$sg,na.rm = T)),][1]
    #best_dl_model_id=as.character(cmod[which(cmod$sg==max(cmod$sg,na.rm = T)),][['Var1']][1])
    
    best_dl_model_id=as.character(cmod$Var1[1])
    
    ### pick the model with highest sg score
    best_dl <- h2o.getModel(best_dl_model_id)
    
    
    best_dl_model_id
    ### print the info for best model 
    catb=dl_grid@summary_table[which(dl_grid@summary_table$model_ids==best_dl_model_id),]
    
    
    cat('best model for ',tar_var,' : ',unlist(catb),'\n')
    
  } # end if.factor 
  else if (is.numeric(df[[tar_var]])) {
    dl_gridperf <- h2o.getGrid(grid_id = "dl_grid"
                               ,sort_by = 'rmse'
                               ,decreasing = FALSE) # smaller the better
  
    best_dl_model_id = dl_gridperf@model_ids[[1]]  
    best_dl <- h2o.getModel(best_dl_model_id)
    
    
    best_dl_model_id
    ### print the info for best model 
    catb=dl_grid@summary_table[which(dl_grid@summary_table$model_ids==best_dl_model_id),]
    
    
    cat('best model for ',tar_var,' : ',unlist(catb),'\n')
  }
  
  
  
  prpos<-which(names(test)==tar_var)  # find tar_var position in test
  testing<-test[-prpos]
  y_pred = h2o.predict(best_dl, newdata = as.h2o(testing))
  
  y_pred3<-as.data.frame(y_pred$predict) # get predicted tar_var values
  
  test_cb=cbind(testKEY$KEY,y_pred3)
  names(test_cb) = c('KEY', tar_var)
  
  train_cb=train[,c('KEY',tar_var)]
  
  imp_a=rbind(train_cb,test_cb) # merge no missing and predicted missing data
  
  pos_tar_in_df=which(names(df)==tar_var)
  imp_b=merge(df[-c(pos_tar_in_df)], imp_a, by ='KEY', all.x=TRUE)
  
  pos_tar_in_miss_ind=which(names(mis_tar_ind)==tar_var)
  imp_c=merge(imp_b, mis_tar_ind[-c(pos_tar_in_miss_ind)], by='KEY', all.x = TRUE)
  imp_c$ind=NULL  # an option in future to add switch for keep missing indicator
  h2o.removeAll()
  return(imp_c)
  
  
  #  RF.model = h2o.randomForest(
  #    y = tar_var,
  #    training_frame = h2o_df_train,
  #    ntrees = tree_setting,
  #   validation_frame = h2o_df_valid,
  #balance_classes=TRUE,
  #nfolds=5,
  #    max_depth = 30,               ## Increase depth, from 20
  #    stopping_rounds = 2,          ##
  #    stopping_tolerance = 1e-2,    ##
  #    score_each_iteration = T,
  #    seed = seedn
  # missing_values_handling="MeanImputation",
  # sparse=TRUE,
  # train_samples_per_iteration = -2
  # )
  
  #summary(RF.model)
  
}



chain_imp_rf = function(data_in, list_imp_var, list_no_imp_var=NULL, n_iter, seedn= seedn, 
                        tree_setting = tree_setting, hyp_grid = hyp_grid, dim0=dim0){
  ### 1st iter
  # 5 variables: pred_ztotchg    FEMALE  bzipinc ASOURCE  RACE   
  dim0=dim0
  
  print("iteration=  1")
  
  kl=c("KEY",list_imp_var)
  dl=list()
  dl[[1]]=data_in
  for ( v in list_imp_var){
    print(paste('impute var=',v))
    posv=which(list_imp_var==v)
    prep_temp1_p1=get_df_wo_var(droplist = list_imp_var[-posv], df = data_in)
    prep_temp1_p2=get_df_w_var(keeplist  = kl[1:posv], df = dl[[posv]])
    prep_temp1=merge(prep_temp1_p1, prep_temp1_p2, by="KEY")
    
    impc=imp_rf_chain(tar_var = v, df = prep_temp1 , seedn = seedn, tree_setting = tree_setting,hyp_grid = hyp_grid, dim0=dim0)
    dl[[posv]]=0
    dl[[posv+1]]=impc
    
  }
  
  
  cyc1=dl[[length(list_imp_var)+1]]
  
  
  iter_cyc=list()
  iter_cyc[[1]]=cyc1
  
  ### iter >=2
  
  for (it in 2:n_iter){
    seedn= seedn+it
    print(paste('iteration: ',it))
    print(paste('seed num: ',seedn))
    
    chain_cyc=list()
    chain_cyc[[1]]=iter_cyc[[it-1]]
    
    #indv_chain_cyc=list()
    #indv_chain_cyc[[1]]=chain_cyc[[1]]
    for (vo in list_imp_var){
      posvo=which(list_imp_var==vo)
      recyc = reset_miss(var = vo, df.true = chain_cyc[[posvo]], df.misind = missind)
      recyc_temp1_p1=get_df_wo_var(droplist = c(paste(vo,'_misind',sep = '')), df = recyc)
      recyc_temp1_p2=get_df_wo_var(droplist = c(vo), df = chain_cyc[[posvo]])
      recyc_temp1=merge(recyc_temp1_p1, recyc_temp1_p2, by="KEY")
      
      recyc_temp2=imp_rf_chain(tar_var = vo, df = recyc_temp1, seedn = seedn, tree_setting = tree_setting,hyp_grid = hyp_grid, dim0=dim0)
      
      chain_cyc[[posvo+1]]=recyc_temp2
      chain_cyc[[posvo]]=0
    }
    iter_cyc[[it]]=chain_cyc[[length(list_imp_var)+1]]
    iter_cyc[[it-1]]=0
  }
  end_cyc=iter_cyc[[n_iter]]
  return(end_cyc)
}  

### RF gen m independent df 
mice_imp_rf=function(df_in, m, var_to_imp_list, n_iter, tree_setting = tree_setting, hyp_grid = hyp_grid){
  imp_mice=list()
  for(mt in 1:m){
    print(paste('working on dataset m=', mt))
    dim0=dim(df_in)[2]-1
    seedn=(9527+u)*mt
    imp_mice[[mt]]=chain_imp_rf(data_in = df_in, list_imp_var = var_to_imp_list, n_iter = n_iter, seedn= seedn, tree_setting=tree_setting,hyp_grid = hyp_grid, dim0=dim0)
  }
  return(imp_mice)
}

### RF gen m independent df one at a time
one_mice_imp_rf=function(df_in, mt, var_to_imp_list, n_iter, tree_setting = tree_setting, hyp_grid = hyp_grid){
  print(paste('working on dataset m=', mt))
  dim0=dim(df_in)[2]-1
  seedn=(9527+u)*mt
  imp_a_mice=chain_imp_rf(data_in = df_in, list_imp_var = var_to_imp_list, n_iter = n_iter, seedn= seedn, tree_setting=tree_setting,hyp_grid = hyp_grid, dim0=dim0)
  return(imp_a_mice)
}

### get distribution for var of interests
get_var_dist = function(tar_var, df){
  if (is.numeric(df[[tar_var]])) {
    dist=list(1,'AUTO',FALSE)
  } else if (is.factor(df[[tar_var]])) {
    if (length(levels(data_simm[[tar_var]]))==2) {
      dist=list(2,'bernoulli',TRUE)
    } else if (length(levels(data_simm[[tar_var]]))>2 ) {
      dist=list(length(levels(data_simm[[tar_var]])), 'multinomial', TRUE)
    }
    
  }
  return(dist)
}




imp_nn_chain=function(tar_var, df, seedn, hyp_grid=hyp_grid, dim0=dim0){
  ###
  ### impute tar_var
  ###
  library(h2o)
  h2o.init(nthreads = -1)
  
  var_dist= get_var_dist(tar_var = tar_var, df = df)
  
  
  print(paste('impute var=', tar_var))
  print(paste('use seed num: ',seedn))
  mis_tar_ind<-df[,c("KEY",tar_var)]
  mis_tar_ind$ind<-missind[[tar_var]]
  
  train<-df[!is.na(df[[tar_var]] ),]
  test<-df[is.na(df[[tar_var]] ),]
  
  pos_KEY_in_train=which(names(train)=='KEY')
  training<-train[,-c(pos_KEY_in_train)] # create training without KEY 
  
  testKEY=test[,c('KEY', tar_var)]# KEY for test data
  
  h2o_train=as.h2o(training)
  
  #h2o_split=h2o.splitFrame(h2o_train, ratios = .8, seed = seedn)
  #h2o_df_train=h2o_split[[1]]
  #h2o_df_valid=h2o_split[[2]]
  
  
  ### get default node number
  dim= dim0
  hnodes1= ceiling((dim + var_dist[[1]] )/2)
  hnodes2= ceiling((dim + var_dist[[1]] )/3*2)
  #print(paste("dim0 s3", dim0))
  
  #  if (is.null(hidden)){
  #    hnodes= ceiling((dim + var_dist[[1]] )/2)
  #    hlayer= c(hnodes,hnodes,hnodes)
  #  } else {
  #    hlayer=hidden
  #  }
  
  ###################### Not auto adapt yet
  
  #  if (is.null(hidden)){
  #    if (tar_var %in% c('ASOURCE','FEMALE'))
  #    {
  #      cat('NN for ',tar_var)
  #      hlayer=c(30,90,30)
  #    } else if (tar_var %in% c('RACE','ZIPINC_QRTL')) {
  #      cat('NN for ',tar_var)
  #      hlayer=c(40,80,40)
  #   } else if (tar_var %in% c('binTOTCHG')) {
  #     cat('NN for ',tar_var)
  #     hlayer=c(40,80,20)
  #  } else if (tar_var %in% c('zTOTCHG')) {
  #    cat('NN for ',tar_var)
  #    hnodes= ceiling((dim + var_dist[[1]] )/2)
  #    hlayer= c(hnodes,hnodes,hnodes)
  #   }  
  # }
  
  if (is.null(hyp_grid)){
    # activation_opt <- c("Rectifier", "RectifierWithDropout", "Maxout", "MaxoutWithDropout", "Tanh", "TanhWithDropout")
    activation_opt <- c("Rectifier", "RectifierWithDropout",  "Tanh", "TanhWithDropout")
    l1_opt <- c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1)
    l2_opt <- c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1)
    hidden_opt <-list(c(hnodes1,hnodes1,hnodes1)
                      #,c(hnodes1,hnodes1)
                      ,c(hnodes1,5*hnodes1,2*hnodes1)
                      ,c(hnodes2,hnodes2,hnodes2)
                      ,c(round(hnodes2*5,-1), 3*round(hnodes2*5,-1), round(hnodes2*3,-1))
                      # ,c(30,90,30), c(40,80,40), c(40,80,20)
    )
    
    ### sample_factors
    #if (tar_var=='ASOURCE') {
    #  sample_factor=list(c(0.01, 0.6, 1),c(0.01,0.8,1))
    #} else if (tar_var=='FEMALE') {
    #  sample_factor=list(c( 1, 0.1),c(0.1,1))
    #}  else if (tar_var=='RACE') {
    #  sample_factor=list(c(.02,1,.7),c(.02,1,.75))
    #}  else if (tar_var=='ZIPINC_QRTL') {
    #  sample_factor=list(c(1.0, 1.0, 0.2, 0.15),c(1.0, 1.0, 0.2, 0.1) )
    #}  else if (tar_var=='binTOTCHG') {
    #  sample_factor=list(c(1,1))
    #}  
    
    sfa=1/prop.table(table(df[[tar_var]]))
    sfb=round(sfa/max(sfa),1)
    sample_factors= list(c(as.matrix(sfb)))
    
    
    
    
    
    
    hyper_params <- list(activation = activation_opt,
                         l1 = l1_opt,
                         l2 = l2_opt,
                         class_sampling_factors =sample_factors,
                         hidden = hidden_opt
                         
    )
    
    search_criteria <- list(strategy = "RandomDiscrete", seed=seedn,
                            max_runtime_secs = 60)
    
    
    
    
    
    
    ### need to complete select by different conditions
    select_by=c('accuracy','r2','logloss','err','rmse')
    
    gl=list()
    
    gl[[1]]=hyper_params
    gl[[2]]=search_criteria
    gl[[3]]=select_by
    
    dl_grid <- h2o.grid("deeplearning", y = tar_var,
                        grid_id = "dl_grid",
                        training_frame = h2o_train,
                        #validation_frame = h2o_df_valid,
                        nfolds=5,
                        seed = seedn,
                        #hidden = c(10,10),
                        balance_classes= TRUE,
                        hyper_params = gl[[1]],
                        search_criteria = gl[[2]])
    
    
  } else if (is.null(hyp_grid)==FALSE) {
    
    hyp_grid=hyp_grid
  }
  
  
  #### consider 5 metrics to select best model
  lmod1=list()
  lmod2=list()
  lsort1=c('accuracy','r2') # larger the better
  lsort2=c('logloss','err','rmse') # smaller the better
  lsort=gl[[3]]
  
  for (ls in 1:length(lsort)){
    
    lsd=lsort[ls]
    if (lsd %in% lsort1) {
      d_opt=TRUE
    } else if (lsd %in% lsort2){
      d_opt=FALSE
    }
    
    dl_gridperf <- h2o.getGrid(grid_id = "dl_grid"
                               ,sort_by = lsd
                               ,decreasing = d_opt) # smaller the better
    lmod1[[ls]]=dl_gridperf@model_ids[1]
    lmod2[[ls]]=dl_gridperf@model_ids[2]
    
  }# end loop ls
  
  
  ### scoring the model from different metrics with weight 1 for the 1st, 0.5 for the 2nd 
  cmod1=as.data.frame(table(unlist(lmod1)))
  cmod2=as.data.frame(table(unlist(lmod2)))
  cmod1$sg=cmod1$Freq
  cmod2$sg=cmod2$Freq*0.5
  
  cmod=merge(cmod1,cmod2,by = 'Var1',all = TRUE)
  cmod$sg.x= cmod$sg.x %>% replace_na( 0)
  cmod$sg.y= cmod$sg.y %>% replace_na( 0)
  cmod$sg=(cmod$sg.x+cmod$sg.y)
  
  cmod=cmod[order(-cmod$sg.x,-cmod$sg.y),]
  cmod
  cmod$Var1[1]
  
  #cmod[which(cmod$sg==max(cmod$sg,na.rm = T)),][1]
  #best_dl_model_id=as.character(cmod[which(cmod$sg==max(cmod$sg,na.rm = T)),][['Var1']][1])
  
  best_dl_model_id=as.character(cmod$Var1[1])
  
  ### pick the model with highest sg score
  best_dl <- h2o.getModel(best_dl_model_id)
  
  
  best_dl_model_id
  ### print the info for best model 
  catb=dl_grid@summary_table[which(dl_grid@summary_table$model_ids==best_dl_model_id),]
  
  cat('best model for ',tar_var,' : ',unlist(catb),'\n')
  
  
  
  
  prpos<-which(names(test)==tar_var)  # find tar_var position in test
  testing<-test[-prpos]
  y_pred = h2o.predict(best_dl, newdata = as.h2o(testing))
  
  y_pred3<-as.data.frame(y_pred$predict) # get predicted tar_var values
  
  test_cb=cbind(testKEY$KEY,y_pred3)
  names(test_cb) = c('KEY', tar_var)
  
  train_cb=train[,c('KEY',tar_var)]
  
  imp_a=rbind(train_cb,test_cb) # merge no missing and predicted missing data
  
  pos_tar_in_df=which(names(df)==tar_var)
  imp_b=merge(df[-c(pos_tar_in_df)], imp_a, by ='KEY', all.x=TRUE)
  
  pos_tar_in_miss_ind=which(names(mis_tar_ind)==tar_var)
  imp_c=merge(imp_b, mis_tar_ind[-c(pos_tar_in_miss_ind)], by='KEY', all.x = TRUE)
  imp_c$ind=NULL  # an option in future to add switch for keep missing indicator
  h2o.removeAll()
  return(imp_c)
  
} # end function imp_nn_chain




chain_imp_nn = function(data_in, list_imp_var, list_no_imp_var=NULL, n_iter, seedn= seedn, hyp_grid = hyp_grid, dim0=dim0){
  ### 1st iter
  # 5 variables: pred_ztotchg    FEMALE  bzipinc ASOURCE  RACE   
  
  dim0=dim0
  #print(paste('dim0 s2: ', dim0))
  
  print("iteration=  1")
  
  kl=c("KEY",list_imp_var)
  dl=list()
  dl[[1]]=data_in
  for ( v in list_imp_var){
    print(paste('impute var=',v))
    posv=which(list_imp_var==v)
    prep_temp1_p1=get_df_wo_var(droplist = list_imp_var[-posv], df = data_in)
    prep_temp1_p2=get_df_w_var(keeplist  = kl[1:posv], df = dl[[posv]])
    prep_temp1=merge(prep_temp1_p1, prep_temp1_p2, by="KEY")
    
    impc=imp_nn_chain(tar_var = v, df = prep_temp1 , seedn = seedn, hyp_grid = hyp_grid, dim0=dim0)
    dl[[posv]]=0
    dl[[posv+1]]=impc
    
  }
  
  
  cyc1=dl[[length(list_imp_var)+1]]
  
  
  iter_cyc=list()
  iter_cyc[[1]]=cyc1
  
  ### iter >=2
  
  for (it in 2:n_iter){
    seedn= seedn+it
    print(paste('iteration: ',it))
    print(paste('seed num: ',seedn))
    
    chain_cyc=list()
    chain_cyc[[1]]=iter_cyc[[it-1]]
    
    #indv_chain_cyc=list()
    #indv_chain_cyc[[1]]=chain_cyc[[1]]
    for (vo in list_imp_var){
      posvo=which(list_imp_var==vo)
      recyc = reset_miss(var = vo, df.true = chain_cyc[[posvo]], df.misind = missind)
      recyc_temp1_p1=get_df_wo_var(droplist = c(paste(vo,'_misind',sep = '')), df = recyc)
      recyc_temp1_p2=get_df_wo_var(droplist = c(vo), df = chain_cyc[[posvo]])
      recyc_temp1=merge(recyc_temp1_p1, recyc_temp1_p2, by="KEY")
      
      recyc_temp2=imp_nn_chain(tar_var = vo, df = recyc_temp1, seedn = seedn, hyp_grid = hyp_grid, dim0=dim0)
      
      chain_cyc[[posvo+1]]=recyc_temp2
      chain_cyc[[posvo]]=0
    }
    iter_cyc[[it]]=chain_cyc[[length(list_imp_var)+1]]
    iter_cyc[[it-1]]=0
  }
  end_cyc=iter_cyc[[n_iter]]
  return(end_cyc)
}  






### NN gen m independent df 
mice_imp_nn=function(df_in, m, var_to_imp_list, n_iter, hyp_grid=hyp_grid){
  imp_mice=list()
  for(mt in 1:m){
    print(paste('working on dataset m=', mt))
    dim0=dim(df_in)[2]-1
    seedn=(9527+u)*mt
    imp_mice[[mt]]=chain_imp_nn(data_in = df_in, list_imp_var = var_to_imp_list, n_iter = n_iter, seedn= seedn, hyp_grid=hyp_grid)
  }
  return(imp_mice)
}


### NN gen m independent df one at a time
one_mice_imp_nn=function(df_in, mt, var_to_imp_list, n_iter, hyp_grid=hyp_grid){
  print(paste('working on dataset m=', mt))
  dim0=dim(df_in)[2]-1
  #print(paste('dim0 s1: ', dim0))
  seedn=(9527+u)*mt
  imp_a_mice=chain_imp_nn(data_in = df_in, list_imp_var = var_to_imp_list, n_iter = n_iter, seedn= seedn, hyp_grid=hyp_grid, dim0 = dim0)
  return(imp_a_mice)
}



### return df only contains var of interests
get_df_eval <- function(var,df.true, df.misind, df.imp ){
  df_cm1=as.data.frame(cbind( df.true$KEY, df.true[[var]], df.misind[[var]]) )
  names(df_cm1)=c("KEY", var, 'misind' )
  df_cm2=df.imp[c('KEY', var)]
  names(df_cm2)=c('KEY','pred')
  df_cm= merge(df_cm1, df_cm2, by='KEY')
  df_var=df_cm[which(df_cm$misind==TRUE),]
  return(df_var)
}

### return list of CM and acc for categorical vars after imp
get_acc_var <- function(var, df_eval){
  c0=squareTable(df_eval[[var]], df_eval$pred)
  acc=sum(diag(c0))/sum(c0)
  res_lst=list(c0, acc)
  names(res_lst)=c('CM', 'Acc')
  return(res_lst)
}

### return list of CM and acc for continuous vars after imp
get_rmsd_var <- function(var, df_eval){
  df_eval$diff2=(df_eval[[var]] - df_eval$pred)^2
  rmsd=sqrt(mean(df_eval$diff2))
  return(rmsd)
}

### return list of acc values for list of categorical vars
acc_list <- function(df, dft, dfm, var_to_acc_list){
  list_acc=list()
  for (va in var_to_acc_list){
    posva=which(var_to_acc_list==va)
    df_v=get_df_eval(va, df.true = dft, df.misind = dfm, df.imp = df)
    c0=squareTable(df_v[[va]], df_v$pred)
    acc0=sum(diag(c0))/sum(c0)
    list_acc[[posva]]=acc0
    names(list_acc[[posva]])=va
  }
  return(list_acc)
}




### return list of acc values for list of continuous vars
rmsd_list <- function(df, dft, dfm, var_to_acc_list){
  list_rmsd=list()
  for (va in var_to_acc_list){
    posva=which(var_to_acc_list==va)
    rmsd0 = get_rmsd_var(va, df_eval= df)
    list_rmsd[[posva]]=rmsd0
    names(list_rmsd[[posva]])=va
  }
  return(list_rmsd)
}






















###### impute using xgb
imp_xgb_chain=function(tar_var, df, seedn, tree_setting){
  ###
  ### impute tar_var
  ###
  library(h2o)
  h2o.init(nthreads = -1)
  var_dist= get_var_dist(tar_var = tar_var, df = df)
  
  
  print(paste('impute var=', tar_var))
  print(paste('use seed num: ',seedn))
  mis_tar_ind<-df[,c("KEY",tar_var)]
  mis_tar_ind$ind<-missind[[tar_var]]
  
  train<-df[!is.na(df[[tar_var]] ),]
  test<-df[is.na(df[[tar_var]] ),]
  
  pos_KEY_in_train=which(names(train)=='KEY')
  training<-train[,-c(pos_KEY_in_train)] # create training without KEY 
  
  testKEY=test[,c('KEY', tar_var)]# KEY for test data
  
  h2o_train=as.h2o(training)
  
  #h2o_split=h2o.splitFrame(h2o_train, ratios = .8, seed = seedn)
  h2o_split=h2o.splitFrame(h2o_train, ratios = .7, seed = seedn)
  
  h2o_df_train=h2o_split[[1]]
  h2o_df_valid=h2o_split[[2]]
  
  xgb.model = h2o.xgboost(
    y = tar_var,
    training_frame = h2o_df_train,
    
    validation_frame = h2o_df_valid,
    #balance_classes=TRUE,
    #nfolds=5,
    
    seed = seedn
    
    
    #,stopping_rounds = 3
    ,stopping_rounds = 20
    ,stopping_tolerance = 0.0001
    ,stopping_metric = "AUTO"
    ,distribution = var_dist[[2]]
    ,score_tree_interval = 1
    ,max_depth = 15
    #,learn_rate=0.1
    ,learn_rate=0.001
    ,ntrees=tree_setting
    ,subsample = 0.75
    ,colsample_bytree = 0.75
    #,subsample = 0.7
    #,colsample_bytree = 0.7
    #,tree_method = "hist"
    #,grow_policy = "lossguide"
    ,booster = "gbtree"
    #,gamma = 0.0
    ,gamma = 100
    
  )
  
  #summary(xgb.model)
  
  prpos<-which(names(test)==tar_var)  # find tar_var position in test
  testing<-test[-prpos]
  xgb.y_pred = h2o.predict(xgb.model, newdata = as.h2o(testing))
  
  #summary(xgb.y_pred)
  #npred<-length(test[[tar_var]])  # number of zTOTCHG predicted
  #xgb.y_pred2<-as.vector(xgb.y_pred)[1:npred] # predicted zTOTCHG
  
  
  xgb.y_pred3<-as.data.frame(xgb.y_pred$predict)
  
  test_cb=cbind(testKEY$KEY,xgb.y_pred3)
  names(test_cb) = c('KEY', tar_var)
  
  train_cb=train[,c('KEY',tar_var)]
  
  imp_a=rbind(train_cb,test_cb)
  
  pos_tar_in_df=which(names(df)==tar_var)
  imp_b=merge(df[-c(pos_tar_in_df)], imp_a, by ='KEY', all.x=TRUE)
  
  pos_tar_in_miss_ind=which(names(mis_tar_ind)==tar_var)
  imp_c=merge(imp_b, mis_tar_ind[-c(pos_tar_in_miss_ind)], by='KEY', all.x = TRUE)
  imp_c$ind=NULL
  h2o.removeAll()
  return(imp_c)
}



chain_imp_xgb = function(data_in, list_imp_var, list_no_imp_var=NULL, n_iter, seedn= seedn, tree_setting = tree_setting){
  ### 1st iter
  # 5 variables: pred_ztotchg    FEMALE  bzipinc ASOURCE  RACE   
  
  print("iteration=  1")
  
  kl=c("KEY",list_imp_var)
  dl=list()
  dl[[1]]=data_in
  for ( v in list_imp_var){
    print(paste('impute var=',v))
    posv=which(list_imp_var==v)
    prep_temp1_p1=get_df_wo_var(droplist = list_imp_var[-posv], df = data_in)
    prep_temp1_p2=get_df_w_var(keeplist  = kl[1:posv], df = dl[[posv]])
    prep_temp1=merge(prep_temp1_p1, prep_temp1_p2, by="KEY")
    
    impc=imp_xgb_chain(tar_var = v, df = prep_temp1 , seedn = seedn, tree_setting = tree_setting)
    dl[[posv]]=0
    dl[[posv+1]]=impc
    
  }
  
  
  cyc1=dl[[length(list_imp_var)+1]]
  
  
  iter_cyc=list()
  iter_cyc[[1]]=cyc1
  
  ### iter >=2
  
  for (it in 2:n_iter){
    seedn= seedn+it
    print(paste('iteration: ',it))
    print(paste('seed num: ',seedn))
    
    chain_cyc=list()
    chain_cyc[[1]]=iter_cyc[[it-1]]
    
    #indv_chain_cyc=list()
    #indv_chain_cyc[[1]]=chain_cyc[[1]]
    for (vo in list_imp_var){
      posvo=which(list_imp_var==vo)
      recyc = reset_miss(var = vo, df.true = chain_cyc[[posvo]], df.misind = missind)
      recyc_temp1_p1=get_df_wo_var(droplist = c(paste(vo,'_misind',sep = '')), df = recyc)
      recyc_temp1_p2=get_df_wo_var(droplist = c(vo), df = chain_cyc[[posvo]])
      recyc_temp1=merge(recyc_temp1_p1, recyc_temp1_p2, by="KEY")
      
      recyc_temp2=imp_xgb_chain(tar_var = vo, df = recyc_temp1, seedn = seedn, tree_setting = tree_setting)
      
      chain_cyc[[posvo+1]]=recyc_temp2
      chain_cyc[[posvo]]=0
    }
    iter_cyc[[it]]=chain_cyc[[length(list_imp_var)+1]]
    iter_cyc[[it-1]]=0
  }
  end_cyc=iter_cyc[[n_iter]]
  return(end_cyc)
}  

### xgb gen m independent df 
mice_imp_xgb=function(df_in, m, var_to_imp_list, n_iter, tree_setting){
  imp_mice=list()
  for(mt in 1:m){
    print(paste('working on dataset m=', mt))
    seedn=(9527+u)*mt
    imp_mice[[mt]]=chain_imp_xgb(data_in = df_in, list_imp_var = var_to_imp_list, n_iter = n_iter, seedn= seedn, tree_setting=tree_setting)
  }
  return(imp_mice)
}

### xgb gen m independent df one at a time
one_mice_imp_xgb=function(df_in, mt, var_to_imp_list, n_iter, tree_setting){
  print(paste('working on dataset m=', mt))
  seedn=(9527+u)*mt
  imp_a_mice=chain_imp_xgb(data_in = df_in, list_imp_var = var_to_imp_list, n_iter = n_iter, seedn= seedn, tree_setting=tree_setting)
  return(imp_a_mice)
}




























###### impute using DENOISE AUTOENCODER
imp_dae_chain=function(tar_var, df, seedn, hyp_grid=NULL, dim0=dim0){
  ###
  ### impute tar_var
  ###
  library(h2o)
  h2o.init(nthreads = -1)
  
  var_dist= get_var_dist(tar_var = tar_var, df = df)
  
  
  print(paste('impute var=', tar_var))
  print(paste('use seed num: ',seedn))
  mis_tar_ind<-df[,c("KEY",tar_var)]
  mis_tar_ind$ind<-missind[[tar_var]]
  
  train<-df[!is.na(df[[tar_var]] ),]
  test<-df[is.na(df[[tar_var]] ),]
  
  pos_KEY_in_train=which(names(train)=='KEY')
  training<-train[,-c(pos_KEY_in_train)] # create training without KEY 
  
  testKEY=test[,c('KEY', tar_var)]# KEY for test data
  
  h2o_train=as.h2o(training)
  
  h2o_split=h2o.splitFrame(h2o_train, ratios = .6, seed = seedn)
  h2o_df_train=h2o_split[[1]]
  h2o_df_valid=h2o_split[[2]]
  
  
  ### get default node number
  dim= dim0
  hnodes1= ceiling((dim + var_dist[[1]] )/2)
  hnodes2= ceiling((dim + var_dist[[1]] )/3*2)
  hnodes3= dim + var_dist[[1]]
  #print(paste("dim0 s3", dim0))
  
  #  if (is.null(hidden)){
  #    hnodes= ceiling((dim + var_dist[[1]] )/2)
  #    hlayer= c(hnodes,hnodes,hnodes)
  #  } else {
  #    hlayer=hidden
  #  }
  
  ###################### Not auto adapt yet
  
  #  if (is.null(hidden)){
  #    if (tar_var %in% c('ASOURCE','FEMALE'))
  #    {
  #      cat('NN for ',tar_var)
  #      hlayer=c(30,90,30)
  #    } else if (tar_var %in% c('RACE','ZIPINC_QRTL')) {
  #      cat('NN for ',tar_var)
  #      hlayer=c(40,80,40)
  #   } else if (tar_var %in% c('binTOTCHG')) {
  #     cat('NN for ',tar_var)
  #     hlayer=c(40,80,20)
  #  } else if (tar_var %in% c('zTOTCHG')) {
  #    cat('NN for ',tar_var)
  #    hnodes= ceiling((dim + var_dist[[1]] )/2)
  #    hlayer= c(hnodes,hnodes,hnodes)
  #   }  
  # }
  
  if (is.null(hyp_grid)){
    activation_opt <- c("Rectifier", "RectifierWithDropout", "Tanh", "TanhWithDropout")
    #l1_opt <- c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1)
    #l2_opt <- c(0, 0.00001, 0.0001, 0.001, 0.01, 0.1)
    l1_opt <- c(0, 0.001, 0.01, 0.1)
    l2_opt <- c(0, 0.001, 0.01, 0.1)
    hidden_opt <-list(#c(hnodes1,hnodes1,hnodes1)
      #,c(hnodes1,hnodes1)
      c(hnodes1,5*hnodes1,2*hnodes1)
      #,c(hnodes2,hnodes2,hnodes2)
      ,c(round(hnodes2*5,-1), 3*round(hnodes2*5,-1), round(hnodes2*3,-1))
      ,c(hnodes3,hnodes3+7,hnodes3+14,hnodes3+21,hnodes3+14,hnodes3+7)
      # ,c(30,90,30), c(40,80,40), c(40,80,20)
    )
    
    elastic_mv_rt_opt=c(.9,.95)
    elastic_reglar_opt=c(0,.0001,.001,.01)
    epoch_opt=c(1000,5000,13000)
    
    
    hyper_params <- list(activation = activation_opt,
                         l1 = l1_opt,
                         l2 = l2_opt,
                         #class_sampling_factors =sample_factors,
                         hidden = hidden_opt,
                         elastic_averaging = TRUE,
                         elastic_averaging_moving_rate=elastic_mv_rt_opt,
                         elastic_averaging_regularization = elastic_reglar_opt
                         
    )
    
    search_criteria <- list(strategy = "RandomDiscrete", seed=seedn,
                            max_runtime_secs = 120)
    
    
    
    ### need to complete select by different conditions
    select_by=c('accuracy','r2','logloss','err','rmse')
    
    gl=list()
    
    gl[[1]]=hyper_params
    gl[[2]]=search_criteria
    gl[[3]]=select_by
    
    dl_grid <- h2o.grid("deeplearning", x=names(h2o_df_train),
                        grid_id = "dl_grid",
                        training_frame = h2o_df_train,
                        validation_frame = h2o_df_valid,
                        #nfolds=5,
                        seed = seedn,
                        autoencoder=T,
                        input_dropout_ratio = 0.2,
                        #hidden = c(10,10),
                        #balance_classes= TRUE,
                        #sparse = TRUE,
                        hyper_params = gl[[1]],
                        search_criteria = gl[[2]])
    
    
  } else if (is.null(hyp_grid)==FALSE) {
    
    hyp_grid=hyp_grid
  }
  
  
  dl_gridperf <- h2o.getGrid(grid_id = "dl_grid"
                             ,sort_by = 'rmse'
                             ,decreasing = FALSE)
  
  
  best_dl_model_id=dl_gridperf@model_ids[[1]]
  
  ### pick the model with highest sg score
  best_dl <- h2o.getModel(best_dl_model_id)
  
  
  best_dl_model_id
  ### print the info for best model 
  catb=dl_grid@summary_table[which(dl_grid@summary_table$model_ids==best_dl_model_id),]
  
  cat('best model for ',tar_var,' : ',unlist(catb),'\n')
  
  
  
  
  #prpos<-which(names(test)==tar_var)  # find tar_var position in test
  #testing<-test[-prpos]
  testing<-test
  y_pred = h2o.predict(best_dl, newdata = as.h2o(testing))
  
  
  y_pred2<-as.data.frame(y_pred) # get predicted tar_var values
  
  summary(y_pred2)
  
  rn=names(y_pred2)
  krn=grep(pattern = tar_var, x=rn, value = TRUE)
  
  
  y_pred2=as.data.frame(y_pred2[,krn])
  head(y_pred2)
  
  
  #y_pred3=round(y_pred2)
  #table(y_pred3)
  
  y_pred2$predlab=colnames(y_pred2)[apply(y_pred2,1,which.max)]
  
  #y_pred2$pred=apply(y_pred2,1,function(x) which(x==max(x)))   ## return the column label of max value in a row, which is the predicted level
  y_pred2$pred= as.numeric(gsub("[^\\d]+", "", y_pred2$predlab, perl=TRUE))  ## return the column label of max value in a row, which is the predicted level
  
  head(y_pred2)
  
  y_pred3=as.data.frame(y_pred2[,c('pred')])
  
  test_cb=cbind(testKEY$KEY,y_pred3)
  names(test_cb) = c('KEY', tar_var)
  
  train_cb=train[,c('KEY',tar_var)]
  
  imp_a=rbind(train_cb,test_cb) # merge no missing and predicted missing data
  
  pos_tar_in_df=which(names(df)==tar_var)
  imp_b=merge(df[-c(pos_tar_in_df)], imp_a, by ='KEY', all.x=TRUE)   # merge the imputed tar_var value to df with tar_var missing
  
  pos_tar_in_miss_ind=which(names(mis_tar_ind)==tar_var)
  imp_c=merge(imp_b, mis_tar_ind[-c(pos_tar_in_miss_ind)], by='KEY', all.x = TRUE)
  imp_c$ind=NULL  # an option in future to add switch for keep missing indicator
  
  h2o.removeAll()
  return(imp_c)
  
} # end function imp_dae_chain




chain_imp_dae = function(data_in, list_imp_var, list_no_imp_var=NULL, n_iter, seedn= seedn, hyp_grid = hyp_grid, dim0=dim0){
  ### 1st iter
  # 5 variables: pred_ztotchg    FEMALE  bzipinc ASOURCE  RACE   
  
  dim0=dim0
  #print(paste('dim0 s2: ', dim0))
  
  print("iteration=  1")
  
  kl=c("KEY",list_imp_var)
  dl=list()
  dl[[1]]=data_in
  for ( v in list_imp_var){
    print(paste('impute var=',v))
    posv=which(list_imp_var==v)
    prep_temp1_p1=get_df_wo_var(droplist = list_imp_var[-posv], df = data_in)
    prep_temp1_p2=get_df_w_var(keeplist  = kl[1:posv], df = dl[[posv]])
    prep_temp1=merge(prep_temp1_p1, prep_temp1_p2, by="KEY")
    
    impc=imp_dae_chain(tar_var = v, df = prep_temp1 , seedn = seedn, hyp_grid = hyp_grid, dim0=dim0)
    dl[[posv]]=0
    dl[[posv+1]]=impc
    
  }
  
  
  cyc1=dl[[length(list_imp_var)+1]]
  
  
  iter_cyc=list()
  iter_cyc[[1]]=cyc1
  
  ### iter >=2
  
  for (it in 2:n_iter){
    seedn= seedn+it
    print(paste('iteration: ',it))
    print(paste('seed num: ',seedn))
    
    chain_cyc=list()
    chain_cyc[[1]]=iter_cyc[[it-1]]
    
    #indv_chain_cyc=list()
    #indv_chain_cyc[[1]]=chain_cyc[[1]]
    for (vo in list_imp_var){
      posvo=which(list_imp_var==vo)
      recyc = reset_miss(var = vo, df.true = chain_cyc[[posvo]], df.misind = missind)
      recyc_temp1_p1=get_df_wo_var(droplist = c(paste(vo,'_misind',sep = '')), df = recyc)
      recyc_temp1_p2=get_df_wo_var(droplist = c(vo), df = chain_cyc[[posvo]])
      recyc_temp1=merge(recyc_temp1_p1, recyc_temp1_p2, by="KEY")
      
      recyc_temp2=imp_dae_chain(tar_var = vo, df = recyc_temp1, seedn = seedn, hyp_grid = hyp_grid, dim0=dim0)
      
      chain_cyc[[posvo+1]]=recyc_temp2
      chain_cyc[[posvo]]=0
    }
    iter_cyc[[it]]=chain_cyc[[length(list_imp_var)+1]]
    iter_cyc[[it-1]]=0
  }
  end_cyc=iter_cyc[[n_iter]]
  return(end_cyc)
}  






### DAE gen m independent df 
mice_imp_dae=function(df_in, m, var_to_imp_list, n_iter, hyp_grid=NULL){
  imp_mice=list()
  for(mt in 1:m){
    print(paste('working on dataset m=', mt))
    dim0=dim(df_in)[2]-1
    seedn=(9527+u)*mt
    imp_mice[[mt]]=chain_imp_dae(data_in = df_in, list_imp_var = var_to_imp_list, n_iter = n_iter, seedn= seedn, hyp_grid=hyp_grid)
  }
  return(imp_mice)
}


### DAE gen m independent df one at a time
one_mice_imp_dae=function(df_in, mt, var_to_imp_list, n_iter, hyp_grid=NULL){
  print(paste('working on dataset m=', mt))
  dim0=dim(df_in)[2]-1
  #print(paste('dim0 s1: ', dim0))
  seedn=(9527+u)*mt
  imp_a_mice=chain_imp_dae(data_in = df_in, list_imp_var = var_to_imp_list, n_iter = n_iter, seedn= seedn, hyp_grid=hyp_grid, dim0 = dim0)
  return(imp_a_mice)
}


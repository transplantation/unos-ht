### In this file, we included all functions we wrote and used in the
### article: A Two-Stage Machine Learning Approach to Predict Heart 
### Transplantation Survival Probabilities over Time with a Monotonic 
### Probability Constraint


## cat_changer() function re-groups values in a categorical variable
## data_set: a data object that contains the varaible that we want to re-group
## var: the variable that we want to re-group
## val_old: some old levels in the variable
## val_new: some new levels in the variable
## val_old and val_new should be ONE-TO-ONE corresponds to each other
## Other levels that are not specified in val_old are put into "OTHER"
## Return value: a data object after re-groupping levels in the given variable
cat_changer <- function(data_set,var,val_old,val_new){
  temp_var <- dplyr::pull(data_set,var)
  cum_index <-c()
  for (i in 1:length(val_old)){
    index <- which(data_set[,var]==val_old[i])
    temp_var[index] <- val_new[i]
    cum_index <- c(cum_index,index)
  }
  na.index <- which(is.na(data_set[,var]))
  temp_var[-sort(c(cum_index, na.index))] <- "OTHER"
  data_set[,var] <- temp_var
  return(data_set)
}


## detect_terms() detects if the given term is in a vector
## x: a vector where we want to detect a given term
## term: a given term we want to detect
## Return Value: "YES" or "NO", "YES" means the given term appears at least once in the vector; "NO" means it doesn't appear in the vector
detect_terms <- function(x,term){
  result <- sapply(x, function(y) any(gsub(" ", "", strsplit(y, ",")[[1]])==term))
  if (all(is.na(result))) return(NA)
  result <- ifelse(is.na(result), FALSE, result)
  if (any(result==TRUE)) return("Y")
  return("N")
}


## col_missing_function() counts number of nas in each column in a given data object
## input_data: a array, including a matrix, dataframe
## Return Values: a data frame object that reconds percentage of missing values in each column in the input data
col_missing_function <- function(input_data){
  # first, we count nas in each column
  na_count_col <- apply(input_data, 2, function(y) length(which(is.na(y)))/nrow(input_data)) 
  # we saved the object into a data_frame
  na_count_col <- data.frame(na_count_col) 
  return(na_count_col)
}


## dummy_maker() performs One-hot encoding (creates dummy variables) for categorical variables in a given data object
## input_data: a data object
## char_var: a vector of all independent categorical variables in the data
## Return Values: a data object after One-hot encoding is applied
dummy_maker <- function(input_data,char_var){
  for (i in 1:ncol(input_data)){
    if(names(input_data[i]) %in% char_var){
      # Use createDummyvars() function to create dummy variables for each categorical variable
      # The definition of createDummyvars() can be found in this file
      temp <- createDummyvars(input_data[i])
      names(temp) <- paste(names(input_data[i]),levels(as.factor(input_data[,i])),sep="_")
      input_data <- cbind(input_data,temp)
      input_data[,ncol(input_data)]<-NULL}
  }
  # We removed the dependent dummy variable in each categorical variable
  input_data <- input_data[-which(names(input_data) %in% char_var)]
  return(input_data)
}


## createDummyvars() function creates dummy variables for a given categorical variable
## data0: a column (vector, categoric) in the data
## Return Values: a data object after one-hot encoding is applied to the given variable
createDummyvars <- function(data0){
  dd<-as.data.frame(table(data0))
  dum_data<-as.data.frame(matrix(0, ncol = nrow(dd), nrow = nrow(data0)))
  names(dum_data) <- dd[,1]
  for(i in 1:ncol(dum_data)){
    dum_data[which(data0==names(dum_data)[i]),i]<-1
    dum_data[i]<-as.factor(dum_data[,i])
  }
  return(dum_data)
}


## class_generator_bino() function creates binary response variable
## gstatus: GRAFT FAILED (1=YES)
## gtime: GRAFT LIFESPAN-Days From Transplant to Failure/Death/Last Follow-Up
## p_unit: time period, the unit is year
## predict_length: the length of days, 365 here if the unit is 1 year 
## Return Value: a vector (our response/dependent variable), 0: death, 1:survival, NA: unknow status
class_generator_bino <- function(gstatus,gtime,p_unit,predict_length){
  p_unit <- as.numeric(p_unit)
  predict_length <- as.numeric(predict_length)
  if(gtime < p_unit*predict_length){
    if(is.na(gstatus)){return(NA)}else{
      if(gstatus==0){return(NA)}
      if(gstatus==1){return(0)}  # death
    }
  }else{
    return(1)  # survival
  }
}


## FFS_bin() function performs Fast Feature Selection (FFS) for a given data object
## ii: a number that indicates which response variable (corresponds to the time point) should be used
## In our study, ii=0: first month, ii=1: 1st year, ii=2: 2nd year, ..., up to ii=10: 10th year 
## df: a data object (this should be the finalized data after one-hot encoding is applied)
## ids: a list object that we save corresponding ids for all traing and test data objects
## exclud: a vector of variables excluded when finding important variables
## seed: random seed used in the function 
## Return Value: a vector of all variables selected, saved it as a list
FFS_bin <- function(ii,df,ids,exclud,seed=110){
  library(Biocomb)
  set.seed(seed)
  df1 <- df[df$ID %in% ids[[paste("ID_train",ii, sep="")]],c(names(df)[
    !names(df) %in% exclud],paste0("year",ii))]
  df1[,paste0("year",ii)]<-as.factor(df1[,paste0("year",ii)])
  disc <- "MDL"
  threshold <- 0.001
  attrs.nominal <- numeric()
  FF_vars <- select.fast.filter(df1, disc.method=disc, threshold=threshold, attrs.nominal=attrs.nominal)
  FF_vars$Information.Gain<-NULL
  FF_vars$NumberFeature<-NULL
  names(FF_vars) <- "variables"
  newlist <- list(FF_vars)
  names(newlist) <- paste("Year", ii, sep="")
  # colnames(FFSV)<-c("variables")
  #return(FF_vars)
  return(newlist)
}


## Lasso_bin() function performs LASSO Feature Selection for a given data object
## ii: a number that indicates which response variable (corresponds to the time point) should be used
## In our study, ii=0: first month, ii=1: 1st year, ii=2: 2nd year, ..., up to ii=10: 10th year 
## df: a data object (this should be the finalized data after one-hot encoding is applied)
## ids: a list object that we save corresponding ids for all traing and test data objects
## exclud: a vector of variables excluded when finding important variables
## seed: random seed used in the function 
## folds, trace, alpha: parameters used in cv.glmnet() function, details can be found in the package: glmnet
## Return Value: a vector of all variables selected, saved it as a list
Lasso_bin <- function(ii,df,ids,exclud,folds=5,trace=F,alpha=1,seed=110){
  set.seed(seed)
  yvar <- paste("year", ii, sep="")
  df1 <- df[df$ID %in% ids[[paste("ID_train",ii,sep="")]],c(names(df)[
    !names(df) %in% exclud],yvar)]
  
  dff <- df1[!names(df1) %in% yvar]
  
  for(i in 1:ncol(dff)){
    dff[i] <- as.numeric(dff[,i])
  }
  x <- data.matrix(dff)
  glmnet1 <- glmnet::cv.glmnet(x=x,y=as.factor(df1[,yvar]),type.measure='auc',nfolds=folds,alpha=alpha, family="binomial")
  co <- coef(glmnet1,s = "lambda.1se")
  inds <- which(co[,1]!=0)
  variables <- row.names(co)[inds]
  variables <- as.data.frame(variables[!(variables %in% '(Intercept)')])
  colnames(variables) <- c("variables")
  newlist <- list(variables)
  names(newlist) <- paste("Year", ii, sep="")
  return(newlist)
}


## RF_bin() function performs Feature Selection using Randon Forest Algorithm in the package: Boruta for a given data object
## ii: a number that indicates which response variable (corresponds to the time point) should be used
## In our study, ii=0: first month, ii=1: 1st year, ii=2: 2nd year, ..., up to ii=10: 10th year 
## df: a data object (this should be the finalized data after one-hot encoding is applied)
## ids: a list object that we save corresponding ids for all traing and test data objects
## exclud: a vector of variables excluded when finding important variables
## seed: random seed used in the function 
## Return Value: a vector of all variables selected, saved it as a list
RF_bin <- function(ii, df, ids, exclud, seed=110){
  set.seed(seed)
  TARGET <- paste("year", ii, sep="")
  df1 <- df[df$ID %in% ids[[paste("ID_train",ii,sep="")]],c(names(df)[
    !names(df) %in% exclud],TARGET)]
  
  dataset <- df1[!names(df1) %in% TARGET]
  dataset$TARGET <- as.factor(df1[,TARGET])
  
  Random_Forrest.train <- Boruta::Boruta(TARGET~., data = dataset, doTrace = 2)
  Random_Forrest_fs <- as.data.frame(Random_Forrest.train$finalDecision)
  names(Random_Forrest_fs)[1] <- paste("criteria")
  Random_Forrest_imp <- as.data.frame.table(subset(Random_Forrest_fs, Random_Forrest_fs[,1] == "Confirmed"))
  names(Random_Forrest_imp)[1] <- paste("variables")
  names(Random_Forrest_imp)[2] <- c("version")
  R_Forrest <- as.data.frame(Random_Forrest_imp$variables)
  colnames(R_Forrest) <- c("variables")
  newlist <- list(R_Forrest)
  names(newlist) <- paste("Year", ii, sep="")
  return(newlist)
}


## all_bin() combines all varaibles selected from FFS, LASSO, and Random Forest
## x: time, x=0: 1st month, x=1: 1st year, ...x=10: 10th year
## df: a list object where all variables from each selection algorithm are saved
## Return Value: a vector of all variables combined, saved it as a list
all_bin <- function(x, df){
  all_bind <- data.frame(variables = union(union(df$FFS[[x]][[1]], df$LASSO[[x]][[1]]), df$RF[[x]][[1]]))
  newlist <- list(all_bind)
  names(newlist) <- paste("Year", (x-1), sep="")
  return(newlist)
}


## pred_func() performs down-sampling with replacement on the training data
## and trains models given a machine learning algorithm and the down-sampling samples
## it also evaluates the performance measures: AUC, Sensitivity, Specificity, Accuracy on the holdout set
## This function is used with parSapply() function in the package: snow
## All functions used in this function should be imported locally since multi-cores are used
## ii: the ii-th down-sampling sample is used
## traindata: the training data object
## hold_out: the hold_out / test data object
## TARGET: the response / dependent variable in the data
## formul: formula used in the machine learning algorithm, refer to the function train() in the package caret
## var_numeric: a vector of names for numerical variables
## assigned_seed: set a random seed
## methods_input="log": machine learning algorithm used
## methods_input can be one of the following: "log", "rf_bag", "gbm_boost", "cart_bag"
## fold_no,repeat_no: parameters used in trainControl() function in the package: caret
## The default fold_no is 5 and repeat_no is 3
## Return Values: a list contains the following:
## 1. Performance measures from the holdout object
## 2. Predicted Survival Probabilities for the patients in the holdout object
## 3. AUC value for the down-sampling sample
pred_func <- function(ii,traindata,hold_out,TARGET,formul,var_numeric,assigned_seed,methods_input="log",fold_no,repeat_no){
  library(caret)
  library(AUC)
  library(MASS)
  
  ## RUS_func() performs down-sampling with replacement algorithm
  ## input_data: a data frame object (should be the training data)
  ## TARGET: the response / dependent variable in the data
  ## Return Value: a down-sampling sample
  RUS_func <- function(input_data,TARGET){
    Train_Two <- input_data[ which(input_data[TARGET]=="Two"), ]
    Train_One <- input_data[ which(input_data[TARGET]=="One"), ]
    if(nrow(Train_Two)<=nrow(Train_One)){
      sample_size<-nrow(Train_Two)
      Train_One <- Train_One[sample(nrow(Train_One), sample_size, replace=T), ]
      Train_Two <- Train_Two[sample(nrow(Train_Two), sample_size, replace=T), ]
    }else{
      sample_size<-nrow(Train_One)
      Train_One <- Train_One[sample(nrow(Train_One), sample_size, replace=T), ]
      Train_Two <- Train_Two[sample(nrow(Train_Two), sample_size, replace=T), ]
    }
    input_data<-rbind(Train_One,Train_Two)
    
    return(input_data)
  }
  
  set.seed((assigned_seed+ii))
  
  traindata[TARGET] <- as.factor(ifelse(traindata[TARGET]==0, "One", "Two"))
  hold_out[TARGET] <- as.factor(ifelse(hold_out[TARGET]==0, "One", "Two"))
  
  ## we obtain the bootstrap down-sampling data using RUS_func() function and use it to train the model 
  traindata <- RUS_func(traindata,TARGET)
  
  ## we scale numerical variables in the training data and holdout data
  for (i in 1:ncol(traindata)){
    if (colnames(traindata)[i]%in%var_numeric){
      temp_normalized <- scale(traindata[,i])
      traindata[,i] <- as.numeric(temp_normalized)
      hold_out[,i] <- (hold_out[,i]-attributes(temp_normalized)$`scaled:center`)/attributes(temp_normalized)$`scaled:scale`
    }
  }
  
  control_ <- trainControl(method = "repeatedcv", number=fold_no,  
                           repeats = repeat_no, classProbs = TRUE)
  
  if(methods_input=="log"){
    result_model <- train(formul, data=traindata, method="glm", family="binomial",
                          trControl = control_, metric="Accuracy")
  }
  
  
  if(methods_input=="rf_bag"){
    #mtry <- sqrt(ncol(traindata))
    #tunegrid <- expand.grid(.mtry=mtry)
    result_model <- train(formul,  data=traindata, method="rf",
                          trControl = control_, tuneLength=10, metric="Accuracy")
  }
  
  if(methods_input=="gbm_boost"){
    result_model <- train(formul, data=traindata, method="gbm",
                          trControl = control_, tuneLength=10, metric="Accuracy")
  }
  
  if(methods_input=="cart_bag"){
    result_model <- train(formul, data=traindata, method="treebag", family="binomial",
                          trControl = control_, tuneLength=10, metric="Accuracy")
  }
  
  resul_raw <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(hold_out)))
  colnames(resul_raw) <- c("TARGET", methods_input, "Probability")
  resul_raw$TARGET <- hold_out[TARGET]
  
  train_raw <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(traindata)))
  colnames(train_raw) <- c("TARGET", methods_input, "Probability")
  train_raw$TARGET <- traindata[TARGET]
  
  resul_pred_perf<-as.data.frame(matrix(NA, ncol = 1, nrow = 4))
  colnames(resul_pred_perf)<-c(methods_input)
  rownames(resul_pred_perf)<-c("auc","sen","spec","accu")
  train_auc <- NA
  
  train_raw$Probability <- predict(result_model, newdata=traindata, type="prob")[,2]
  train_raw[methods_input] <- predict(result_model, newdata=traindata, type="raw")
  train_auc <- AUC::auc(roc(train_raw$Probability, traindata[,TARGET]))
  resul_raw[methods_input] <- predict(result_model, newdata=hold_out, type="raw")
  resul_raw$Probability <- predict(result_model, newdata=hold_out, type="prob")[,2]
  resul_pred_perf[1,1] <- AUC::auc(roc(resul_raw$Probability,hold_out[,TARGET]))
  resul_pred_perf[2,1] <- caret::sensitivity(resul_raw[,methods_input],hold_out[,TARGET])
  resul_pred_perf[3,1] <- caret::specificity(resul_raw[,methods_input],hold_out[,TARGET])
  resul_pred_perf[4,1] <- (as.data.frame(confusionMatrix(resul_raw[,methods_input], hold_out[,TARGET])$overall))[1,]
  
  return(list(Performance=resul_pred_perf, Predicted=resul_raw, AUC=train_auc))
}


#The following two ________________________________________________________________________

# check_range: check if numerical values are in the corresponding range
check_range <- function(df, vars, Num_Range){
  out_range <- rep(NA, 13)
  for (i in 1:length(vars)){
    m <- min(df[,vars[i]])
    M <- max(df[,vars[i]])
    out_range[i] <- ifelse((m<Num_Range[i,1]|M>Num_Range[i,2]), 1, 0)
  }
  return(out_range)
}
#________________________________________________________________________

# check_levels: check if categorical values are in our data
check_levels <- function(df, vars, Cat_Levels){
  CAT_match <- rep(NA, length(vars))
  for (i in 1:length(vars)){
    temp_level <- levels(as.factor(df[,vars[i]]))
    CAT_match[i] <- ifelse(all(temp_level%in%Cat_Levels[[i]]), 0, 1)
  }
  return(CAT_match)
}


#________________________________________________________________________

# creating dummy variables for input data (this is for shinny app)
dummy_maker_app<-function(input_data, char_var){
  for (i in 1:ncol(input_data)){
    if(names(input_data[i]) %in% char_var){
      temp<-createDummyvars(input_data[i])
      names(temp)<-paste(names(input_data[i]),levels(as.factor(input_data[,i])),sep="_")
      
      input_data<-cbind(input_data,temp)
      if (length(levels(as.factor(input_data[,i])))!=1){
        input_data[,ncol(input_data)]<-NULL
      }else{
        input_data[,ncol(input_data)]<-factor(input_data[,ncol(input_data)], levels = c(0,1))
      }
    }
  }
  input_data<-input_data[-which(names(input_data) %in% char_var)]
  return(input_data)
}

#________________________________________________________________________

# change all "U" to Unknown (this is for shinny app)
find_U <- function(x){
  col_index <- which(x=="U")
  if (length(col_index)!=0) x[col_index] <- "UNKOWN"
  return(x)
}

# finding mode for categorical variables

cahrmode <- function(x) {
  x<-x[complete.cases(x)]
  uniq_val <- unique(x)
  return(uniq_val[which.max(tabulate(match(x, uniq_val)))])
}

#==================================================================
# a function for finding proportion of missings in the rows
row_missing_function <- function(input_data){
  na_count_row <- apply(input_data, 1, function(y) length(which(is.na(y)))/ncol(input_data)) # counting nas in each row
  na_count_row <- data.frame(na_count_row) #transforming this into a data_frame
  return(na_count_row)
}
# a function for finding proportion of missings in the columns
col_missing_function <- function(input_data){
  na_count_col <- apply(input_data, 2, function(y) length(which(is.na(y)))/nrow(input_data)) # counting nas in each column
  na_count_col <- data.frame(na_count_col) #transforming this into a data_frame
  return(na_count_col)
}
#==================================================================

# function for cleaning tables
# Try_data<-heart.df.cleaned
# kept_opt<-c("FUNC_STAT_TRR","FUNC_STAT_TCR")
# ID<-"ID"
# kept_col<-c("year1")
# col2row_emp<-1
table_cleaner<-function(Try_data,col2row_emp,ID,kept_col,kept_opt){
  
  data0<-Try_data
  Try_data<-data0[complete.cases(data0[kept_col]),]
  
  max_col<-0
  max_row<-0
  Try_data_temp <<- Try_data
  
  for(i in 1:100000){
    Try_data<-Try_data_temp
    gc()
    count_col <- col_missing_function(Try_data)
    gc()
    count_row <- row_missing_function(Try_data)
    max_col<-max(count_col)
    max_row<-max(count_row)
    max_emp<- max(max_col,max_row)
    
    if(max_col==0){break()}
    if(max_row==0){break()}
    
    if (max_emp==max_row)
    {
      if((nrow(Try_data))>3){
        a<-which(count_row$na_count_row==max_row)
        
        Try_data_temp <<- Try_data[-a,]
        print(c("nrow",nrow(Try_data_temp)))
      }
    }
    if(max_col>=max_row*col2row_emp){
      b<-names(Try_data[which(count_col$na_count_col ==max_col)])
      
      if (max_emp==max_col){
        if((ncol(Try_data))>2){
          #in the next lines, I try to find which columns are the emptiest and then take out the one that I don't want to be deleted
          
          f<-which(count_col$na_count_col==max_col)
          print(paste("iteration: ",i," // dropped column name: ",names(Try_data[f]),sep=""))
          Try_data_temp<<- Try_data[, -f]
        }
        
      }
      
      print(c("ncol",ncol(Try_data_temp)))
    }else{
      if((nrow(Try_data))>3){
        a<-which(count_row$na_count_row==max_row)
        
        Try_data_temp <<- Try_data[-a,]
        
      }
      print(c("nrow",nrow(Try_data_temp)))
    }
    
  }
  
  
  input_object<-list()
  
  print("Data Cleaning is Done! Exit the tool. Run the tool again, load the cleaned file for Data Analysis")
  out_object<-list()
  out_object$col_names<-names(Try_data_temp)
  out_object$row_ID<-Try_data_temp[ID]
  
  if(sum(names(data0) %in% kept_opt)>0){
    temp<-data0[out_object$row_ID[,],which(names(data0) %in% kept_opt)]
    if(sum(names(data0) %in% kept_opt)==1){
      temp<-as.data.frame(temp)
      names(temp)<-names(data0[kept_opt])}
    
    data_inc<-cbind(Try_data_temp,temp)
    data_inc<-data_inc[complete.cases(data_inc[kept_opt]),]
    out_object$col_names_inc<-names(data_inc)
    out_object$row_ID_inc<-data_inc[ID]
  }else{
    out_object$col_names_inc<-out_object$col_names
    out_object$row_ID_inc<-out_object$row_ID
  }
  out_object$data<-Try_data_temp
  return(out_object)
}
#==================================================================
# fast feature selection filter method

FFS_features<-function(ii,data,dependent,exclud, seed=110){
  
  set.seed(seed)
  exclud<- c(exclud,dependent)
  df1<-data[[ii]]
  
  data_use<-cbind(df1[!names(df1)%in%exclud],df1[[dependent]])
  
  disc<-"MDL"
  threshold=0.001
  attrs.nominal=numeric()
  FF_vars=Biocomb::select.fast.filter(data_use, disc.method=disc, threshold=threshold,
                                      attrs.nominal=attrs.nominal)
  
  FF_vars$Information.Gain<-NULL
  FF_vars$NumberFeature<-NULL
  names(FF_vars) <- "variables"
  newlist <- list(FF_vars)
  names(newlist) <- paste("scenario", ii, sep="")
  # colnames(FFSV)<-c("variables")
  #return(FF_vars)
  return(newlist)
}

#==================================================================

# feaature selection using LASSO algorithm
LASSO_features <- function(ii,data,dependent,exclud,folds=5,trace=F,alpha=1,seed=110){
  set.seed(seed)
  excluds<-c(exclud,dependent)
  dataset<-data[[ii]]
  dff<-dataset[!names(dataset)%in%excluds]
  
  for(i in 1:ncol(dff)){
    dff[i] <- as.numeric(dff[,i])
  }
  x <- data.matrix(dff)
  glmnet1 <- glmnet::cv.glmnet(x=x,y=as.factor(dataset[,dependent]),type.measure='auc',nfolds=folds,alpha=alpha, family="binomial")
  co <- coef(glmnet1,s = "lambda.1se")
  inds <- which(co[,1]!=0)
  variables <- row.names(co)[inds]
  variables <- as.data.frame(variables[!(variables %in% '(Intercept)')])
  colnames(variables) <- c("variables")
  newlist <- list(variables)
  names(newlist) <- paste("scenario", ii, sep="")
  return(newlist)
}


#==================================================================

# data<-dfs_dum[[1]]
# dependent<-"year1"
# exclud<-"ID"

# feaature selection using Random Forrest algorithm
RF_features <- function(ii,data,dependent,exclud, seed=110){
  set.seed(seed)
  exclud<- c(exclud,dependent)
  df1<-data[[ii]]
  
  dataset<-df1[!names(df1)%in%exclud]
  
  dataset$TARGET <- as.factor(df1[,dependent])
  
  Random_Forrest.train <- Boruta::Boruta(TARGET~., data = dataset, doTrace = 2)
  Random_Forrest_fs <- as.data.frame(Random_Forrest.train$finalDecision)
  names(Random_Forrest_fs)[1] <- paste("criteria")
  Random_Forrest_imp <- as.data.frame.table(subset(Random_Forrest_fs, Random_Forrest_fs[,1] == "Confirmed"))
  names(Random_Forrest_imp)[1] <- paste("variables")
  names(Random_Forrest_imp)[2] <- c("version")
  R_Forrest <- as.data.frame(Random_Forrest_imp$variables)
  colnames(R_Forrest) <- c("variables")
  newlist <- list(R_Forrest)
  names(newlist) <- paste("scenario", ii, sep="")
  return(newlist)
}




## encode_cat() retuns the data after the assigned encoding method for
## categorical variables
## df: data; method = c("numeric", "factor")
encode_cat <- function(df, cat_vars, method){
  cat_index <- which(colnames(df)%in%cat_vars)
  if (method=="numeric"){
    df[,cat_index] <- apply(df[,cat_index], 2, function(x) as.numeric(as.factor(x)))
  }else if (method=="factor"){
    variables <- colnames(df)[cat_index]
    df <- dummy_maker(df, variables)
  }
  return(df)
}


## testdata_index() creates IDs for test data
## n: number of rows from the entire dataset
## seed: random seed used in the function, the default is 2019
holdout_index <- function(IDs, seed=2019){
  n<-length(IDs)
  set.seed <- seed
  index <- list()
  index0 <- c()
  
  for (i in 1:5){
    index[[i]] <- sample(setdiff(IDs, index0), floor(n/5))
    index0 <- c(index0, index[[i]])
  }
  return(index)
} 

## select_vars() does the variable selection for the data
## three variable selection methods are available: FFS, LASSO, Random Forest
## ii: the index to find the iith training data (we have 5 disjoint holdout datasets)
## df: the whole dataset
## y: name of the response variable
## method: FFS or LASSO or RF
## ID: the ID for the test data set 
## seed: random seed used in the function, the default is 201905

select_vars <- function(ii, df, y, method, ID, seed=201905){
  set.seed <- seed
  df <- df[!df$ID%in%ID[[ii]],]
  df <- df[,colnames(df)!="ID"]
  X <- df[,colnames(df)!=y]
  if (method=="FFS"){
    X$y <- df[,y]
    disc <- "MDL"
    threshold <- 0.001
    attrs.nominal <- numeric()
    result_temp <- Biocomb::select.fast.filter(X, disc.method=disc, threshold=threshold, attrs.nominal=attrs.nominal)
    Imp_vars <- as.character(result_temp$Biomarker)
  }else if (method=="LASSO"){
    '%ni%' <- Negate('%in%')
    X <- data.matrix(apply(X, 2, as.numeric))
    glmnet1 <- glmnet::cv.glmnet(x=X,y=as.factor(df[,y]),type.measure='auc', family="binomial")
    co <- coef(glmnet1,s = "lambda.1se")
    co <- as.matrix(co[-1,])
    Imp_vars <- row.names(co)[which((co[,1]!=0))]
    Imp_vars <- Imp_vars[Imp_vars %ni% '(Intercept)']
  }else if (method=="RF"){
    library(party)
    library(varImp)
    # formul<- as.formula(paste0(y,"~."))
    # cf1 <- Boruta::Boruta(formul, data = df, doTrace = 2)
    # find_vars <- cf1$finalDecision[which(cf1$finalDecision!="Rejected")]
    # Imp_vars <- names(find_vars)
    # for parallel processing I have to load the functions
    source("https://raw.githubusercontent.com/transplantation/DoE/Hamid/DoE_paper_functions.R")
    y<-as.factor(df[[y]])
    Imp_vars <- var.sel.r2vim(X, y, no.runs = 1,ntree = 500, random_seed = seed)
    Imp_vars<-Imp_vars$var
  }
  return(list(Imp_vars))
}


## modeling() utilitizes a machine learning model in the caret package using cv with 5 folds
## The available models are glm, nnet, svmRadial, rf, gbm, earth, rpart, xgbTree, naive_bayes,
## treebag, lda, ranger
## it also evaluates the performance measures: AUC, Sensitivity, Specificity, Accuracy on the holdout set
## This function is used with parSapply() function in the package: snow
## All functions used in this function should be imported locally since multi-cores are used
## ii: the ii-th data indicator for different row in Index_matrix
## All_df: a list object contains all data with difference scenarios (16 components)
## TARGET: name of the response / dependent variable in the data
## Index_matrix: a matrix I(i,j) that helps to indicate which scenairo (j) and the corresponding 
## train and holdout data (i) and features used (j,i)
## assigned_seed: set a random seed
## methods_input="glm" : default
## Return Values: a list contains the following:
## 1. Performance measures from the holdout object
## 2. Predicted Survival Probabilities for the patients in the holdout object


model_data<-function(aa,All_df,TARGET,Index_matrix,Index_test,features,methods_input="none", 
                     sampling.method=c("none", "down", "up", "rose", "smote"),iteration=5 ,
                     data_experiments,assigned_seed1=2019){
  gc()
  
  library(pacman) # needs to be installed first
  # p_load is equivalent to combining both install.packages() and library()
  p_load(caret,AUC,MASS,ROSE,DMwR,snow,ranger,parallel)
  
  All_df1<-All_df[aa]
  All_data.scenario1<-aa
  
  
  modeling <- function(ii,All_df2,TARGET2,Index_matrix2,Index_test2,features2,methods_input2="none", 
                       sampling.method2=c("none", "down", "up", "rose", "smote"),All_data.scenario2 ,assigned_seed2=2019){
    gc()
    library(pacman) # needs to be installed first
    # p_load is equivalent to combining both install.packages() and library()
    p_load(caret,AUC,MASS,ROSE,DMwR,snow,ranger,parallel)
    
    test_no <- Index_matrix2[ii,1]
    impute_no <- Index_matrix2[ii,2]
    df <- All_df2[[impute_no]]
    index_t<- which(df$ID%in%Index_test2[[test_no]])
    #df <- df[,c(features[[impute_no]][[test_no]], TARGET2)]
    df <- df[,which(names(df)%in% c((features2[[impute_no[1]]][test_no[1]])[[1]],as.character(TARGET2[1])) )]
    hold_out <- df[index_t,]
    traindata <- df[-index_t,]
    traindata$ID <- NULL
    fold2<-ii
    
    set.seed(assigned_seed2+ii)
    
    # 1: survival; 0: death
    traindata[,as.character(TARGET2[1])] <- as.factor(ifelse(traindata[,as.character(TARGET2[1])]==0, "Death", "Survival"))
    hold_out[,as.character(TARGET2[1])] <- as.factor(ifelse(hold_out[,as.character(TARGET2[1])]==0, "Death", "Survival"))
    
    formul2<- as.formula(paste0(as.character(TARGET2[1]),"~."))
    
    trainers<-function(iii,formul3,traindata3,hold_out3,methods_input3="none",All_data.scenario3=All_data.scenario2,
                       fold3,TARGET3,assigned_seed3=assigned_seed2){
      gc()
      set.seed(assigned_seed3)
      
      library(pacman) # needs to be installed first
      # p_load is equivalent to combining both install.packages() and library()
      p_load(caret,AUC,MASS,ROSE,DMwR,snow,ranger)
      
      sampling.method<-iii
      if(iii=="none"){sampling.method<-NULL}
      formul<-formul3
      
      
      # I used 5 fold cross validation 
      control_setting <- caret::trainControl(method = "cv", number=5, sampling=sampling.method , summaryFunction = twoClassSummary, 
                                             search="random", classProbs = TRUE, selectionFunction="tolerance")
      
      if (methods_input3%in% c("glm", "nnet", "svmRadial")){
        result_model <- train(formul, data=traindata3, method=methods_input3, family="binomial",
                              trControl = control_setting, metric="ROC")
      }else if (methods_input3%in%c("rf", "gbm", "earth", "rpart", "xgbTree", "naive_bayes","xgbDART" ,"ranger")){
        result_model <- train(formul,  data=traindata3, method=methods_input3,
                              trControl = control_setting, tuneLength=10, metric="ROC")
      }else if(methods_input3=="lda"){
        result_model <- train(formul, data=traindata3, method=methods_input3, preProcess="pca", preProcOptions = list(method="BoxCox"),
                              trControl = control_setting, metric="ROC")
      }else if(methods_input3=="treebag"){
        result_model <- train(formul, data=traindata3, method=methods_input3, family="binomial",
                              trControl = control_setting, tuneLength=10, metric="ROC")
      }
      
      resul_raw <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(hold_out)))
      colnames(resul_raw) <- c("TARGET", methods_input3, "Probability")
      resul_raw$TARGET <- hold_out[as.character(TARGET3[1])]
      
      train_raw <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(traindata3)))
      colnames(train_raw) <- c("TARGET", methods_input3, "Probability")
      train_raw$TARGET <- traindata3[as.character(TARGET3[1])]
      
      resul_pred_perf<-as.data.frame(matrix(NA, ncol = 1, nrow = 4))
      colnames(resul_pred_perf)<-c(methods_input3)
      rownames(resul_pred_perf)<-c("auc","sen","spec","accu")
      train_auc <- NA
      
      #if(class(try(varImp(result_model),silent = TRUE))!="try-error"){
      train_raw$Probability <- predict(result_model, newdata=traindata3, type="prob")[,2]
      train_raw[methods_input3] <- predict(result_model, newdata=traindata3, type="raw")
      #train_auc <- AUC::auc(roc(train_raw$Probability, traindata[,TARGET]))
      resul_raw[methods_input3] <- predict(result_model, newdata=hold_out, type="raw")
      resul_raw$Probability <- predict(result_model, newdata=hold_out, type="prob")[,2]
      resul_pred_perf[1,1] <- AUC::auc(roc(resul_raw$Probability,hold_out[,as.character(TARGET3[1])]))
      resul_pred_perf[2,1] <- caret::sensitivity(resul_raw[,methods_input3],hold_out[,as.character(TARGET3[1])])
      resul_pred_perf[3,1] <- caret::specificity(resul_raw[,methods_input3],hold_out[,as.character(TARGET3[1])])
      resul_pred_perf[4,1] <- (as.data.frame(confusionMatrix(resul_raw[,methods_input3], hold_out[,as.character(TARGET3[1])])$overall))[1,]
      #}
      
      
      # putting summary of the experiment in here
      experiment.summary<-as.data.frame(matrix(0, ncol = 5, nrow = 1))
      names(experiment.summary)<-c("data_scenario","training_algorithm","resampling_method","fold","TARGET")
      experiment.summary$data_scenario<-All_data.scenario3
      experiment.summary$training_algorithm<-methods_input3
      experiment.summary$resampling_method<-iii
      experiment.summary$fold<-fold3
      experiment.summary$TARGET<-as.character(TARGET3[1])
      
      return(list(experiment.summary=experiment.summary, Performance=resul_pred_perf, Predicted=resul_raw
      ))
      
    }
    
    
    
    cl <- makeCluster(length(sampling.method2), type="SOCK")
    modeling_result <- parSapply(cl, sampling.method2, trainers ,formul3=formul2,traindata3=traindata,hold_out3=hold_out,
                                 methods_input3=methods_input2,fold3=fold2,TARGET3=TARGET2,assigned_seed3=assigned_seed2) 
    stopCluster(cl)
    
    return(modeling_result)
    
    
  }
  
  cl.main <- makeCluster(length(data_experiments), type="SOCK")
  modeling_result.main <- parSapply(cl.main, c(1:iteration), modeling ,All_df2=All_df1,TARGET2=TARGET,Index_matrix2=Index_matrix,
                                    Index_test2=Index_test,features2=features,methods_input2=methods_input, 
                                    sampling.method2=sampling.method, All_data.scenario2=All_data.scenario1, assigned_seed2=2019) 
  stopCluster(cl.main)
  
  return(modeling_result.main)
  
}

#' Variable selection using recurrent relative variable importance (r2VIM).
#'
#' Generates several random forests using all variables and different random
#' number seeds. For each run, the importance score is divided by the (absolute)
#' minimal importance score (relative importance scores). Variables are selected
#' if the minimal relative importance score is >= factor.
#'
#' Note: This function is a reimplementation of the R package \code{RFVarSelGWAS}.
#'
#' @inheritParams wrapper.rf
#' @param no.runs number of random forests to be generated
#' @param factor minimal relative importance score for a variable to be selected
#'
#' @return List with the following components:
#'   \itemize{
#'   \item \code{info} data.frame
#'   with information for each variable
#'   \itemize{
#'   \item vim.run.x = original variable importance (VIM) in run x
#'   \item rel.vim.run.x = relative VIM in run x
#'   \item rel.vim.min = minimal relative VIM over all runs
#'   \item rel.vim.med = median relative VIM over all runs
#'   \item selected = variable has been selected
#'   }
#'   \item \code{var} vector of selected variables
#'   }
#'
#'  @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # select variables
#' res = var.sel.r2vim(x = data[, -1], y = data[, 1], no.runs = 5, factor = 1)
#' res$var
#'
#' @export

var.sel.r2vim <- function(x, y, no.runs = 1, factor = 1, ntree = 500, 
                          mtry.prop = 0.2, nodesize.prop = 0.1,
                          no.threads = 1, method = "ranger", 
                          type = "classification",random_seed =seed) {
  set.seed(random_seed)
  
  ## importance for each run
  imp.all = NULL
  for (r in 1:no.runs) {
    print(paste("run", r))
    rf = wrapper.rf(x = x, y = y,
                    ntree = ntree, mtry.prop = mtry.prop, nodesize.prop = nodesize.prop, no.threads = no.threads,
                    method = method, type = type)
    imp.all = cbind(imp.all, get.vim(rf))
  }
  
  ## factors
  min.global = min(imp.all)
  if (min.global >= 0) {
    stop("Global minimal importance score is not negative!")
  }
  no.neg.min = 0
  fac = matrix(nrow = nrow(imp.all), ncol = ncol(imp.all),
               dimnames = dimnames(imp.all))
  for (i in 1:ncol(imp.all)) {
    x = imp.all[,i]
    min = min(x)
    if (min >= 0) {
      no.neg.min = no.neg.min + 1
      fac[,i] = x / abs(min.global)
    } else {
      fac[, i] = x / abs(min)
    }
  }
  if (no.neg.min > 0) {
    print(paste(no.neg.min, "runs with no negative importance score!"))
  }
  fac.min = apply(fac, 1, min)
  fac.med = apply(fac, 1, median)
  
  ## select variables
  ind.sel = as.numeric(fac.min >= factor)
  
  ## info about variables
  info = data.frame(imp.all, fac, fac.min, fac.med, ind.sel)
  colnames(info) = c(paste("vim.run.", 1:no.runs, sep = ""),
                     paste("rel.vim.run.", 1:no.runs, sep = ""),
                     "rel.vim.min", "rel.vim.median", "selected")
  return(list(info = info, var = sort(rownames(info)[info$selected == 1])))
}

#' Wrapper function to call random forests function.
#'
#' Provides an interface to different parallel implementations of the random
#' forest algorithm. Currently, only the \code{ranger} package is
#' supported.
#'
#' @param x matrix or data.frame of predictor variables with variables in
#'   columns and samples in rows (Note: missing values are not allowed).
#' @param y vector with values of phenotype variable (Note: will be converted to factor if
#'   classification mode is used).
#' @param ntree number of trees.
#' @param mtry.prop proportion of variables that should be used at each split.
#' @param nodesize.prop proportion of minimal number of samples in terminal
#'   nodes.
#' @param no.threads number of threads used for parallel execution.
#' @param method implementation to be used ("ranger").
#' @param type mode of prediction ("regression", "classification" or "probability").
#' @param ... further arguments needed for \code{\link[relVarId]{holdout.rf}} function only.
#'
#' @return An object of class \code{\link[ranger]{ranger}}.
#'
#' @import methods stats
#'
#' @export
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # regression
#' wrapper.rf(x = data[, -1], y = data[, 1],
#'            type = "regression", method = "ranger")

wrapper.rf <- function(x, y, ntree = 100, mtry.prop = 0.2, nodesize.prop = 0.1, no.threads = 1,
                       method = "ranger", type = "regression", ...) {
  
  ## check data
  if (length(y) != nrow(x)) {
    stop("length of y and number of rows in x are different")
  }
  
  if (any(is.na(x))) {
    stop("missing values are not allowed")
  }
  
  if (type %in% c("probability", "regression") & (is.character(y) | is.factor(y))) {
    stop("only numeric y allowed for probability or regression mode")
  }
  
  ## set global parameters
  nodesize = floor(nodesize.prop * nrow(x))
  mtry = floor(mtry.prop * ncol(x))
  if (mtry == 0) mtry = 1
  
  if (type == "classification") {
    #    print("in classification")
    y = as.factor(y)
  }
  
  ## run RF
  if (method == "ranger") {
    if (type == "probability") {
      y = as.factor(y)
      prob = TRUE
    } else {
      prob = FALSE
    }
    
    rf = ranger::ranger(data = data.frame(y, x),
                        dependent.variable.name = "y",
                        probability = prob,
                        importance = "permutation", scale.permutation.importance = FALSE,
                        num.trees = ntree,
                        mtry = mtry,
                        min.node.size = nodesize,
                        num.threads = no.threads,
                        write.forest = TRUE,
                        ...)
  } else {
    stop(paste("method", method, "undefined. Use 'ranger'."))
  }
  
  return(rf)
}


#' Error calculation.
#'
#' Calculates errors by comparing predictions with the true values. For
#' regression and probability mode, it will give root mean squared error (rmse) and
#' pseudo R-squared (rsq). For classification mode, overall accuracy (acc), overall
#' error (err), Matthews correlation coefficient (mcc), sensitivity (sens) and
#' specificity (spec) are returned.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}
#' @param true vector with true value for each sample
#' @param test.set matrix or data.frame of predictor variables for test set with variables in
#'   columns and samples in rows (Note: missing values are not allowed)
#' @inheritParams wrapper.rf
#'
#' @return numeric vector with two elements for regression and probability estimation (rmse, rsq) and
#' five elements for classification (acc, err, mcc, sens, spec)
#'
#' @export
#'
#' @examples
#' # simulate toy data set
#' data = simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
#'
#' # random forest
#' rf = wrapper.rf(x = data[, -1], y = data[, 1],
#'                 type = "regression")
#'
#' # error
#' calculate.error(rf = rf, true = data[, 1])

calculate.error <- function(rf, true, test.set = NULL) {
  
  if (is(rf, "ranger")) {
    if (!is.null(test.set)) {
      pred = predict(rf, data = test.set)$predictions
    } else {
      pred = rf$predictions
    }
    if (rf$treetype == "Probability estimation") {
      pred = pred[, 2]
    }
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  
  if ((is(rf, "randomForest") && rf$type == "classification") |
      (is(rf, "ranger") && rf$treetype == "Classification")) {
    conf.matrix = table(pred = pred, true = true)
    tp = conf.matrix[2, 2]
    tn = conf.matrix[1, 1]
    fn = conf.matrix[2, 1]
    fp = conf.matrix[1, 2]
    
    ## accuracy
    acc = (tp + tn) / sum(conf.matrix)
    
    ## Matthews correlation coefficient
    mcc = (tp * tn - fp * fn) /
      sqrt( (tp + fn) * (tn + fp) * (tp + fp) * (tn + fn))
    
    ## sensitivity
    sens = tp / (tp + fn)
    
    ## specificity
    spec = tn / (fp + tn)
    
    error = c(err = 1 - acc, acc = acc, mcc = mcc, sens = sens, spec = spec)
  } else {
    mse = sum((pred - true)^2, na.rm = TRUE) / sum(!is.na(pred))
    
    ## pseudo R-squared uses sum of squared differences divided by n instead of variance!
    v = sum((true - mean(true))^2) / length(true)
    rsq = 1 - mse/v
    error = c(rmse = sqrt(mse), rsq = rsq)
  }
  
  return(error)
}


#' Get variable importance.
#'
#' Extracts variable importance depending on class of random forest object.
#'
#' @param rf Object of class \code{\link[ranger]{ranger}}
#'
#' @return numeric vector with importance value for each variable (in original order)
#'
#' @export

get.vim <- function(rf) {
  if (is(rf, "ranger")) {
    vim = ranger::importance(rf)
  } else {
    stop(paste("rf needs to be of class ranger"))
  }
  return(vim)
}
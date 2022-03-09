### Custom functions 

detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}



#map2_df with adding cohort info and subjectID at once

cohort_rbindlist<-function(list){
  list %>% 
    map2_df(.,names(.),~mutate(.x,cohort=.y)) %>% 
    mutate(subjectID=paste(cohort,mouseID,sep="."))
}

## get tidied up data frame - winner/loser format 
get_df_new_ethogram<-function(dfx){
  #Add Timestamp, Hours, Minutes, Seconds

  dfx1<-dfx %>% 
    filter(!is.na(Timestamp)) %>% 
    mutate(date = as.Date(Timestamp, format="%m/%d/%Y %H:%M")) %>% 
    mutate(day = as.numeric(date - min(date))+1,
           hour = strptime(Timestamp, format="%m/%d/%Y %H:%M") %>% lubridate::hour(.)-11,
           minute = strptime(Timestamp, format="%m/%d/%Y %H:%M") %>% lubridate::minute(.)
           ) %>% 
    arrange(day,hour,minute) #arrange by date and time
  
  #Remove Duplicates/Ties that are not Start/End
  dfx1$`Animal 1 - Initiator` <- as.character(dfx1$`Animal 1 - Initiator`)
  dfx1$`Animal 2`<- as.character(dfx1$`Animal 2`)
  dfx1 <- dfx1[(dfx1$`Animal 1 - Initiator`!=dfx1$`Animal 2`) | (dfx1$`Start/End`=="End") | (dfx1$`Start/End`=="Start"),]  
  
  # Add observation sample number:
  dfx2 <- dfx1 %>% 
    mutate(`Start/End`=ifelse(is.na(`Start/End`),"",`Start/End`)) %>% 
    mutate(uniqueobs = cumsum(`Start/End`=="Start"))
  
  # Add time variable 
  dfx3<-dfx2 %>% 
    mutate(time=hour*60 + minute )  
  
  #get rid of start/end
  dfx4<-dfx3 %>%
    filter(`Start/End`=="")
  
  dfx5a<-dfx4 %>% filter(result==1)
  dfx5b<-dfx4 %>% filter(result==-1) %>% mutate(result=result*(-1))
  
  colnames(dfx5a)[c(2,4)]<-c("winner","loser")
  colnames(dfx5b)[c(2,4)]<-c("loser","winner")  
  
  dfx6 <- rbind(dfx5a,dfx5b) %>% 
    arrange(day,hour,minute) 
  
  return(dfx6)  
}

get_response_freq<-function(df){
  df1<-df %>%
    filter(!is.na(Timestamp)) %>% 
    mutate(date = as.Date(Timestamp, format="%m/%d/%Y %H:%M")) %>% 
    mutate(day = as.numeric(date - min(date))+1)%>% 
    arrange(Timestamp) %>% #arrange by date and time
    select(day,winner,loser,`Animal 2 - Behavior`) 
  
  colnames(df1)<-c("GD","aggressorID","mouseID","response")
  
  df2<-df1 %>% 
    filter(!is.na(response)) %>% 
    mutate(preem_flee = ifelse(str_detect(response,"Preemptive")==T,"Flee",NA),
           flee = ifelse(str_detect(response,"Flee")==T,"Flee",NA), #it distinguishes from preemtive flee because of capital F
           sub_pos = ifelse(str_detect(response,"posture")==T,"Sub",NA),
           freeze = ifelse(str_detect(response,"Freeze")==T,"Freeze",NA)) %>% 
    select(-response)
  
  df3<-df2 %>% 
    gather(col,response,4:7) %>% 
    select(GD,aggressorID,mouseID,response) %>% 
    filter(!is.na(response)) %>% 
    mutate(GD=as.numeric(GD)) %>% 
    arrange(GD)
  
  return(df3)}



##get win/loss matrix with new ethogram
get_wl_matrix_new<-function(df){
mylevs = unique(c(as.character(df[,1]),as.character(df[,2])))

df[,1] <- factor(df[,1], levels=mylevs)
df[,2] <- factor(df[,2], levels=mylevs)
m1 = stats::xtabs(result ~ ., data = df)
m1 <- m1[order(rownames(m1)), order(colnames(m1))]
return(m1)
}


# Calculate despotism
despotism <- function(x) {
  rev(sort(round(100*(rowSums(x)/sum(x)),2)))
}


# Calculate p-value from steepness test
getStp.pval <- function(x){
  a <- steepness::steeptest(x,rep=1000)
  return( sum(a$Stpsim>=a$Stp)/1000 )
}


# Make dataframe out of David's Scores
get_dsdf <- function(d){
  d.df <- data.frame(rank=1:length(d), ds = rev(d[order(d)]))
  return(d.df)
}



get_glickos <- function(aa, cval=3){
  PlayerRatings::glicko(aa %>% 
      mutate(event=row_number()) %>% 
      select(event,winner,loser,result),
    history=T,plot=F,cval=cval)
}




#count how many times each animal initiated social/aggressive interactions

initiation<-function(df){
  dfx<-df %>% 
    count(`Animal 1 - Initiator`) %>% 
    rename(mouseID=`Animal 1 - Initiator`,
           initiation= n) #I feel like an idiot not knowing rename function so far!
  
  dfy<-df %>% 
    filter(result>0) %>% 
    count(`Animal 1 - Initiator`) %>% 
    rename(mouseID=`Animal 1 - Initiator`,
           agg_initiation= n) 
  
  dfz<-left_join(dfx,dfy)
  return(dfz)
}

#get the dataframe 1/0 for each behavior for each recorded observation data point 
each_behavior<-function(df){
  
  df0 <- df %>% 
    filter(is.na(`Start/End`))
  
  df1a <- df0 %>% select(Timestamp, `Animal 1 - Initiator`, `Animal 1 - Behavior`)
  df1b <- df0 %>% select(Timestamp, `Animal 2`, `Animal 2 - Behavior`)
  
  colnames(df1a)<-c("Timestamp","mouseID","behavior")
  colnames(df1b)<-c("Timestamp","mouseID","behavior")
  
  df2<-rbind(df1a,df1b) %>% 
    arrange(Timestamp)
  
  df3 <- df2 %>% 
    filter(!is.na(behavior)) %>% 
    mutate(lunge = ifelse(str_detect(behavior,"Lunge")==T,1,0),  
           fight = ifelse(str_detect(behavior,"Fighting")==T,1,0),
           chase = ifelse(str_detect(behavior,"Chase")==T,1,0),
           mount = ifelse(str_detect(behavior,"Mount")==T,1,0),
           preem_flee = ifelse(str_detect(behavior,"Preemptive")==T,1,0),
           flee = ifelse(str_detect(behavior,"Flee")==T,1,0), #it distinguishes from preemtive flee because of capital F
           sub_pos = ifelse(str_detect(behavior,"posture")==T,1,0),
           freeze = ifelse(str_detect(behavior,"Freeze")==T,1,0),
           sniff = ifelse(str_detect(behavior,"Sniff")==T,1,0),
           allogroom = ifelse(str_detect(behavior,"Allo")==T,1,0),
           walk = ifelse(str_detect(behavior,"Walk")==T,1,0),
           none = ifelse(str_detect(behavior,"N/A")==T,1,0)
    ) 
  return(df3)
}

#summarize each behavior data frame 
behavior_freq<-function(df){
  df %>% 
    select(-Timestamp,-behavior) %>% 
    group_by(mouseID) %>%
    summarize_all(~sum(.))
  
}

#calculate percentage of each behavior each individual showed during entire group housing period
behavior_perc<-function(df){
  df1<- df %>% 
    mutate(allbehav = select(.,2:13) %>% rowSums()) %>% 
    mutate_at(.funs = funs(perc = ./allbehav), .vars = vars(lunge:none)) %>% 
    select(mouseID,lunge_perc:none_perc)
}


#calculate percentage of three subordinate behavior each individual showed during entire group housing period
sub_behavior_perc<-function(df){
  df1<- df %>% 
    mutate(subbehav = select(.,preem_flee:freeze) %>% rowSums()) %>% 
    mutate_at(.funs = funs(perc = ./subbehav), .vars = vars(preem_flee:freeze)) %>% 
    select(mouseID,preem_flee_perc:freeze_perc)
}



contests <- function(df,a,b){
  df[(df$winner==a & df$loser==b)|(df$winner==b & df$loser==a),]
}

# Individual dominance uncertainty - Perc package 
get_uncertainty <- function(d){
  mat<-get_wl_matrix(as.data.frame(d[,2:3]))
  davids <- ds(mat)
  ord <- names(davids)[rev(order(davids))]
  mat1<-mat[ord,ord]
  #rownames(mat1)<-colnames(mat1)<-LETTERS[1:10]
  mat1c <-  as.conflictmat(mat1)
  mat1dp <- conductance(mat1c, maxLength = 3)
  out1<-valueConverter(mat1dp$p.hat) # certainty dyadic
  out2<-individualDomProb(mat1dp$p.hat)
  return(list(mat = mat1, dyads=out1, indivs=out2))
}



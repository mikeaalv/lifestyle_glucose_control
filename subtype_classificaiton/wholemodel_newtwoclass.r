# the lasso model for classify prediabtes groups
rm(list=ls())
options(warn=1)
options(stringsAsFactors=FALSE)
options(digits=15)
require(stringr)
require(magrittr)
require(ggplot2)
require(dplyr)
require(tidyr)
require(glmnet)
require(rms)
require(readxl)
require(reshape2)
comp="/Users/yuewu/"
proj_dir=paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_dataset/habit_run_02062025/")
resdir=paste0(proj_dir,"res/result_refine2/")
setwd(resdir)
# all day to avg day
# tab_allday=read.table("/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_dataset/habit_run_02062025/data/all_diet_sleep_PA_features_alldays.csv",header=TRUE,sep=",")
# cols=colnames(tab_allday)
# rmcol=cols[str_detect(string=cols,pattern=fixed(".."))]
# rmcol=c(rmcol,"date","cycle","day")
# tab_allday=tab_allday[,!cols%in%rmcol]
# tab_avgday<-tab_allday%>%group_by(subject.id)%>%dplyr::summarize(across(where(is.numeric),~ mean(.x,na.rm=TRUE),.names="{.col}"))%>%as.data.frame()
# save(tab_avgday,file="/Users/yuewu/Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_dataset/habit_run_02062025/data/all_diet_sleep_PA_features_avgdays_newpart.RData")
# feature list
response_class_list=c("a1c_t2d_status","sspg_status","DI_status","IE_status","HIR_status","adiposeIR_status")
response_cont_list=c("a1c_bl","sspg","modified_DI","ie","hepatic_IR","ffa_avg")
feat_general_list=c("age","sex","ethnicity","alcohol","exercise","education","bmi_bl")#features from meta data, demographic
# habtis related features
meal_energy=c("Meal2_Br_energy","Meal3_L_energy","Meal4_L_energy","Meal5_D_energy")#d_timing
carb_source=c("Carb_Fruits","Carb_Rice","Carb_Starchy_Vegetables","Carb_Pasta_Noodles","Carb_Non_Starchy_Vegetables","Carb_Bread_Products","Carb_Coffee","Carb_Legume_Products","Carb_Sweets","Carb_Snacks")#d_composition
eatsele_oth=c("energy","fib","chol","sat.fat","trans.fat","sodium","VitD","Ca","Fe","K","carb","pro","fat","meal.num","fasting.period","fft_t","sugar")
otherfeat=c("latency.dur","sleep.dur","waso.dur","wake_t","daily_steps","move.over.sed","active.hours")#s,s,s,s,a,a,a,D_S,d_s
habits_feat=c("subject.id",meal_energy,carb_source,eatsele_oth,otherfeat)
# new meta table
metatest_tab=read.table(paste0(proj_dir,"data/combined_meta_data_all_final3.csv"),header=TRUE,sep=",")
## consistent columns 
indnew=metatest_tab[,"cohort"]=="validation"
metatest_tab[indnew,"sex"]=ifelse(metatest_tab[indnew,"sex"]=="Female",0,1)
metatest_tab[indnew,"ethnicity"]=ifelse(metatest_tab[indnew,"ethnicity"]=="White",1,0)
grp_sets=list(training=metatest_tab[!indnew,"subject.id"],validation=c(22,74,77,79,83,85,87,90,95,100))
for(feat in c(response_cont_list,feat_general_list)){
    metatest_tab[,feat]=as.numeric(metatest_tab[,feat])
}
# features training
load(paste0(proj_dir,"data/diet_sleep_PA_features_Avgdays_old.RData"))
habit_mean=diet_sleep_PA_features_Avgdays[,habits_feat]
feattab=as.data.frame(habit_mean)
feattab[,"subject.id"]<-feattab[,"subject.id"]%>%str_remove(string=.,pattern="43883-")%>%as.numeric()
data_full=merge(feattab,metatest_tab[,c("subject.id",feat_general_list,response_cont_list,response_class_list)],by="subject.id")
data_full=data_full[data_full[,"subject.id"]%in%grp_sets[["training"]],]
# use mean to replace na for each features
xcollab=setdiff(colnames(data_full),c("subject.id",response_class_list,response_cont_list))
listscale=list()#use to scale validation 
for(featname in xcollab){
    locvec=data_full[,featname]
    locvec[is.na(locvec)]=mean(locvec,na.rm=TRUE)
    locvec=scale(locvec)
    data_full[,featname]=locvec
    listscale[[featname]]=c(attr(locvec,"scaled:center"),attr(locvec,"scaled:scale"))
}
# response training
# assign the only T2D to preDM
data_full[str_which(string=data_full[,"a1c_t2d_status"],pattern=fixed("T2D")),"a1c_t2d_status"]="preDM"
# factor classes
levelsorder=list("a1c_t2d_status"=c("Normal","preDM"),
                "sspg_status"=c("IS","IR"),
                "DI_status"=c("BC_normal","BC_dys"),#,,"BC_inter"
                "IE_status"=c("IE_normal","IE_dys"),#"IE_inter",
                "HIR_status"=c("H_IS","H_IR"),
                "adiposeIR_status"=c("adipose_IS","adipose_IR"))#"Adipocyte_Intermediate",
# put NA and factor
for(featname in response_class_list){
    locvec=data_full[,featname]
    locvec[str_which(string=locvec,pattern=fixed("Unknown"))]=NA
    # print(levels(locvec))
    data_full[,featname]=factor(locvec,levels=levelsorder[[featname]])
}
# scale continuous response
listscaleresp=list()
for(featname in response_cont_list){
    locvec=scale(data_full[,featname])
    data_full[,featname]=locvec
    listscaleresp[[featname]]=c(attr(locvec,"scaled:center"),attr(locvec,"scaled:scale"))
}
# feature validaiton 
load(paste0(proj_dir,"data/all_diet_sleep_PA_features_avgdays_newpart.RData"))
tab_avgday[,"move.over.sed"]=tab_avgday[,"move.min"]/tab_avgday[,"sedentary.min"]
habit_mean=tab_avgday[,habits_feat]
feattab=as.data.frame(habit_mean)
feattab[,"subject.id"]<-feattab[,"subject.id"]%>%str_remove(string=.,pattern="43883-")%>%as.numeric()
data_full_val=merge(feattab,metatest_tab[,c("subject.id",feat_general_list,response_cont_list,response_class_list)],by="subject.id")
data_full_val=data_full_val[data_full_val[,"subject.id"]%in%grp_sets[["validation"]],]
# use mean to replace na for each features
xcollab=setdiff(colnames(data_full_val),c("subject.id",response_class_list,response_cont_list))
for(featname in xcollab){
    locvec=data_full_val[,featname]
    locvec[is.na(locvec)]=mean(locvec,na.rm=TRUE)
    locvec=scale(locvec,center=listscale[[featname]][[1]],scale=listscale[[featname]][[2]])
    data_full_val[,featname]=locvec
}
# response validation
# assign the only T2D to preDM
data_full_val[str_which(string=data_full_val[,"a1c_t2d_status"],pattern=fixed("T2D")),"a1c_t2d_status"]="preDM"
# put NA and factor
for(featname in response_class_list){
    locvec=data_full_val[,featname]
    locvec[str_which(string=locvec,pattern=fixed("Unknown"))]=NA
    # print(levels(locvec))
    data_full_val[,featname]=factor(locvec,levels=levelsorder[[featname]])
}
# scale continuous response
for(featname in response_cont_list){
    data_full_val[,featname]=scale(data_full_val[,featname],center=listscaleresp[[featname]][[1]],scale=listscaleresp[[featname]][[2]])
}
# 
lambda_best_coll=c()
miss_class_error_coll=c()
baselineest_coll=c() #just guess the majority: the proportion of the majority class
baselineest_coll_val=c()
miss_class_error_coll_val=c()
npara_coll=c()
var_select=list()
model_coll=list()
var_select_each=list()
varexp=c()
response_class_list_clean=response_class_list
for(res_class in response_class_list_clean){
    x=as.matrix(data_full[,xcollab])
    y=data_full[,res_class]
    nafilter=!is.na(y)
    x=x[nafilter,]
    y=y[nafilter]
    # validition 
    x_val=as.matrix(data_full_val[,xcollab])
    y_val=data_full_val[,res_class]
    nafilter=!is.na(y_val)
    x_val=x_val[nafilter,]
    y_val=y_val[nafilter]
    set.seed(101)
    perfs<-NULL
    modelfamily="binomial"
    # leave one out
    cv <- cv.glmnet(x,y=y,alpha=1,family=modelfamily,type.measure="class",nfolds=length(y)) #,
    lambda.min=cv$lambda.min
    minind=which(cv$lambda==lambda.min)
    if(cv$nzero[minind]==0){
        ind_equal=sort(which(cv$cvm==min(cv$cvm)))
        minind=ind_equal[length(ind_equal)]
        lambda.min=cv$lambda[minind]
    }
    # 
    lambda_best_coll=c(lambda_best_coll,lambda.min)
    npara_coll=c(npara_coll,cv$nzero[minind])
    coefmat=coef(cv,s=lambda.min)
    # selected features
    selectedvar=xcollab[coefmat@i[-1]]
    constant=rep(1,times=dim(x)[1])
    xsele=x[,selectedvar,drop=FALSE]#add constant col as glmnet doesn't work with one feature
    if(dim(xsele)[2]==1){
        xsele=cbind(xsele,constant)
    }
    # calculate performance by just rerun with the selected feature set
    cv <- cv.glmnet(xsele,y=y,alpha=1,family=modelfamily,type.measure="class",nfolds=length(y),lambda=c(0,1),maxit=1000000) #,
    zerolambdain=which(cv$lambda==0)
    perfs<-cv$cvm
    miss_class_error_coll=c(miss_class_error_coll,perfs[zerolambdain])
    varexp=c(varexp,cv$glmnet.fit$dev.ratio[zerolambdain])
    # baseline perf
    ytab=table(y)
    baselineest_coll=c(baselineest_coll,1-max(ytab)/sum(ytab))
    # validation baseline perf
    ytab=table(y_val)
    baselineest_coll_val=c(baselineest_coll_val,1-max(ytab)/sum(ytab))
    # predictions and for saving of final coef estimation
    mfit<-glmnet(xsele,y,family=modelfamily,alpha=1,lambda=0,maxit=1000000)
    # validation
    y_pred_val=predict(mfit,newx=x_val[,selectedvar,drop=FALSE],type="class")
    miss_class_error_coll_val=c(miss_class_error_coll_val,sum(y_val!=y_pred_val)/length(y_pred_val))
    # only works for binomial (2 classes)
    # logistic regression for EACH feature TO calcuate p value 
    temp_sin_tab=c()
    for(feat in selectedvar){
        datadf=data.frame(x=xsele[,feat],class=y)
        glm_fit=glm(class~x,data=datadf,family=binomial)
        glm_fit_sum=summary(glm_fit)
        temp_sin_tab=rbind(temp_sin_tab,glm_fit_sum$coefficients["x",c(1,2,4)])
    }
    colnames(temp_sin_tab)=c("val","ste","pvalue")
    rownames(temp_sin_tab)=selectedvar
    model_coll[[res_class]]=mfit
    # print(1-sum(predict(mfit,xtest[,selectedvar],type="class")==ytest)/5)
    var_select[[res_class]]=coef(mfit)
    var_select_each[[res_class]]=temp_sin_tab
}
datsumary=data.frame(classes=response_class_list_clean,lambda=lambda_best_coll,nzero=npara_coll,var_exp=varexp,Miss_class_err=miss_class_error_coll,Miss_class_err_val=miss_class_error_coll_val,baseline=baselineest_coll,baseline_validation=baselineest_coll_val)
var_select_class=var_select
save(datsumary,var_select_class,var_select_each,model_coll,file=paste0(resdir,"lasso_fit_class.RData"))
write.table(format(datsumary,digits=2),file="perf_class.tsv",sep="\t",row.names=FALSE)
# visualize targetted features
## t.test for all selected features for the classification
ttest_tab=c()
for(classname in names(var_select)){
    subsetsele=var_select[[classname]]@Dimnames[[1]][-1]
    subsetsele=subsetsele[subsetsele!="constant"]
    x=as.matrix(data_full[,subsetsele,drop=FALSE])
    y=data_full[,classname]
    nafilter=!is.na(y)
    x=x[nafilter,,drop=FALSE]
    y=y[nafilter]
    for(featele in subsetsele){
        aovtest=aov(val~type,data.frame(val=x[,featele],type=y))
        ttest_tab=rbind(ttest_tab,data.frame(typename=classname,feature=featele,pval=summary(aovtest)[[1]][["Pr(>F)"]][1])) 
    }
}
# plot box plot 
ttest_tab_sig=ttest_tab[ttest_tab[,"pval"]<0.05,]
for(rowi in seq(nrow(ttest_tab_sig))){
    subtab_p=ttest_tab_sig[rowi,]
    data_full_sub=data_full[!is.na(data_full[,subtab_p[,"typename"]]),]
    p<-ggplot(data_full_sub,aes_string(x=subtab_p[,"typename"],y=subtab_p[,"feature"])) + geom_boxplot()
    ggsave(paste0(resdir,subtab_p[,"feature"],"_diff_btw_",subtab_p[,"typename"],".pdf"))
}
for(rowi in seq(nrow(ttest_tab_sig))){
    subtab_p=ttest_tab_sig[rowi,]
    data_full_sub=data_full_val[!is.na(data_full_val[,subtab_p[,"typename"]]),]
    p<-ggplot(data_full_sub,aes_string(x=subtab_p[,"typename"],y=subtab_p[,"feature"])) + geom_boxplot()
    ggsave(paste0(resdir,subtab_p[,"feature"],"_diff_btw_",subtab_p[,"typename"],"valid.pdf"))
}
# plotting features 
color_group=list("diet"="#A6B438","activity"="#F16667","sleep"="#BD72AF","demographic"="#C0C0C0")
color_feature_vec=vector(length=length(xcollab))
names(color_feature_vec)=xcollab
color_feature_vec[feat_general_list]=color_group[["demographic"]]
color_feature_vec[c(meal_energy,carb_source,eatsele_oth)]=color_group[["diet"]]
color_feature_vec[c("latency.dur","sleep.dur","waso.dur","wake_t","last.to.bed","wake.to.first")]=color_group[["sleep"]]# set all diet_sleep features sleep features
color_feature_vec[c("daily_steps","move.over.sed","active.hours")]=color_group[["activity"]]
color_feature_vec[c("exercise")]=color_group[["activity"]]
varlabtab=read_excel(paste0(comp,"Library/CloudStorage/Box-Box/Yue Wu's Files/cgm_dataset/heyjun_data_run/Variables_Labeling.xlsx"),col_names=FALSE)
label_vec=varlabtab[[2]]
names(label_vec)=varlabtab[[1]]
for(res_class in response_class_list_clean){
    coefmat=var_select_class[[res_class]]
    names=coefmat@Dimnames[[1]]
    vals=coefmat@x
    coef_df=data.frame(coef=names[-1],val=vals[-1])
    reordind=order(coef_df[,"val"])
    coef_df=coef_df[reordind,]
    featname=coef_df[,"coef"]
    coef_df[,"coef"]=ifelse(featname%in%names(label_vec),label_vec[featname],featname)
    coef_df[,"coef"]=factor(coef_df[,"coef"],levels=coef_df[,"coef"])
    lablval=paste0("error rates ",format(datsumary[datsumary[,"classes"]==res_class,"Miss_class_err"],digits=2,nsmall=2))
    p<-ggplot(coef_df,aes(x=coef,y=val))+geom_bar(stat="identity",fill=color_feature_vec[featname])+coord_flip()+labs(x="",y="Coefficients")+ggtitle(res_class)+theme_bw()+theme(panel.background=element_blank(),panel.grid=element_blank())+annotate("text",x= Inf,y=0,label=lablval,vjust=1,hjust=1)
    ggsave(paste0(resdir,"coef_",res_class,".pdf"),width=5,height=5)
}
# each with separate logistic regresssion 
for(res_class in response_class_list_clean){
    coefmat=var_select_each[[res_class]]
    coef_df=as.data.frame(coefmat)
    reordind=order(coef_df[,"val"])
    coef_df=coef_df[reordind,]
    featname=rownames(coef_df)
    coef_df[,"coef"]=ifelse(featname%in%names(label_vec),label_vec[featname],featname)
    coef_df[,"coef"]=factor(coef_df[,"coef"],levels=coef_df[,"coef"])
    lablval=paste0("error rates ",format(datsumary[datsumary[,"classes"]==res_class,"Miss_class_err"],digits=2,nsmall=2))
    p<-ggplot(coef_df,aes(x=coef,y=val))+geom_bar(stat="identity",fill=color_feature_vec[featname]) + geom_errorbar(aes(ymin=val-ste*2,ymax=val+ste*2))+coord_flip()+labs(x="",y="Coefficients")+ggtitle(res_class)+theme_bw()+theme(panel.background=element_blank(),panel.grid=element_blank())+annotate("text",x= Inf,y=0,label=lablval,vjust=1,hjust=1)
    ggsave(paste0(resdir,"coef_",res_class,"_sep_logisticreg_errorbar.pdf"),width=5,height=5)
    p<-ggplot(coef_df,aes(x=coef,y=val))+geom_bar(stat="identity",fill=color_feature_vec[featname])+geom_text(aes(label=ifelse(pvalue<0.05,"*","")),position=position_dodge(width =.9),size=20/.pt)+coord_flip()+labs(x="",y="Coefficients")+ggtitle(res_class)+theme_bw()+theme(panel.background=element_blank(),panel.grid=element_blank())+annotate("text",x= Inf,y=0,label=lablval,vjust=1,hjust=1)
    ggsave(paste0(resdir,"coef_",res_class,"_sep_logisticreg_pval.pdf"),width=5,height=5)
}
# plot error rates
datsumary2=datsumary
colnames(datsumary2)[5]="Lasso"
datsumary_l=melt(datsumary2[,c("classes","Lasso","baseline")],id=c("classes"))
p<-ggplot(datsumary_l,aes(x=classes,y=value,fill=variable))+geom_bar(stat="identity",position=position_dodge())+scale_y_continuous(limits=c(-0.1,1))+labs(x="subtypes",y="miss classification error")+theme_bw()+theme(panel.background=element_blank(),panel.grid=element_blank())
ggsave(paste0(resdir,"coef_bar.pdf"),width=5,height=5)
# continuous variable as response 
lambda_best_coll=c()
corr_coll=c()
npara_coll=c()
var_select=list()
model_coll=list()
varexp=c()
response_cont_list_clean=response_cont_list[c(-3,-4,-6)]#no modified_DI as the problem from small sample size
for(res_class in response_cont_list_clean){
    x=as.matrix(data_full[,xcollab])
    y=data_full[,res_class]
    nafilter=!is.na(y)
    x=x[nafilter,]
    y=y[nafilter]
    set.seed(101)
    perfs<-NULL
    modelfamily="gaussian"
    # leave one out
    cv <- cv.glmnet(x,y=y,alpha=1,family=modelfamily,type.measure="mse",nfolds=length(y)) #,
    lambda.min=cv$lambda.min
    minind=which(cv$lambda==lambda.min)
    if(cv$nzero[minind]==0){
        ind_equal=sort(which(cv$cvm==min(cv$cvm)))
        minind=ind_equal[length(ind_equal)]
        lambda.min=cv$lambda[minind]
    }
    # 
    lambda_best_coll=c(lambda_best_coll,lambda.min)
    npara_coll=c(npara_coll,cv$nzero[minind])
    coefmat=coef(cv,s=lambda.min)
    selectedvar=xcollab[coefmat@i[-1]]
    xsele=x[,selectedvar,drop=FALSE]
    y_coll=c()
    for(i in seq(length(y))){
        mfit1<-glmnet(xsele[-i,],y[-i],family=modelfamily,alpha=1,lambda=0)
        y_coll=c(y_coll,predict(mfit1,newx=xsele[c(1,i),])[2])
    }
    corr_coll=c(corr_coll,cor(y_coll,y))
    # 
    cv <- cv.glmnet(xsele,y=y,alpha=1,family=modelfamily,type.measure="mse",nfolds=length(y),lambda=c(0,1),keep=TRUE) #,
    zerolambdain=which(cv$lambda==0)
    varexp=c(varexp,cv$glmnet.fit$dev.ratio[zerolambdain])
    # 
    mfit<-glmnet(xsele,y,family=modelfamily,alpha=1,lambda=0)
    model_coll[[res_class]]=mfit
    # print(1-sum(predict(mfit,xtest[,selectedvar],type="class")==ytest)/5)
    var_select[[res_class]]=coef(mfit)
}
datsumary=data.frame(classes=response_cont_list_clean,lambda=lambda_best_coll,nzero=npara_coll,var_exp=varexp,corr=corr_coll)
var_select_cont=var_select
save(datsumary,var_select_cont,model_coll,file=paste0(resdir,"lasso_fit_cont.RData"))
# plotting features 
for(res_class in response_cont_list_clean){
    coefmat=var_select_cont[[res_class]]
    names=coefmat@Dimnames[[1]]
    vals=coefmat@x
    coef_df=data.frame(coef=names[-1],val=vals[-1])
    reordind=order(coef_df[,"val"])
    coef_df=coef_df[reordind,]
    featname=coef_df[,"coef"]
    coef_df[,"coef"]=factor(coef_df[,"coef"],levels=coef_df[,"coef"])
    lablval=paste0("correlation ",format(datsumary[datsumary[,"classes"]==res_class,"corr"],digits=2,nsmall=2))
    p<-ggplot(coef_df,aes(x=coef,y=val)) + geom_bar(stat="identity",fill=color_feature_vec[featname]) + coord_flip() + labs(x="",y="Coefficients") + ggtitle(res_class) + theme_bw() + theme(panel.background=element_blank(),panel.grid=element_blank()) + annotate("text",x= Inf,y= 0,label=lablval,vjust=1,hjust=1)
    ggsave(paste0(resdir,"coef_",res_class,".reg.pdf"))
}

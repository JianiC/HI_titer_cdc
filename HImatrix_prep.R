## preapare HI titer data
## filter from influenza subtype and collection year
source("setup.R")

HI_titer_raw<-read.csv("uga_titers_2010_onward_H3.csv")
HI_titer_raw %>% 
  filter(ag_ha_type=="H3")%>%
  mutate(ag_collection_date2=as_date(ag_collection_date,format="%m/%d/%y"))%>%
  filter(ag_collection_date2>="2010-9-1" & ag_collection_date2 <= "2020-5-1")-> uga_titer_2020_H3

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  group_by(ag_strain_name)%>%
  summarise(sr_count=n())->ag_srcount
  
uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  group_by(sr_strain_name)%>%
  summarise(count=n())->sr_strain_count

#######################################
## combine sr_strain_name

#write.csv(sr_strain_count,"sr_strain_name.csv",row.names=FALSE)
#write.csv(ag_srcount,"ag_strain_name.csv",row.names=FALSE)
sr_strain<-read.csv("sr_strain_name.csv")
ag_strain<-read.csv("ag_strain_name.csv")
## translate sr_strain_name
sr_strain%>%
  select(-count)%>%
  merge(uga_titer_2020_H3,by=c("sr_strain_name","sr_strain_name"))%>%
  merge(ag_strain,by=c("ag_strain_name","ag_strain_name"))->uga_titer_2020_H3





## all sequence matrix
uga_titer_2020_H3%>%
  mutate(ag_strain_name2=gsub(" ","_",ag_strain_name2))%>%
  distinct(ag_epi_isolate_id,sr_epi_isolate_id,.keep_all= TRUE)%>%
  select(ag_epi_isolate_id,sr_epi_isolate_id,titer_value)%>%
  spread(key = "sr_epi_isolate_id", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->test_matrix

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_epi_isolate_id,.keep_all= TRUE)%>%
  select(ag_epi_isolate_id,sr_strain_name,titer_value)%>%
  group_by(ag_strain_name)%>%
  summarise(count=n())%>%
  ggplot(aes(x=count))+
  geom_histogram()+
  scale_x_continuous(limits = c(0,50,2))

## get strain name 


uga_titer_2020_H3%>%
  distinct(ag_epi_isolate_id)->uga_titer_strain_2020_H3



write.table(uga_titer_2020_H3$ag_strain_taxa,"HI_titer_strain.txt",sep="\t",row.names=FALSE)


write.csv(test_matrix,"2010_2020_raw_titer.csv",row.names = FALSE)

################################################################
## filter the sequence with high accuracy

uga_titer_2020_H3%>%
  mutate(ag_strain_name2=gsub(" ","_",ag_strain_name2))%>%
  merge(ag_srcount,by=c("ag_strain_name","ag_strain_name"))%>%
  mutate(ag_collection_year=year(ag_collection_date2))%>%
  distinct(ag_epi_isolate_id,.keep_all= TRUE)%>%
  #arrange(desc(titer_value,sr_count))%>%
  group_by(ag_collection_year)%>%
  slice_sample(n=30)->uga_titer_strain_filteryear

uga_titer_2020_H3%>%
  mutate(ag_strain_name2=gsub(" ","_",ag_strain_name2))%>%
  distinct(ag_epi_isolate_id,sr_epi_isolate_id,.keep_all= TRUE)%>%
  filter(ag_epi_isolate_id %in% uga_titer_strain_filteryear$ag_epi_isolate_id)%>%
  select(ag_epi_isolate_id,sr_epi_isolate_id,titer_value)%>%
  spread(key = "sr_epi_isolate_id", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->HI_titer_filter_matrix

## get sequence to test

names(HI_titer_filter_matrix)[1] <- ""
write.csv(HI_titer_filter_matrix,"2010_2020_raw_titer_filter.csv",row.names = FALSE)
write.table(uga_titer_strain_filteryear$ag_strain_taxa,"HI_titer_yearfilter_strain.txt",sep="\t",row.names=FALSE)

## from test 1

rm_sr<-c('A/Alaska/01/2021',
         'A/Arizona/45/2018',
         'A/Brisbane/190/2017',
         'A/California/213/2016',
         'A/California/55/2020',
         'A/Cambodia/E0826360/2020',
         'A/Darwin/6/2021',
         'A/Darwin/9/2021',
         'A/Delaware/01/2021',
         'A/Greece/4/2017',
         'A/Honduras/0188/2017',
         'A/Hong Kong/681/2018',
         'A/Kentucky/04/2012',
         'A/Maryland/02/2021',
         'A/Maryland/07/2018',
         'A/Michigan/173/2020',
         'A/Netherlands/10260/2018',
         'A/New York/46/2018',
         'A/Norway/3275/2018',
         'A/Singapore/GP2646/2016',
         'A/Sydney/22/2018',
         'A/Sydney/53/2019',
         'A/Texas/68/2017',
         'A/Togo/771/2020',
         'A/Uruguay/43/2017',
         'A/Washington/04/2017',
         'A/Wisconsin/02/2021',
         'A/Alaska/240/2015',
         'A/Hong Kong/2286/2017',
         'A/Hong Kong/5738/2014',
         'A/Maryland/3/2014')
## read the strain name that contains the sequence data
ag_seq<-scan("filter_ag_strain.txt", what="", sep="\n")




rm_ag<-c(
         'A/Virginia/03/2020')


## further filter the dataframe to test

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name2,.keep_all= TRUE)%>%
  filter(ag_strain_name %in% uga_titer_strain_filteryear$ag_strain_name)%>%
  #filter(ag_strain_name %in% ag_seq)%>%
  filter(!(sr_strain_name %in% rm_sr))%>%
  filter(!(ag_strain_name %in% rm_ag))%>%
  select(ag_strain_name,sr_strain_name2,titer_value)%>%
  spread(key = "sr_strain_name2", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->HI_titer_filter_matrix2

names(HI_titer_filter_matrix2)[1] <- ""
write.csv(HI_titer_filter_matrix2,"2010_2020_raw_titer_filter2.csv",row.names = FALSE)
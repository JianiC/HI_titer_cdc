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

## all sequence matrix
uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  spread(key = "sr_strain_name", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->test_matrix

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  group_by(ag_strain_name)%>%
  summarise(count=n())%>%
  ggplot(aes(x=count))+
  geom_histogram()+
  scale_x_continuous(limits = c(0,50,2))

## get strain name 

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,.keep_all= TRUE)%>%
  mutate(ag_strain_taxa=paste0(ag_epi_isolate_id,"|",ag_strain_name))->uga_titer_strain_2020_H3

write.table(uga_titer_2020_H3$ag_strain_taxa,"HI_titer_strain.txt",sep="\t",row.names=FALSE)

write.csv(test_matrix,"2010_2020_raw_titer.csv",row.names = FALSE)

################################################################
## filter the sequence with high accuracy

uga_titer_strain_2020_H3%>%
  merge(ag_srcount,by=c("ag_strain_name","ag_strain_name"))%>%
  mutate(ag_collection_year=year(ag_collection_date2))%>%
  arrange(desc(titer_value,sr_count))%>%
  group_by(ag_collection_year)%>%
  slice(1:40)->uga_titer_strain_filteryear

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  filter(ag_strain_name %in% uga_titer_strain_filteryear$ag_strain_name)%>%
  
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  spread(key = "sr_strain_name", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->HI_titer_filter_matrix

## get sequence to test


write.csv(HI_titer_filter_matrix,"2010_2020_raw_titer_filter.csv",row.names = FALSE)
write.table(uga_titer_strain_filteryear$ag_strain_taxa,"HI_titer_yearfilter_strain.txt",sep="\t",row.names=FALSE)

## from test 1

rm_sr<-c('A/Brisbane/01/2018',
              'A/Brisbane/1/2018 X-311A',
              'A/Brisbane/29/2017',
              'A/Greece/4/2017',
              'A/Hawaii/22/2012 X-225A',
              'A/Hong Kong/5738/2014',
              'A/Hong Kong/681/2018',
              'A/Idaho/33/2016 X-293',
              'A/Kansas/14/2017 CBER-22B',
              'A/New Jersey/26/2014',
              'A/New York/46/2018',
              'A/Newcastle/82/2018',
              'A/North Carolina/13/2014 X-255A',
              'A/Norway/3275/2018',
              'A/Norway/3806/2016 CBER-08 C1.4',
              'A/Singapore/GP2646/2016',
              'A/Switzerland/9715293/2013 FA2337',
              'A/Switzerland/9715293/2013 FB2119',
              'A/Sydney/53/2019',
              'A/Uruguay/43/2017',
              'A/Victoria/361/2011 X-217',
              'A/Victoria/361/2011 X-217A',
              'A/Victoria/505/2013',
              'A/Washington/04/2017',
              'A/Washington/16/2017 X-301',
              'A/Washington/16/2017 X-301A',
              'A/Wisconsin/95/2018',
              'A/Abu Dhabi/240/2018 X-325',
              'A/Alaska/232/2015 X-289',
              'A/Alaska/240/2015',
              'A/American Samoa/4786/2013',
              'A/Arizona/45/2018',
              'A/Delaware/32/2016',
              'A/Delaware/41/2019',
              'A/Hawaii/22/2012 X-225',
              'A/Hong Kong/4801/2014 CDC-LV15A',
              'A/Kansas/14/2017 X-327A',
              'A/Maryland/3/2014',
              'A/Montana/18/2019',
              'A/Nebraska/19/2015',
              'A/Netherlands/10260/2018',
              'A/Ohio/02/2012 X-221',
              'A/Palau/6759/2014 X-245',
              'A/Palau/6759/2014 X-245A',
              'A/SOUTH AUSTRALIA/2/2019-CDC-LV26A',
              'A/Switzerland/8060/2017',
              'A/Texas/50/2012 X-223',
              'A/Texas/68/2017',
              'A/Wisconsin/19/2017')

rm_ag<-c('A/Alaska/01/2011',
         'A/Alaska/38/2018',
         'A/Argentina/139/2014',
         'A/Brazil/7920/2012',
         'A/California/16/2011',
         'A/Florida/39/2014',
         'A/Idaho/20/2012',
         'A/Laos/348/2013',
         'A/Michigan/482/2019',
         'A/Virginia/03/2020')


## further filter the dataframe to test

uga_titer_2020_H3%>%
  mutate(ag_strain_name=gsub(" ","_",ag_strain_name))%>%
  distinct(ag_strain_name,sr_strain_name,.keep_all= TRUE)%>%
  filter(ag_strain_name %in% uga_titer_strain_filteryear$ag_strain_name)%>%
  filter(!(ag_strain_name %in% rm_ag))%>%
  filter(!(sr_strain_name %in% rm_sr))%>%
  select(ag_strain_name,sr_strain_name,titer_value)%>%
  spread(key = "sr_strain_name", value = "titer_value")%>%
  mutate_all(~replace(., is.na(.), "*"))->HI_titer_filter_matrix2


write.csv(HI_titer_filter_matrix2,"2010_2020_raw_titer_filter2.csv",row.names = FALSE)
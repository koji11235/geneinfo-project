import urllib.request
from bs4 import BeautifulSoup
import pandas as pd
import time
import lxml.html
import copy
import json
import sqlite3
import xmltodict
import ssl
import numpy as np

def fetch_gwascatalog(chromosome,GenomeLocation):
    GWASCatalog_link=f"https://www.ebi.ac.uk/gwas/rest/api/snpLocation/{chromosome}:{GenomeLocation}"
    print(GWASCatalog_link)
    data = urllib.request.urlopen(GWASCatalog_link)
    #jsonの整形出力
    data = json.load(data)
    #print(json.dumps(data["_embedded"]["singleNucleotidePolymorphisms"][0], indent=2))
    #print(data)
    rsid_list=[]
    #アクセスしたChr:GenomeLocationに情報があれば、返ってくるファイルには"_embedded"が入ってるみたいだからそれを確認してる
    #3つくらいでしか試してないから全対応してるか分からん
    if "_embedded" in data.keys():
        for content in data["_embedded"]["singleNucleotidePolymorphisms"]:
            #print(content)
            #print(content['rsId'])
            rsid_list.append(content['rsId'])
            #print("rsid: ",content["rsid"])

        GWAS_df=pd.DataFrame(columns=["rsid","PubmedID","Study_Accession","First_Author","Publication_Date","Journal","Title","Reported_Trait","Discovery_Sample_Number_And_Ancestry","Traits"])
        for rsid in rsid_list:
            PubmedID=""
            First_Author=""
            Study_Accession=""
            Publication_Date=""
            Journal=""
            Title=""
            Reported_Trait=""
            Traits=[]#Traitsだけ複数の可能性があるからリスト使ってる
            Discovery_Sample_Number_And_Ancestry=""
            Replication_sample_number_and_ancestry=""
            #Association_count=""
            #Summary_statistics=""
            GWAS_Data_Dict={"rsid":rsid,"PubmedID":PubmedID,"Study_Accession":Study_Accession,"First_Author":First_Author,"Publication_Date":Publication_Date,"Journal":Journal,"Title":Title,"Reported_Trait":Reported_Trait,"Discovery_Sample_Number_And_Ancestry":Discovery_Sample_Number_And_Ancestry,"Traits":Traits}

            RSID_study_link="https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/studies".format(rsid=rsid)
            data = urllib.request.urlopen(RSID_study_link)
            data = json.load(data)
            print(RSID_study_link)
            #print(data)
            if not data["_embedded"]["studies"]==[]:
                Study_Accession=data["_embedded"]["studies"][0]["accessionId"]
                First_Author=data["_embedded"]["studies"][0]["publicationInfo"]["author"]["fullname"]
                PubmedID=data["_embedded"]["studies"][0]["publicationInfo"]["pubmedId"]
                Publication_Date=data["_embedded"]["studies"][0]["publicationInfo"]["publicationDate"]
                Journal=data["_embedded"]["studies"][0]["publicationInfo"]["publication"]
                Title=data["_embedded"]["studies"][0]["publicationInfo"]["title"]
                Reported_Trait=data["_embedded"]["studies"][0]["diseaseTrait"]["trait"]
                Discovery_Sample_Number_And_Ancestry=str(data["_embedded"]["studies"][0]["ancestries"][0]["numberOfIndividuals"])+" "+data["_embedded"]["studies"][0]["ancestries"][0]["ancestralGroups"][0]["ancestralGroup"]
                """
                print("Study_Accession",Study_Accession)
                print("First_Author",First_Author)
                print("PubmedID",PubmedID)
                print("Publication_Date",Publication_Date)
                print("Journal",Journal)
                print("Title",Title)
                print("Reported_Trait",Reported_Trait)
                print("Discovery_Sample_Number_And_Ancestry",Discovery_Sample_Number_And_Ancestry)
                """
                Accession_Traits_link="https://www.ebi.ac.uk/gwas/rest/api/studies/{Study_Accession}/efoTraits".format(Study_Accession=Study_Accession)
                data = urllib.request.urlopen(Accession_Traits_link)
                data = json.load(data)
                #print(data)
                for trait_data in data["_embedded"]["efoTraits"]:
                    Traits.append(trait_data["trait"])
                #print(Traits)
            GWAS_Data_Dict={"rsid":rsid,"PubmedID":PubmedID,"Study_Accession":Study_Accession,"First_Author":First_Author,"Publication_Date":Publication_Date,"Journal":Journal,"Title":Title,"Reported_Trait":Reported_Trait,"Discovery_Sample_Number_And_Ancestry":Discovery_Sample_Number_And_Ancestry,"Traits":Traits}
            GWAS_df=GWAS_df.append(GWAS_Data_Dict,ignore_index=True)
            #print(GWAS_df.values.tolist())
            chosen_column=["rsid","PubmedID","Study_Accession","Reported_Trait"]
        return GWAS_df[chosen_column]
    else:
        result=pd.DataFrame(columns=["rsid","PubmedID","Study_Accession","First_Author","Publication_Date","Journal","Title","Reported_Trait","Discovery_Sample_Number_And_Ancestry","Traits"])
        chosen_column=["rsid","PubmedID","Study_Accession","Reported_Trait"]
        return result[chosen_column]
"""
chromosome='9'
GenomeLocation='136881912-136926607'
GWAS_df=fetch_gwascatalog(chromosome,GenomeLocation)
print(GWAS_df)#Empty DataFrameが返ってくる（TRAF2にGWAS情報はない）

chromosome='5'
GenomeLocation='157142933-157255191'
GWAS_df=fetch_gwascatalog(chromosome,GenomeLocation)
print(GWAS_df)
"""
chromosome='5'
GenomeLocation='157142933-157255191'
GWAS_df=fetch_gwascatalog(chromosome,GenomeLocation)
print(GWAS_df)

#CSF1R 1436
chromosome='5'
GenomeLocation='150053291-150113372 '
GWAS_df=fetch_gwascatalog(chromosome,GenomeLocation)
print(GWAS_df)
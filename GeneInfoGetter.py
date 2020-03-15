#インターネットから情報取得するための関数のコレクション
#変数名の命名方法がlower_separatedとCapitalどちらもあって統一されてない
#TODO:変数名の命名方法をlower_separatedに統一する

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

#geneidentifier(http://biogps.org/ext/symatlasbar/?geneid={})から
#GeneName,Uniprot_ID, HomogeneID, Homologene_linkを抽出して辞書型で返す
#TODO:まだ三つくらいのGeneIDでしかテストしてない
def fetch_geneidentifier(geneid):
    GeneIdentifier_link="http://biogps.org/ext/symatlasbar/?geneid={}&hidespecies=1".format(geneid)
    data = urllib.request.urlopen(GeneIdentifier_link)
    html = data.read()
    soup=BeautifulSoup(html, "lxml")
    Symbol=soup.find_all("tr")[0].find(class_="annval").text
    Description=soup.find_all("tr")[1].find(class_="annval").text
    Aliases=soup.find_all("tr")[3].find(class_="annval").text.strip()
    
    
    #XXX:このサイクルでUniprot_ID、Homologene_link、GenomeLocationの三つを一気に取得して
    #GenomeLocationが取れたらBreakしてるけどfor文を三つに分けて一つずつ取得した方が確実だと思う。
    #現状これで動いてるけど
    Uniprot_ID=""
    OMIM_ID=""
    OMIM_link=""
    HomogeneID=""
    Homologene_link=""
    chromosome=""
    GenomeLocation=""

    for a in soup.find_all("a"):
        if "http://www.uniprot.org/uniprot/" in a.get("href"):
            Uniprot_ID=a.text
        
        if "http://omim.org/entry/" in a.get("href"):
            OMIM_ID=a.text
            OMIM_link=a.get("href")

        if "http://www.ncbi.nlm.nih.gov/homologene/" in a.get("href"):
            HomogeneID=a.text
            Homologene_link=a.get("href")
        
        if "http://genome.ucsc.edu/cgi-bin/hgTracks?db=" in a.get("href"):
            GenomeLocation=a.text
            chromosome=GenomeLocation.split(":")[0].strip("chr")
            GenomeLocation=GenomeLocation.split(":")[1]
            break
    #print("Uniprot_ID:",Uniprot_ID)
    #print("HomogeneID",HomogeneID)
    #print(Homologene_link)
    result={"Symbol":Symbol,\
            "Description":Description,\
            "Aliases":Aliases,\
            "Uniprot_ID":Uniprot_ID,\
            "OMIM_ID":OMIM_ID,\
            "OMIM_link":OMIM_link,\
            "HomogeneID":HomogeneID,\
            "Homologene_link":Homologene_link,\
            "chromosome":chromosome,\
            "GenomeLocation":GenomeLocation}
    return result

"""
geneid=1#3702
result=fetch_geneidentifier(geneid)
print(result)
"""

#関数名はfetch_sequenceとかのがわかりやすいか
def fetch_uniprot(uniprotid):
    Uniprot_link=f"https://www.uniprot.org/uniprot/{uniprotid}.xml"
    data = urllib.request.urlopen(Uniprot_link)
    xml=data.read()
    soup=BeautifulSoup(xml, "xml")
    sequence=soup.find_all("sequence")[-1].text.strip()
    #print("sequence: ",sequence)
    if soup.find_all("comment")!=[]:
        function=soup.find_all("comment")[0].text.strip()
    else:
        function=""
    result={"sequence":sequence,"function": function}
    return result

"""
uniprotid="Q12933"
uniprotid="Q08881"
sequence=fetch_uniprot(uniprotid)
print(sequence)
"""

def fetch_allignment(homologeneid):
    Homologene_alignment_link=f"https://www.ncbi.nlm.nih.gov/homologene/{homologeneid}?report=alignmentscores"
    data = urllib.request.urlopen(Homologene_alignment_link)
    xml=data.read()
    soup=BeautifulSoup(xml, "xml")
    alignment_table=soup.find("table",class_="alignment_table").find_all("tr")
    #print(alignment_table)
    for idx, content in enumerate(alignment_table):
        #print(content)
        #print(content.find_all(class_="top_species"))
        if content.find_all(class_="top_species")!=[]:
            if content.find(class_="top_species").text=="H.sapiens":
                H_sapiens_begin=copy.deepcopy(idx)
                break
    #print(H_sapiens_begin)

    for idx, content in enumerate(alignment_table[H_sapiens_begin+1:]):
        #print(content)
        #print(content.find_all(class_="top_species"))
        if content.find_all(class_="top_species")!=[]:
            H_sapiens_end=copy.deepcopy(idx)+H_sapiens_begin+1
            break
    #print(H_sapiens_end)
    #print(alignment_table[H_sapiens_begin:H_sapiens_end])
    Identity_Protein_Mouse=""
    Identity_DNA_Mouse=""
    Identity_Protein_Rat=""
    Identity_DNA_Rat=""
    for content in alignment_table[H_sapiens_begin+1:H_sapiens_end]:
        if content.find(class_="second_species").text=="vs. M.musculus":
            Identity_Protein_Mouse=content.find(class_="alignment_table_col3").text
            Identity_DNA_Mouse=content.find(class_="alignment_table_col4").text
        if content.find(class_="second_species").text=="vs. R.norvegicus":
            Identity_Protein_Rat=content.find(class_="alignment_table_col3").text
            Identity_DNA_Rat=content.find(class_="alignment_table_col4").text
    #print("Identity_Protein_Mouse;",Identity_Protein_Mouse)
    #print("Identity_DNA_Mouse;",Identity_DNA_Mouse)
    #print("Identity_Protein_Rat;",Identity_Protein_Rat)
    #print("Identity_DNA_Rat;",Identity_DNA_Rat)
    
    #TODO:各々の数字がstrで返される。intにした方がいいと思う。
    result={"Identity_Protein_Mouse":Identity_Protein_Mouse,\
            "Identity_DNA_Mouse":Identity_DNA_Mouse,\
            "Identity_Protein_Rat":Identity_Protein_Rat,\
            "Identity_DNA_Rat":Identity_DNA_Rat}
    return result
"""
homologeneid=4138
fa_result=fetch_allignment(homologeneid)
print(fa_result)
"""

def fetch_gwascatalog(chromosome,GenomeLocation):
    context = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
    GWASCatalog_link=f"https://www.ebi.ac.uk/gwas/rest/api/snpLocation/{chromosome}:{GenomeLocation}"
    try:
        servererror=False
        data = urllib.request.urlopen(GWASCatalog_link,context=context)
        #jsonの整形出力
        data = json.load(data)
        #print(json.dumps(data["_embedded"]["singleNucleotidePolymorphisms"][0], indent=2))
        #print(data)
        rsid_list=[]
        #アクセスしたChr:GenomeLocationに情報があれば、返ってくるファイルには"_embedded"が入ってるみたいだからそれを確認してる
        #3つくらいでしか試してないから全対応してるか分からん
        if "_embedded" in data.keys():
            for content in data["_embedded"]["singleNucleotidePolymorphisms"]:
                rsid_list.append(content['rsId'])
                
            GWAS_df=pd.DataFrame(columns=["rsid","PubmedID","Study_Accession","First_Author","Publication_Date","Journal","Title","Reported_Trait","Discovery_Sample_Number_And_Ancestry","Traits"])
            for rsid in rsid_list[:10]:
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
                GWAS_Data_Dict={"rsid":rsid,"PubmedID":PubmedID,"Study_Accession":Study_Accession,"First_Author":First_Author,"Publication_Date":Publication_Date,"Journal":Journal,"Title":Title,"Reported_Trait":Reported_Trait,"Discovery_Sample_Number_And_Ancestry":Discovery_Sample_Number_And_Ancestry,"Traits":Traits}

                RSID_study_link="https://www.ebi.ac.uk/gwas/rest/api/singleNucleotidePolymorphisms/{rsid}/studies".format(rsid=rsid)
                data = urllib.request.urlopen(RSID_study_link)
                data = json.load(data)
                print(data)
                print(data["_embedded"]["studies"])
                print(data["_embedded"]["studies"]==[])
                if data["_embedded"]["studies"]==[]:
                    continue
                Study_Accession=data["_embedded"]["studies"][0]["accessionId"]
                First_Author=data["_embedded"]["studies"][0]["publicationInfo"]["author"]["fullname"]
                PubmedID=data["_embedded"]["studies"][0]["publicationInfo"]["pubmedId"]
                Publication_Date=data["_embedded"]["studies"][0]["publicationInfo"]["publicationDate"]
                Journal=data["_embedded"]["studies"][0]["publicationInfo"]["publication"]
                Title=data["_embedded"]["studies"][0]["publicationInfo"]["title"]
                Reported_Trait=data["_embedded"]["studies"][0]["diseaseTrait"]["trait"]
                Discovery_Sample_Number_And_Ancestry=str(data["_embedded"]["studies"][0]["ancestries"][0]["numberOfIndividuals"])+" "+data["_embedded"]["studies"][0]["ancestries"][0]["ancestralGroups"][0]["ancestralGroup"]
                Accession_Traits_link="https://www.ebi.ac.uk/gwas/rest/api/studies/{Study_Accession}/efoTraits".format(Study_Accession=Study_Accession)
                data = urllib.request.urlopen(Accession_Traits_link)
                data = json.load(data)
                for trait_data in data["_embedded"]["efoTraits"]:
                    Traits.append(trait_data["trait"])
                GWAS_Data_Dict={"rsid":rsid,"PubmedID":PubmedID,"Study_Accession":Study_Accession,"First_Author":First_Author,"Publication_Date":Publication_Date,"Journal":Journal,"Title":Title,"Reported_Trait":Reported_Trait,"Discovery_Sample_Number_And_Ancestry":Discovery_Sample_Number_And_Ancestry,"Traits":Traits}
                GWAS_df=GWAS_df.append(GWAS_Data_Dict,ignore_index=True)
            chosen_column=["rsid","PubmedID","Study_Accession","Reported_Trait"]
            truncated=False
            if len(GWAS_df.index)<len(rsid_list):
                print(len(GWAS_df.index),len(rsid_list))
                truncated=True
            return GWAS_df[chosen_column], truncated, servererror
        else:
            truncated=False
            result=pd.DataFrame(columns=["rsid","PubmedID","Study_Accession","First_Author","Publication_Date","Journal","Title","Reported_Trait","Discovery_Sample_Number_And_Ancestry","Traits"])
            chosen_column=["rsid","PubmedID","Study_Accession","Reported_Trait"]
            return result[chosen_column], truncated, servererror

    except urllib.error.HTTPError as e:
        servererror=True
        truncated=False
        
        print("urllib.error.HTTPError: ",e)
        result=pd.DataFrame(columns=["rsid","PubmedID","Study_Accession","First_Author","Publication_Date","Journal","Title","Reported_Trait","Discovery_Sample_Number_And_Ancestry","Traits"])
        chosen_column=["rsid","PubmedID","Study_Accession","Reported_Trait"]
        return result[chosen_column], truncated, servererror

"""
chromosome='9'
GenomeLocation='136881912-136926607'
GWAS_df=fetch_gwascatalog(chromosome,GenomeLocation)
print(GWAS_df)#Empty DataFrameが返ってくる（TRAF2にGWAS情報はない）

chromosome='5'
GenomeLocation='157142933-157255191'
GWAS_df=fetch_gwascatalog(chromosome,GenomeLocation)
print(GWAS_df)

chromosome='6'
GenomeLocation='132743870-132763459'
GWAS_df=fetch_gwascatalog(chromosome,GenomeLocation)
print(GWAS_df)
"""





def get_mgi_result(geneid):
    sqlite_path = './db/MGI.db'

    connection = sqlite3.connect(sqlite_path)
    #connection.row_factory = sqlite3.Row
    cursor = connection.cursor()
    res1 = cursor.execute('SELECT "MGI Marker Accession ID" FROM HMD_HumanPhenotype WHERE "Human Entrez Gene ID"=?',(geneid,))
    res1_tmp=res1.fetchone()
    mgi_acc_id=""
    print(res1_tmp)
    print(res1_tmp is not None)
    if res1_tmp is not None:
        mgi_acc_id=res1_tmp[0]

    res2 = cursor.execute('SELECT "Allelic Composition","Mammalian Phenotype ID","PubMed ID" FROM MGI_GenePheno WHERE "MGI Marker Accession ID"=?',(mgi_acc_id,))
    #print(res2.fetchall())
    allele_list=res2.fetchall()
    table=[]
    for allele in allele_list:
        mpid=allele[1]
        res3 = cursor.execute('SELECT "phenotype_term" FROM VOC_MammalianPhenotype WHERE "Mammalian Phenotype ID"=?',(mpid,))
        table.append((allele[0],allele[1],res3.fetchone()[0],allele[2]))
        #print(allele,res3.fetchone())
    result={"mgi_acc_id":mgi_acc_id,"table":table}
    return result


def add_params2url(url,params):
    for k,i in params.items():
        url+=k+"="+i+"&"
    return url.strip("&")

def element2dfrow(element):
    result_dict=xmltodict.parse(str(element))
    result_dict=result_dict["element"]
    if type(result_dict["bodySystems"]["element"])!=str:
        result_dict["bodySystems"]=",".join(result_dict["bodySystems"]["element"])
    else:
        result_dict["bodySystems"]=result_dict["bodySystems"]["element"]
    if type(result_dict["standardDeviation"])!=str:
        result_dict["standardDeviation"]=np.nan
        result_dict["SDNull"]=True
    else:
        result_dict["SDNull"]=False
    #print(pd.DataFrame(data=[list(result_dict.values())],columns=list(result_dict.keys())))
    print(list(result_dict.keys()))
    return pd.DataFrame(data=[list(result_dict.values())],columns=list(result_dict.keys()))


"""
def fetch_BaseSpace(Symbol):
    #context = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
    context = ssl._create_unverified_context()
    #ssl._create_default_https_context = ssl._create_unverified_context
    
    url = 'https://japantobacco.ussc.informatics.illumina.com/c/nbapi/bodyatlas.api?'
    params = {
    "q": Symbol,
    "apikey": "5defa1ede2cf3a21b233315f07b0163b",
    "v": "0",
    "bodyatlastype": "CELL_TYPE",
    "bodyatlasview": "RANK",
    "source": "1",
    "fmt": "xml"
    }

    url=add_params2url(url,params)
    html = urllib.request.urlopen(url,context=context)
    soup = BeautifulSoup(html, "xml")
    elements=soup.find("result").contents
    result_df=pd.DataFrame()

    for e in elements[:20]:
        row=element2dfrow(e)
        #print(row)
        result_df=result_df.append(row, ignore_index=True)
    #print(result_df)
    result_df["standardDeviation"]=result_df["standardDeviation"].fillna(0)
    return result_df

"""
def fetch_BaseSpace(GeneSymbol):
    #https://japantobacco.ussc.informatics.illumina.com/c/search/ba/CCR6?type=gene&id=5043#view=rank
    #BodyAtlas: TISSUE, view: RANK, source: Array
    #https://japantobacco.ussc.informatics.illumina.com/c/search/ba/CCR6?type=gene&id=5043#source=2&view=rank&winId=27285--31003
    #BodyAtlas: TISSUE, view: RANK, source: RNA-seq Based (2)
    context = ssl._create_unverified_context()
    ssl._create_default_https_context = ssl._create_unverified_context
    param_list=[{"bodyatlastype": "CELL_TYPE","source": "1","param_name":"CELL_TYPE"},\
                    {"bodyatlastype": "TISSUE","source": "1","param_name":"TISSUE_Array_Based"},\
                    {"bodyatlastype": "TISSUE","source": "2","param_name":"TISSUE_RNA-seq_Based"}\
                ]
    BaseSpase_result={}
    for p in param_list:
        url = 'https://japantobacco.ussc.informatics.illumina.com/c/nbapi/bodyatlas.api?'
        params = {
        "q": GeneSymbol,
        "apikey": "5defa1ede2cf3a21b233315f07b0163b",
        "v": "0",
        "bodyatlastype": p["bodyatlastype"],
        "bodyatlasview": "RANK",
        "source": p["source"],
        "fmt": "xml"
        }

        url=add_params2url(url,params)
        print(url)
        html = urllib.request.urlopen(url)
        soup = BeautifulSoup(html, "xml")
        if soup.find("error"):
            print(soup.find("error").text)
            BaseSpase_result[p["param_name"]]=pd.DataFrame(columns=['bodySystems', 'conceptId', 'conceptLabel', 'controlExpression', 'direction', 'standardDeviation', 'tissueExpression', 'SDNull'])
            continue
        elements=soup.find("result").contents
        result_df=pd.DataFrame()

        #Only top 20 results are shown
        for e in elements[:20]:
            row=element2dfrow(e)
            #print(row)
            result_df=result_df.append(row, ignore_index=True)
        #print(result_df)
        result_df["standardDeviation"]=result_df["standardDeviation"].fillna(0)
        BaseSpase_result[p["param_name"]]=copy.deepcopy(result_df)
    return BaseSpase_result
"""
GeneSymbol="TRAF2"
BaseSpase_result=fetch_BaseSpace(GeneSymbol)
print(BaseSpase_result)

GeneSymbol="ABRAXAS2"
BaseSpase_result=fetch_BaseSpace(GeneSymbol)

print(BaseSpase_result)
"""
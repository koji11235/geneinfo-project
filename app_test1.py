# -*- coding: utf-8 -*-
from flask import Flask, render_template, request

from GeneInfoGetter import fetch_geneidentifier
from GeneInfoGetter import fetch_uniprot
from GeneInfoGetter import fetch_allignment
from GeneInfoGetter import fetch_gwascatalog
from GeneInfoGetter import get_mgi_result
from GeneInfoGetter import fetch_BaseSpace

import mygene


def searchterm2geneid(search_term):
    if search_term.isdigit():
        return search_term
    mg = mygene.MyGeneInfo()
    candidates=mg.query(search_term, size=5)
    #print(candidates)
    candidate_GeneID=[]
    for c in candidates["hits"]:
        if c["taxid"]==9606:
            candidate_GeneID.append({"entrezgene":c["entrezgene"],"symbol":c["symbol"],"score":c["_score"]})
    print(candidate_GeneID)
    return candidate_GeneID[0]["entrezgene"]

app = Flask(__name__)

@app.route("/")
def index():
    return render_template('index.html')

@app.route("/result")
def result():
    search_term = request.args.get('geneid', '')
    geneid=searchterm2geneid(search_term)
    print(geneid)
    fg_result=fetch_geneidentifier(geneid)
    Symbol=fg_result["Symbol"]
    fullname=fg_result["Description"]
    Aliases=fg_result["Aliases"]
    uniprotid=fg_result["Uniprot_ID"]
    OMIM_ID=fg_result["OMIM_ID"]
    OMIM_link=fg_result["OMIM_link"]
    homologeneid=fg_result["HomogeneID"]
    #homologene_link=fg_result["Homologene_link"]
    chromosome=fg_result["chromosome"]
    GenomeLocation=fg_result["GenomeLocation"]

    uniprot_result=fetch_uniprot(uniprotid)
    function=uniprot_result["function"]
    #sequence=uniprot_result["sequence"]

    fa_result=fetch_allignment(homologeneid)
    Identity_Protein_Mouse=fa_result["Identity_Protein_Mouse"]
    Identity_DNA_Mouse=fa_result["Identity_DNA_Mouse"]
    Identity_Protein_Rat=fa_result["Identity_Protein_Rat"]
    Identity_DNA_Rat=fa_result["Identity_DNA_Rat"]

    GWAS_df=fetch_gwascatalog(chromosome,GenomeLocation)

    MGI_result=get_mgi_result(geneid)
    mgi_acc_id=MGI_result["mgi_acc_id"]
    MGI_table=MGI_result["table"]

    BaseSpase_result=fetch_BaseSpace(Symbol)
    BaseSpase_result_len=14
    for param_name,BaseSpase_result_df in BaseSpase_result.items():
        BaseSpase_result[param_name]=BaseSpase_result_df.iloc[0:BaseSpase_result_len,:]
    
    #print(type(BaseSpase_result_df.iloc[0,3]))


    return render_template("index.html",
                            geneid=geneid,\
                            Symbol=Symbol,\
                            fullname=fullname,\
                            Aliases=Aliases,\
                            uniprotid=uniprotid,\
                            OMIM_ID=OMIM_ID,\
                            OMIM_link=OMIM_link,\
                            homologeneid=homologeneid,\
                            function=function,\
                            #sequence=sequence,\
                            Identity_Protein_Mouse=Identity_Protein_Mouse,\
                            Identity_DNA_Mouse=Identity_DNA_Mouse,\
                            Identity_Protein_Rat=Identity_Protein_Rat,\
                            Identity_DNA_Rat=Identity_DNA_Rat,\
                            GWAS_catalog=GWAS_df,\
                            mgi_acc_id=mgi_acc_id,\
                            MGI_table=MGI_table,\
                            BaseSpase_result=BaseSpase_result,\
                            BaseSpase_result_len=BaseSpase_result_len
                            )


if __name__ == "__main__":
    app.run(debug=True)

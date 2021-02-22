# -*- coding: utf-8 -*-
from flask import Flask, render_template, request
from GeneInfo import GeneInfo

app = Flask(__name__)

@app.route("/")
def index():
    return render_template('index.html')

@app.route("/result")
def result():
    entrez_id = request.args.get('entrez_id', '')
    #print(geneid)
    geneinfo = GeneInfo(entrez_id=entrez_id)
    geneinfo.fetch_mygene()

    geneinfo.fetch_uniprot()
    geneinfo.fetch_homologene()

    geneinfo.fetch_basespace()
    
    geneinfo.fetch_opentarget_genetic_association()
    geneinfo.get_mgi_result()

    geneinfo.fetch_opentarget_drug()

    # create charts -----------------------------------
    chart_ot_clinical_trial = geneinfo.plot_ot_clinical_trial_pie()
    chart_ot_drug2disease = geneinfo.plot_ot_drug2disease()
    chart_domain = geneinfo.plot_domains()
    

    return render_template(
        "index.html", 
        entrez_id=entrez_id,
        Symbol=geneinfo.symbol,
        fullname=geneinfo.fullname,
        Aliases=', '.join(geneinfo.aliases),
        ids=geneinfo.ids,
        links=geneinfo.links,
        function=geneinfo.uniprot_result['function'],
        #sequence=sequence,
        sequence_homologies=geneinfo.sequence_homologies,
        opentarget_genetic_association=geneinfo.opentarget_genetic_association_table,
        mgi_table=geneinfo.mgi_result_table,
        basespase_result=geneinfo.basespase_result,
        opentarget_drug_table=geneinfo.opentarget_drug_table,
        # figs -------------------------------------------------------
        chart_ot_clinical_trial=chart_ot_clinical_trial,
        chart_ot_drug2disease=chart_ot_drug2disease,
        chart_domain=chart_domain,
        is_domainplot_omitted=geneinfo.is_domainplot_omitted,
        MAX_PLOT_NUM_DOMAINS=geneinfo.MAX_PLOT_NUM_DOMAINS
    )


if __name__ == "__main__":
    app.run(debug=True)

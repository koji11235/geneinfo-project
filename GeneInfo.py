from bs4 import BeautifulSoup
import pandas as pd
import numpy as np
import time
import lxml.html
import copy
import json
import sqlite3
import xmltodict
#import ssl
import mygene
import requests

import plotly
import plotly.graph_objects as go
import plotly.express as px
import networkx as nx
from matplotlib import cm


from math import log10
from math import floor



class GeneInfo(object):
    def __init__(self, entrez_id):
        api_keys_path = 'api_keys.csv'
        self.mgi_db_path = './db/MGI.db'

        api_key_df = pd.read_csv(api_keys_path, header=None, index_col=0)
        self.api_keys = dict(zip(api_key_df.index.values, api_key_df[1]))

        self.ids = {}
        self.links = {}
        self.ids["entrez"] = entrez_id

        self.symbol = None
        self.fullname = None
        self.aliases = None



        self.uniprot_result = {}
        self.alignment = {}
        self.basespase_result = {}

        self.opentarget_genetic_association_table = pd.DataFrame() #column namesはoutputと統一させる
        self.opentarget_drug_table = pd.DataFrame() #column namesはoutputと統一させる
        self.mgi_result_table = pd.DataFrame() #column namesはoutputと統一させる

        # 情報を取得してきたかどうかのフラグだけど必要ない？
        # self._is_geneidentifier_fetched = False
        # self._is_uniprot_fetched = False

        

    # def fetch_geneidentifier(self):
    #     pass

    def fetch_mygene(self):
        fields = [
            'symbol', 'name', 'alias',
            'ensembl.gene',
            'uniprot.Swiss-Prot',
            'MIM',
            'homologene.id',
            'genomic_pos',
            'exons',
            #'interpro' #domain information
        ]
        self.links['mygene'] = f"http://mygene.info/v3/gene/{self.ids['entrez']}?fields={','.join(fields)}"
        mygene_response = requests.get(self.links['mygene']).json()

        self.symbol = mygene_response['symbol']
        self.fullname = mygene_response['name']
        self.aliases = mygene_response['alias']
        self.ids['ensembl'] = mygene_response['ensembl']['gene']

        self.ids['uniprot'] = mygene_response['uniprot']['Swiss-Prot']
        self.ids['omim'] = mygene_response['MIM']
        self.links['omim'] = f"https://www.omim.org/entry/{self.ids['omim']}"
        self.ids['homologene'] = mygene_response['homologene']['id']
        self.links['homologene'] = f"https://www.ncbi.nlm.nih.gov/homologene/{self.ids['homologene']}"
        self.genome_location = mygene_response['genomic_pos']
        self.exons = mygene_response['exons']


    def fetch_geneidentifier(self):
        """
        Deprecated. Please use fetch mygene instead.
        retrieve basic gene informations from geneidentifier.

        必要情報：Entrez Gene ID (GeneInfoオブジェクト生成時に引数として渡す)
        """
        self.links['gene_identifier'] = f"http://biogps.org/ext/symatlasbar/?geneid={self.ids['entrez']}&hidespecies=1"
        geneidentifier_xml = requests.get(self.links['gene_identifier']).text
        geneidentifier_soup = BeautifulSoup(geneidentifier_xml, "lxml")
        self.symbol = geneidentifier_soup.find_all("tr")[0].find(class_="annval").text
        self.description = geneidentifier_soup.find_all("tr")[1].find(class_="annval").text
        self.aliases = geneidentifier_soup.find_all("tr")[3].find(class_="annval").text.strip()

        self.ids['uniprot'] = None
        self.ids['omim'] = None
        self.links['omim'] = None
        self.ids['homologen'] = None
        self.links['homologene'] = None
        self.chromosome = None
        self.genome_location = None

        for a in geneidentifier_soup.find_all("a"):
            if "http://www.uniprot.org/uniprot/" in a.get("href"):
                self.ids['uniprot']=a.text
            
            if "http://omim.org/entry/" in a.get("href"):
                self.ids['omim'] = a.text
                self.links['omim'] = a.get("href")

            if "http://www.ncbi.nlm.nih.gov/homologene/" in a.get("href"):
                self.ids['homologene'] = a.text
                self.links['homologene'] = a.get("href")
            
            if "http://genome.ucsc.edu/cgi-bin/hgTracks?db=" in a.get("href"):
                genome_location = a.text
                self.chromosome = genome_location.split(":")[0].strip("chr")
                self.genome_location=  genome_location.split(":")[1]
                break

        self._is_geneidentifier_fetched = True

    # previous ver -------------------------------------------------------------------------------
    # def fetch_uniprot(self):
    #     """
    #     UniprotからFunction情報（とアミノ酸配列）の取得
    #     必要情報：Uniprot ID（fetch_geneidentifierで取得可能）
    #     """
    #     if self.ids['uniprot'] is None:
    #         print("Uniprot ID is not defined. Please run .fetch_geneidentifier() method first.")
    #         raise ValueError 

    #     self.links['uniprot'] = f"https://www.uniprot.org/uniprot/{self.ids['uniprot']}.xml"
    #     data = urllib.request.urlopen(self.links['uniprot'])
    #     xml = data.read()
    #     soup = BeautifulSoup(xml, "xml")
    #     self.sequence = soup.find_all("sequence")[-1].text.strip()
    #     #print("sequence: ",sequence)
    #     if soup.find_all("comment") != []:
    #         self.function = soup.find_all("comment")[0].text.strip()
    #     else:
    #         self.function = ""

    #     self._is_uniprot_fetched = True

    def fetch_uniprot(self):
        self.links['uniprot'] = f"https://www.ebi.ac.uk/proteins/api/proteins/{self.ids['uniprot']}"
        uniprot_response = requests.get(self.links['uniprot'], headers={ "Accept" : "application/json"}).json()
        
        self.uniprot_result['function'] = uniprot_response['comments'][0]['text'][0]['value']
        self.uniprot_result['sequence'] = uniprot_response['sequence']['sequence']
        self.uniprot_result['domains'] = [feature for feature in uniprot_response['features'] if feature['category'] == 'DOMAINS_AND_SITES']
        

    def fetch_homologene(self):
        self.links['homologene'] = f"https://www.ncbi.nlm.nih.gov/homologene/{self.ids['homologene']}?report=alignmentscores"
        # data = urllib.request.urlopen(self.links['homologene'])
        # xml = data.read()
        # soup = BeautifulSoup(xml, "xml")
        print(self.links['homologene'])

        homologene_xml = requests.get(self.links['homologene']).text
        homologene_soup = BeautifulSoup(homologene_xml, "xml")

        alignment_table = homologene_soup.find("table",class_="alignment_table").find_all("tr")
        #print(alignment_table)
        for idx, content in enumerate(alignment_table):

            #print(content)
            #print(content.find_all(class_="top_species"))
            if content.find_all(class_="top_species")!=[]:
                if content.find(class_="top_species").text=="H.sapiens":
                    H_sapiens_begin = copy.deepcopy(idx)
                    break
        #print(H_sapiens_begin)

        for idx, content in enumerate(alignment_table[H_sapiens_begin+1:]):
            #print(content)
            #print(content.find_all(class_="top_species"))
            if content.find_all(class_="top_species")!=[]:
                H_sapiens_end = copy.deepcopy(idx) + H_sapiens_begin + 1
                break
        #print(H_sapiens_end)
        #print(alignment_table[H_sapiens_begin:H_sapiens_end])
        identity_protein_mouse = ""
        identity_dna_mouse = ""
        identity_protein_rat = ""
        identity_dna_rat = ""
        for content in alignment_table[H_sapiens_begin + 1: H_sapiens_end]:
            if content.find(class_="second_species").text == "vs. M.musculus":
                identity_protein_mouse = content.find(class_="alignment_table_col3").text
                identity_dna_mouse = content.find(class_="alignment_table_col4").text
            if content.find(class_="second_species").text == "vs. R.norvegicus":
                identity_protein_rat = content.find(class_="alignment_table_col3").text
                identity_dna_rat = content.find(class_="alignment_table_col4").text
        #print("identity_protein_mouse;",identity_protein_mouse)
        #print("identity_dna_mouse;",identity_dna_mouse)
        #print("identity_protein_rat;",identity_protein_rat)
        #print("identity_dna_rat;",identity_dna_rat)

        self.sequence_homologies = {
            "identity_protein_mouse": identity_protein_mouse,
            "identity_dna_mouse": identity_dna_mouse,
            "identity_protein_rat": identity_protein_rat,
            "identity_dna_rat": identity_dna_rat
        }

    def fetch_gwascatalog(self):
        """
        Deprecated. Please use fetch_opentarget_genetic_association instead.
        """
        # opentarget APIで書き直す
        pass

    def fetch_opentarget_genetic_association(self):

        def round_sig(x, sig=2):
            return np.around(x, sig-int(floor(log10(abs(x))))-1)

        self.links['opentarget_genetic_association'] = f"https://platform-api.opentargets.io/v3/platform/public/evidence/filter?target={self.ids['ensembl']}&datasource=ot_genetics_portal&size=1000"
        ot_genetic_association = requests.get(self.links['opentarget_genetic_association']).json()

        ids = []
        source_links = []
        rs_ids = []
        reported_traits = []
        efo_ids = []
        efo_links = []
        study_links = []
        study_articles = []
        study_article_authors = []
        study_article_years = []
        variant2disease_resource_score_types = []
        variant2disease_resource_score_values = []
        variant_types = []
        variant_type_links = []
        gene2variant_resource_score_types = []
        gene2variant_resource_scores = []

        for data in ot_genetic_association["data"]:
            ids.append(data['variant']['id'])
            source_links.append(data['variant'].get('source_link', ""))
            rs_ids.append(data['variant'].get('rs_id'))
            reported_traits.append(data['disease'].get('reported_trait', data['disease'].get('name',''))) #genetic disorder のときreported traitsでなくてnameになる
            efo_ids.append(data['disease']['id'])
            efo_links.append(data['disease']['efo_info']['efo_id'])
            study_links.append(data['evidence']['variant2disease'].get('study_link'))

            study_articles.append(data['evidence']['variant2disease']['provenance_type'].get('literature', {}).get('references',[{'lit_id':""}])[0]['lit_id'])
            study_article_authors.append(data['evidence']['variant2disease']['provenance_type'].get('literature', {}).get('references',[{'author':""}])[0].get('author'))
            study_article_years.append(data['evidence']['variant2disease']['provenance_type'].get('literature', {}).get('references',[{'year':""}])[0].get('year'))

            variant2disease_resource_score_types.append(data['evidence']['variant2disease']['resource_score']['type'])
            variant2disease_resource_score_values.append(data['evidence']['variant2disease']['resource_score']['value'])
            variant_types.append(data['evidence']['gene2variant'].get('consequence_code', 'others')) # SO_0001583 misssense variantのときconsequence_codeがない
            variant_type_links.append(data['evidence']['gene2variant']['functional_consequence']) 
            gene2variant_resource_score_types.append(data['evidence']['gene2variant']['resource_score']['type'])
            gene2variant_resource_scores.append(data['evidence']['gene2variant']['resource_score']['value'])

        opentarget_genetic_association_columns = [
            'ids','source_links','rs_ids',
            'reported_traits',
            'efo_ids', 'efo_links',
            'study_links', 'study_articles', 'study_article_authors', 'study_article_years',
            'variant2disease_resource_score_types', 'variant2disease_resource_score_values',
            'variant_types', 'variant_type_links',
            'gene2variant_resource_score_types', 'gene2variant_resource_scores'
        ]
        self.opentarget_genetic_association_table = pd.DataFrame({
            'ids': ids,
            'source_links': source_links,
            'rs_ids': rs_ids,
            'reported_traits': reported_traits,
            'efo_ids': efo_ids,
            'efo_links': efo_links,
            'study_links': study_links,
            'study_articles': study_articles,
            'study_article_authors': study_article_authors,
            'study_article_years': study_article_years,
            'variant2disease_resource_score_types': variant2disease_resource_score_types,
            'variant2disease_resource_score_values': [round_sig(score, sig=2) for score in variant2disease_resource_score_values],
            'variant_types': variant_types,
            'variant_type_links': variant_type_links,
            'gene2variant_resource_score_types': gene2variant_resource_score_types,
            'gene2variant_resource_scores': [round_sig(score, sig=2) for score in gene2variant_resource_scores]},
            columns=opentarget_genetic_association_columns
        )




    def get_mgi_result(self):
        # create connection -----------------------------------------------------------------------
        #sqlite_path = './db/MGI.db'
        connection = sqlite3.connect(self.mgi_db_path)
        #connection.row_factory = sqlite3.Row
        cursor = connection.cursor()

        # get mgi accession id -----------------------------------------------------------------------
        res1 = cursor.execute('SELECT "MGI Marker Accession ID" FROM HMD_HumanPhenotype WHERE "Human Entrez Gene ID"=?', (self.ids['entrez'],))
        res1_tmp = res1.fetchone()
        mgi_acc_id = ""
        print(res1_tmp)
        print(res1_tmp is not None)
        if res1_tmp is not None:
            mgi_acc_id = res1_tmp[0]


        # create allele table -----------------------------------------------------------------------
        res2 = cursor.execute('SELECT "Allelic Composition","Mammalian Phenotype ID","PubMed ID" FROM MGI_GenePheno WHERE "MGI Marker Accession ID"=?',(mgi_acc_id,))
        #print(res2.fetchall())
        allele_list = res2.fetchall()
        table = []
        for allele in allele_list:
            mpid = allele[1]
            res3 = cursor.execute('SELECT "phenotype_term" FROM VOC_MammalianPhenotype WHERE "Mammalian Phenotype ID"=?', (mpid,))
            table.append((allele[0], allele[1], res3.fetchone()[0], allele[2]))
            #print(allele,res3.fetchone())
        #result={"mgi_acc_id": mgi_acc_id, "table": table}
        # store the results in object -----------------------------------------------------------------------
        self.ids['mgi'] = mgi_acc_id
        self.links['mgi'] = f"http://www.informatics.jax.org/marker/{mgi_acc_id}"
        self.mgi_result_table = pd.DataFrame(table, columns=['alele_composition', 'mammalian_phenotype_id', 'phenotype', 'pubmed_id'])
        self.mgi_result_table['mammalian_phenotype_link'] = 'http://www.informatics.jax.org/vocab/mp_ontology/' + self.mgi_result_table['mammalian_phenotype_id']
        self.mgi_result_table['pubmed_link'] = 'https://www.ncbi.nlm.nih.gov/pubmed/' + self.mgi_result_table['pubmed_id']





    def fetch_basespace(self):
        """
        Fetch gene expression data from BaseSpace by using API.
        This function retrieve information by three different conditions(parameters).
            1. gene expression by cell type (param_name: "CELL_TYPE")
            2. gene expression by tissue (Array Based) (param_name: "TISSUE_Array_Based")
            3. gene expression by tissue (RNA-seq Based) (param_name: "TISSUE_RNA-seq_Based")

        args:

        return: 
            dictionary of the results.
            You can access it by self.basespase_result
            Its keys are names of conditions. ("CELL_TYPE", "TISSUE_Array_Based", "TISSUE_RNA-seq_Based")
            Its values are dataframes which contain the gene expression data of each condition.

        """

        def add_params2url(url,params):
            for k,i in params.items():
                url += k + "=" + i + "&"
            return url.strip("&")

        def element2dfrow(element):
            result_dict=xmltodict.parse(str(element))
            result_dict=result_dict["element"]
            if type(result_dict["bodySystems"]["element"]) != str:
                result_dict["bodySystems"] = ",".join(result_dict["bodySystems"]["element"])
            else:
                result_dict["bodySystems"] = result_dict["bodySystems"]["element"]
            if type(result_dict["standardDeviation"]) != str:
                result_dict["standardDeviation"] = np.nan
                result_dict["SDNull"] = True
            else:
                result_dict["SDNull"] = False
            #print(pd.DataFrame(data=[list(result_dict.values())],columns=list(result_dict.keys())))
            #print(list(result_dict.keys()))
            return pd.DataFrame(data=[list(result_dict.values())], columns=list(result_dict.keys()))
        
        # ssl._create_default_https_context = ssl._create_unverified_context

        param_list = [{"bodyatlastype": "CELL_TYPE", "source": "1", "param_name": "CELL_TYPE"},
                        {"bodyatlastype": "TISSUE", "source": "1", "param_name": "TISSUE_ARRAY_BASED"},
                        {"bodyatlastype": "TISSUE", "source": "2", "param_name": "TISSUE_RNA-SEQ_BASED"}
                    ]
        show_categories = 12 # How many categories are shown in page
        
        for p in param_list:
            url = 'https://japantobacco.ussc.informatics.illumina.com/c/nbapi/bodyatlas.api?'
            params = {
            "q": self.symbol,
            "apikey": self.api_keys["basespace"],
            "v": "0",
            "bodyatlastype": p["bodyatlastype"],
            "bodyatlasview": "RANK",
            "source": p["source"],
            "fmt": "xml"
            }

            url = add_params2url(url,params)
            #print(url)
            #html = urllib.request.urlopen(url)
            basespace_xml = requests.get(url).text
            basespace_soup = BeautifulSoup(basespace_xml, "xml")

            result_df = pd.DataFrame(
                    columns=['bodySystems', 'conceptId', 'conceptLabel', 
                    'controlExpression', 'direction', 'standardDeviation', 
                    'tissueExpression', 'SDNull'])

            if basespace_soup.find("error"):
                print(basespace_soup.find("error").text)
                continue

            elements = basespace_soup.find("result").contents
            #Only top 20 results are shown
            for e in elements[:20]:
                row = element2dfrow(e)
                #print(row)
                result_df = result_df.append(row, ignore_index=True)
            #print(result_df)
            result_df["standardDeviation"] = result_df["standardDeviation"].fillna(0)
            self.basespase_result[p["param_name"]] = result_df.iloc[0: show_categories, :]

    def fetch_geneexpressionatlas(self):
        #
        pass

    def fetch_opentarget_drug(self):
        self.links['opentarget_drug'] = f"https://platform-api.opentargets.io/v3/platform/public/evidence/filter?target={self.ids['ensembl']}&datatype=known_drug&size=1000"
        ot_drug_data = requests.get(self.links['opentarget_drug']).json()
        
        # create empty lists ----------------------------------------
        drug = []
        drug_url = []

        disease = []
        disease_id = []
        disease_url = []

        phase = []
        status = []

        mechanism_of_action = []
        drug_molecule_type = []
        action_type = []
        activity = []

        target = []
        target_class = []
        evidence_curated_from = []
        clinical_trial_url = []

        # append contents ----------------------------------------
        for data in ot_drug_data["data"]:
            #print(data["disease"]['efo_info']['therapeutic_area'])
            drug.append(data["drug"]["molecule_name"])
            drug_url.append(data["unique_association_fields"]["chembl_molecule"])

            disease.append(data["disease"]["efo_info"]["label"])
            disease_id.append(data["disease"]["id"])
            disease_url.append(data["disease"]['efo_info']['efo_id'])

            phase.append(data['evidence']['drug2clinic']['clinical_trial_phase']['label'])
            status.append(data['evidence']['drug2clinic'].get('status', ''))

            mechanism_of_action.append(data['evidence']['target2drug']['mechanism_of_action'])
            drug_molecule_type.append(data["drug"]['molecule_type'])
            action_type.append(data["evidence"]['target2drug']['action_type'])
            activity.append(data["target"]['activity'])

            target.append(data['target']['gene_info']['symbol'])
            target_class.append(data['target']['target_class'][0])
            evidence_curated_from.append(data['evidence']['drug2clinic']['urls'][0]['nice_name'])
            clinical_trial_url.append(data['evidence']['drug2clinic']['urls'][0]['url'])

        # create result dataframe ----------------------------------------
        opentarget_drug_table_columns = [
            'Drug', 'Drug URL', 
            'Disease', 'Disease ID', 'Disease URL', 
            'Phase', 'Status', 
            'Mechanism of action', 'Drug Molecule Type', 'Action type', 'Activity', 
            'Target', 'Target class', 
            'Evidence curated from', 'Clinical trial url']

        self.opentarget_drug_table = pd.DataFrame(
                data={
                    'Drug': drug,
                    'Drug URL': drug_url, 
                    'Disease': disease, 
                    'Disease ID': disease_id, 
                    'Disease URL': disease_url, 
                    'Phase': phase,
                    'Status': status, 
                    'Mechanism of action': mechanism_of_action, 
                    'Drug Molecule Type': drug_molecule_type,
                    'Action type': action_type, 
                    'Activity': activity, 
                    'Target': target, 
                    'Target class': target_class,
                    'Evidence curated from': evidence_curated_from, 
                    'Clinical trial url': clinical_trial_url},
                columns=opentarget_drug_table_columns
        )
        # remove dupulicates by 'Target','Chemble id','Disease'
        self.opentarget_drug_table.drop_duplicates(subset=['Target','Drug URL','Disease'], inplace=True)



    def plot_ot_clinical_trial_pie(self):
        drug_data = self.opentarget_drug_table

        d_order = {'Phase 0': 0, 'Phase I': 1, 'Phase II': 2, 'Phase III': 3, 'Phase IV': 4}

        phase_df = pd.DataFrame(drug_data['Phase'].value_counts())
        phase_df['order'] = phase_df.index.map(d_order)
        phase_df = phase_df.sort_values('order')
        #print(phase_df)

        layout = go.Layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            xaxis = go.layout.XAxis(
                title = '',
                showticklabels=False),
            yaxis = go.layout.YAxis(
                title = '',
                showticklabels=False)
        )

        px.colors.sequential.Blues
        fig = go.Figure(
            data=[
                go.Pie(
                    labels=phase_df.index, 
                    values=phase_df['Phase'], 
                    hole=.6, 
                    marker_colors=px.colors.sequential.Blues[3:8],
                    sort=False
                )
            ],
            layout=layout
        )

        fig.update_traces(
        textposition='outside',
        textinfo='label+value',
        showlegend=False,
        ).update_layout(
            annotations=[
                dict(
                    x=0.5, y=0.5,
                    xanchor= 'center',
                    yanchor='middle',
                    text= '<b>Clinical Trials</b>',
                    font=dict(family="Helvetica", size=16),
                    showarrow=False
                )
            ],

        )
        
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    def plot_ot_drug2disease(self):
        drug_data = self.opentarget_drug_table
        drug2disease = drug_data.loc[:,['Drug', 'Disease']]
        #print(drug2disease)
        drug_nodes = drug2disease.Drug.unique().tolist()
        disease_nodes = drug2disease.Disease.unique().tolist()

        G = nx.Graph()
        for node in drug_nodes + disease_nodes:
            G.add_node(node)

        #pos = nx.spring_layout(G, k=0.3, seed=1)

        for _, row in drug2disease.iterrows():
            G.add_edge(row.Drug, row.Disease)

        #top = nx.bipartite.sets(G)[0]
        #pos = nx.bipartite_layout(G, top)
        shells = [drug_nodes, disease_nodes]
        pos = nx.shell_layout(G, shells)
        for node in G.nodes():
            G.nodes[node]["pos"] = pos[node]

        node_x = []
        node_y = []
        for n in G.nodes():
            x, y = G.nodes[n]["pos"]
            node_x.append(x)
            node_y.append(y)

        edge_x = []
        edge_y = []
        for e in G.edges():
            x0, y0 = G.nodes[e[0]]["pos"]
            x1, y1 = G.nodes[e[1]]["pos"]
            edge_x.append(x0)
            edge_y.append(y0)
            edge_x.append(x1)
            edge_y.append(y1)
            edge_x.append(None)
            edge_y.append(None)
        text = drug_nodes + disease_nodes

        nodes = go.Scatter(
            x=node_x,
            y=node_y,
            mode="markers+text",
            text=text,
            textposition="middle center",
            #marker=dict(size=30, line=dict(width=2)),
            marker=dict(
                size=[40] * len(drug_nodes) + [20] * len(disease_nodes),
                color=[px.colors.sequential.Blues[4]] * len(drug_nodes) + [px.colors.sequential.Blues[3]] * len(disease_nodes),
                opacity=1,
            ),
            hovertemplate="%{text}<extra></extra>",
            showlegend=False
        )

        px.colors.sequential.Blues[3:8]
        # '#fb9f3a' '#f0f921'

        edges = go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line=dict(
                width=2, 
                color='rgba(50,50,50,0.3)', 
            ),
            showlegend=False
        )
        layout = go.Layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            xaxis = go.layout.XAxis(
                title = '',
                showticklabels=False),
            yaxis = go.layout.YAxis(
                title = '',
                showticklabels=False)
        )
        fig = go.Figure(data=[edges, nodes], layout=layout)
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    def plot_domains(self):
        keys = ['DOMAIN', 'REGION', 'MOTIF', 'SITE', 'COILED', 'ZN_FING', 'REPEAT', 'BINDING', 'DNA_BIND', 'METAL', 'CA_BIND', 'NP_BIND', 'ACT_SITE']
        colormap = cm.get_cmap('tab10', 12)
        colors = [
            {'linecolor': f"rgba({colormap(i)[0]*255},{colormap(i)[1]*255},{colormap(i)[2]*255},1)",
            'fillcolor': f"rgba({colormap(i)[0]*255},{colormap(i)[1]*255},{colormap(i)[2]*255},0.5)"} for i in range(len(keys))
        ]
        domain2color = dict(zip(keys, colors))

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=[0, len(self.uniprot_result['sequence'])],
            y=[0, 0],
            mode="lines",
            line=dict(
                width=2, 
                color='rgba(50,50,50,0.7)', 
                ),
            showlegend=False
            )
        )

        for i, domain in enumerate(self.uniprot_result['domains']):
            begin = int(domain['begin'])
            end = int(domain['end'])
            #annotation_y = 1 if (i ％ 2) == 1　else -1
            fig.add_trace(
                go.Scatter(
                    x=[begin,begin,end,end,begin], 
                    y=[-1,1,1,-1,-1], 
                    fill="toself",
                    mode='lines',
                    name='',
                    text=f"{domain['type']}<br>{domain.get('description','')}<br>begin: {begin}<br>end: {end}",
                    opacity=1,
                    fillcolor=domain2color[domain['type']]['fillcolor'],
                    line=dict(
                        color=domain2color[domain['type']]['linecolor'],
                        width=3,
                    ),
                    showlegend=False
                )
            ).add_annotation(
                    x=(begin + end)/2,
                    y=-1.3,
                    xref="x",
                    yref="y",
                    text=f"{domain['type']}",
                    showarrow=False,
            )

        fig.update_layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            xaxis = go.layout.XAxis(
                title = '',
                showticklabels=True
            ),
            yaxis = go.layout.YAxis(
                title = '',
                showticklabels=False,
                range=[-1.5,1.5]
            ),
            width=1000,
            height=150,
            margin=dict(
                l=0,
                r=0,
                b=0,
                t=10,
                pad=10
            ),
        )
        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
 


                


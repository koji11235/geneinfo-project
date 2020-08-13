from GeneInfo import GeneInfo
import numpy
import sys
import json
import urllib.request
import ssl
import mygene
import pandas as pd

entrez_id = "3716"
print("create instance")
geneinfo = GeneInfo(entrez_id=entrez_id)
print("instance created")
print("\n\n\n")

# geneinfo.fetch_geneidentifier()
# print(geneinfo.ids['uniprot'])

print("mygene","------------------------------------")
geneinfo.fetch_mygene()
#print(geneinfo.symbol)
print(geneinfo.aliases)

#print("genome location", geneinfo.genome_location)
print("\n\n\n")

# print("uniprot","------------------------------------")
# geneinfo.fetch_uniprot()
# print(geneinfo.uniprot_result)
# print("\n\n\n")

# print("alignment","------------------------------------")
# geneinfo.fetch_alignment()
# print(geneinfo.alignment)
# print("\n\n\n")


# print("basespace","------------------------------------")
# geneinfo.fetch_basespace()
# print(geneinfo.basespase_result)
# print(len(geneinfo.basespase_result['CELL_TYPE']))
# print("\n\n\n")

# print("opentarget drug","------------------------------------")
# geneinfo.fetch_opentarget_drug()
# print(geneinfo.opentarget_drug_table)
# print("\n\n\n")

# print("opentarget genetic_association","------------------------------------")
# geneinfo.fetch_opentarget_genetic_association()
# print(geneinfo.opentarget_genetic_association_table)
# #self.links['gene_identifier'] = f"http://biogps.org/ext/symatlasbar/?geneid={self.ids["entrez"]}&hidespecies=1"
# print("\n\n\n")

# print("mgi","------------------------------------")
# geneinfo.get_mgi_result()
# print(geneinfo.mgi_result_table)

# print(geneinfo.opentarget_drug_table)
# geneinfo.opentarget_drug_table.to_csv("ot_drug.csv")

import mygene

mg = mygene.MyGeneInfo()

search_term="dock2"
candidates=mg.query(search_term, size=5)
candidate_GeneID=[]

for c in candidates["hits"]:
    if c["taxid"]==9606:
        candidate_GeneID.append({"entrezgene":c["entrezgene"],"symbol":c["symbol"],"score":c["_score"]})
print(candidate_GeneID)


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

print(searchterm2geneid("dock2"))
print(searchterm2geneid("1704"))
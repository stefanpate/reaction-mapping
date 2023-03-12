import pythoncyc

pwy_ids = ['PWY1A0-6325']
starters = ['|PROPIONYL-COA|', '|MALONYL-COA|', '|ACETYL-COA|', '|METHYL-MALONYL-COA|', '|D-METHYL-MALONYL-COA|']

def check_for_starter(pwy_frame):
    '''
    Searches a given metacyc pathway frame
    for starters in an attribute called
    'reaction_layout'.
    '''
    rxn_layout = pwy_frame.reaction_layout
    for rxn in rxn_layout:
        for elt in rxn:
            if type(elt) == str:
                if elt in starters:
                    return True
            elif type(elt) == list:
                for sub_elt in elt:
                    if sub_elt in starters:
                        return True
    return False

meta = pythoncyc.select_organism('meta')
pwy_frames = [meta[elt] for elt in pwy_ids]

rxn_ids = []
for elt in pwy_frames:
    found_starter = False
    found_starter = check_for_starter(elt)
    rxn_ids += elt.reaction_list
    if not found_starter:
        for link_list in elt.pathway_links.values():
            for link in link_list:
                if 'PWY' in link:
                    predecessor_pathway = meta[link]
                    found_starter = check_for_starter(predecessor_pathway)
                    if found_starter:
                        rxn_ids += predecessor_pathway.reaction_list
print(pwy_ids)
print(rxn_ids)
    

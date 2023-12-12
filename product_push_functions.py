# product push functions
def get_reaction_information(rxn, model):
    rxn_df=model.reactions.get_by_id(rxn)
    reversibility=rxn_df.reversibility
    rxn_reactants, rxn_products=rxn_df.reactants, rxn_df.products
    rxn_reactant_names=[i.id for i in rxn_reactants]
    rxn_product_names=[i.id for i in rxn_products]
    all_metabolites=rxn_reactants+rxn_products
    all_metabolite_names=rxn_reactant_names+rxn_product_names
    rxn_object={'df':rxn_df, 'rev':reversibility, 'r':rxn_reactants, 'p':rxn_products, 'rn': rxn_reactant_names, 'pn':rxn_product_names, 'm':all_metabolites, 'mn':all_metabolite_names}
    return rxn_object

# conducts product push for all mets in a reaction (based on reversibility)
def f1_product_push(rxn, model, graph, f1_rev=False, f2_rev=False):
    rxn_object=get_reaction_information(rxn, model)
    reversibility=rxn_object['rev']
    mets=rxn_object['pn']
    if (f1_rev and reversibility):
        mets=rxn_object['mn']
    rxn_neighbors=[n for n in graph.neighbors(rxn)]
    rxn_p_push={}
    for m in mets:
        p_push=[]
        for nxt_rxn in rxn_neighbors:
            p_push.append(f2_product_push(m, nxt_rxn, model, f2_rev=f2_rev))
        rxn_p_push[m]=sum(p_push)
    return rxn_p_push
        
# check how many times the given product appears in the next reaction
def f2_product_push(p, rxn, model, f2_rev=False):
    rxn_object=get_reaction_information(rxn, model)
    reversibility=rxn_object['rev']
    mets=rxn_object['rn']
    if (f2_rev and reversibility):
        mets=rxn_object['mn']
    p_push=sum([1 if (m==p) else 0 for m in mets])
    return p_push

def product_push_main(lofo_dict, model, graph, func, mets_to_ignore=[], f1_rev=False, f2_rev=False, ignore_zero_push=False, return_func_df=False):
    product_push_data={}
    product_push_func={}
    for rxn, obj in lofo_dict.items():
        rxn_p_push=f1_product_push(rxn, model, graph, f1_rev=f1_rev, f2_rev=f2_rev)
        if (len(mets_to_ignore)!=0 and len(rxn_p_push)!=0):
                rxn_p_push={key: value for key, value in rxn_p_push.items() if key not in mets_to_ignore}
        if (ignore_zero_push and len(rxn_p_push)!=0):
            rxn_p_push={key: value for key, value in rxn_p_push.items() if value!=0}
        if (len(rxn_p_push)==0):
            rxn_p_push_func=0
        else:
            rxn_p_push_func=func(list(rxn_p_push.values()))
        
        product_push_data[rxn]=rxn_p_push
        product_push_func[rxn]=rxn_p_push_func
        
    if (return_func_df):
        product_push_func=pd.DataFrame.from_dict(product_push_func, orient='index')
        product_push_func.columns=["p_push"]

    return product_push_data, product_push_func
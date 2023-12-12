# product sharing functions
def product_sharing(lofo_dict, graph, all_reactants_dict, all_products_dict, mets_to_ignore=[], s_func=min, pp_func=sum, return_pushin_p=False, ignore_zero_sharing=False):
    data={}
    num_reactants=[]
    num_products=[]
    pushin_p={}
    product_push_data, product_push_func=product_push_main(fba_ecoli_lofo, ecoli_model, ecoli_dg, pp_func, mets_to_ignore=mets_to_ignore, f1_rev=rev1, f2_rev=rev2, ignore_zero_push=ignore_zero_push, return_func_df=False)
    pushin_p=product_push_data
    for rxn, obj in zip(fba_ecoli_lofo.keys(), fba_ecoli_lofo.values()):
        loop_data=[]
        rxn_df=ecoli_model.reactions.get_by_id(rxn)
        reversibility=rxn_df.reversibility
        rxn_reactants=rxn_df.reactants
        rxn_products=rxn_df.products
        rxn_reactant_names=[i.id for i in rxn_reactants]
        rxn_product_names=[i.id for i in rxn_products]

        # get num reactants and products for each reaction
        num_reactants.append(len(rxn_reactant_names))
        num_products.append(len(rxn_product_names))

        # should be using all reactants or products - should we split by essentiality like below?
        r_sharing=[]
        p_sharing=[]
        for r in rxn_reactant_names:
            r_sharing.append(all_reactants_dict[r])

        for p in rxn_product_names:
            p_sharing.append(all_products_dict[p]) 

        # take min or max? or sum?
        loop_data.append(func(r_sharing)) # min, max, or sum
        loop_data.append(func(p_sharing)) # min, max, or sum

        # append len of each
        loop_data.append(len(r_sharing))
        loop_data.append(len(p_sharing))

        # append reversibility
        loop_data.append(reversibility)

        # taking len() function as a metric for reaction sharing does not work
        r_sharing=[len(i.reactions) for i in rxn_reactants]
        p_sharing=[len(i.reactions) for i in rxn_products]
        loop_data.append(func(r_sharing)) # min, max, or sum
        loop_data.append(func(p_sharing)) # min, max, or sum
        

        if (obj < 2):
            loop_data.append(1)
        else:
            loop_data.append(0)


        # add final data
        data[rxn]=loop_data
        
    if (return_pushin_p):
        return data, pushin_p
    return data # pushin_p

# run it
def get_sharing_data(fba_ecoli_lofo, ecoli_dg, all_reactants_dict, all_products_dict, mets_to_ignore, func=min, pp_func=sum, ignore_mets=True, rev1=True, rev2=True, return_pushin_p=False, ignore_zero_push=False):
    data={}
    num_reactants=[]
    num_products=[]
    pushin_p={}
#     product_push_data, product_push_func=product_push_main(fba_ecoli_lofo, ecoli_model, ecoli_dg, pp_func, mets_to_ignore=mets_to_ignore, f1_rev=rev1, f2_rev=rev2, ignore_zero_push=ignore_zero_push, return_func_df=False)
#     pushin_p=product_push_data
    for rxn, obj in zip(fba_ecoli_lofo.keys(), fba_ecoli_lofo.values()):
        loop_data=[]
        rxn_df=ecoli_model.reactions.get_by_id(rxn)
        reversibility=rxn_df.reversibility
        rxn_reactants=rxn_df.reactants
        rxn_products=rxn_df.products
        rxn_reactant_names=[i.id for i in rxn_reactants]
        rxn_product_names=[i.id for i in rxn_products]

        # get num reactants and products for each reaction
        num_reactants.append(len(rxn_reactant_names))
        num_products.append(len(rxn_product_names))

        # should be using all reactants or products - should we split by essentiality like below?
        r_sharing=[]
        p_sharing=[]
        for r in rxn_reactant_names:
            r_sharing.append(all_reactants_dict[r])

        for p in rxn_product_names:
            p_sharing.append(all_products_dict[p]) 

        # take min or max? or sum?
        loop_data.append(func(r_sharing)) # min, max, or sum
        loop_data.append(func(p_sharing)) # min, max, or sum

        # append len of each
        loop_data.append(len(r_sharing))
        loop_data.append(len(p_sharing))

        # append reversibility
        loop_data.append(reversibility)

        # taking len() function as a metric for reaction sharing does not work
        r_sharing=[len(i.reactions) for i in rxn_reactants]
        p_sharing=[len(i.reactions) for i in rxn_products]
        loop_data.append(func(r_sharing)) # min, max, or sum
        loop_data.append(func(p_sharing)) # min, max, or sum
        
        ##### START PRODUCT PUSH CODE
        # see how products of selected rxn are used as reactants in next rxn
        # a.index(min(a))
        node_out_neighbors=[n for n in ecoli_dg.neighbors(rxn)]
        product_pass_dict={}
        nxt_react=[]

        if (rev1):
            if (reversibility): # if the reaction is reversible, use both the products and the reactants to feed into the next reaction
                mets_to_consider=rxn_product_names+rxn_reactant_names
            else:
                mets_to_consider=rxn_product_names
        else:
            mets_to_consider=rxn_product_names

        for p in mets_to_consider:
            p_share_count=0
            
            # NEW LOCATION
#             if (ignore_mets):
#                 if p in mets_to_ignore.keys():
#                     continue

            for next_rxn in node_out_neighbors:
                next_rxn_df=ecoli_model.reactions.get_by_id(next_rxn)
                next_reversibility=next_rxn_df.reversibility
                next_rxn_reactants=next_rxn_df.reactants
                next_rxn_products=next_rxn_df.products
                next_rxn_reactant_names=[i.id for i in next_rxn_reactants]
                next_rxn_product_names=[i.id for i in next_rxn_products]
                nxt_react+=next_rxn_product_names
                if (len(nxt_react)==0):
                    print('no next products')

                # second reversibility clause
                if (rev2):
                    if (next_reversibility):
                        next_mets_to_consider=next_rxn_reactant_names+next_rxn_product_names
                    else:
                        next_mets_to_consider=next_rxn_reactant_names
                else:
                    next_mets_to_consider=next_rxn_reactant_names

                # ORIGINAL LOCATION
                if (ignore_mets):
                    if p in mets_to_ignore.keys():
                        continue
                        
                if p in next_mets_to_consider:
                    p_share_count+=1

            # NEW
            if (ignore_zero_push):
                if (p_share_count==0):
                    continue
                else:
                    product_pass_dict[p]=p_share_count
            else:
                product_pass_dict[p]=p_share_count

        pushin_p[rxn]=product_pass_dict
        
        # get total product push forward
        product_push=pp_func(list(product_pass_dict.values())) # min, max, or sum
        ##### END PRODUCT PUSH CODE
        
        # NEW LINE ADDED
#         product_push=product_push_func[rxn]
        
        # append
        loop_data.append(product_push)

        if (obj < 2):
    #         for r in rxn_reactant_names:
    #             r_sharing.append(all_reactants_ess_dict[r])
    #         for p in rxn_products_names:
    #             p_sharing.append(all_products_ess_dict[p])

            loop_data.append(1)
        else:
    #         for r in rxn_reactant_names:
    #             r_sharing.append(all_reactants_noness_dict[r])
    #         for p in rxn_products_names:
    #             p_sharing.append(all_products_noness_dict[p])

            loop_data.append(0)


        # add final data
        data[rxn]=loop_data
        
    if (return_pushin_p):
        return data, pushin_p
    return data # pushin_p
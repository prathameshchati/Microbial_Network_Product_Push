def dict_isolate(dic, val):
    match_keys=[]
    for key, value in dic.items():
        if (value==val):
            match_keys.append(key)
    return match_keys

def df_isolate(df, col, val):
    return df[df[col]==val]

def log_columns(df, cols, base=2):
    for col in cols:
        if (base==2):
            df[f"log_{col}"]=np.log(df[col]+1)
        elif (base==10):
            df[f"log_{col}"]=np.log10(df[col]+1)
        else:
            print("Choose valid base (2 or 10)")
    return df

from sklearn import preprocessing
def normalize(data, lower, upper):
    data_norm=[lower + (upper - lower) * x for x in data]
    return data_norm

def adj_columns(df, cols, prefix):
    for col in cols:
        d=np.array(df[col])
        df[f"{prefix}_{col}"]=preprocessing.normalize([d]).tolist()[0]
    return df

def z_score(df, cols, prefix):
    for col in cols:
        d=np.array(df[col])
        df[f"{prefix}_{col}"]=stats.zscore(d).tolist()
    return df

def scale(df, cols, lb=0, ub=1):
    for col in cols:
        df[f"scale_{col}"]=scale_list(list(df[col]), lb, ub)
    return df

def scale_number(unscaled, to_min, to_max, from_min, from_max):
    return (to_max-to_min)*(unscaled-from_min)/(from_max-from_min)+to_min

def scale_list(l, to_min, to_max):
    return [scale_number(i, to_min, to_max, min(l), max(l)) for i in l]

def replace_index(df, replace_pairs):
    for r1, r2 in replace_pairs:
        df.index=df.index.str.replace(r1,r2)
    return df

def replace_column(df, replace_pairs):
    for r1, r2 in replace_pairs:
        df.columns=df.columns.str.replace(r1,r2)
    return df

def sort_dict(d):
    return {k: v for k, v in sorted(d.items(), key=lambda item: item[1])}

def sorted_counter(data):
    return sort_dict(Counter(data))

def sharing_data_processing(data_df, base=2):
    # format
    data_df=pd.DataFrame.from_dict(data_df).T
    data_df.columns=["r", "p", "num_r", "num_p", "rev", "r_all", "p_all", "p_push", "ess"]

    # merge
    data_df=pd.concat([data_df, cent_metrics_df], axis=1) # was degree_df
    
    # modify columns 
    mod_cols=["r", "p", "k_in", "k_out", "k", "r_all", "p_all", "p_push"]
    data_df[mod_cols]=data_df[mod_cols].astype(float)

    # log
    data_df=log_columns(data_df, mod_cols, base=base)

    # norm
    data_df=adj_columns(data_df, mod_cols, "norm")

    # norm
    data_df=z_score(data_df, mod_cols, "z")

    # scale
    data_df=scale(data_df, mod_cols)
    
    # reaction degree
    data_df_sources=data_df.loc[sources.index]
    data_df_sinks=data_df.loc[sinks.index]
    
    # reaction essentiality
    data_df_ess=data_df[data_df["ess"]==1]
    data_df_noness=data_df[data_df["ess"]==0]
    
    data_df=data_df.sort_values(by="r", ascending=False)
    
#     data_df["log_p_push_k_out"]=data_df["log_p_push"]/data_df["log_k_out"]

    return data_df, data_df_sources, data_df_sinks, data_df_ess, data_df_noness

def sharing_data_processing_2(data_df, data_df_sources, data_df_sinks, rm_sources=True, rm_sinks=True, base=2):

    # remove source and sink nodes
    if (rm_sources):
        data_df_f=data_df.drop(data_df_sources.index)
    else:
        data_df_f=data_df.copy()
        
    if (rm_sinks):
        data_df_f=data_df_f.drop(data_df_sinks.index)
        
    data_df_f["pp_kout"]=data_df_f["p_push"]/data_df_f["k_out"]
    data_df_f=log_columns(data_df_f, ["pp_kout"], base=base)

    # get ess and non ess 
    data_df_f_ess=data_df_f[data_df_f["ess"]==1]
    data_df_f_noness=data_df_f[data_df_f["ess"]==0]
    
    return data_df_f, data_df_f_ess, data_df_f_noness
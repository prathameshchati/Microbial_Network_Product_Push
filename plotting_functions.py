# assign bins for histogram
def assign_bins(df, col, bins, manual_bins=[], remove_outliers=False, q=0.99, z_filter=False, z_t=2, z_col='z_r', return_hist=False):
    
    # remove outliers (https://stackoverflow.com/questions/23199796/detect-and-exclude-outliers-in-a-pandas-dataframe)
    if (remove_outliers):
        q = df[col].quantile(q)
        df=df[df[col] < q]
#         df=adj_columns(df, ["r", "p", "k_in", "k_out", "k"], "norm") # RENORMALIZE

    if (z_filter):
        df=df[df[z_col]<z_t]
        
    col_data=list(df[col])
    
    # make sure min of data is not zero
    if (min(col_data)==0):
        print("Warning, minimum of data is zero.")
        
    if (len(manual_bins)==0):
        hist,bin_edges=np.histogram(col_data, bins=bins)
    else:
        hist,bin_edges=np.histogram(col_data, bins=manual_bins)
        
    hist=list(hist)
    bin_edges=list(bin_edges)
    
    # not needed
#     bin_edges.insert(0, 0)
#     max_p1=max(col_data)+1
#     bin_edges.append(max_p1)
    
    # get bin edge pairs
    bin_edge_pairs={}
    for i in range(len(bin_edges)):
        if (i>=len(bin_edges)-1):
            break
        pair=(bin_edges[i], bin_edges[i+1])
        bin_edge_pairs[pair]=f"g{i}"        
        
    # bound the data
    df_hist=pd.DataFrame()
    for bin_pair, label in zip(bin_edge_pairs.keys(), bin_edge_pairs.values()):
        df_bounded=df.loc[df[col].between(bin_pair[0], bin_pair[1])]
        df_bounded["hist_label"]=label
        df_hist=pd.concat([df_hist, df_bounded])
        
    if (return_hist):
        return df_hist, hist, bin_edges, bin_edge_pairs
    else:
        return df_hist
    
# computes fraction of essential reactions
def compute_fraction_essential_hist(df_hist, hist, bin_edges, bin_edge_pairs, compute_ratio=False):
    df_ess=pd.DataFrame()
    df_ess["hist_label"]=list(bin_edge_pairs.values())
    df_ess["hist_counts"]=hist
    ess_counts=[]
    noness_counts=[]
    ess_fracs=[]
    bin_centers=[]
    left_bin=[]
    right_bin=[]
    for pair, label in zip(bin_edge_pairs.keys(), bin_edge_pairs.values()):
        bin_centers.append(np.mean(pair))
        left_bin.append(pair[0])
        right_bin.append(pair[1])
        df_hist_label=df_hist[df_hist["hist_label"]==label]
        num_ess=len(df_hist_label[df_hist_label["ess"]==1])
        num_noness=len(df_hist_label[df_hist_label["ess"]==0])
        if (compute_ratio):
            if (num_noness==0):
                frac_ess=0
            else:
                frac_ess=num_ess/num_noness
        else:
            if ((num_ess+num_noness)==0):
                frac_ess=0
            else:
                frac_ess=num_ess/(num_ess+num_noness)
        
        # append data
        ess_counts.append(num_ess)
        noness_counts.append(num_noness)
        ess_fracs.append(frac_ess)
        
    # update df
    df_ess["num_ess"]=ess_counts
    df_ess["num_noness"]=noness_counts
    df_ess["frac_ess"]=ess_fracs
    df_ess["left_bin"]=left_bin
    df_ess["bin_centers"]=bin_centers
    df_ess["right_bin"]=right_bin
    
    # return 
    return df_ess

def create_bins(lb, ub, s, mb=1):
    if (mb==1):
        manual_bins=list(np.arange(lb, ub, s))
        print(list(np.arange(lb, ub, s)))
        print(len(list(np.arange(lb, ub, s))))
    else: 
        manual_bins=[]
    return manual_bins

def abline(slope, intercept):
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--', color='red')

def sk_lr(df, x_label, y_label):
    x=np.array(df[x_label]).reshape((-1, 1))
    y=np.array(df[y_label])
    model = LinearRegression()
    model.fit(x, y.reshape((-1, 1)))
    r_sq=model.score(x, y)
    print(f"coefficient of determination: {r_sq}")
    
def sp_lr(df, x_label, y_label):
    x=np.array(df[x_label]) # .reshape((-1, 1))
    y=np.array(df[y_label])
    slope, intercept, r_value, p_value, std_err = linregress(y, x)
    print("Slope: ", slope)
    print("Intercept: ", intercept)
    print("r-value: ", r_value)
    print("p-value: ", p_value)
    print("standard error: ", std_err)
    return slope, intercept, r_value, p_value, std_err
    
# plot histogram 
def plot_histogram(data_df_hist, data_df_ess_frac, xl, yl, subset=[], frac_ess=True):
    if (len(subset)!=0):
        data_df_hist=data_df_hist[data_df_hist["hist_label"].isin(subset)]
        data_df_ess_frac=data_df_ess_frac[data_df_ess_frac["hist_label"].isin(subset)]
    if (frac_ess):
        df_to_use=data_df_ess_frac.copy()
    else:
        df_to_use=data_df_hist.copy()
    x=df_to_use[xl]
    y=df_to_use[yl]
    sk_lr(df_to_use, xl, yl)
    plt.scatter(x, y)  
    plt.xlabel(xl)
    plt.ylabel(yl)
#     plt.show()
    
def group_ranges(b1, b2):
    subset=[]
    for i in range(b1,b2+1):
        subset.append(f"g{str(i)}")
    return subset

# def get_mean()
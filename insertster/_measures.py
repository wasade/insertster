def make_f_beta(beta):
    """Create a f beta function

    Parameters
    ----------
    beta : float
        The beta to use where a beta of 1 is the f1-score or F-measure

    Returns
    -------
    function
        A function to compute the f_beta score
    """
    beta_2 = beta**2
    coeff = (1 + beta_2)
    def f(global_, local_, node):
        """Compute the f-measure

        Parameters
        ----------
        global_ : np.array
            All of the scores for a given query
        local_ : np.array
            The scores for the query at the current node
        node : skbio.TreeNode
            The current node being evaluated
        """
        p = len(global_) / len(local_)
        r = len(local_) / node.ntips

        return coeff * (p * r) / ((beta_2 * p) + r)
    return f


f1_measure = make_f_beta(1.0)
fhalf_measure = make_f_beta(0.5)
f2_measure = make_f_beta(2.0)


available_score_functions = {'f1-measure': f1_measure,
                             'fhalf_measure': fhalf_measure,
                             'f2_measure': f2_measure}

def parse_newick(newick_string):
    """
    Parse a Newick string into a dictionary representation.
    
    Returns:
        dict: Tree structure where keys are node IDs and values are lists of children
        dict: Node metadata (names, branch lengths, etc.)
        int: Root node ID
    """
    
    # Remove whitespace and trailing semicolon
    newick = newick_string.strip().rstrip(';')
    
    # Global variables for the parser
    node_counter = 0
    node_children = {}
    node_label = {}
    branch_length = {}
    
    def get_next_node_id():
        nonlocal node_counter
        node_id = node_counter
        node_counter += 1
        return node_id
    
    def parse_node(s, pos=0):
        """
        Recursively parse a node from position 'pos' in string 's'.
        
        Returns:
            tuple: (node_id, new_position)
        """
        nonlocal node_children, node_label, branch_length
        
        current_node_id = get_next_node_id()
        node_children[current_node_id] = []  # Initialize children list
        node_label[current_node_id] = None
        branch_length[current_node_id] = None
        
        # Skip whitespace
        while pos < len(s) and s[pos].isspace(): pos += 1
        
        # Case 1: Internal node - starts with '('
        if pos < len(s) and s[pos] == '(':
            pos += 1  # Skip opening parenthesis
            
            # Parse children until we hit the closing parenthesis
            while pos < len(s) and s[pos] != ')':
                # Skip whitespace and commas
                while pos < len(s) and s[pos] in ' \t,': pos += 1
                
                if pos < len(s) and s[pos] != ')':
                    # Recursively parse child
                    child_id, pos = parse_node(s, pos)
                    node_children[current_node_id].append(child_id)
            
            if pos < len(s) and s[pos] == ')': pos += 1  # Skip closing parenthesis
        
        # Case 2 & 3: Parse node label and branch length (for both leaf and internal nodes)
        # Node label comes first
        label_start = pos
        while pos < len(s) and s[pos] not in '[:,();': pos += 1
        
        if pos > label_start:
            node_label[current_node_id] = s[label_start:pos].strip()
        
        # Branch length comes after ':'
        if pos < len(s) and s[pos] == ':':
            pos += 1  # Skip ':'
            
            # Parse branch length
            length_start = pos
            while pos < len(s) and s[pos] not in '[,();': pos += 1
            
            if pos > length_start:
                try:
                    branch_length[current_node_id] = float(s[length_start:pos])
                except:
                    raise ValueError(f"{length_start:pos} not a valid branch length")
                    #pass  # Invalid branch length, keep as None
        
        #progress past anything included beyond branch length
        while pos < len(s) and s[pos] not in ',();': pos += 1
        
        return current_node_id, pos
    
    # Start parsing from the beginning
    root_id, _ = parse_node(newick, 0)
    
    #now figure out which nodes are leaves and which are parents
    n_nodes = get_next_node_id()
    parents = []
    leaves = []
    for id in range(n_nodes):
        if node_children[id] == []: leaves.append(id)
        else: parents.append(id)
    
    return node_children, node_label, branch_length, leaves, parents, root_id


from sticcs import sticcs


# another function to convert a newick string to a sticcs tree.
# This object uses numeric leaf IDs from 0 to N-1.
# These can be passed to the function as a dictionary.
# Otherwise the function attempts to convert leaf IDs to integers
# Internal nodes don't need the same level of consistency, so those will
# be kept in the order i which they were found, but renumbered so they come after the leaf indices
def newick_to_sticcs_Tree(newick_string, leaf_idx_dict=None, interval=None):
    
    #first parse newick and number nodes and leaves by order of finding them
    node_children, node_label, branch_length, leaves, parents, root_id = parse_newick(newick_string)
    
    n_leaves = len(leaves)
    
    #Dictionary to convert leaf number to new index starting from zero
    if leaf_idx_dict is not None:
        new_leaf_IDs = leaf_idx_dict.values()
    
    else:
            #assume leaves are integers
        try:
            leaf_labels = [node_label[leaf_number] for leaf_number in leaves]
            new_leaf_IDs = [int(label) for label in leaf_labels]
            leaf_idx_dict = dict(zip(leaf_labels, new_leaf_IDs))
        
        except: raise ValueError("Leaf labels are not numbers. Please provide a leaf_idx_dict to link leaf labels to a numeric value")
    
    assert len(new_leaf_IDs) == n_leaves and max(new_leaf_IDs) == n_leaves-1, "Please provide a leaf_idx_dict with consecutive indices from 0 to n_leaves-1"
    
    #now assign the new ID to each leaf number (remember leaf numbers in the parsed newick are simply in the order they were encountered)
    new_node_idx = dict(zip(leaves, [leaf_idx_dict[node_label[leaf_number]] for leaf_number in leaves]))
    
    #the new parents are just n_leaves along from where they were
    new_parents = [id+n_leaves for id in parents]
    
    new_node_idx.update(dict(zip(parents, new_parents)))
    
    #now update objects for the new sticcs Tree object
    new_node_children = dict()
    for item in node_children.items():
        new_node_children[new_node_idx[item[0]]] = [new_node_idx[x] for x in item[1]]
    
    new_node_parent = dict()
    for item in new_node_children.items():
        for child in item[1]:
            new_node_parent[child] = item[0]
    
    new_node_label = dict()
    for item in node_label.items():
        new_node_label[new_node_idx[item[0]]] = item[1]
    
    objects = {"leaves":sorted(new_leaf_IDs), "root":n_leaves, "parents":new_parents,
               "node_children": new_node_children, "node_parent":new_node_parent,
               "interval":interval}
    
    return (sticcs.Tree(objects=objects), new_node_label)


def sim_test_newick_parser(n=12, reps=5):
    import msprime
    from twisst2 import TopologySummary
    
    for i in range(reps):
        ts = msprime.sim_ancestry(n, ploidy=1)
        t = ts.first()
        sim_tree_summary = TopologySummary.TopologySummary(t)
        sim_tree_ID = sim_tree_summary.get_topology_ID()
        sim_tree_newick = t.as_newick(include_branch_lengths=False)
        print("Simulated tree:", sim_tree_newick)
        
        #now read it
        parsed_tree, node_labels = newick_to_sticcs_Tree(sim_tree_newick, leaf_idx_dict= dict([("n"+str(i), i) for i in range(n)]))
        
        parsed_tree_summary = TopologySummary.TopologySummary(parsed_tree)
        parsed_tree_ID = parsed_tree_summary.get_topology_ID()
        parsed_tree_ID == sim_tree_ID
        print("Parsed tree:   ", parsed_tree.as_newick(node_labels=node_labels))
        
        print("Match:", parsed_tree_ID == sim_tree_ID)




# Test the parser
if __name__ == "__main__":
    
    #Manual test
    #newick_string = '((n1,n4),((n0,n3),(n2,n5)));'
    #parsed_tree, node_labels = newick_to_sticcs_Tree(newick_string, leaf_idx_dict={"n0":0, "n1":1, "n2":2, "n3":3, "n4":4, "n5": 5})
    #print(parsed_tree.as_newick(node_labels=node_labels))
    
    #Auto test multiple times
    sim_test_newick_parser()

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
def newick_to_sticcs_Tree(newick_string, leaf_idx_dict=None, allow_additional_leaves=False, interval=None):
    
    #first parse newick and number nodes and leaves by order of finding them
    node_children, node_label, branch_length, leaves, parents, root_id = parse_newick(newick_string)
    #nodes that come out of parse_newick are integers, but they are in the order they were encountered in reading the string, so they are meaningless
    #What matters is the node_labels, as these are what those nodes were called (i.e. leaf names, but internal nodes can theoretically have labels too)
    #Now the sticcs tree needs numeric leaves, but these need to start from zero - so they are NOT the same as the nodes that come from parse_newick()
    #Instead each node_label needs to map to an integer between zero and n_leaves
    #This can be provided in the leaf_idx_dict
    #If not, all node_labels need to be integers in the newick string
    
    n_leaves = len(leaves)
    
    #get the label for each leaf
    leaf_labels = [node_label[leaf_number] for leaf_number in leaves]
    assert len(set(leaf_labels)) == n_leaves, "Leaf labels must all be unique."
    
    #now, if there is no dictionary provided, we assume leaves are already numbered 0:n-1
    if not leaf_idx_dict:
        if not set(leaf_labels) == set([str(i) for i in range(n_leaves)]):
            raise ValueError("Leaf labels not consecutive integers. Please provide a leaf_idx_dict with consecutive indices from 0 to n_leaves-1")
        
        _leaf_idx_dict_ = dict([(label, int(label)) for label in leaf_labels])
    
    else:
        #We already have the dictionary to convert leaf number to new index
        _leaf_idx_dict_ = {}
        _leaf_idx_dict_.update(leaf_idx_dict)
        assert set(_leaf_idx_dict_.values()) == set(range(len(_leaf_idx_dict_))), "leaf_idx_dict must provide consecutive indices from zero"
        #Check that all leaves are represented
        unrepresented_leaves = [label for label in leaf_labels if label not in leaf_idx_dict]
        if len(unrepresented_leaves) > 0:
            if allow_additional_leaves:
                i = len(_leaf_idx_dict_)
                for leaf_label in sorted(unrepresented_leaves):
                    _leaf_idx_dict_[leaf_label] = i
                    i += 1
            else:
                raise ValueError("Some leaves are not listed in leaf_idx_dict. Set allow_additional_leaves=True if you want to risk it.")
    
    new_leaf_IDs = sorted(_leaf_idx_dict_.values())
    
    #now assign the new ID to each leaf number (remember leaf numbers in the parsed newick are simply in the order they were encountered)
    new_node_idx = dict(zip(leaves, [_leaf_idx_dict_[node_label[leaf_number]] for leaf_number in leaves]))
    
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
    
    objects = {"leaves":new_leaf_IDs, "root":n_leaves, "parents":new_parents,
               "node_children": new_node_children, "node_parent":new_node_parent,
               "interval":interval}
    
    return (sticcs.Tree(objects=objects), new_node_label)



# a class that just holds a bunch of trees, but has a few attributes and methods that resemble a tskit treesequence object
# This allos us to import newick trees and treat the list as a treesequence for a few things (like my quartet distance calculation)
class TreeList:
    def __init__(self, trees):
        self.num_trees = len(trees)
        self.tree_list = trees
    
    def trees(self):
        for tree in self.tree_list:
            yield tree


def parse_newick_file(newickfile, leaf_names=None):
    
    #if leaves are not integers, need to link each leaf to its index
    leaf_idx_dict = dict(zip(leaf_names, range(len(leaf_names)))) if leaf_names is not None else None
    
    trees = []
    
    i = 1
    
    for line in newickfile:
        tree, node_labels = newick_to_sticcs_Tree(line.strip(), leaf_idx_dict=leaf_idx_dict, allow_additional_leaves=True, interval = (i, i))
        trees.append(tree)
        i+=1
    
    return TreeList(trees)


def parse_argweaver_smc(argfile, leaf_names=None):
    #annoyingly, the trees have numeric leaves, but these are NOT in the order of the haps
    #The order if given by the first line, so we can figure out which hap each leaf points to
    
    #get the order of haps from argweaver output
    #The position of each number in this list links it to a leaf (leaves are numeric) 
    leaf_names_reordered = argfile.readline().split()[1:]
    
    n_leaves = len(leaf_names_reordered)
    
    _leaf_names = leaf_names if leaf_names is not None else [str(i) for i in range(n_leaves)]
    
    #link each leaf name to its real index. This is the index that will be used for topology weighting
    leaf_name_to_idx_dict = dict([(_leaf_names[i], i) for i in range(n_leaves)])
    
    #link the leaf number in the current file (which is arbitrary) to its real index
    leaf_idx_dict = dict([(str(i), leaf_name_to_idx_dict[leaf_names_reordered[i]]) for i in range(n_leaves)])
    
    trees = []
    
    chrom, chrom_start, chrom_len = argfile.readline().split()[1:]
    for line in argfile:
        if line.startswith("TREE"):
            elements = line.split()
            interval=(int(elements[1]), int(elements[2]),)
            tree_newick = elements[3]
            tree, node_labels = newick_to_sticcs_Tree(tree_newick, leaf_idx_dict=leaf_idx_dict, allow_additional_leaves=True, interval=interval)
            trees.append(tree)
    
    return TreeList(trees)



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
    
    #newick_string = '((1,4),((0,3),(2,5)));'
    #parsed_tree, node_labels = newick_to_sticcs_Tree(newick_string)
    #print(parsed_tree.as_newick())
    
    #Auto test multiple times
    sim_test_newick_parser()


#! /usr/bin/env python3

from anytree import AnyNode, RenderTree, PostOrderIter, LevelOrderGroupIter
from anytree.render import AsciiStyle
from xml.dom.minidom import parse
import sys
import re

#
# This is aweful spaghetti-code and I am sorry for this!!
#

def check_position_coverage(poly, data, all_parent_positions=[]):
    def return_uncovered():
        return 0,0,0
    # ignore insertions (maybe do that later...)
    if '.' in poly:
        return return_uncovered()
    #haplogroup-info
    pos = re.search('[0-9]+', poly).group()
    base = re.search('[A-Zd]', poly).group()
    # pileup-info
    try:
        pile = data[pos] # e.g. CCTC
        cov = len(pile) # coverage
        if cov == 0:
            return return_uncovered()
        target = pile.count(base) # on target
        perc = round((target/cov)*100, 2)
        if poly.endswith("!"): #a remutation
            # check, if the mutation is already represented in the branch leading here
            # e.g. 16311T in L3 --> 16311T! in several U subgroups
            if pos+base in all_parent_positions:
                return return_uncovered()

    except KeyError:
        return return_uncovered()
    ##
    ## compare the poly with the mpileup_data to return coverage
    ##
    ##
    return perc, target, cov

prefix = sys.argv[3]

# open XML file
with open(sys.argv[1]) as xml_file:
    tree = parse(xml_file)

# open pileup file
with open(sys.argv[2]) as pileup_file:
    pileup_data = {}
    for line in pileup_file:
        cols = line.split('\t')
        sequence = cols[4]
        quality = cols[5]
        good_bases = ''.join([base for base, q in zip(sequence, quality) if q != '!'])
        # import position -> bases
        pileup_data[cols[1]] = good_bases.upper()

# create the first tree
node = AnyNode(id='mtMRCA', parent=None, data=[])
name_node_dict = {'mtMRCA':node}
node_parent_polys = {'mtMRCA':set()}


# walk through the XML and fill the tree
for haplogroup in tree.getElementsByTagName('haplogroup'):
    name = haplogroup.getAttribute('name')
    
    # get the haplogroup name from the XML parent
    parent_node = haplogroup.parentNode
    if parent_node.nodeType == parent_node.ELEMENT_NODE and parent_node.tagName == 'haplogroup':
        parent = parent_node.getAttribute('name')
    else:
        continue

    # overwrite parent_node with the anytree Node element
    parent_node = name_node_dict[parent]
    
    # get all positions that are checked in the same branch
    # this is to double check if a position was aleady found within the same branch
    # e.g L3: 16311T --> U1b3: 16311T!
    all_parent_positions = node_parent_polys[parent_node.id]

    data = []
    found_covered = 0
    found_target = 0
    for child in haplogroup.childNodes: #child here means XML childs --> for parsing the POLYs
        if child.nodeType == child.ELEMENT_NODE and child.tagName == 'details':
        # Get data for each haplogroup-defining position 
        # Get the 'poly' elements directly under the 'details' element
            polys = []
            raw_polys = []
            for poly in child.getElementsByTagName('poly'):
                # extract position from XML
                poly_data = poly.firstChild.data
                raw_polys.append(poly_data)
                # now parse the positions and calculate cpverage
                perc, target, cov = check_position_coverage(poly_data, pileup_data, all_parent_positions)
                poly_string = f"{poly_data} ({perc:.2f}% {target}/{cov})"
                if cov > 0:
                    found_covered += 1
                if perc > 10:
                    found_target += 1
                # save for later 
                polys.append(poly_string)
            data.append(polys)
            # now update the all_parent_positions to store them in the dict for all descendent haplogroups
            tmp = all_parent_positions.copy()
            tmp.update(raw_polys)
            node_parent_polys.update({name: tmp})
    #add node-data to the tree
    data.append([found_target, found_covered])
    # and add node-data for later overwriting of the stats for each node (including leaves)
    data.append([found_target, found_covered])
    # add penalty-value for later updating
    data.append(1)
    #create the tree
    tmp = AnyNode(id=name, parent=parent_node, data=data)
    name_node_dict.update({name:tmp})   

# now summarize stats on each node
for hap in PostOrderIter(node):
    #walk from the leaves up and count the number of successful (leave) haplogroups on each node:
    #thats what the PostOrderIter does
    total_target = sum(x.data[2][0] for x in hap.children)
    total_cov = sum(x.data[2][1] for x in hap.children)
    if total_cov == 0:
        continue
    if len(hap.data) < 3:
        hap.data=[total_target, total_cov]
        continue # back to root
    # introduce a penalty for nodes in a branch that are
    # covered, but dont have the target allel
    #   Tree                       Node   Position  Penalty
    #  +-- L1c3b'c                 1/20   0/1       4 (Distance to Tip)
    #      |-- L1c3b               1/22   0/2       3 (Distance to Tip)
    #      |   |-- L1c3b1          1/4    0/3       2 (Distance to Tip)
    #      |   |   |-- L1c3b1a     1/2    1/2       1 (Tip)
    #      |   |   +-- L1c3b1b     0/2    0/2       1 (Tip)
    #      |   +-- L1c3b2          0/7    0/7       1 (Tip)
    #      +-- L1c3c               0/9    0/9       1 (Tip)
    #
    # So get the minimum Penalty from the children (if it has node-coverage)
    try:
        child_penalty = min(x.data[3] for x in hap.children if x.data[2][0] != 0)
    except:
        child_penalty = 1 # no counts on the nodes will be filtered out anyways later...
    if  hap.data[1][0] == 0:
        penalty = child_penalty + 1
    else:
        penalty = 1

    hap.data[3] = penalty
    hap.data[2] = [total_target, total_cov]

def print_header(file):
    print('PhyloTree\tNodeCoverage\tPositionCoverage\tPenalty\tReadCoverage', file=file)

def print_line(row, file):
    print('\t'.join(
            [
                f"{row.pre.rstrip()} {row.node.id}",
                f"{row.node.data[2][0]}/{row.node.data[2][1]}",
                f"{row.node.data[1][0]}/{row.node.data[1][1]}",
                f"{row.node.data[3]}",
                f"{'; '.join(row.node.data[0])}"
            ]
        ), file=file
    )

# and print summary files
rendered = list(RenderTree(node, style=AsciiStyle))

# First, print stats for every haplogroup!
with open(f"{prefix}.01_stats_all_groups.tsv", 'w') as tree1:
    print_header(tree1)
    for row in rendered:
        if len(row.node.data) > 2: #any node in the middle
            print_line(row, tree1)
        else: # root
            print(f"{row.pre.rstrip()} {row.node.id}\t{row.node.data[0]}/{row.node.data[1]}\t\t", file=tree1)

#Then, print the BEST path
# create a new tree that only consists of the best nodes
best_node = AnyNode(id='mtMRCA', parent=None, data=node.data)

def get_best_child(node, best_parent):
    #The penalty is _really_ strong e.g. Penalty of 2 -> cut node-count in half 
    order = sorted(node.children, key=lambda x: x.data[2][0]/x.data[3])
    try:
        best_child = order[-1]
        tmp = AnyNode(id=best_child.id, parent=best_parent, data=best_child.data)
        return get_best_child(order[-1], tmp)
    except IndexError: #reached the tip
        pass

get_best_child(node, best_node)

best_path = list(RenderTree(best_node, style=AsciiStyle))

with open(f"{prefix}.03_stats_best_path.tsv", 'w') as tree3:    
    print_header(tree3)
    for row in best_path:
        try:
            if row.node.data[2][0] == 0:
                continue #the last tip is random and could be empty...
            print_line(row, tree3)
        except IndexError: # root
            print(f"{row.pre.rstrip()} {row.node.id}\t{row.node.data[0]}/{row.node.data[1]}\t\t", file=tree3)

#Then, prune the tree, so that:
# - nodes with 0 coverage are removed
# - nodes with penalty 3 or more are removed
# 
# iterate through the tree and set all parents to None if these conditions apply

for hap in PostOrderIter(node):
    #walk from the leaves up and count the number of successful (leave) haplogroups on each node:
    #thats what the PostOrderIter does
    if len(hap.data) < 3:
        continue # root

    total_coverage = hap.data[2][0]
    penalty = hap.data[3]

    if total_coverage == 0 or penalty >= 4:
        hap.parent = None

prune_rendered = list(RenderTree(node, style=AsciiStyle))
with open(f"{prefix}.02_stats_all_groups_with_coverage.tsv", 'w') as tree2:    
    print_header(tree2)
    for row in prune_rendered:
        try:
            test = row.node.data[2][0]
            print_line(row, tree2)
        except IndexError: # root
            print(f"{row.pre.rstrip()} {row.node.id}\t{row.node.data[0]}/{row.node.data[1]}\t\t", file=tree2)

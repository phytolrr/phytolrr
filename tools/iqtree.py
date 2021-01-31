# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019 phytolrr.com. All rights reserved.

import re
import logging


class Atrribute(object):
    def __init__(self):
        self.label = None
        self.color = None


class Node(object):
    def __init__(self):
        self.bootstrap_support = None
        self.attrs = Atrribute()

    def add_attr(self, key, value):
        if key == 'label':
            self.attrs.label = value
        elif key == 'color':
            self.attrs.color = value
        else:
            logging.error("Unexpected attr " + key)


class Leaf(Node):
    def __init__(self):
        super().__init__()
        self.seq_id = None


class Container(Node):
    def __init__(self):
        super().__init__()
        self.children = []
        self.nodes = []

    def add_child(self, container):
        self.children.append(container)
        return container

    def add_node(self, node):
        self.nodes.append(node)
        return node


class Tree(object):
    def __init__(self, name):
        self.name = name
        self.root_node = None


SEQ_ID_RE = re.compile(r"'?[\w\.-]+'?")
CONTAINER_BEGIN = re.compile(r"\(")
ATTR_BEGIN = re.compile(r"\[")
ATTR_END = re.compile(r"\]")
BS_RE = re.compile(r":(?P<bs>[\d\.E-]+)")


def parse_attributes(tree_buf, start_pos, node):
    start_pos += 1
    end_ma = ATTR_END.search(tree_buf, start_pos)
    end_i = end_ma.start()
    for attr in tree_buf[start_pos:end_i].split(','):
        key, value = attr.split('=')
        key = key.replace('&', '').replace('!', '')
        node.add_attr(key, value)
    return end_i + 1


def parse_bs(tree_buf, start_pos, node):
    ma = BS_RE.match(tree_buf, start_pos)
    if ma is not None:
        node.bootstrap_support = ma.group('bs')
        start_pos += ma.end() - ma.start()
    return start_pos


def parse_node(tree_buf, start_pos):
    ma = SEQ_ID_RE.match(tree_buf, start_pos)
    seq_id = ma.group().replace("'", '')
    node = Node()
    node.seq_id = seq_id
    start_pos += ma.end() - ma.start()

    ma = ATTR_BEGIN.match(tree_buf, start_pos)
    if ma is not None:
        start_pos = parse_attributes(tree_buf, start_pos, node)

    start_pos = parse_bs(tree_buf, start_pos, node)

    if tree_buf[start_pos] == ',':
        start_pos += 1

    logging.debug(str.format("Match seq {}", node))
    return node, start_pos


def parse_container(tree_buf, start_pos):
    container = Container()
    start_pos += 1
    while start_pos < len(tree_buf):
        if SEQ_ID_RE.match(tree_buf, start_pos):
            node, start_pos = parse_node(tree_buf, start_pos)
            container.add_node(node)
            continue
        if CONTAINER_BEGIN.match(tree_buf, start_pos):
            child_container, start_pos = parse_container(tree_buf, start_pos)
            container.add_child(child_container)
            continue
        # find the brackets of the container
        if tree_buf[start_pos] == ')':
            start_pos += 1
            break
        logging.error(str.format("Unexpected at pos {}, {}", start_pos, tree_buf[start_pos:start_pos+30]))
        start_pos += 1

    if ATTR_BEGIN.match(tree_buf, start_pos):
        start_pos = parse_attributes(tree_buf, start_pos, container)
    if BS_RE.match(tree_buf, start_pos):
        start_pos = parse_bs(tree_buf, start_pos, container)
    if tree_buf[start_pos] == ',':
        start_pos += 1

    return container, start_pos


def parse_tree(tree_buf):
    containers = []
    i = 0
    while i < len(tree_buf):
        if tree_buf[i] == ';':
            logging.info("Found ';', the process of parsing tree is goinng to end")
            i += 1
            continue
        container, i = parse_container(tree_buf, i)
        containers.append(container)
    if len(containers) != 1:
        logging.error(str.format("More than one root containers were found({})", len(containers)))
        return None
    return containers[0]


def parse_tree_by_buf(fb):
    IN_TREE = False
    for line in fb.splitlines():
        if line.strip() == 'begin trees;':
            IN_TREE = True
            continue
        if line.strip() == 'end;':
            IN_TREE = False
            continue
        if IN_TREE:
            eles = line.split()
            assert len(eles) == 5
            assert eles[0] == 'tree'
            assert eles[2] == '='
            assert eles[3] == '[&R]'
            tree = Tree(eles[1])
            tree_buf = eles[4]
            tree.root_node = parse_tree(tree_buf)
            return tree


def parse_tree_by_file(file_path):
    with open(file_path, 'r') as f:
        return parse_tree_by_buf(f.read())


def _group_by_colors_recursively(container, colors_to_seq_ids):
    for child_node in container.nodes:
        color = child_node.attrs.color
        if color is not None:
            if color not in colors_to_seq_ids:
                colors_to_seq_ids[color] = []
            colors_to_seq_ids[color].append(child_node.seq_id)
    for child_container in container.children:
        _group_by_colors_recursively(child_container, colors_to_seq_ids)


def group_by_colors(tree):
    colors_to_seq_ids = {}
    _group_by_colors_recursively(tree.root_node, colors_to_seq_ids)
    return colors_to_seq_ids


def main_test():
    with open("../test_files/trees/pepr1_figtree.tree", 'r') as f :
        fb = f.read()
    tree = parse_tree_by_buf(fb)
    return group_by_colors(tree)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    ret = main_test()

    for key, value in ret.items():
        print(str.format("color {}, count {}", key, len(value)))
        for seq_id in value:
            print(seq_id)
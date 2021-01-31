# THIS FILE IS PART OF phytolrr.com PROJECT.
# Copyright 2019 phytolrr.com. All rights reserved.

def get_names_from_file(gene_file_path):
    names = []
    with open(gene_file_path, "r") as f:
        for line in f.readlines():
            for name in line.split():
                stripped_name = name.strip()
                if len(stripped_name) == 0:
                    continue
                else:
                    names.append(stripped_name)
    return names
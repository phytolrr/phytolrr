def nsites_entity_output(nsites):
    output = []
    for nsite in nsites:
        output.append({
            'offset': nsite.start_pos,
            'ntype': ['S', 'T'][nsite.ntype]
        })
    output.sort(key=lambda k: k['offset'])
    return output
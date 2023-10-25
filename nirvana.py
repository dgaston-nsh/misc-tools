import gzip
import json

def parseNirvana(file):
    is_header_line = True
    is_position_line = False
    is_gene_line = False

    gene_section_line = '],"genes":['
    end_line = ']}'

    data = dict()

    data['header'] = ''
    data['samples'] = []
    data['positions'] = []
    data['genes'] = []

    with gzip.open(file, 'rt') as f:
        data['position_count'] = 0
        data['gene_count'] = 0

        for line in f:
            trimmed_line = line.strip()
            if is_header_line:
                data['header'] = trimmed_line[10:-14]
                header_dict = json.loads(data['header'])
                data['samples'] = header_dict['samples']
                is_header_line = False
                is_position_line = True
                continue
            if trimmed_line == gene_section_line:
                is_gene_line = True
                is_position_line = False
                continue
            elif trimmed_line == end_line:
                break
            else:
                if is_position_line:
                    position_dict = json.loads(trimmed_line.rstrip(','))
                    data['positions'].append(position_dict)
                    data['position_count'] += 1
                if is_gene_line:
                    data['genes'].append(trimmed_line.rstrip(','))
                    data['gene_count'] += 1
    return data

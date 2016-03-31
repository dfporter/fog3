import re
import argparse
import sys
import config
import collections


def add_wb_name_to_gtf(lib):
    fname = lib['gtf_raw']
    gtf = {}
    header = "0\t1\t2\t3\t4\t5\t6\t7\t8\tgene_name"
    header += "\ttranscript_id\ttranscript_name\texon_number\tgene_id\tbiotype\n"
    gene_id_to_transcript_id = collections.defaultdict(list)
    gtf_lines = []
    with open(fname, 'r') as f:
        for line in f:
            gene_id = re.match('.*gene_id "([^"]+)";.*', line)
            if gene_id is None: print line
            gene_id = gene_id.group(1).partition('.')[0]
            gene_name = re.match('.*gene_name "([^"]+)";.*', line).group(1)
            biotype = re.match('.*biotype "([^"]+)";.*', line)
            biotype = biotype.group(1)
            transcript_name_m = re.match('.*transcript_name "([^"]+)";.*', line)
            transcript_id_m = re.match('.*transcript_id "([^"]+)";.*', line)
            exon_number_m = re.match('.*exon_number "([^"]+)";.*', line)
            if exon_number_m is not None:
                exon_number = exon_number_m.group(1)
            else: exon_number = '.'
            if transcript_name_m is not None:
                transcript_name = transcript_name_m.group(1)
            else: transcript_name = '.'
            if transcript_id_m is not None:
                transcript_id = transcript_id_m.group(1)
                gene_id_to_transcript_id[gene_id].append(transcript_id)
            else:
                transcript_id = '.'
            gtf_lines.append({
                    'gene_id': gene_id, 'line': line,
                    'transcript_name': transcript_name,
                    'transcript_id': transcript_id, 'exon_number': exon_number,
                    'biotype': biotype, 'gene_name': gene_name})
            if gene_name in gtf:
                if gene_id == gtf[gene_name]['gene_id']:
                    gtf[gene_name]['lines'].append(line)
                    # if gtf[gene_name]['transcript_name'] == '.'
                else:
                    print "Two gene IDs for the same gene name:",
                    print "%s: %s vs %s" % (gene_name, gene_id, gtf[gene_name]['gene_id'])
            else:
                gtf[gene_name] = {
                    'gene_id': gene_id, 'lines': [line],
                    'transcript_name': transcript_name,
                    'transcript_id': transcript_id, 'exon_number': exon_number,
                    'biotype': biotype, 'gene_name': gene_name}
    with open(lib['gtf'], 'w') as f:
        f.write(header)
        # for gene in gtf:
            # for line in gtf[gene]['lines']:
        for _line in gtf_lines:
            f.write(_line['line'].rstrip('\n') + '\t' + '\t'.join([
                _line['gene_name'], _line['transcript_id'],
                _line['transcript_name'], _line['exon_number'],
                _line['gene_id'], _line['biotype']
                ]) + '\n'
            )
                # f.write(line.rstrip('\n') + '\t' + '\t'.join([
                #     gtf[gene]['gene_name'], gtf[gene]['transcript_id'],
                #     gtf[gene]['transcript_name'], gtf[gene]['exon_number'],
                #     gtf[gene]['gene_id'], gtf[gene]['biotype']
                #     ]) + '\n'
                # )

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Converts a raw gtf file into a table with \
        a header and gene_id and biotype columns.''')
    parser.add_argument(
        '-c', '--config_ini', default='config.ini',
        hel'config.ini file.'
    )
    args = parser.parse_args()
    lib = config.config(args.config_ini)
    add_wb_name_to_gtf(lib)
